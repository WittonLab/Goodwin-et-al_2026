import numpy as np
import matplotlib.pyplot as plt 
import json
import pyhomogeneity as hg
import pymannkendall as mk
import os
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from sklearn.neighbors import NearestNeighbors
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import WhiteKernel, Matern
from scipy.ndimage import gaussian_filter1d
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def plot_manifold(manifold, laps, ax, switch, y_mean):
    """3D visualisation of the manifold."""
    sc = ax.scatter(
        manifold[:, 0],
        manifold[:, 1],
        manifold[:, 2],
        c=laps,
        cmap="magma",
        s=12,
        alpha=0.9,
    )
    ax.grid(True, color="0.9", linewidth=0.3)
    x = y_mean[:, 0]
    y = y_mean[:, 1]
    z = y_mean[:, 2]
    t = np.arange(len(x))
    points = np.array([x, y, z]).T.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = Line3DCollection(
    segments,
    cmap="magma",
    linewidth=2
    )
    lc.set_array(t[:-1])
    ax.add_collection3d(lc)

    ax.set_xlabel("Psi1", fontsize=9)
    ax.set_ylabel("Psi2", fontsize=9)
    ax.set_zlabel("Psi3", labelpad=12, fontsize=9)
    

    cbar = plt.colorbar(sc, ax=ax, pad=0.15, label="Laps")
    cbar.ax.invert_yaxis()
    cbar.ax.axhline(switch, color="fuchsia", label="Switch")

def fix_diffusion_signs(coords, ref_idx=0):
    coords = coords.copy()
    n_components = coords.shape[1]
    for k in range(n_components):
        if coords[ref_idx, k] < 0:
            coords[:, k] *= -1
    return coords

def diffusion_map(X, k, sigma=2.0, n_components=3, skip_first=True):
    """
    Compute a diffusion map from a correlation matrix.

    Parameters
    ----------
    X : np.ndarray
       A correlation matrix
    k : int
       Number of neighbours
    sigma : float
       Kernel density
    n_components : int
       Number of dimensions in the output manifold
    skip_first : bool
       Whether or not to skip the first (trivial) eigenvector
    """

    # Build a symmetric weighted distance graph
    nbrs = NearestNeighbors(n_neighbors=k).fit(X)
    knn_graph = nbrs.kneighbors_graph(X, mode="distance")
    W = 0.5 * (knn_graph + knn_graph.T)
    W = W.toarray()
    K = np.exp(-(W**2) / (2 * sigma**2))
    K[W == 0] = 0 

    plt.axis("off")
    plt.savefig("network.svg", bbox_inches="tight")
    plt.close()

    # Build a markovian transition matrix
    q = np.sum(K, axis=1)
    Q_inv = np.diag(1.0 / q)
    K_tilde = Q_inv @ K @ Q_inv 
    d = np.sum(K_tilde, axis=1)
    d_sqrt_inv = np.diag(1.0 / np.sqrt(d))
    P = (np.diag(1.0 / d) @ K_tilde)
    
    # Derive diffusion coordinates
    M_s = d_sqrt_inv @ K_tilde @ d_sqrt_inv 
    eigvals, eigvecs_s = eigh(M_s)
    idx = np.argsort(eigvals)[::-1]
    eigvals = np.clip(eigvals[idx], 0, 1)
    eigvecs_s = eigvecs_s[:, idx]
    psi = d_sqrt_inv @ eigvecs_s
    start = 1 if skip_first else 0
    if n_components is None:
        n_components = psi.shape[1] - start
    lam = eigvals[start:start+n_components]
    gaps = np.abs(np.diff(eigvals[1:]))
    Psi = psi[:, start:start+n_components]
    diffusion_coords = Psi * (lam[None, :])
    # This step ensures sign consistency across different experiments
    # Just makes plots a little bit easier to immediately compare
    Phi = fix_diffusion_signs(diffusion_coords) 
    return lam, Phi, P

def compute_l_norm(Phi):
    # Compute total step length
    dPhi = Phi[1:] - Phi[:-1]
    step_lengths = np.linalg.norm(dPhi, axis=1)
    L_total = np.sum(step_lengths)

    # Normalize by embedding energy
    Phi_center = Phi - Phi.mean(axis=0, keepdims=True)
    scale = np.mean(np.sum(Phi_center**2, axis=1))  # mean squared radius
    L_norm = L_total / (np.sqrt(scale) + 1e-10)
    return L_norm, L_total, scale

def compute_manifold(cor_mat, switch, n_dims=2):
    """Compute import manifold features"""
    num_laps = cor_mat.shape[0]
    k = int(np.sqrt(num_laps) * 2)
    lam, Phi, P = diffusion_map(cor_mat, k)
    L_norm_pre, L_total_pre, _ = compute_l_norm(Phi[:switch])
    L_norm_post, L_total_post, _ = compute_l_norm(Phi[switch:])
    manifold = Phi[:, :n_dims+1]  
    spectral_gap = lam[0] - lam[1]
    return manifold, L_norm_pre, L_norm_post, spectral_gap, L_total_pre, L_total_post

def split_manifold(manifold, switch):
    """Segment manifold into pre and post switch"""
    pre_manifold = manifold[:switch]
    post_manifold = manifold[switch:]
    return pre_manifold, post_manifold

def compute_convergence(post_manifold, fraction=0.15, sigma=1.5):
    """Compute the convergence of the post-switch manifold from some initial centroid."""
    centroid_point = int(np.floor(fraction * post_manifold.shape[0]))
    initial_centroid = np.median(
        post_manifold[:centroid_point],
        axis=0
    )
    convergence_dist = np.linalg.norm(post_manifold - initial_centroid, axis=1)
    return convergence_dist

def qscale(x, lo=5, hi=95, eps=1e-9):
    a, b = np.percentile(x, [lo, hi])
    return np.clip((x - a) / (b - a + eps), 0, 1)

def rqa(X, mouse, experiment, upsample_factor=8):
    """
    Run recurrence quantification analysis and estimate stability onset.

    Parameters
    ----------
    X : np.ndarray
        Shape (T, features)
    mouse : str
        Output folder label.
    experiment : str
        Output folder label.
    upsample_factor : int, default=8
        Temporal interpolation factor.

    Returns
    -------
    int
        Estimated stability onset in original time units.
    """
    T_orig = X.shape[0]
    W = max(3, round(T_orig / 6))

    X_interp = interpolate_timeseries(X, upsample_factor)
    T = X_interp.shape[0]
    W_scaled = int(W * upsample_factor)

    R = build_recurrence_matrix(X_interp, target_density=10)

    rr, det, lam = sliding_rqa_metrics(R, W_scaled)

    score = build_stability_score(
        rr=rr,
        det=det,
        lam=lam,
        upsample_factor=upsample_factor,
    )

    T_stable_upsampled = detect_stability_onset(
        score,
        upsample_factor=upsample_factor,
    )

    T_stable_real = int(T_stable_upsampled / upsample_factor)

    save_recurrence_plot(
        R=R,
        T=T,
        upsample_factor=upsample_factor,
        T_stable_upsampled=T_stable_upsampled,
        T_stable_real=T_stable_real,
        mouse=mouse,
        experiment=experiment,
    )

    return T_stable_real

def interpolate_timeseries(X, upsample_factor):
    """Cubic spline temporal upsampling."""
    T_orig = X.shape[0]

    t_orig = np.arange(T_orig)
    t_new = np.linspace(0, T_orig - 1, num=T_orig * upsample_factor)

    cs = CubicSpline(t_orig, X, axis=0)
    return cs(t_new)

def pairwise_distance_matrix(X):
    """Euclidean pairwise distance matrix."""
    s = np.sum(X * X, axis=1)
    D2 = s[:, None] + s[None, :] - 2.0 * (X @ X.T)
    np.maximum(D2, 0.0, out=D2)
    return np.sqrt(D2)


def recurrence_threshold(R, target_density=10):
    """Percentile threshold excluding Theiler window."""
    T = R.shape[0]
    theiler = max(1, int(0.05 * T))

    np.fill_diagonal(R, np.inf)

    i, j = np.triu_indices(T, k=theiler + 1)
    dists = R[i, j]

    return np.percentile(dists, target_density)


def build_recurrence_matrix(X, target_density=10):
    """Binary recurrence matrix."""
    R = pairwise_distance_matrix(X)
    tol = recurrence_threshold(R, target_density)
    R = (R < tol).astype(float)
    np.fill_diagonal(R, 0)
    return R

def run_lengths(binary_vec):
    """Lengths of contiguous 1-runs."""
    d = np.diff(np.pad(binary_vec.astype(np.int8), (1, 1)))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    return ends - starts

def laminarity_metric(R_sub):
    """LAM from vertical lines length >= 2."""
    lam_sum = 0

    for j in range(R_sub.shape[1]):
        lengths = run_lengths(R_sub[:, j])
        lam_sum += lengths[lengths >= 2].sum()

    total = np.sum(R_sub)
    return lam_sum / total if total > 0 else 0.0

def determinism_metric(R_sub):
    """DET from diagonal lines length >= 2."""
    n = R_sub.shape[0]
    offsets = list(range(-(n - 2), 0)) + list(range(1, n - 1))

    diag_sum = 0
    for k in offsets:
        diag = np.diagonal(R_sub, offset=k)
        lengths = run_lengths(diag)
        diag_sum += lengths[lengths >= 2].sum()

    total = np.sum(R_sub)
    return diag_sum / total if total > 0 else 0.0

def recurrence_rate(R_sub):
    """RR."""
    return np.sum(R_sub) / (R_sub.shape[0] ** 2)

def sliding_rqa_metrics(R, window):
    """Compute RR, DET, LAM over sliding windows."""
    rr, det, lam = [], [], []

    for t in range(0, R.shape[0] - window):
        R_sub = R[t:t + window, t:t + window]

        total = np.sum(R_sub)
        rr_val = recurrence_rate(R_sub)

        if total == 0:
            rr.append(rr_val)
            det.append(0.0)
            lam.append(0.0)
            continue

        rr.append(rr_val)
        det.append(determinism_metric(R_sub))
        lam.append(laminarity_metric(R_sub))

    return np.array(rr), np.array(det), np.array(lam)

def build_stability_score(rr, det, lam, upsample_factor):
    """Normalized composite score of RQA metrics."""
    warmup = 3 * upsample_factor

    rr_n = qscale(rr[warmup:])
    det_n = qscale(det[warmup:])
    lam_n = qscale(lam[warmup:])

    structure = (det_n + lam_n) / 2
    score = (rr_n + structure) / 2

    return gaussian_filter1d(score, 1.5 * upsample_factor)


def detect_stability_onset(score, upsample_factor, alpha=0.05, max_iter=10):
    """
    Iterative Pettitt + Mann-Kendall detection of stability onset.
    """
    S = np.asarray(score).ravel()

    start = 0
    count = 1

    for _ in range(max_iter):
        segment = S[start:]

        if len(segment) < (3 * upsample_factor):
            
            break

        pettitt = hg.pettitt_test(segment)

        if pettitt.p >= (alpha / count):
            break

        candidate = start + pettitt.cp
        tail = S[candidate:]

        mk_result = mk.hamed_rao_modification_test(tail)

        if mk_result.trend == "no trend":
            start = candidate
            break

        start = candidate
        count += 1

    return start

def save_recurrence_plot(
    R,
    T,
    upsample_factor,
    T_stable_upsampled,
    T_stable_real,
    mouse,
    experiment,
):
    """Save recurrence matrix figure."""
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)

    ax.imshow(R, cmap="magma", interpolation="nearest", aspect="equal")

    ax.set_xlabel("Post-Switch Lap")
    ax.set_ylabel("Post-Switch Lap")
    ax.set_title("Recurrence Plot of Post-Switch Manifold")

    n_ticks = 6
    tick_pos = np.linspace(0, T - 1, n_ticks).astype(int)
    tick_lab = np.round(tick_pos / upsample_factor).astype(int)

    ax.set_xticks(tick_pos)
    ax.set_yticks(tick_pos)
    ax.set_xticklabels(tick_lab)
    ax.set_yticklabels(tick_lab)

    ax.axvline(T_stable_upsampled, color="black", lw=2.2, ls="--", alpha=0.6)
    ax.axvline(
        T_stable_upsampled,
        color="white",
        lw=1.2,
        ls="--",
        label=f"Stability Onset = {T_stable_real}",
    )

    ax.axhline(T_stable_upsampled, color="black", lw=2.2, ls="--", alpha=0.6)
    ax.axhline(T_stable_upsampled, color="white", lw=1.2, ls="--")

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.tick_params(length=2, width=0.8, labelsize=8)
    ax.legend(framealpha=0.7)

    fig.tight_layout(pad=0.5)

    outdir = f"AB2026_RSC/{mouse}/{experiment}"
    os.makedirs(outdir, exist_ok=True)

    plt.savefig(f"{outdir}/recurrence_matrix.svg")
    plt.close()

def save_metrics(post_manifold,
                 t_star,   
                 mouse, 
                 experiment, 
                 L_norm_post, 
                 L_norm_pre, 
                 spectral_gap, 
                 L_total_pre, 
                 L_total_post):

    post_steps = np.linalg.norm(np.diff(post_manifold, axis=0), axis=1)
    path_length = np.sum(post_steps)
    displacement = np.linalg.norm(post_manifold[-1] - post_manifold[0])
    tortuosity = path_length / displacement

    metrics =  {
        "relative_drift_speed_pre": str(L_norm_pre),
        "drift_speed_pre": str(L_total_pre),
        "relative_drift_speed_post": str(L_norm_post),
        "drift_speed_post": str(L_total_post),
        "spectral_gap": str(spectral_gap),
        "tortuosity": str(tortuosity),
        "t_star": str(t_star),
        "t_star_percent": t_star / post_manifold.shape[0],
    }
    with open(f'AB2026_RSC/{mouse}/{experiment}/manifold_results.json', 'w+') as fp:
        json.dump(metrics, fp, indent=2)

def manifold_analysis(cor_mat, switch, mouse, experiment):
    manifold, L_norm_pre, L_norm_post, spectral_gap, L_total_pre, L_total_post = compute_manifold(cor_mat, switch)
    pre_manifold, post_manifold = split_manifold(manifold, switch)
    kernel = 1.0 * Matern(length_scale=20.0, nu=1.5) + WhiteKernel(noise_level=1.0)
    gpr = GaussianProcessRegressor(kernel=kernel,
                                   optimizer='fmin_l_bfgs_b',
                                   alpha=0,
                                   normalize_y=True,
                                   n_restarts_optimizer=10)
    x = np.arange(cor_mat.shape[0]).reshape(-1, 1)
    y = manifold
    gpr.fit(x, y)
    x_smooth = np.linspace(x.min(), x.max(), 100).reshape(-1, 1)
    y_mean, y_std = gpr.predict(x_smooth, return_std=True)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    plot_manifold(manifold, np.arange(cor_mat.shape[0]), ax, switch, y_mean)
    plt.tight_layout()
    plt.savefig(f"AB2026_RSC/{mouse}/{experiment}/manifold.svg")
    plt.close()
    t_star = rqa(post_manifold, mouse, experiment)
    save_metrics(post_manifold, t_star, mouse, experiment, L_norm_post, L_norm_pre, spectral_gap, L_total_pre, L_total_post)
    return t_star


    

    