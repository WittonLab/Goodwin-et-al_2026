"""
Compute place-cell and behavioural correlation matrices and estimate
manifold stabilisation times. 

Usage
-----
python analyse_manifold.py DATA_ROOT OUTPUT_ROOT

Example
-------
python analyse_manifold.py RC_DGoodwin figures
"""

from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.ndimage import gaussian_filter

from manifold_analysis import manifold_analysis

MICE = [
    "DG31",
    "DG32",
    "DG33",
    "DG35",
]

EXPERIMENTS = [
    "250409",
    "250417",
    "250424",
    "250513",
]


EXPERIMENT_LABELS = {
    "250409": "Novel Environment",
    "250417": "Unexpected Uncertainty",
    "250424": "Expected Uncertainty",
    "250513": "Uncertainty Interaction",
}

def generate_lap_indices(lap_starts, total_length):
    """
    Generate sample indices belonging to each lap.
    """

    lap_starts = lap_starts[:-1]

    return [
        np.arange(
            lap_starts[i],
            lap_starts[i + 1]
            if i < len(lap_starts) - 1
            else total_length,
        )
        for i in range(len(lap_starts))
    ]


def compute_lap_ratemaps(
    calcium,
    coordinates,
    lap_indices,
    bin_size=5,
    sigma=2,
):
    """
    Compute spatial firing maps for every lap.
    """

    track_start = np.min(coordinates)
    track_end = np.max(coordinates)

    n_bins = int(
        np.round(
            (track_end - track_start)
            / bin_size
        )
    )

    bins = np.linspace(
        track_start,
        track_start + n_bins * bin_size,
        n_bins + 1,
    )

    n_cells = calcium.shape[1]

    ratemaps = []

    for lap in lap_indices:

        pos = coordinates[lap]
        traces = calcium[lap]

        occupancy, _ = np.histogram(
            pos,
            bins=bins,
        )

        occupancy = occupancy.astype(float)

        occupancy = gaussian_filter(
            occupancy,
            sigma=sigma,
            mode="nearest",
        )

        occupancy[occupancy == 0] = np.nan

        lap_map = np.zeros(
            (n_bins, n_cells)
        )

        for cell in range(n_cells):

            activity, _ = np.histogram(
                pos,
                bins=bins,
                weights=traces[:, cell],
            )

            activity = gaussian_filter(
                activity.astype(float),
                sigma=sigma,
                mode="nearest",
            )

            lap_map[:, cell] = np.nan_to_num(
                activity / occupancy,
                nan=0,
            )

        ratemaps.append(lap_map)

    return np.stack(ratemaps)


def normalise_ratemaps(ratemaps):
    """
    Normalise each cell by maximum activity.
    """

    peaks = ratemaps.max(axis=(0, 1))

    peaks[peaks == 0] = 1

    return ratemaps / peaks[None, None, :]


def compute_session_correlation(ratemaps):
    """
    Compute lap-by-lap population vector correlation.
    """

    flattened = ratemaps.reshape(
        ratemaps.shape[0],
        -1,
    )

    return np.nan_to_num(
        np.corrcoef(flattened)
    )


def compute_binned_speed_and_occupancy(coordinates, speeds, lap_idx, bin_size=5, min_occupancy=3):
    """Compute binned speeds and occupancy per lap."""
    track_start, track_end = np.min(coordinates), np.max(coordinates)
    num_bins = int(np.round((track_end - track_start) / bin_size))
    bin_edges = np.linspace(track_start, track_start + num_bins * bin_size, num_bins + 1)
    
    bin_idx = np.digitize(coordinates, bin_edges) - 1
    valid = (bin_idx >= 0) & (bin_idx < num_bins)

    num_laps = lap_idx.max() + 1
    speed_bin = np.zeros((num_laps, num_bins))
    occupancy = np.zeros((num_laps, num_bins))
    for t in range(len(coordinates)):
        if not valid[t]:
            continue
        l = lap_idx[t]
        b = bin_idx[t]
        speed_bin[l, b] += speeds[t]
        occupancy[l, b] += 1
    return speed_bin, occupancy

def compute_binned_licks(coordinates, licks, lap_idx, bin_size=5):

    track_start, track_end = np.min(coordinates), np.max(coordinates)

    bin_edges = np.arange(track_start, track_end + bin_size, bin_size)
    num_bins = len(bin_edges) - 1

    bin_idx = np.digitize(coordinates, bin_edges) - 1
    valid = (bin_idx >= 0) & (bin_idx < num_bins)

    num_laps = int(lap_idx.max()) + 1

    licks_bin = np.zeros((num_laps, num_bins))
    occupancy = np.zeros((num_laps, num_bins))

    np.add.at(licks_bin, (lap_idx[valid], bin_idx[valid]), licks[valid])
    np.add.at(occupancy, (lap_idx[valid], bin_idx[valid]), 1)
    return licks_bin


def zscore_rows(matrix):

    return (
        matrix - matrix.mean(axis=1, keepdims=True)
    ) / (
        matrix.std(axis=1, keepdims=True)
        + 1e-8
    )

def plot_correlation(
    matrix,
    output,
    title,
    switch,
    stability,
):

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    plt.figure(figsize=(5, 4))

    ax = sns.heatmap(
        matrix,
        cmap="magma",
        rasterized=True,
    )

    ax.axvline(
        switch,
        color="#00BFC4",
        linestyle="--",
        linewidth=2,
    )

    ax.axhline(
        switch,
        color="#00BFC4",
        linestyle="--",
        linewidth=2,
    )

    ax.axvline(
        stability,
        color="#00FF00",
        linestyle="--",
        linewidth=2,
    )

    ax.axhline(
        stability,
        color="#00FF00",
        linestyle="--",
        linewidth=2,
    )

    plt.title(title)

    plt.tight_layout()

    plt.savefig(
        output,
        format="svg",
    )

    plt.close()

def analyse_session(
    session_dir,
    output_root,
):

    coordinates = np.loadtxt(
        session_dir / "coords.csv"
    )

    speeds = np.loadtxt(
        session_dir / "speeds.csv"
    )

    licks = np.loadtxt(
        session_dir / "licks.csv"
    )

    laps = np.loadtxt(
        session_dir / "laps.csv",
        dtype=int,
    )

    lap_starts = np.loadtxt(
        session_dir / "lap_starts.csv",
        dtype=int,
    )

    switch_frame = int(
        np.loadtxt(
            session_dir / "switch_point.txt"
        )
    )

    calcium = pd.read_csv(
        session_dir / "Fc3_cleaned.csv"
    ).values


    experiment = session_dir.parent.name
    mouse = session_dir.name


    switch_lap = laps[switch_frame]


    if experiment != "250409":
        switch_lap += 1


    lap_indices = generate_lap_indices(
        lap_starts,
        len(coordinates),
    )


    ratemaps = compute_lap_ratemaps(
        calcium,
        coordinates,
        lap_indices,
    )

    ratemaps = normalise_ratemaps(
        ratemaps
    )

    pvc = compute_session_correlation(
        ratemaps
    )


    speed_matrix, _ = (
        compute_binned_speed_and_occupancy(
            coordinates,
            speeds,
            laps,
        )
    )

    lick_matrix = compute_binned_licks(
        coordinates,
        licks,
        laps,
    )


    combined = np.concatenate(
        [
            zscore_rows(lick_matrix),
            zscore_rows(speed_matrix),
        ],
        axis=1,
    )

    combined /= np.linalg.norm(
        combined,
        axis=1,
        keepdims=True,
    )


    behavioural_similarity = combined @ combined.T

    # To analyse behaviour, pass behvaioural_similarity as the first argument to manifold_analysis.
    # Make similar amendments to other functions that take "pvc" as an argument.
    stability = (
        manifold_analysis(
            pvc,
            switch_lap,
            mouse,
            experiment,
        )
        + switch_lap
    )


    output = (
        output_root /
        mouse /
        experiment /
        "pvc_matrix.svg"
    )


    plot_correlation(
        pvc,
        output,
        EXPERIMENT_LABELS[experiment],
        switch_lap,
        stability,
    )

def main():

    parser = argparse.ArgumentParser(
        description="Analyse manifold stability."
    )

    parser.add_argument(
        "data_root",
        type=Path,
        help="Input dataset directory.",
    )

    parser.add_argument(
        "output_root",
        type=Path,
        help="Output directory.",
    )

    args = parser.parse_args()


    for experiment in EXPERIMENTS:

        for mouse in MICE:

            session = (
                args.data_root /
                experiment /
                mouse
            )

            if not session.exists():
                print(
                    f"Skipping missing session: {session}"
                )
                continue

            print(
                f"Processing {experiment}/{mouse}"
            )

            analyse_session(
                session,
                args.output_root,
            )


if __name__ == "__main__":
    main()