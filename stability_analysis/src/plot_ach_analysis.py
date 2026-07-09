import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.optimize import curve_fit
from scipy.stats import spearmanr


EXPERIMENT_FILES = {
    "EU": "pooledSwitch_0_10cm_s_weightedBins2(N2B).csv",
    "Novelty": "pooledSwitch_0_10cm_s_weightedBins2(Teleport1).csv",
    "UI": "pooledSwitch_0_10cm_s_weightedBins2(B2N).csv",
    "UU": "pooledSwitch_0_10cm_s_weightedBins2(N2F).csv",
}

COLOURS = {
    "EU": "cyan",
    "UU": "magenta",
    "UI": "darkblue",
    "Novelty": "darkgoldenrod",
}


def single_exponential(x, A, B, C):
    """
    Single exponential model.

    y = A * exp(B*x) + C
    """
    return A * np.exp(B * x) + C


def load_dff_values(data_dir):
    """
    Extract late post-switch dFF values from experiment files.

    Parameters
    ----------
    data_dir : Path
        Directory containing experiment CSV files.

    Returns
    -------
    means : ndarray
        Mean late post-switch dFF values.
    experiments : ndarray
        Experiment identity for each mouse.
    """

    means = []
    experiments = []

    for experiment, filename in EXPERIMENT_FILES.items():

        file_path = data_dir / filename

        if not file_path.exists():
            print(f"Skipping missing file: {file_path}")
            continue

        df = pd.read_csv(
            file_path,
            header=None,
        )

        current_mouse = None

        for _, row in df.iterrows():

            identifier = str(row.iloc[0]).strip()

            if identifier.startswith("DG"):

                current_mouse = identifier

            elif identifier == "dFF" and current_mouse:

                values = pd.to_numeric(
                    row.iloc[1:],
                    errors="coerce",
                )

                late_post = (
                    values.iloc[22:27]
                    .dropna()
                    .to_numpy()
                )

                if late_post.size:
                    means.append(
                        np.mean(late_post)
                    )

                    experiments.append(
                        experiment
                    )

    return (
        np.asarray(means),
        np.asarray(experiments),
    )


def fit_exponential(x, y):
    """
    Fit a single exponential relationship.
    """

    initial_guess = [
        y.max() - y.min(),
        -1,
        y.min(),
    ]

    parameters, _ = curve_fit(
        single_exponential,
        x,
        y,
        p0=initial_guess,
        maxfev=10000,
    )

    return parameters


def plot_relationship(
    means,
    stabilities,
    experiments,
    parameters,
    output_file,
):
    """
    Plot dFF-stability relationship.
    """

    fig, ax = plt.subplots(
        figsize=(5, 5)
    )

    for experiment in EXPERIMENT_FILES:

        mask = experiments == experiment

        if not np.any(mask):
            continue

        if experiment == "Novelty":

            ax.scatter(
                means[mask],
                stabilities[mask],
                facecolors="white",
                edgecolors="black",
                linewidths=1.2,
                s=50,
                label="Novelty (excluded)",
            )

        else:

            ax.scatter(
                means[mask],
                stabilities[mask],
                color=COLOURS[experiment],
                s=50,
                label=experiment,
            )


    x_fit = np.linspace(
        means.min(),
        means.max(),
        300,
    )

    ax.plot(
        x_fit,
        single_exponential(
            x_fit,
            *parameters,
        ),
        color="black",
        linewidth=2,
        label="Exponential fit",
    )

    ax.set_xlabel(
        "Mean late post-switch - baseline dFF"
    )

    ax.set_ylabel(
        "Stability"
    )


    legend_order = [
        "EU",
        "UU",
        "UI",
        "Novelty (excluded)",
        "Exponential fit",
    ]

    handles, labels = ax.get_legend_handles_labels()

    ordered = [
        (h, l)
        for l in legend_order
        for h, label in zip(handles, labels)
        if label == l
    ]

    if ordered:
        handles, labels = zip(*ordered)

        ax.legend(
            handles,
            labels,
            title="Experiment",
        )


    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)

    ax.tick_params(
        width=1.5
    )

    fig.tight_layout()

    fig.savefig(
        output_file,
        bbox_inches="tight",
    )

    plt.close(fig)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Analyse relationship between "
            "post-switch activity and manifold stability."
        )
    )

    parser.add_argument(
        "data_dir",
        type=Path,
        help="Directory containing dFF CSV files.",
    )

    parser.add_argument(
        "stability_file",
        type=Path,
        help="Path to stabs.txt.",
    )

    parser.add_argument(
        "output_file",
        type=Path,
        help="Output figure filename.",
    )

    args = parser.parse_args()


    means, experiments = load_dff_values(
        args.data_dir
    )

    stabilities = np.loadtxt(
        args.stability_file
    )

    if len(means) != len(stabilities):
        raise ValueError(
            "Number of dFF measurements does not "
            "match number of stability values."
        )


    analysis_mask = (
        experiments != "Novelty"
    )

    x = means[analysis_mask]
    y = stabilities[analysis_mask]


    rho, p_value = spearmanr(
        x,
        y,
    )

    print(
        f"Spearman rho = {rho:.3f}"
    )

    print(
        f"p = {p_value:.4f}"
    )


    parameters = fit_exponential(
        x,
        y,
    )

    print(
        "Exponential parameters:",
        parameters,
    )


    plot_relationship(
        means,
        stabilities,
        experiments,
        parameters,
        args.output_file,
    )


if __name__ == "__main__":
    main()