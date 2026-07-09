import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pingouin as pg
import starbars

parser = argparse.ArgumentParser(
    description="Generate manifold analysis figures."
)

parser.add_argument(
    "data_dir",
    type=Path,
    help="Path to the analysis directory containing the mouse folders."
)

args = parser.parse_args()

DATA_ROOT = args.data_dir

MICE = ["DG31", "DG32", "DG33", "DG35"]

EXPERIMENTS = [
    "250409",
    "250417",
    "250424",
    "250513",
]

LABELS = {
    "250409": "Env",
    "250417": "EU",
    "250424": "UU",
    "250513": "UI",
}

COLOURS = {
    "250409": "darkgoldenrod",
    "250417": "#5493A0",
    "250424": "#956480",
    "250513": "#458d50",
}



stability_points = {exp: [] for exp in EXPERIMENTS}
tortuosities = {exp: [] for exp in EXPERIMENTS}

for experiment in EXPERIMENTS:
    for mouse in MICE:

        folder = DATA_ROOT / mouse / experiment

        with open(folder / "manifold_results.json") as f:
            results = json.load(f)

        cka_score = np.loadtxt(folder / "cka_score.txt")
        shuffle_stabs = np.loadtxt(folder / "shuffle.txt")

        stability_points[experiment].append(float(results["t_star"]))
        tortuosities[experiment].append(float(results["tortuosity"]))




def repeated_measures_statistics(data):
    """Run repeated-measures ANOVA and Bonferroni-corrected post hoc tests."""

    df = pd.DataFrame(
        {k: v for k, v in data.items() if k != "250409"}
    )

    df["mouse"] = np.arange(len(df))

    df_long = df.melt(
        id_vars="mouse",
        var_name="experiment",
        value_name="value",
    )

    anova = pg.rm_anova(
        data=df_long,
        dv="value",
        within="experiment",
        subject="mouse",
        detailed=True,
    )

    posthoc = pg.pairwise_tests(
        data=df_long,
        dv="value",
        within="experiment",
        subject="mouse",
        padjust="bonf",
        effsize="cohen",
    )

    return anova, posthoc

def plot_metric(
    data,
    output_name,
    title,
    ylabel,
    ylim=None,
    outlier_idx=None,
):

    keep = np.ones(len(MICE), dtype=bool)

    if outlier_idx is not None:
        keep[outlier_idx] = False

    labels = list(data.keys())
    x = np.arange(len(labels))

    jitter = np.random.normal(0, 0.05, keep.size)

    fig, ax = plt.subplots(figsize=(5, 3))

    # Mean ± SEM
    for i, key in enumerate(labels):

        values = np.asarray(data[key])[keep]

        mean = values.mean()
        sem = values.std(ddof=1) / np.sqrt(values.size)

        ax.errorbar(
            i,
            mean,
            yerr=sem,
            fmt="o",
            markersize=8,
            color=COLOURS[key],
            ecolor=COLOURS[key],
            capsize=4,
            linewidth=1.5,
            zorder=4,
        )

    # Individual animals
    for i, key in enumerate(labels):

        values = np.asarray(data[key])

        ax.scatter(
            np.full(keep.sum(), i) + jitter[keep],
            values[keep],
            s=28,
            color="black",
            zorder=3,
        )

    # Connecting lines
    Y = np.column_stack([data[k] for k in labels])[keep, 1:]
    X = x[1:] + jitter[keep, None]

    for xi, yi in zip(X, Y):
        ax.plot(
            xi,
            yi,
            color="black",
            alpha=0.6,
            linewidth=0.9,
            zorder=2,
        )

    # Formatting
    ax.set_xticks(x)
    ax.set_xticklabels(
        [LABELS[k] for k in labels],
        rotation=30,
        ha="right",
    )

    ax.set_ylabel(ylabel)
    ax.set_title(title)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if ylim is not None:
        ax.set_ylim(*ylim)

    # Statistics
    anova, posthoc = repeated_measures_statistics(data)

    print(anova)
    print(posthoc)
    p1 = posthoc.query("A == '250417' and B == '250424'")["p_corr"].values[0]
    p2 = posthoc.query("A == '250417' and B == '250513'")["p_corr"].values[0]

    if output_name == "stability_points":

        starbars.draw_annotation(
            [
                (2, 3, p1),
                (2, 4, p2),
            ],
            ax=ax,
        )

    fig.tight_layout()

    fig.savefig(DATA_ROOT / f"{output_name}.svg")

    plt.close(fig)


plot_metric(
    tortuosities,
    output_name="tortuosity",
    title="Manifold Tortuosity",
    ylabel="Path length / displacement",
)

plot_metric(
    stability_points,
    output_name="stability_points",
    title="Stabilisation Time",
    ylabel="Laps after switch",
    ylim=(0, 40),
)