import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from manifold_analysis import manifold_analysis
from main import compute_binned_speed_and_occupancy, compute_binned_licks


MICE = [
    "DG19", "DG20", "DG21", "DG22",
    "DG23", "DG24", "DG25", "DG26",
    "DG27", "DG28", "DG29", "DG30",
]

EXPERIMENTS = [
    "EU",
    "NovelEnv",
    "UI",
    "UU",
]


def zscore_rows(A):
    """Z-score each row independently."""
    return (
        A - A.mean(axis=1, keepdims=True)
    ) / (
        A.std(axis=1, keepdims=True) + 1e-8
    )


def process_session(input_dir, output_dir, mouse, experiment):
    """
    Process one mouse/experiment session.

    Returns
    -------
    list
        Stability deviations from the switch point.
    """

    stabs = []

    behaviour_file = input_dir / experiment / mouse / "GRAB_behaviour_switch_export.mat"
    print
    if not behaviour_file.exists():
        return stabs

    session_dir = input_dir / experiment / mouse 

    with h5py.File(behaviour_file, "r") as beh:

        vel = np.array(beh.get("VdsSm"))[0]
        pos = np.array(beh.get("Y_ds"))[0]
        licks = np.array(beh.get("licks_ds"))[0]
        switch = int(np.array(beh.get("switch_lap"))[0, 0])

    lap_info = pd.read_csv(
        session_dir / "lap_level_switch_export.csv"
    )

    n_samples = len(pos)

    lap_start_idx = lap_info["lap_start_idx"].values

    lap_indices = np.searchsorted(
        lap_start_idx,
        np.arange(n_samples),
        side="right",
    ) - 1

    binned_speeds, occ = compute_binned_speed_and_occupancy(
        pos,
        vel,
        lap_indices,
    )

    binned_licks = compute_binned_licks(
        pos,
        licks,
        lap_indices,
    )

    V = pd.DataFrame(data=binned_speeds).T.corr().values

    L = (
        pd.DataFrame(data=binned_licks)
        .T
        .corr()
        .fillna(0)
        .values
    )

    output_session = output_dir / mouse / experiment
    output_session.mkdir(
        parents=True,
        exist_ok=True,
    )

    sns.heatmap(V, cmap="magma")
    plt.savefig(output_session / "speeds.svg")
    plt.close()

    sns.heatmap(L, cmap="magma")
    plt.savefig(output_session / "licks.svg")
    plt.close()

    Xn = zscore_rows(binned_licks)
    Yn = zscore_rows(binned_speeds)

    Z = np.concatenate(
        [Xn, Yn],
        axis=1,
    )

    Z = Z / np.linalg.norm(
        Z,
        axis=1,
        keepdims=True,
    )

    C = Z @ Z.T
    stab = (
        manifold_analysis(
            C,
            switch,
            mouse,
            experiment,
        )
        + switch
    )

    if stab != switch:
        stabs.append(stab - switch)

    return stabs


def main():

    parser = argparse.ArgumentParser(
        description="Run manifold analysis on GRAB behavioural data."
    )

    parser.add_argument(
        "input_dir",
        type=Path,
        help="Input directory containing mouse/experiment folders.",
    )

    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory where analysis outputs will be saved.",
    )

    args = parser.parse_args()

    stabs = []

    for experiment in EXPERIMENTS:
        for mouse in MICE:

            print(f"Processing {experiment}/{mouse}")

            stabs.extend(
                process_session(
                    args.input_dir,
                    args.output_dir,
                    mouse,
                    experiment,
                )
            )

    np.savetxt(
        args.output_dir / "stabs.txt",
        stabs,
    )


if __name__ == "__main__":
    main()