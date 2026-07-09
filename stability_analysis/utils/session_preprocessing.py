import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat


MICE = ["DG31", "DG32", "DG33", "DG35", "DG36"]
EXPERIMENTS = ["250409", "250417", "250424", "250513"]


def resample_to_imaging_frames(signal, behaviour_time, n_frames):
    """
    Resample a behavioural signal to match the number of imaging frames.

    Parameters
    ----------
    signal : ndarray
        Behavioural signal sampled at behavioural timestamps.
    behaviour_time : ndarray
        Behavioural timestamps.
    n_frames : int
        Number of imaging frames.

    Returns
    -------
    ndarray
        Behavioural signal aligned to imaging frames.
    """
    behaviour_time = behaviour_time - behaviour_time[0]
    frame_times = np.linspace(0, behaviour_time[-1], n_frames)

    indices = np.searchsorted(behaviour_time, frame_times, side="left")
    indices = np.clip(indices, 0, len(signal) - 1)

    return signal[indices] - np.min(signal)


def process_session(session_dir):
    """Process one mouse/session pair."""

    mouse = session_dir.name
    experiment = session_dir.parent.name

    pre_behaviour = loadmat(session_dir / f"{mouse}_{experiment}_pre_tyv.mat")
    post_behaviour = loadmat(session_dir / f"{mouse}_{experiment}_post_tyv.mat")

    pre_cells = pd.read_csv(session_dir / "output_data_pre" / "Fc3_cleaned.csv")
    post_cells = pd.read_csv(session_dir / "output_data_post" / "Fc3_cleaned.csv")

    coords_pre = resample_to_imaging_frames(
        np.asarray(pre_behaviour["y_pre"]).squeeze(),
        np.asarray(pre_behaviour["t_pre"]).squeeze(),
        len(pre_cells),
    )

    coords_post = resample_to_imaging_frames(
        np.asarray(post_behaviour["y_post"]).squeeze(),
        np.asarray(post_behaviour["t_post"]).squeeze(),
        len(post_cells),
    )

    speed_pre = resample_to_imaging_frames(
        np.asarray(pre_behaviour["v_pre"]).squeeze(),
        np.asarray(pre_behaviour["t_pre"]).squeeze(),
        len(pre_cells),
    )

    speed_post = resample_to_imaging_frames(
        np.asarray(post_behaviour["v_post"]).squeeze(),
        np.asarray(post_behaviour["t_post"]).squeeze(),
        len(post_cells),
    )

    licks_pre = resample_to_imaging_frames(
        np.asarray(pre_behaviour["licks_pre"]).squeeze(),
        np.asarray(pre_behaviour["t_pre"]).squeeze(),
        len(pre_cells),
    )

    licks_post = resample_to_imaging_frames(
        np.asarray(post_behaviour["licks_post"]).squeeze(),
        np.asarray(post_behaviour["t_post"]).squeeze(),
        len(post_cells),
    )

    np.savetxt(
        session_dir / "coords.csv",
        np.concatenate([coords_pre, coords_post]),
    )

    np.savetxt(
        session_dir / "speeds.csv",
        np.concatenate([speed_pre, speed_post]),
    )

    np.savetxt(
        session_dir / "licks.csv",
        np.concatenate([licks_pre, licks_post]),
    )

    np.savetxt(
        session_dir / "switch_point.txt",
        [len(coords_pre)],
        fmt="%d",
    )


def main():
    parser = argparse.ArgumentParser(
        description="Prepare behavioural data for manifold analysis."
    )

    parser.add_argument(
        "data_root",
        type=Path,
        help="Root directory containing experiment folders.",
    )

    args = parser.parse_args()

    for experiment in EXPERIMENTS:
        for mouse in MICE:
            session_dir = args.data_root / experiment / mouse

            if not session_dir.exists():
                print(f"Skipping missing session: {session_dir}")
                continue

            print(f"Processing {experiment} / {mouse}")
            process_session(session_dir)


if __name__ == "__main__":
    main()