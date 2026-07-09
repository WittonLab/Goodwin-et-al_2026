import argparse
from pathlib import Path

import numpy as np


MICE = ["DG31", "DG32", "DG33", "DG35", "DG36"]
EXPERIMENTS = ["250409", "250417", "250424", "250513"]


def compute_laps(coords, track_length=280):
    """
    Compute lap numbers and lap start indices from position data.

    Parameters
    ----------
    coords : ndarray
        One-dimensional position trace.
    track_length : float, optional
        Length of the linear track.

    Returns
    -------
    lap_numbers : ndarray
        Lap index for every sample.
    lap_starts : ndarray
        Sample indices at which each lap begins.
    """
    coords = np.mod(coords, track_length)

    wrap_threshold = -track_length / 2
    wrap_indices = np.where(np.diff(coords) < wrap_threshold)[0] + 1

    lap_starts = np.concatenate(([0], wrap_indices, [len(coords)]))

    lap_numbers = np.zeros(len(coords), dtype=int)
    lap_numbers[wrap_indices] = 1
    lap_numbers = np.cumsum(lap_numbers)

    return lap_numbers, lap_starts


def process_session(session_dir, track_length):
    """Generate lap labels for a single recording session."""

    coords = np.loadtxt(session_dir / "coords.csv")

    lap_numbers, lap_starts = compute_laps(coords, track_length)

    np.savetxt(session_dir / "laps.csv", lap_numbers, fmt="%d")
    np.savetxt(session_dir / "lap_starts.csv", lap_starts, fmt="%d")


def main():
    parser = argparse.ArgumentParser(
        description="Generate lap labels from behavioural position data."
    )

    parser.add_argument(
        "data_root",
        type=Path,
        help="Root directory containing experiment folders.",
    )

    parser.add_argument(
        "--track-length",
        type=float,
        default=280.0,
        help="Track length used for lap detection (default: 280).",
    )

    args = parser.parse_args()

    for experiment in EXPERIMENTS:
        for mouse in MICE:
            session_dir = args.data_root / experiment / mouse

            if not session_dir.exists():
                print(f"Skipping missing session: {session_dir}")
                continue

            print(f"Processing {experiment} / {mouse}")
            process_session(session_dir, args.track_length)


if __name__ == "__main__":
    main()