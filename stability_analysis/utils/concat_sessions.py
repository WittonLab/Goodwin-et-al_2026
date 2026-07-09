import argparse
from pathlib import Path

import pandas as pd


MICE = [
    "DG31",
    "DG32",
    "DG33",
    "DG35",
    "DG36",
]

EXPERIMENTS = ["250409", "250417", "250424", "250513"]


def merge_cell_traces(session_dir):
    """
    Concatenate pre- and post-switch calcium imaging traces.

    Parameters
    ----------
    session_dir : pathlib.Path
        Recording session directory containing pre/ and post/ folders.
    """

    pre_path = session_dir / "output_data_pre" / "Fc3_cleaned.csv"
    post_path = session_dir / "output_data_post" / "Fc3_cleaned.csv"

    if not pre_path.exists():
        print(f"Missing pre-switch data: {pre_path}")
        return

    if not post_path.exists():
        print(f"Missing post-switch data: {post_path}")
        return

    pre_cells = pd.read_csv(pre_path)
    post_cells = pd.read_csv(post_path)

    full_cells = pd.concat(
        [pre_cells, post_cells],
        axis=0,
        ignore_index=True,
    )

    full_cells.to_csv(
        session_dir / "Fc3_cleaned.csv",
        index=False,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Merge pre- and post-switch calcium imaging traces."
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

            print(f"Processing {experiment}/{mouse}")

            merge_cell_traces(session_dir)


if __name__ == "__main__":
    main()