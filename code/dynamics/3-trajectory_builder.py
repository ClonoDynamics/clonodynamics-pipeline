"""
3-trajectory_builder.py

Build clonotype trajectories in long format from a per-clone denoised table.

Purpose
-------
This script converts a denoised per-clone table (one or more rows per clonotype,
subject, and time point) into a standardized longitudinal trajectory table suitable
for downstream transition and dynamics analyses.

It performs only structural formatting and aggregation:
- selects one frequency column (via --freq-col)
- keeps the required identifiers: subject, time, aaSeqCDR3
- aggregates duplicate rows at the same subject/clone/time
- computes log-transformed frequency as:
      log_freq = ln(freq + epsilon)
- preserves the column `observable` if present
  (otherwise creates it and sets it to False)

This script intentionally applies:
- NO filtering of clonotypes
- NO model-selection criteria
- NO dropout correction
- NO transition building

Input
-----
CSV file such as:
  per_clone_denoised_*.csv

Required columns:
  - subject
  - time
  - aaSeqCDR3

Required frequency column:
  - one column specified with --freq-col
    (default: freq_geo_x)

Optional column:
  - observable

Behavior
--------
For each unique combination:
    (subject, aaSeqCDR3, time)

the script keeps one row only.
If duplicate rows are present, it aggregates them as:
  - freq       -> mean
  - observable -> max

Outputs
-------
Written to --outdir:
  - trajectories_long.csv
      Long-format trajectory table with:
      subject, aaSeqCDR3, time, freq, observable, log_freq

  - included_clones.csv
      Unique list of included (subject, aaSeqCDR3) pairs

  - trajectory_build_report.md
      Small report summarizing input, parameters, and output sizes

Typical use
-----------
1) Standard run with default frequency column:
   python 2-1_trajectory_builder.py \
       --per-clone results/per_clone_denoised.csv \
       --outdir results/trajectories_build

2) Use a different frequency column:
   python 2-1_trajectory_builder.py \
       --per-clone results/per_clone_denoised.csv \
       --outdir results/trajectories_build \
       --freq-col freq_geo_y

3) Use a different epsilon for log transform:
   python 2-1_trajectory_builder.py \
       --per-clone results/per_clone_denoised.csv \
       --outdir results/trajectories_build \
       --epsilon 1e-12

Notes
-----
- Frequencies must be numeric and >= 0.
- Rows with non-finite or negative frequency are discarded.
- If `observable` is present, it is converted robustly to boolean.
- This script prepares trajectories only; transition construction is performed later.

Examples
python ./code/dynamics/2-trajectory_builder.py \
    --per-clone ./results/2-noise_model/p_005/per_clone_denoised_subject.csv \
    --outdir ./results/3-trajectories/p_005
"""

from __future__ import annotations

import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd


REQUIRED = {"subject", "time", "aaSeqCDR3"}


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build trajectories_long.csv from per_clone_denoised_*.csv (no filtering).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--per-clone", type=Path, required=True, help="Path to per_clone_denoised_*.csv")
    p.add_argument("--outdir", type=Path, default=Path("results/trajectories_build"), help="Output directory")

    p.add_argument(
        "--freq-col",
        type=str,
        default="freq_geo_x",
        help="Frequency column to use as freq (e.g., freq_geo_x).",
    )
    p.add_argument(
        "--epsilon",
        type=float,
        default=1e-15,
        help="Small constant added before log transform: log_freq = ln(freq + epsilon).",
    )
    return p


def main() -> None:
    args = build_argparser().parse_args()
    outdir = ensure_dir(args.outdir)

    df = pd.read_csv(args.per_clone)

    missing = sorted(REQUIRED - set(df.columns))
    if missing:
        raise ValueError(f"Input per-clone file missing required columns: {missing}")

    if args.freq_col not in df.columns:
        raise ValueError(f"--freq-col '{args.freq_col}' not found. Available columns: {list(df.columns)}")

    # Coerce types
    df["subject"] = pd.to_numeric(df["subject"], errors="raise").astype(int)
    df["time"] = pd.to_numeric(df["time"], errors="raise")
    df["aaSeqCDR3"] = df["aaSeqCDR3"].astype(str)

    # Frequency
    freq = pd.to_numeric(df[args.freq_col], errors="coerce")
    df = df[np.isfinite(freq) & (freq >= 0)].copy()
    df["freq"] = freq.loc[df.index].to_numpy(dtype=float)

    # Keep observable if present, else create it (False)
    if "observable" in df.columns:
        # robust boolean coercion
        if df["observable"].dtype != bool:
            df["observable"] = df["observable"].astype(str).str.strip().str.lower().isin(["true", "t", "1", "yes", "y"])
        else:
            df["observable"] = df["observable"].astype(bool)
    else:
        df["observable"] = False

    # One row per subject/clone/time: if duplicates exist, aggregate
    agg = (
        df.groupby(["subject", "aaSeqCDR3", "time"], as_index=False)
        .agg(
            freq=("freq", "mean"),
            observable=("observable", "max"),
        )
    )
    agg["log_freq"] = np.log(agg["freq"].to_numpy(dtype=float) + float(args.epsilon))

    # Included clones = all unique (subject, clone) seen after aggregation
    included = (
        agg[["subject", "aaSeqCDR3"]]
        .drop_duplicates()
        .sort_values(["subject", "aaSeqCDR3"])
        .reset_index(drop=True)
    )

    traj_long = agg.sort_values(["subject", "aaSeqCDR3", "time"]).reset_index(drop=True)

    # Outputs
    traj_fp = outdir / "trajectories_long.csv"
    inc_fp = outdir / "included_clones.csv"
    report_fp = outdir / "trajectory_build_report.md"

    traj_long.to_csv(traj_fp, index=False)
    included.to_csv(inc_fp, index=False)

    overall = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input_per_clone": str(Path(args.per_clone)),
        "freq_col": str(args.freq_col),
        "epsilon": float(args.epsilon),
        "n_rows_in_after_freq_filter": int(len(df)),
        "n_rows_after_agg": int(len(agg)),
        "n_subjects": int(agg["subject"].nunique()) if len(agg) else 0,
        "n_unique_clones": int(included.shape[0]),
        "n_rows_trajectories_long": int(traj_long.shape[0]),
    }

    with open(report_fp, "w", encoding="utf-8") as f:
        f.write("# Trajectory build report\n\n")
        for k, v in overall.items():
            f.write(f"- {k}: {v}\n")
        f.write("\n## Outputs\n")
        f.write(f"- trajectories_long: `{traj_fp.name}`\n")
        f.write(f"- included_clones: `{inc_fp.name}`\n")

    print("Saved:")
    print(" -", traj_fp)
    print(" -", inc_fp)
    print(" -", report_fp)


if __name__ == "__main__":
    main()
