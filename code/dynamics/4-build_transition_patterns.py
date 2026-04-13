#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_transition_patterns.py

Construct transition patterns from long-format clonotype trajectories.

Default behavior:
- build ONLY the all-pairs dataset (`all`)
- build the adjacent/consecutive-observed dataset (`adj`) only if explicitly requested
  at launch with `--build_adj`.

Input long table columns (defaults):
  subject, aaSeqCDR3, time, freq, observable, log_freq

We define x(t) in log-frequency space:
  x = log_freq if provided and numeric, else x = ln(freq) with a floor.

Each row is a timepoint for a clonotype trajectory.

----------------------------------------------------------------------
PATTERNS
----------------------------------------------------------------------
(1) all_pairs ("all"):
    For each (subject, clone), include ALL endpoint pairs among observed timepoints:
      for any i<j in the sorted observed rows:
        (t_i -> t_j), dt = t_j - t_i
    Optional constraint: dt <= max_dt (default 6; can be disabled by max_dt<=0).

(2) consecutive_pairs ("adj"):
    For each (subject, clone), include transitions between consecutive observed rows
    after sorting by time:
      (t_i -> t_{i+1}), dt = t_{i+1} - t_i
    This is "adjacent in observed sequence", NOT necessarily dt=1.

Endpoint observation classes:
  obs_class ∈ {TT, TF, FT, FF} from observable at endpoints.

----------------------------------------------------------------------
WEIGHTS (optional)
----------------------------------------------------------------------
If --add_weights, we emit:
  w = w_clone * w_dt * w_class

  w_clone: 1 / n_transitions_for_that_(subject,clone) within the *pattern dataset*
  w_dt (optional via --balance_dt): 1 / n_transitions_for_that_dt
  w_class: TT=1, TF/FT=gamma_cross, FF=gamma_ff (default 0.25, 0.1)

Each component can be normalized to mean~1 via --normalize_weight_components.

----------------------------------------------------------------------
OUTPUTS
----------------------------------------------------------------------
Writes into --outdir:

  transitions_all.csv
  counts_by_subject_class_dt_all.csv
  transition_pattern_report.md

And, only if `--build_adj` is requested:

  transitions_adj.csv
  counts_by_subject_class_dt_adj.csv

Each transitions CSV includes:
  subject, aaSeqCDR3, t0, t1, dt, x0, x1, dx, obs0, obs1, obs_class,
  pattern, n_steps, w_clone, w_dt, w_class, w   (weights columns only if requested)

----------------------------------------------------------------------
EXAMPLE
----------------------------------------------------------------------
Default (build ONLY all-pairs transitions):
python build_transition_patterns.py \
  --input trajectories_long.csv \
  --outdir results/transitions_patterns \
  --max_dt 6

Build all-pairs + adjacent/consecutive-observed transitions:
python build_transition_patterns.py \
  --input trajectories_long.csv \
  --outdir results/transitions_patterns \
  --max_dt 6 \
  --build_adj \
  --add_weights --balance_dt \
  --gamma_cross 0.25 --gamma_ff 0.1 \
  --normalize_weight_components
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


# -----------------------------
# IO / utility
# -----------------------------
def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def parse_bool_like(series: pd.Series) -> pd.Series:
    """
    Convert a column containing boolean-like values to pandas boolean.
    Unrecognized values become NA.
    """
    if pd.api.types.is_bool_dtype(series):
        return series.astype("boolean")

    if pd.api.types.is_numeric_dtype(series):
        s = series.copy()
        out = pd.Series(pd.NA, index=s.index, dtype="boolean")
        out[s == 1] = True
        out[s == 0] = False
        return out

    s = series.astype(str).str.strip().str.lower()
    true_set = {"true", "t", "1", "yes", "y"}
    false_set = {"false", "f", "0", "no", "n"}

    out = pd.Series(pd.NA, index=series.index, dtype="boolean")
    out[s.isin(true_set)] = True
    out[s.isin(false_set)] = False
    return out


def safe_ln_freq(freq: np.ndarray, min_freq: float) -> np.ndarray:
    f = np.asarray(freq, float)
    f = np.where(np.isfinite(f), f, np.nan)
    f = np.where(f > 0, f, np.nan)
    if min_freq is not None and min_freq > 0:
        f = np.where(np.isnan(f), min_freq, np.maximum(f, min_freq))
    return np.log(f)


def obs_class_from_bool(o0: np.ndarray, o1: np.ndarray) -> np.ndarray:
    o0 = np.asarray(o0, bool)
    o1 = np.asarray(o1, bool)
    out = np.empty(o0.shape[0], dtype=object)
    out[(o0) & (o1)] = "TT"
    out[(o0) & (~o1)] = "TF"
    out[(~o0) & (o1)] = "FT"
    out[(~o0) & (~o1)] = "FF"
    return out


def summarize_counts(df: pd.DataFrame, subject_col: str, pattern_name: str, outdir: Path) -> None:
    out_path = outdir / f"counts_by_subject_class_dt_{pattern_name}.csv"
    if df.empty:
        pd.DataFrame(columns=[subject_col, "obs_class", "dt", "n"]).to_csv(out_path, index=False)
        return

    g = (
        df.groupby([subject_col, "obs_class", "dt"], observed=True)
          .size()
          .reset_index(name="n")
          .sort_values([subject_col, "obs_class", "dt"])
    )
    g.to_csv(out_path, index=False)


# -----------------------------
# Weights
# -----------------------------
def add_weights(
    df: pd.DataFrame,
    subject_col: str,
    clone_col: str,
    dt_col: str = "dt",
    class_col: str = "obs_class",
    balance_dt: bool = True,
    class_weights: Optional[Dict[str, float]] = None,
    normalize_components: bool = True,
) -> pd.DataFrame:
    """
    Add weight components and total weight:
      w = w_clone * w_dt * w_class
    """
    out = df.copy()
    if out.empty:
        out["w_clone"] = np.nan
        out["w_dt"] = np.nan
        out["w_class"] = np.nan
        out["w"] = np.nan
        return out

    # w_clone
    n_per_clone = out.groupby([subject_col, clone_col], observed=True).size().rename("n_clone")
    out = out.merge(n_per_clone.reset_index(), on=[subject_col, clone_col], how="left")
    out["w_clone"] = 1.0 / out["n_clone"].astype(float)

    # w_dt
    if balance_dt:
        n_per_dt = out.groupby(dt_col, observed=True).size().rename("n_dt")
        out = out.merge(n_per_dt.reset_index(), on=[dt_col], how="left")
        out["w_dt"] = 1.0 / out["n_dt"].astype(float)
    else:
        out["w_dt"] = 1.0

    # w_class
    if class_weights is None:
        out["w_class"] = 1.0
    else:
        out["w_class"] = out[class_col].map(class_weights).astype(float).fillna(1.0)

    if normalize_components:
        for c in ["w_clone", "w_dt", "w_class"]:
            m = np.isfinite(out[c]) & (out[c] > 0)
            if m.any():
                out[c] = out[c] / float(out.loc[m, c].mean())

    out["w"] = out["w_clone"] * out["w_dt"] * out["w_class"]

    for c in ["n_clone", "n_dt"]:
        if c in out.columns:
            out = out.drop(columns=[c])

    return out


# -----------------------------
# Builders (ONLY TWO)
# -----------------------------
def build_all_pairs(
    df: pd.DataFrame,
    subject_col: str,
    clone_col: str,
    time_col: str,
    max_dt: int,
) -> pd.DataFrame:
    """
    Pattern 1: ALL endpoint pairs among observed rows for each (subject, clone).
    Optionally restrict to dt<=max_dt if max_dt>0.
    """
    rows: List[dict] = []
    df = df.sort_values([subject_col, clone_col, time_col]).copy()

    for (subj, clone), g in df.groupby([subject_col, clone_col], sort=False, observed=True):
        g = g.sort_values(time_col)
        t = g[time_col].astype(int).to_numpy()
        x = g["x"].to_numpy(float)
        o = np.array([bool(v) for v in g["obsT"].to_numpy()])

        n = len(g)
        if n < 2:
            continue

        for i in range(n - 1):
            t0 = int(t[i])
            o0 = bool(o[i])
            for j in range(i + 1, n):
                t1 = int(t[j])
                dt = int(t1 - t0)
                if dt <= 0:
                    continue
                if max_dt and max_dt > 0 and dt > max_dt:
                    # because t is sorted increasing, dt will only grow with j
                    break

                o1 = bool(o[j])
                rows.append({
                    subject_col: subj,
                    clone_col: clone,
                    "t0": t0,
                    "t1": int(t1),
                    "dt": dt,
                    "x0": float(x[i]),
                    "x1": float(x[j]),
                    "dx": float(x[j] - x[i]),
                    "obs0": "T" if o0 else "F",
                    "obs1": "T" if o1 else "F",
                    "obs_class": obs_class_from_bool(np.array([o0]), np.array([o1]))[0],
                    "pattern": "all",
                    "n_steps": dt,
                })

    out = pd.DataFrame(rows)
    if out.empty:
        keep = [subject_col, clone_col, "t0", "t1", "dt", "x0", "x1", "dx",
                "obs0", "obs1", "obs_class", "pattern", "n_steps"]
        return pd.DataFrame(columns=keep)

    return out.sort_values([subject_col, clone_col, "t0", "dt", "t1"]).reset_index(drop=True)


def build_consecutive_observed(
    df: pd.DataFrame,
    subject_col: str,
    clone_col: str,
    time_col: str,
    max_dt: int,
) -> pd.DataFrame:
    """
    Pattern 2: transitions between consecutive observed rows for each (subject, clone),
    i.e. (row i -> row i+1) after sorting by time. Not necessarily dt=1.
    Optionally restrict to dt<=max_dt if max_dt>0.
    """
    df = df.sort_values([subject_col, clone_col, time_col]).copy()
    g = df.groupby([subject_col, clone_col], sort=False, observed=True)

    left = df.rename(columns={time_col: "t0", "x": "x0", "obsT": "obs0T"}).copy()
    right = g[[time_col, "x", "obsT"]].shift(-1).rename(
        columns={time_col: "t1", "x": "x1", "obsT": "obs1T"}
    )

    tmp = pd.concat([left.reset_index(drop=True), right.reset_index(drop=True)], axis=1)
    tmp = tmp.dropna(subset=["t1", "x1", "obs1T"]).copy()

    tmp["t0"] = tmp["t0"].astype(int)
    tmp["t1"] = tmp["t1"].astype(int)
    tmp["dt"] = (tmp["t1"] - tmp["t0"]).astype(int)
    tmp = tmp[tmp["dt"] > 0].copy()

    if max_dt and max_dt > 0:
        tmp = tmp[tmp["dt"] <= int(max_dt)].copy()

    tmp["dx"] = tmp["x1"] - tmp["x0"]

    tmp["obs0"] = np.where(tmp["obs0T"].astype(bool), "T", "F")
    tmp["obs1"] = np.where(tmp["obs1T"].astype(bool), "T", "F")
    tmp["obs_class"] = obs_class_from_bool(
        tmp["obs0T"].astype(bool).to_numpy(),
        tmp["obs1T"].astype(bool).to_numpy(),
    )

    tmp["pattern"] = "adj"
    tmp["n_steps"] = tmp["dt"].astype(int)

    keep = [subject_col, clone_col, "t0", "t1", "dt", "x0", "x1", "dx",
            "obs0", "obs1", "obs_class", "pattern", "n_steps"]
    return tmp[keep].reset_index(drop=True)


# -----------------------------
# CLI
# -----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build transition patterns from trajectories. By default builds only all endpoint pairs (all); build adjacent/consecutive-observed pairs (adj) only with --build_adj."
    )
    p.add_argument("--input", required=True, help="Input trajectories long CSV")
    p.add_argument("--outdir", required=True, help="Output directory")

    # columns
    p.add_argument("--subject_col", default="subject")
    p.add_argument("--clone_col", default="aaSeqCDR3")
    p.add_argument("--time_col", default="time")

    p.add_argument("--freq_col", default="freq",
                   help="Frequency column used if --log_freq_col not present/usable.")
    p.add_argument("--log_freq_col", default="log_freq",
                   help="If present, use as x=log_freq (assumed natural log by default).")
    p.add_argument("--observable_col", default="observable",
                   help="Boolean-like observable column (True/False, 1/0, T/F, etc.).")

    # parameters
    p.add_argument("--max_dt", type=int, default=6,
                   help="Maximum dt to generate. Use 0 or negative to disable dt filtering (keep all dt).")
    p.add_argument("--min_freq", type=float, default=1e-12,
                   help="Floor before log if computing from freq.")
    p.add_argument("--build_adj", action="store_true",
                   help="Also build adjacent/consecutive-observed transitions. By default only the all-pairs dataset is produced.")

    # weights
    p.add_argument("--add_weights", action="store_true",
                   help="Add weights w_clone, w_dt, w_class and total w to each output.")
    p.add_argument("--balance_dt", action="store_true",
                   help="If adding weights, include w_dt to balance dt groups (recommended for all pattern).")
    p.add_argument("--gamma_cross", type=float, default=0.25,
                   help="Class weight for TF/FT when adding weights (diagnostic-only).")
    p.add_argument("--gamma_ff", type=float, default=0.1,
                   help="Class weight for FF when adding weights (diagnostic-only).")
    p.add_argument("--normalize_weight_components", action="store_true",
                   help="Normalize each weight component to mean~1 (recommended).")

    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = ensure_dir(Path(args.outdir))

    # Load
    df = pd.read_csv(args.input)
    required = {args.subject_col, args.clone_col, args.time_col, args.observable_col}
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns in input: {missing}")

    df = df.copy()
    df[args.time_col] = df[args.time_col].astype(int)

    # Parse observable
    obs_parsed = parse_bool_like(df[args.observable_col])
    n_na_obs = int(obs_parsed.isna().sum())
    if n_na_obs > 0:
        obs_parsed = obs_parsed.fillna(False)  # conservative
    df["obsT"] = obs_parsed.astype(bool)

    # Choose x (log-frequency)
    if args.log_freq_col in df.columns:
        x = pd.to_numeric(df[args.log_freq_col], errors="coerce")
        if x.isna().any() and args.freq_col in df.columns:
            x_fallback = safe_ln_freq(df[args.freq_col].to_numpy(), args.min_freq)
            x = x.fillna(pd.Series(x_fallback, index=df.index))
        df["x"] = x
        x_source = args.log_freq_col
    else:
        if args.freq_col not in df.columns:
            raise SystemExit(f"Neither '{args.log_freq_col}' present nor '{args.freq_col}' available.")
        df["x"] = safe_ln_freq(df[args.freq_col].to_numpy(), args.min_freq)
        x_source = f"ln({args.freq_col})"

    # Drop non-finite x
    n_in = len(df)
    bad_x = ~np.isfinite(df["x"].to_numpy(float))
    n_bad_x = int(bad_x.sum())
    if n_bad_x > 0:
        df = df.loc[~bad_x].copy()

    # Class weights
    class_w = {
        "TT": 1.0,
        "TF": float(args.gamma_cross),
        "FT": float(args.gamma_cross),
        "FF": float(args.gamma_ff),
    }

    max_dt = int(args.max_dt)

    print(f"[INFO] Input rows: {n_in} | After x-QC: {len(df)} | Dropped non-finite x: {n_bad_x}")
    print(f"[INFO] Observable parse: NA->False count = {n_na_obs}")
    print(f"[INFO] x source: {x_source}")
    print(f"[INFO] max_dt filter: {max_dt} (<=0 means no filtering)")

    # Pattern 1: all
    all_df = build_all_pairs(df, args.subject_col, args.clone_col, args.time_col, max_dt)
    if args.add_weights:
        all_df = add_weights(
            all_df,
            subject_col=args.subject_col,
            clone_col=args.clone_col,
            balance_dt=bool(args.balance_dt),
            class_weights=class_w,
            normalize_components=bool(args.normalize_weight_components),
        )
    all_df.to_csv(outdir / "transitions_all.csv", index=False)
    summarize_counts(all_df, args.subject_col, "all", outdir)

    # Pattern 2: adj (consecutive observed) -- optional
    adj_df = None
    if args.build_adj:
        adj_df = build_consecutive_observed(df, args.subject_col, args.clone_col, args.time_col, max_dt)
        if args.add_weights:
            # adj is still multilag; w_dt is optional; by default you probably want it OFF here
            adj_df = add_weights(
                adj_df,
                subject_col=args.subject_col,
                clone_col=args.clone_col,
                balance_dt=False,
                class_weights=class_w,
                normalize_components=bool(args.normalize_weight_components),
            )
        adj_df.to_csv(outdir / "transitions_adj.csv", index=False)
        summarize_counts(adj_df, args.subject_col, "adj", outdir)

    # Markdown report
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report: List[str] = []
    report.append("# Transition pattern report (2 patterns)\n")
    report.append(f"Generated: {now}\n")
    report.append("## Input\n")
    report.append(f"- input: `{args.input}`\n")
    report.append(f"- rows (original): {n_in}\n")
    report.append(f"- rows (after x QC): {len(df)} (dropped non-finite x: {n_bad_x})\n")
    report.append(f"- subject_col: `{args.subject_col}`\n")
    report.append(f"- clone_col: `{args.clone_col}`\n")
    report.append(f"- time_col: `{args.time_col}`\n")
    report.append(f"- observable_col: `{args.observable_col}` (parsed; NA->False: {n_na_obs})\n")
    report.append(f"- x source: `{x_source}`\n")
    report.append(f"- max_dt: {max_dt} (<=0 means no filtering)\n")

    report.append("\n## Patterns generated\n")
    report.append(f"- all: {len(all_df)} rows (all endpoint pairs among observed times)\n")
    if args.build_adj:
        report.append(f"- adj: {len(adj_df)} rows (consecutive observed rows per trajectory)\n")
    else:
        report.append("- adj: not generated (use `--build_adj` to request it)\n")

    report.append("\n## Weights\n")
    report.append(f"- add_weights: {bool(args.add_weights)}\n")
    if args.add_weights:
        report.append("  - w = w_clone * w_dt * w_class\n")
        report.append(f"  - balance_dt applied to ALL pattern: {bool(args.balance_dt)}\n")
        report.append(f"  - class weights: TT=1, TF/FT={args.gamma_cross}, FF={args.gamma_ff}\n")
        report.append(f"  - normalize_weight_components: {bool(args.normalize_weight_components)}\n")
        report.append("  - Note: if ADJ is requested, w_dt is disabled by default for that pattern (balance_dt=False).\n")

    report.append("\n## Count tables\n")
    report.append("- counts_by_subject_class_dt_all.csv\n")
    if args.build_adj:
        report.append("- counts_by_subject_class_dt_adj.csv\n")

    report.append("\n## Notes\n")
    if args.build_adj:
        report.append("- `adj` means consecutive *observed* timepoints, not necessarily dt=1.\n")
        report.append("- `all` includes adj as a subset (plus non-adj pairs).\n")
    else:
        report.append("- Default mode builds only the `all` dataset. Use `--build_adj` to also write adjacent/consecutive-observed transitions.\n")

    (outdir / "transition_pattern_report.md").write_text("\n".join(report), encoding="utf-8")

    print(f"[OK] Wrote transition datasets to: {outdir.resolve()}")
    print("     - transitions_all.csv")
    if args.build_adj:
        print("     - transitions_adj.csv")


if __name__ == "__main__":
    main()