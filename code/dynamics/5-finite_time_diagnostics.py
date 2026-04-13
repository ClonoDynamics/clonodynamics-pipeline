#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
finite_time_diagnostics.py
=========================

Unified analysis script for:
  1) local model-agnostic diagnostics
  2) finite-time drift/diffusion summaries

All outputs are standardized to NATURAL LOG space (ln).

Modes
-----
--mode local
    Compute local diagnostics:
      - drift_by_dt_xbins.csv
      - msd_by_dt.csv
      - msd_by_dt_xbins.csv
      - msd_by_dt_freqbins.csv
      - dt_freqbin_counts.csv

--mode finite_time
    Compute finite-time summaries:
      - finite_time_summary.csv
      - finite_time_by_dt_xbins_boot.csv

--mode all
    Run both blocks.

Input
-----
Required columns:
  subject, aaSeqCDR3, dt, x0, x1, dx, obs_class

Optional:
  w

Log base
--------
Input can be in ln or log10:
  --x_base ln
  --x_base log10

Internally and in all outputs:
  ln_f0, ln_f1, ln_fgeo, dx are in ln space.

  Examples
For weighted transitions:
  python ./code/dynamics/5-finite_time_diagnostics.py \
    --mode all \
    --x_base ln \
    --n_boot 300 \
    --bootstrap_reps_ft 300 \
    --min_rows_per_bin 200 \
    --nbins_x 6 \
    --weight_col w \
    --dt_max 5 \
    --use_wls_slope \
--input ./results/4-transitions/p_02/transitions_all.csv \
--out_dir ./results/5-finite_time_diagnostics/TT_w   

For unweighted tranistions:
    python ./code/dynamics/finite_time_diagnostics.py \
    --mode all \
    --x_base ln \
    --n_boot 300 \
    --bootstrap_reps_ft 300 \
    --min_rows_per_bin 200 \
    --nbins_x 6 \
    --dt_max 5 \
    --use_wls_slope \
    --input ./results/4-transitions/p_02/transitions_all.csv \
    --out_dir ./results/5-finite_time_diagnostics/TT 
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

LN10 = float(np.log(10.0))


# =========================================================
# CLI
# =========================================================
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Unified local + finite-time model-agnostic analysis in ln space."
    )

    ap.add_argument("--input", required=True, help="Transitions CSV")
    ap.add_argument("--out_dir", required=True, help="Output directory")

    ap.add_argument(
        "--mode",
        choices=["local", "finite_time", "all"],
        default="all",
        help="Which analysis block to run."
    )

    ap.add_argument(
        "--x_base",
        choices=["ln", "log10"],
        default="ln",
        help="Input base of x0/x1/dx. Outputs are always in ln."
    )

    ap.add_argument("--csv_sep", default=",", help="CSV separator")
    ap.add_argument("--min_dt", type=float, default=0.0, help="Keep dt > min_dt")
    ap.add_argument("--dt_max", type=int, default=5, help="Keep dt <= dt_max")
    ap.add_argument("--dt_values", default="", help="Optional comma-separated dt values to keep exactly")

    ap.add_argument("--only_TT", action="store_true", help="Keep only obs_class == tt_value")
    ap.add_argument("--tt_value", default="TT", help="Value of obs_class representing TT")

    ap.add_argument(
        "--classes",
        default="",
        help="Optional comma-separated obs_class values to keep. Overrides --only_TT if given."
    )

    ap.add_argument(
        "--weight_col",
        default="w",
        help="Weight column name, or '1' for unweighted."
    )

    ap.add_argument(
        "--subsample_per_dt",
        type=int,
        default=0,
        help="If >0, subsample N rows per dt before analysis."
    )

    # Binning
    ap.add_argument("--nbins_x", type=int, default=30, help="Bins along ln_f0")
    ap.add_argument("--nbins_f", type=int, default=6, help="Bins along ln_fgeo")
    ap.add_argument("--min_rows_per_bin", type=int, default=200, help="Minimum rows per grouped estimate")

    ap.add_argument(
        "--binning_mode",
        choices=["per_dataset_quantile", "fixed_edges", "save_edges_only"],
        default="per_dataset_quantile",
        help="Binning strategy."
    )
    ap.add_argument("--x_edges_csv", default="", help="Fixed x edges CSV")
    ap.add_argument("--f_edges_csv", default="", help="Fixed f edges CSV")
    ap.add_argument("--save_x_edges_name", default="xbin_edges.csv", help="Saved x edges filename")
    ap.add_argument("--save_f_edges_name", default="fbin_edges.csv", help="Saved f edges filename")

    # Bootstrap
    ap.add_argument("--n_boot", type=int, default=300, help="Bootstrap reps for local diagnostics")
    ap.add_argument("--seed", type=int, default=123, help="Random seed")

    # Finite-time specific
    ap.add_argument("--window_bins", type=int, default=5, help="Bins used for local slope fit")
    ap.add_argument(
        "--crossing_mode",
        choices=["first", "closest_to_zero"],
        default="first",
        help="How to choose zero crossing in median(dx|x0)"
    )
    ap.add_argument("--mad_to_sigma", type=float, default=1.4826, help="MAD -> sigma factor")
    ap.add_argument(
        "--robust_center_mode",
        choices=["median", "weighted_median_if_possible"],
        default="weighted_median_if_possible",
        help="Center used for MAD"
    )
    ap.add_argument("--bootstrap_reps_ft", type=int, default=300, help="Bootstrap reps for finite-time block")
    ap.add_argument("--bootstrap_ci", type=float, default=0.95, help="Finite-time CI level")
    ap.add_argument("--cluster_cols", default="subject,aaSeqCDR3", help="Bootstrap cluster columns")
    ap.add_argument("--use_wls_slope", action="store_true", help="Use weighted local slope fit if possible")
    ap.add_argument("--include_total_dt", action="store_true", help="Include TOTAL rows in dt_freqbin_counts.csv")

    ap.add_argument("--debug", action="store_true", help="Print debug info")
    return ap


# =========================================================
# General helpers
# =========================================================
def parse_list_str(s: str) -> List[str]:
    return [x.strip() for x in (s or "").split(",") if x.strip()]


def parse_list_int(s: str) -> List[int]:
    return [int(x.strip()) for x in (s or "").split(",") if x.strip()]


def finite_mask(*arrs: np.ndarray) -> np.ndarray:
    m = np.ones(len(arrs[0]), dtype=bool)
    for a in arrs:
        m &= np.isfinite(a)
    return m


def read_and_prepare(path: Path, x_base: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep)
    required = {"subject", "aaSeqCDR3", "dt", "x0", "x1", "dx", "obs_class"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    d = df.copy()
    for c in ["dt", "x0", "x1", "dx"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    d = d.dropna(subset=["dt", "x0", "x1", "dx"]).copy()
    d = d[np.isfinite(d["dt"]) & np.isfinite(d["x0"]) & np.isfinite(d["x1"]) & np.isfinite(d["dx"])].copy()

    if x_base == "ln":
        d["ln_f0"] = d["x0"]
        d["ln_f1"] = d["x1"]
    else:
        d["ln_f0"] = d["x0"] * LN10
        d["ln_f1"] = d["x1"] * LN10
        d["dx"] = d["dx"] * LN10

    d["ln_fgeo"] = 0.5 * (d["ln_f0"] + d["ln_f1"])
    d["dt"] = d["dt"].astype(int)
    return d


def filter_dataframe(
    d: pd.DataFrame,
    min_dt: float,
    dt_max: int,
    dt_values: List[int],
    only_TT: bool,
    tt_value: str,
    classes: List[str],
) -> pd.DataFrame:
    out = d.copy()

    if classes:
        out = out[out["obs_class"].astype(str).isin(classes)].copy()
    elif only_TT:
        out = out[out["obs_class"].astype(str) == str(tt_value)].copy()

    out = out[out["dt"] > float(min_dt)].copy()
    out = out[out["dt"] <= int(dt_max)].copy()

    if dt_values:
        out = out[out["dt"].isin(dt_values)].copy()

    return out.reset_index(drop=True)


def subsample_by_dt(d: pd.DataFrame, n_per_dt: int, seed: int) -> pd.DataFrame:
    if n_per_dt <= 0:
        return d
    rng = np.random.default_rng(seed)
    parts = []
    for dt, g in d.groupby("dt", sort=True):
        if len(g) <= n_per_dt:
            parts.append(g)
        else:
            idx = rng.choice(g.index.to_numpy(), size=int(n_per_dt), replace=False)
            parts.append(g.loc[idx])
    return pd.concat(parts, axis=0).reset_index(drop=True)


def get_weights(df: pd.DataFrame, weight_col: str) -> np.ndarray:
    if weight_col.strip() == "1":
        return np.ones(len(df), dtype=float)
    if weight_col not in df.columns:
        raise ValueError(f"weight_col '{weight_col}' not found.")
    w = pd.to_numeric(df[weight_col], errors="coerce").to_numpy(float)
    return np.where(np.isfinite(w) & (w > 0), w, 0.0)


# =========================================================
# Binning
# =========================================================
def compute_quantile_edges(values: np.ndarray, nbins: int) -> np.ndarray:
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise RuntimeError("No finite values to compute edges.")
    edges = np.unique(np.quantile(x, np.linspace(0, 1, nbins + 1)))
    if len(edges) < 3:
        raise RuntimeError("Not enough distinct values to form bins.")
    return edges.astype(float)


def save_edges_csv(edges: np.ndarray, out_csv: Path, value_col: str) -> None:
    pd.DataFrame({value_col: np.asarray(edges, dtype=float)}).to_csv(out_csv, index=False)


def load_edges_csv(path: Path) -> np.ndarray:
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"Edges CSV is empty: {path}")
    col = df.columns[0]
    edges = pd.to_numeric(df[col], errors="coerce").dropna().to_numpy(float)
    edges = np.unique(edges)
    if len(edges) < 3:
        raise ValueError(f"Need at least 3 unique edges in {path}")
    if not np.all(np.diff(edges) > 0):
        raise ValueError(f"Edges must be strictly increasing in {path}")
    return edges


def assign_bins_from_edges(d: pd.DataFrame, col: str, edges: np.ndarray, prefix: str) -> pd.DataFrame:
    bins = pd.cut(d[col], bins=edges, include_lowest=True, duplicates="drop")
    out = d.copy()
    out[f"{prefix}_bin"] = bins
    cats = out[f"{prefix}_bin"].cat.categories

    meta = pd.DataFrame({
        f"{prefix}_bin": cats,
        f"{prefix}_center": [(c.left + c.right) / 2.0 for c in cats],
        f"{prefix}_lo": [c.left for c in cats],
        f"{prefix}_hi": [c.right for c in cats],
        f"{prefix}_id": np.arange(1, len(cats) + 1),
    })
    out = out.merge(meta, on=f"{prefix}_bin", how="left")
    return out


def define_edges(
    d: pd.DataFrame,
    binning_mode: str,
    nbins_x: int,
    nbins_f: int,
    x_edges_csv: str,
    f_edges_csv: str,
) -> Tuple[np.ndarray, np.ndarray]:
    if binning_mode == "per_dataset_quantile":
        x_edges = compute_quantile_edges(d["ln_f0"].to_numpy(float), nbins_x)
        f_edges = compute_quantile_edges(d["ln_fgeo"].to_numpy(float), nbins_f)
        return x_edges, f_edges

    if binning_mode == "fixed_edges":
        if not x_edges_csv or not f_edges_csv:
            raise ValueError("fixed_edges mode requires both --x_edges_csv and --f_edges_csv")
        return load_edges_csv(Path(x_edges_csv)), load_edges_csv(Path(f_edges_csv))

    raise ValueError(f"Unsupported edge definition mode: {binning_mode}")


# =========================================================
# Weighted statistics
# =========================================================
def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(m):
        return np.nan
    x, w = x[m], w[m]
    return float(np.sum(w * x) / np.sum(w))


def weighted_prob_positive(x: np.ndarray, w: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(m):
        return np.nan
    x, w = x[m], w[m]
    return float(np.sum(w * (x > 0)) / np.sum(w))


def weighted_quantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(m):
        return np.nan
    x, w = x[m], w[m]
    order = np.argsort(x)
    x, w = x[order], w[order]
    cw = np.cumsum(w)
    cutoff = q * cw[-1]
    idx = np.searchsorted(cw, cutoff, side="left")
    return float(x[min(idx, len(x) - 1)])


def weighted_median(x: np.ndarray, w: np.ndarray) -> float:
    return weighted_quantile(x, w, 0.5)


def bootstrap_group(
    x: np.ndarray,
    w: np.ndarray,
    stat: str,
    n_boot: int,
    rng: np.random.Generator,
) -> Tuple[float, float, float]:
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[m], w[m]
    n = x.size
    if n == 0:
        return np.nan, np.nan, np.nan

    if stat == "wmedian":
        fn = weighted_median
    elif stat == "ppos":
        fn = weighted_prob_positive
    elif stat == "wmean":
        fn = weighted_mean
    else:
        raise ValueError(f"Unknown stat: {stat}")

    point = fn(x, w)
    if n < 50 or n_boot <= 0:
        return float(point), np.nan, np.nan

    boots = np.empty(int(n_boot), dtype=float)
    for i in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        boots[i] = fn(x[idx], w[idx])

    lo, hi = np.quantile(boots, [0.025, 0.975])
    return float(point), float(lo), float(hi)


# =========================================================
# Local block
# =========================================================
def compute_dt_freqbin_counts(
    d: pd.DataFrame,
    include_total_dt: bool,
    x_edges: np.ndarray,
    f_edges: np.ndarray,
) -> pd.DataFrame:
    parts = []

    if include_total_dt:
        dt_tot = d.groupby("dt", sort=True).size().rename("n").reset_index()
        dt_tot["bin_type"] = "TOTAL"
        dt_tot["bin_id"] = np.nan
        dt_tot["bin_lo"] = np.nan
        dt_tot["bin_hi"] = np.nan
        dt_tot["bin_center"] = np.nan
        parts.append(dt_tot[["bin_type", "dt", "bin_id", "bin_lo", "bin_hi", "bin_center", "n"]])

    dx = assign_bins_from_edges(d, col="ln_f0", edges=x_edges, prefix="x").dropna(subset=["x_id"]).copy()
    if not dx.empty:
        dx["x_id"] = dx["x_id"].astype(int)
        x_counts = (
            dx.groupby(["dt", "x_id", "x_lo", "x_hi", "x_center"], sort=True)
            .size().rename("n").reset_index()
        )
        x_counts["bin_type"] = "ln_f0"
        x_counts = x_counts.rename(columns={
            "x_id": "bin_id", "x_lo": "bin_lo", "x_hi": "bin_hi", "x_center": "bin_center"
        })
        parts.append(x_counts[["bin_type", "dt", "bin_id", "bin_lo", "bin_hi", "bin_center", "n"]])

    df = assign_bins_from_edges(d, col="ln_fgeo", edges=f_edges, prefix="f").dropna(subset=["f_id"]).copy()
    if not df.empty:
        df["f_id"] = df["f_id"].astype(int)
        f_counts = (
            df.groupby(["dt", "f_id", "f_lo", "f_hi", "f_center"], sort=True)
            .size().rename("n").reset_index()
        )
        f_counts["bin_type"] = "ln_fgeo"
        f_counts = f_counts.rename(columns={
            "f_id": "bin_id", "f_lo": "bin_lo", "f_hi": "bin_hi", "f_center": "bin_center"
        })
        parts.append(f_counts[["bin_type", "dt", "bin_id", "bin_lo", "bin_hi", "bin_center", "n"]])

    if not parts:
        return pd.DataFrame(columns=["bin_type", "dt", "bin_id", "bin_lo", "bin_hi", "bin_center", "n"])

    out = pd.concat(parts, axis=0, ignore_index=True)
    return out.sort_values(["bin_type", "dt", "bin_id"], kind="mergesort").reset_index(drop=True)


def compute_drift_tables(
    d: pd.DataFrame,
    weight_col: str,
    x_edges: np.ndarray,
    min_rows_per_bin: int,
    n_boot: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    d2 = assign_bins_from_edges(d, col="ln_f0", edges=x_edges, prefix="x").dropna(subset=["x_id"]).copy()
    if d2.empty:
        return pd.DataFrame(columns=[
            "dt", "x_id", "x_center", "x_lo", "x_hi", "n",
            "median_dx", "median_lo", "median_hi",
            "ppos", "ppos_lo", "ppos_hi"
        ])
    d2["x_id"] = d2["x_id"].astype(int)

    rows = []
    for (dt, x_id), g in d2.groupby(["dt", "x_id"], sort=True):
        n = int(len(g))
        x_center = float(g["x_center"].iloc[0])
        x_lo = float(g["x_lo"].iloc[0])
        x_hi = float(g["x_hi"].iloc[0])

        if n < min_rows_per_bin:
            rows.append({
                "dt": int(dt), "x_id": int(x_id),
                "x_center": x_center, "x_lo": x_lo, "x_hi": x_hi, "n": n,
                "median_dx": np.nan, "median_lo": np.nan, "median_hi": np.nan,
                "ppos": np.nan, "ppos_lo": np.nan, "ppos_hi": np.nan,
            })
            continue

        dx = g["dx"].to_numpy(float)
        wg = get_weights(g, weight_col)

        med, med_lo, med_hi = bootstrap_group(dx, wg, "wmedian", n_boot, rng)
        ppos, ppos_lo, ppos_hi = bootstrap_group(dx, wg, "ppos", n_boot, rng)

        rows.append({
            "dt": int(dt), "x_id": int(x_id),
            "x_center": x_center, "x_lo": x_lo, "x_hi": x_hi, "n": n,
            "median_dx": med, "median_lo": med_lo, "median_hi": med_hi,
            "ppos": ppos, "ppos_lo": ppos_lo, "ppos_hi": ppos_hi,
        })

    return pd.DataFrame(rows).sort_values(["dt", "x_center"]).reset_index(drop=True)


def compute_msd_by_dt(d: pd.DataFrame, weight_col: str, n_boot: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for dt, g in d.groupby("dt", sort=True):
        dx2 = g["dx"].to_numpy(float) ** 2
        wg = get_weights(g, weight_col)
        msd, lo, hi = bootstrap_group(dx2, wg, "wmean", n_boot, rng)
        rows.append({"dt": int(dt), "n": int(len(g)), "msd": msd, "msd_lo": lo, "msd_hi": hi})
    return pd.DataFrame(rows).sort_values("dt").reset_index(drop=True)


def compute_msd_by_dt_xbins(
    d: pd.DataFrame,
    weight_col: str,
    x_edges: np.ndarray,
    min_rows_per_bin: int,
    n_boot: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    d2 = assign_bins_from_edges(d, col="ln_f0", edges=x_edges, prefix="x").dropna(subset=["x_id"]).copy()
    if d2.empty:
        return pd.DataFrame(columns=["dt", "x_id", "x_center", "x_lo", "x_hi", "n", "msd", "msd_lo", "msd_hi"])
    d2["x_id"] = d2["x_id"].astype(int)

    rows = []
    for (dt, x_id), g in d2.groupby(["dt", "x_id"], sort=True):
        n = int(len(g))
        x_center = float(g["x_center"].iloc[0])
        x_lo = float(g["x_lo"].iloc[0])
        x_hi = float(g["x_hi"].iloc[0])

        if n < min_rows_per_bin:
            rows.append({
                "dt": int(dt), "x_id": int(x_id),
                "x_center": x_center, "x_lo": x_lo, "x_hi": x_hi, "n": n,
                "msd": np.nan, "msd_lo": np.nan, "msd_hi": np.nan,
            })
            continue

        dx2 = g["dx"].to_numpy(float) ** 2
        wg = get_weights(g, weight_col)
        msd, lo, hi = bootstrap_group(dx2, wg, "wmean", n_boot, rng)

        rows.append({
            "dt": int(dt), "x_id": int(x_id),
            "x_center": x_center, "x_lo": x_lo, "x_hi": x_hi, "n": n,
            "msd": msd, "msd_lo": lo, "msd_hi": hi,
        })

    return pd.DataFrame(rows).sort_values(["dt", "x_center"]).reset_index(drop=True)


def compute_msd_by_dt_freqbins(
    d: pd.DataFrame,
    weight_col: str,
    f_edges: np.ndarray,
    min_rows_per_bin: int,
    n_boot: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    d2 = assign_bins_from_edges(d, col="ln_fgeo", edges=f_edges, prefix="f").dropna(subset=["f_id"]).copy()
    if d2.empty:
        return pd.DataFrame(columns=["f_id", "dt", "f_center", "f_lo", "f_hi", "n", "msd", "msd_lo", "msd_hi"])
    d2["f_id"] = d2["f_id"].astype(int)

    rows = []
    for (f_id, dt), g in d2.groupby(["f_id", "dt"], sort=True):
        n = int(len(g))
        f_center = float(g["f_center"].iloc[0])
        f_lo = float(g["f_lo"].iloc[0])
        f_hi = float(g["f_hi"].iloc[0])

        if n < min_rows_per_bin:
            rows.append({
                "f_id": int(f_id), "dt": int(dt),
                "f_center": f_center, "f_lo": f_lo, "f_hi": f_hi, "n": n,
                "msd": np.nan, "msd_lo": np.nan, "msd_hi": np.nan,
            })
            continue

        dx2 = g["dx"].to_numpy(float) ** 2
        wg = get_weights(g, weight_col)
        msd, lo, hi = bootstrap_group(dx2, wg, "wmean", n_boot, rng)

        rows.append({
            "f_id": int(f_id), "dt": int(dt),
            "f_center": f_center, "f_lo": f_lo, "f_hi": f_hi, "n": n,
            "msd": msd, "msd_lo": lo, "msd_hi": hi,
        })

    return pd.DataFrame(rows).sort_values(["f_center", "dt"]).reset_index(drop=True)


# =========================================================
# Finite-time block
# =========================================================
def usable_weights(w: Optional[np.ndarray]) -> bool:
    if w is None:
        return False
    w = np.asarray(w, float)
    return bool(np.any(np.isfinite(w) & (w > 0)))


def weighted_mad(values: np.ndarray, weights: Optional[np.ndarray], center: Optional[float] = None) -> float:
    v = np.asarray(values, float)
    if weights is None:
        v = v[np.isfinite(v)]
        if v.size == 0:
            return np.nan
        c = float(np.median(v)) if center is None else float(center)
        return float(np.median(np.abs(v - c)))

    w = np.asarray(weights, float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[m], w[m]
    if v.size == 0:
        return np.nan
    c = weighted_median(v, w) if center is None else float(center)
    return float(weighted_median(np.abs(v - c), w))


def mad_to_sigma(mad: float, c: float) -> float:
    if not np.isfinite(mad):
        return np.nan
    return float(c * mad)


def find_crossings(x: np.ndarray, y: np.ndarray) -> List[int]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 2:
        return []
    s = np.sign(y)
    s[s == 0] = 1.0
    out = []
    for i in range(y.size - 1):
        if s[i] * s[i + 1] < 0:
            out.append(i)
    return out


def choose_crossing_idx(x: np.ndarray, y: np.ndarray, mode: str) -> Optional[int]:
    idxs = find_crossings(x, y)
    if not idxs:
        return None
    if mode == "first":
        return idxs[0]
    best_i = idxs[0]
    best = np.inf
    for i in idxs:
        score = min(abs(y[i]), abs(y[i + 1]))
        if score < best:
            best = score
            best_i = i
    return best_i


def interpolate_zero(x0: float, y0: float, x1: float, y1: float) -> float:
    den = y1 - y0
    if not np.isfinite(den) or abs(den) < 1e-15:
        return np.nan
    t = -y0 / den
    return float(x0 + t * (x1 - x0))


def pick_window(n: int, center: int, w: int) -> np.ndarray:
    w = max(2, int(w))
    half = w // 2
    lo = max(0, center - half)
    hi = min(n, center + half + 1)
    idx = np.arange(lo, hi)
    while idx.size < w and (lo > 0 or hi < n):
        if lo > 0:
            lo -= 1
        elif hi < n:
            hi += 1
        idx = np.arange(lo, hi)
    return idx


def slope_fit(x: np.ndarray, y: np.ndarray, w: Optional[np.ndarray], use_wls: bool) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if w is not None:
        w = np.asarray(w, float)
        m &= np.isfinite(w) & (w > 0)

    x, y = x[m], y[m]
    w = w[m] if w is not None else None

    if x.size < 2:
        return np.nan

    X = np.column_stack([np.ones(x.size), x])
    if w is not None and use_wls:
        sw = np.sqrt(w)
        beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    else:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(beta[1])


def tau_drift(dt: float, slope: float) -> float:
    if not np.isfinite(dt) or not np.isfinite(slope):
        return np.nan
    if slope >= 0:
        return np.nan
    if abs(slope) < 1e-12:
        return np.inf
    return float(dt / (-slope))


def build_moments_curve(
    df: pd.DataFrame,
    dt_val: int,
    weight_col: str,
    x_edges: np.ndarray,
    mad_sigma_c: float,
    robust_center_mode: str,
    min_rows_per_bin: int,
) -> pd.DataFrame:
    d2 = assign_bins_from_edges(df, col="ln_f0", edges=x_edges, prefix="x").dropna(subset=["x_id"]).copy()
    if d2.empty:
        return pd.DataFrame(columns=[
            "bin", "x_center", "n", "w_sum",
            "median_dx", "mad_dx", "sigma_dx", "var_mad_dx", "D_mad"
        ])
    d2["x_id"] = d2["x_id"].astype(int)

    rows = []
    for x_id, g in d2.groupby("x_id", sort=True):
        n = int(len(g))
        if n < min_rows_per_bin:
            continue

        xk = g["ln_f0"].to_numpy(float)
        dxk = g["dx"].to_numpy(float)
        wk = get_weights(g, weight_col)
        wk_use = wk if usable_weights(wk) else None

        med_dx = weighted_median(dxk, wk_use) if wk_use is not None else float(np.median(dxk))

        if robust_center_mode == "median" or wk_use is None:
            center = float(np.median(dxk))
            mad = float(np.median(np.abs(dxk - center)))
        else:
            center = float(weighted_median(dxk, wk_use))
            mad = float(weighted_mad(dxk, wk_use, center=center))

        sigma = mad_to_sigma(mad, mad_sigma_c)
        var_mad = float(sigma * sigma) if np.isfinite(sigma) else np.nan
        D = float(var_mad / (2.0 * dt_val)) if np.isfinite(var_mad) and dt_val > 0 else np.nan

        rows.append({
            "bin": int(x_id),
            "x_center": float(np.median(xk)),
            "n": n,
            "w_sum": float(np.sum(wk_use)) if wk_use is not None else float(n),
            "median_dx": float(med_dx),
            "mad_dx": float(mad),
            "sigma_dx": float(sigma),
            "var_mad_dx": float(var_mad),
            "D_mad": float(D),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("x_center").reset_index(drop=True)


def compute_xstar_slope(
    curve: pd.DataFrame,
    crossing_mode: str,
    window_bins: int,
    use_wls: bool,
) -> Dict[str, float]:
    if curve is None or curve.empty or curve.shape[0] < 2:
        return {"x_star": np.nan, "slope": np.nan, "n_win": 0}

    c = curve.sort_values("x_center")
    x = c["x_center"].to_numpy(float)
    y = c["median_dx"].to_numpy(float)
    w = c["w_sum"].to_numpy(float) if "w_sum" in c.columns else None

    i = choose_crossing_idx(x, y, crossing_mode)
    if i is None:
        return {"x_star": np.nan, "slope": np.nan, "n_win": 0}

    xs = interpolate_zero(x[i], y[i], x[i + 1], y[i + 1])
    if not np.isfinite(xs):
        return {"x_star": np.nan, "slope": np.nan, "n_win": 0}

    nearest = int(np.argmin(np.abs(x - xs)))
    idx = pick_window(len(x), nearest, window_bins)
    slope = slope_fit(x[idx], y[idx], w[idx] if w is not None else None, use_wls)
    return {"x_star": float(xs), "slope": float(slope), "n_win": int(idx.size)}


def pooled_var_mad(
    df: pd.DataFrame,
    weight_col: str,
    mad_sigma_c: float,
    robust_center_mode: str,
) -> float:
    dx = pd.to_numeric(df["dx"], errors="coerce").to_numpy(float)
    dx = dx[np.isfinite(dx)]
    if dx.size == 0:
        return np.nan

    if weight_col == "1" or weight_col not in df.columns:
        center = float(np.median(dx))
        mad = float(np.median(np.abs(dx - center)))
        sigma = mad_to_sigma(mad, mad_sigma_c)
        return float(sigma * sigma) if np.isfinite(sigma) else np.nan

    w = pd.to_numeric(df[weight_col], errors="coerce").to_numpy(float)
    n = min(dx.size, w.size)
    dx2 = dx[:n]
    w2 = w[:n]
    m = np.isfinite(dx2) & np.isfinite(w2)
    dx2, w2 = dx2[m], w2[m]

    if (not usable_weights(w2)) or robust_center_mode == "median":
        center = float(np.median(dx2))
        mad = float(np.median(np.abs(dx2 - center)))
    else:
        center = float(weighted_median(dx2, w2))
        mad = float(weighted_mad(dx2, w2, center=center))

    sigma = mad_to_sigma(mad, mad_sigma_c)
    return float(sigma * sigma) if np.isfinite(sigma) else np.nan


def make_cluster_key(df: pd.DataFrame, cols: List[str]) -> pd.Series:
    if len(cols) == 1:
        return df[cols[0]].astype(str).fillna("")
    s = df[cols[0]].astype(str).fillna("")
    for c in cols[1:]:
        s = s.str.cat(df[c].astype(str).fillna(""), sep="||")
    return s


def build_cluster_index_map(cluster_ids: np.ndarray) -> Dict[str, np.ndarray]:
    uniq = np.unique(cluster_ids)
    return {str(c): np.where(cluster_ids == c)[0] for c in uniq}


def sample_cluster_rows(rng: np.random.Generator, cluster_map: Dict[str, np.ndarray]) -> np.ndarray:
    keys = list(cluster_map.keys())
    if not keys:
        return np.array([], dtype=int)
    sampled = rng.choice(keys, size=len(keys), replace=True)
    return np.concatenate([cluster_map[k] for k in sampled]) if len(sampled) else np.array([], dtype=int)


# =========================================================
# Main blocks
# =========================================================
def run_local_block(
    d: pd.DataFrame,
    out_dir: Path,
    args,
    x_edges: np.ndarray,
    f_edges: np.ndarray,
) -> Dict[str, str]:
    counts = compute_dt_freqbin_counts(d, args.include_total_dt, x_edges, f_edges)
    drift = compute_drift_tables(d, args.weight_col, x_edges, args.min_rows_per_bin, args.n_boot, args.seed)
    msd = compute_msd_by_dt(d, args.weight_col, args.n_boot, args.seed)
    msd_x = compute_msd_by_dt_xbins(d, args.weight_col, x_edges, args.min_rows_per_bin, args.n_boot, args.seed)
    msd_f = compute_msd_by_dt_freqbins(d, args.weight_col, f_edges, args.min_rows_per_bin, args.n_boot, args.seed)

    counts.to_csv(out_dir / "dt_freqbin_counts.csv", index=False)
    drift.to_csv(out_dir / "drift_by_dt_xbins.csv", index=False)
    msd.to_csv(out_dir / "msd_by_dt.csv", index=False)
    msd_x.to_csv(out_dir / "msd_by_dt_xbins.csv", index=False)
    msd_f.to_csv(out_dir / "msd_by_dt_freqbins.csv", index=False)

    return {
        "dt_freqbin_counts.csv": "counts by dt and bins",
        "drift_by_dt_xbins.csv": "median_dx and P(dx>0) by dt,x-bin",
        "msd_by_dt.csv": "global MSD by dt",
        "msd_by_dt_xbins.csv": "MSD by dt,x-bin",
        "msd_by_dt_freqbins.csv": "MSD by dt,f-bin",
    }


def run_finite_time_block(
    d: pd.DataFrame,
    out_dir: Path,
    args,
    x_edges: np.ndarray,
) -> Dict[str, str]:
    cluster_cols = parse_list_str(args.cluster_cols)
    for c in cluster_cols:
        if c not in d.columns:
            raise ValueError(f"cluster column '{c}' not found")

    df = d.copy()
    df["_cluster_key"] = make_cluster_key(df, cluster_cols)

    point_rows = []
    curves_all = []
    edges_by_dt: Dict[int, np.ndarray] = {}

    dt_values = sorted(df["dt"].unique().tolist())
    for dt in dt_values:
        ddt = df[df["dt"] == int(dt)].copy()
        if ddt.empty:
            continue

        curve = build_moments_curve(
            ddt,
            dt_val=int(dt),
            weight_col=args.weight_col,
            x_edges=x_edges,
            mad_sigma_c=args.mad_to_sigma,
            robust_center_mode=args.robust_center_mode,
            min_rows_per_bin=args.min_rows_per_bin,
        )

        if curve.empty or curve.shape[0] < 2:
            continue

        edges_by_dt[int(dt)] = x_edges.copy()

        curve = curve.copy()
        curve.insert(0, "dt", int(dt))
        curves_all.append(curve)

        xs = compute_xstar_slope(
            curve,
            crossing_mode=args.crossing_mode,
            window_bins=args.window_bins,
            use_wls=bool(args.use_wls_slope and args.weight_col != "1"),
        )

        var_pool = pooled_var_mad(
            ddt,
            weight_col=args.weight_col,
            mad_sigma_c=args.mad_to_sigma,
            robust_center_mode=args.robust_center_mode,
        )
        D_pool = float(var_pool / (2.0 * float(dt))) if np.isfinite(var_pool) and dt > 0 else np.nan

        point_rows.append({
            "dt": int(dt),
            "x_star": xs["x_star"],
            "slope": xs["slope"],
            "abs_slope": abs(xs["slope"]) if np.isfinite(xs["slope"]) else np.nan,
            "tau_drift": tau_drift(float(dt), float(xs["slope"])),
            "n_win": int(xs["n_win"]),
            "n_rows": int(len(ddt)),
            "n_bins_kept": int(curve.shape[0]),
            "var_mad_pooled": float(var_pool),
            "D_mad_pooled": float(D_pool),
        })

    if not point_rows:
        pd.DataFrame().to_csv(out_dir / "finite_time_summary.csv", index=False)
        pd.DataFrame().to_csv(out_dir / "finite_time_by_dt_xbins_boot.csv", index=False)
        return {
            "finite_time_summary.csv": "empty (no valid finite-time bins)",
            "finite_time_by_dt_xbins_boot.csv": "empty (no valid finite-time bins)",
        }

    summary_df = pd.DataFrame(point_rows).sort_values("dt").reset_index(drop=True)
    curves_df = pd.concat(curves_all, ignore_index=True)

    for c in ["median_lo", "median_hi", "sigma_lo", "sigma_hi", "var_mad_lo", "var_mad_hi"]:
        if c not in curves_df.columns:
            curves_df[c] = np.nan

    B = int(args.bootstrap_reps_ft)
    ci = float(args.bootstrap_ci)
    alpha = 0.5 * (1.0 - ci)
    qlo = alpha
    qhi = 1.0 - alpha
    rng = np.random.default_rng(int(args.seed))

    boot_records = []

    for dt in summary_df["dt"].tolist():
        ddt = df[df["dt"] == int(dt)].copy()
        cluster_map = build_cluster_index_map(ddt["_cluster_key"].to_numpy(str))
        point_curve = curves_df[curves_df["dt"] == int(dt)].copy()
        kept_bins = point_curve["bin"].to_numpy(int)

        boot_med_by_bin = {int(b): [] for b in kept_bins.tolist()}
        boot_sig_by_bin = {int(b): [] for b in kept_bins.tolist()}
        boot_var_by_bin = {int(b): [] for b in kept_bins.tolist()}

        for rep in range(B):
            ridx = sample_cluster_rows(rng, cluster_map)
            if ridx.size == 0:
                boot_records.append({
                    "dt": int(dt), "rep": rep,
                    "x_star": np.nan, "slope": np.nan, "tau_drift": np.nan,
                    "var_mad_pooled": np.nan, "D_mad_pooled": np.nan,
                })
                continue

            bs = ddt.iloc[ridx].copy()

            curve_b = build_moments_curve(
                bs,
                dt_val=int(dt),
                weight_col=args.weight_col,
                x_edges=x_edges,
                mad_sigma_c=args.mad_to_sigma,
                robust_center_mode=args.robust_center_mode,
                min_rows_per_bin=args.min_rows_per_bin,
            )

            if curve_b.shape[0] >= 2:
                xs_b = compute_xstar_slope(
                    curve_b,
                    crossing_mode=args.crossing_mode,
                    window_bins=args.window_bins,
                    use_wls=bool(args.use_wls_slope and args.weight_col != "1"),
                )
                slope_b = float(xs_b["slope"])
                tau_b = tau_drift(float(dt), slope_b)
            else:
                xs_b = {"x_star": np.nan}
                slope_b = np.nan
                tau_b = np.nan

            var_pool_b = pooled_var_mad(
                bs,
                weight_col=args.weight_col,
                mad_sigma_c=args.mad_to_sigma,
                robust_center_mode=args.robust_center_mode,
            )
            D_pool_b = float(var_pool_b / (2.0 * float(dt))) if np.isfinite(var_pool_b) and dt > 0 else np.nan

            boot_records.append({
                "dt": int(dt), "rep": rep,
                "x_star": xs_b["x_star"], "slope": slope_b, "tau_drift": tau_b,
                "var_mad_pooled": float(var_pool_b), "D_mad_pooled": float(D_pool_b),
            })

            if not curve_b.empty:
                for _, r in curve_b.iterrows():
                    k = int(r["bin"])
                    if k in boot_med_by_bin:
                        boot_med_by_bin[k].append(float(r["median_dx"]))
                        boot_sig_by_bin[k].append(float(r["sigma_dx"]))
                        boot_var_by_bin[k].append(float(r["var_mad_dx"]))

        def qq(arr: List[float], q: float) -> float:
            a = np.asarray(arr, float)
            a = a[np.isfinite(a)]
            return float(np.quantile(a, q)) if a.size else np.nan

        idx_dt = curves_df.index[curves_df["dt"] == int(dt)].to_numpy()
        point_curve = curves_df.loc[idx_dt].copy()

        curves_df.loc[idx_dt, "median_lo"] = [qq(boot_med_by_bin.get(int(k), []), qlo) for k in point_curve["bin"]]
        curves_df.loc[idx_dt, "median_hi"] = [qq(boot_med_by_bin.get(int(k), []), qhi) for k in point_curve["bin"]]
        curves_df.loc[idx_dt, "sigma_lo"] = [qq(boot_sig_by_bin.get(int(k), []), qlo) for k in point_curve["bin"]]
        curves_df.loc[idx_dt, "sigma_hi"] = [qq(boot_sig_by_bin.get(int(k), []), qhi) for k in point_curve["bin"]]
        curves_df.loc[idx_dt, "var_mad_lo"] = [qq(boot_var_by_bin.get(int(k), []), qlo) for k in point_curve["bin"]]
        curves_df.loc[idx_dt, "var_mad_hi"] = [qq(boot_var_by_bin.get(int(k), []), qhi) for k in point_curve["bin"]]

    boot_df = pd.DataFrame(boot_records)

    ci_rows = []
    for dt in summary_df["dt"].tolist():
        bd = boot_df[boot_df["dt"] == int(dt)]

        def q(arr: np.ndarray, qv: float) -> float:
            a = np.asarray(arr, float)
            a = a[np.isfinite(a)]
            return float(np.quantile(a, qv)) if a.size else np.nan

        xstar = bd["x_star"].to_numpy(float)
        slope = bd["slope"].to_numpy(float)
        tau = bd["tau_drift"].to_numpy(float)
        vpool = bd["var_mad_pooled"].to_numpy(float)
        Dpool = bd["D_mad_pooled"].to_numpy(float)

        ci_rows.append({
            "dt": int(dt),
            "x_star_lo": q(xstar, qlo),
            "x_star_hi": q(xstar, qhi),
            "slope_lo": q(slope, qlo),
            "slope_hi": q(slope, qhi),
            "tau_drift_lo": q(tau, qlo),
            "tau_drift_hi": q(tau, qhi),
            "var_mad_pooled_lo": q(vpool, qlo),
            "var_mad_pooled_hi": q(vpool, qhi),
            "D_mad_pooled_lo": q(Dpool, qlo),
            "D_mad_pooled_hi": q(Dpool, qhi),
            "n_boot_valid_xstar": int(np.isfinite(xstar).sum()),
            "n_boot_valid_slope": int(np.isfinite(slope).sum()),
            "n_boot_valid_varpool": int(np.isfinite(vpool).sum()),
        })

    ci_df = pd.DataFrame(ci_rows)
    summary_df = summary_df.merge(ci_df, on="dt", how="left")
    summary_df["has_crossing"] = np.isfinite(summary_df["x_star"])
    summary_df["mean_reverting"] = np.isfinite(summary_df["slope"]) & (summary_df["slope"] < 0)

    summary_df.to_csv(out_dir / "finite_time_summary.csv", index=False)
    curves_df.to_csv(out_dir / "finite_time_by_dt_xbins_boot.csv", index=False)

    return {
        "finite_time_summary.csv": "finite-time pooled summaries by dt",
        "finite_time_by_dt_xbins_boot.csv": "finite-time binned curves with bootstrap bands",
    }


# =========================================================
# Main
# =========================================================
def main() -> None:
    args = build_argparser().parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    dt_values = parse_list_int(args.dt_values)
    classes = parse_list_str(args.classes)

    d = read_and_prepare(Path(args.input), x_base=args.x_base, sep=args.csv_sep)
    d = filter_dataframe(
        d,
        min_dt=args.min_dt,
        dt_max=args.dt_max,
        dt_values=dt_values,
        only_TT=args.only_TT,
        tt_value=args.tt_value,
        classes=classes,
    )

    if args.subsample_per_dt > 0:
        d = subsample_by_dt(d, args.subsample_per_dt, args.seed)

    if d.empty:
        raise ValueError("No rows remain after filtering.")

    if args.binning_mode == "save_edges_only":
        x_edges = compute_quantile_edges(d["ln_f0"].to_numpy(float), args.nbins_x)
        f_edges = compute_quantile_edges(d["ln_fgeo"].to_numpy(float), args.nbins_f)

        save_edges_csv(x_edges, out_dir / args.save_x_edges_name, "x_edge_ln")
        save_edges_csv(f_edges, out_dir / args.save_f_edges_name, "f_edge_ln")

        manifest = {
            "input": str(Path(args.input).resolve()),
            "mode": "save_edges_only",
            "x_base_in": args.x_base,
            "output_base": "ln",
            "n_rows_used": int(len(d)),
            "outputs": [args.save_x_edges_name, args.save_f_edges_name, "manifest.json"],
        }
        (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        print(f"[SAVED EDGES] {out_dir}")
        print(f"  - {args.save_x_edges_name}")
        print(f"  - {args.save_f_edges_name}")
        print("  - manifest.json")
        return

    x_edges, f_edges = define_edges(
        d=d,
        binning_mode=args.binning_mode,
        nbins_x=args.nbins_x,
        nbins_f=args.nbins_f,
        x_edges_csv=args.x_edges_csv,
        f_edges_csv=args.f_edges_csv,
    )

    save_edges_csv(x_edges, out_dir / args.save_x_edges_name, "x_edge_ln")
    save_edges_csv(f_edges, out_dir / args.save_f_edges_name, "f_edge_ln")

    produced = {}

    if args.mode in {"local", "all"}:
        produced.update(run_local_block(d, out_dir, args, x_edges, f_edges))

    if args.mode in {"finite_time", "all"}:
        produced.update(run_finite_time_block(d, out_dir, args, x_edges))

    manifest = {
        "input": str(Path(args.input).resolve()),
        "out_dir": str(out_dir),
        "mode": args.mode,
        "x_base_in": args.x_base,
        "output_base": "ln",
        "filters": {
            "min_dt": float(args.min_dt),
            "dt_max": int(args.dt_max),
            "dt_values": dt_values,
            "only_TT": bool(args.only_TT),
            "tt_value": args.tt_value,
            "classes": classes,
            "subsample_per_dt": int(args.subsample_per_dt),
        },
        "weights": {"weight_col": args.weight_col},
        "binning": {
            "binning_mode": args.binning_mode,
            "x_var": "ln_f0",
            "f_var": "ln_fgeo",
            "nbins_x": int(args.nbins_x),
            "nbins_f": int(args.nbins_f),
            "min_rows_per_bin": int(args.min_rows_per_bin),
            "x_edges_saved": args.save_x_edges_name,
            "f_edges_saved": args.save_f_edges_name,
            "x_edges_csv_in": str(Path(args.x_edges_csv).resolve()) if args.x_edges_csv else "",
            "f_edges_csv_in": str(Path(args.f_edges_csv).resolve()) if args.f_edges_csv else "",
        },
        "local_bootstrap": {
            "n_boot": int(args.n_boot),
            "seed": int(args.seed),
        },
        "finite_time": {
            "window_bins": int(args.window_bins),
            "crossing_mode": args.crossing_mode,
            "mad_to_sigma": float(args.mad_to_sigma),
            "robust_center_mode": args.robust_center_mode,
            "bootstrap_reps_ft": int(args.bootstrap_reps_ft),
            "bootstrap_ci": float(args.bootstrap_ci),
            "cluster_cols": parse_list_str(args.cluster_cols),
            "use_wls_slope": bool(args.use_wls_slope),
        },
        "n_rows_used": int(len(d)),
        "outputs": list(produced.keys()) + [args.save_x_edges_name, args.save_f_edges_name, "manifest.json"],
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    if args.debug:
        print("[DEBUG] rows used:", len(d))
        print("[DEBUG] dt values:", sorted(d["dt"].unique().tolist()))
        print("[DEBUG] x_edges ln:", x_edges)
        print("[DEBUG] f_edges ln:", f_edges)

    print(f"[SAVED DATA] {out_dir}")
    for f in manifest["outputs"]:
        print(f"  - {f}")


if __name__ == "__main__":
    main()