#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12-longitudinal_vs_pseudo_fluctuation_baseline.py
=================================================

Cross-dataset technical-baseline comparison for finite-time clonotype
fluctuations.

Scientific role
---------------
Step 12 compares genuine longitudinal fluctuation magnitude with the
pseudo-longitudinal technical baseline after the displacement representation and
conditioning geometry have been fixed upstream.

The analysis combines:

    Step 11
        genuine longitudinal fluctuation dynamics
        -> representation: latent-frequency displacement (`dx_latent`)
        -> conditioning: transition-centred latent abundance (`xstar_latent`)
        -> primary transition class: TT
        -> abundance-resolved Var(dx | xstar, dt), MSD and MAD

    Step 8
        pseudo-longitudinal representation assessment
        -> fixed downstream displacement representation: latent_frequency
        -> fixed conditioning geometry: xstar_latent
        -> Step-8 absolute xstar bin edges used for the pseudo reference

    pseudo Step 4
        latent trajectories used to construct directed pseudo-time pair
        sufficient statistics, unless a previously generated sufficient table is
        supplied directly.

Step 12 asks how the one-interval genuine longitudinal fluctuation magnitude
compares with the fluctuation magnitude expected from the pseudo-longitudinal
technical reference.

Primary estimand
----------------
    conditioning coordinate : xstar_latent
    displacement            : dx_latent
    transition class        : TT
    longitudinal lag        : dt = 1 (default)

Metrics
-------
    PRIMARY
        sample Var(dx_latent | xstar_latent)

    COMPLEMENTARY
        MSD = E(dx_latent^2 | xstar_latent)

    DESCRIPTIVE
        MAD(dx_latent | xstar_latent)

Var and MSD are evaluated across arbitrary pseudo-time configurations using
additive sufficient statistics. MAD is retained as a descriptive comparison
against the fixed pooled pseudo reference from Step 8. It is deliberately not
assigned a configuration-level permutation P value because exact pooled MAD is
not additive across pseudo-time pairs.

Absolute-abundance comparison rule
----------------------------------
Longitudinal and pseudo datasets retain their native absolute xstar abundance
geometry.

No percentile normalization, rank transformation or within-dataset abundance
normalization is applied. No common clonotype-level bin edges are imposed across
datasets before estimation.

Direct quantitative comparison is restricted to the absolute-abundance support
represented by both datasets. Already-estimated longitudinal and pseudo curves
are interpolated only within this shared support onto a common evaluation grid.
Underlying clonotypes are never re-binned across datasets, and no extrapolation
outside native support is used.

Pseudo configuration null
-------------------------
Pseudo-time orderings are generated independently within each pseudo subject.
For Var and MSD, a configuration and the configuration obtained by reversing
ALL subject-specific orders are exactly redundant because:

    xstar is endpoint symmetric;
    dx changes sign globally;
    dx^2 is unchanged;
    sample variance is unchanged.

Only one member of each global order/reverse equivalence class is therefore
retained.

Inputs
------
--longitudinal-dir
    Output directory from:

        11-longitudinal_fluctuation_dynamics.py

    Required primary curve:

        10_binned_dynamics_long.csv

    or its parquet equivalent.

--pseudo-representation-dir
    Output directory from:

        8-pseudo_longitudinal_representation_assessment.py

    Required:

        02_conditioning_bin_edges.csv

    Optional fixed pooled pseudo reference used for descriptive MAD:

        03_pooled_representation_curves.csv

Pseudo pair source: provide exactly ONE of

--pseudo-trajectories
    Pseudo Step-4 trajectory table:

        latent_trajectories_long.parquet

    normally produced by:

        4-latent_trajectory_construction.py

    The script builds directed TT pseudo-pair sufficient statistics once and
    then evaluates the pseudo-time configuration ensemble.

--pseudo-sufficient
    Previously generated:

        02_pseudo_directed_pair_sufficient.csv

    or parquet equivalent. This is the fast rerun path and avoids rebuilding
    directed-pair sufficient statistics.

Main outputs
------------
    00_run_config.json
    00_input_manifest.csv
    01_longitudinal_dt1_native_curve.csv
    02_pseudo_directed_pair_sufficient.csv/parquet
    03_pseudo_unique_orderings.csv
    04_pseudo_configuration_curves.csv
    05_pseudo_native_envelope.csv
    06_shared_absolute_support_summary.csv
    07_shared_grid_comparison.csv
    08_pseudo_shared_grid_configuration_metrics.csv
    09_empirical_tests.csv
    10_comparison_summary.csv
    11_pseudo_fixed_reference_curve.csv          [when available]
    12_fluctuation_baseline_report.json

Interpretive boundary
---------------------
The pseudo-longitudinal reference is an empirical technical fluctuation baseline.
Differences labelled "excess" are descriptive longitudinal-minus-pseudo
quantities. They are NOT interpreted as pure biological variance, pure
biological noise, or an exact technical/biological variance decomposition.

The pseudo baseline is NOT subtracted from the genuine longitudinal Step-11
fluctuation curves before Step-13 temporal-scaling analysis. Step 12 is therefore
a cross-dataset comparison layer, not a preprocessing correction for Step 13.

Pipeline position
-----------------
    Step 11 + Step 8 + pseudo Step 4  ->  Step 12

    Step 11                           ->  Step 13

Step 12 does not computationally feed Step 13.

Python compatibility
--------------------
Python >= 3.9
Required: numpy, pandas
Required only when --pseudo-trajectories is used: polars
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SCRIPT_VERSION = "v2-step12-current-provenance-2026-08-25"
PRIMARY_REPRESENTATION = "latent"
PRIMARY_CONDITIONING = "xstar"
PRIMARY_MODE = "TT"
PSEUDO_REFERENCE_REPRESENTATION = "latent_frequency"


# =============================================================================
# Generic helpers
# =============================================================================


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(path: Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if not np.isfinite(value):
            return None
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, np.ndarray):
        return [json_safe(v) for v in value.tolist()]
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError("Input table not found: {}".format(path))
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix in {".parquet", ".pq"}:
        try:
            return pd.read_parquet(path)
        except Exception:
            try:
                import polars as pl
                return pl.read_parquet(str(path)).to_pandas()
            except Exception as exc:
                raise RuntimeError(
                    "Could not read parquet {}. Install pyarrow or polars. {}".format(
                        path, exc
                    )
                )
    raise ValueError("Unsupported table format: {}".format(path))


def write_csv(df: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, index=False)
    print("[WRITE] {}".format(path))


def write_parquet_best_effort(df: pd.DataFrame, path: Path) -> bool:
    ensure_dir(path.parent)
    try:
        df.to_parquet(path, index=False)
        print("[WRITE] {}".format(path))
        return True
    except Exception:
        try:
            import polars as pl
            pl.from_pandas(df).write_parquet(str(path), compression="zstd")
            print("[WRITE] {}".format(path))
            return True
        except Exception as exc:
            eprint("[WARN] parquet write skipped for {}: {}".format(path, exc))
            return False


def find_table(directory: Path, stem: str) -> Path:
    for suffix in [".csv", ".parquet", ".pq"]:
        p = Path(directory) / (stem + suffix)
        if p.exists():
            return p
    raise FileNotFoundError(
        "Could not find {}.[csv|parquet|pq] in {}".format(stem, directory)
    )


def require_columns(df: pd.DataFrame, required: Sequence[str], label: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            "{} is missing columns {}. Available: {}".format(
                label, missing, list(df.columns)
            )
        )


def finite_numeric(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    return x[np.isfinite(x)]


def interpolate_curve(
    x: np.ndarray,
    y: np.ndarray,
    grid: np.ndarray,
) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    grid = np.asarray(grid, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if x.size < 2:
        return np.full(grid.shape, np.nan, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    ux, inv = np.unique(x, return_inverse=True)
    if ux.size != x.size:
        sums = np.zeros(ux.size, dtype=float)
        counts = np.zeros(ux.size, dtype=float)
        for i, g in enumerate(inv):
            if np.isfinite(y[i]):
                sums[g] += y[i]
                counts[g] += 1.0
        good = counts > 0
        x = ux[good]
        y = sums[good] / counts[good]
    if x.size < 2:
        return np.full(grid.shape, np.nan, dtype=float)
    out = np.full(grid.shape, np.nan, dtype=float)
    inside = (grid >= np.min(x)) & (grid <= np.max(x))
    out[inside] = np.interp(grid[inside], x, y)
    return out


def empirical_p_upper(null: np.ndarray, observed: float) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if null.size == 0 or not np.isfinite(observed):
        return float("nan")
    return float((1.0 + np.sum(null >= observed)) / (null.size + 1.0))


def empirical_p_lower(null: np.ndarray, observed: float) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if null.size == 0 or not np.isfinite(observed):
        return float("nan")
    return float((1.0 + np.sum(null <= observed)) / (null.size + 1.0))


def empirical_p_two_sided_about_median(null: np.ndarray, observed: float) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if null.size == 0 or not np.isfinite(observed):
        return float("nan")
    med = float(np.median(null))
    obs_dev = abs(float(observed) - med)
    p = (1.0 + np.sum(np.abs(null - med) >= obs_dev)) / (null.size + 1.0)
    return float(min(1.0, p))


# =============================================================================
# Longitudinal Step-11 input
# =============================================================================


def load_longitudinal_curve(
    longitudinal_dir: Path,
    dt_value: int,
) -> Tuple[pd.DataFrame, Path]:
    path = find_table(longitudinal_dir, "10_binned_dynamics_long")
    df = read_table(path)
    require_columns(
        df,
        [
            "representation", "conditioning", "mode", "dt", "bin_id", "bin_mid",
            "n", "n_subjects", "var_dx", "msd", "mad_dx",
        ],
        "Step-11 binned dynamics",
    )

    rep = df["representation"].astype(str).str.lower()
    cond = df["conditioning"].astype(str).str.lower()
    mode = df["mode"].astype(str).str.upper()
    dt = pd.to_numeric(df["dt"], errors="coerce")

    out = df.loc[
        (rep == PRIMARY_REPRESENTATION)
        & (cond == PRIMARY_CONDITIONING)
        & (mode == PRIMARY_MODE)
        & np.isfinite(dt)
        & (np.abs(dt - float(dt_value)) < 1e-12)
    ].copy()

    if out.empty:
        raise ValueError(
            "No Step-11 rows for representation={}, conditioning={}, mode={}, dt={}".format(
                PRIMARY_REPRESENTATION,
                PRIMARY_CONDITIONING,
                PRIMARY_MODE,
                dt_value,
            )
        )

    out = out.sort_values("bin_mid").reset_index(drop=True)
    return out, path


# =============================================================================
# Pseudo Step-8 fixed xstar bins and descriptive reference
# =============================================================================


def load_pseudo_edges(pseudo_representation_dir: Path) -> Tuple[np.ndarray, pd.DataFrame, Path]:
    path = Path(pseudo_representation_dir) / "02_conditioning_bin_edges.csv"
    edge = read_table(path)
    require_columns(edge, ["bin_id", "x_lo", "x_hi", "x_center"], str(path))
    edge = edge.copy()
    for c in ["bin_id", "x_lo", "x_hi", "x_center"]:
        edge[c] = pd.to_numeric(edge[c], errors="coerce")
    edge = edge.loc[
        np.isfinite(edge["bin_id"])
        & np.isfinite(edge["x_lo"])
        & np.isfinite(edge["x_hi"])
        & np.isfinite(edge["x_center"])
    ].sort_values("x_lo").reset_index(drop=True)
    if edge.empty:
        raise ValueError("Pseudo xstar edge table is empty after numeric filtering")

    if not np.all(edge["x_hi"].to_numpy() > edge["x_lo"].to_numpy()):
        raise ValueError("Pseudo bin edges contain non-positive widths")

    # Require contiguous bins up to floating tolerance.
    if len(edge) > 1:
        gap = edge["x_lo"].to_numpy()[1:] - edge["x_hi"].to_numpy()[:-1]
        if np.max(np.abs(gap)) > 1e-7:
            raise ValueError("Pseudo xstar bins are not contiguous")

    edges = np.asarray(
        [float(edge.loc[0, "x_lo"])] + edge["x_hi"].astype(float).tolist(),
        dtype=float,
    )
    if np.any(np.diff(edges) <= 0):
        raise ValueError("Pseudo xstar edges are not strictly increasing")

    edge = edge.rename(columns={"x_lo": "bin_left", "x_hi": "bin_right", "x_center": "bin_mid"})
    edge["bin_index"] = np.arange(len(edge), dtype=int)
    edge["bin_id"] = edge["bin_id"].astype(int)
    return edges, edge[["bin_index", "bin_id", "bin_left", "bin_mid", "bin_right"]], path


def load_pseudo_fixed_reference(pseudo_representation_dir: Path) -> Tuple[pd.DataFrame, Optional[Path]]:
    path = Path(pseudo_representation_dir) / "03_pooled_representation_curves.csv"
    if not path.exists():
        return pd.DataFrame(), None
    df = read_table(path)
    required = ["representation", "bin_id", "x_center", "var_dx", "msd", "mad_dx"]
    require_columns(df, required, str(path))
    out = df.loc[
        df["representation"].astype(str) == PSEUDO_REFERENCE_REPRESENTATION
    ].copy()
    if "meets_min_n" in out.columns:
        flag = out["meets_min_n"]
        if not pd.api.types.is_bool_dtype(flag):
            flag = flag.astype(str).str.lower().isin(["true", "1", "yes", "y"])
        out = out.loc[flag].copy()
    out = out.sort_values("x_center").reset_index(drop=True)
    return out, path


# =============================================================================
# Pseudo directed-pair construction
# =============================================================================


def assign_bins(x: np.ndarray, edges: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.full(x.shape, -1, dtype=np.int16)
    finite = np.isfinite(x)
    if not np.any(finite):
        return out
    vals = x[finite]
    ids = np.digitize(vals, edges[1:-1], right=False)
    ids = np.clip(ids, 0, len(edges) - 2)
    inside = (vals >= edges[0]) & (vals <= edges[-1])
    out_vals = np.where(inside, ids, -1).astype(np.int16)
    out[finite] = out_vals
    return out


def summarize_pair_arrays(
    xstar: np.ndarray,
    dx: np.ndarray,
    edges: np.ndarray,
) -> pd.DataFrame:
    xstar = np.asarray(xstar, dtype=float)
    dx = np.asarray(dx, dtype=float)
    bins = assign_bins(xstar, edges)
    valid = (bins >= 0) & np.isfinite(xstar) & np.isfinite(dx)
    bins = bins[valid].astype(np.int64)
    xv = xstar[valid]
    dv = dx[valid]
    n_bins = len(edges) - 1

    if bins.size == 0:
        return pd.DataFrame(
            {
                "bin_index": np.arange(n_bins, dtype=int),
                "n": np.zeros(n_bins, dtype=np.int64),
                "sum_x": np.zeros(n_bins, dtype=float),
                "sum_dx": np.zeros(n_bins, dtype=float),
                "sum_dx2": np.zeros(n_bins, dtype=float),
            }
        )

    return pd.DataFrame(
        {
            "bin_index": np.arange(n_bins, dtype=int),
            "n": np.bincount(bins, minlength=n_bins).astype(np.int64),
            "sum_x": np.bincount(bins, weights=xv, minlength=n_bins).astype(float),
            "sum_dx": np.bincount(bins, weights=dv, minlength=n_bins).astype(float),
            "sum_dx2": np.bincount(bins, weights=dv * dv, minlength=n_bins).astype(float),
        }
    )


def normalize_bool_expr(pl, column: str):
    s = pl.col(column)
    return (
        pl.when(s.cast(pl.Boolean, strict=False).is_not_null())
        .then(s.cast(pl.Boolean, strict=False).fill_null(False))
        .otherwise(
            s.cast(pl.String, strict=False)
            .str.strip_chars()
            .str.to_lowercase()
            .is_in(["true", "t", "1", "yes", "y", "observable", "detected", "detectable"])
            .fill_null(False)
        )
    )


def collect_streaming(lf):
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)


def discover_pseudo_times(trajectories: Path) -> Dict[int, List[int]]:
    try:
        import polars as pl
    except Exception as exc:
        raise RuntimeError("polars is required with --pseudo-trajectories") from exc

    schema = set(pl.scan_parquet(str(trajectories)).collect_schema().names())
    for col in ["subject", "time"]:
        if col not in schema:
            raise ValueError("Pseudo trajectory table missing column {!r}".format(col))

    df = collect_streaming(
        pl.scan_parquet(str(trajectories))
        .select(
            pl.col("subject").cast(pl.Int64, strict=False),
            pl.col("time").cast(pl.Int64, strict=False),
        )
        .drop_nulls()
        .unique()
        .sort(["subject", "time"])
    )
    out: Dict[int, List[int]] = {}
    for row in df.iter_rows(named=True):
        out.setdefault(int(row["subject"]), []).append(int(row["time"]))
    for subject in list(out):
        out[subject] = sorted(set(out[subject]))
    return out


def build_pseudo_directed_pair_sufficient(
    trajectories: Path,
    edges: np.ndarray,
    edge_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[int, List[int]]]:
    try:
        import polars as pl
    except Exception as exc:
        raise RuntimeError("polars is required with --pseudo-trajectories") from exc

    trajectories = Path(trajectories)
    if not trajectories.exists():
        raise FileNotFoundError("Pseudo trajectories not found: {}".format(trajectories))

    schema = set(pl.scan_parquet(str(trajectories)).collect_schema().names())
    required = [
        "subject", "time", "aaSeqCDR3", "observable", "x_latent", "x_latent_sd"
    ]
    missing = [c for c in required if c not in schema]
    if missing:
        raise ValueError(
            "Pseudo trajectory table missing required columns {}".format(missing)
        )

    times_by_subject = discover_pseudo_times(trajectories)
    outputs: List[pd.DataFrame] = []

    total_pairs = sum(len(v) * (len(v) - 1) for v in times_by_subject.values())
    pair_counter = 0

    for subject in sorted(times_by_subject):
        times = times_by_subject[subject]
        if len(times) < 2:
            continue
        print("[PSEUDO] subject {} times {}".format(subject, times))

        subject_df = collect_streaming(
            pl.scan_parquet(str(trajectories))
            .filter(pl.col("subject").cast(pl.Int64, strict=False) == int(subject))
            .select([pl.col(c) for c in required])
            .with_columns(
                [
                    pl.col("time").cast(pl.Int64, strict=False),
                    pl.col("aaSeqCDR3").cast(pl.String, strict=False),
                    normalize_bool_expr(pl, "observable").alias("observable"),
                    pl.col("x_latent").cast(pl.Float64, strict=False),
                    pl.col("x_latent_sd").cast(pl.Float64, strict=False),
                ]
            )
        )

        frames: Dict[int, object] = {}
        for t in times:
            frame = subject_df.filter(pl.col("time") == int(t)).select(
                ["aaSeqCDR3", "observable", "x_latent", "x_latent_sd"]
            )
            # Trajectory tables should be unique by subject x time x clonotype.
            if frame["aaSeqCDR3"].n_unique() != frame.height:
                raise ValueError(
                    "Duplicate aaSeqCDR3 rows for pseudo subject={} time={}".format(
                        subject, t
                    )
                )
            frames[int(t)] = frame

        for t0 in times:
            for t1 in times:
                if int(t0) == int(t1):
                    continue
                pair_counter += 1
                print(
                    "  [pair {}/{}] S{} {}->{}".format(
                        pair_counter, total_pairs, subject, t0, t1
                    )
                )

                left = frames[int(t0)].rename(
                    {
                        "observable": "observable_t0",
                        "x_latent": "x0",
                        "x_latent_sd": "sd0",
                    }
                )
                right = frames[int(t1)].rename(
                    {
                        "observable": "observable_t1",
                        "x_latent": "x1",
                        "x_latent_sd": "sd1",
                    }
                )
                joined = left.join(right, on="aaSeqCDR3", how="inner")
                joined = joined.filter(
                    pl.col("observable_t0").fill_null(False)
                    & pl.col("observable_t1").fill_null(False)
                )

                if joined.height == 0:
                    pair = summarize_pair_arrays(
                        np.asarray([], dtype=float),
                        np.asarray([], dtype=float),
                        edges,
                    )
                else:
                    x0 = joined["x0"].to_numpy().astype(float, copy=False)
                    x1 = joined["x1"].to_numpy().astype(float, copy=False)
                    sd0 = joined["sd0"].to_numpy().astype(float, copy=False)
                    sd1 = joined["sd1"].to_numpy().astype(float, copy=False)
                    dx = x1 - x0

                    xstar = 0.5 * (x0 + x1)
                    good_sd = (
                        np.isfinite(sd0)
                        & np.isfinite(sd1)
                        & (sd0 > 0)
                        & (sd1 > 0)
                    )
                    if np.any(good_sd):
                        w0 = 1.0 / (sd0[good_sd] ** 2 + 1e-12)
                        w1 = 1.0 / (sd1[good_sd] ** 2 + 1e-12)
                        xstar[good_sd] = (
                            w0 * x0[good_sd] + w1 * x1[good_sd]
                        ) / (w0 + w1)

                    pair = summarize_pair_arrays(xstar, dx, edges)

                pair.insert(0, "t1", int(t1))
                pair.insert(0, "t0", int(t0))
                pair.insert(0, "subject", int(subject))
                outputs.append(pair)

        del frames, subject_df

    if not outputs:
        raise ValueError("No directed pseudo pairs were constructed")

    out = pd.concat(outputs, ignore_index=True)
    out = out.merge(edge_df, on="bin_index", how="left")
    out = out[
        [
            "subject", "t0", "t1", "bin_index", "bin_id",
            "bin_left", "bin_mid", "bin_right",
            "n", "sum_x", "sum_dx", "sum_dx2",
        ]
    ].sort_values(["subject", "t0", "t1", "bin_index"])
    return out.reset_index(drop=True), times_by_subject


def infer_times_from_sufficient(sufficient: pd.DataFrame) -> Dict[int, List[int]]:
    require_columns(sufficient, ["subject", "t0", "t1"], "pseudo sufficient table")
    out: Dict[int, set] = {}
    for row in sufficient[["subject", "t0", "t1"]].drop_duplicates().itertuples(index=False):
        s = int(row.subject)
        out.setdefault(s, set()).add(int(row.t0))
        out.setdefault(s, set()).add(int(row.t1))
    return {s: sorted(v) for s, v in sorted(out.items())}


# =============================================================================
# Dense sufficient-statistics index and pseudo configurations
# =============================================================================


class PairSufficientIndex:
    STAT_COLUMNS = ["n", "sum_x", "sum_dx", "sum_dx2"]

    def __init__(self, sufficient: pd.DataFrame, edge_df: pd.DataFrame):
        required = {
            "subject", "t0", "t1", "bin_index", "n", "sum_x", "sum_dx", "sum_dx2"
        }
        missing = required - set(sufficient.columns)
        if missing:
            raise ValueError(
                "Pseudo sufficient table missing columns {}".format(sorted(missing))
            )
        self.edge_df = edge_df.sort_values("bin_index").reset_index(drop=True).copy()
        self.n_bins = len(self.edge_df)
        self.data: Dict[Tuple[int, int, int], np.ndarray] = {}

        for key, g in sufficient.groupby(["subject", "t0", "t1"], sort=False):
            subject, t0, t1 = key
            arr = np.zeros((self.n_bins, 4), dtype=np.float64)
            bins = pd.to_numeric(g["bin_index"], errors="coerce").to_numpy(dtype=float)
            good = np.isfinite(bins)
            bins_i = bins[good].astype(int)
            ok = (bins_i >= 0) & (bins_i < self.n_bins)
            bins_i = bins_i[ok]
            if bins_i.size:
                vals = g.loc[good, self.STAT_COLUMNS].to_numpy(dtype=float)[ok]
                arr[bins_i, :] = vals
            self.data[(int(subject), int(t0), int(t1))] = arr

    def pair(self, subject: int, t0: int, t1: int) -> np.ndarray:
        arr = self.data.get((int(subject), int(t0), int(t1)))
        if arr is None:
            return np.zeros((self.n_bins, 4), dtype=np.float64)
        return arr

    def subject_order(self, subject: int, order: Sequence[int]) -> np.ndarray:
        total = np.zeros((self.n_bins, 4), dtype=np.float64)
        for t0, t1 in zip(order[:-1], order[1:]):
            total += self.pair(int(subject), int(t0), int(t1))
        return total

    def cohort_order(
        self,
        orders_by_subject: Dict[int, Sequence[int]],
    ) -> Tuple[np.ndarray, np.ndarray]:
        total = np.zeros((self.n_bins, 4), dtype=np.float64)
        subject_presence = np.zeros((self.n_bins,), dtype=np.int64)
        for subject, order in orders_by_subject.items():
            s = self.subject_order(int(subject), order)
            total += s
            subject_presence += (s[:, 0] > 0).astype(np.int64)
        return total, subject_presence


def configuration_key(orders_by_subject: Dict[int, Sequence[int]]) -> Tuple:
    direct = tuple(
        (int(s), tuple(int(t) for t in orders_by_subject[s]))
        for s in sorted(orders_by_subject)
    )
    reverse = tuple(
        (int(s), tuple(int(t) for t in reversed(orders_by_subject[s])))
        for s in sorted(orders_by_subject)
    )
    return min(direct, reverse)


def generate_unique_orderings(
    times_by_subject: Dict[int, Sequence[int]],
    n_configurations: int,
    seed: int,
) -> List[Dict[int, Tuple[int, ...]]]:
    if n_configurations < 1:
        raise ValueError("n_configurations must be >=1")
    for subject, times in times_by_subject.items():
        if len(times) < 2:
            raise ValueError("Pseudo subject {} has fewer than 2 times".format(subject))

    rng = np.random.default_rng(int(seed))
    seen = set()
    out: List[Dict[int, Tuple[int, ...]]] = []
    attempts = 0
    max_attempts = max(10000, int(n_configurations) * 100)

    while len(out) < int(n_configurations) and attempts < max_attempts:
        attempts += 1
        cfg: Dict[int, Tuple[int, ...]] = {}
        for subject in sorted(times_by_subject):
            times = np.asarray(list(times_by_subject[subject]), dtype=int)
            perm = tuple(int(v) for v in rng.permutation(times).tolist())
            cfg[int(subject)] = perm
        key = configuration_key(cfg)
        if key in seen:
            continue
        seen.add(key)
        out.append(cfg)

    if len(out) < int(n_configurations):
        raise RuntimeError(
            "Could generate only {} unique global order/reverse equivalence classes "
            "after {} attempts; requested {}".format(
                len(out), attempts, n_configurations
            )
        )
    return out


def sufficient_to_curve(
    stats: np.ndarray,
    n_subjects: np.ndarray,
    edge_df: pd.DataFrame,
    min_n: int,
    min_subjects: int,
) -> pd.DataFrame:
    stats = np.asarray(stats, dtype=float)
    n = stats[:, 0]
    sum_x = stats[:, 1]
    sum_dx = stats[:, 2]
    sum_dx2 = stats[:, 3]

    x_mean = np.full(n.shape, np.nan, dtype=float)
    mean_dx = np.full(n.shape, np.nan, dtype=float)
    msd = np.full(n.shape, np.nan, dtype=float)
    var_dx = np.full(n.shape, np.nan, dtype=float)

    positive = n > 0
    x_mean[positive] = sum_x[positive] / n[positive]
    mean_dx[positive] = sum_dx[positive] / n[positive]
    msd[positive] = sum_dx2[positive] / n[positive]

    enough = n > 1
    if np.any(enough):
        var = (
            sum_dx2[enough] - (sum_dx[enough] ** 2) / n[enough]
        ) / (n[enough] - 1.0)
        var_dx[enough] = np.maximum(var, 0.0)

    out = edge_df.copy()
    out["n"] = np.rint(n).astype(np.int64)
    out["n_subjects"] = np.asarray(n_subjects, dtype=np.int64)
    out["x_mean"] = x_mean
    out["mean_dx"] = mean_dx
    out["msd"] = msd
    out["var_dx"] = var_dx
    out["valid_bin"] = (
        (out["n"] >= int(min_n))
        & (out["n_subjects"] >= int(min_subjects))
        & np.isfinite(out["var_dx"])
        & np.isfinite(out["msd"])
    )
    return out


def evaluate_pseudo_configurations(
    index: PairSufficientIndex,
    configurations: Sequence[Dict[int, Sequence[int]]],
    min_n: int,
    min_subjects: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    curve_rows: List[pd.DataFrame] = []
    order_rows: List[Dict[str, object]] = []

    for config_id, cfg in enumerate(configurations):
        stats, presence = index.cohort_order(cfg)
        curve = sufficient_to_curve(
            stats,
            presence,
            index.edge_df,
            min_n=min_n,
            min_subjects=min_subjects,
        )
        curve.insert(0, "configuration_id", int(config_id))
        curve_rows.append(curve)

        order_rows.append(
            {
                "configuration_id": int(config_id),
                "orders_by_subject_json": json.dumps(
                    {str(s): list(cfg[s]) for s in sorted(cfg)}, sort_keys=True
                ),
                "global_reverse_redundant": True,
            }
        )

        if (config_id + 1) % 100 == 0 or config_id == 0:
            print(
                "[PSEUDO] evaluated configuration {}/{}".format(
                    config_id + 1, len(configurations)
                )
            )

    curves = pd.concat(curve_rows, ignore_index=True)
    orders = pd.DataFrame(order_rows)
    return curves, orders


# =============================================================================
# Pseudo native envelope and shared-support comparison
# =============================================================================


def quantile_safe(x: pd.Series, q: float) -> float:
    a = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return float(np.quantile(a, q))


def build_pseudo_native_envelope(
    config_curves: pd.DataFrame,
    n_configurations: int,
    min_config_fraction: float,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for bin_index, g in config_curves.groupby("bin_index", sort=True):
        g = g.copy()
        valid = g["valid_bin"].astype(bool)
        valid_g = g.loc[valid]
        n_valid = int(valid_g["configuration_id"].nunique())
        frac = float(n_valid / float(n_configurations)) if n_configurations else 0.0
        first = g.iloc[0]
        row: Dict[str, object] = {
            "bin_index": int(bin_index),
            "bin_id": int(first["bin_id"]),
            "bin_left": float(first["bin_left"]),
            "bin_mid": float(first["bin_mid"]),
            "bin_right": float(first["bin_right"]),
            "n_configurations_valid": n_valid,
            "configuration_valid_fraction": frac,
            "stable_bin": bool(frac >= float(min_config_fraction)),
            "median_n": quantile_safe(valid_g["n"], 0.5),
            "median_n_subjects": quantile_safe(valid_g["n_subjects"], 0.5),
        }
        for metric in ["var_dx", "msd", "mean_dx"]:
            for q, suffix in [
                (0.005, "q005"),
                (0.025, "q025"),
                (0.25, "q250"),
                (0.5, "median"),
                (0.75, "q750"),
                (0.975, "q975"),
                (0.995, "q995"),
            ]:
                row["{}_{}".format(metric, suffix)] = quantile_safe(valid_g[metric], q)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("bin_mid").reset_index(drop=True)


def determine_shared_support(
    longitudinal: pd.DataFrame,
    pseudo_envelope: pd.DataFrame,
) -> Dict[str, float]:
    lx = pd.to_numeric(longitudinal["bin_mid"], errors="coerce").to_numpy(dtype=float)
    lx = lx[np.isfinite(lx)]
    p = pseudo_envelope.loc[pseudo_envelope["stable_bin"].astype(bool)].copy()
    px = pd.to_numeric(p["bin_mid"], errors="coerce").to_numpy(dtype=float)
    px = px[np.isfinite(px)]
    if lx.size < 2:
        raise ValueError("Longitudinal dt=1 curve has fewer than 2 finite xstar anchors")
    if px.size < 2:
        raise ValueError("Pseudo stable envelope has fewer than 2 finite xstar anchors")

    long_min = float(np.min(lx))
    long_max = float(np.max(lx))
    pseudo_min = float(np.min(px))
    pseudo_max = float(np.max(px))
    shared_min = max(long_min, pseudo_min)
    shared_max = min(long_max, pseudo_max)
    if not shared_max > shared_min:
        raise ValueError(
            "No shared absolute xstar support: longitudinal=[{},{}], pseudo=[{},{}]".format(
                long_min, long_max, pseudo_min, pseudo_max
            )
        )
    return {
        "longitudinal_min": long_min,
        "longitudinal_max": long_max,
        "pseudo_stable_min": pseudo_min,
        "pseudo_stable_max": pseudo_max,
        "shared_min": shared_min,
        "shared_max": shared_max,
        "shared_width": shared_max - shared_min,
    }


def interpolate_longitudinal_to_grid(
    longitudinal: pd.DataFrame,
    grid: np.ndarray,
) -> Dict[str, np.ndarray]:
    x = pd.to_numeric(longitudinal["bin_mid"], errors="coerce").to_numpy(dtype=float)
    out: Dict[str, np.ndarray] = {}
    for metric in ["var_dx", "msd", "mad_dx", "mean_dx"]:
        y = pd.to_numeric(longitudinal[metric], errors="coerce").to_numpy(dtype=float)
        out[metric] = interpolate_curve(x, y, grid)

    ci_map = {
        "var_dx": ("var_dx_ci_low", "var_dx_ci_high"),
        "msd": ("msd_ci_low", "msd_ci_high"),
        "mean_dx": ("mean_dx_ci_low", "mean_dx_ci_high"),
    }
    for metric, (low_col, high_col) in ci_map.items():
        if low_col in longitudinal.columns:
            out["{}_ci_low".format(metric)] = interpolate_curve(
                x,
                pd.to_numeric(longitudinal[low_col], errors="coerce").to_numpy(dtype=float),
                grid,
            )
        else:
            out["{}_ci_low".format(metric)] = np.full(grid.shape, np.nan)
        if high_col in longitudinal.columns:
            out["{}_ci_high".format(metric)] = interpolate_curve(
                x,
                pd.to_numeric(longitudinal[high_col], errors="coerce").to_numpy(dtype=float),
                grid,
            )
        else:
            out["{}_ci_high".format(metric)] = np.full(grid.shape, np.nan)
    return out


def interpolate_pseudo_configs_to_grid(
    config_curves: pd.DataFrame,
    grid: np.ndarray,
    min_grid_fraction: float,
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    rows: List[Dict[str, object]] = []
    grid_values: Dict[str, List[np.ndarray]] = {"var_dx": [], "msd": []}

    for config_id, g in config_curves.groupby("configuration_id", sort=True):
        d = g.loc[g["valid_bin"].astype(bool)].sort_values("bin_mid")
        x = pd.to_numeric(d["bin_mid"], errors="coerce").to_numpy(dtype=float)
        record: Dict[str, object] = {"configuration_id": int(config_id)}
        keep_configuration = True

        for metric in ["var_dx", "msd"]:
            y = pd.to_numeric(d[metric], errors="coerce").to_numpy(dtype=float)
            values = interpolate_curve(x, y, grid)
            finite = np.isfinite(values)
            frac = float(np.mean(finite)) if values.size else 0.0
            record["{}_grid_fraction".format(metric)] = frac
            if frac < float(min_grid_fraction):
                keep_configuration = False
            record["mean_{}_shared".format(metric)] = (
                float(np.nanmean(values)) if np.any(finite) else float("nan")
            )
            grid_values[metric].append(values)

        record["valid_for_integrated_test"] = bool(keep_configuration)
        rows.append(record)

    metrics = pd.DataFrame(rows)
    arrays = {
        metric: np.vstack(values) if values else np.empty((0, len(grid)))
        for metric, values in grid_values.items()
    }
    return metrics, arrays


def interpolate_pseudo_reference_mad(
    pseudo_reference: pd.DataFrame,
    grid: np.ndarray,
) -> np.ndarray:
    if pseudo_reference.empty:
        return np.full(grid.shape, np.nan, dtype=float)
    x = pd.to_numeric(pseudo_reference["x_center"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(pseudo_reference["mad_dx"], errors="coerce").to_numpy(dtype=float)
    return interpolate_curve(x, y, grid)


def build_shared_grid_table(
    grid: np.ndarray,
    longitudinal_interp: Dict[str, np.ndarray],
    pseudo_arrays: Dict[str, np.ndarray],
    pseudo_reference_mad: np.ndarray,
) -> pd.DataFrame:
    out = pd.DataFrame({"xstar": grid})

    for metric in ["var_dx", "msd"]:
        long_values = longitudinal_interp[metric]
        arr = pseudo_arrays[metric]
        if arr.shape[0] == 0:
            q025 = med = q975 = np.full(grid.shape, np.nan)
            n_valid = np.zeros(grid.shape, dtype=int)
        else:
            q025 = np.nanquantile(arr, 0.025, axis=0)
            med = np.nanquantile(arr, 0.5, axis=0)
            q975 = np.nanquantile(arr, 0.975, axis=0)
            n_valid = np.sum(np.isfinite(arr), axis=0)

        out["longitudinal_{}".format(metric)] = long_values
        out["longitudinal_{}_ci_low".format(metric)] = longitudinal_interp.get(
            "{}_ci_low".format(metric), np.full(grid.shape, np.nan)
        )
        out["longitudinal_{}_ci_high".format(metric)] = longitudinal_interp.get(
            "{}_ci_high".format(metric), np.full(grid.shape, np.nan)
        )
        out["pseudo_{}_median".format(metric)] = med
        out["pseudo_{}_q025".format(metric)] = q025
        out["pseudo_{}_q975".format(metric)] = q975
        out["pseudo_{}_n_configurations".format(metric)] = n_valid
        out["excess_{}".format(metric)] = long_values - med
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = long_values / med
        ratio[~np.isfinite(ratio)] = np.nan
        out["longitudinal_to_pseudo_{}_ratio".format(metric)] = ratio
        out["longitudinal_outside_pseudo_{}_95".format(metric)] = (
            (long_values < q025) | (long_values > q975)
        ) & np.isfinite(long_values) & np.isfinite(q025) & np.isfinite(q975)

    out["longitudinal_mad_dx"] = longitudinal_interp["mad_dx"]
    out["pseudo_mad_dx_fixed_reference"] = pseudo_reference_mad
    out["excess_mad_dx_descriptive"] = (
        out["longitudinal_mad_dx"] - out["pseudo_mad_dx_fixed_reference"]
    )
    return out


def build_empirical_tests(
    shared_grid: pd.DataFrame,
    config_metrics: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    valid_cfg = config_metrics.loc[config_metrics["valid_for_integrated_test"].astype(bool)].copy()

    for metric in ["var_dx", "msd"]:
        long = pd.to_numeric(
            shared_grid["longitudinal_{}".format(metric)], errors="coerce"
        ).to_numpy(dtype=float)
        observed = float(np.nanmean(long))
        null = pd.to_numeric(
            valid_cfg["mean_{}_shared".format(metric)], errors="coerce"
        ).to_numpy(dtype=float)
        null = null[np.isfinite(null)]
        if null.size == 0:
            summary = [float("nan")] * 7
        else:
            summary = [
                float(np.mean(null)),
                float(np.std(null, ddof=1)) if null.size > 1 else 0.0,
                float(np.quantile(null, 0.005)),
                float(np.quantile(null, 0.025)),
                float(np.median(null)),
                float(np.quantile(null, 0.975)),
                float(np.quantile(null, 0.995)),
            ]

        outside = shared_grid[
            "longitudinal_outside_pseudo_{}_95".format(metric)
        ].astype(bool)
        rows.append(
            {
                "metric": metric,
                "longitudinal_mean_over_shared_support": observed,
                "pseudo_n_configurations": int(null.size),
                "pseudo_mean": summary[0],
                "pseudo_sd": summary[1],
                "pseudo_q005": summary[2],
                "pseudo_q025": summary[3],
                "pseudo_median": summary[4],
                "pseudo_q975": summary[5],
                "pseudo_q995": summary[6],
                "empirical_p_upper": empirical_p_upper(null, observed),
                "empirical_p_lower": empirical_p_lower(null, observed),
                "empirical_p_two_sided_about_median": empirical_p_two_sided_about_median(
                    null, observed
                ),
                "n_shared_grid_points": int(len(shared_grid)),
                "n_longitudinal_outside_pseudo_95": int(outside.sum()),
                "fraction_longitudinal_outside_pseudo_95": float(outside.mean()),
            }
        )
    return pd.DataFrame(rows)


# =============================================================================
# CLI and main
# =============================================================================


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Step 12: compare genuine longitudinal latent xstar-conditioned "
            "fluctuation magnitude with the pseudo-longitudinal technical baseline "
            "on shared absolute abundance support."
        ),
    )

    ap.add_argument(
        "--longitudinal-dir",
        required=True,
        help=(
            "Step-11 output directory from "
            "11-longitudinal_fluctuation_dynamics.py; must contain "
            "10_binned_dynamics_long.csv/parquet."
        ),
    )
    ap.add_argument(
        "--pseudo-representation-dir",
        required=True,
        help=(
            "Step-8 output directory from "
            "8-pseudo_longitudinal_representation_assessment.py; must contain "
            "02_conditioning_bin_edges.csv. "
            "03_pooled_representation_curves.csv is optional for descriptive MAD."
        ),
    )
    ap.add_argument(
        "--outdir",
        required=True,
        help="Step-12 output directory.",
    )

    pseudo = ap.add_mutually_exclusive_group(required=True)
    pseudo.add_argument(
        "--pseudo-trajectories",
        default=None,
        help=(
            "Pseudo Step-4 latent_trajectories_long.parquet produced by "
            "4-latent_trajectory_construction.py. Step 12 builds directed TT "
            "pseudo-pair sufficient statistics from this table."
        ),
    )
    pseudo.add_argument(
        "--pseudo-sufficient",
        default=None,
        help=(
            "Precomputed 02_pseudo_directed_pair_sufficient.csv/parquet from a "
            "previous Step-12 run. This is the fast rerun path."
        ),
    )

    ap.add_argument(
        "--longitudinal-dt",
        type=int,
        default=1,
        help=(
            "Longitudinal Step-11 lag used for the cross-dataset comparison. "
            "The manuscript-primary comparison uses dt=1."
        ),
    )
    ap.add_argument(
        "--n-configurations",
        type=int,
        default=2000,
        help=(
            "Number of pseudo-time configurations sampled after global "
            "order/reversal deduplication."
        ),
    )
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument(
        "--pseudo-min-n",
        type=int,
        default=200,
        help=(
            "Minimum pseudo observations in a configuration x abundance bin "
            "for that bin-level pseudo estimate."
        ),
    )
    ap.add_argument(
        "--pseudo-min-subjects",
        type=int,
        default=3,
        help=(
            "Minimum pseudo subjects contributing to a configuration x "
            "abundance bin."
        ),
    )
    ap.add_argument(
        "--min-config-fraction",
        type=float,
        default=0.95,
        help=(
            "Fraction of pseudo configurations required for a native pseudo bin "
            "to define stable pseudo support."
        ),
    )
    ap.add_argument(
        "--min-grid-fraction",
        type=float,
        default=0.95,
        help=(
            "Minimum fraction of shared-grid points required for a pseudo "
            "configuration to enter integrated Var/MSD tests."
        ),
    )
    ap.add_argument(
        "--grid-points",
        type=int,
        default=101,
        help=(
            "Number of evaluation points on the shared absolute xstar support. "
            "This is an interpolation grid, not a new clonotype-level binning."
        ),
    )
    ap.add_argument(
        "--write-pseudo-sufficient-parquet",
        action="store_true",
        default=True,
        help="Write 02_pseudo_directed_pair_sufficient.parquet in addition to CSV.",
    )
    ap.add_argument(
        "--no-write-pseudo-sufficient-parquet",
        action="store_false",
        dest="write_pseudo_sufficient_parquet",
        help="Disable the optional parquet copy of pseudo sufficient statistics.",
    )
    return ap


def main() -> int:
    args = build_parser().parse_args()

    if not (0 < args.min_config_fraction <= 1):
        raise ValueError("--min-config-fraction must be in (0,1]")
    if not (0 < args.min_grid_fraction <= 1):
        raise ValueError("--min-grid-fraction must be in (0,1]")
    if args.grid_points < 21:
        raise ValueError("--grid-points must be >=21")

    longitudinal_dir = Path(args.longitudinal_dir)
    pseudo_rep_dir = Path(args.pseudo_representation_dir)
    outdir = ensure_dir(Path(args.outdir))

    print("[INFO] {} | {}".format(Path(__file__).name, SCRIPT_VERSION))
    print("[INFO] Primary: sample Var(dx_latent | xstar_latent), TT, longitudinal dt={}".format(args.longitudinal_dt))
    print("[INFO] Complementary: MSD; descriptive: MAD")
    print("[INFO] Percentile/rank abundance normalization: DISABLED")
    print("[INFO] Cross-dataset re-binning of clonotypes: DISABLED")

    longitudinal, longitudinal_path = load_longitudinal_curve(
        longitudinal_dir, args.longitudinal_dt
    )
    edges, edge_df, edge_path = load_pseudo_edges(pseudo_rep_dir)
    pseudo_reference, pseudo_reference_path = load_pseudo_fixed_reference(pseudo_rep_dir)

    # Pseudo directed-pair sufficient statistics.
    sufficient_source = None
    if args.pseudo_sufficient:
        sufficient_path = Path(args.pseudo_sufficient)
        sufficient = read_table(sufficient_path)
        times_by_subject = infer_times_from_sufficient(sufficient)
        sufficient_source = str(sufficient_path)
        print("[INFO] Reusing pseudo sufficient statistics: {}".format(sufficient_path))
    else:
        trajectories_path = Path(args.pseudo_trajectories)
        sufficient, times_by_subject = build_pseudo_directed_pair_sufficient(
            trajectories_path, edges, edge_df
        )
        sufficient_source = str(trajectories_path)

    # If a reused sufficient table came from the same edge table, merge edge metadata
    # only when those columns are absent.
    if "bin_index" not in sufficient.columns:
        raise ValueError("Pseudo sufficient table lacks bin_index")
    for col in ["bin_id", "bin_left", "bin_mid", "bin_right"]:
        if col not in sufficient.columns:
            sufficient = sufficient.merge(
                edge_df[["bin_index", col]], on="bin_index", how="left"
            )

    write_csv(sufficient, outdir / "02_pseudo_directed_pair_sufficient.csv")
    if args.write_pseudo_sufficient_parquet:
        write_parquet_best_effort(
            sufficient, outdir / "02_pseudo_directed_pair_sufficient.parquet"
        )

    index = PairSufficientIndex(sufficient, edge_df)
    configurations = generate_unique_orderings(
        times_by_subject,
        n_configurations=args.n_configurations,
        seed=args.seed,
    )
    config_curves, order_table = evaluate_pseudo_configurations(
        index,
        configurations,
        min_n=args.pseudo_min_n,
        min_subjects=args.pseudo_min_subjects,
    )
    write_csv(order_table, outdir / "03_pseudo_unique_orderings.csv")
    write_csv(config_curves, outdir / "04_pseudo_configuration_curves.csv")

    envelope = build_pseudo_native_envelope(
        config_curves,
        n_configurations=args.n_configurations,
        min_config_fraction=args.min_config_fraction,
    )
    write_csv(envelope, outdir / "05_pseudo_native_envelope.csv")

    support = determine_shared_support(longitudinal, envelope)
    print(
        "[INFO] longitudinal xstar support: [{:.4f}, {:.4f}]".format(
            support["longitudinal_min"], support["longitudinal_max"]
        )
    )
    print(
        "[INFO] pseudo stable xstar support: [{:.4f}, {:.4f}]".format(
            support["pseudo_stable_min"], support["pseudo_stable_max"]
        )
    )
    print(
        "[INFO] shared absolute xstar support: [{:.4f}, {:.4f}]".format(
            support["shared_min"], support["shared_max"]
        )
    )

    support_table = pd.DataFrame(
        [
            {
                "dataset": "longitudinal_dt{}".format(args.longitudinal_dt),
                "native_support_min": support["longitudinal_min"],
                "native_support_max": support["longitudinal_max"],
                "shared_support_min": support["shared_min"],
                "shared_support_max": support["shared_max"],
                "shared_support_width": support["shared_width"],
            },
            {
                "dataset": "pseudo_stable",
                "native_support_min": support["pseudo_stable_min"],
                "native_support_max": support["pseudo_stable_max"],
                "shared_support_min": support["shared_min"],
                "shared_support_max": support["shared_max"],
                "shared_support_width": support["shared_width"],
            },
        ]
    )
    write_csv(support_table, outdir / "06_shared_absolute_support_summary.csv")

    grid = np.linspace(
        support["shared_min"], support["shared_max"], int(args.grid_points)
    )
    long_interp = interpolate_longitudinal_to_grid(longitudinal, grid)
    config_metrics, pseudo_arrays = interpolate_pseudo_configs_to_grid(
        config_curves,
        grid,
        min_grid_fraction=args.min_grid_fraction,
    )
    write_csv(
        config_metrics,
        outdir / "08_pseudo_shared_grid_configuration_metrics.csv",
    )

    pseudo_mad = interpolate_pseudo_reference_mad(pseudo_reference, grid)
    shared_grid = build_shared_grid_table(
        grid,
        long_interp,
        pseudo_arrays,
        pseudo_mad,
    )
    write_csv(shared_grid, outdir / "07_shared_grid_comparison.csv")

    tests = build_empirical_tests(shared_grid, config_metrics)
    write_csv(tests, outdir / "09_empirical_tests.csv")

    # Native longitudinal curve is copied/normalized for plotting convenience.
    long_native_cols = [
        c for c in [
            "dataset_label", "representation", "conditioning", "mode", "dt",
            "bin_id", "bin_left", "bin_mid", "bin_right", "n", "n_subjects",
            "var_dx", "var_dx_ci_low", "var_dx_ci_high",
            "msd", "msd_ci_low", "msd_ci_high", "mad_dx", "mean_dx",
        ] if c in longitudinal.columns
    ]
    write_csv(
        longitudinal[long_native_cols],
        outdir / "01_longitudinal_dt1_native_curve.csv",
    )

    if not pseudo_reference.empty:
        write_csv(
            pseudo_reference,
            outdir / "11_pseudo_fixed_reference_curve.csv",
        )

    summary: Dict[str, object] = {
        "longitudinal_dt": int(args.longitudinal_dt),
        "n_pseudo_subjects": int(len(times_by_subject)),
        "pseudo_times_by_subject": times_by_subject,
        "n_pseudo_configurations": int(args.n_configurations),
        "global_reversal_deduplicated": True,
        "shared_support_min": support["shared_min"],
        "shared_support_max": support["shared_max"],
        "shared_support_width": support["shared_width"],
        "n_shared_grid_points": int(args.grid_points),
        "mad_permutation_inference": False,
        "mad_role": "descriptive_fixed_pseudo_reference_only",
    }
    for _, row in tests.iterrows():
        metric = str(row["metric"])
        summary["{}_longitudinal_mean_shared".format(metric)] = row[
            "longitudinal_mean_over_shared_support"
        ]
        summary["{}_pseudo_median_shared".format(metric)] = row["pseudo_median"]
        summary["{}_empirical_p_upper".format(metric)] = row["empirical_p_upper"]
        summary["{}_empirical_p_two_sided".format(metric)] = row[
            "empirical_p_two_sided_about_median"
        ]
        summary["{}_fraction_grid_outside_pseudo95".format(metric)] = row[
            "fraction_longitudinal_outside_pseudo_95"
        ]
    write_csv(pd.DataFrame([summary]), outdir / "10_comparison_summary.csv")

    manifest_rows = [
        {"role": "longitudinal_step11", "path": str(longitudinal_path)},
        {"role": "pseudo_xstar_bins_step8", "path": str(edge_path)},
        {"role": "pseudo_pair_source", "path": sufficient_source},
    ]
    if pseudo_reference_path is not None:
        manifest_rows.append(
            {"role": "pseudo_fixed_reference_step8", "path": str(pseudo_reference_path)}
        )
    write_csv(pd.DataFrame(manifest_rows), outdir / "00_input_manifest.csv")

    config_json = {
        "step": 12,
        "script": "12-longitudinal_vs_pseudo_fluctuation_baseline.py",
        "script_version": SCRIPT_VERSION,
        "longitudinal_source_step": 11,
        "pseudo_representation_source_step": 8,
        "pseudo_trajectory_source_step": 4,
        "primary_estimand": "sample Var(dx_latent | xstar_latent), TT",
        "complementary_metric": "MSD = E(dx_latent^2 | xstar_latent)",
        "descriptive_metric": "MAD(dx_latent | xstar_latent)",
        "longitudinal_dir": str(longitudinal_dir),
        "pseudo_representation_dir": str(pseudo_rep_dir),
        "pseudo_trajectories": args.pseudo_trajectories,
        "pseudo_sufficient": args.pseudo_sufficient,
        "longitudinal_dt": int(args.longitudinal_dt),
        "n_configurations": int(args.n_configurations),
        "seed": int(args.seed),
        "pseudo_min_n": int(args.pseudo_min_n),
        "pseudo_min_subjects": int(args.pseudo_min_subjects),
        "min_config_fraction": float(args.min_config_fraction),
        "min_grid_fraction": float(args.min_grid_fraction),
        "grid_points": int(args.grid_points),
        "pseudo_times_by_subject": times_by_subject,
        "global_reverse_equivalence": (
            "one representative retained because reversing all subject orders flips dx globally "
            "while leaving xstar, MSD and sample variance unchanged"
        ),
        "cross_dataset_abundance_rule": "native absolute xstar bins; comparison only on shared absolute support",
        "percentile_or_rank_normalization": False,
        "cross_dataset_rebinning": False,
        "mad_inference": "descriptive fixed Step-8 pseudo reference; no configuration-level permutation P value",
        "temporal_scaling_rule": "Step 12 does not feed Step 13; do not subtract the pseudo baseline from Step-11 curves before Step-13 temporal scaling",
    }
    with open(outdir / "00_run_config.json", "w", encoding="utf-8") as fh:
        json.dump(json_safe(config_json), fh, indent=2)
    print("[WRITE] {}".format(outdir / "00_run_config.json"))

    report = {
        "step": 12,
        "script": "12-longitudinal_vs_pseudo_fluctuation_baseline.py",
        "script_version": SCRIPT_VERSION,
        "analysis": "longitudinal_vs_pseudo_fluctuation_baseline",
        "summary": summary,
        "interpretive_boundary": (
            "Pseudo values estimate a technical fluctuation baseline on the shared absolute xstar support. "
            "Differences are descriptive excess-over-technical-baseline quantities and are not treated as pure biological variance."
        ),
    }
    with open(outdir / "12_fluctuation_baseline_report.json", "w", encoding="utf-8") as fh:
        json.dump(json_safe(report), fh, indent=2)
    print("[WRITE] {}".format(outdir / "12_fluctuation_baseline_report.json"))

    print("[RESULT] empirical shared-support tests")
    for _, row in tests.iterrows():
        print(
            "  {}: longitudinal={:.6g}, pseudo median={:.6g}, upper P={:.6g}, two-sided P={:.6g}".format(
                row["metric"],
                row["longitudinal_mean_over_shared_support"],
                row["pseudo_median"],
                row["empirical_p_upper"],
                row["empirical_p_two_sided_about_median"],
            )
        )
    print("[DONE] {}".format(outdir))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("[ERROR] {}".format(exc), file=sys.stderr)
        raise
