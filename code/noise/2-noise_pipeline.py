#!/usr/bin/env python3
"""
noise_pipeline.py

Streamlined end-to-end pipeline for technical replicate noise modelling,
empirical denoising, and trajectory-ready observability outputs.

This version removes:
    - the internal power-law validation section
    - the internal QC computation section

It keeps:
    1) pairwise NB-priorPL noise-model fitting on technical replicates
    2) null-pool definition (clone / subject / global)
    3) empirical denoising via pooled logP distributions
    4) per-clonotype observability / trajectory-ready outputs
    5) optional alpha-sensitivity diagnostics
    6) aggregated null-parameter reporting

QC handling in this version
---------------------------
Pairs are filtered using an external QC table passed at launch with --qc-csv.
The QC table must contain at least:
    - subject
    - time
    - QC_pass
where QC_pass is interpreted as True/False.

Only replicate pairs with QC_pass == True are processed from the longitudinal
repertoire folder. Pairs failing QC are excluded from fitting and from all
subsequent downstream steps.

Required input
--------------
A directory of replicate files named by default like:
    <subject>_<time>-1
    <subject>_<time>-2
with at least the columns:
    aaSeqCDR3, readCount

Main outputs
------------
results_dir/
  null_params_by_pair_std.csv
  null_params_aggregated.csv
  logP_cutoffs_<null>.csv
  per_clone_denoised_<null>.csv
  trajectory_observability_<null>.csv
  trajectory_observability_counts_<null>.csv
  sensitivity_alpha_overall_<null>.csv                 [optional]
  sensitivity_alpha_by_subject_<null>.csv              [optional]
  sensitivity_alpha_pairwise_stability_<null>.csv      [optional]
  included_pairs_from_qc.csv
  per_clone_long_logP/
      null_per_clone_logP_std<subject>_<time>.csv

Requirements
------------
- Python 3.9+
- pandas, numpy
- an importable module noiseK_nb exposing fit_noiseK_nb_powerlaw(...)
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def add_prefix(name: str, prefix: str) -> str:
    if not prefix:
        return name
    pfx = str(prefix)
    if not (pfx.endswith("_") or pfx.endswith("-")):
        pfx += "_"
    return pfx + name


def empirical_cdf(sorted_arr: np.ndarray, x: float) -> float:
    if sorted_arr.size == 0:
        return float("nan")
    return float(np.searchsorted(sorted_arr, x, side="right") / sorted_arr.size)


# ---------------------------------------------------------------------
# Noise model fitting
# ---------------------------------------------------------------------
PATTERN_DEFAULT = r"^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])$"


def read_rep_file(fp: Path, sep: str = "\t") -> pd.DataFrame:
    df = pd.read_csv(fp, sep=sep)
    required = {"aaSeqCDR3", "readCount"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{fp.name}: missing columns: {missing}. Expected: {sorted(required)}")
    return df[["aaSeqCDR3", "readCount"]].copy()


def merge_two_reps(fp1: Path, fp2: Path, sep: str = "\t", key: str = "aaSeqCDR3"):
    d1 = read_rep_file(fp1, sep=sep).rename(columns={"readCount": "count_rep1"})
    d2 = read_rep_file(fp2, sep=sep).rename(columns={"readCount": "count_rep2"})

    merged = d1.merge(d2, on=key, how="outer")
    merged["count_rep1"] = merged["count_rep1"].fillna(0).astype(int)
    merged["count_rep2"] = merged["count_rep2"].fillna(0).astype(int)

    depth_rep1 = int(merged["count_rep1"].sum())
    depth_rep2 = int(merged["count_rep2"].sum())
    depths = [depth_rep1, depth_rep2]
    count_cols = ["count_rep1", "count_rep2"]
    return merged, depths, count_cols, depth_rep1, depth_rep2


def index_pairs(data_dir: Path, pattern: re.Pattern) -> List[Tuple[int, int, Path, Path]]:
    all_files = [fp for fp in data_dir.iterdir() if fp.is_file()]
    index: Dict[Tuple[int, int, int], Path] = {}
    for fp in all_files:
        m = pattern.match(fp.name)
        if not m:
            continue
        subject = int(m.group("subject"))
        time = int(m.group("time"))
        replica = int(m.group("replica"))
        index[(subject, time, replica)] = fp

    pairs: List[Tuple[int, int, Path, Path]] = []
    seen_st = sorted({(s, t) for (s, t, r) in index.keys()})
    for subject, time in seen_st:
        fp1 = index.get((subject, time, 1))
        fp2 = index.get((subject, time, 2))
        if fp1 and fp2:
            pairs.append((subject, time, fp1, fp2))
    return pairs


def parse_bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s.fillna(False)
    norm = s.astype(str).str.strip().str.lower()
    truthy = {"true", "1", "yes", "y", "pass"}
    return norm.isin(truthy)


def load_qc_pass_pairs(qc_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(qc_csv)
    required = {"subject", "time", "QC_pass"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"QC CSV missing required columns: {sorted(missing)}. "
            f"Expected at least: {sorted(required)}"
        )

    out = df[["subject", "time", "QC_pass"]].copy()
    out["subject"] = pd.to_numeric(out["subject"], errors="raise").astype(int)
    out["time"] = pd.to_numeric(out["time"], errors="raise").astype(int)
    out["QC_pass"] = parse_bool_series(out["QC_pass"])
    out = out.drop_duplicates(subset=["subject", "time"], keep="first")
    return out


def build_allowed_pair_set(qc_df: pd.DataFrame) -> Set[Tuple[int, int]]:
    keep = qc_df.loc[qc_df["QC_pass"] == True, ["subject", "time"]].copy()
    return set(map(tuple, keep.to_records(index=False)))


def run_noise_fit(
    data_dir: Path,
    results_dir: Path,
    perclone_dir: Path,
    file_sep: str,
    pattern_str: str,
    gamma_init: float,
    k_init: float,
    grid_size: int,
    allowed_pairs: Set[Tuple[int, int]],
    out_prefix: str = "",
) -> Path:
    try:
        from noiseK_nb import fit_noiseK_nb_powerlaw
    except Exception as e:
        raise ImportError(
            "Unable to import noiseK_nb.fit_noiseK_nb_powerlaw. "
            "Make sure noiseK_nb.py (or the package) is on the PYTHONPATH or in the same folder."
        ) from e

    ensure_dir(results_dir)
    ensure_dir(perclone_dir)

    pattern = re.compile(pattern_str)
    all_pairs = index_pairs(data_dir, pattern)
    selected_pairs = [(s, t, fp1, fp2) for (s, t, fp1, fp2) in all_pairs if (s, t) in allowed_pairs]

    print(f"Found {len(all_pairs)} complete pairs (replicate 1 & 2) in {data_dir}")
    print(f"Keeping {len(selected_pairs)} pairs with QC_pass == True from external QC table")

    rows = []
    for subject, time, fp1, fp2 in selected_pairs:
        pair_id = f"{subject}_{time}"
        print(f"\n--- Processing {pair_id}: {fp1.name} + {fp2.name}")
        try:
            merged, depths, count_cols, depth_rep1, depth_rep2 = merge_two_reps(fp1, fp2, sep=file_sep)

            merged["freq_rep1"] = merged["count_rep1"] / (depth_rep1 if depth_rep1 > 0 else 1)
            merged["freq_rep2"] = merged["count_rep2"] / (depth_rep2 if depth_rep2 > 0 else 1)

            f1 = merged["freq_rep1"].to_numpy(dtype=float)
            f2 = merged["freq_rep2"].to_numpy(dtype=float)
            c1 = merged["count_rep1"].to_numpy(dtype=int)
            c2 = merged["count_rep2"].to_numpy(dtype=int)
            freq_geo = np.zeros(len(merged), dtype=float)
            both = (c1 > 0) & (c2 > 0)
            only1 = (c1 > 0) & (c2 == 0)
            only2 = (c2 > 0) & (c1 == 0)
            freq_geo[both] = np.sqrt(f1[both] * f2[both])
            freq_geo[only1] = f1[only1]
            freq_geo[only2] = f2[only2]
            merged["freq_geo"] = freq_geo

            fit = fit_noiseK_nb_powerlaw(
                merged,
                depths,
                count_cols,
                gamma_init=gamma_init,
                k_init=k_init,
                grid_size=grid_size,
            )

            perclone_out = perclone_dir / f"{add_prefix('null_per_clone_logP_std', out_prefix)}{pair_id}.csv"
            per_clone = fit.per_clone_logP.copy()

            if "aaSeqCDR3" in per_clone.columns:
                key_col = "aaSeqCDR3"
            elif "cloneId" in per_clone.columns:
                key_col = "cloneId"
            elif "cloneID" in per_clone.columns:
                key_col = "cloneID"
            else:
                key_col = per_clone.columns[0]

            extra_cols = ["count_rep1", "count_rep2", "freq_rep1", "freq_rep2", "freq_geo"]
            extra = merged[[key_col] + [c for c in extra_cols if c in merged.columns]].copy() if key_col in merged.columns else None
            if extra is not None:
                per_clone = per_clone.merge(extra, on=key_col, how="left")

            per_clone["depth_rep1"] = depth_rep1
            per_clone["depth_rep2"] = depth_rep2
            per_clone["subject"] = subject
            per_clone["time"] = time
            per_clone["pair_id"] = pair_id
            per_clone.to_csv(perclone_out, index=False)

            n_clones_rep1 = int((merged["count_rep1"] > 0).sum())
            n_clones_rep2 = int((merged["count_rep2"] > 0).sum())
            inter = int(((merged["count_rep1"] > 0) & (merged["count_rep2"] > 0)).sum())
            union = int(((merged["count_rep1"] > 0) | (merged["count_rep2"] > 0)).sum())
            overlap_jaccard = float(inter / union) if union > 0 else float("nan")
            depth_ratio = float(max(depth_rep1, depth_rep2) / max(1, min(depth_rep1, depth_rep2)))

            row = {
                "subject": subject,
                "time": time,
                "file_rep1": fp1.name,
                "file_rep2": fp2.name,
                "depth_rep1": depth_rep1,
                "depth_rep2": depth_rep2,
                "depth_ratio": depth_ratio,
                "n_clones_rep1": n_clones_rep1,
                "n_clones_rep2": n_clones_rep2,
                "overlap_jaccard": overlap_jaccard,
                "success": bool(getattr(fit, "success", True)),
                "message": str(getattr(fit, "message", "")),
                "log_likelihood": float(getattr(fit, "log_likelihood", np.nan)),
                "QC_pass": True,
            }

            params = getattr(fit, "params", {})
            if isinstance(params, dict):
                for k, v in params.items():
                    row[k] = v

            rows.append(row)

        except Exception as e:
            eprint(f"ERROR in {pair_id}: {e}")
            rows.append({
                "subject": subject,
                "time": time,
                "file_rep1": fp1.name,
                "file_rep2": fp2.name,
                "depth_rep1": None,
                "depth_rep2": None,
                "success": False,
                "message": f"ERROR: {e}",
                "QC_pass": True,
            })

    out_table = results_dir / add_prefix("null_params_by_pair_std.csv", out_prefix)
    pd.DataFrame(rows).to_csv(out_table, index=False)
    print(f"\nSaved parameter table: {out_table}")
    print(f"Saved per-clone logP CSVs in: {perclone_dir}")
    return out_table


# ---------------------------------------------------------------------
# Null choice + denoising + aggregations
# ---------------------------------------------------------------------
def infer_clone_and_logp_cols(perclone: pd.DataFrame) -> Tuple[str, str]:
    clone_col_candidates = ["aaSeqCDR3", "cloneId", "cloneID", "clonotype", "cdr3"]
    logp_col_candidates = ["logP", "logp", "log_prob", "log_probability", "logP_noise", "logPnull"]

    clone_col = next((c for c in clone_col_candidates if c in perclone.columns), None)
    if clone_col is None:
        raise ValueError(f"Cannot find a clonotype column among {clone_col_candidates}. Columns: {list(perclone.columns)}")

    logp_col = next((c for c in logp_col_candidates if c in perclone.columns), None)
    if logp_col is None:
        raise ValueError(f"Cannot find a logP column among {logp_col_candidates}. Columns: {list(perclone.columns)}")

    return clone_col, logp_col


def load_perclone_with_metadata(perclone_dir: Path, qc_df: pd.DataFrame) -> pd.DataFrame:
    files = sorted(perclone_dir.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files in {perclone_dir}")

    qc_key = qc_df.set_index(["subject", "time"])["QC_pass"].to_dict()

    rows = []
    pair_pat = re.compile(r"(?P<sid>\d+_\d+)")
    for fp in files:
        m = pair_pat.search(fp.stem)
        if not m:
            continue
        pair_id = m.group("sid")
        subject, time = map(int, pair_id.split("_"))
        if not bool(qc_key.get((subject, time), False)):
            continue
        d = pd.read_csv(fp)
        d["subject"] = subject
        d["time"] = time
        d["pair_id"] = pair_id
        d["QC_pass_pair"] = True
        rows.append(d)

    if not rows:
        raise ValueError(f"Could not parse any QC-pass subject_time from filenames in {perclone_dir}")
    return pd.concat(rows, ignore_index=True)


def apply_pool_and_cutoff(
    perclone: pd.DataFrame,
    qc_df: pd.DataFrame,
    choice: str,
    alpha: float,
    min_qc_pairs_per_subject: int,
    logp_col: str,
    results_dir: Path,
    out_prefix: str = "",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    clean = perclone[perclone["QC_pass_pair"] == True].copy()
    if clean.empty:
        raise ValueError("No QC-pass pairs available to build the null.")

    if choice == "clone":
        clean["pool_id"] = clean["pair_id"]
        perclone["pool_id"] = perclone["pair_id"]
    elif choice == "subject":
        clean["pool_id"] = clean["subject"].astype(str)
        perclone["pool_id"] = perclone["subject"].astype(str)
    elif choice == "global":
        clean["pool_id"] = "GLOBAL"
        perclone["pool_id"] = "GLOBAL"
    else:
        raise ValueError("choice must be clone/subject/global")

    if choice == "subject":
        qc_pass_counts = (
            qc_df[qc_df["QC_pass"] == True]
            .groupby("subject")
            .size()
            .rename("n_qc_pairs")
            .to_dict()
        )
        perclone["n_qc_pairs_subject"] = perclone["subject"].map(qc_pass_counts).fillna(0).astype(int)
        clean["n_qc_pairs_subject"] = clean["subject"].map(qc_pass_counts).fillna(0).astype(int)

        perclone["pool_id"] = np.where(
            perclone["n_qc_pairs_subject"] >= min_qc_pairs_per_subject,
            perclone["subject"].astype(str),
            "GLOBAL",
        )
        clean["pool_id"] = np.where(
            clean["n_qc_pairs_subject"] >= min_qc_pairs_per_subject,
            clean["subject"].astype(str),
            "GLOBAL",
        )

    cutoffs = (
        clean.groupby("pool_id")[logp_col]
        .quantile(alpha)
        .rename("logP_cutoff")
        .reset_index()
    )
    cutoff_map = dict(zip(cutoffs["pool_id"], cutoffs["logP_cutoff"]))
    perclone["logP_cutoff"] = perclone["pool_id"].map(cutoff_map)

    if choice != "global":
        global_cut = float(clean[logp_col].quantile(alpha))
        perclone["logP_cutoff"] = perclone["logP_cutoff"].fillna(global_cut)

    perclone["observable_cutoff"] = perclone[logp_col] >= perclone["logP_cutoff"]

    ensure_dir(results_dir)
    cutoffs_out = results_dir / add_prefix(f"logP_cutoffs_{choice}.csv", out_prefix)
    cutoffs.to_csv(cutoffs_out, index=False)
    print("Saved cutoffs:", cutoffs_out)
    return perclone, cutoffs


def denoise_empirical(perclone: pd.DataFrame, logp_col: str) -> pd.DataFrame:
    clean = perclone[perclone["QC_pass_pair"] == True].copy()
    if clean.empty:
        raise ValueError("No QC-pass pairs: cannot compute empirical CDFs.")

    pool_arrays: Dict[str, np.ndarray] = {pid: np.sort(g[logp_col].to_numpy(dtype=float)) for pid, g in clean.groupby("pool_id")}
    if "GLOBAL" not in pool_arrays:
        pool_arrays["GLOBAL"] = np.sort(clean[logp_col].to_numpy(dtype=float))

    p_emp_low = []
    for pid, x in zip(perclone["pool_id"].astype(str).to_numpy(), perclone[logp_col].to_numpy(dtype=float)):
        arr = pool_arrays.get(pid, pool_arrays["GLOBAL"])
        p_emp_low.append(empirical_cdf(arr, float(x)))

    perclone["p_emp_low"] = p_emp_low
    return perclone


def aggregate_null_params(
    null_params_df: pd.DataFrame,
    qc_df: pd.DataFrame,
    choice: str,
    results_dir: Path,
    out_prefix: str = "",
) -> Optional[Path]:
    allowed = qc_df.loc[qc_df["QC_pass"] == True, ["subject", "time"]].copy()
    allowed["_key"] = allowed["subject"].astype(str) + "_" + allowed["time"].astype(str)

    params_clean = null_params_df.copy()
    params_clean["_key"] = params_clean["subject"].astype(str) + "_" + params_clean["time"].astype(str)
    params_clean = params_clean[params_clean["_key"].isin(set(allowed["_key"]))].copy()
    params_clean = params_clean[params_clean.get("success", False) == True].copy()

    if params_clean.empty:
        raise ValueError("No QC-pass successful pairs: cannot aggregate parameters.")

    candidates = ["gamma", "fmin", "k", "K", "N_total_mean", "N_total_min", "N_total_max", "log_likelihood"]
    cols = [c for c in candidates if c in params_clean.columns]
    if not cols:
        print("Note: no aggregatable parameter columns found. Skipping aggregation.")
        return None

    ensure_dir(results_dir)
    out = results_dir / add_prefix("null_params_aggregated.csv", out_prefix)

    if choice == "subject":
        agg = params_clean.groupby("subject")[cols].median(numeric_only=True).reset_index()
        agg.to_csv(out, index=False)
        print("Saved aggregated parameters (per subject):", out)
        return out

    med = params_clean[cols].median(numeric_only=True)
    agg = pd.DataFrame([med])
    agg.insert(0, "pool_id", "GLOBAL")
    agg.insert(1, "aggregation_level", choice)
    agg["n_qc_pairs"] = int(len(params_clean))
    agg.to_csv(out, index=False)
    if choice == "global":
        print("Saved aggregated parameters (global):", out)
    else:
        print("Saved aggregated parameters (global median; clone mode):", out)
    return out


def alpha_sensitivity_reports(
    perclone: pd.DataFrame,
    clone_col: str,
    alphas: List[float],
    results_dir: Path,
    null_choice: str,
    out_prefix: str = "",
) -> Tuple[Path, Path, Path]:
    if "p_value" not in perclone.columns:
        raise ValueError("perclone must contain 'p_value' (run denoising first).")

    ensure_dir(results_dir)

    overall_rows = []
    bysub_rows = []
    traj_counts_by_alpha = {}

    for a in alphas:
        obs = perclone["p_value"] < float(a)
        traj = (
            perclone.assign(_obs=obs)
            .groupby(["subject", clone_col, "time"], as_index=False)
            .agg(observable=("_obs", "max"))
        )
        traj_counts = (
            traj.groupby(["subject", clone_col], as_index=False)
            .agg(
                n_times_observable=("observable", "sum"),
                n_times_total=("time", "nunique"),
            )
        )
        traj_counts_by_alpha[a] = traj_counts

        overall_rows.append({
            "alpha": float(a),
            "n_rows": int(len(perclone)),
            "frac_observable_rows": float(obs.mean()),
            "n_unique_clonotypes": int(perclone[clone_col].nunique()),
        })

        tmp = traj_counts.groupby("subject").agg(
            n_clonotypes=("n_times_observable", "size"),
            frac_observable_clonotypes=("n_times_observable", lambda x: float((x > 0).mean())),
            median_n_times_observable=("n_times_observable", "median"),
            mean_n_times_observable=("n_times_observable", "mean"),
        ).reset_index()
        tmp["alpha"] = float(a)
        bysub_rows.append(tmp)

    overall = pd.DataFrame(overall_rows)
    bysub = pd.concat(bysub_rows, ignore_index=True)

    out_overall = results_dir / add_prefix(f"sensitivity_alpha_overall_{null_choice}.csv", out_prefix)
    out_bysub = results_dir / add_prefix(f"sensitivity_alpha_by_subject_{null_choice}.csv", out_prefix)
    overall.to_csv(out_overall, index=False)
    bysub.to_csv(out_bysub, index=False)

    al = sorted(traj_counts_by_alpha.keys())
    stab_rows = []
    for i in range(len(al)):
        for j in range(i + 1, len(al)):
            a1, a2 = al[i], al[j]
            d1 = traj_counts_by_alpha[a1].copy()
            d2 = traj_counts_by_alpha[a2].copy()
            d1["key"] = d1["subject"].astype(str) + "||" + d1[clone_col].astype(str)
            d2["key"] = d2["subject"].astype(str) + "||" + d2[clone_col].astype(str)
            m = d1[["key", "n_times_observable"]].merge(
                d2[["key", "n_times_observable"]],
                on="key",
                how="outer",
                suffixes=("_a1", "_a2"),
            ).fillna(0)
            x = m["n_times_observable_a1"].to_numpy(dtype=float)
            y = m["n_times_observable_a2"].to_numpy(dtype=float)
            corr = float(np.corrcoef(x, y)[0, 1]) if (np.std(x) > 0 and np.std(y) > 0) else float("nan")
            stab_rows.append({
                "alpha_1": float(a1),
                "alpha_2": float(a2),
                "n_keys": int(len(m)),
                "fraction_identical_n_times": float((x == y).mean()),
                "pearson_corr_n_times": corr,
            })

    stab = pd.DataFrame(stab_rows)
    out_stab = results_dir / add_prefix(f"sensitivity_alpha_pairwise_stability_{null_choice}.csv", out_prefix)
    stab.to_csv(out_stab, index=False)

    return out_overall, out_bysub, out_stab


def trajectory_outputs(perclone: pd.DataFrame, clone_col: str, choice: str, results_dir: Path, out_prefix: str = "") -> Tuple[Path, Path]:
    traj = (
        perclone.groupby(["subject", clone_col, "time"], as_index=False)
        .agg(observable=("observable", "max"), p_emp_low=("p_emp_low", "min"))
    )

    traj_counts = (
        traj.groupby(["subject", clone_col], as_index=False)
        .agg(
            n_times_observable=("observable", "sum"),
            n_times_total=("time", "nunique"),
        )
    )

    out1 = results_dir / add_prefix(f"trajectory_observability_{choice}.csv", out_prefix)
    out2 = results_dir / add_prefix(f"trajectory_observability_counts_{choice}.csv", out_prefix)
    traj.to_csv(out1, index=False)
    traj_counts.to_csv(out2, index=False)

    print("Saved trajectory-ready outputs:")
    print(" -", out1)
    print(" -", out2)
    return out1, out2


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Noise pipeline: external QC filter -> NB-priorPL fit -> null choice -> denoising -> trajectory-ready outputs"
    )
    p.add_argument("--data-dir", type=Path, required=True,
                   help="Folder with replicate files named like '1_1-1' and '1_1-2'.")
    p.add_argument("--qc-csv", type=Path, required=True,
                   help="Path to external QC table (e.g. A2_pair_gamma_qc.csv) containing subject, time, QC_pass.")
    p.add_argument("--results-dir", type=Path, default=Path("results"),
                   help="Main output folder (default: ./results)")
    p.add_argument("--out-prefix", type=str, default="",
                   help="Optional prefix for all output files.")
    p.add_argument("--file-sep", type=str, default="\t",
                   help="Input file separator (default: tab).")
    p.add_argument("--pattern", type=str, default=PATTERN_DEFAULT,
                   help=f"Regex for replicate filenames (default: {PATTERN_DEFAULT})")

    # Fit
    p.add_argument("--gamma-init", type=float, default=1.6)
    p.add_argument("--k-init", type=float, default=50.0)
    p.add_argument("--grid-size", type=int, default=500)

    # Null choice + alpha
    p.add_argument("--null", dest="null_choice", choices=["clone", "subject", "global"], required=True,
                   help="Null model choice: clone (=per pair), subject (=pool per subject), global (=global pool).")
    p.add_argument("--alpha", type=float, default=0.01,
                   help="Alpha threshold for observability (default 0.01).")
    p.add_argument("--min-qc-pairs-per-subject", type=int, default=3,
                   help="If null=subject, minimum number of QC-pass pairs required to use the subject pool.")
    p.add_argument(
        "--tail",
        choices=["low", "high"],
        default="low",
        help=(
            "Which tail defines rare events under the null. "
            "low -> p_value = p_emp_low; high -> p_value = 1 - p_emp_low. "
            "Final observable is always defined as (p_value < alpha)."
        ),
    )

    # Skip / rerun controls
    p.add_argument("--skip-fit", action="store_true",
                   help="Skip fitting and reuse existing null_params_by_pair_std.csv and per-clone CSVs.")
    p.add_argument("--skip-alpha-sensitivity", action="store_true",
                   help="Skip alpha-sensitivity CSV summaries.")

    return p


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    args = build_argparser().parse_args()

    data_dir: Path = args.data_dir
    qc_csv: Path = args.qc_csv
    results_dir: Path = args.results_dir
    out_prefix: str = str(args.out_prefix or "")
    perclone_dir: Path = results_dir / "per_clone_long_logP"

    if not data_dir.exists():
        raise FileNotFoundError(f"data-dir not found: {data_dir}")
    if not qc_csv.exists():
        raise FileNotFoundError(f"qc-csv not found: {qc_csv}")
    if not (0.0 < args.alpha < 1.0):
        raise ValueError("alpha must be between 0 and 1")

    ensure_dir(results_dir)
    ensure_dir(perclone_dir)

    qc_df = load_qc_pass_pairs(qc_csv)
    allowed_pairs = build_allowed_pair_set(qc_df)
    if len(allowed_pairs) == 0:
        raise ValueError("The QC CSV does not contain any pair with QC_pass == True")

    included_pairs_out = results_dir / add_prefix("included_pairs_from_qc.csv", out_prefix)
    qc_df.loc[qc_df["QC_pass"] == True].sort_values(["subject", "time"]).to_csv(included_pairs_out, index=False)
    print(f"Loaded QC table: {qc_csv}")
    print(f"QC-pass pairs: {len(allowed_pairs)}")
    print(f"Saved included pairs list: {included_pairs_out}")

    null_params_csv = results_dir / add_prefix("null_params_by_pair_std.csv", out_prefix)

    # 1) Fit only QC-pass pairs
    if args.skip_fit:
        if not null_params_csv.exists():
            raise FileNotFoundError(f"--skip-fit requested but missing: {null_params_csv}")
        print(f"[skip-fit] Using {null_params_csv}")
    else:
        run_noise_fit(
            data_dir=data_dir,
            results_dir=results_dir,
            perclone_dir=perclone_dir,
            file_sep=args.file_sep,
            pattern_str=args.pattern,
            gamma_init=args.gamma_init,
            k_init=args.k_init,
            grid_size=args.grid_size,
            allowed_pairs=allowed_pairs,
            out_prefix=out_prefix,
        )

    # 2) load per-clone + metadata, restricted to QC-pass pairs from external table
    perclone = load_perclone_with_metadata(
        perclone_dir=perclone_dir,
        qc_df=qc_df,
    )
    print("Total per-clone rows:", len(perclone))
    clone_col, logp_col = infer_clone_and_logp_cols(perclone)
    print("Using clone_col =", clone_col, "| logp_col =", logp_col)

    # 3) null pool + cutoff
    perclone, _cutoffs = apply_pool_and_cutoff(
        perclone=perclone,
        qc_df=qc_df,
        choice=args.null_choice,
        alpha=args.alpha,
        min_qc_pairs_per_subject=args.min_qc_pairs_per_subject,
        logp_col=logp_col,
        results_dir=results_dir,
        out_prefix=out_prefix,
    )

    # 4) empirical denoising consistent with the chosen null
    perclone = denoise_empirical(perclone=perclone, logp_col=logp_col)
    perclone["p_emp_low"] = pd.to_numeric(perclone["p_emp_low"], errors="raise")
    perclone["p_emp_high"] = 1.0 - perclone["p_emp_low"]
    perclone["p_value"] = perclone["p_emp_low"] if args.tail == "low" else perclone["p_emp_high"]
    perclone["observable"] = perclone["p_value"] < float(args.alpha)

    bad = perclone[(perclone["p_value"] >= float(args.alpha)) & (perclone["observable"] == True)]
    if len(bad) > 0:
        cols = [c for c in ["subject", "time", clone_col, logp_col, "p_emp_low", "p_emp_high", "p_value", "observable"] if c in bad.columns]
        raise RuntimeError(
            f"Incoherent observability detected: {len(bad)} rows with p_value>=alpha but observable=True.\n"
            f"Examples:\n{bad[cols].head(5).to_string(index=False)}"
        )

    den_out = results_dir / add_prefix(f"per_clone_denoised_{args.null_choice}.csv", out_prefix)
    perclone.to_csv(den_out, index=False)
    print("Saved per-clone denoising:", den_out)

    # 5) alpha sensitivity (optional)
    if not args.skip_alpha_sensitivity:
        sens_alphas = [1e-2, 5e-3, 1e-3]
        out_overall, out_bysub, out_stab = alpha_sensitivity_reports(
            perclone=perclone,
            clone_col=clone_col,
            alphas=sens_alphas,
            results_dir=results_dir,
            null_choice=args.null_choice,
            out_prefix=out_prefix,
        )
        print("Saved alpha-sensitivity diagnostics:")
        print(" -", out_overall)
        print(" -", out_bysub)
        print(" -", out_stab)
    else:
        print("[skip-alpha-sensitivity] No alpha-sensitivity CSVs were written.")

    # 6) aggregate null parameters (reporting/downstream)
    null_params_df = pd.read_csv(null_params_csv)
    aggregate_null_params(
        null_params_df=null_params_df,
        qc_df=qc_df,
        choice=args.null_choice,
        results_dir=results_dir,
        out_prefix=out_prefix,
    )

    # 7) trajectory-ready outputs
    trajectory_outputs(
        perclone=perclone,
        clone_col=clone_col,
        choice=args.null_choice,
        results_dir=results_dir,
        out_prefix=out_prefix,
    )

    print("\nDONE.")


if __name__ == "__main__":
    main()
