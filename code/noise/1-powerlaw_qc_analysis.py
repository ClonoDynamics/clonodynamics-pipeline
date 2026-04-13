#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
powerlaw_qc_analysis.py

STEP A — Analysis-only script for:
(A1) empirical heavy-tail validation and tail-model comparison
(A2) technical QC based on replicate concordance of the power-law exponent gamma

No figures are produced here.
Only analysis tables are written.

Conceptual workflow
-------------------
A1. Heavy-tail / tail-model validation
    - verifies that clone-size distributions are heavy-tailed
    - compares power-law vs truncated log-normal in the high-frequency regime
    - scans tail thresholds defined by *quantiles*, not by raw xmin values
    - writes tables supporting figure panels analogous to:
        A) heavy-tailed structure
        B) distribution of per-repertoire mean Δ log-likelihood
        C) stability of model preference across tail-threshold quantiles

A2. Replicate QC based on gamma concordance
    - fits a power law separately to each original replicate
    - computes pairwise discordance:
          delta_gamma = |gamma_rep1 - gamma_rep2|
    - defines PASS / FAIL using ONLY gamma discordance
    - no WARN class in this version

Important QC rule
-----------------
PASS / FAIL is based only on gamma outliers:
- FAIL if one replicate fit failed
- FAIL if |z_delta_gamma| > gamma_mad_cutoff
- PASS otherwise

Other quantities such as xmin, KS, depth ratio, tail size and tail fraction are
retained in the output tables as diagnostics, but do not enter the final QC rule.

Inputs
------
A directory of repertoire files named by default like:
    <subject>_<time>-1
    <subject>_<time>-2
with at least:
    aaSeqCDR3, readCount

Outputs
-------
<out_dir>/
  A1_repertoire_tail_fits.csv
  A1_repertoire_threshold_scan.csv
  A1_repertoire_model_compare.csv
  A1_pair_tail_preference_summary.csv

  A2_replicate_powerlaw_fits.csv
  A2_pair_gamma_qc.csv
  A2_gamma_qc_thresholds.csv
  A2_summary_overall.csv
  A2_summary_by_subject.csv

Example
-------
python stepA_powerlaw_qc_analysis.py \
  --data_dir ./longitudinali \
  --out_dir ./results_stepA \
  --file_sep '\\t'

If filenames end with .tsv:
python stepA_powerlaw_qc_analysis.py \
  --data_dir ./longitudinali \
  --out_dir ./results_stepA \
  --file_sep '\\t' \
  --pattern '^(?P<subject>\\d+)_(?P<time>\\d+)-(?P<replica>[12])\\.tsv$'
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


PATTERN_DEFAULT = r"^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])$"


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------
def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def robust_mad(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    med = np.median(x)
    return float(np.median(np.abs(x - med)))


def robust_z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = robust_mad(x)
    if not np.isfinite(mad) or mad < 1e-12:
        return np.full_like(x, np.nan, dtype=float)
    return (x - med) / mad


def depth_ratio(d1: int, d2: int) -> float:
    d1 = int(d1)
    d2 = int(d2)
    if d1 <= 0 or d2 <= 0:
        return np.nan
    return float(max(d1, d2) / min(d1, d2))


# ---------------------------------------------------------------------
# Power-law fit (continuous Clauset-style)
# ---------------------------------------------------------------------
@dataclass
class PowerLawFitResult:
    success: bool
    message: str
    n_total: int
    n_tail: int
    alpha_hat: float
    xmin_hat: float
    ks_distance: float
    tail_fraction: float


def validate_positive_1d(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x > 0]
    return np.sort(x)


def continuous_alpha_mle(x_tail: np.ndarray, xmin: float) -> float:
    x_tail = np.asarray(x_tail, dtype=float)
    if x_tail.size == 0:
        return np.nan
    denom = np.sum(np.log(x_tail / xmin))
    if denom <= 0:
        return np.nan
    return 1.0 + x_tail.size / denom


def powerlaw_cdf_continuous(x: np.ndarray, xmin: float, alpha: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if alpha <= 1.0 or xmin <= 0:
        return np.full_like(x, np.nan, dtype=float)
    return 1.0 - np.power(x / xmin, 1.0 - alpha)


def ks_distance_continuous(x_tail_sorted: np.ndarray, xmin: float, alpha: float) -> float:
    n = x_tail_sorted.size
    if n == 0:
        return np.nan
    F_emp = np.arange(1, n + 1, dtype=float) / n
    F_model = powerlaw_cdf_continuous(x_tail_sorted, xmin, alpha)
    return float(np.max(np.abs(F_emp - F_model)))


def estimate_xmin_ks_continuous(
    x: np.ndarray,
    min_tail_size: int = 50,
    candidate_step: int = 1,
) -> Tuple[float, float, float, int]:
    x = validate_positive_1d(x)
    if x.size < min_tail_size:
        return np.nan, np.nan, np.nan, 0

    unique_vals = np.unique(x)
    candidates = unique_vals[::max(1, candidate_step)]

    best_xmin = np.nan
    best_alpha = np.nan
    best_ks = np.inf
    best_n_tail = 0

    for xmin in candidates:
        x_tail = x[x >= xmin]
        n_tail = x_tail.size
        if n_tail < min_tail_size:
            continue

        alpha_hat = continuous_alpha_mle(x_tail, xmin)
        if not np.isfinite(alpha_hat) or alpha_hat <= 1.0:
            continue

        ks = ks_distance_continuous(x_tail, xmin, alpha_hat)
        if ks < best_ks:
            best_xmin = float(xmin)
            best_alpha = float(alpha_hat)
            best_ks = float(ks)
            best_n_tail = int(n_tail)

    return best_xmin, best_alpha, best_ks, best_n_tail


def fit_powerlaw_on_frequencies(
    freqs: np.ndarray,
    min_tail_size: int = 50,
    candidate_step: int = 1,
) -> PowerLawFitResult:
    x = validate_positive_1d(freqs)
    if x.size < min_tail_size:
        return PowerLawFitResult(
            success=False,
            message=f"Not enough positive frequencies ({x.size})",
            n_total=int(x.size),
            n_tail=0,
            alpha_hat=np.nan,
            xmin_hat=np.nan,
            ks_distance=np.nan,
            tail_fraction=np.nan,
        )

    xmin_hat, alpha_hat, ks_hat, n_tail = estimate_xmin_ks_continuous(
        x=x,
        min_tail_size=min_tail_size,
        candidate_step=candidate_step,
    )

    if not np.isfinite(xmin_hat):
        return PowerLawFitResult(
            success=False,
            message="Failed KS-based xmin estimation",
            n_total=int(x.size),
            n_tail=0,
            alpha_hat=np.nan,
            xmin_hat=np.nan,
            ks_distance=np.nan,
            tail_fraction=np.nan,
        )

    return PowerLawFitResult(
        success=True,
        message="ok",
        n_total=int(x.size),
        n_tail=int(n_tail),
        alpha_hat=float(alpha_hat),
        xmin_hat=float(xmin_hat),
        ks_distance=float(ks_hat),
        tail_fraction=float(n_tail / max(len(x), 1)),
    )


# ---------------------------------------------------------------------
# Tail-model comparison: power-law vs truncated log-normal
# ---------------------------------------------------------------------
def _normal_logpdf(z: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    sigma = max(float(sigma), 1e-12)
    return -0.5 * np.log(2.0 * np.pi) - np.log(sigma) - 0.5 * ((z - mu) / sigma) ** 2


def powerlaw_tail_loglik(x_tail: np.ndarray, xmin: float, alpha: float) -> float:
    """
    Continuous Pareto density on x >= xmin:
        p(x) = (alpha - 1) * xmin^(alpha - 1) * x^(-alpha)
    """
    x_tail = np.asarray(x_tail, dtype=float)
    if x_tail.size == 0 or xmin <= 0 or alpha <= 1:
        return np.nan
    ll = (
        x_tail.size * np.log(alpha - 1.0)
        + (alpha - 1.0) * x_tail.size * np.log(xmin)
        - alpha * np.sum(np.log(x_tail))
    )
    return float(ll)


def truncated_lognormal_tail_loglik(x_tail: np.ndarray, xmin: float) -> Tuple[float, float, float]:
    """
    Truncated log-normal on x >= xmin.
    Uses plug-in estimates mu,sigma in log-space and includes truncation normalization.
    """
    x_tail = np.asarray(x_tail, dtype=float)
    if x_tail.size < 3 or xmin <= 0:
        return np.nan, np.nan, np.nan

    z = np.log(x_tail)
    mu = float(np.mean(z))
    sigma = float(np.std(z, ddof=1)) if z.size > 1 else np.nan
    if not np.isfinite(sigma) or sigma <= 1e-12:
        return np.nan, np.nan, np.nan

    # Standard normal CDF approximation via erf
    from math import erf, sqrt
    zmin = np.log(xmin)
    Phi = 0.5 * (1.0 + erf((zmin - mu) / (sigma * sqrt(2.0))))
    surv = max(1.0 - Phi, 1e-300)

    ll = np.sum(_normal_logpdf(z, mu, sigma) - np.log(x_tail)) - x_tail.size * np.log(surv)
    return float(ll), mu, sigma


def threshold_scan_tail_models(
    freqs: np.ndarray,
    subject: int,
    time: int,
    replica: int,
    threshold_quantiles: np.ndarray,
    min_tail_size: int,
) -> pd.DataFrame:
    """
    For each threshold quantile q:
      xmin = quantile(freqs, q)
      tail = freqs[freqs >= xmin]
    Then compares power-law vs truncated log-normal on the same tail.

    Output includes:
      - threshold_quantile
      - xmin_threshold
      - model log-likelihoods
      - mean delta log-likelihood per clonotype
    """
    x = validate_positive_1d(freqs)
    rows: List[Dict[str, float]] = []
    if x.size == 0:
        return pd.DataFrame()

    for q in threshold_quantiles:
        q = float(q)
        xmin = float(np.quantile(x, q))
        x_tail = x[x >= xmin]
        n_tail = int(x_tail.size)
        if n_tail < min_tail_size:
            continue

        alpha_hat = continuous_alpha_mle(x_tail, xmin)
        ll_pl = powerlaw_tail_loglik(x_tail, xmin, alpha_hat) if np.isfinite(alpha_hat) and alpha_hat > 1 else np.nan
        ks_pl = ks_distance_continuous(x_tail, xmin, alpha_hat) if np.isfinite(alpha_hat) and alpha_hat > 1 else np.nan

        ll_tln, mu_tln, sigma_tln = truncated_lognormal_tail_loglik(x_tail, xmin)

        delta_ll = ll_pl - ll_tln if np.isfinite(ll_pl) and np.isfinite(ll_tln) else np.nan
        mean_delta_ll = delta_ll / n_tail if np.isfinite(delta_ll) and n_tail > 0 else np.nan

        rows.append({
            "subject": subject,
            "time": time,
            "replica": replica,
            "pair_id": f"{subject}_{time}",
            "threshold_quantile": q,
            "xmin_threshold": xmin,
            "n_total": int(x.size),
            "n_tail": n_tail,
            "tail_fraction": float(n_tail / max(len(x), 1)),

            "pl_alpha_hat": float(alpha_hat) if np.isfinite(alpha_hat) else np.nan,
            "pl_gamma_hat": float(alpha_hat) if np.isfinite(alpha_hat) else np.nan,
            "pl_loglik": float(ll_pl) if np.isfinite(ll_pl) else np.nan,
            "pl_mean_loglik_per_clonotype": float(ll_pl / n_tail) if np.isfinite(ll_pl) and n_tail > 0 else np.nan,
            "pl_ks": float(ks_pl) if np.isfinite(ks_pl) else np.nan,

            "tln_loglik": float(ll_tln) if np.isfinite(ll_tln) else np.nan,
            "tln_mean_loglik_per_clonotype": float(ll_tln / n_tail) if np.isfinite(ll_tln) and n_tail > 0 else np.nan,
            "tln_mu": float(mu_tln) if np.isfinite(mu_tln) else np.nan,
            "tln_sigma": float(sigma_tln) if np.isfinite(sigma_tln) else np.nan,

            "delta_loglik_pl_minus_tln": float(delta_ll) if np.isfinite(delta_ll) else np.nan,
            "mean_delta_loglik_pl_minus_tln": float(mean_delta_ll) if np.isfinite(mean_delta_ll) else np.nan,
            "pl_better": bool(delta_ll > 0) if np.isfinite(delta_ll) else False,
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# File indexing and per-replicate processing
# ---------------------------------------------------------------------
def read_rep_file(fp: Path, sep: str) -> pd.DataFrame:
    df = pd.read_csv(fp, sep=sep)
    needed = {"aaSeqCDR3", "readCount"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"{fp.name}: missing columns {missing}. Expected at least {sorted(needed)}")

    out = df[["aaSeqCDR3", "readCount"]].copy()
    out["readCount"] = pd.to_numeric(out["readCount"], errors="coerce")
    out = out[np.isfinite(out["readCount"])].copy()
    out = out[out["readCount"] > 0].copy()
    return out


def index_pairs(data_dir: Path, pattern_str: str) -> List[Tuple[int, int, Path, Path]]:
    pattern = re.compile(pattern_str)
    idx: Dict[Tuple[int, int, int], Path] = {}

    for fp in data_dir.iterdir():
        if not fp.is_file():
            continue
        m = pattern.match(fp.name)
        if not m:
            continue
        s = int(m.group("subject"))
        t = int(m.group("time"))
        r = int(m.group("replica"))
        idx[(s, t, r)] = fp

    pairs: List[Tuple[int, int, Path, Path]] = []
    seen = sorted({(s, t) for (s, t, _) in idx.keys()})
    for s, t in seen:
        fp1 = idx.get((s, t, 1))
        fp2 = idx.get((s, t, 2))
        if fp1 is not None and fp2 is not None:
            pairs.append((s, t, fp1, fp2))
    return pairs


def fit_one_replicate(
    fp: Path,
    file_sep: str,
    min_tail_size: int,
    candidate_step: int,
) -> Tuple[PowerLawFitResult, int, np.ndarray]:
    df = read_rep_file(fp, sep=file_sep)
    depth = int(df["readCount"].sum())
    if depth <= 0:
        return (
            PowerLawFitResult(False, "Zero depth", 0, 0, np.nan, np.nan, np.nan, np.nan),
            0,
            np.array([], dtype=float),
        )

    freqs = (df["readCount"] / depth).to_numpy(dtype=float)
    fit = fit_powerlaw_on_frequencies(
        freqs=freqs,
        min_tail_size=min_tail_size,
        candidate_step=candidate_step,
    )
    return fit, depth, freqs


# ---------------------------------------------------------------------
# QC logic: PASS/FAIL only, based only on gamma outliers
# ---------------------------------------------------------------------
def assign_gamma_only_qc(pair_df: pd.DataFrame, gamma_mad_cutoff: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = pair_df.copy()
    d["QC_traffic"] = "PASS"
    d["QC_reason"] = "pass"

    d["z_delta_gamma"] = robust_z(d["delta_gamma"].to_numpy(dtype=float))
    d["gamma_outlier"] = pd.to_numeric(d["z_delta_gamma"], errors="coerce").abs() > float(gamma_mad_cutoff)

    fit_fail = (~d["fit_success_rep1"].astype(bool)) | (~d["fit_success_rep2"].astype(bool))
    fail_mask = fit_fail | d["gamma_outlier"].fillna(False)

    d.loc[fail_mask, "QC_traffic"] = "FAIL"

    for idx in d.index[fail_mask]:
        reasons = []
        if bool(fit_fail.loc[idx]):
            reasons.append("fit_failed")
        if bool(d.loc[idx, "gamma_outlier"]) if pd.notna(d.loc[idx, "gamma_outlier"]) else False:
            reasons.append("delta_gamma")
        d.at[idx, "QC_reason"] = "fail:" + ",".join(reasons) if reasons else "fail"

    d["QC_pass_main"] = d["QC_traffic"] == "PASS"
    d["QC_pass"] = d["QC_pass_main"]

    thresholds = pd.DataFrame([{
        "gamma_mad_cutoff": float(gamma_mad_cutoff),
        "pass_fail_rule": "FAIL if one replicate fit failed or |z_delta_gamma| > gamma_mad_cutoff; PASS otherwise",
    }])

    return d, thresholds


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="STEP A analysis-only script: heavy-tail validation + gamma-only QC.")
    ap.add_argument("--data_dir", type=Path, required=True, help="Directory with replicate files.")
    ap.add_argument("--out_dir", type=Path, required=True, help="Directory for analysis outputs.")
    ap.add_argument("--file_sep", type=str, default="\t", help="Input file separator. Default: tab.")
    ap.add_argument("--pattern", type=str, default=PATTERN_DEFAULT, help="Filename regex pattern.")

    # replicate-level power-law fit
    ap.add_argument("--min_tail_size", type=int, default=50, help="Minimum tail size for power-law fitting.")
    ap.add_argument("--candidate_step", type=int, default=1, help="Step for xmin candidate thinning.")

    # threshold scan
    ap.add_argument(
        "--threshold_quantiles",
        type=str,
        default="0.70,0.75,0.80,0.85,0.90,0.95",
        help="Comma-separated threshold quantiles for tail-model comparison.",
    )
    ap.add_argument("--scan_min_tail_size", type=int, default=50, help="Minimum tail size for threshold-scan comparisons.")

    # gamma-only QC
    ap.add_argument("--gamma_mad_cutoff", type=float, default=3.5, help="Robust z cutoff on delta_gamma for FAIL.")
    return ap


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    args = build_argparser().parse_args()

    if not args.data_dir.exists():
        raise FileNotFoundError(f"data_dir not found: {args.data_dir}")

    out_dir = ensure_dir(args.out_dir)

    threshold_quantiles = np.array(
        [float(s.strip()) for s in str(args.threshold_quantiles).split(",") if s.strip()],
        dtype=float,
    )
    if threshold_quantiles.size == 0 or np.any((threshold_quantiles <= 0) | (threshold_quantiles >= 1)):
        raise ValueError("threshold_quantiles must contain values strictly between 0 and 1")

    pairs = index_pairs(args.data_dir, args.pattern)
    print(f"Found {len(pairs)} complete replicate pairs.")

    rep_rows: List[Dict[str, object]] = []
    pair_rows: List[Dict[str, object]] = []
    scan_frames: List[pd.DataFrame] = []

    for subject, time, fp1, fp2 in pairs:
        pair_id = f"{subject}_{time}"
        print(f"Processing {pair_id}: {fp1.name} | {fp2.name}")

        fit1, depth1, freqs1 = fit_one_replicate(
            fp=fp1,
            file_sep=args.file_sep,
            min_tail_size=args.min_tail_size,
            candidate_step=args.candidate_step,
        )
        fit2, depth2, freqs2 = fit_one_replicate(
            fp=fp2,
            file_sep=args.file_sep,
            min_tail_size=args.min_tail_size,
            candidate_step=args.candidate_step,
        )

        rep_rows.extend([
            {
                "subject": subject,
                "time": time,
                "replica": 1,
                "pair_id": pair_id,
                "file": fp1.name,
                "depth": depth1,
                "fit_success": fit1.success,
                "fit_message": fit1.message,
                "n_total": fit1.n_total,
                "n_tail": fit1.n_tail,
                "tail_fraction": fit1.tail_fraction,
                "gamma_hat": fit1.alpha_hat,
                "xmin_hat": fit1.xmin_hat,
                "ks_distance": fit1.ks_distance,
            },
            {
                "subject": subject,
                "time": time,
                "replica": 2,
                "pair_id": pair_id,
                "file": fp2.name,
                "depth": depth2,
                "fit_success": fit2.success,
                "fit_message": fit2.message,
                "n_total": fit2.n_total,
                "n_tail": fit2.n_tail,
                "tail_fraction": fit2.tail_fraction,
                "gamma_hat": fit2.alpha_hat,
                "xmin_hat": fit2.xmin_hat,
                "ks_distance": fit2.ks_distance,
            },
        ])

        gamma1 = fit1.alpha_hat if fit1.success else np.nan
        gamma2 = fit2.alpha_hat if fit2.success else np.nan
        xmin1 = fit1.xmin_hat if fit1.success else np.nan
        xmin2 = fit2.xmin_hat if fit2.success else np.nan

        pair_rows.append({
            "subject": subject,
            "time": time,
            "pair_id": pair_id,
            "file_rep1": fp1.name,
            "file_rep2": fp2.name,
            "fit_success_rep1": fit1.success,
            "fit_success_rep2": fit2.success,

            "depth_rep1": depth1,
            "depth_rep2": depth2,
            "depth_ratio": depth_ratio(depth1, depth2),

            "gamma_rep1": gamma1,
            "gamma_rep2": gamma2,
            "delta_gamma": abs(gamma1 - gamma2) if (np.isfinite(gamma1) and np.isfinite(gamma2)) else np.nan,
            "gamma_mean": np.nanmean([gamma1, gamma2]) if (np.isfinite(gamma1) or np.isfinite(gamma2)) else np.nan,

            "xmin_rep1": xmin1,
            "xmin_rep2": xmin2,
            "delta_xmin": abs(xmin1 - xmin2) if (np.isfinite(xmin1) and np.isfinite(xmin2)) else np.nan,

            "ks_rep1": fit1.ks_distance if fit1.success else np.nan,
            "ks_rep2": fit2.ks_distance if fit2.success else np.nan,

            "n_tail_rep1": fit1.n_tail if fit1.success else np.nan,
            "n_tail_rep2": fit2.n_tail if fit2.success else np.nan,
            "tail_fraction_rep1": fit1.tail_fraction if fit1.success else np.nan,
            "tail_fraction_rep2": fit2.tail_fraction if fit2.success else np.nan,
        })

        # A1 threshold scan per replicate
        for replica, freqs in [(1, freqs1), (2, freqs2)]:
            freqs_pos = validate_positive_1d(freqs)
            if freqs_pos.size < args.scan_min_tail_size:
                continue

            scan_df = threshold_scan_tail_models(
                freqs=freqs_pos,
                subject=subject,
                time=time,
                replica=replica,
                threshold_quantiles=threshold_quantiles,
                min_tail_size=int(args.scan_min_tail_size),
            )
            if not scan_df.empty:
                scan_frames.append(scan_df)

    rep_df = pd.DataFrame(rep_rows)
    pair_df = pd.DataFrame(pair_rows)

    # -------------------- A2 QC --------------------
    pair_qc, thresholds_qc = assign_gamma_only_qc(pair_df=pair_df, gamma_mad_cutoff=args.gamma_mad_cutoff)

    overall = pd.DataFrame([{
        "n_pairs": int(len(pair_qc)),
        "n_pass": int((pair_qc["QC_traffic"] == "PASS").sum()),
        "n_fail": int((pair_qc["QC_traffic"] == "FAIL").sum()),
        "median_delta_gamma": float(np.nanmedian(pair_qc["delta_gamma"].to_numpy(dtype=float))) if len(pair_qc) else np.nan,
        "median_delta_xmin": float(np.nanmedian(pair_qc["delta_xmin"].to_numpy(dtype=float))) if len(pair_qc) else np.nan,
        "median_depth_ratio": float(np.nanmedian(pair_qc["depth_ratio"].to_numpy(dtype=float))) if len(pair_qc) else np.nan,
        "median_gamma_mean": float(np.nanmedian(pair_qc["gamma_mean"].to_numpy(dtype=float))) if len(pair_qc) else np.nan,
    }])

    by_subject = (
        pair_qc.groupby("subject")
        .agg(
            n_pairs=("pair_id", "size"),
            n_pass=("QC_traffic", lambda x: int((x == "PASS").sum())),
            n_fail=("QC_traffic", lambda x: int((x == "FAIL").sum())),
            median_delta_gamma=("delta_gamma", "median"),
            median_delta_xmin=("delta_xmin", "median"),
            median_depth_ratio=("depth_ratio", "median"),
            median_gamma_mean=("gamma_mean", "median"),
        )
        .reset_index()
    )

    # -------------------- A1 heavy-tail outputs --------------------
    scan_df = pd.concat(scan_frames, ignore_index=True) if scan_frames else pd.DataFrame()

    rep_tail_fits = rep_df.rename(columns={
        "gamma_hat": "pl_gamma_hat",
        "xmin_hat": "pl_xmin_hat",
        "ks_distance": "pl_ks_distance",
        "n_tail": "pl_n_tail",
        "tail_fraction": "pl_tail_fraction",
        "fit_success": "pl_fit_success",
        "fit_message": "pl_fit_message",
    })

    if not scan_df.empty:
        rep_compare = (
            scan_df.groupby(["subject", "time", "replica", "pair_id"], as_index=False)
            .agg(
                n_thresholds=("threshold_quantile", "size"),
                mean_delta_loglik=("mean_delta_loglik_pl_minus_tln", "mean"),
                median_delta_loglik=("mean_delta_loglik_pl_minus_tln", "median"),
                frac_thresholds_pl_better=("pl_better", lambda x: float(pd.Series(x).mean())),
                median_tail_fraction=("tail_fraction", "median"),
            )
        )

        pair_tail_pref = (
            rep_compare.groupby(["subject", "time", "pair_id"], as_index=False)
            .agg(
                mean_median_delta_loglik=("median_delta_loglik", "mean"),
                mean_frac_thresholds_pl_better=("frac_thresholds_pl_better", "mean"),
                min_frac_thresholds_pl_better=("frac_thresholds_pl_better", "min"),
                n_replicates=("replica", "size"),
            )
        )
    else:
        rep_compare = pd.DataFrame(columns=[
            "subject", "time", "replica", "pair_id", "n_thresholds",
            "mean_delta_loglik", "median_delta_loglik",
            "frac_thresholds_pl_better", "median_tail_fraction"
        ])
        pair_tail_pref = pd.DataFrame(columns=[
            "subject", "time", "pair_id",
            "mean_median_delta_loglik",
            "mean_frac_thresholds_pl_better",
            "min_frac_thresholds_pl_better",
            "n_replicates"
        ])

    # -------------------- Save --------------------
    # A1
    rep_tail_fits.to_csv(out_dir / "A1_repertoire_tail_fits.csv", index=False)
    scan_df.to_csv(out_dir / "A1_repertoire_threshold_scan.csv", index=False)
    rep_compare.to_csv(out_dir / "A1_repertoire_model_compare.csv", index=False)
    pair_tail_pref.to_csv(out_dir / "A1_pair_tail_preference_summary.csv", index=False)

    # A2
    rep_df.to_csv(out_dir / "A2_replicate_powerlaw_fits.csv", index=False)
    pair_qc.to_csv(out_dir / "A2_pair_gamma_qc.csv", index=False)
    thresholds_qc.to_csv(out_dir / "A2_gamma_qc_thresholds.csv", index=False)
    overall.to_csv(out_dir / "A2_summary_overall.csv", index=False)
    by_subject.to_csv(out_dir / "A2_summary_by_subject.csv", index=False)

    print("\nSaved heavy-tail / tail-model outputs:")
    print(" -", out_dir / "A1_repertoire_tail_fits.csv")
    print(" -", out_dir / "A1_repertoire_threshold_scan.csv")
    print(" -", out_dir / "A1_repertoire_model_compare.csv")
    print(" -", out_dir / "A1_pair_tail_preference_summary.csv")

    print("\nSaved gamma-only QC outputs:")
    print(" -", out_dir / "A2_replicate_powerlaw_fits.csv")
    print(" -", out_dir / "A2_pair_gamma_qc.csv")
    print(" -", out_dir / "A2_gamma_qc_thresholds.csv")
    print(" -", out_dir / "A2_summary_overall.csv")
    print(" -", out_dir / "A2_summary_by_subject.csv")


if __name__ == "__main__":
    main()