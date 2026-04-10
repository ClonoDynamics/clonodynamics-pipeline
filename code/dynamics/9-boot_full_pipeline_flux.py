#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
09_bootstrap_full_pipeline_flux.py

Full-pipeline bootstrap for Fokker–Planck flux/source diagnostics.

WHAT THIS SCRIPT DOES
---------------------
For each bootstrap replicate, the script:

  1) cluster-resamples transitions
  2) refits the drift function b(x) for each selected model
  3) infers diffusion D(x) either by:
       a) recomputing residuals from the bootstrapped transitions
          using the bootstrapped drift fit (recommended), or
       b) resampling a precomputed residual table
  4) cluster-resamples trajectories and estimates the empirical density p_emp(x)
     with optional observable handling: ignore / filter / weight
  5) computes the probability flux:
         J(x) = b(x) p(x) - d/dx [ D(x) p(x) ]
  6) computes the source proxy:
         S(x) = dJ/dx

SUPPORTED DRIFT MODELS
----------------------
  - OU_linear
  - hinge_plateau_smooth

WHY THIS SCRIPT IS USEFUL
-------------------------
This is a full-pipeline uncertainty propagation step.
It does not only bootstrap the final flux/source curves, but also propagates
uncertainty from:
  - drift refitting
  - diffusion estimation
  - empirical density estimation

This makes it more rigorous than a partial bootstrap that keeps some inferred
quantities fixed.

KEY DESIGN CHOICE
-----------------
Outputs are written into a dedicated run directory:

    <out_dir>/<run_name>/

with simple, stable filenames such as:
    J_summary.csv
    S_summary.csv
    global_summary.csv
    flux_J.png
    source_S.png
    report.txt

This avoids excessively long filenames and makes downstream plotting, figure
assembly, and manuscript work much easier.

EXAMPLE
-------
python 09_bootstrap_full_pipeline_flux.py \
  --transitions ./results/3-transitions/p_02/transitions_all.csv \
  --trajectories ./results/4-trajectories/p_02/trajectories_long.csv \
  --out_dir ./results/9-fp-bootstrap \
  --run_name p02_TT_hinge_filter \
  --dt 1 \
  --obs_class ALL \
  --models OU_linear,hinge_plateau_smooth \
  --weight_col w \
  --observable_mode filter \
  --bootstrap_reps 500 \
  --cluster_cols subject,aaSeqCDR3

OPTIONAL RESIDUAL TABLE
-----------------------
If you already have a precomputed residual table, you may pass:
  --residuals diffusion_residuals_long.csv

and choose:
  --diffusion_source precomputed_residuals

Otherwise, the recommended default is:
  --diffusion_source recompute_from_transitions
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# basic utilities
# ============================================================

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def gaussian_smooth(y: np.ndarray, sigma_pts: float) -> np.ndarray:
    if sigma_pts is None or sigma_pts <= 0:
        return np.asarray(y, dtype=float).copy()

    y = np.asarray(y, dtype=float)
    r = int(max(3, np.ceil(4 * sigma_pts)))
    xs = np.arange(-r, r + 1)
    w = np.exp(-0.5 * (xs / sigma_pts) ** 2)
    w /= w.sum()
    return np.convolve(y, w, mode="same")


def trapz_norm(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    z = np.trapezoid(y, x)
    if not np.isfinite(z) or z <= 0:
        raise ValueError("Normalization failed.")
    return y / z


def summarize_ci(arr: np.ndarray, ci: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    alpha = 0.5 * (1.0 - ci)
    mean = np.mean(arr, axis=0)
    lo = np.quantile(arr, alpha, axis=0)
    hi = np.quantile(arr, 1.0 - alpha, axis=0)
    return mean, lo, hi


def scalar_ci(arr: np.ndarray, ci: float) -> Tuple[float, float, float]:
    alpha = 0.5 * (1.0 - ci)
    return (
        float(np.mean(arr)),
        float(np.quantile(arr, alpha)),
        float(np.quantile(arr, 1.0 - alpha)),
    )


def robust_bulk_mask(
    xgrid: np.ndarray,
    p_emp: np.ndarray,
    trim_frac: float,
    min_density_frac: float,
) -> np.ndarray:
    g = len(xgrid)
    k = int(np.floor(trim_frac * g))
    mask = np.ones(g, dtype=bool)

    if k > 0:
        mask[:k] = False
        mask[-k:] = False

    if min_density_frac > 0:
        pmax = np.nanmax(p_emp)
        if np.isfinite(pmax) and pmax > 0:
            mask &= (p_emp >= min_density_frac * pmax)

    return mask


def parse_true_mask(series: pd.Series, true_values=("true", "1", "t", "yes", "y")) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.to_numpy(dtype=bool)

    s = series.astype(str).str.strip().str.lower()
    tv = [v.strip().lower() for v in true_values]
    return s.isin(tv).to_numpy(dtype=bool)


def combine_optional_weights(
    base: Optional[np.ndarray],
    extra: Optional[np.ndarray],
) -> Optional[np.ndarray]:
    if base is None and extra is None:
        return None
    if base is None:
        return np.asarray(extra, dtype=float)
    if extra is None:
        return np.asarray(base, dtype=float)
    return np.asarray(base, dtype=float) * np.asarray(extra, dtype=float)


def parse_obs_class_selector(obs_class_arg: str) -> Tuple[bool, List[str]]:
    """
    Parse the --obs_class selector.

    Returns
    -------
    use_all_obs_classes : bool
        True if no filtering on observability class should be applied.
    selected_obs_classes : list[str]
        Explicit classes to retain when use_all_obs_classes is False.

    Accepted examples
    -----------------
      --obs_class TT
      --obs_class TT,TF,FT,FF
      --obs_class ALL
      --obs_class *
      --obs_class all
    """
    raw = str(obs_class_arg).strip()
    if raw == "":
        return True, []

    parts = [p.strip() for p in raw.split(",") if p.strip()]
    if len(parts) == 0:
        return True, []

    parts_upper = [p.upper() for p in parts]
    if any(p in {"ALL", "*"} for p in parts_upper):
        return True, []

    return False, parts_upper


# ============================================================
# density estimation on grid
# ============================================================

def density_on_grid_from_x(
    x: np.ndarray,
    xgrid: np.ndarray,
    smooth_sigma_pts: float = 0.0,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    xgrid = np.asarray(xgrid, dtype=float)

    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        good = np.isfinite(x) & np.isfinite(weights) & (weights > 0)
        x = x[good]
        weights = weights[good]
    else:
        good = np.isfinite(x)
        x = x[good]
        weights = None

    if x.size == 0:
        return np.zeros_like(xgrid)

    edges = np.zeros(len(xgrid) + 1, dtype=float)
    edges[1:-1] = 0.5 * (xgrid[:-1] + xgrid[1:])
    edges[0] = xgrid[0] - (edges[1] - xgrid[0])
    edges[-1] = xgrid[-1] + (xgrid[-1] - edges[-2])

    counts, _ = np.histogram(x, bins=edges, weights=weights)
    dx = np.diff(edges)
    total = counts.sum()

    if total <= 0 or not np.isfinite(total):
        p = np.zeros_like(xgrid)
    else:
        p = counts / (total * dx)

    if smooth_sigma_pts and smooth_sigma_pts > 0:
        p = gaussian_smooth(p, smooth_sigma_pts)
        p = np.clip(p, 0, np.inf)
        p = trapz_norm(p, xgrid)

    return p


# ============================================================
# cluster bootstrap helpers
# ============================================================

def build_cluster_indices(df: pd.DataFrame, cluster_cols: List[str]) -> List[np.ndarray]:
    cols = [c for c in cluster_cols if c in df.columns]

    if len(cols) == 0:
        if "aaSeqCDR3" in df.columns:
            cols = ["aaSeqCDR3"]
        else:
            df = df.copy()
            df["_row_id"] = np.arange(len(df))
            cols = ["_row_id"]

    g = df.groupby(cols, sort=False).indices
    return [np.asarray(v, dtype=int) for v in g.values()]


def resample_df_by_clusters(
    df: pd.DataFrame,
    cluster_cols: List[str],
    rng: np.random.Generator,
) -> pd.DataFrame:
    if len(df) == 0:
        return df.copy()

    cluster_indices = build_cluster_indices(df, cluster_cols)
    n_clusters = len(cluster_indices)
    draw = rng.integers(0, n_clusters, size=n_clusters)
    idx = np.concatenate([cluster_indices[i] for i in draw])
    return df.iloc[idx].copy()


# ============================================================
# integrals for global summaries
# ============================================================

def integrate_abs(y: np.ndarray, x: np.ndarray, mask: np.ndarray) -> float:
    yy = np.where(mask, y, 0.0)
    return float(np.trapezoid(np.abs(yy), x))


def integrate_sq(y: np.ndarray, x: np.ndarray, mask: np.ndarray) -> float:
    yy = np.where(mask, y, 0.0)
    return float(np.trapezoid(yy * yy, x))


def integrate_signed(y: np.ndarray, x: np.ndarray, mask: np.ndarray) -> float:
    yy = np.where(mask, y, 0.0)
    return float(np.trapezoid(yy, x))


def bulk_length(x: np.ndarray, mask: np.ndarray) -> float:
    xx = x[mask]
    if len(xx) < 2:
        return float("nan")
    return float(xx.max() - xx.min())


# ============================================================
# plotting helper
# ============================================================

def save_curve_plot(
    x: np.ndarray,
    mean: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    mask: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
    out_png: Path,
    fig_w: float = 5.0,
    fig_h: float = 3.5,
    dpi: int = 300,
) -> None:
    plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    plt.axhline(0, linewidth=1)
    plt.plot(x[mask], mean[mask], linewidth=1.6)
    plt.fill_between(x[mask], lo[mask], hi[mask], alpha=0.2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.close()


# ============================================================
# drift fitting
# ============================================================

def binned_medians(
    x: np.ndarray,
    y: np.ndarray,
    nbins: int,
    weights: Optional[np.ndarray] = None,
    min_rows_per_bin: int = 30,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    good = np.isfinite(x) & np.isfinite(y)
    if weights is not None:
        good &= np.isfinite(weights) & (weights > 0)

    x = x[good]
    y = y[good]

    if len(x) < max(min_rows_per_bin * 2, 20):
        raise ValueError("Not enough rows for binned medians.")

    qs = np.linspace(0, 1, nbins + 1)
    edges = np.quantile(x, qs)
    edges = np.unique(edges)

    if len(edges) < 4:
        raise ValueError("Too few unique bin edges.")

    xm, ym, nm = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x <= hi if hi == edges[-1] else x < hi)
        if m.sum() < min_rows_per_bin:
            continue
        xm.append(np.median(x[m]))
        ym.append(np.median(y[m]))
        nm.append(int(m.sum()))

    if len(xm) < 4:
        raise ValueError("Too few populated bins.")

    return np.asarray(xm), np.asarray(ym), np.asarray(nm)


def fit_ou_linear(xm: np.ndarray, ym: np.ndarray, nm: np.ndarray) -> Dict[str, float]:
    coef = np.polyfit(xm, ym, deg=1, w=np.sqrt(nm))
    slope, intercept = coef[0], coef[1]
    tau = (-1.0 / slope) if np.isfinite(slope) and slope < 0 else np.nan
    x_star = (-intercept / slope) if np.isfinite(slope) and slope != 0 else np.nan
    return {
        "intercept": float(intercept),
        "slope": float(slope),
        "tau": float(tau),
        "x_star": float(x_star),
    }


def eval_ou_linear(xgrid: np.ndarray, pars: Dict[str, float]) -> np.ndarray:
    return pars["intercept"] + pars["slope"] * xgrid


def fit_hinge_plateau_smooth(xm: np.ndarray, ym: np.ndarray, nm: np.ndarray) -> Dict[str, float]:
    x_candidates = np.quantile(xm, np.linspace(0.2, 0.8, 31))
    best = None

    for xstar in x_candidates:
        z = np.minimum(xm, xstar)
        coef = np.polyfit(z, ym, deg=1, w=np.sqrt(nm))
        slope, intercept = coef[0], coef[1]
        yhat = intercept + slope * z
        sse = np.sum(nm * (ym - yhat) ** 2)

        if (best is None) or (sse < best["sse"]):
            best = {
                "intercept": float(intercept),
                "slope": float(slope),
                "xstar": float(xstar),
                "sse": float(sse),
            }

    if best is None:
        raise ValueError("Hinge fit failed.")

    tau = (-1.0 / best["slope"]) if np.isfinite(best["slope"]) and best["slope"] < 0 else np.nan
    best["tau"] = float(tau)
    return best


def eval_hinge_plateau_smooth(
    xgrid: np.ndarray,
    pars: Dict[str, float],
    smooth_sigma_pts: float = 1.0,
) -> np.ndarray:
    z = np.minimum(xgrid, pars["xstar"])
    y = pars["intercept"] + pars["slope"] * z
    if smooth_sigma_pts and smooth_sigma_pts > 0:
        y = gaussian_smooth(y, smooth_sigma_pts)
    return y


def fit_drift_model(
    df: pd.DataFrame,
    model: str,
    xcol: str,
    ycol: str,
    weight_col: Optional[str],
    nbins_x: int,
    min_rows_per_bin: int,
    hinge_smooth_sigma_pts: float,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    x = df[xcol].to_numpy(dtype=float)
    y = df[ycol].to_numpy(dtype=float)
    weights = df[weight_col].to_numpy(dtype=float) if (weight_col and weight_col in df.columns) else None

    xm, ym, nm = binned_medians(
        x, y,
        nbins=nbins_x,
        weights=weights,
        min_rows_per_bin=min_rows_per_bin
    )
    xgrid = np.linspace(np.nanmin(xm), np.nanmax(xm), 800)

    if model == "OU_linear":
        pars = fit_ou_linear(xm, ym, nm)
        b = eval_ou_linear(xgrid, pars)

    elif model == "hinge_plateau_smooth":
        pars = fit_hinge_plateau_smooth(xm, ym, nm)
        b = eval_hinge_plateau_smooth(xgrid, pars, smooth_sigma_pts=hinge_smooth_sigma_pts)

    else:
        raise ValueError(f"Unsupported model: {model}")

    return xgrid, b, pars


# ============================================================
# diffusion fitting
# ============================================================

def fit_diffusion_from_residuals(
    df: pd.DataFrame,
    xgrid: np.ndarray,
    xcol: str,
    resid_col: str,
    dt: float,
    weight_col: Optional[str],
    nbins_x: int,
    min_rows_per_bin: int,
    smooth_sigma_pts: float,
    estimator: str = "var_mad",
) -> Tuple[np.ndarray, pd.DataFrame]:
    x = df[xcol].to_numpy(dtype=float)
    r = df[resid_col].to_numpy(dtype=float)
    w = df[weight_col].to_numpy(dtype=float) if (weight_col and weight_col in df.columns) else None

    good = np.isfinite(x) & np.isfinite(r)
    if w is not None:
        good &= np.isfinite(w) & (w > 0)

    x = x[good]
    r = r[good]

    if len(x) < max(min_rows_per_bin * 2, 30):
        raise ValueError("Not enough residual rows for diffusion fit.")

    qs = np.linspace(0, 1, nbins_x + 1)
    edges = np.quantile(x, qs)
    edges = np.unique(edges)

    if len(edges) < 4:
        raise ValueError("Too few unique diffusion bin edges.")

    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x <= hi if hi == edges[-1] else x < hi)
        if m.sum() < min_rows_per_bin:
            continue

        xm = float(np.median(x[m]))
        rr = r[m]

        if estimator == "var_mad":
            mad = np.median(np.abs(rr - np.median(rr)))
            var_est = (1.4826 * mad) ** 2
        elif estimator == "mean_resid2":
            var_est = float(np.mean(rr ** 2))
        else:
            raise ValueError("Unsupported estimator.")

        rows.append({
            "x_mid": xm,
            "var_est": float(var_est),
            "n": int(m.sum()),
        })

    ddf = pd.DataFrame(rows).sort_values("x_mid")
    if len(ddf) < 4:
        raise ValueError("Too few diffusion bins after filtering.")

    d_est = ddf["var_est"].to_numpy(dtype=float) / (2.0 * dt)
    d_est = np.clip(d_est, 1e-12, np.inf)

    d_on = np.interp(
        xgrid,
        ddf["x_mid"].to_numpy(dtype=float),
        d_est,
        left=d_est[0],
        right=d_est[-1],
    )
    d_on = gaussian_smooth(d_on, smooth_sigma_pts)
    d_on = np.clip(d_on, 1e-12, np.inf)

    ddf["D_est"] = d_est
    return d_on, ddf


# ============================================================
# trajectory observability handling
# ============================================================

def prepare_trajectory_table(
    traj: pd.DataFrame,
    traj_x_col: str,
    weight_col: Optional[str],
    observable_col: str,
    observable_mode: str,
    observable_threshold: float,
    observable_true_values: Tuple[str, ...],
) -> Tuple[pd.DataFrame, Optional[np.ndarray]]:
    out = traj.copy()

    if traj_x_col not in out.columns:
        if traj_x_col == "log_freq" and "freq" in out.columns:
            out["log_freq"] = np.log(np.clip(out["freq"].to_numpy(dtype=float), 1e-300, np.inf))
        else:
            raise ValueError(f"Trajectories missing '{traj_x_col}' and no fallback available.")

    out = out[np.isfinite(out[traj_x_col].to_numpy(dtype=float))].copy()
    if len(out) == 0:
        raise ValueError("No finite trajectory x values available.")

    base_weights = None
    if weight_col and weight_col in out.columns:
        base_weights = out[weight_col].to_numpy(dtype=float)

    observable_weights = None

    if observable_mode != "ignore":
        if observable_col not in out.columns:
            raise ValueError(
                f"observable_mode={observable_mode} but column '{observable_col}' is missing from trajectories."
            )

        obs = out[observable_col]

        if observable_mode == "filter":
            if pd.api.types.is_bool_dtype(obs):
                mask = obs.to_numpy(dtype=bool)
            elif pd.api.types.is_numeric_dtype(obs):
                mask = obs.to_numpy(dtype=float) >= float(observable_threshold)
            else:
                mask = parse_true_mask(obs, true_values=observable_true_values)

            out = out.loc[mask].copy()
            if len(out) == 0:
                raise ValueError("No trajectory rows left after observable filtering.")

            if weight_col and weight_col in out.columns:
                base_weights = out[weight_col].to_numpy(dtype=float)
            else:
                base_weights = None

        elif observable_mode == "weight":
            if pd.api.types.is_bool_dtype(obs):
                observable_weights = obs.astype(float).to_numpy()
            elif pd.api.types.is_numeric_dtype(obs):
                observable_weights = np.clip(obs.to_numpy(dtype=float), 0.0, np.inf)
            else:
                observable_weights = parse_true_mask(obs, true_values=observable_true_values).astype(float)

            keep = np.isfinite(observable_weights) & (observable_weights > 0)
            out = out.loc[keep].copy()
            observable_weights = observable_weights[keep]

            if len(out) == 0:
                raise ValueError("No trajectory rows left after observable weighting.")

            if weight_col and weight_col in out.columns:
                base_weights = out[weight_col].to_numpy(dtype=float)
            else:
                base_weights = None

    combined_weights = combine_optional_weights(base_weights, observable_weights)
    return out, combined_weights


# ============================================================
# metadata and IO helpers
# ============================================================

def build_run_dir(out_dir: Path, run_name: str, dt: float, obs_class: str, model_label: str) -> Path:
    rn = run_name.strip() if run_name else f"dt{dt:g}_{obs_class}_{model_label}"
    return ensure_dir(out_dir / rn)


def write_metadata_json(out_json: Path, metadata: Dict) -> None:
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


# ============================================================
# main
# ============================================================

def main() -> None:
    ap = argparse.ArgumentParser(description="Full-pipeline bootstrap for FP flux/source diagnostics.")

    # inputs
    ap.add_argument("--transitions", required=True)
    ap.add_argument("--residuals", default=None, help="Optional precomputed residual table.")
    ap.add_argument("--trajectories", required=True)

    # outputs
    ap.add_argument("--out_dir", required=True, help="Parent directory for run outputs.")
    ap.add_argument("--run_name", default="", help="Short run label. If empty, a default is built automatically.")

    # filters
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument(
        "--obs_class",
        default="TT",
        help=(
            "Observability class selector for transitions/residuals. "
            "Examples: TT ; TT,TF ; TT,TF,FT,FF ; ALL. "
            "Use ALL (or *) to disable filtering on observability class."
        ),
    )
    ap.add_argument("--models", default="OU_linear,hinge_plateau_smooth")
    ap.add_argument("--weight_col", default=None)
    ap.add_argument("--cluster_cols", default="subject,aaSeqCDR3")

    # column names
    ap.add_argument("--x0_col", default="x0")
    ap.add_argument("--dx_col", default="dx")
    ap.add_argument("--resid_col", default="resid")
    ap.add_argument("--traj_x_col", default="log_freq")
    ap.add_argument("--model_col_resid", default="model")
    ap.add_argument("--dt_col", default="dt")
    ap.add_argument("--obs_class_col", default="obs_class")

    # diffusion source
    ap.add_argument(
        "--diffusion_source",
        default="recompute_from_transitions",
        choices=["recompute_from_transitions", "precomputed_residuals", "auto"],
        help=(
            "How to infer D(x): recompute residuals from bootstrapped transitions "
            "(recommended), use precomputed residual table, or auto."
        ),
    )

    # trajectory observability handling
    ap.add_argument("--observable_col", default="observable")
    ap.add_argument(
        "--observable_mode",
        default="ignore",
        choices=["ignore", "filter", "weight"],
        help="ignore=do nothing; filter=keep observable rows only; weight=weight p_emp by observable.",
    )
    ap.add_argument("--observable_threshold", type=float, default=0.5)
    ap.add_argument("--observable_true_values", default="true,1,t,yes,y")

    # fitting params
    ap.add_argument("--nbins_x_drift", type=int, default=32)
    ap.add_argument("--nbins_x_diff", type=int, default=32)
    ap.add_argument("--min_rows_per_bin", type=int, default=200)
    ap.add_argument("--hinge_smooth_sigma_pts", type=float, default=1.0)
    ap.add_argument("--smooth_p_sigma_pts", type=float, default=2.0)
    ap.add_argument("--smooth_D_sigma_pts", type=float, default=2.0)
    ap.add_argument("--smooth_J_sigma_pts", type=float, default=2.0)
    ap.add_argument("--diffusion_estimator", default="var_mad", choices=["var_mad", "mean_resid2"])

    # bootstrap
    ap.add_argument("--bootstrap_reps", type=int, default=200)
    ap.add_argument("--bootstrap_seed", type=int, default=123)
    ap.add_argument("--ci", type=float, default=0.95)

    # bulk
    ap.add_argument("--bulk_trim_frac", type=float, default=0.03)
    ap.add_argument("--bulk_min_density_frac", type=float, default=0.01)

    # plot
    ap.add_argument("--fig_w", type=float, default=5.0)
    ap.add_argument("--fig_h", type=float, default=3.5)
    ap.add_argument("--dpi", type=int, default=300)

    args = ap.parse_args()

    out_dir = ensure_dir(Path(args.out_dir))
    rng = np.random.default_rng(args.bootstrap_seed)
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    cluster_cols = [c.strip() for c in args.cluster_cols.split(",") if c.strip()]
    true_values = tuple(v.strip() for v in args.observable_true_values.split(",") if v.strip())
    use_all_obs_classes, selected_obs_classes = parse_obs_class_selector(args.obs_class)

    # --------------------------------------------------------
    # load inputs
    # --------------------------------------------------------
    trans = pd.read_csv(args.transitions)
    traj_raw = pd.read_csv(args.trajectories)
    resid = pd.read_csv(args.residuals) if args.residuals else None

    for req in [args.dt_col, args.obs_class_col, args.x0_col, args.dx_col]:
        if req not in trans.columns:
            raise ValueError(f"transitions: missing required column '{req}'")

    trans_mask = (
        np.isfinite(trans[args.dt_col].to_numpy(dtype=float))
        & (np.abs(trans[args.dt_col].to_numpy(dtype=float) - args.dt) <= 1e-12)
    )
    if not use_all_obs_classes:
        trans_mask &= trans[args.obs_class_col].astype(str).str.upper().isin(selected_obs_classes)

    trans = trans.loc[trans_mask].copy()

    if len(trans) == 0:
        raise ValueError("No transition rows left after filtering.")

    if resid is not None:
        for req in [args.dt_col, args.obs_class_col, args.x0_col, args.resid_col, args.model_col_resid]:
            if req not in resid.columns:
                raise ValueError(f"residuals: missing required column '{req}'")

        resid_mask = (
            np.isfinite(resid[args.dt_col].to_numpy(dtype=float))
            & (np.abs(resid[args.dt_col].to_numpy(dtype=float) - args.dt) <= 1e-12)
        )
        if not use_all_obs_classes:
            resid_mask &= resid[args.obs_class_col].astype(str).str.upper().isin(selected_obs_classes)

        resid = resid.loc[resid_mask].copy()

    if args.diffusion_source == "precomputed_residuals" and resid is None:
        raise ValueError("--diffusion_source precomputed_residuals requires --residuals.")

    traj, traj_base_weights = prepare_trajectory_table(
        traj=traj_raw,
        traj_x_col=args.traj_x_col,
        weight_col=args.weight_col,
        observable_col=args.observable_col,
        observable_mode=args.observable_mode,
        observable_threshold=args.observable_threshold,
        observable_true_values=true_values,
    )

    # --------------------------------------------------------
    # baseline grid from transitions
    # --------------------------------------------------------
    x_full = trans[args.x0_col].to_numpy(dtype=float)
    good_x = np.isfinite(x_full)

    if good_x.sum() < 20:
        raise ValueError("Too few finite x0 values in transitions.")

    xgrid_global = np.linspace(
        np.nanquantile(x_full[good_x], 0.01),
        np.nanquantile(x_full[good_x], 0.99),
        800
    )

    # --------------------------------------------------------
    # baseline p_emp for bulk mask
    # --------------------------------------------------------
    p_base = density_on_grid_from_x(
        traj[args.traj_x_col].to_numpy(dtype=float),
        xgrid_global,
        smooth_sigma_pts=args.smooth_p_sigma_pts,
        weights=traj_base_weights,
    )

    bulk_mask = robust_bulk_mask(
        xgrid_global,
        p_base,
        trim_frac=args.bulk_trim_frac,
        min_density_frac=args.bulk_min_density_frac,
    )

    if bulk_mask.sum() < 10:
        raise ValueError("Bulk mask too small.")

    l_bulk = bulk_length(xgrid_global, bulk_mask)
    if not np.isfinite(l_bulk) or l_bulk <= 0:
        raise ValueError("Invalid bulk length.")

    # --------------------------------------------------------
    # analyze each selected model
    # --------------------------------------------------------
    for model in models:
        use_precomputed_residuals = False
        residual_mode_effective = args.diffusion_source
        resid_model = None

        if args.diffusion_source == "precomputed_residuals":
            resid_model = resid.loc[resid[args.model_col_resid].astype(str) == str(model)].copy() if resid is not None else None
            if resid_model is None or len(resid_model) == 0:
                raise ValueError(f"No residual rows for model '{model}'.")
            use_precomputed_residuals = True

        elif args.diffusion_source == "auto":
            resid_model = resid.loc[resid[args.model_col_resid].astype(str) == str(model)].copy() if resid is not None else None
            if resid_model is not None and len(resid_model) > 0:
                use_precomputed_residuals = True
                residual_mode_effective = "precomputed_residuals"
            else:
                use_precomputed_residuals = False
                residual_mode_effective = "recompute_from_transitions"

        else:
            use_precomputed_residuals = False
            residual_mode_effective = "recompute_from_transitions"

        model_short = "ou" if model == "OU_linear" else ("hinge" if model == "hinge_plateau_smooth" else model)
        obs_class_label = "ALL" if use_all_obs_classes else "-".join(selected_obs_classes)
        run_dir = build_run_dir(
            out_dir=out_dir,
            run_name=(f"{args.run_name}_{model_short}" if args.run_name else ""),
            dt=args.dt,
            obs_class=obs_class_label,
            model_label=model_short,
        )

        # fixed output paths
        j_reps_path = run_dir / "J_boot_reps.npy"
        s_reps_path = run_dir / "S_boot_reps.npy"
        j_summary_path = run_dir / "J_summary.csv"
        s_summary_path = run_dir / "S_summary.csv"
        global_summary_path = run_dir / "global_summary.csv"
        drift_summary_path = run_dir / "drift_summary.csv"
        diffusion_summary_path = run_dir / "diffusion_summary.csv"
        drift_params_path = run_dir / "drift_params.csv"
        diffusion_bins_path = run_dir / "diffusion_bins.csv"
        flux_png_path = run_dir / "flux_J.png"
        source_png_path = run_dir / "source_S.png"
        report_path = run_dir / "report.txt"
        metadata_path = run_dir / "run_metadata.json"

        # bootstrap storage
        j_reps = np.zeros((args.bootstrap_reps, len(xgrid_global)), dtype=float)
        s_reps = np.zeros((args.bootstrap_reps, len(xgrid_global)), dtype=float)
        b_reps = np.zeros((args.bootstrap_reps, len(xgrid_global)), dtype=float)
        d_reps = np.zeros((args.bootstrap_reps, len(xgrid_global)), dtype=float)

        global_metrics = {
            "int_abs_J": np.zeros(args.bootstrap_reps),
            "int_sq_J": np.zeros(args.bootstrap_reps),
            "int_signed_J": np.zeros(args.bootstrap_reps),
            "mean_signed_J": np.zeros(args.bootstrap_reps),
            "int_abs_S": np.zeros(args.bootstrap_reps),
            "int_sq_S": np.zeros(args.bootstrap_reps),
            "p_emp_int": np.zeros(args.bootstrap_reps),
        }

        drift_pars_rows: List[Dict[str, float]] = []
        diffusion_bin_rows: List[pd.DataFrame] = []

        # ----------------------------------------------------
        # bootstrap loop
        # ----------------------------------------------------
        for r in range(args.bootstrap_reps):
            # 1) resample transitions
            trans_b = resample_df_by_clusters(trans, cluster_cols, rng)

            # 2) refit drift
            xgrid_fit, b_fit, pars = fit_drift_model(
                trans_b,
                model=model,
                xcol=args.x0_col,
                ycol=args.dx_col,
                weight_col=args.weight_col,
                nbins_x=args.nbins_x_drift,
                min_rows_per_bin=args.min_rows_per_bin,
                hinge_smooth_sigma_pts=args.hinge_smooth_sigma_pts,
            )

            b_on = np.interp(xgrid_global, xgrid_fit, b_fit, left=b_fit[0], right=b_fit[-1])
            b_reps[r, :] = b_on

            pars_row = {"rep": r, "model": model}
            pars_row.update(pars)
            drift_pars_rows.append(pars_row)

            # 3) diffusion
            if use_precomputed_residuals:
                resid_b = resample_df_by_clusters(resid_model, cluster_cols, rng).copy()
            else:
                resid_b = resample_df_by_clusters(trans_b, cluster_cols, rng).copy()
                b_rows = np.interp(
                    resid_b[args.x0_col].to_numpy(dtype=float),
                    xgrid_global,
                    b_on,
                    left=b_on[0],
                    right=b_on[-1],
                )
                resid_b[args.resid_col] = resid_b[args.dx_col].to_numpy(dtype=float) - b_rows

            d_on, ddf = fit_diffusion_from_residuals(
                resid_b,
                xgrid=xgrid_global,
                xcol=args.x0_col,
                resid_col=args.resid_col,
                dt=args.dt,
                weight_col=args.weight_col,
                nbins_x=args.nbins_x_diff,
                min_rows_per_bin=args.min_rows_per_bin,
                smooth_sigma_pts=args.smooth_D_sigma_pts,
                estimator=args.diffusion_estimator,
            )

            d_reps[r, :] = d_on
            dtmp = ddf.copy()
            dtmp["rep"] = r
            dtmp["model"] = model
            dtmp["diffusion_source"] = residual_mode_effective
            diffusion_bin_rows.append(dtmp)

            # 4) trajectories -> p_emp
            traj_b = resample_df_by_clusters(traj, cluster_cols, rng)

            traj_weights = None
            if args.weight_col and args.weight_col in traj_b.columns:
                traj_weights = traj_b[args.weight_col].to_numpy(dtype=float)

            if args.observable_mode == "weight":
                obs = traj_b[args.observable_col]
                if pd.api.types.is_bool_dtype(obs):
                    obs_w = obs.astype(float).to_numpy()
                elif pd.api.types.is_numeric_dtype(obs):
                    obs_w = np.clip(obs.to_numpy(dtype=float), 0.0, np.inf)
                else:
                    obs_w = parse_true_mask(obs, true_values=true_values).astype(float)

                traj_weights = combine_optional_weights(traj_weights, obs_w)

            p_emp = density_on_grid_from_x(
                traj_b[args.traj_x_col].to_numpy(dtype=float),
                xgrid_global,
                smooth_sigma_pts=args.smooth_p_sigma_pts,
                weights=traj_weights,
            )

            global_metrics["p_emp_int"][r] = float(np.trapezoid(p_emp, xgrid_global))

            # 5) J(x)
            dp = d_on * p_emp
            d_dp_dx = np.gradient(dp, xgrid_global)
            j = b_on * p_emp - d_dp_dx
            j = gaussian_smooth(j, args.smooth_J_sigma_pts)

            # 6) S(x)
            s = np.gradient(j, xgrid_global)
            s = gaussian_smooth(s, args.smooth_J_sigma_pts)

            j_reps[r, :] = j
            s_reps[r, :] = s

            global_metrics["int_abs_J"][r] = integrate_abs(j, xgrid_global, bulk_mask)
            global_metrics["int_sq_J"][r] = integrate_sq(j, xgrid_global, bulk_mask)
            global_metrics["int_signed_J"][r] = integrate_signed(j, xgrid_global, bulk_mask)
            global_metrics["mean_signed_J"][r] = global_metrics["int_signed_J"][r] / l_bulk
            global_metrics["int_abs_S"][r] = integrate_abs(s, xgrid_global, bulk_mask)
            global_metrics["int_sq_S"][r] = integrate_sq(s, xgrid_global, bulk_mask)

        # ----------------------------------------------------
        # summarize bootstrap outputs
        # ----------------------------------------------------
        j_mean, j_lo, j_hi = summarize_ci(j_reps, args.ci)
        s_mean, s_lo, s_hi = summarize_ci(s_reps, args.ci)
        b_mean, b_lo, b_hi = summarize_ci(b_reps, args.ci)
        d_mean, d_lo, d_hi = summarize_ci(d_reps, args.ci)

        j_df = pd.DataFrame({
            "x": xgrid_global[bulk_mask],
            "J_mean": j_mean[bulk_mask],
            "J_lo": j_lo[bulk_mask],
            "J_hi": j_hi[bulk_mask],
        })

        s_df = pd.DataFrame({
            "x": xgrid_global[bulk_mask],
            "S_mean": s_mean[bulk_mask],
            "S_lo": s_lo[bulk_mask],
            "S_hi": s_hi[bulk_mask],
        })

        b_df = pd.DataFrame({
            "x": xgrid_global,
            "b_mean": b_mean,
            "b_lo": b_lo,
            "b_hi": b_hi,
        })

        d_df = pd.DataFrame({
            "x": xgrid_global,
            "D_mean": d_mean,
            "D_lo": d_lo,
            "D_hi": d_hi,
        })

        glob_rows = []
        for key, vals in global_metrics.items():
            m, lo, hi = scalar_ci(vals, args.ci)
            glob_rows.append({
                "metric": key,
                "mean": m,
                "lo": lo,
                "hi": hi
            })
        glob_df = pd.DataFrame(glob_rows)

        # ----------------------------------------------------
        # save arrays and tables
        # ----------------------------------------------------
        np.save(j_reps_path, j_reps)
        np.save(s_reps_path, s_reps)

        j_df.to_csv(j_summary_path, index=False)
        s_df.to_csv(s_summary_path, index=False)
        glob_df.to_csv(global_summary_path, index=False)
        b_df.to_csv(drift_summary_path, index=False)
        d_df.to_csv(diffusion_summary_path, index=False)
        pd.DataFrame(drift_pars_rows).to_csv(drift_params_path, index=False)
        pd.concat(diffusion_bin_rows, ignore_index=True).to_csv(diffusion_bins_path, index=False)

        # ----------------------------------------------------
        # save plots
        # ----------------------------------------------------
        save_curve_plot(
            x=xgrid_global,
            mean=j_mean,
            lo=j_lo,
            hi=j_hi,
            mask=bulk_mask,
            xlabel="x = ln f",
            ylabel="J(x)",
            title=f"Flux J(x) | {model}",
            out_png=flux_png_path,
            fig_w=args.fig_w,
            fig_h=args.fig_h,
            dpi=args.dpi,
        )

        save_curve_plot(
            x=xgrid_global,
            mean=s_mean,
            lo=s_lo,
            hi=s_hi,
            mask=bulk_mask,
            xlabel="x = ln f",
            ylabel="S(x) = dJ/dx",
            title=f"Source S(x) | {model}",
            out_png=source_png_path,
            fig_w=args.fig_w,
            fig_h=args.fig_h,
            dpi=args.dpi,
        )

        # ----------------------------------------------------
        # report
        # ----------------------------------------------------
        signed_row = glob_df.loc[glob_df["metric"] == "int_signed_J"].iloc[0]
        includes0 = (signed_row["lo"] <= 0.0) and (signed_row["hi"] >= 0.0)

        with open(report_path, "w", encoding="utf-8") as f:
            f.write("FULL-PIPELINE BOOTSTRAP REPORT\n")
            f.write("================================\n\n")
            f.write(f"model: {model}\n")
            f.write(f"dt: {args.dt}\n")
            f.write(f"obs_class_requested: {args.obs_class}\n")
            f.write(f"obs_class_effective: {obs_class_label}\n")
            f.write(f"bootstrap_reps: {args.bootstrap_reps}\n")
            f.write(f"cluster_cols: {cluster_cols}\n")
            f.write(f"diffusion_source_requested: {args.diffusion_source}\n")
            f.write(f"diffusion_source_effective: {residual_mode_effective}\n")
            f.write(f"observable_col: {args.observable_col}\n")
            f.write(f"observable_mode: {args.observable_mode}\n")
            f.write(f"observable_threshold: {args.observable_threshold}\n")
            f.write(f"bulk_gridpoints: {int(bulk_mask.sum())}/{len(xgrid_global)}\n")
            f.write(f"bulk_x_range: [{xgrid_global[bulk_mask].min():.6g}, {xgrid_global[bulk_mask].max():.6g}]\n")
            f.write(f"int_signed_J_CI_includes_0: {includes0}\n\n")

            f.write("GLOBAL METRICS\n")
            f.write("------------------------------\n")
            for _, row in glob_df.iterrows():
                f.write(f"{row['metric']}: {row['mean']:.6g}  {row['lo']:.6g}  {row['hi']:.6g}\n")

        # ----------------------------------------------------
        # metadata json
        # ----------------------------------------------------
        metadata = {
            "inputs": {
                "transitions": str(Path(args.transitions).resolve()),
                "trajectories": str(Path(args.trajectories).resolve()),
                "residuals": str(Path(args.residuals).resolve()) if args.residuals else None,
            },
            "run": {
                "model": model,
                "dt": args.dt,
                "obs_class_requested": args.obs_class,
                "obs_class_effective": obs_class_label,
                "obs_class_filter_disabled": bool(use_all_obs_classes),
                "selected_obs_classes": selected_obs_classes,
                "run_dir": str(run_dir.resolve()),
            },
            "filters": {
                "weight_col": args.weight_col,
                "cluster_cols": cluster_cols,
                "observable_col": args.observable_col,
                "observable_mode": args.observable_mode,
                "observable_threshold": args.observable_threshold,
            },
            "diffusion": {
                "source_requested": args.diffusion_source,
                "source_effective": residual_mode_effective,
                "estimator": args.diffusion_estimator,
            },
            "bootstrap": {
                "bootstrap_reps": args.bootstrap_reps,
                "bootstrap_seed": args.bootstrap_seed,
                "ci": args.ci,
            },
            "grid": {
                "n_grid": int(len(xgrid_global)),
                "bulk_gridpoints": int(bulk_mask.sum()),
                "bulk_trim_frac": args.bulk_trim_frac,
                "bulk_min_density_frac": args.bulk_min_density_frac,
            },
            "fitting": {
                "nbins_x_drift": args.nbins_x_drift,
                "nbins_x_diff": args.nbins_x_diff,
                "min_rows_per_bin": args.min_rows_per_bin,
                "hinge_smooth_sigma_pts": args.hinge_smooth_sigma_pts,
                "smooth_p_sigma_pts": args.smooth_p_sigma_pts,
                "smooth_D_sigma_pts": args.smooth_D_sigma_pts,
                "smooth_J_sigma_pts": args.smooth_J_sigma_pts,
            },
            "outputs": {
                "J_boot_reps": j_reps_path.name,
                "S_boot_reps": s_reps_path.name,
                "J_summary": j_summary_path.name,
                "S_summary": s_summary_path.name,
                "global_summary": global_summary_path.name,
                "drift_summary": drift_summary_path.name,
                "diffusion_summary": diffusion_summary_path.name,
                "drift_params": drift_params_path.name,
                "diffusion_bins": diffusion_bins_path.name,
                "flux_png": flux_png_path.name,
                "source_png": source_png_path.name,
                "report": report_path.name,
            },
        }
        write_metadata_json(metadata_path, metadata)

        print(f"[DONE] model={model}")
        print(f"  run_dir: {run_dir}")
        print(f"  saved: {j_summary_path.name}")
        print(f"  saved: {s_summary_path.name}")
        print(f"  saved: {global_summary_path.name}")
        print(f"  saved: {drift_summary_path.name}")
        print(f"  saved: {diffusion_summary_path.name}")
        print(f"  saved: {report_path.name}")

    print("\nAll selected models completed.")


if __name__ == "__main__":
    main()