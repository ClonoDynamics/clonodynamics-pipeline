#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
7-fit_drift_models_analysis_and_residuals.py
ROBUST v3.0 - multi-model drift fitting + residual-based diffusion summaries

PURPOSE
-------
This is the ANALYSIS-ONLY stage of the finite-time drift / diffusion workflow.

The script operates on transition-level data at a fixed lag dt and performs
four conceptually distinct tasks:

1) FIT MULTIPLE DRIFT MODELS
   It fits one or more candidate drift functions to the conditional mean
   displacement dx as a function of initial log-frequency x0.

   Supported models:
   - OU_linear
   - tanh_saturating
   - rational_saturating
   - hinge_plateau_smooth

   These fits are intended to characterize the deterministic / mean component
   of finite-time clonotype dynamics.

2) COMPUTE MODEL-SPECIFIC RESIDUALS
   For each fitted model, it computes:
       resid_model  = dx - dx_hat_model
       resid2_model = resid_model^2

   These residuals represent the stochastic component left after subtraction
   of the fitted mean drift.

3) BUILD EMPIRICAL, BINNED RESIDUAL / DIFFUSION SUMMARIES
   For each fitted model, residuals are summarized in x0-bins, producing
   robust and non-robust measures of residual spread:
       - resid_med       : median residual
       - mad             : median absolute deviation around resid_med
       - sigma_mad       : 1.4826 * mad
       - var_mad         : sigma_mad^2
       - mean_resid2     : mean(resid^2)

   In addition, the script writes diffusion-scaled versions of these quantities:
       - D_var_mad       = var_mad / (2*dt)
       - D_mean_resid2   = mean_resid2 / (2*dt)
       - D_mad2          = mad^2 / (2*dt)

   IMPORTANT:
   These are still BINNED, EMPIRICAL summaries.
   They are NOT yet the final smooth function D(x) used in the Fokker–Planck step.

4) EXPORT ALL REQUIRED DOWNSTREAM DATASETS
   The outputs are meant to support:
   - model comparison and plotting
   - residual diagnostics
   - later construction of a smooth D(x) in the FP diagnostics script

WORKFLOW POSITION
-----------------
This script should be run BEFORE the FP diagnostics builder.

Recommended workflow:
   transitions CSV
      -> this script
         -> drift_function_grid.csv
         -> fitted_values.csv
         -> diffusion_residuals_long.csv
         -> diffusion_binned.csv
      -> FP diagnostics builder
         -> selects one model (typically hinge_plateau_smooth)
         -> selects one diffusion proxy (typically D_var_mad)
         -> interpolates / smooths / floors D(x)
         -> computes p_closed, J, S

WHY DIFFUSION IS NOT FINALIZED HERE
-----------------------------------
The role of this script is to produce a transparent empirical description of
residual variability after drift subtraction.

It deliberately does NOT:
   - smooth D(x)
   - impose a positivity floor for final FP usage
   - choose the final operational diffusion function

Those choices belong downstream, because they are part of the FP model-building
step, not of the raw residual-summary step.

MAIN OUTPUTS
------------
1) fit_summary.csv
   Parameter estimates and information criteria for all fitted drift models.

2) model_compare.csv
   Delta-AIC / delta-BIC relative to a baseline model.

3) fitted_values.csv
   Original transitions plus model predictions and residuals.

4) drift_function_grid.csv
   Common x-grid with fitted drift curves b_model(x).

5) binned_medians.csv
   Binned medians for:
      - dx
      - residuals of each model
   Intended mainly for plotting the drift fit and residual-bias diagnostics.

6) diffusion_residuals_long.csv
   Long-format residual dataset, one row per transition × model.
   Intended for downstream analyses or custom diagnostics.

7) diffusion_binned.csv
   Binned empirical residual-spread summaries for each model.
   This is the key bridge to the FP stage.

OUTPUT INTERPRETATION
---------------------
- drift_function_grid.csv provides b(x)
- diffusion_binned.csv provides empirical ingredients from which D(x) can be built
- the final D(x) used in FP must be selected and regularized downstream

ASSUMPTIONS
-----------
Required input columns:
   - dt
   - x0
   - dx

Optional / commonly used columns:
   - obs_class
   - subject
   - aaSeqCDR3
   - time0, time1
   - one weight column (e.g. w, w_clone, w_class, ...)

EXAMPLE
-------
python ./code/dynamics/7-fit_drift_models_analysis_and_residuals.py \
  --csv ./results/3-transitions/p_01/transitions_all.csv \
  --data_dir ./results/5-drift-multimodel/p_01_dt1_TT_w/data \
  --dt 1 \
  --filter_class TT \
  --weight_col w_clone \
  --weight_cap_q 0.999 \
  --grid_n 500 \
  --grid_qlo 0.02 \
  --grid_qhi 0.98 \
  --nbins_plot 20 \
  --min_per_bin 200 \
  --nbins_diff 20 \
  --min_per_bin_diff 200 \
  --bootstrap \ 
  --cluster_cols subject,aaSeqCDR3 \ 
  --seed 123 \
  --models OU_linear,tanh_saturating,rational_saturating,hinge_plateau_smooth

"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Tuple, Dict, List

import numpy as np
import pandas as pd
from scipy.optimize import least_squares


# ============================================================
# IO / utilities
# ============================================================

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def finite_mask(*arrs: np.ndarray) -> np.ndarray:
    m = None
    for a in arrs:
        a = np.asarray(a, float)
        if a.ndim == 1:
            mm = np.isfinite(a)
        elif a.ndim == 2:
            mm = np.isfinite(a).all(axis=1)
        else:
            raise ValueError(f"finite_mask supports only 1D/2D arrays, got ndim={a.ndim}")
        m = mm if m is None else (m & mm)
    return m if m is not None else np.array([], dtype=bool)


def standardize_safe(x: np.ndarray) -> Tuple[np.ndarray, float, float]:
    x = np.asarray(x, float)
    m = np.isfinite(x)
    if not m.any():
        return x * np.nan, np.nan, np.nan
    xx = x[m]
    mu = float(np.mean(xx))
    sd = float(np.std(xx, ddof=1)) if xx.size > 1 else 1.0
    if (not np.isfinite(sd)) or sd <= 0:
        sd = 1.0
    z = (x - mu) / sd
    return z, mu, sd


def aic_bic(sse: float, n: int, k: int) -> Tuple[float, float]:
    if (not np.isfinite(sse)) or sse <= 0 or n <= 0:
        return np.nan, np.nan
    aic = n * np.log(sse / n) + 2 * k
    bic = n * np.log(sse / n) + np.log(n) * k
    return float(aic), float(bic)


def quantile_ci(x: np.ndarray, ci: float) -> Tuple[float, float]:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan
    lo = (1 - ci) / 2
    hi = 1 - lo
    return float(np.quantile(x, lo)), float(np.quantile(x, hi))


def sanitize_weights_for_fit(w: np.ndarray, cap_q: float = 0.999) -> np.ndarray:
    ww = np.asarray(w, float)
    ww = np.where(np.isfinite(ww) & (ww > 0), ww, 0.0)
    if not np.any(ww > 0):
        return ww

    q = float(cap_q)
    q = min(max(q, 0.90), 0.999999)
    wcap = np.quantile(ww[ww > 0], q)
    if np.isfinite(wcap) and wcap > 0:
        ww = np.clip(ww, 0.0, wcap)

    mean_w = np.mean(ww[ww > 0])
    if np.isfinite(mean_w) and mean_w > 0:
        ww = ww / mean_w

    return ww


# ============================================================
# Binning helpers
# ============================================================

def binned_median(x0: np.ndarray, y: np.ndarray, nbins: int, min_per_bin: int) -> pd.DataFrame:
    x0 = np.asarray(x0, float)
    y = np.asarray(y, float)
    m = finite_mask(x0, y)
    x0, y = x0[m], y[m]

    cols = ["bin", "x0_lo", "x0_hi", "x0_mid", "y_med", "n"]
    if x0.size == 0:
        return pd.DataFrame(columns=cols)

    edges = np.quantile(x0, np.linspace(0, 1, nbins + 1))
    edges = np.unique(edges)
    if edges.size < 2:
        return pd.DataFrame(columns=cols)

    nb = edges.size - 1
    rows = []
    for i in range(nb):
        lo, hi = edges[i], edges[i + 1]
        sel = (x0 >= lo) & (x0 <= hi) if i == nb - 1 else (x0 >= lo) & (x0 < hi)
        n = int(sel.sum())
        if n < int(min_per_bin):
            continue
        rows.append({
            "bin": int(i),
            "x0_lo": float(lo),
            "x0_hi": float(hi),
            "x0_mid": float(np.median(x0[sel])),
            "y_med": float(np.median(y[sel])),
            "n": n,
        })
    return pd.DataFrame(rows)


def _weighted_mean(x: np.ndarray, w: Optional[np.ndarray]) -> float:
    x = np.asarray(x, float)
    if w is None:
        return float(np.mean(x)) if x.size else np.nan
    w = np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not m.any():
        return np.nan
    return float(np.sum(w[m] * x[m]) / np.sum(w[m]))


def _weighted_median(x: np.ndarray, w: Optional[np.ndarray]) -> float:
    x = np.asarray(x, float)
    if w is None:
        return float(np.median(x)) if x.size else np.nan
    w = np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not m.any():
        return np.nan
    x2 = x[m]
    w2 = w[m]
    order = np.argsort(x2)
    x2 = x2[order]
    w2 = w2[order]
    cw = np.cumsum(w2) / np.sum(w2)
    return float(x2[np.searchsorted(cw, 0.5)])


def binned_diffusion(
    x0: np.ndarray,
    resid: np.ndarray,
    nbins: int,
    min_per_bin: int,
    dt: float,
    w: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Bin by x0 quantiles and estimate conditional residual / diffusion summaries.

    Returned columns:
      - resid_med
      - mad
      - sigma_mad
      - var_mad
      - mean_resid2
      - D_var_mad
      - D_mean_resid2
      - D_mad2

    Here D_* are only bin-wise empirical summaries:
        D = Var(resid | x0-bin) / (2*dt)

    They are NOT smoothed or regularized.
    """
    x0 = np.asarray(x0, float)
    resid = np.asarray(resid, float)
    m = finite_mask(x0, resid)
    if w is not None:
        w = np.asarray(w, float)
        m = m & np.isfinite(w) & (w >= 0)

    x0 = x0[m]
    resid = resid[m]
    ww = None if w is None else w[m]

    cols = [
        "bin", "x0_lo", "x0_hi", "x0_mid", "n",
        "resid_med", "mad", "sigma_mad", "var_mad", "mean_resid2",
        "D_var_mad", "D_mean_resid2", "D_mad2"
    ]
    if x0.size == 0:
        return pd.DataFrame(columns=cols)

    edges = np.quantile(x0, np.linspace(0, 1, nbins + 1))
    edges = np.unique(edges)
    if edges.size < 2:
        return pd.DataFrame(columns=cols)

    nb = edges.size - 1
    rows = []
    for i in range(nb):
        lo, hi = edges[i], edges[i + 1]
        sel = (x0 >= lo) & (x0 <= hi) if i == nb - 1 else (x0 >= lo) & (x0 < hi)
        n = int(sel.sum())
        if n < int(min_per_bin):
            continue

        r = resid[sel]
        wsel = None if ww is None else ww[sel]

        r_med = _weighted_median(r, wsel)
        mad = _weighted_median(np.abs(r - r_med), wsel)
        sigma_mad = 1.4826 * mad if np.isfinite(mad) else np.nan
        var_mad = float(sigma_mad * sigma_mad) if np.isfinite(sigma_mad) else np.nan
        mean_resid2 = _weighted_mean(r * r, wsel)

        D_var_mad = var_mad / (2.0 * float(dt)) if np.isfinite(var_mad) else np.nan
        D_mean_resid2 = mean_resid2 / (2.0 * float(dt)) if np.isfinite(mean_resid2) else np.nan
        D_mad2 = (mad * mad) / (2.0 * float(dt)) if np.isfinite(mad) else np.nan

        rows.append({
            "bin": int(i),
            "x0_lo": float(lo),
            "x0_hi": float(hi),
            "x0_mid": float(np.median(x0[sel])),
            "n": n,
            "resid_med": float(r_med),
            "mad": float(mad),
            "sigma_mad": float(sigma_mad) if np.isfinite(sigma_mad) else np.nan,
            "var_mad": float(var_mad) if np.isfinite(var_mad) else np.nan,
            "mean_resid2": float(mean_resid2) if np.isfinite(mean_resid2) else np.nan,
            "D_var_mad": float(D_var_mad) if np.isfinite(D_var_mad) else np.nan,
            "D_mean_resid2": float(D_mean_resid2) if np.isfinite(D_mean_resid2) else np.nan,
            "D_mad2": float(D_mad2) if np.isfinite(D_mad2) else np.nan,
        })

    return pd.DataFrame(rows)


def cluster_bootstrap_samples(df: pd.DataFrame, cluster_cols: List[str], reps: int, seed: int) -> List[pd.DataFrame]:
    rng = np.random.default_rng(seed)
    clusters = df[cluster_cols].drop_duplicates().reset_index(drop=True)
    idx = np.arange(len(clusters))
    boots: List[pd.DataFrame] = []
    for _ in range(int(reps)):
        take = rng.choice(idx, size=len(idx), replace=True)
        sel = clusters.iloc[take]
        boots.append(df.merge(sel, on=cluster_cols, how="inner"))
    return boots


# ============================================================
# Drift model functions
# ============================================================

def ou_predict(x: np.ndarray, a: float, b: float) -> np.ndarray:
    return a + b * x


def tanh_predict(x: np.ndarray, A: float, x_star: float, Delta: float) -> np.ndarray:
    return A * np.tanh((x_star - x) / Delta)


def rational_predict(x: np.ndarray, A: float, x_star: float, Delta: float) -> np.ndarray:
    z = (x_star - x) / Delta
    return A * (z / np.sqrt(1.0 + z * z))


def sigmoid(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, float)
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


def hinge_plateau_predict(x: np.ndarray, a: float, b: float, A: float, x_c: float, Delta: float) -> np.ndarray:
    """
    Smooth blend between linear (a+b*x) and plateau (-A) across x_c with width Delta.
    For x << x_c, s≈1 -> linear; for x >> x_c, s≈0 -> plateau.
    """
    s = 1.0 - sigmoid((x - x_c) / Delta)
    return (a + b * x) * s + (-A) * (1.0 - s)


# ============================================================
# Fit helpers
# ============================================================

def fit_ou_linear_standardized(
    x0: np.ndarray,
    dx: np.ndarray,
    w: Optional[np.ndarray],
    weight_cap_q: float,
) -> Dict[str, float]:
    x0 = np.asarray(x0, float)
    dx = np.asarray(dx, float)

    m = finite_mask(x0, dx)
    if w is not None:
        w = np.asarray(w, float)
        m = m & np.isfinite(w) & (w >= 0)

    if not m.any():
        return {"ok": 0, "a": np.nan, "b": np.nan, "x_star": np.nan, "tau": np.nan, "sse": np.inf, "n_used": 0}

    x = x0[m]
    y = dx[m]
    ww = None if w is None else w[m]

    z, mu, sd = standardize_safe(x)
    mm = finite_mask(z, y)
    z, y = z[mm], y[mm]
    if ww is not None:
        ww = ww[mm]

    if z.size < 10:
        return {"ok": 0, "a": np.nan, "b": np.nan, "x_star": np.nan, "tau": np.nan, "sse": np.inf, "n_used": int(z.size)}

    X = np.column_stack([np.ones_like(z), z])

    try:
        if ww is None:
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            yhat = X @ beta
            resid = y - yhat
            sse = float(np.sum(resid**2))
        else:
            ww2 = sanitize_weights_for_fit(ww, cap_q=weight_cap_q)
            if not np.any(ww2 > 0):
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                yhat = X @ beta
                resid = y - yhat
                sse = float(np.sum(resid**2))
            else:
                W = np.sqrt(ww2)[:, None]
                Xw = X * W
                yw = y * W[:, 0]
                beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
                yhat = X @ beta
                resid = y - yhat
                sse = float(np.sum(ww2 * resid**2))

        a_z, b_z = float(beta[0]), float(beta[1])
        if not (np.isfinite(mu) and np.isfinite(sd) and sd > 0 and np.isfinite(a_z) and np.isfinite(b_z)):
            return {"ok": 0, "a": np.nan, "b": np.nan, "x_star": np.nan, "tau": np.nan, "sse": np.inf, "n_used": int(z.size)}

        b_raw = b_z / sd
        a_raw = a_z - b_z * (mu / sd)
        x_star = (-a_raw / b_raw) if np.isfinite(b_raw) and b_raw != 0 else np.nan
        tau = (-1.0 / b_raw) if np.isfinite(b_raw) and b_raw < 0 else np.nan

        if not (np.isfinite(sse) and sse > 0):
            return {"ok": 0, "a": np.nan, "b": np.nan, "x_star": np.nan, "tau": np.nan, "sse": np.inf, "n_used": int(z.size)}

        return {
            "ok": 1,
            "a": float(a_raw),
            "b": float(b_raw),
            "x_star": float(x_star),
            "tau": float(tau),
            "sse": float(sse),
            "n_used": int(z.size),
        }
    except Exception:
        return {"ok": 0, "a": np.nan, "b": np.nan, "x_star": np.nan, "tau": np.nan, "sse": np.inf, "n_used": int(z.size)}


def _prep_xyw(x0: np.ndarray, dx: np.ndarray, w: Optional[np.ndarray], weight_cap_q: float):
    x0 = np.asarray(x0, float)
    dx = np.asarray(dx, float)
    m = finite_mask(x0, dx)
    if w is not None:
        w = np.asarray(w, float)
        m = m & np.isfinite(w) & (w >= 0)
    x = x0[m]
    y = dx[m]
    if x.size == 0:
        return x, y, None
    if w is None:
        ww2 = np.ones_like(y)
    else:
        ww2 = sanitize_weights_for_fit(w[m], cap_q=weight_cap_q)
        if not np.any(ww2 > 0):
            ww2 = np.ones_like(y)
    return x, y, ww2


def fit_tanh_saturating(x0, dx, w, weight_cap_q) -> Dict[str, float]:
    x, y, ww2 = _prep_xyw(x0, dx, w, weight_cap_q)
    if x.size < 30:
        return {"ok": 0, "A": np.nan, "x_star": np.nan, "Delta": np.nan, "sse": np.inf, "nfev": np.nan, "status": np.nan, "message": ""}

    sqrtw = np.sqrt(ww2)
    A0 = float(np.nanpercentile(np.abs(y), 75))
    if (not np.isfinite(A0)) or A0 <= 0:
        A0 = float(np.std(y)) if np.std(y) > 0 else 0.1

    idx = np.argsort(np.abs(y))
    k0 = max(50, int(0.02 * len(y)))
    xstar0 = float(np.median(x[idx[:k0]])) if len(idx) >= k0 else float(np.median(x))

    Delta0 = float(np.std(x))
    if (not np.isfinite(Delta0)) or Delta0 <= 0:
        Delta0 = 1.0

    p0 = np.array([A0, xstar0, Delta0], float)
    lb = np.array([-np.inf, -np.inf, 1e-6], float)
    ub = np.array([ np.inf,  np.inf,  np.inf], float)

    def resid(p):
        A, x_star, Delta = p
        yhat = tanh_predict(x, A, x_star, Delta)
        return sqrtw * (y - yhat)

    try:
        res = least_squares(resid, x0=p0, bounds=(lb, ub), max_nfev=20000, loss="soft_l1", f_scale=1.0)
        if (not res.success) or (not np.all(np.isfinite(res.x))):
            return {"ok": 0, "A": np.nan, "x_star": np.nan, "Delta": np.nan, "sse": np.inf,
                    "nfev": float(res.nfev), "status": float(res.status), "message": str(res.message)}
        A_hat, x_star_hat, Delta_hat = map(float, res.x)
        if not (np.isfinite(Delta_hat) and Delta_hat > 0):
            return {"ok": 0, "A": np.nan, "x_star": np.nan, "Delta": np.nan, "sse": np.inf,
                    "nfev": float(res.nfev), "status": float(res.status), "message": str(res.message)}
        r = y - tanh_predict(x, A_hat, x_star_hat, Delta_hat)
        sse = float(np.sum(ww2 * r * r))
        slope_local = (-A_hat / Delta_hat)
        tau_local = (Delta_hat / A_hat) if np.isfinite(A_hat) and A_hat != 0 else np.nan
        return {
            "ok": 1, "A": A_hat, "x_star": x_star_hat, "Delta": Delta_hat,
            "sse": sse, "slope_local": slope_local, "tau_local": tau_local,
            "nfev": int(res.nfev), "status": int(res.status), "message": str(res.message),
            "A0": float(A0), "xstar0": float(xstar0), "Delta0": float(Delta0),
        }
    except Exception:
        return {"ok": 0, "A": np.nan, "x_star": np.nan, "Delta": np.nan, "sse": np.inf, "nfev": np.nan, "status": np.nan, "message": ""}


def fit_rational_saturating(x0, dx, w, weight_cap_q) -> Dict[str, float]:
    x, y, ww2 = _prep_xyw(x0, dx, w, weight_cap_q)
    if x.size < 30:
        return {"ok": 0, "A": np.nan, "x_star": np.nan, "Delta": np.nan, "sse": np.inf, "nfev": np.nan, "status": np.nan, "message": ""}

    sqrtw = np.sqrt(ww2)
    A0 = float(np.nanpercentile(np.abs(y), 75))
    if (not np.isfinite(A0)) or A0 <= 0:
        A0 = float(np.std(y)) if np.std(y) > 0 else 0.1

    idx = np.argsort(np.abs(y))
    k0 = max(50, int(0.02 * len(y)))
    xstar0 = float(np.median(x[idx[:k0]])) if len(idx) >= k0 else float(np.median(x))

    Delta0 = float(np.std(x))
    if (not np.isfinite(Delta0)) or Delta0 <= 0:
        Delta0 = 1.0

    p0 = np.array([A0, xstar0, Delta0], float)
    lb = np.array([-np.inf, -np.inf, 1e-6], float)
    ub = np.array([ np.inf,  np.inf,  np.inf], float)

    def resid(p):
        A, x_star, Delta = p
        yhat = rational_predict(x, A, x_star, Delta)
        return sqrtw * (y - yhat)

    try:
        res = least_squares(resid, x0=p0, bounds=(lb, ub), max_nfev=20000, loss="soft_l1", f_scale=1.0)
        if (not res.success) or (not np.all(np.isfinite(res.x))):
            return {"ok": 0, "A": np.nan, "x_star": np.nan, "Delta": np.nan, "sse": np.inf,
                    "nfev": float(res.nfev), "status": float(res.status), "message": str(res.message)}
        A_hat, x_star_hat, Delta_hat = map(float, res.x)
        if not (np.isfinite(Delta_hat) and Delta_hat > 0):
            return {"ok": 0, "A": np.nan, "x_star": np.nan, "Delta": np.nan, "sse": np.inf,
                    "nfev": float(res.nfev), "status": float(res.status), "message": str(res.message)}
        r = y - rational_predict(x, A_hat, x_star_hat, Delta_hat)
        sse = float(np.sum(ww2 * r * r))
        slope_local = (-A_hat / Delta_hat)
        tau_local = (Delta_hat / A_hat) if np.isfinite(A_hat) and A_hat != 0 else np.nan
        return {
            "ok": 1, "A": A_hat, "x_star": x_star_hat, "Delta": Delta_hat,
            "sse": sse, "slope_local": slope_local, "tau_local": tau_local,
            "nfev": int(res.nfev), "status": int(res.status), "message": str(res.message),
            "A0": float(A0), "xstar0": float(xstar0), "Delta0": float(Delta0),
        }
    except Exception:
        return {"ok": 0, "A": np.nan, "x_star": np.nan, "Delta": np.nan, "sse": np.inf, "nfev": np.nan, "status": np.nan, "message": ""}


def fit_hinge_plateau_smooth(x0, dx, w, weight_cap_q) -> Dict[str, float]:
    x, y, ww2 = _prep_xyw(x0, dx, w, weight_cap_q)
    if x.size < 50:
        return {"ok": 0, "a": np.nan, "b": np.nan, "A": np.nan, "x_c": np.nan, "Delta": np.nan,
                "sse": np.inf, "nfev": np.nan, "status": np.nan, "message": ""}

    sqrtw = np.sqrt(ww2)

    ou0 = fit_ou_linear_standardized(x, y, w=ww2, weight_cap_q=1.0)
    a0 = float(ou0["a"]) if ou0.get("ok", 0) == 1 else float(np.median(y))
    b0 = float(ou0["b"]) if ou0.get("ok", 0) == 1 else -0.1

    x_q = float(np.quantile(x, 0.90))
    tail = y[x >= x_q]
    A0 = float(-np.median(tail)) if tail.size > 20 else float(np.nanpercentile(np.abs(y), 75))
    if (not np.isfinite(A0)) or A0 <= 0:
        A0 = 0.2

    x_c0 = float(np.quantile(x, 0.80))
    Delta0 = float(np.std(x) * 0.2)
    if (not np.isfinite(Delta0)) or Delta0 <= 0:
        Delta0 = 0.3

    p0 = np.array([a0, b0, A0, x_c0, Delta0], float)
    lb = np.array([-np.inf, -np.inf, 1e-6, -np.inf, 1e-6], float)
    ub = np.array([ np.inf,  np.inf,  np.inf,  np.inf,  np.inf], float)

    def resid(p):
        a, b, A, x_c, Delta = p
        yhat = hinge_plateau_predict(x, a, b, A, x_c, Delta)
        return sqrtw * (y - yhat)

    try:
        res = least_squares(resid, x0=p0, bounds=(lb, ub), max_nfev=30000, loss="soft_l1", f_scale=1.0)
        if (not res.success) or (not np.all(np.isfinite(res.x))):
            return {"ok": 0, "a": np.nan, "b": np.nan, "A": np.nan, "x_c": np.nan, "Delta": np.nan,
                    "sse": np.inf, "nfev": float(res.nfev), "status": float(res.status), "message": str(res.message)}
        a_hat, b_hat, A_hat, x_c_hat, Delta_hat = map(float, res.x)
        if not (np.isfinite(Delta_hat) and Delta_hat > 0 and np.isfinite(A_hat) and A_hat > 0):
            return {"ok": 0, "a": np.nan, "b": np.nan, "A": np.nan, "x_c": np.nan, "Delta": np.nan,
                    "sse": np.inf, "nfev": float(res.nfev), "status": float(res.status), "message": str(res.message)}
        r = y - hinge_plateau_predict(x, a_hat, b_hat, A_hat, x_c_hat, Delta_hat)
        sse = float(np.sum(ww2 * r * r))
        tau_local = (-1.0 / b_hat) if np.isfinite(b_hat) and b_hat < 0 else np.nan
        return {
            "ok": 1, "a": a_hat, "b": b_hat, "A": A_hat, "x_c": x_c_hat, "Delta": Delta_hat,
            "sse": sse, "tau_from_b": float(tau_local),
            "nfev": int(res.nfev), "status": int(res.status), "message": str(res.message),
            "a0": float(a0), "b0": float(b0), "A0": float(A0), "x_c0": float(x_c0), "Delta0": float(Delta0),
        }
    except Exception:
        return {"ok": 0, "a": np.nan, "b": np.nan, "A": np.nan, "x_c": np.nan, "Delta": np.nan,
                "sse": np.inf, "nfev": np.nan, "status": np.nan, "message": ""}


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Transitions CSV")
    ap.add_argument("--data_dir", required=True, help="Directory for DATA outputs (CSVs/JSON).")

    ap.add_argument("--dt", type=int, default=1)
    ap.add_argument("--filter_class", default="", help="Value of obs_class to keep; empty disables.")
    ap.add_argument("--weight_col", default="", help="Optional weight column; empty => unweighted.")
    ap.add_argument("--weight_cap_q", type=float, default=0.999, help="Quantile cap for weights.")

    ap.add_argument("--grid_n", type=int, default=500, help="Number of points for drift grid.")
    ap.add_argument("--grid_qlo", type=float, default=0.02, help="Lower quantile for x-grid.")
    ap.add_argument("--grid_qhi", type=float, default=0.98, help="Upper quantile for x-grid.")

    ap.add_argument("--nbins_plot", type=int, default=20, help="Binning for binned_medians.csv.")
    ap.add_argument("--min_per_bin", type=int, default=200)

    ap.add_argument("--nbins_diff", type=int, default=0, help="Binning for diffusion_binned.csv (0 => use nbins_plot).")
    ap.add_argument("--min_per_bin_diff", type=int, default=0, help="Min per bin for diffusion_binned.csv (0 => use min_per_bin).")

    ap.add_argument("--bootstrap", type=int, default=0)
    ap.add_argument("--bootstrap_ci", type=float, default=0.95)
    ap.add_argument("--cluster_cols", default="subject,aaSeqCDR3")
    ap.add_argument("--seed", type=int, default=123)

    ap.add_argument(
        "--models",
        default="OU_linear,tanh_saturating,rational_saturating,hinge_plateau_smooth",
        help="Comma-separated list of models to fit."
    )
    return ap.parse_args()


# ============================================================
# Main
# ============================================================

def main() -> None:
    args = parse_args()
    data_out = ensure_dir(Path(args.data_dir))

    df = pd.read_csv(args.csv)
    need = {"dt", "x0", "dx"}
    miss = sorted(list(need - set(df.columns)))
    if miss:
        raise SystemExit(f"Missing required columns: {miss}")

    df = df[df["dt"] == int(args.dt)].copy()

    if args.filter_class:
        if "obs_class" not in df.columns:
            raise SystemExit("filter_class requested but obs_class missing in CSV.")
        df = df[df["obs_class"].astype(str) == str(args.filter_class)].copy()

    df["x0"] = pd.to_numeric(df["x0"], errors="coerce")
    df["dx"] = pd.to_numeric(df["dx"], errors="coerce")
    m = np.isfinite(df["x0"].to_numpy(float)) & np.isfinite(df["dx"].to_numpy(float))
    df = df[m].copy()
    if df.empty:
        raise SystemExit("No data after filtering dt/class and finite x0/dx.")

    weight_col = args.weight_col.strip()
    w = None
    if weight_col:
        if weight_col not in df.columns:
            raise SystemExit(f"weight_col '{weight_col}' not found in CSV.")
        df[weight_col] = pd.to_numeric(df[weight_col], errors="coerce")
        w = df[weight_col].to_numpy(float)

    x0 = df["x0"].to_numpy(float)
    dx = df["dx"].to_numpy(float)
    n = int(len(df))

    models = [s.strip() for s in args.models.split(",") if s.strip()]
    allowed = {"OU_linear", "tanh_saturating", "rational_saturating", "hinge_plateau_smooth"}
    bad = [m for m in models if m not in allowed]
    if bad:
        raise SystemExit(f"Unknown models: {bad}. Allowed: {sorted(allowed)}")

    fits: Dict[str, Dict[str, float]] = {}

    if "OU_linear" in models:
        fits["OU_linear"] = fit_ou_linear_standardized(x0, dx, w=w, weight_cap_q=float(args.weight_cap_q))
        if fits["OU_linear"].get("ok", 0) != 1:
            raise SystemExit("OU fit failed.")

    if "tanh_saturating" in models:
        fits["tanh_saturating"] = fit_tanh_saturating(x0, dx, w=w, weight_cap_q=float(args.weight_cap_q))
        if fits["tanh_saturating"].get("ok", 0) != 1:
            raise SystemExit("tanh_saturating fit failed.")

    if "rational_saturating" in models:
        fits["rational_saturating"] = fit_rational_saturating(x0, dx, w=w, weight_cap_q=float(args.weight_cap_q))
        if fits["rational_saturating"].get("ok", 0) != 1:
            raise SystemExit("rational_saturating fit failed.")

    if "hinge_plateau_smooth" in models:
        fits["hinge_plateau_smooth"] = fit_hinge_plateau_smooth(x0, dx, w=w, weight_cap_q=float(args.weight_cap_q))
        if fits["hinge_plateau_smooth"].get("ok", 0) != 1:
            raise SystemExit("hinge_plateau_smooth fit failed.")

    # --------------------------------------------------------
    # Fit summary + model comparison
    # --------------------------------------------------------
    fit_rows = []
    for name, f in fits.items():
        if f.get("ok", 0) != 1:
            continue

        k = {"OU_linear": 2, "tanh_saturating": 3, "rational_saturating": 3, "hinge_plateau_smooth": 5}[name]
        aic, bic = aic_bic(float(f["sse"]), n=n, k=k)

        row = {
            "model": name,
            "n": n,
            "k": k,
            "sse": float(f["sse"]),
            "aic": float(aic),
            "bic": float(bic),
        }

        if name == "OU_linear":
            row.update({
                "a": float(f["a"]),
                "b": float(f["b"]),
                "x_star": float(f["x_star"]),
                "tau": float(f["tau"]),
            })
        elif name in ("tanh_saturating", "rational_saturating"):
            row.update({
                "A": float(f["A"]),
                "x_star": float(f["x_star"]),
                "Delta": float(f["Delta"]),
                "slope_local": float(f.get("slope_local", np.nan)),
                "tau_local": float(f.get("tau_local", np.nan)),
                "nfev": f.get("nfev", np.nan),
                "status": f.get("status", np.nan),
            })
        elif name == "hinge_plateau_smooth":
            row.update({
                "a": float(f["a"]),
                "b": float(f["b"]),
                "A": float(f["A"]),
                "x_c": float(f["x_c"]),
                "Delta": float(f["Delta"]),
                "tau_from_b": float(f.get("tau_from_b", np.nan)),
                "nfev": f.get("nfev", np.nan),
                "status": f.get("status", np.nan),
            })

        row["weight_cap_q"] = float(args.weight_cap_q) if weight_col else np.nan
        fit_rows.append(row)

    fit_summary = pd.DataFrame(fit_rows).sort_values("aic", ascending=True)
    fit_summary.to_csv(data_out / "fit_summary.csv", index=False)

    if "OU_linear" in fits:
        base = "OU_linear"
    else:
        base = str(fit_summary.iloc[0]["model"]) if not fit_summary.empty else ""

    base_aic = float(fit_summary.loc[fit_summary["model"] == base, "aic"].iloc[0]) if base else np.nan
    base_bic = float(fit_summary.loc[fit_summary["model"] == base, "bic"].iloc[0]) if base else np.nan

    comp_rows = []
    for _, r in fit_summary.iterrows():
        comp_rows.append({
            "dt": int(args.dt),
            "obs_class": str(args.filter_class) if args.filter_class else "ALL",
            "weight_col": weight_col if weight_col else "unweighted",
            "weight_cap_q": float(args.weight_cap_q) if weight_col else np.nan,
            "n": n,
            "baseline": base,
            "model": str(r["model"]),
            "delta_aic": float(r["aic"] - base_aic) if np.isfinite(base_aic) else np.nan,
            "delta_bic": float(r["bic"] - base_bic) if np.isfinite(base_bic) else np.nan,
        })
    model_compare = pd.DataFrame(comp_rows)
    model_compare.to_csv(data_out / "model_compare.csv", index=False)

    # --------------------------------------------------------
    # Fitted values per transition
    # --------------------------------------------------------
    out = df.copy()

    if "OU_linear" in fits:
        out["dx_hat_OU_linear"] = ou_predict(x0, fits["OU_linear"]["a"], fits["OU_linear"]["b"])
        out["resid_OU_linear"] = out["dx"] - out["dx_hat_OU_linear"]
        out["resid2_OU_linear"] = out["resid_OU_linear"].to_numpy(float) ** 2

    if "tanh_saturating" in fits:
        f = fits["tanh_saturating"]
        out["dx_hat_tanh_saturating"] = tanh_predict(x0, f["A"], f["x_star"], f["Delta"])
        out["resid_tanh_saturating"] = out["dx"] - out["dx_hat_tanh_saturating"]
        out["resid2_tanh_saturating"] = out["resid_tanh_saturating"].to_numpy(float) ** 2

    if "rational_saturating" in fits:
        f = fits["rational_saturating"]
        out["dx_hat_rational_saturating"] = rational_predict(x0, f["A"], f["x_star"], f["Delta"])
        out["resid_rational_saturating"] = out["dx"] - out["dx_hat_rational_saturating"]
        out["resid2_rational_saturating"] = out["resid_rational_saturating"].to_numpy(float) ** 2

    if "hinge_plateau_smooth" in fits:
        f = fits["hinge_plateau_smooth"]
        out["dx_hat_hinge_plateau_smooth"] = hinge_plateau_predict(x0, f["a"], f["b"], f["A"], f["x_c"], f["Delta"])
        out["resid_hinge_plateau_smooth"] = out["dx"] - out["dx_hat_hinge_plateau_smooth"]
        out["resid2_hinge_plateau_smooth"] = out["resid_hinge_plateau_smooth"].to_numpy(float) ** 2

    out.to_csv(data_out / "fitted_values.csv", index=False)

    # --------------------------------------------------------
    # Drift function grid
    # --------------------------------------------------------
    if not (0.0 <= args.grid_qlo < args.grid_qhi <= 1.0):
        raise SystemExit("Invalid grid quantiles (need 0<=qlo<qhi<=1).")
    if int(args.grid_n) < 20:
        raise SystemExit("grid_n too small (use >=20).")

    x_lo = float(np.quantile(x0, float(args.grid_qlo)))
    x_hi = float(np.quantile(x0, float(args.grid_qhi)))
    if not (np.isfinite(x_lo) and np.isfinite(x_hi) and x_hi > x_lo):
        raise SystemExit("Failed to build x-grid from quantiles (check x0 after filters).")

    xs = np.linspace(x_lo, x_hi, int(args.grid_n))
    drift_grid = pd.DataFrame({
        "x": xs.astype(float),
        "dt": int(args.dt),
        "obs_class": str(args.filter_class) if args.filter_class else "ALL",
        "weight_col": weight_col if weight_col else "unweighted",
        "weight_cap_q": float(args.weight_cap_q) if weight_col else np.nan,
        "n": int(n),
        "grid_qlo": float(args.grid_qlo),
        "grid_qhi": float(args.grid_qhi),
        "grid_n": int(args.grid_n),
    })

    if "OU_linear" in fits:
        drift_grid["b_OU_linear"] = ou_predict(xs, fits["OU_linear"]["a"], fits["OU_linear"]["b"])
    if "tanh_saturating" in fits:
        f = fits["tanh_saturating"]
        drift_grid["b_tanh_saturating"] = tanh_predict(xs, f["A"], f["x_star"], f["Delta"])
    if "rational_saturating" in fits:
        f = fits["rational_saturating"]
        drift_grid["b_rational_saturating"] = rational_predict(xs, f["A"], f["x_star"], f["Delta"])
    if "hinge_plateau_smooth" in fits:
        f = fits["hinge_plateau_smooth"]
        drift_grid["b_hinge_plateau_smooth"] = hinge_plateau_predict(xs, f["a"], f["b"], f["A"], f["x_c"], f["Delta"])

    drift_grid.to_csv(data_out / "drift_function_grid.csv", index=False)

    # --------------------------------------------------------
    # Binned medians for plotting
    # --------------------------------------------------------
    frames = []
    b_dx = binned_median(x0, dx, nbins=int(args.nbins_plot), min_per_bin=int(args.min_per_bin))
    b_dx["series"] = "dx"
    frames.append(b_dx)

    for name in fits.keys():
        resid_col = f"resid_{name}"
        if resid_col in out.columns:
            br = binned_median(
                x0,
                out[resid_col].to_numpy(float),
                nbins=int(args.nbins_plot),
                min_per_bin=int(args.min_per_bin)
            )
            br["series"] = resid_col
            frames.append(br)

    bmed = pd.concat(frames, ignore_index=True)
    bmed["dt"] = int(args.dt)
    bmed["obs_class"] = str(args.filter_class) if args.filter_class else "ALL"
    bmed["weight_col"] = weight_col if weight_col else "unweighted"
    bmed.to_csv(data_out / "binned_medians.csv", index=False)

    # --------------------------------------------------------
    # Residual-based diffusion datasets
    # --------------------------------------------------------
    w_sanit = None
    if weight_col:
        w_sanit = sanitize_weights_for_fit(out[weight_col].to_numpy(float), cap_q=float(args.weight_cap_q))

    # Long-format residuals
    long_frames = []
    base_cols = ["dt", "x0", "dx"]
    for c in ["subject", "aaSeqCDR3", "time0", "time1", "obs_class"]:
        if c in out.columns:
            base_cols.append(c)

    if weight_col:
        out["_w_sanit"] = w_sanit
        base_cols.append(weight_col)
        base_cols.append("_w_sanit")

    for name in fits.keys():
        rcol = f"resid_{name}"
        pcol = f"dx_hat_{name}" if f"dx_hat_{name}" in out.columns else None
        if rcol not in out.columns:
            continue

        tmp = out[base_cols].copy()
        tmp["model"] = name
        if pcol is not None:
            tmp["dx_hat"] = out[pcol].to_numpy(float)
        tmp["resid"] = out[rcol].to_numpy(float)
        tmp["resid2"] = tmp["resid"].to_numpy(float) ** 2
        long_frames.append(tmp)

    diffusion_long = pd.concat(long_frames, ignore_index=True) if long_frames else pd.DataFrame()
    diffusion_long.to_csv(data_out / "diffusion_residuals_long.csv", index=False)

    # Binned diffusion summaries
    nbins_diff = int(args.nbins_diff) if int(args.nbins_diff) > 0 else int(args.nbins_plot)
    min_per_bin_diff = int(args.min_per_bin_diff) if int(args.min_per_bin_diff) > 0 else int(args.min_per_bin)

    bdiff_frames = []
    for name in fits.keys():
        rcol = f"resid_{name}"
        if rcol not in out.columns:
            continue

        bd = binned_diffusion(
            x0=out["x0"].to_numpy(float),
            resid=out[rcol].to_numpy(float),
            nbins=nbins_diff,
            min_per_bin=min_per_bin_diff,
            dt=float(args.dt),
            w=(w_sanit if weight_col else None),
        )
        bd["model"] = name
        bd["dt"] = int(args.dt)
        bd["obs_class"] = str(args.filter_class) if args.filter_class else "ALL"
        bd["weight_col"] = weight_col if weight_col else "unweighted"
        bd["weight_cap_q"] = float(args.weight_cap_q) if weight_col else np.nan
        bd["nbins_diff"] = nbins_diff
        bd["min_per_bin_diff"] = min_per_bin_diff
        bdiff_frames.append(bd)

    diffusion_binned = pd.concat(bdiff_frames, ignore_index=True) if bdiff_frames else pd.DataFrame()
    diffusion_binned.to_csv(data_out / "diffusion_binned.csv", index=False)

    # --------------------------------------------------------
    # Bootstrap
    # --------------------------------------------------------
    kept = 0
    skipped = 0
    reps_df = pd.DataFrame()

    if int(args.bootstrap) > 0:
        cluster_cols = [c.strip() for c in args.cluster_cols.split(",") if c.strip()]
        if not set(cluster_cols).issubset(df.columns):
            raise SystemExit(f"cluster_cols not present: {cluster_cols}")

        boots = cluster_bootstrap_samples(df, cluster_cols, reps=int(args.bootstrap), seed=int(args.seed))
        rows = []

        for i, bdf in enumerate(boots):
            x0b = pd.to_numeric(bdf["x0"], errors="coerce").to_numpy(float)
            dxb = pd.to_numeric(bdf["dx"], errors="coerce").to_numpy(float)
            mb = np.isfinite(x0b) & np.isfinite(dxb)
            x0b, dxb = x0b[mb], dxb[mb]
            if x0b.size < 500:
                skipped += 1
                continue

            wb = None
            if weight_col:
                wb = pd.to_numeric(bdf.loc[mb, weight_col], errors="coerce").to_numpy(float)

            bf: Dict[str, Dict[str, float]] = {}
            ok = True

            if "OU_linear" in models:
                bf["OU_linear"] = fit_ou_linear_standardized(x0b, dxb, w=wb, weight_cap_q=float(args.weight_cap_q))
                ok = ok and (bf["OU_linear"].get("ok", 0) == 1)

            if "tanh_saturating" in models:
                bf["tanh_saturating"] = fit_tanh_saturating(x0b, dxb, w=wb, weight_cap_q=float(args.weight_cap_q))
                ok = ok and (bf["tanh_saturating"].get("ok", 0) == 1)

            if "rational_saturating" in models:
                bf["rational_saturating"] = fit_rational_saturating(x0b, dxb, w=wb, weight_cap_q=float(args.weight_cap_q))
                ok = ok and (bf["rational_saturating"].get("ok", 0) == 1)

            if "hinge_plateau_smooth" in models:
                bf["hinge_plateau_smooth"] = fit_hinge_plateau_smooth(x0b, dxb, w=wb, weight_cap_q=float(args.weight_cap_q))
                ok = ok and (bf["hinge_plateau_smooth"].get("ok", 0) == 1)

            if not ok:
                skipped += 1
                continue

            nb = int(len(x0b))
            row = {"rep": int(i), "n": nb}

            ic = {}
            for name, f in bf.items():
                k = {"OU_linear": 2, "tanh_saturating": 3, "rational_saturating": 3, "hinge_plateau_smooth": 5}[name]
                aic, bic = aic_bic(float(f["sse"]), n=nb, k=k)
                ic[name] = {"aic": aic, "bic": bic}

            if "OU_linear" in bf:
                row.update({
                    "ou_a": float(bf["OU_linear"]["a"]),
                    "ou_b": float(bf["OU_linear"]["b"]),
                    "ou_x_star": float(bf["OU_linear"]["x_star"]),
                    "ou_tau": float(bf["OU_linear"]["tau"]),
                    "ou_sse": float(bf["OU_linear"]["sse"]),
                    "ou_aic": float(ic["OU_linear"]["aic"]),
                    "ou_bic": float(ic["OU_linear"]["bic"]),
                })

            for nm in ("tanh_saturating", "rational_saturating"):
                if nm in bf:
                    pfx = "tanh" if nm == "tanh_saturating" else "rat"
                    row.update({
                        f"{pfx}_A": float(bf[nm]["A"]),
                        f"{pfx}_Delta": float(bf[nm]["Delta"]),
                        f"{pfx}_x_star": float(bf[nm]["x_star"]),
                        f"{pfx}_tau_local": float(bf[nm].get("tau_local", np.nan)),
                        f"{pfx}_sse": float(bf[nm]["sse"]),
                        f"{pfx}_aic": float(ic[nm]["aic"]),
                        f"{pfx}_bic": float(ic[nm]["bic"]),
                    })

            if "hinge_plateau_smooth" in bf:
                row.update({
                    "hinge_a": float(bf["hinge_plateau_smooth"]["a"]),
                    "hinge_b": float(bf["hinge_plateau_smooth"]["b"]),
                    "hinge_A": float(bf["hinge_plateau_smooth"]["A"]),
                    "hinge_x_c": float(bf["hinge_plateau_smooth"]["x_c"]),
                    "hinge_Delta": float(bf["hinge_plateau_smooth"]["Delta"]),
                    "hinge_tau_from_b": float(bf["hinge_plateau_smooth"].get("tau_from_b", np.nan)),
                    "hinge_sse": float(bf["hinge_plateau_smooth"]["sse"]),
                    "hinge_aic": float(ic["hinge_plateau_smooth"]["aic"]),
                    "hinge_bic": float(ic["hinge_plateau_smooth"]["bic"]),
                })

            if "OU_linear" in ic:
                for nm in ic.keys():
                    if nm == "OU_linear":
                        continue
                    key = "tanh" if nm == "tanh_saturating" else ("rat" if nm == "rational_saturating" else "hinge")
                    row[f"delta_aic_{key}"] = float(ic[nm]["aic"] - ic["OU_linear"]["aic"])
                    row[f"delta_bic_{key}"] = float(ic[nm]["bic"] - ic["OU_linear"]["bic"])

            rows.append(row)
            kept += 1

        reps_df = pd.DataFrame(rows)
        reps_df.to_csv(data_out / "bootstrap_reps.csv", index=False)

        if not reps_df.empty:
            ci = float(args.bootstrap_ci)

            def summarize(series: pd.Series, name: str) -> Dict[str, float]:
                arr = series.to_numpy(float)
                med = float(np.nanmedian(arr)) if np.isfinite(arr).any() else np.nan
                lo, hi = quantile_ci(arr, ci=ci)
                return {"param": name, "median": med, "ci_lo": lo, "ci_hi": hi, "n_valid": int(np.isfinite(arr).sum())}

            summ_rows = []
            if "ou_tau" in reps_df.columns:
                summ_rows.append(summarize(reps_df["ou_tau"], "ou_tau"))

            for col, label in [
                ("delta_aic_tanh", "delta_aic (tanh - ou)"),
                ("delta_bic_tanh", "delta_bic (tanh - ou)"),
                ("delta_aic_rat", "delta_aic (rational - ou)"),
                ("delta_bic_rat", "delta_bic (rational - ou)"),
                ("delta_aic_hinge", "delta_aic (hinge - ou)"),
                ("delta_bic_hinge", "delta_bic (hinge - ou)"),
            ]:
                if col in reps_df.columns:
                    summ_rows.append(summarize(reps_df[col], label))

            boot_summary = pd.DataFrame(summ_rows)
            boot_summary.to_csv(data_out / "bootstrap_summary.csv", index=False)

    # --------------------------------------------------------
    # JSON summary
    # --------------------------------------------------------
    meta = {
        "input_csv": str(Path(args.csv).resolve()),
        "data_dir": str(data_out.resolve()),
        "filters": {
            "dt": int(args.dt),
            "obs_class": str(args.filter_class) if args.filter_class else "ALL",
        },
        "weight": {
            "col": weight_col if weight_col else "",
            "cap_q": float(args.weight_cap_q) if weight_col else None,
        },
        "models": models,
        "n": int(n),
        "grid": {
            "n": int(args.grid_n),
            "qlo": float(args.grid_qlo),
            "qhi": float(args.grid_qhi),
        },
        "plot_binning": {
            "nbins_plot": int(args.nbins_plot),
            "min_per_bin": int(args.min_per_bin),
        },
        "diffusion_binning": {
            "nbins_diff": nbins_diff,
            "min_per_bin_diff": min_per_bin_diff,
        },
        "bootstrap": {
            "requested": int(args.bootstrap),
            "kept": int(kept),
            "skipped": int(skipped),
            "ci": float(args.bootstrap_ci),
            "seed": int(args.seed),
            "cluster_cols": str(args.cluster_cols),
        },
        "info_criteria_note": (
            "AIC/BIC computed from SSE under a Gaussian residual approximation; "
            "used for relative comparison across drift models."
        ),
        "diffusion_note": (
            "diffusion_residuals_long.csv contains per-transition residuals. "
            "diffusion_binned.csv contains unsmoothed bin-wise conditional residual "
            "summaries and derived diffusion-scale proxies (D_var_mad, D_mean_resid2, D_mad2). "
            "Final D(x) for FP must be selected, interpolated, smoothed, and floored downstream."
        ),
    }

    (data_out / "fit_summary.json").write_text(
        json.dumps(
            {
                "meta": meta,
                "fit_summary": fit_summary.to_dict(orient="records"),
                "model_compare": model_compare.to_dict(orient="records"),
            },
            indent=2
        ),
        encoding="utf-8",
    )

    title_cls = (args.filter_class if args.filter_class else "ALL")
    title_w = (weight_col if weight_col else "unweighted")

    print("✅ Completed ANALYSIS (multi-model drift fits + residual-based diffusion summaries).")
    print(f"[DATA] n={n} dt={args.dt} class={title_cls} weight={title_w}")
    print(f"[MODELS] {', '.join(models)}")
    if weight_col:
        print(f"[WEIGHT] cap_q={float(args.weight_cap_q)} (clipped + rescaled to mean=1)")
    if int(args.bootstrap) > 0:
        print(f"[BOOT] kept={kept} skipped={skipped} requested={int(args.bootstrap)}")
    print(f"[PATH] data_dir: {data_out.resolve()}")
    print("[OUT] fit_summary.csv")
    print("[OUT] model_compare.csv")
    print("[OUT] fitted_values.csv")
    print("[OUT] drift_function_grid.csv")
    print("[OUT] binned_medians.csv")
    print("[OUT] diffusion_residuals_long.csv")
    print("[OUT] diffusion_binned.csv")


if __name__ == "__main__":
    main()