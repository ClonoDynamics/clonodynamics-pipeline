#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
finite_time_fit.py

OU finite-time validation (ANALYSIS ONLY):
- Read finite_time_summary.csv (from finite_time_drift_diffusion.py)
- Calibrate OU from dt_ref:
    lambda_hat ≈ -slope(dt_ref)   (requires slope<0)
    D_hat from Var_OU(dt_ref)=var_mad_pooled(dt_ref)
- Predict Var_OU(dt) and D_app(dt)=Var_OU(dt)/(2 dt)
- Optional MC prediction bands using CI bounds at dt_ref (uniform in CI intervals)
- Save:
    out_dir/ou_pred_vs_empirical.csv
    out_dir/ou_fit_meta.json

NO PLOTTING here. Use 4b_plot_ou_finite_time_validation.py.

Example:
python ./code/dynamics/finite_time_fit.py \
  --summary_csv ./results/finite_time_moments_TT_consecutive/finite_time_summary.csv \
  --out_dir ./results/finite_time/ou_validation \
  --dt_ref 1 \
  --mc_reps 20000 --pred_ci 0.95 --mc_seed 123
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd


# -------------------------
# CLI
# -------------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="OU finite-time validation (analysis only).")
    ap.add_argument("--summary_csv", required=True, help="Path to finite_time_summary.csv")
    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--csv_sep", default=",", help="CSV separator for input (default ',')")

    ap.add_argument("--dt_ref", type=int, default=1, help="Reference dt (weeks) used to calibrate lambda and D (default 1)")
    ap.add_argument("--dt_values", default="", help="Optional comma list of dt to include (default: all in file)")

    # uncertainty propagation
    ap.add_argument("--mc_reps", type=int, default=20000,
                    help="Monte Carlo reps for OU prediction bands using CI bounds (default 20000; set 0 to disable)")
    ap.add_argument("--mc_seed", type=int, default=123)
    ap.add_argument("--pred_ci", type=float, default=0.95, help="CI level for prediction band (default 0.95)")
    ap.add_argument("--min_mc_ok", type=int, default=200,
                    help="Minimum number of valid MC samples required to compute bands (default 200).")

    # behavior
    ap.add_argument("--fail_on_nonneg_slope", action="store_true",
                    help="If set, raise error when slope(dt_ref) >= 0. Otherwise predictions are skipped but CSV/JSON still saved.")
    ap.add_argument("--quiet", action="store_true", help="Reduce stdout logs.")

    return ap


# -------------------------
# Helpers
# -------------------------
def parse_dt_values(s: str) -> Optional[np.ndarray]:
    s = (s or "").strip()
    if not s:
        return None
    return np.array([int(x.strip()) for x in s.split(",") if x.strip()], dtype=int)


def has_cols(df: pd.DataFrame, cols) -> bool:
    return all(c in df.columns for c in cols)


def ou_var(dt: np.ndarray, lam: float, D: float) -> np.ndarray:
    dt = np.asarray(dt, float)
    if not np.isfinite(lam) or not np.isfinite(D) or lam <= 0:
        return np.full(dt.shape, np.nan, float)
    return (D / lam) * (1.0 - np.exp(-2.0 * lam * dt))


def ou_D_app(dt: np.ndarray, lam: float, D: float) -> np.ndarray:
    v = ou_var(dt, lam, D)
    dt = np.asarray(dt, float)
    out = np.full(dt.shape, np.nan, float)
    m = np.isfinite(v) & (dt > 0)
    out[m] = v[m] / (2.0 * dt[m])
    return out


def qband(x: np.ndarray, ci: float) -> Tuple[float, float]:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (np.nan, np.nan)
    alpha = 0.5 * (1.0 - ci)
    return (float(np.quantile(x, alpha)), float(np.quantile(x, 1.0 - alpha)))


def log(msg: str, quiet: bool) -> None:
    if not quiet:
        print(msg)


# -------------------------
# Main
# -------------------------
def main() -> None:
    args = build_argparser().parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(Path(args.summary_csv).resolve(), sep=args.csv_sep)
    if "dt" not in df.columns:
        raise ValueError("Missing required column 'dt' in summary CSV.")

    df["dt"] = pd.to_numeric(df["dt"], errors="coerce")
    df = df[np.isfinite(df["dt"])].copy()
    df["dt"] = df["dt"].astype(int)
    df = df.sort_values("dt").reset_index(drop=True)

    dt_keep = parse_dt_values(args.dt_values)
    if dt_keep is not None:
        df = df[df["dt"].isin(dt_keep)].copy()
        df = df.sort_values("dt").reset_index(drop=True)

    # Required empirical columns for the validation
    if "var_mad_pooled" not in df.columns:
        raise ValueError("Missing 'var_mad_pooled' in summary CSV (needed for OU validation).")
    if "slope" not in df.columns:
        raise ValueError("Missing 'slope' in summary CSV (needed for OU validation).")

    # Ensure D_mad_pooled exists (plot script expects it). If missing, compute from var_mad_pooled.
    if "D_mad_pooled" not in df.columns:
        df["D_mad_pooled"] = np.where(df["dt"] > 0, df["var_mad_pooled"].astype(float) / (2.0 * df["dt"].astype(float)), np.nan)

    # Reference row
    ref = df[df["dt"] == int(args.dt_ref)]
    if ref.empty:
        raise ValueError(f"dt_ref={args.dt_ref} not found in summary CSV.")
    ref = ref.iloc[0]

    slope_ref = float(ref["slope"]) if np.isfinite(ref["slope"]) else np.nan
    var_ref = float(ref["var_mad_pooled"]) if np.isfinite(ref["var_mad_pooled"]) else np.nan
    lam_hat = -slope_ref if np.isfinite(slope_ref) else np.nan

    # Prepare base output frame
    out = df.copy()

    # Add CI cols if missing (keeps downstream robust)
    for c in ["var_mad_pooled_lo", "var_mad_pooled_hi", "D_mad_pooled_lo", "D_mad_pooled_hi"]:
        if c not in out.columns:
            out[c] = np.nan

    # If D CI missing but var CI exists, derive D CI from var CI (simple transformation)
    if (out["D_mad_pooled_lo"].isna().all() and out["D_mad_pooled_hi"].isna().all()
            and out["var_mad_pooled_lo"].notna().any() and out["var_mad_pooled_hi"].notna().any()):
        out["D_mad_pooled_lo"] = out["var_mad_pooled_lo"].astype(float) / (2.0 * out["dt"].astype(float))
        out["D_mad_pooled_hi"] = out["var_mad_pooled_hi"].astype(float) / (2.0 * out["dt"].astype(float))

    if not np.isfinite(lam_hat) or lam_hat <= 0:
        msg = (f"Cannot calibrate OU: slope(dt_ref={args.dt_ref}) must be finite and negative. "
               f"Got slope={slope_ref}.")
        out["lambda_hat_from_dtref"] = np.nan
        out["D_hat_from_dtref"] = np.nan
        out["var_ou_pred"] = np.nan
        out["var_ou_pred_lo"] = np.nan
        out["var_ou_pred_hi"] = np.nan
        out["D_ou_pred"] = np.nan
        out["D_ou_pred_lo"] = np.nan
        out["D_ou_pred_hi"] = np.nan

        out_path = out_dir / "ou_pred_vs_empirical.csv"
        out.to_csv(out_path, index=False)

        meta = {"error": msg, "dt_ref": int(args.dt_ref), "slope_ref": slope_ref, "var_ref": var_ref}
        (out_dir / "ou_fit_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

        if args.fail_on_nonneg_slope:
            raise ValueError(msg)
        log(f"[WARN] {msg}", args.quiet)
        log(f"[SAVED] {out_path}", args.quiet)
        return

    denom = 1.0 - np.exp(-2.0 * lam_hat * float(args.dt_ref))
    if not np.isfinite(var_ref) or var_ref < 0 or denom <= 0:
        raise ValueError(f"Cannot solve for D_hat: var_ref={var_ref}, denom={denom} (check inputs).")
    D_hat = lam_hat * var_ref / denom

    dts = out["dt"].to_numpy(int)
    var_pred = ou_var(dts, lam_hat, D_hat)
    D_app_pred = ou_D_app(dts, lam_hat, D_hat)

    # Optional prediction bands via Monte Carlo from CI bounds at dt_ref
    var_pred_lo = np.full_like(var_pred, np.nan, float)
    var_pred_hi = np.full_like(var_pred, np.nan, float)
    D_pred_lo = np.full_like(D_app_pred, np.nan, float)
    D_pred_hi = np.full_like(D_app_pred, np.nan, float)

    can_mc = (
        args.mc_reps > 0
        and has_cols(out, ["slope_lo", "slope_hi"])
        and (has_cols(out, ["var_mad_pooled_lo", "var_mad_pooled_hi"]))
        and (out["dt"] == int(args.dt_ref)).any()
    )

    used_mc = False
    if can_mc:
        rng = np.random.default_rng(int(args.mc_seed))
        B = int(args.mc_reps)
        ci = float(args.pred_ci)

        refrow = out[out["dt"] == int(args.dt_ref)].iloc[0]
        s_lo = float(refrow["slope_lo"])
        s_hi = float(refrow["slope_hi"])
        v_lo = float(refrow["var_mad_pooled_lo"])
        v_hi = float(refrow["var_mad_pooled_hi"])

        # Uniform sampling within CI at dt_ref (assumption-light)
        slope_s = rng.uniform(s_lo, s_hi, size=B)
        var_s = rng.uniform(max(0.0, v_lo), max(0.0, v_hi), size=B)

        lam_s = -slope_s
        ok = np.isfinite(lam_s) & (lam_s > 0) & np.isfinite(var_s) & (var_s >= 0)

        denom_s = 1.0 - np.exp(-2.0 * lam_s * float(args.dt_ref))
        ok = ok & np.isfinite(denom_s) & (denom_s > 0)

        lam_s = lam_s[ok]
        var_s = var_s[ok]
        denom_s = denom_s[ok]

        if lam_s.size >= int(args.min_mc_ok):
            D_s = lam_s * var_s / denom_s
            used_mc = True

            # Vectorized per-dt computation of bands
            for i, dt in enumerate(dts):
                v = (D_s / lam_s) * (1.0 - np.exp(-2.0 * lam_s * float(dt)))
                dapp = v / (2.0 * float(dt))

                var_pred_lo[i], var_pred_hi[i] = qband(v, ci)
                D_pred_lo[i], D_pred_hi[i] = qband(dapp, ci)

    out["lambda_hat_from_dtref"] = lam_hat
    out["D_hat_from_dtref"] = D_hat

    out["var_ou_pred"] = var_pred
    out["var_ou_pred_lo"] = var_pred_lo
    out["var_ou_pred_hi"] = var_pred_hi

    out["D_ou_pred"] = D_app_pred
    out["D_ou_pred_lo"] = D_pred_lo
    out["D_ou_pred_hi"] = D_pred_hi

    out_path = out_dir / "ou_pred_vs_empirical.csv"
    out.to_csv(out_path, index=False)

    meta = {
        "input": str(Path(args.summary_csv).resolve()),
        "dt_ref": int(args.dt_ref),
        "slope_ref": slope_ref,
        "lambda_hat": lam_hat,
        "var_ref": var_ref,
        "D_hat": D_hat,
        "used_mc_bands": bool(used_mc),
        "mc_reps_requested": int(args.mc_reps),
        "mc_seed": int(args.mc_seed),
        "pred_ci": float(args.pred_ci),
        "min_mc_ok": int(args.min_mc_ok),
        "notes": [
            "lambda_hat ≈ -slope(dt_ref) (requires slope<0).",
            "D_hat is solved from Var_OU(dt_ref) = var_mad_pooled(dt_ref).",
            "Prediction bands (if enabled) use uniform sampling within CI bounds at dt_ref.",
        ],
    }
    (out_dir / "ou_fit_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    log(f"[SAVED] {out_path}", args.quiet)
    log(f"[SAVED] {out_dir / 'ou_fit_meta.json'}", args.quiet)
    log(f"[OU CALIBRATION] dt_ref={args.dt_ref}  lambda_hat={lam_hat:.6g}  D_hat={D_hat:.6g}  (from var_ref={var_ref:.6g})", args.quiet)
    if can_mc and not used_mc:
        log("[WARN] MC bands requested but not enough valid samples / missing CI bounds at dt_ref.", args.quiet)


if __name__ == "__main__":
    main()