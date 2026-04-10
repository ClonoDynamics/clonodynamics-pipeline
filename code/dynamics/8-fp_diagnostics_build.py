#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-fp_diagnostics_build_obsclass_all.py
ROBUST v3.1 - build closed-vs-open FP diagnostic grids in log-frequency space

Key update versus v3.0
----------------------
This version adds explicit support for running the FP diagnostic build on:
  - a single observability class (e.g. TT)
  - a comma-separated subset of classes (e.g. TT,TF,FT)
  - all classes together (ALL or *)

This applies consistently to both:
  - drift_function_grid.csv filtering
  - diffusion_binned.csv filtering

Importantly, this is distinct from trajectory-level observability handling:
  --traj_observable_only still only affects the construction of p_emp(x)
  from trajectories_long.csv, not the transition-class filtering itself.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def pick_col(df: pd.DataFrame, candidates, ctx: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"[{ctx}] None of {candidates} found. Available columns: {list(df.columns)}")


def gaussian_smooth(y: np.ndarray, sigma_pts: float) -> np.ndarray:
    y = np.asarray(y, float)
    if sigma_pts is None or sigma_pts <= 0:
        return y.copy()
    r = int(max(3, np.ceil(4 * sigma_pts)))
    xs = np.arange(-r, r + 1)
    w = np.exp(-0.5 * (xs / sigma_pts) ** 2)
    w /= w.sum()
    ypad = np.pad(y, pad_width=r, mode="edge")
    out = np.convolve(ypad, w, mode="same")
    return out[r:-r]


def trapz_norm(p: np.ndarray, x: np.ndarray) -> np.ndarray:
    z = np.trapezoid(p, x)
    if not np.isfinite(z) or z <= 0:
        raise ValueError("Normalization failed: integral <= 0 or non-finite.")
    return p / z


def density_on_grid(x: np.ndarray, xgrid: np.ndarray, smooth_sigma_pts: float = 0.0) -> np.ndarray:
    xgrid = np.asarray(xgrid, float)
    x = np.asarray(x, float)
    edges = np.zeros(len(xgrid) + 1)
    edges[1:-1] = 0.5 * (xgrid[:-1] + xgrid[1:])
    edges[0] = xgrid[0] - (edges[1] - xgrid[0])
    edges[-1] = xgrid[-1] + (xgrid[-1] - edges[-2])
    counts, _ = np.histogram(x, bins=edges)
    dx = np.diff(edges)
    tot = counts.sum()
    if tot <= 0:
        raise ValueError("No samples available to build empirical density p_emp.")
    p = counts / (tot * dx)
    if smooth_sigma_pts and smooth_sigma_pts > 0:
        p = gaussian_smooth(p, smooth_sigma_pts)
        p = np.clip(p, 0, np.inf)
        p = trapz_norm(p, xgrid)
    return p


def build_closed_stationary(xgrid: np.ndarray, b: np.ndarray, D: np.ndarray) -> np.ndarray:
    xgrid = np.asarray(xgrid, float)
    b = np.asarray(b, float)
    D = np.asarray(D, float)
    Dsafe = np.clip(D, 1e-300, np.inf)
    integrand = b / Dsafe
    I = np.zeros_like(xgrid)
    I[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(xgrid))
    p_unn = (1.0 / Dsafe) * np.exp(I - np.max(I))
    return trapz_norm(p_unn, xgrid)


def compute_js_ks(xgrid: np.ndarray, p_emp: np.ndarray, p_closed: np.ndarray) -> tuple[float, float]:
    dx = np.gradient(xgrid)
    P = p_emp * dx
    Q = p_closed * dx
    P = P / P.sum()
    Q = Q / Q.sum()
    M = 0.5 * (P + Q)

    def KL(A, B):
        A = np.clip(A, 1e-300, 1.0)
        B = np.clip(B, 1e-300, 1.0)
        return float(np.sum(A * np.log(A / B)))

    JS = 0.5 * KL(P, M) + 0.5 * KL(Q, M)
    KS = float(np.max(np.abs(np.cumsum(P) - np.cumsum(Q))))
    return JS, KS


def parse_truthy_series(s: pd.Series) -> np.ndarray:
    if s.dtype == bool:
        return s.to_numpy()
    vals = s.astype(str).str.strip().str.lower()
    return vals.isin(["true", "1", "t", "yes", "y"]).to_numpy()


def filter_exact(df: pd.DataFrame, col: str, val, ctx: str) -> pd.DataFrame:
    if val is None:
        return df
    if col not in df.columns:
        raise ValueError(f"[{ctx}] Requested filter {col}={val} but column '{col}' not found.")
    return df.loc[df[col] == val].copy()


def filter_float(df: pd.DataFrame, col: str, val: float, ctx: str, tol: float = 1e-12) -> pd.DataFrame:
    if val is None:
        return df
    if col not in df.columns:
        raise ValueError(f"[{ctx}] Requested filter {col}={val} but column '{col}' not found.")
    c = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
    mask = np.isfinite(c) & (np.abs(c - val) <= tol)
    return df.loc[mask].copy()


def parse_obs_class_spec(spec: str) -> tuple[str, list[str] | None]:
    raw = str(spec).strip()
    if raw == "" or raw.upper() == "ALL" or raw == "*":
        return "ALL", None
    vals = [v.strip() for v in raw.split(",") if v.strip()]
    vals_norm = []
    for v in vals:
        u = v.upper()
        if u in {"ALL", "*"}:
            return "ALL", None
        vals_norm.append(u)
    vals_norm = list(dict.fromkeys(vals_norm))
    if not vals_norm:
        return "ALL", None
    return (vals_norm[0], vals_norm) if len(vals_norm) == 1 else ("MULTI", vals_norm)


def filter_obs_classes(df: pd.DataFrame, col: str, spec: str, ctx: str) -> tuple[pd.DataFrame, list[str] | None]:
    mode, classes = parse_obs_class_spec(spec)
    if col not in df.columns:
        raise ValueError(f"[{ctx}] Requested obs_class filtering but column '{col}' not found.")
    if classes is None:
        return df.copy(), None
    ser = df[col].astype(str).str.strip().str.upper()
    out = df.loc[ser.isin(classes)].copy()
    return out, classes


def extract_diffusion_from_binned(diff: pd.DataFrame, diffusion_source: str, dt: float) -> tuple[np.ndarray, str]:
    src = str(diffusion_source)
    if src == "var_mad":
        if "var_mad" not in diff.columns:
            raise ValueError("Requested diffusion_source='var_mad' but column not found.")
        arr = pd.to_numeric(diff["var_mad"], errors="coerce").to_numpy(float)
        return arr / (2.0 * float(dt)), "var_mad"
    if src == "mean_resid2":
        if "mean_resid2" not in diff.columns:
            raise ValueError("Requested diffusion_source='mean_resid2' but column not found.")
        arr = pd.to_numeric(diff["mean_resid2"], errors="coerce").to_numpy(float)
        return arr / (2.0 * float(dt)), "mean_resid2"
    if src == "mad2":
        if "mad" not in diff.columns:
            raise ValueError("Requested diffusion_source='mad2' but column 'mad' not found.")
        arr = pd.to_numeric(diff["mad"], errors="coerce").to_numpy(float)
        return (arr ** 2) / (2.0 * float(dt)), "mad^2"
    if src in ("D_var_mad", "D_mean_resid2", "D_mad2"):
        if src not in diff.columns:
            raise ValueError(f"Requested diffusion_source='{src}' but column not found.")
        arr = pd.to_numeric(diff[src], errors="coerce").to_numpy(float)
        return arr, src
    raise ValueError(
        f"Unknown diffusion_source='{src}'. Allowed: var_mad, mean_resid2, mad2, D_var_mad, D_mean_resid2, D_mad2."
    )


def main():
    ap = argparse.ArgumentParser(
        description="Build FP closed-vs-open diagnostic grids from fitted drift and residual-based diffusion summaries."
    )
    ap.add_argument("--drift_grid", required=True, help="drift_function_grid.csv")
    ap.add_argument("--diffusion_binned", required=True, help="diffusion_binned.csv")
    ap.add_argument("--trajectories", required=True, help="trajectories_long.csv")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--out_prefix", default="fp", help="Short prefix for output files (e.g. 'hinge_main').")
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument(
        "--obs_class",
        default="TT",
        help=(
            "Observability-class filter for drift/diffusion tables. "
            "Examples: 'TT', 'TT,TF,FT,FF', 'ALL', '*'."
        ),
    )
    ap.add_argument(
        "--model",
        default="hinge_plateau_smooth",
        choices=["OU_linear", "tanh_saturating", "rational_saturating", "hinge_plateau_smooth"],
        help="Which model's drift and diffusion slice to use.",
    )
    ap.add_argument("--weight_col", default=None, help="Optional exact filter, e.g. 'w_clone'.")
    ap.add_argument("--weight_cap_q", type=float, default=None, help="Optional exact float filter.")
    ap.add_argument(
        "--diffusion_source",
        default="D_var_mad",
        choices=["var_mad", "mean_resid2", "mad2", "D_var_mad", "D_mean_resid2", "D_mad2"],
        help="Which diffusion proxy from diffusion_binned.csv to use.",
    )
    ap.add_argument("--smooth_D_sigma_pts", type=float, default=2.0)
    ap.add_argument("--D_floor", type=float, default=1e-4)
    ap.add_argument("--traj_observable_only", action="store_true")
    ap.add_argument("--traj_observable_col", default="observable")
    ap.add_argument("--x_col_traj", default="")
    ap.add_argument("--smooth_p_sigma_pts", type=float, default=2.0)
    ap.add_argument("--smooth_J_sigma_pts", type=float, default=2.0)
    args = ap.parse_args()

    out = Path(args.out_dir)
    ensure_dir(out)

    drift = pd.read_csv(args.drift_grid)
    drift = filter_float(drift, "dt", args.dt, ctx="drift_grid")
    drift, obs_classes_effective = filter_obs_classes(drift, "obs_class", args.obs_class, ctx="drift_grid")
    if args.weight_col is not None:
        drift = filter_exact(drift, "weight_col", args.weight_col, ctx="drift_grid")
    if args.weight_cap_q is not None:
        drift = filter_float(drift, "weight_cap_q", args.weight_cap_q, ctx="drift_grid")
    if len(drift) == 0:
        raise ValueError("drift_grid is empty after filtering. Check --dt/--obs_class/--weight_col/--weight_cap_q.")

    xcol = pick_col(drift, ["x", "x_mid", "x0_mid", "xgrid", "grid_x"], "drift_x")
    bcol_map = {
        "OU_linear": "b_OU_linear",
        "tanh_saturating": "b_tanh_saturating",
        "rational_saturating": "b_rational_saturating",
        "hinge_plateau_smooth": "b_hinge_plateau_smooth",
    }
    bcol = bcol_map[args.model]
    if bcol not in drift.columns:
        raise ValueError(f"Drift column '{bcol}' not found. Available: {list(drift.columns)}")

    drift = drift[[xcol, bcol]].dropna().copy()
    drift[xcol] = pd.to_numeric(drift[xcol], errors="coerce")
    drift[bcol] = pd.to_numeric(drift[bcol], errors="coerce")
    drift = drift[np.isfinite(drift[xcol]) & np.isfinite(drift[bcol])].sort_values(xcol)
    xgrid = drift[xcol].to_numpy(float)
    b = drift[bcol].to_numpy(float)
    if len(xgrid) < 5:
        raise ValueError("Drift grid too small after filtering; need enough points to compute derivatives.")

    diff = pd.read_csv(args.diffusion_binned)
    diff = filter_float(diff, "dt", args.dt, ctx="diffusion_binned")
    diff, diff_obs_classes_effective = filter_obs_classes(diff, "obs_class", args.obs_class, ctx="diffusion_binned")
    diff = filter_exact(diff, "model", args.model, ctx="diffusion_binned")
    if args.weight_col is not None:
        diff = filter_exact(diff, "weight_col", args.weight_col, ctx="diffusion_binned")
    if args.weight_cap_q is not None:
        diff = filter_float(diff, "weight_cap_q", args.weight_cap_q, ctx="diffusion_binned")
    if len(diff) == 0:
        raise ValueError("diffusion_binned is empty after filtering. Check --dt/--obs_class/--model/--weight_*.")

    diff_xcol = pick_col(diff, ["x0_mid", "x_mid", "x", "x_center"], "diff_x")
    D_bins, D_var_source = extract_diffusion_from_binned(diff, args.diffusion_source, dt=float(args.dt))
    diff2 = diff[[diff_xcol]].copy()
    diff2[diff_xcol] = pd.to_numeric(diff2[diff_xcol], errors="coerce")
    diff2["D_bin"] = pd.to_numeric(D_bins, errors="coerce")
    diff2 = diff2[np.isfinite(diff2[diff_xcol]) & np.isfinite(diff2["D_bin"])].sort_values(diff_xcol)
    if len(diff2) < 2:
        raise ValueError("Not enough finite diffusion points after filtering and source selection.")

    D_raw_interp = np.interp(
        xgrid,
        diff2[diff_xcol].to_numpy(float),
        diff2["D_bin"].to_numpy(float),
        left=np.nan,
        right=np.nan,
    )
    if np.any(~np.isfinite(D_raw_interp)):
        good = np.isfinite(D_raw_interp)
        if good.sum() < 2:
            raise ValueError("Not enough finite diffusion points to interpolate D(x) onto drift xgrid.")
        D_raw_interp[~good] = np.interp(xgrid[~good], xgrid[good], D_raw_interp[good])

    D_smooth = gaussian_smooth(D_raw_interp, float(args.smooth_D_sigma_pts))
    D_smooth = np.where(np.isfinite(D_smooth), D_smooth, np.nan)
    D_final = np.clip(D_smooth, float(args.D_floor), np.inf)

    traj = pd.read_csv(args.trajectories)
    if args.traj_observable_only:
        if args.traj_observable_col not in traj.columns:
            raise ValueError(f"--traj_observable_only but '{args.traj_observable_col}' not found in trajectories.")
        mask = parse_truthy_series(traj[args.traj_observable_col])
        traj = traj.loc[mask].copy()

    if args.x_col_traj:
        if args.x_col_traj not in traj.columns:
            raise ValueError(f"--x_col_traj '{args.x_col_traj}' not found in trajectories.")
        xtraj = pd.to_numeric(traj[args.x_col_traj], errors="coerce").to_numpy(float)
    else:
        if "log_freq" in traj.columns:
            xtraj = pd.to_numeric(traj["log_freq"], errors="coerce").to_numpy(float)
        elif "x" in traj.columns:
            xtraj = pd.to_numeric(traj["x"], errors="coerce").to_numpy(float)
        elif "freq" in traj.columns:
            freq = pd.to_numeric(traj["freq"], errors="coerce").to_numpy(float)
            xtraj = np.log(np.clip(freq, 1e-300, np.inf))
        else:
            raise ValueError("trajectories must have one of: log_freq, x, freq (or pass --x_col_traj).")

    xtraj = xtraj[np.isfinite(xtraj)]
    if xtraj.size == 0:
        raise ValueError("No finite x values in trajectories after filtering.")

    p_emp = density_on_grid(xtraj, xgrid, smooth_sigma_pts=float(args.smooth_p_sigma_pts))
    p_closed = build_closed_stationary(xgrid, b=b, D=D_final)

    Dp = D_final * p_emp
    dDp_dx = np.gradient(Dp, xgrid)
    J_raw = b * p_emp - dDp_dx
    J_smooth = gaussian_smooth(J_raw, float(args.smooth_J_sigma_pts))
    S_raw = np.gradient(J_smooth, xgrid)
    S_smooth = gaussian_smooth(S_raw, float(args.smooth_J_sigma_pts))
    ratio = p_emp / np.clip(p_closed, 1e-300, np.inf)

    tag = str(args.out_prefix).strip() or "fp"
    grid = pd.DataFrame({
        "x": xgrid,
        "b": b,
        "D_raw_interp": D_raw_interp,
        "D_smooth": D_smooth,
        "D_final": D_final,
        "D_var_source": D_var_source,
        "p_emp": p_emp,
        "p_closed": p_closed,
        "ratio_pemp_over_pclosed": ratio,
        "J_raw": J_raw,
        "J_smooth": J_smooth,
        "S_raw": S_raw,
        "S_smooth": S_smooth,
    })
    grid_path = out / f"{tag}.grid.csv"
    grid.to_csv(grid_path, index=False)

    JS, KS = compute_js_ks(xgrid, p_emp, p_closed)
    summary = {
        "tag": tag,
        "drift_grid": str(args.drift_grid),
        "diffusion_binned": str(args.diffusion_binned),
        "trajectories": str(args.trajectories),
        "dt": float(args.dt),
        "obs_class_requested": args.obs_class,
        "obs_classes_effective": obs_classes_effective,
        "diffusion_obs_classes_effective": diff_obs_classes_effective,
        "model": args.model,
        "weight_col": args.weight_col,
        "weight_cap_q": args.weight_cap_q,
        "diffusion_source": args.diffusion_source,
        "D_var_source": D_var_source,
        "smooth_D_sigma_pts": float(args.smooth_D_sigma_pts),
        "D_floor": float(args.D_floor),
        "traj_observable_only": bool(args.traj_observable_only),
        "traj_observable_col": args.traj_observable_col,
        "x_col_traj": args.x_col_traj,
        "smooth_p_sigma_pts": float(args.smooth_p_sigma_pts),
        "smooth_J_sigma_pts": float(args.smooth_J_sigma_pts),
        "JS": float(JS),
        "KS": float(KS),
        "grid_csv": str(grid_path),
    }
    summary_json = out / f"{tag}.summary.json"
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    summary_txt = out / f"{tag}.summary.txt"
    with open(summary_txt, "w", encoding="utf-8") as f:
        for k, v in summary.items():
            if isinstance(v, dict) or isinstance(v, list):
                f.write(f"{k}={json.dumps(v)}\n")
            else:
                f.write(f"{k}={v}\n")

    print("✅ Completed FP diagnostic build.")
    print(f"[SAVED] {grid_path}")
    print(f"[SAVED] {summary_json}")
    print(f"[SAVED] {summary_txt}")
    print(f"[INFO] model={args.model}")
    print(f"[INFO] obs_class_requested={args.obs_class}")
    print(f"[INFO] obs_classes_effective={obs_classes_effective}")
    print(f"[INFO] diffusion_source={args.diffusion_source} -> stored as {D_var_source}")


if __name__ == "__main__":
    main()
