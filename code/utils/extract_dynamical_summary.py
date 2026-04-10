#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extract_dynamical_summary_multi.py

Multi-condition extraction of dynamical summary metrics with single-condition and
comparative overlay plots for drift, diffusion, and flux.

What the script does
--------------------
For each condition listed in a manifest CSV, the script reads:
- a transitions CSV
- a drift CSV
- a flux CSV
- optionally, a diffusion CSV

Then it computes, for each condition:
- tau from the linear relation dx ~ slope * x0 over a central x0 range
- bootstrap CI for tau using cluster bootstrap by (subject, aaSeqCDR3)
- equilibrium point x* and f* from the zero-crossing of the drift curve
- fraction of large changes, defined by |dx| > threshold
- fraction of total |J(x)| concentrated below a chosen low-frequency threshold

Outputs
-------
1) summary_all_conditions.csv with one row per condition
2) Single-condition plots:
   - drift_equilibrium__COND.ext
   - diffusion__COND.ext        [if diffusion is available]
   - flux__COND.ext
3) Overlay plots across all conditions:
   - overlay_drift.ext
   - overlay_diffusion.ext      [only if at least one diffusion file is available]
   - overlay_flux.ext
4) Summary bar plots:
   - summary_tau_boot_median.ext
   - summary_x_star.ext
   - summary_percent_large_change.ext
   - summary_percent_flux_low_freq.ext

Manifest format
---------------
The manifest must be a CSV with at least these columns:

condition,transitions,drift,flux

Optionally, it can also contain:

diffusion

Example:
condition,transitions,drift,flux,diffusion
unweighted,./res/u/transitions.csv,./res/u/drift.csv,./res/u/flux.csv,./res/u/diffusion.csv
weight_w,./res/w/transitions.csv,./res/w/drift.csv,./res/w/flux.csv,./res/w/diffusion.csv
weight_clone,./res/c/transitions.csv,./res/c/drift.csv,./res/c/flux.csv,

If the diffusion column is absent, or if a row has an empty diffusion path,
that condition is still processed, but diffusion plots are skipped for it.

Expected columns in the input files
-----------------------------------
Transitions file:
- dt
- obs_class
- x0
- dx
- optionally subject, aaSeqCDR3

Drift file:
- one x column, default: x
- one drift column, default: b_hinge_plateau_smooth

Flux file:
- one x column, default: x
- one flux column, default: J_smooth

Diffusion file:
- one x column, default: x
- one diffusion column, default: D

Usage example
-------------
python extract_dynamical_summary_multi.py \
  --manifest ./manifest_conditions.csv \
  --out_data_dir ./results/dynamical_summary_multi \
  --out_fig_dir ./figures/dynamical_summary_multi \
  --dt 1 \
  --obs_class TT \
  --drift_col b_hinge_plateau_smooth \
  --x_col_drift x \
  --x_col_flux x \
  --J_col J_smooth \
  --x_col_diff x \
  --D_col D \
  --n_boot 300 \
  --dx_threshold log2 \
  --low_freq_threshold -9.0 \
  --linear_qlo 0.2 \
  --linear_qhi 0.8 \
  --fig_w 5.0 \
  --fig_h 3.6 \
  --dpi 300 \
  --fmt png \
  --font Arial \
  --font_size 10

Notes
-----
- The script uses matplotlib default colors.
- Diffusion plots are optional and depend on the manifest.
- Overlay plots are meant for direct comparison of conditions on the same axes.
"""

import argparse
import os
from typing import Dict, Any, Optional, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def parse_dx_threshold(val: str) -> float:
    return np.log(2) if val == "log2" else float(val)


def fit_tau(x: pd.Series, dx: pd.Series):
    if len(x) < 30:
        return np.nan, np.nan
    slope, _ = np.polyfit(x, dx, 1)
    tau = -1.0 / slope if slope < 0 else np.nan
    return slope, tau


def bootstrap_tau(df: pd.DataFrame, n_boot: int, seed: int):
    rng = np.random.default_rng(seed)

    if not {"subject", "aaSeqCDR3"}.issubset(df.columns):
        return np.nan, np.nan, np.nan

    clusters = df[["subject", "aaSeqCDR3"]].drop_duplicates().reset_index(drop=True)
    if len(clusters) == 0:
        return np.nan, np.nan, np.nan

    taus = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(clusters), size=len(clusters))
        sampled = clusters.iloc[idx].copy()
        merged = pd.merge(sampled, df, on=["subject", "aaSeqCDR3"], how="left")
        _, tau = fit_tau(merged["x0"], merged["dx"])
        taus.append(tau)

    taus = np.asarray(taus, dtype=float)
    taus = taus[np.isfinite(taus)]
    if len(taus) == 0:
        return np.nan, np.nan, np.nan

    return np.median(taus), np.percentile(taus, 2.5), np.percentile(taus, 97.5)


def estimate_equilibrium(drift: pd.DataFrame, x_col: str, b_col: str):
    x = drift[x_col].to_numpy(dtype=float)
    b = drift[b_col].to_numpy(dtype=float)

    mask = np.isfinite(x) & np.isfinite(b)
    x = x[mask]
    b = b[mask]

    if len(x) < 2:
        return np.nan, np.nan

    sign = np.sign(b)
    idx = np.where(sign[:-1] * sign[1:] < 0)[0]
    if len(idx) == 0:
        return np.nan, np.nan

    i = idx[0]
    x1, x2 = x[i], x[i + 1]
    b1, b2 = b[i], b[i + 1]

    if b2 == b1:
        return np.nan, np.nan

    x_star = x1 - b1 * (x2 - x1) / (b2 - b1)
    f_star = np.exp(x_star)
    return x_star, f_star


def load_optional_csv(path_value: Any) -> Optional[pd.DataFrame]:
    if path_value is None:
        return None
    if pd.isna(path_value):
        return None
    path = str(path_value).strip()
    if path == "":
        return None
    return pd.read_csv(path)


def load_and_summarize_condition(row: pd.Series, args) -> Dict[str, Any]:
    condition = str(row["condition"])
    trans = pd.read_csv(row["transitions"])
    drift = pd.read_csv(row["drift"])
    flux = pd.read_csv(row["flux"])
    diffusion = load_optional_csv(row["diffusion"]) if "diffusion" in row.index else None

    if "dt" not in trans.columns or "obs_class" not in trans.columns:
        raise ValueError(f"Condition {condition}: transitions file must contain 'dt' and 'obs_class'.")

    trans = trans[(trans["dt"] == args.dt) & (trans["obs_class"] == args.obs_class)].copy()

    if "subject" not in trans.columns:
        trans["subject"] = "all"
    if "aaSeqCDR3" not in trans.columns:
        trans["aaSeqCDR3"] = trans.index.astype(str)

    needed_trans_cols = {"x0", "dx"}
    if not needed_trans_cols.issubset(trans.columns):
        raise ValueError(f"Condition {condition}: transitions file must contain {sorted(needed_trans_cols)}.")

    if args.x_col_drift not in drift.columns or args.drift_col not in drift.columns:
        raise ValueError(
            f"Condition {condition}: drift file must contain '{args.x_col_drift}' and '{args.drift_col}'."
        )

    if args.x_col_flux not in flux.columns or args.J_col not in flux.columns:
        raise ValueError(
            f"Condition {condition}: flux file must contain '{args.x_col_flux}' and '{args.J_col}'."
        )

    if diffusion is not None:
        if args.x_col_diff not in diffusion.columns or args.D_col not in diffusion.columns:
            raise ValueError(
                f"Condition {condition}: diffusion file must contain '{args.x_col_diff}' and '{args.D_col}'."
            )

    qlo = trans["x0"].quantile(args.linear_qlo)
    qhi = trans["x0"].quantile(args.linear_qhi)
    sub = trans[(trans["x0"] >= qlo) & (trans["x0"] <= qhi)].copy()

    slope, tau = fit_tau(sub["x0"], sub["dx"])
    tau_med, tau_lo, tau_hi = bootstrap_tau(sub, args.n_boot, args.seed)

    x_star, f_star = estimate_equilibrium(drift, args.x_col_drift, args.drift_col)

    dx_thr = parse_dx_threshold(args.dx_threshold)
    trans["large"] = np.abs(trans["dx"]) > dx_thr
    frac_large = float(trans["large"].mean()) if len(trans) else np.nan

    x_flux = flux[args.x_col_flux].to_numpy(dtype=float)
    J = flux[args.J_col].to_numpy(dtype=float)
    J_abs = np.abs(J)
    J_sum = np.nansum(J_abs)
    if J_sum > 0:
        frac_low = np.nansum(J_abs[x_flux < args.low_freq_threshold]) / J_sum
    else:
        frac_low = np.nan

    summary = {
        "condition": condition,
        "n_transitions": int(len(trans)),
        "n_transitions_tau_fit": int(len(sub)),
        "slope_dx_vs_x0": slope,
        "tau_weeks": tau,
        "tau_boot_median": tau_med,
        "tau_ci_lo": tau_lo,
        "tau_ci_hi": tau_hi,
        "x_star": x_star,
        "f_star": f_star,
        "percent_large_change": frac_large * 100 if np.isfinite(frac_large) else np.nan,
        "percent_flux_low_freq": frac_low * 100 if np.isfinite(frac_low) else np.nan,
        "qlo_x0": qlo,
        "qhi_x0": qhi,
        "has_diffusion": diffusion is not None,
    }

    return {
        "condition": condition,
        "trans": trans,
        "drift": drift,
        "flux": flux,
        "diffusion": diffusion,
        "summary": summary,
    }


def apply_style(args) -> None:
    plt.rcParams.update({
        "font.family": args.font,
        "font.size": args.font_size,
        "axes.titlesize": args.font_size,
        "axes.labelsize": args.font_size,
        "xtick.labelsize": args.font_size,
        "ytick.labelsize": args.font_size,
        "legend.fontsize": max(args.font_size - 1, 7),
    })


def savefig(path: str, args) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=args.dpi, bbox_inches="tight")
    plt.close()


def plot_single_drift(res: Dict[str, Any], args, out_dir: str) -> None:
    condition = res["condition"]
    drift = res["drift"]
    x_star = res["summary"]["x_star"]

    plt.figure(figsize=(args.fig_w, args.fig_h))
    plt.plot(drift[args.x_col_drift], drift[args.drift_col], label=condition)
    plt.axhline(0, linestyle="--", linewidth=1)
    if np.isfinite(x_star):
        plt.axvline(x_star, linestyle="--", linewidth=1, label=f"x* = {x_star:.2f}")
    plt.xlabel("log-frequency")
    plt.ylabel("drift")
    plt.title(f"Drift: {condition}")
    plt.legend(frameon=False)
    savefig(os.path.join(out_dir, f"drift_equilibrium__{condition}.{args.fmt}"), args)


def plot_single_diffusion(res: Dict[str, Any], args, out_dir: str) -> None:
    condition = res["condition"]
    diffusion = res["diffusion"]
    if diffusion is None:
        return

    plt.figure(figsize=(args.fig_w, args.fig_h))
    plt.plot(diffusion[args.x_col_diff], diffusion[args.D_col], label=condition)
    plt.xlabel("log-frequency")
    plt.ylabel(args.D_col)
    plt.title(f"Diffusion: {condition}")
    plt.legend(frameon=False)
    savefig(os.path.join(out_dir, f"diffusion__{condition}.{args.fmt}"), args)


def plot_single_flux(res: Dict[str, Any], args, out_dir: str) -> None:
    condition = res["condition"]
    flux = res["flux"]

    plt.figure(figsize=(args.fig_w, args.fig_h))
    plt.plot(flux[args.x_col_flux], flux[args.J_col], label=condition)
    plt.axvline(args.low_freq_threshold, linestyle="--", linewidth=1,
                label=f"x = {args.low_freq_threshold:g}")
    plt.xlabel("log-frequency")
    plt.ylabel(args.J_col)
    plt.title(f"Flux: {condition}")
    plt.legend(frameon=False)
    savefig(os.path.join(out_dir, f"flux__{condition}.{args.fmt}"), args)


def plot_overlay_drift(results: List[Dict[str, Any]], args, out_dir: str) -> None:
    plt.figure(figsize=(args.fig_w, args.fig_h))
    for res in results:
        drift = res["drift"]
        condition = res["condition"]
        x_star = res["summary"]["x_star"]
        plt.plot(drift[args.x_col_drift], drift[args.drift_col], label=condition)
        if np.isfinite(x_star):
            plt.axvline(x_star, linestyle="--", linewidth=0.8)
    plt.axhline(0, linestyle="--", linewidth=1)
    plt.xlabel("log-frequency")
    plt.ylabel("drift")
    plt.title("Drift overlay")
    plt.legend(frameon=False)
    savefig(os.path.join(out_dir, f"overlay_drift.{args.fmt}"), args)


def plot_overlay_diffusion(results: List[Dict[str, Any]], args, out_dir: str) -> None:
    valid = [res for res in results if res["diffusion"] is not None]
    if not valid:
        return

    plt.figure(figsize=(args.fig_w, args.fig_h))
    for res in valid:
        diffusion = res["diffusion"]
        condition = res["condition"]
        plt.plot(diffusion[args.x_col_diff], diffusion[args.D_col], label=condition)
    plt.xlabel("log-frequency")
    plt.ylabel(args.D_col)
    plt.title("Diffusion overlay")
    plt.legend(frameon=False)
    savefig(os.path.join(out_dir, f"overlay_diffusion.{args.fmt}"), args)


def plot_overlay_flux(results: List[Dict[str, Any]], args, out_dir: str) -> None:
    plt.figure(figsize=(args.fig_w, args.fig_h))
    for res in results:
        flux = res["flux"]
        condition = res["condition"]
        plt.plot(flux[args.x_col_flux], flux[args.J_col], label=condition)
    plt.axvline(args.low_freq_threshold, linestyle="--", linewidth=1)
    plt.xlabel("log-frequency")
    plt.ylabel(args.J_col)
    plt.title("Flux overlay")
    plt.legend(frameon=False)
    savefig(os.path.join(out_dir, f"overlay_flux.{args.fmt}"), args)


def plot_summary_bars(summary_df: pd.DataFrame, args, out_dir: str) -> None:
    metrics = [
        ("tau_boot_median", "Tau (weeks)"),
        ("x_star", "Equilibrium x*"),
        ("percent_large_change", "% large changes"),
        ("percent_flux_low_freq", "% |J| at low freq"),
    ]

    for col, ylabel in metrics:
        vals = summary_df[col].to_numpy(dtype=float)
        labels = summary_df["condition"].tolist()
        fig_w = max(args.fig_w, 0.8 * len(labels) + 2.0)

        plt.figure(figsize=(fig_w, args.fig_h))
        xpos = np.arange(len(labels))
        plt.bar(xpos, vals)
        plt.xticks(xpos, labels, rotation=45, ha="right")
        plt.ylabel(ylabel)
        plt.title(ylabel)
        savefig(os.path.join(out_dir, f"summary_{col}.{args.fmt}"), args)


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("--manifest", required=True,
                        help="CSV with columns: condition,transitions,drift,flux[,diffusion]")
    parser.add_argument("--out_data_dir", required=True)
    parser.add_argument("--out_fig_dir", required=True)

    parser.add_argument("--dt", type=float, default=1)
    parser.add_argument("--obs_class", default="TT")

    parser.add_argument("--drift_col", default="b_hinge_plateau_smooth")
    parser.add_argument("--x_col_drift", default="x")

    parser.add_argument("--x_col_flux", default="x")
    parser.add_argument("--J_col", default="J_smooth")

    parser.add_argument("--x_col_diff", default="x")
    parser.add_argument("--D_col", default="D")

    parser.add_argument("--n_boot", type=int, default=300)
    parser.add_argument("--dx_threshold", default="log2")
    parser.add_argument("--low_freq_threshold", type=float, default=-9.0)
    parser.add_argument("--linear_qlo", type=float, default=0.2)
    parser.add_argument("--linear_qhi", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=123)

    parser.add_argument("--fig_w", type=float, default=5.0)
    parser.add_argument("--fig_h", type=float, default=3.6)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--fmt", default="png", choices=["png", "pdf", "svg"])
    parser.add_argument("--font", default="Arial")
    parser.add_argument("--font_size", type=int, default=10)

    args = parser.parse_args()

    ensure_dir(args.out_data_dir)
    ensure_dir(args.out_fig_dir)
    ensure_dir(os.path.join(args.out_fig_dir, "single"))
    ensure_dir(os.path.join(args.out_fig_dir, "overlay"))
    ensure_dir(os.path.join(args.out_fig_dir, "summary"))

    apply_style(args)

    manifest = pd.read_csv(args.manifest)
    needed_cols = {"condition", "transitions", "drift", "flux"}
    if not needed_cols.issubset(manifest.columns):
        raise ValueError(f"Manifest must contain columns: {sorted(needed_cols)}")

    results = []
    for _, row in manifest.iterrows():
        results.append(load_and_summarize_condition(row, args))

    summary_df = pd.DataFrame([r["summary"] for r in results])
    summary_df.to_csv(os.path.join(args.out_data_dir, "summary_all_conditions.csv"), index=False)

    for res in results:
        plot_single_drift(res, args, os.path.join(args.out_fig_dir, "single"))
        plot_single_diffusion(res, args, os.path.join(args.out_fig_dir, "single"))
        plot_single_flux(res, args, os.path.join(args.out_fig_dir, "single"))

    plot_overlay_drift(results, args, os.path.join(args.out_fig_dir, "overlay"))
    plot_overlay_diffusion(results, args, os.path.join(args.out_fig_dir, "overlay"))
    plot_overlay_flux(results, args, os.path.join(args.out_fig_dir, "overlay"))
    plot_summary_bars(summary_df, args, os.path.join(args.out_fig_dir, "summary"))

    print("\n===== DONE =====")
    print(f"Conditions processed: {len(results)}")
    print(f"Summary: {os.path.join(args.out_data_dir, 'summary_all_conditions.csv')}")
    print(f"Figures: {args.out_fig_dir}")


if __name__ == "__main__":
    main()
