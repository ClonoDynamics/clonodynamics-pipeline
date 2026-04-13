#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_fit_drift_models_plot.py  (ROBUST v2.3 - multi-model, plot-only + residual spread QC)

Reads outputs from the analysis script:
  - drift_function_grid.csv   [required]
  - binned_medians.csv        [required]
  - diffusion_binned.csv      [optional but recommended]
  - bootstrap_reps.csv        [optional]
  - fit_summary.csv           [optional, used to infer best model if not passed]

Outputs
-------
Main figures
- fig_drift_compare.png
- fig_resid_median_by_bin.png
- fig_bootstrap_overlay.png                       [if bootstrap_reps.csv exists]

Residual-spread control figures
- fig_resid_spread_by_bin_var_mad.png            [if diffusion_binned.csv exists]
- fig_resid_spread_by_bin_mean_resid2.png        [if diffusion_binned.csv exists]
- fig_resid_spread_by_bin_mad.png                [if diffusion_binned.csv exists]
- fig_resid_spread_best_model.png                [if diffusion_binned.csv exists]
- fig_resid_spread_bin_counts.png                [if diffusion_binned.csv exists]

Rationale
---------
These extra figures are intended to check whether the residual spread pattern is:
1) robust to the statistic used (var_mad vs mean_resid2 vs mad),
2) driven by one specific model or shared across models,
3) potentially distorted by unstable bin counts.

Example
-------
python ./code/plotting/plot_drift_models_plot.py \
  --data_dir ./results/7-drift-multimodel/TT_dt1_unweighted/data \
  --fig_dir  ./figures/7-drift-multimodel/TT_dt1_unweighted \
  --dpi 300 --fig_w 8 --fig_h 3 \
  --legend_outside --legend_loc "upper left"
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Callable, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# IO / utils
# -----------------------------
def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def _as_float_col(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        raise SystemExit(f"Missing required column '{col}' in input.")
    return pd.to_numeric(df[col], errors="coerce").to_numpy(float)


def _finite_xy(x: np.ndarray, y: np.ndarray):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def _sort_xy(x: np.ndarray, y: np.ndarray):
    if len(x) == 0:
        return x, y
    order = np.argsort(x)
    return x[order], y[order]


# -----------------------------
# Model predictions (bootstrap overlay)
# -----------------------------
def ou_predict(x: np.ndarray, a: float, b: float) -> np.ndarray:
    return a + b * x


def tanh_predict(x: np.ndarray, A: float, x_star: float, Delta: float) -> np.ndarray:
    return A * np.tanh((x_star - x) / Delta)


def rational_predict(x: np.ndarray, A: float, x_star: float, Delta: float) -> np.ndarray:
    z = (x_star - x) / Delta
    return A * (z / np.sqrt(1.0 + z * z))


def _sigmoid(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, float)
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


def hinge_plateau_predict(x: np.ndarray, a: float, b: float, A: float, x_c: float, Delta: float) -> np.ndarray:
    s = 1.0 - _sigmoid((x - x_c) / Delta)
    return (a + b * x) * s + (-A) * (1.0 - s)


# -----------------------------
# Plot aesthetics
# -----------------------------
def apply_rcparams(font_family: str, font_size: float, axes_w: float, tick_w: float, tick_len: float) -> None:
    plt.rcParams.update({
        "font.family": font_family,
        "font.size": font_size,
        "axes.linewidth": axes_w,
        "xtick.major.width": tick_w,
        "ytick.major.width": tick_w,
        "xtick.major.size": tick_len,
        "ytick.major.size": tick_len,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def _maybe_set_limits(ax, xlim_lo, xlim_hi, ylim_lo, ylim_hi):
    if np.isfinite(xlim_lo) or np.isfinite(xlim_hi):
        lo = xlim_lo if np.isfinite(xlim_lo) else None
        hi = xlim_hi if np.isfinite(xlim_hi) else None
        ax.set_xlim(lo, hi)
    if np.isfinite(ylim_lo) or np.isfinite(ylim_hi):
        lo = ylim_lo if np.isfinite(ylim_lo) else None
        hi = ylim_hi if np.isfinite(ylim_hi) else None
        ax.set_ylim(lo, hi)


# -----------------------------
# Legend / save helpers
# -----------------------------
def apply_legend(ax: plt.Axes, args: argparse.Namespace) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    if getattr(args, "legend_outside", False):
        ax.legend(
            loc=args.legend_loc,
            bbox_to_anchor=(args.legend_bbox_x, args.legend_bbox_y),
            frameon=bool(args.legend_frame),
            fontsize=args.legend_fontsize,
            borderaxespad=0.0,
        )
    else:
        ax.legend(
            loc=args.legend_loc,
            frameon=bool(args.legend_frame),
            fontsize=args.legend_fontsize,
        )


def apply_layout(fig: plt.Figure, args: argparse.Namespace) -> None:
    if getattr(args, "legend_outside", False):
        fig.tight_layout(rect=[0.0, 0.0, float(args.legend_pad_right), 1.0])
    else:
        fig.tight_layout()


def save_png(fig: plt.Figure, path: Path, args: argparse.Namespace) -> None:
    if getattr(args, "legend_outside", False):
        fig.savefig(path, dpi=int(args.dpi), bbox_inches="tight", pad_inches=float(args.legend_pad_inches))
    else:
        fig.savefig(path, dpi=int(args.dpi))


def save_pdf(fig: plt.Figure, path: Path, args: argparse.Namespace) -> None:
    if getattr(args, "legend_outside", False):
        fig.savefig(path, bbox_inches="tight", pad_inches=float(args.legend_pad_inches))
    else:
        fig.savefig(path)


# -----------------------------
# Discovery helpers
# -----------------------------
def discover_b_cols(drift_grid: pd.DataFrame) -> List[str]:
    return [c for c in drift_grid.columns if str(c).startswith("b_")]


def discover_resid_series(bmed: pd.DataFrame) -> List[str]:
    if "series" not in bmed.columns:
        return []
    vals = bmed["series"].astype(str).unique().tolist()
    return sorted([s for s in vals if s.startswith("resid_")])


def discover_diff_models(diff_binned: pd.DataFrame) -> List[str]:
    if "model" not in diff_binned.columns:
        return []
    vals = diff_binned["model"].astype(str).dropna().unique().tolist()
    return sorted(vals)


def pretty_label(name: str) -> str:
    x = str(name)
    if x.startswith("b_"):
        x = x[2:]
    if x.startswith("resid_"):
        x = x[6:]
    return x.replace("_", " ")


def parse_show_list(arg: str) -> Optional[List[str]]:
    s = (arg or "").strip()
    if not s:
        return None
    return [t.strip() for t in s.split(",") if t.strip()]


def infer_best_model(data_dir: Path, fallback: Optional[str] = None) -> Optional[str]:
    fit_path = data_dir / "fit_summary.csv"
    if fit_path.exists():
        try:
            fs = pd.read_csv(fit_path)
            if {"model", "aic"}.issubset(fs.columns) and not fs.empty:
                fs = fs.copy()
                fs["aic"] = pd.to_numeric(fs["aic"], errors="coerce")
                fs = fs[np.isfinite(fs["aic"])]
                if not fs.empty:
                    fs = fs.sort_values("aic", ascending=True)
                    return str(fs.iloc[0]["model"])
        except Exception:
            pass
    return fallback


# -----------------------------
# Bootstrap mapping
# -----------------------------
def build_bootstrap_predictors(reps_df: pd.DataFrame) -> Dict[str, Callable[[pd.Series, np.ndarray], np.ndarray]]:
    preds: Dict[str, Callable[[pd.Series, np.ndarray], np.ndarray]] = {}
    cols = set(reps_df.columns.astype(str))

    if {"ou_a", "ou_b"}.issubset(cols):
        preds["OU_linear"] = lambda r, xs: ou_predict(xs, float(r["ou_a"]), float(r["ou_b"]))
    elif {"OU_a", "OU_b"}.issubset(cols):
        preds["OU_linear"] = lambda r, xs: ou_predict(xs, float(r["OU_a"]), float(r["OU_b"]))

    if {"log_A", "log_x_star", "log_Delta"}.issubset(cols):
        preds["tanh_saturating"] = lambda r, xs: tanh_predict(xs, float(r["log_A"]), float(r["log_x_star"]), float(r["log_Delta"]))
    if {"tanh_A", "tanh_x_star", "tanh_Delta"}.issubset(cols):
        preds["tanh_saturating"] = lambda r, xs: tanh_predict(xs, float(r["tanh_A"]), float(r["tanh_x_star"]), float(r["tanh_Delta"]))

    if {"rat_A", "rat_x_star", "rat_Delta"}.issubset(cols):
        preds["rational_saturating"] = lambda r, xs: rational_predict(xs, float(r["rat_A"]), float(r["rat_x_star"]), float(r["rat_Delta"]))

    if {"hinge_a", "hinge_b", "hinge_A", "hinge_x_c", "hinge_Delta"}.issubset(cols):
        preds["hinge_plateau_smooth"] = lambda r, xs: hinge_plateau_predict(
            xs, float(r["hinge_a"]), float(r["hinge_b"]), float(r["hinge_A"]), float(r["hinge_x_c"]), float(r["hinge_Delta"])
        )

    return preds


# -----------------------------
# Generic plot helper
# -----------------------------
def lineplot_by_group(
    df: pd.DataFrame,
    group_col: str,
    x_col: str,
    y_col: str,
    labels: Optional[List[str]],
    xlabel: str,
    ylabel: str,
    title: Optional[str],
    fig_path_png: Path,
    args: argparse.Namespace,
    fig_path_pdf: Optional[Path] = None,
    add_zero_hline: bool = False,
    xlim_lo: float = np.nan,
    xlim_hi: float = np.nan,
    ylim_lo: float = np.nan,
    ylim_hi: float = np.nan,
) -> bool:
    if df.empty:
        return False
    if group_col not in df.columns or x_col not in df.columns or y_col not in df.columns:
        return False

    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(1, 1, 1)

    plotted = False
    groups = labels if labels is not None else sorted(df[group_col].astype(str).dropna().unique().tolist())

    for g in groups:
        part = df[df[group_col].astype(str) == str(g)].copy()
        if part.empty:
            continue
        x = _as_float_col(part, x_col)
        y = _as_float_col(part, y_col)
        x, y = _finite_xy(x, y)
        x, y = _sort_xy(x, y)
        if len(x) == 0:
            continue
        ax.plot(x, y, linewidth=args.line_w, label=pretty_label(g))
        plotted = True

    if add_zero_hline:
        ax.axhline(0, linewidth=args.axline_w)

    if not plotted:
        ax.text(0.5, 0.5, "No finite data to plot", ha="center", va="center")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    _maybe_set_limits(ax, xlim_lo, xlim_hi, ylim_lo, ylim_hi)

    apply_legend(ax, args)
    if title and not args.no_titles:
        ax.set_title(title)

    apply_layout(fig, args)
    save_png(fig, fig_path_png, args)
    if fig_path_pdf is not None and args.pdf:
        save_pdf(fig, fig_path_pdf, args)
    plt.close(fig)
    return plotted


# -----------------------------
# CLI
# -----------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--fig_dir", required=True)

    # Typography
    ap.add_argument("--font_family", default="Arial")
    ap.add_argument("--font_size", type=float, default=10.0)

    # Geometry
    ap.add_argument("--fig_w", type=float, default=7.2)
    ap.add_argument("--fig_h", type=float, default=4.0)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--pdf", action="store_true")

    # Style
    ap.add_argument("--line_w", type=float, default=2.0)
    ap.add_argument("--axline_w", type=float, default=1.0)
    ap.add_argument("--tick_w", type=float, default=1.0)
    ap.add_argument("--tick_len", type=float, default=4.0)
    ap.add_argument("--axes_w", type=float, default=1.0)
    ap.add_argument("--marker_s", type=float, default=26.0)
    ap.add_argument("--marker_alpha", type=float, default=1.0)

    # Legend
    ap.add_argument("--legend_loc", default="best")
    ap.add_argument("--legend_frame", type=int, default=0)
    ap.add_argument("--legend_fontsize", type=float, default=10.0)
    ap.add_argument("--legend_outside", action="store_true")
    ap.add_argument("--legend_bbox_x", type=float, default=1.02)
    ap.add_argument("--legend_bbox_y", type=float, default=1.0)
    ap.add_argument("--legend_pad_right", type=float, default=0.78)
    ap.add_argument("--legend_pad_inches", type=float, default=0.05)

    # Labels
    ap.add_argument("--xlabel", default="x0")
    ap.add_argument("--ylabel_dx", default="dx")
    ap.add_argument("--ylabel_resid", default="median residual (dx - fit)")
    ap.add_argument("--ylabel_spread_var_mad", default="residual spread (var_mad)")
    ap.add_argument("--ylabel_spread_mean_resid2", default="residual spread (mean_resid2)")
    ap.add_argument("--ylabel_spread_mad", default="residual spread (mad)")
    ap.add_argument("--ylabel_bin_count", default="bin count (n)")

    # Limits
    ap.add_argument("--xlim_lo", type=float, default=np.nan)
    ap.add_argument("--xlim_hi", type=float, default=np.nan)
    ap.add_argument("--ylim_dx_lo", type=float, default=np.nan)
    ap.add_argument("--ylim_dx_hi", type=float, default=np.nan)
    ap.add_argument("--ylim_resid_lo", type=float, default=np.nan)
    ap.add_argument("--ylim_resid_hi", type=float, default=np.nan)
    ap.add_argument("--ylim_spread_lo", type=float, default=np.nan)
    ap.add_argument("--ylim_spread_hi", type=float, default=np.nan)
    ap.add_argument("--ylim_count_lo", type=float, default=np.nan)
    ap.add_argument("--ylim_count_hi", type=float, default=np.nan)

    # Titles
    ap.add_argument("--title_prefix", default="")
    ap.add_argument("--no_titles", action="store_true")

    # Models
    ap.add_argument("--show_models", default="", help="Comma-separated models to show. Empty=all.")
    ap.add_argument("--best_model", default="", help="Optional manual override for best-model residual spread plot.")

    # Bootstrap
    ap.add_argument("--bootstrap_draw", type=int, default=80)
    ap.add_argument("--bootstrap_alpha", type=float, default=0.15)
    ap.add_argument("--bootstrap_line_w", type=float, default=1.0)

    return ap.parse_args()


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    fig_dir = ensure_dir(Path(args.fig_dir))

    apply_rcparams(args.font_family, args.font_size, args.axes_w, args.tick_w, args.tick_len)

    grid_path = data_dir / "drift_function_grid.csv"
    bmed_path = data_dir / "binned_medians.csv"
    diff_path = data_dir / "diffusion_binned.csv"
    boot_path = data_dir / "bootstrap_reps.csv"

    if not grid_path.exists():
        raise SystemExit(f"Missing: {grid_path}")
    if not bmed_path.exists():
        raise SystemExit(f"Missing: {bmed_path}")

    drift_grid = pd.read_csv(grid_path)
    bmed = pd.read_csv(bmed_path)
    diff_binned = pd.read_csv(diff_path) if diff_path.exists() else pd.DataFrame()
    reps_df = pd.read_csv(boot_path) if boot_path.exists() else pd.DataFrame()

    if "x" not in drift_grid.columns:
        raise SystemExit("drift_function_grid.csv must contain column 'x'.")
    if "series" not in bmed.columns:
        raise SystemExit("binned_medians.csv must contain column 'series'.")
    for c in ["x0_mid", "y_med"]:
        if c not in bmed.columns:
            raise SystemExit(f"binned_medians.csv must contain column '{c}'.")

    dt = int(drift_grid["dt"].iloc[0]) if "dt" in drift_grid.columns else None
    obs_class = str(drift_grid["obs_class"].iloc[0]) if "obs_class" in drift_grid.columns else ""
    weight_col = str(drift_grid["weight_col"].iloc[0]) if "weight_col" in drift_grid.columns else ""
    n = int(drift_grid["n"].iloc[0]) if "n" in drift_grid.columns else None

    tcls = obs_class if obs_class else "ALL"
    tw = weight_col if weight_col else "unweighted"
    base_title = f"dt={dt}, class={tcls}, {tw}, n={n}"
    if args.title_prefix.strip():
        base_title = f"{args.title_prefix.strip()}  ({base_title})"

    show = parse_show_list(args.show_models)

    # -----------------------------
    # Drift compare
    # -----------------------------
    b_cols = discover_b_cols(drift_grid)
    if not b_cols:
        raise SystemExit("No drift curves found (need columns starting with 'b_').")
    if show is not None:
        keep = []
        for suf in show:
            c = f"b_{suf}" if not suf.startswith("b_") else suf
            if c in b_cols:
                keep.append(c)
        if not keep:
            raise SystemExit(f"--show_models requested {show} but none found among: {b_cols}")
        b_cols = keep

    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(1, 1, 1)

    bdx = bmed[bmed["series"].astype(str) == "dx"].copy()
    if not bdx.empty:
        x = _as_float_col(bdx, "x0_mid")
        y = _as_float_col(bdx, "y_med")
        x, y = _finite_xy(x, y)
        x, y = _sort_xy(x, y)
        ax.scatter(x, y, s=args.marker_s, alpha=args.marker_alpha, label="binned median(dx)")

    xs = pd.to_numeric(drift_grid["x"], errors="coerce").to_numpy(float)
    for c in b_cols:
        y = pd.to_numeric(drift_grid[c], errors="coerce").to_numpy(float)
        x2, y2 = _finite_xy(xs, y)
        x2, y2 = _sort_xy(x2, y2)
        if len(x2) == 0:
            continue
        ax.plot(x2, y2, linewidth=args.line_w, label=pretty_label(c))

    ax.axhline(0, linewidth=args.axline_w)
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel_dx)
    _maybe_set_limits(ax, args.xlim_lo, args.xlim_hi, args.ylim_dx_lo, args.ylim_dx_hi)
    apply_legend(ax, args)
    if not args.no_titles:
        ax.set_title(f"Finite-time drift fit: {base_title}")
    apply_layout(fig, args)
    save_png(fig, fig_dir / "fig_drift_compare.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "fig_drift_compare.pdf", args)
    plt.close(fig)

    # -----------------------------
    # Median residual by x0-bin
    # -----------------------------
    resid_series = discover_resid_series(bmed)
    if show is not None:
        resid_series = [f"resid_{s}" for s in show if f"resid_{s}" in resid_series]

    _ = lineplot_by_group(
        df=bmed[bmed["series"].astype(str).isin(resid_series)].copy() if resid_series else pd.DataFrame(),
        group_col="series",
        x_col="x0_mid",
        y_col="y_med",
        labels=resid_series if resid_series else None,
        xlabel=args.xlabel,
        ylabel=args.ylabel_resid,
        title=f"Median residual by x0-bin: {base_title}",
        fig_path_png=fig_dir / "fig_resid_median_by_bin.png",
        fig_path_pdf=(fig_dir / "fig_resid_median_by_bin.pdf"),
        args=args,
        add_zero_hline=True,
        xlim_lo=args.xlim_lo,
        xlim_hi=args.xlim_hi,
        ylim_lo=args.ylim_resid_lo,
        ylim_hi=args.ylim_resid_hi,
    )

    # -----------------------------
    # Residual spread controls
    # -----------------------------
    made_spread_plots = False
    if diff_binned.empty:
        print("[NOTE] diffusion_binned.csv missing -> residual spread QC plots skipped.")
    else:
        need_base = {"model", "x0_mid", "n"}
        if not need_base.issubset(diff_binned.columns):
            print("[NOTE] diffusion_binned.csv lacks required columns for spread QC plots.")
        else:
            diff_models = discover_diff_models(diff_binned)
            if show is not None:
                diff_models = [m for m in diff_models if m in show]

            # Plot 1: var_mad
            if "var_mad" in diff_binned.columns:
                ok = lineplot_by_group(
                    df=diff_binned[diff_binned["model"].astype(str).isin(diff_models)].copy(),
                    group_col="model",
                    x_col="x0_mid",
                    y_col="var_mad",
                    labels=diff_models,
                    xlabel=args.xlabel,
                    ylabel=args.ylabel_spread_var_mad,
                    title=f"Residual spread by x0-bin (var_mad): {base_title}",
                    fig_path_png=fig_dir / "fig_resid_spread_by_bin_var_mad.png",
                    fig_path_pdf=(fig_dir / "fig_resid_spread_by_bin_var_mad.pdf"),
                    args=args,
                    add_zero_hline=False,
                    xlim_lo=args.xlim_lo,
                    xlim_hi=args.xlim_hi,
                    ylim_lo=args.ylim_spread_lo,
                    ylim_hi=args.ylim_spread_hi,
                )
                made_spread_plots = made_spread_plots or ok

            # Plot 2: mean_resid2
            if "mean_resid2" in diff_binned.columns:
                ok = lineplot_by_group(
                    df=diff_binned[diff_binned["model"].astype(str).isin(diff_models)].copy(),
                    group_col="model",
                    x_col="x0_mid",
                    y_col="mean_resid2",
                    labels=diff_models,
                    xlabel=args.xlabel,
                    ylabel=args.ylabel_spread_mean_resid2,
                    title=f"Residual spread by x0-bin (mean_resid2): {base_title}",
                    fig_path_png=fig_dir / "fig_resid_spread_by_bin_mean_resid2.png",
                    fig_path_pdf=(fig_dir / "fig_resid_spread_by_bin_mean_resid2.pdf"),
                    args=args,
                    add_zero_hline=False,
                    xlim_lo=args.xlim_lo,
                    xlim_hi=args.xlim_hi,
                    ylim_lo=args.ylim_spread_lo,
                    ylim_hi=args.ylim_spread_hi,
                )
                made_spread_plots = made_spread_plots or ok

            # Plot 3: mad
            if "mad" in diff_binned.columns:
                ok = lineplot_by_group(
                    df=diff_binned[diff_binned["model"].astype(str).isin(diff_models)].copy(),
                    group_col="model",
                    x_col="x0_mid",
                    y_col="mad",
                    labels=diff_models,
                    xlabel=args.xlabel,
                    ylabel=args.ylabel_spread_mad,
                    title=f"Residual spread by x0-bin (mad): {base_title}",
                    fig_path_png=fig_dir / "fig_resid_spread_by_bin_mad.png",
                    fig_path_pdf=(fig_dir / "fig_resid_spread_by_bin_mad.pdf"),
                    args=args,
                    add_zero_hline=False,
                    xlim_lo=args.xlim_lo,
                    xlim_hi=args.xlim_hi,
                    ylim_lo=args.ylim_spread_lo,
                    ylim_hi=args.ylim_spread_hi,
                )
                made_spread_plots = made_spread_plots or ok

            # Plot 4: best model only
            best_model = args.best_model.strip() or infer_best_model(data_dir, fallback=(diff_models[0] if diff_models else None))
            if best_model and best_model in diff_models and "var_mad" in diff_binned.columns:
                one = diff_binned[diff_binned["model"].astype(str) == str(best_model)].copy()
                fig = plt.figure(figsize=(args.fig_w, args.fig_h))
                ax = fig.add_subplot(1, 1, 1)

                x = _as_float_col(one, "x0_mid")
                y = _as_float_col(one, "var_mad")
                nn = _as_float_col(one, "n")
                x, y = _finite_xy(x, y)
                x, y = _sort_xy(x, y)

                # n va riletto e ordinato coerentemente
                one2 = one.copy()
                one2["x0_mid"] = pd.to_numeric(one2["x0_mid"], errors="coerce")
                one2["var_mad"] = pd.to_numeric(one2["var_mad"], errors="coerce")
                one2["n"] = pd.to_numeric(one2["n"], errors="coerce")
                one2 = one2[np.isfinite(one2["x0_mid"]) & np.isfinite(one2["var_mad"])].copy()
                one2 = one2.sort_values("x0_mid")

                ax.plot(one2["x0_mid"].to_numpy(float), one2["var_mad"].to_numpy(float), linewidth=args.line_w, label=pretty_label(best_model))
                ax.scatter(one2["x0_mid"].to_numpy(float), one2["var_mad"].to_numpy(float), s=args.marker_s, alpha=args.marker_alpha)

                ax.set_xlabel(args.xlabel)
                ax.set_ylabel(args.ylabel_spread_var_mad)
                _maybe_set_limits(ax, args.xlim_lo, args.xlim_hi, args.ylim_spread_lo, args.ylim_spread_hi)
                apply_legend(ax, args)
                if not args.no_titles:
                    ax.set_title(f"Residual spread by x0-bin (best model = {pretty_label(best_model)}): {base_title}")

                apply_layout(fig, args)
                save_png(fig, fig_dir / "fig_resid_spread_best_model.png", args)
                if args.pdf:
                    save_pdf(fig, fig_dir / "fig_resid_spread_best_model.pdf", args)
                plt.close(fig)

            # Plot 5: bin counts
            ok = lineplot_by_group(
                df=diff_binned[diff_binned["model"].astype(str).isin(diff_models)].copy(),
                group_col="model",
                x_col="x0_mid",
                y_col="n",
                labels=diff_models,
                xlabel=args.xlabel,
                ylabel=args.ylabel_bin_count,
                title=f"Residual spread bin counts: {base_title}",
                fig_path_png=fig_dir / "fig_resid_spread_bin_counts.png",
                fig_path_pdf=(fig_dir / "fig_resid_spread_bin_counts.pdf"),
                args=args,
                add_zero_hline=False,
                xlim_lo=args.xlim_lo,
                xlim_hi=args.xlim_hi,
                ylim_lo=args.ylim_count_lo,
                ylim_hi=args.ylim_count_hi,
            )
            made_spread_plots = made_spread_plots or ok

    # -----------------------------
    # Bootstrap overlay
    # -----------------------------
    if not reps_df.empty:
        fig = plt.figure(figsize=(args.fig_w, args.fig_h))
        ax = fig.add_subplot(1, 1, 1)

        if not bdx.empty:
            x = _as_float_col(bdx, "x0_mid")
            y = _as_float_col(bdx, "y_med")
            x, y = _finite_xy(x, y)
            x, y = _sort_xy(x, y)
            ax.scatter(x, y, s=args.marker_s, alpha=args.marker_alpha, label="binned median(dx)")

        preds = build_bootstrap_predictors(reps_df)
        if show is not None:
            preds = {k: v for k, v in preds.items() if k in show}

        max_draw = int(min(max(args.bootstrap_draw, 0), len(reps_df)))
        if max_draw > 0 and preds:
            draw_idx = np.linspace(0, len(reps_df) - 1, max_draw).astype(int)
            for j in draw_idx:
                r = reps_df.iloc[j]
                for _, fn in preds.items():
                    try:
                        yb = fn(r, xs)
                        x2, y2 = _finite_xy(xs, yb)
                        x2, y2 = _sort_xy(x2, y2)
                        if len(x2) > 0:
                            ax.plot(x2, y2, linewidth=args.bootstrap_line_w, alpha=args.bootstrap_alpha)
                    except Exception:
                        continue

        for c in b_cols:
            y = pd.to_numeric(drift_grid[c], errors="coerce").to_numpy(float)
            x2, y2 = _finite_xy(xs, y)
            x2, y2 = _sort_xy(x2, y2)
            if len(x2) == 0:
                continue
            ax.plot(x2, y2, linewidth=args.line_w, label=pretty_label(c))

        ax.axhline(0, linewidth=args.axline_w)
        ax.set_xlabel(args.xlabel)
        ax.set_ylabel(args.ylabel_dx)
        _maybe_set_limits(ax, args.xlim_lo, args.xlim_hi, args.ylim_dx_lo, args.ylim_dx_hi)
        apply_legend(ax, args)
        if not args.no_titles:
            ax.set_title(f"Bootstrap overlay: {base_title}")
        apply_layout(fig, args)
        save_png(fig, fig_dir / "fig_bootstrap_overlay.png", args)
        if args.pdf:
            save_pdf(fig, fig_dir / "fig_bootstrap_overlay.pdf", args)
        plt.close(fig)

    # -----------------------------
    # Console summary
    # -----------------------------
    print("✅ Completed SINGLE-FIGURE plotting (multi-model + residual spread QC).")
    print(f"[PATH] fig_dir: {fig_dir.resolve()}")
    print("[OUT]  fig_drift_compare.png")
    print("[OUT]  fig_resid_median_by_bin.png")
    if not diff_binned.empty:
        if "var_mad" in diff_binned.columns:
            print("[OUT]  fig_resid_spread_by_bin_var_mad.png")
        if "mean_resid2" in diff_binned.columns:
            print("[OUT]  fig_resid_spread_by_bin_mean_resid2.png")
        if "mad" in diff_binned.columns:
            print("[OUT]  fig_resid_spread_by_bin_mad.png")
        if "n" in diff_binned.columns:
            print("[OUT]  fig_resid_spread_bin_counts.png")
        best_model = args.best_model.strip() or infer_best_model(data_dir, fallback=None)
        if best_model:
            print(f"[OUT]  fig_resid_spread_best_model.png  [best_model={best_model}]")
    else:
        print("[NOTE] diffusion_binned.csv missing -> no residual spread QC plots.")
    if not reps_df.empty:
        print("[OUT]  fig_bootstrap_overlay.png")
    else:
        print("[NOTE] bootstrap_reps.csv missing or empty -> no bootstrap overlay.")


if __name__ == "__main__":
    main()