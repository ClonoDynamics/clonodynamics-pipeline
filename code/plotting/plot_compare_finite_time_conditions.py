#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_compare_finite_time_conditions.py

Compare pooled finite-time summaries across multiple analysis conditions.

This script supports two output families:
1) pooled finite-time summaries vs dt
2) optional binned finite-time curves vs initial frequency, across multiple dt

Main features
-------------
- stable export geometry for multipanel assembly
- full axis frame
- no bbox_inches="tight" on save
- optional outside shared legend
- optional removal of titles with --no_title
- adaptive margins to reduce clipping when figure size is small
- fixed ylabel x-position for more homogeneous panels
- more conservative right/top margins to avoid clipped legends and titles

Examples
python ./code/plotting/plot_compare_finite_time_conditions.py \
--cond "TT=./results/5-finite_time_diagnostics/TT_dt1-5/unweighted" \
--cond "ALL=./results/5-finite_time_diagnostics/ALL_dt1-5/unweighted" \
--cond "TT_weight=./results/5-finite_time_diagnostics/TT_dt1-5/weight_w" \
--cond "ALL_weight=./results/5-finite_time_diagnostics/ALL_dt1-5/weight_w" \
--dt 1 \
--out_dir ./figures/5-compare_conditions_dt1/figure_SX.png \
--panel_w 4 --panel_h 3 \
--no_title
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =========================================================
# CLI
# =========================================================
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Compare pooled finite-time parameters across multiple conditions."
    )

    ap.add_argument(
        "--cond",
        action="append",
        required=True,
        help='Condition spec: "LABEL=PATH". Repeat as needed.',
    )

    ap.add_argument(
        "--dt_values",
        default="",
        help="Optional comma-separated dt values to keep, e.g. 1,2,3,4,5",
    )

    ap.add_argument(
        "--metrics",
        default="x_star,slope,tau_drift,var_mad_pooled,D_mad_pooled",
        help=(
            "Comma-separated pooled metrics to plot. "
            "Valid: x_star,slope,tau_drift,var_mad_pooled,D_mad_pooled"
        ),
    )

    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--csv_sep", default=",", help="CSV separator. Default ','")

    ap.add_argument("--fmt", choices=["png", "pdf"], default="png", help="Output format")
    ap.add_argument("--dpi", type=int, default=300, help="DPI for PNG")

    ap.add_argument("--font", default="Arial", help="Matplotlib font family")
    ap.add_argument("--font_size", type=float, default=10.0, help="Base font size")

    ap.add_argument("--panel_w", type=float, default=5.0, help="Panel width in inches")
    ap.add_argument("--panel_h", type=float, default=3.5, help="Panel height in inches")

    ap.add_argument("--line_width", type=float, default=1.5, help="Line width")
    ap.add_argument("--marker", default="o", help="Marker style")
    ap.add_argument("--marker_size", type=float, default=3.5, help="Marker size")

    ap.add_argument("--show_ci", action="store_true", help="Draw CI ribbons if available")
    ap.add_argument("--ci_alpha", type=float, default=0.12, help="Alpha for CI ribbons")

    ap.add_argument("--legend_out", action="store_true", help="Place legend outside")
    ap.add_argument("--legend_ncol", type=int, default=1, help="Legend columns")

    ap.add_argument("--title_prefix", default="", help="Optional title prefix")
    ap.add_argument("--suptitle", default="", help="Optional figure super-title")
    ap.add_argument("--no_title", action="store_true", help="Remove panel titles")

    ap.add_argument("--ncols", type=int, default=3, help="Number of columns in pooled multipanel")
    ap.add_argument("--sharex", action="store_true", help="Share x-axis across pooled panels")

    # Fixed-geometry export controls
    ap.add_argument("--left_margin", type=float, default=0.12, help="Base figure left margin")
    ap.add_argument("--right_margin_in", type=float, default=0.97, help="Figure right margin when legend is inside")
    ap.add_argument("--right_margin_out", type=float, default=0.68, help="Figure right margin when legend is outside")
    ap.add_argument("--bottom_margin", type=float, default=0.16, help="Base figure bottom margin")
    ap.add_argument("--top_margin", type=float, default=0.90, help="Base figure top margin")
    ap.add_argument("--wspace", type=float, default=0.35, help="Subplot horizontal spacing")
    ap.add_argument("--hspace", type=float, default=0.35, help="Subplot vertical spacing")
    ap.add_argument("--ylabel_x", type=float, default=-0.10, help="Fixed ylabel x position in axes coordinates")
    ap.add_argument("--legend_x", type=float, default=0.70, help="Legend x anchor in figure coordinates when legend is outside")
    ap.add_argument("--legend_y", type=float, default=0.50, help="Legend y anchor in figure coordinates when legend is outside")

    # Optional: binned finite-time overlays
    ap.add_argument(
        "--also_plot_binned",
        action="store_true",
        help="Also compare binned finite-time curves from finite_time_by_dt_xbins_boot.csv",
    )
    ap.add_argument(
        "--binned_metrics",
        default="median_dx,var_mad_dx",
        help="Comma-separated binned metrics. Valid: median_dx,var_mad_dx,D_mad",
    )
    ap.add_argument("--binned_panel_w", type=float, default=5.0, help="Width of one binned panel")
    ap.add_argument("--binned_panel_h", type=float, default=3.5, help="Height of one binned panel")

    ap.add_argument("--debug", action="store_true", help="Print debug info")
    return ap


# =========================================================
# Metric metadata
# =========================================================
VALID_POOLED_METRICS = {
    "x_star": {
        "ylabel": r"$x^*$",
        "title": r"Equilibrium scale $x^*$ vs $\Delta t$",
        "fname": "compare_xstar_vs_dt",
    },
    "slope": {
        "ylabel": r"slope",
        "title": r"Local drift slope vs $\Delta t$",
        "fname": "compare_slope_vs_dt",
    },
    "tau_drift": {
        "ylabel": r"$\tau_{\mathrm{drift}}$",
        "title": r"Apparent drift timescale vs $\Delta t$",
        "fname": "compare_tau_vs_dt",
    },
    "var_mad_pooled": {
        "ylabel": r"$\mathrm{var}_{MAD}^{pooled}(\Delta x)$",
        "title": r"Pooled robust variance vs $\Delta t$",
        "fname": "compare_varmad_pooled_vs_dt",
    },
    "D_mad_pooled": {
        "ylabel": r"$D_{MAD}^{pooled}$",
        "title": r"Pooled diffusion proxy vs $\Delta t$",
        "fname": "compare_Dmad_pooled_vs_dt",
    },
}

VALID_BINNED_METRICS = {
    "median_dx": {
        "ylabel": r"median$(\Delta x)$",
        "title": r"Finite-time drift by initial frequency",
        "fname": "compare_binned_median_dx_multidt",
        "ci_lo": "median_lo",
        "ci_hi": "median_hi",
        "zero_line": True,
    },
    "var_mad_dx": {
        "ylabel": r"$\mathrm{var}_{MAD}(\Delta x)$",
        "title": r"Finite-time robust variance by initial frequency",
        "fname": "compare_binned_var_mad_dx_multidt",
        "ci_lo": "var_mad_lo",
        "ci_hi": "var_mad_hi",
        "zero_line": False,
    },
    "D_mad": {
        "ylabel": r"$D_{MAD}$",
        "title": r"Finite-time diffusion proxy by initial frequency",
        "fname": "compare_binned_D_mad_multidt",
        "ci_lo": None,
        "ci_hi": None,
        "zero_line": False,
    },
}


# =========================================================
# Style helpers
# =========================================================
def set_matplotlib_defaults(font: str, font_size: float) -> None:
    plt.rcParams["font.family"] = font
    plt.rcParams["font.size"] = font_size
    plt.rcParams["axes.titlesize"] = font_size
    plt.rcParams["axes.labelsize"] = font_size
    plt.rcParams["legend.fontsize"] = font_size
    plt.rcParams["xtick.labelsize"] = font_size
    plt.rcParams["ytick.labelsize"] = font_size
    plt.rcParams["axes.titlepad"] = 6
    plt.rcParams["axes.labelpad"] = 3
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["figure.facecolor"] = "white"


def style_axes(ax: plt.Axes, ylabel_x: float) -> None:
    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(True)
    ax.tick_params(direction="out")
    ax.yaxis.set_label_coords(ylabel_x, 0.5)


def maybe_set_title(ax: plt.Axes, title: str, args) -> None:
    if not args.no_title:
        final_title = f"{args.title_prefix}{title}" if args.title_prefix else title
        ax.set_title(final_title)


def parse_conditions(cond_args: List[str]) -> List[Tuple[str, Path]]:
    out: List[Tuple[str, Path]] = []
    for item in cond_args:
        if "=" not in item:
            raise ValueError(f"Invalid --cond '{item}'. Use LABEL=PATH")
        label, path = item.split("=", 1)
        label = label.strip()
        path = Path(path.strip()).resolve()
        if not label:
            raise ValueError(f"Empty label in --cond '{item}'")
        out.append((label, path))
    return out


def parse_list_str(s: str) -> List[str]:
    return [x.strip() for x in (s or "").split(",") if x.strip()]


def parse_dt_values(s: str) -> List[int]:
    vals = parse_list_str(s)
    return [int(v) for v in vals]


def read_required_csv(path: Path, sep: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path, sep=sep)


def ensure_columns(df: pd.DataFrame, required: List[str], path: Path) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing columns {missing} in {path}\n"
            f"Available columns: {list(df.columns)}"
        )


def normalize_dt_column(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    out = df.copy()
    if "dt" not in out.columns:
        raise ValueError(f"'dt' column not found in {path}")
    out["dt"] = pd.to_numeric(out["dt"], errors="coerce")
    out = out.dropna(subset=["dt"]).copy()
    out["dt"] = out["dt"].astype(int)
    return out


def save_figure(fig: plt.Figure, out_base: Path, fmt: str, dpi: int) -> None:
    # Deliberately avoid bbox_inches="tight" to preserve panel geometry
    if fmt == "png":
        fig.savefig(out_base.with_suffix(".png"), dpi=dpi)
    else:
        fig.savefig(out_base.with_suffix(".pdf"))


def place_shared_legend(fig: plt.Figure, handles, labels, args) -> None:
    if not handles:
        return

    if args.legend_out:
        fig.legend(
            handles,
            labels,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(args.legend_x, args.legend_y),
            bbox_transform=fig.transFigure,
            borderaxespad=0.0,
            ncol=max(1, int(args.legend_ncol)),
        )
    else:
        fig.legend(
            handles,
            labels,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.995),
            bbox_transform=fig.transFigure,
            borderaxespad=0.0,
            ncol=max(1, int(args.legend_ncol)),
        )


def adaptive_margins(
    width_in: float,
    height_in: float,
    legend_out: bool,
    args,
    has_suptitle: bool = False,
) -> Tuple[float, float, float, float]:
    """
    Conservative adaptive margins to avoid clipping of:
    - y labels on small figures
    - x labels on short figures
    - panel titles / suptitle on top
    - outside legends on the right
    """
    if legend_out:
        if width_in < 6.0:
            right = min(args.right_margin_out, 0.72)
        elif width_in < 10.0:
            right = min(args.right_margin_out, 0.76)
        else:
            right = min(args.right_margin_out, 0.80)
    else:
        right = args.right_margin_in

    if width_in < 4.0:
        left = max(args.left_margin, 0.20)
    elif width_in < 5.0:
        left = max(args.left_margin, 0.16)
    else:
        left = max(args.left_margin, 0.12)

    if height_in < 3.0:
        bottom = max(args.bottom_margin, 0.22)
    elif height_in < 3.5:
        bottom = max(args.bottom_margin, 0.19)
    else:
        bottom = max(args.bottom_margin, 0.16)

    if height_in < 3.0:
        top = min(0.84, args.top_margin)
    elif height_in < 3.5:
        top = min(0.86, args.top_margin)
    else:
        top = min(0.88, args.top_margin)

    if not getattr(args, "no_title", False):
        top -= 0.03

    if has_suptitle:
        top -= 0.05

    top = max(top, bottom + 0.18)

    return left, right, bottom, top


def apply_fixed_layout(
    fig: plt.Figure,
    args,
    width_in: float,
    height_in: float,
    has_suptitle: bool = False,
) -> None:
    left, right, bottom, top = adaptive_margins(
        width_in=width_in,
        height_in=height_in,
        legend_out=args.legend_out,
        args=args,
        has_suptitle=has_suptitle,
    )
    fig.subplots_adjust(
        left=left,
        right=right,
        bottom=bottom,
        top=top,
        wspace=args.wspace,
        hspace=args.hspace,
    )


def add_ci_band(ax: plt.Axes, x: np.ndarray, df: pd.DataFrame, lo_col: str, hi_col: str, alpha: float) -> None:
    if {lo_col, hi_col}.issubset(df.columns):
        lo = pd.to_numeric(df[lo_col], errors="coerce").to_numpy(float)
        hi = pd.to_numeric(df[hi_col], errors="coerce").to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(lo) & np.isfinite(hi)
        if ok.any():
            ax.fill_between(x[ok], lo[ok], hi[ok], alpha=alpha)


# =========================================================
# Data loading
# =========================================================
def load_summary(
    cdir: Path,
    csv_sep: str,
    dt_values: List[int],
    metrics: List[str],
) -> pd.DataFrame:
    path = cdir / "finite_time_summary.csv"
    df = read_required_csv(path, csv_sep)
    df = normalize_dt_column(df, path)

    required = ["dt"] + metrics
    ensure_columns(df, required, path)

    if dt_values:
        df = df[df["dt"].isin(dt_values)].copy()

    if df.empty:
        raise ValueError(f"No rows remain in {path} after dt filtering.")

    return df.sort_values("dt").reset_index(drop=True)


def load_binned(
    cdir: Path,
    csv_sep: str,
    dt_values: List[int],
    binned_metrics: List[str],
) -> pd.DataFrame:
    path = cdir / "finite_time_by_dt_xbins_boot.csv"
    df = read_required_csv(path, csv_sep)
    df = normalize_dt_column(df, path)

    required = ["dt", "x_center"] + binned_metrics
    ensure_columns(df, required, path)

    if dt_values:
        df = df[df["dt"].isin(dt_values)].copy()

    if df.empty:
        raise ValueError(f"No rows remain in {path} after dt filtering.")

    return df.sort_values(["dt", "x_center"]).reset_index(drop=True)


# =========================================================
# Pooled plots
# =========================================================
def plot_pooled_metric(
    ax: plt.Axes,
    metric: str,
    summaries: Dict[str, pd.DataFrame],
    args,
) -> Tuple[List, List]:
    handles = []
    labels = []

    meta = VALID_POOLED_METRICS[metric]
    lo_col = f"{metric}_lo"
    hi_col = f"{metric}_hi"

    for label, df in summaries.items():
        x = df["dt"].to_numpy(int)
        y = pd.to_numeric(df[metric], errors="coerce").to_numpy(float)

        ok = np.isfinite(x) & np.isfinite(y)
        line, = ax.plot(
            x[ok],
            y[ok],
            lw=args.line_width,
            marker=args.marker,
            ms=args.marker_size,
            label=label,
        )
        handles.append(line)
        labels.append(label)

        if args.show_ci and {lo_col, hi_col}.issubset(df.columns):
            add_ci_band(ax, x.astype(float), df, lo_col, hi_col, args.ci_alpha)

    ax.set_xlabel(r"$\Delta t$ (weeks)")
    ax.set_ylabel(meta["ylabel"])
    maybe_set_title(ax, meta["title"], args)

    xticks = sorted(np.unique(np.concatenate([df["dt"].to_numpy(int) for df in summaries.values()])))
    ax.set_xticks(xticks)

    vals = []
    for df in summaries.values():
        vals.extend(pd.to_numeric(df[metric], errors="coerce").tolist())
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    if vals.size and np.nanmin(vals) < 0 < np.nanmax(vals):
        ax.axhline(0.0, linewidth=1.0)

    style_axes(ax, args.ylabel_x)
    return handles, labels


def make_pooled_multipanel(
    metrics: List[str],
    summaries: Dict[str, pd.DataFrame],
    out_dir: Path,
    args,
) -> None:
    n = len(metrics)
    ncols = max(1, int(args.ncols))
    nrows = int(np.ceil(n / ncols))

    fig_w = args.panel_w * ncols
    fig_h = args.panel_h * nrows

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(fig_w, fig_h),
        sharex=args.sharex,
        squeeze=False,
    )

    axes_flat = axes.ravel().tolist()
    all_handles = []
    all_labels = []

    for i, metric in enumerate(metrics):
        ax = axes_flat[i]
        handles, labels = plot_pooled_metric(ax, metric, summaries, args)
        if i == 0:
            all_handles = handles
            all_labels = labels

    for j in range(len(metrics), len(axes_flat)):
        axes_flat[j].axis("off")

    has_suptitle = bool(args.suptitle)
    if has_suptitle:
        fig.suptitle(args.suptitle)

    place_shared_legend(fig, all_handles, all_labels, args)
    apply_fixed_layout(fig, args, width_in=fig_w, height_in=fig_h, has_suptitle=has_suptitle)
    save_figure(fig, out_dir / "compare_finite_time_pooled_multipanel", args.fmt, args.dpi)
    plt.close(fig)


def make_pooled_single_figures(
    metrics: List[str],
    summaries: Dict[str, pd.DataFrame],
    out_dir: Path,
    args,
) -> None:
    for metric in metrics:
        fig, ax = plt.subplots(figsize=(args.panel_w, args.panel_h))
        handles, labels = plot_pooled_metric(ax, metric, summaries, args)
        place_shared_legend(fig, handles, labels, args)
        apply_fixed_layout(
            fig,
            args,
            width_in=args.panel_w,
            height_in=args.panel_h,
            has_suptitle=False,
        )
        save_figure(fig, out_dir / VALID_POOLED_METRICS[metric]["fname"], args.fmt, args.dpi)
        plt.close(fig)


# =========================================================
# Binned plots
# =========================================================
def plot_binned_metric_panel(
    ax: plt.Axes,
    metric: str,
    dt: int,
    binned_data: Dict[str, pd.DataFrame],
    args,
) -> Tuple[List, List]:
    handles = []
    labels = []

    meta = VALID_BINNED_METRICS[metric]
    lo_col = meta["ci_lo"]
    hi_col = meta["ci_hi"]

    for label, df in binned_data.items():
        d = df[df["dt"] == int(dt)].copy().sort_values("x_center")
        x = pd.to_numeric(d["x_center"], errors="coerce").to_numpy(float)

        ycol = "D_mad" if metric == "D_mad" else metric
        y = pd.to_numeric(d[ycol], errors="coerce").to_numpy(float)

        ok = np.isfinite(x) & np.isfinite(y)
        line, = ax.plot(
            x[ok],
            y[ok],
            lw=args.line_width,
            marker=args.marker,
            ms=args.marker_size,
            label=label,
        )
        handles.append(line)
        labels.append(label)

        if args.show_ci and lo_col and hi_col and {lo_col, hi_col}.issubset(d.columns):
            add_ci_band(ax, x, d, lo_col, hi_col, args.ci_alpha)

    if meta["zero_line"]:
        ax.axhline(0.0, linewidth=1.0)

    ax.set_xlabel("Initial ln frequency")
    ax.set_ylabel(meta["ylabel"])
    maybe_set_title(ax, rf"$\Delta t = {dt}$", args)
    style_axes(ax, args.ylabel_x)
    return handles, labels


def make_binned_multidt_figure(
    metric: str,
    dt_values: List[int],
    binned_data: Dict[str, pd.DataFrame],
    out_dir: Path,
    args,
) -> None:
    ncols = len(dt_values)
    fig_w = args.binned_panel_w * ncols
    fig_h = args.binned_panel_h

    fig, axes = plt.subplots(
        nrows=1,
        ncols=ncols,
        figsize=(fig_w, fig_h),
        squeeze=False,
    )
    axes = axes[0]

    all_handles = []
    all_labels = []

    for i, dt in enumerate(dt_values):
        handles, labels = plot_binned_metric_panel(axes[i], metric, dt, binned_data, args)
        if i == 0:
            all_handles = handles
            all_labels = labels

    has_suptitle = bool(args.suptitle)
    if has_suptitle:
        fig.suptitle(args.suptitle)

    place_shared_legend(fig, all_handles, all_labels, args)
    apply_fixed_layout(fig, args, width_in=fig_w, height_in=fig_h, has_suptitle=has_suptitle)
    save_figure(fig, out_dir / VALID_BINNED_METRICS[metric]["fname"], args.fmt, args.dpi)
    plt.close(fig)


# =========================================================
# Main
# =========================================================
def main() -> None:
    args = build_argparser().parse_args()
    set_matplotlib_defaults(args.font, args.font_size)

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    conditions = parse_conditions(args.cond)
    dt_values = parse_dt_values(args.dt_values)

    metrics = parse_list_str(args.metrics)
    bad_metrics = [m for m in metrics if m not in VALID_POOLED_METRICS]
    if bad_metrics:
        raise ValueError(f"Invalid pooled metrics: {bad_metrics}")

    binned_metrics = parse_list_str(args.binned_metrics)
    bad_binned = [m for m in binned_metrics if m not in VALID_BINNED_METRICS]
    if bad_binned:
        raise ValueError(f"Invalid binned metrics: {bad_binned}")

    summaries: Dict[str, pd.DataFrame] = {}
    binned_data: Dict[str, pd.DataFrame] = {}

    for label, cdir in conditions:
        summaries[label] = load_summary(
            cdir=cdir,
            csv_sep=args.csv_sep,
            dt_values=dt_values,
            metrics=metrics,
        )

        if args.also_plot_binned:
            needed_binned_cols = list({m for m in binned_metrics if m != "D_mad"} | {"D_mad"})
            binned_data[label] = load_binned(
                cdir=cdir,
                csv_sep=args.csv_sep,
                dt_values=dt_values,
                binned_metrics=needed_binned_cols,
            )

        if args.debug:
            print(f"[DEBUG] {label}")
            print(f"  dir: {cdir}")
            print(f"  summary rows: {len(summaries[label])}")
            print(f"  summary cols: {list(summaries[label].columns)}")
            if args.also_plot_binned:
                print(f"  binned rows: {len(binned_data[label])}")
                print(f"  binned cols: {list(binned_data[label].columns)}")

    make_pooled_multipanel(metrics, summaries, out_dir, args)
    make_pooled_single_figures(metrics, summaries, out_dir, args)

    if args.also_plot_binned:
        dt_plot = dt_values if dt_values else sorted(
            set(np.concatenate([df["dt"].to_numpy(int) for df in binned_data.values()]))
        )
        for metric in binned_metrics:
            make_binned_multidt_figure(metric, dt_plot, binned_data, out_dir, args)

    print(f"[SAVED] {(out_dir / f'compare_finite_time_pooled_multipanel.{args.fmt}').resolve()}")
    for metric in metrics:
        print(f"[SAVED] {(out_dir / (VALID_POOLED_METRICS[metric]['fname'] + '.' + args.fmt)).resolve()}")
    if args.also_plot_binned:
        for metric in binned_metrics:
            print(f"[SAVED] {(out_dir / (VALID_BINNED_METRICS[metric]['fname'] + '.' + args.fmt)).resolve()}")


if __name__ == "__main__":
    main()