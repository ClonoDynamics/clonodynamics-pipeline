#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_compare_dt1.py

Compare dt=1 model-agnostic diagnostics across multiple analysis conditions.

For each condition, the script reads:
  - drift_by_dt_xbins.csv   -> median_dx vs x_center, ppos vs x_center
  - msd_by_dt_xbins.csv     -> msd vs x_center

Then it filters the requested dt (default: 1) and creates three plots:
  1) median_dx vs x_center
  2) msd vs x_center
  3) P(Δx > 0) vs x_center

Main features
-------------
- stable panel geometry for multipanel assembly
- optional CI ribbons
- optional legend outside
- optional removal of titles
- adaptive margins to avoid clipping when figures are small
- ylabel spacing handled with labelpad, not fixed coordinates
- outside legends anchored in figure coordinates

Examples
python ./code/plotting/plot_compare_dt1.py \
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


# -------------------------
# CLI
# -------------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Compare median_dx, MSD, and P(Δx>0) across multiple conditions."
    )

    ap.add_argument(
        "--cond",
        action="append",
        required=True,
        help='Condition spec: "LABEL=PATH". Repeat as needed.',
    )

    ap.add_argument("--dt", type=int, default=1, help="dt value to plot. Default: 1")
    ap.add_argument("--out_dir", required=True, help="Output directory")

    ap.add_argument("--csv_sep", default=",", help="CSV separator. Default: ','")
    ap.add_argument("--fmt", choices=["png", "pdf"], default="png", help="Output format")
    ap.add_argument("--dpi", type=int, default=300, help="DPI for PNG")

    ap.add_argument("--font", default="Arial", help="Matplotlib font family")
    ap.add_argument("--font_size", type=float, default=10.0, help="Base font size")

    ap.add_argument("--fig_w", type=float, default=5.0, help="Figure width in inches")
    ap.add_argument("--fig_h", type=float, default=3.5, help="Figure height in inches")

    ap.add_argument("--line_width", type=float, default=1.5, help="Line width")
    ap.add_argument("--marker", default="o", help="Marker style")
    ap.add_argument("--marker_size", type=float, default=3.5, help="Marker size")

    ap.add_argument("--show_ci", action="store_true", help="Draw CI ribbons if available")
    ap.add_argument("--ci_alpha", type=float, default=0.12, help="Alpha for CI ribbons")

    ap.add_argument("--title_prefix", default="", help="Optional prefix for titles")
    ap.add_argument("--no_title", action="store_true", help="Remove plot titles")

    ap.add_argument("--legend_out", action="store_true", help="Place legend outside")
    ap.add_argument("--legend_ncol", type=int, default=1, help="Legend columns")
    ap.add_argument("--legend_x", type=float, default=0.70, help="Legend x anchor in figure coordinates if outside")
    ap.add_argument("--legend_y", type=float, default=0.50, help="Legend y anchor in figure coordinates if outside")

    # Geometry controls
    ap.add_argument("--left_margin", type=float, default=0.16, help="Base left margin")
    ap.add_argument("--right_margin_in", type=float, default=0.97, help="Right margin if legend is inside")
    ap.add_argument("--right_margin_out", type=float, default=0.68, help="Right margin if legend is outside")
    ap.add_argument("--bottom_margin", type=float, default=0.17, help="Base bottom margin")
    ap.add_argument("--top_margin", type=float, default=0.90, help="Base top margin")

    # Label/title spacing
    ap.add_argument("--xlabel_pad", type=float, default=4.0, help="Padding for x-axis label")
    ap.add_argument("--ylabel_pad", type=float, default=8.0, help="Padding for y-axis label")
    ap.add_argument("--title_pad", type=float, default=8.0, help="Padding for title")

    ap.add_argument("--debug", action="store_true", help="Print debug info")

    return ap


# -------------------------
# Global style
# -------------------------
def set_matplotlib_defaults(font: str, font_size: float, title_pad: float) -> None:
    plt.rcParams["font.family"] = font
    plt.rcParams["font.size"] = font_size
    plt.rcParams["axes.titlesize"] = font_size
    plt.rcParams["axes.labelsize"] = font_size
    plt.rcParams["legend.fontsize"] = font_size
    plt.rcParams["xtick.labelsize"] = font_size
    plt.rcParams["ytick.labelsize"] = font_size

    plt.rcParams["axes.titlepad"] = title_pad
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["figure.facecolor"] = "white"


# -------------------------
# Helpers
# -------------------------
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


def filter_dt(df: pd.DataFrame, dt: int) -> pd.DataFrame:
    d = df.copy()
    d["dt"] = pd.to_numeric(d["dt"], errors="coerce")
    d = d.dropna(subset=["dt"]).copy()
    d["dt"] = d["dt"].astype(int)
    d = d[d["dt"] == int(dt)].copy()
    return d.sort_values("x_center").reset_index(drop=True)


def style_axes(ax: plt.Axes) -> None:
    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(True)
    ax.tick_params(direction="out")


def set_axis_labels(ax: plt.Axes, xlabel: str, ylabel: str, args) -> None:
    ax.set_xlabel(xlabel, labelpad=args.xlabel_pad)
    ax.set_ylabel(ylabel, labelpad=args.ylabel_pad)


def adaptive_margins(args) -> Tuple[float, float, float, float]:
    """
    Conservative adaptive margins to avoid clipping of labels, titles, and legends.
    Returns: left, right, bottom, top
    """
    # Right margin
    if args.legend_out:
        if args.fig_w < 5.0:
            right = min(args.right_margin_out, 0.64)
        elif args.fig_w < 7.0:
            right = min(max(args.right_margin_out, 0.68), 0.72)
        else:
            right = min(max(args.right_margin_out, 0.72), 0.78)
    else:
        right = args.right_margin_in

    # Left margin: generous because ylabel is no longer manually forced
    if args.fig_w < 4.0:
        left = max(args.left_margin, 0.24)
    elif args.fig_w < 5.0:
        left = max(args.left_margin, 0.19)
    else:
        left = max(args.left_margin, 0.16)

    # Bottom margin
    if args.fig_h < 3.0:
        bottom = max(args.bottom_margin, 0.24)
    elif args.fig_h < 3.5:
        bottom = max(args.bottom_margin, 0.20)
    else:
        bottom = max(args.bottom_margin, 0.17)

    # Top margin
    if args.fig_h < 3.0:
        top = min(0.82, args.top_margin)
    elif args.fig_h < 3.5:
        top = min(0.85, args.top_margin)
    else:
        top = min(0.88, args.top_margin)

    if not args.no_title:
        top -= 0.04

    top = max(top, bottom + 0.18)

    return left, right, bottom, top


def apply_fixed_layout(fig: plt.Figure, args) -> None:
    left, right, bottom, top = adaptive_margins(args)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)


def save_figure(fig: plt.Figure, out_base: Path, fmt: str, dpi: int) -> None:
    # Deliberately avoid bbox_inches="tight" to keep panel geometry stable
    if fmt == "png":
        fig.savefig(out_base.with_suffix(".png"), dpi=dpi)
    else:
        fig.savefig(out_base.with_suffix(".pdf"))


def place_legend(
    fig: plt.Figure,
    ax: plt.Axes,
    legend_out: bool,
    ncol: int,
    legend_x: float,
    legend_y: float
) -> None:
    if legend_out:
        fig.legend(
            *ax.get_legend_handles_labels(),
            frameon=False,
            loc="center left",
            bbox_to_anchor=(legend_x, legend_y),
            bbox_transform=fig.transFigure,
            borderaxespad=0.0,
            ncol=max(1, int(ncol)),
            handlelength=2.2,
        )
    else:
        ax.legend(frameon=False)


def maybe_set_title(ax: plt.Axes, title: str, args) -> None:
    if not args.no_title:
        final_title = f"{args.title_prefix}{title}" if args.title_prefix else title
        ax.set_title(final_title, pad=args.title_pad)


def add_ci_band(ax: plt.Axes, x: np.ndarray, df: pd.DataFrame, lo_col: str, hi_col: str, alpha: float) -> None:
    if {lo_col, hi_col}.issubset(df.columns):
        lo = pd.to_numeric(df[lo_col], errors="coerce").to_numpy(float)
        hi = pd.to_numeric(df[hi_col], errors="coerce").to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(lo) & np.isfinite(hi)
        if ok.any():
            ax.fill_between(x[ok], lo[ok], hi[ok], alpha=alpha)


# -------------------------
# Main plotting
# -------------------------
def main() -> None:
    args = build_argparser().parse_args()
    set_matplotlib_defaults(args.font, args.font_size, args.title_pad)

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    conditions = parse_conditions(args.cond)

    drift_data: Dict[str, pd.DataFrame] = {}
    msd_data: Dict[str, pd.DataFrame] = {}

    for label, cdir in conditions:
        drift_path = cdir / "drift_by_dt_xbins.csv"
        msd_path = cdir / "msd_by_dt_xbins.csv"

        drift = read_required_csv(drift_path, args.csv_sep)
        msd = read_required_csv(msd_path, args.csv_sep)

        ensure_columns(drift, ["dt", "x_center", "median_dx", "ppos"], drift_path)
        ensure_columns(msd, ["dt", "x_center", "msd"], msd_path)

        drift = filter_dt(drift, args.dt)
        msd = filter_dt(msd, args.dt)

        if drift.empty:
            raise ValueError(f"No rows with dt={args.dt} in {drift_path}")
        if msd.empty:
            raise ValueError(f"No rows with dt={args.dt} in {msd_path}")

        drift_data[label] = drift
        msd_data[label] = msd

        if args.debug:
            print(f"[DEBUG] {label}")
            print(f"  dir: {cdir}")
            print(f"  drift rows dt={args.dt}: {len(drift)}")
            print(f"  msd rows dt={args.dt}:   {len(msd)}")

    # ---------------------------------
    # Plot 1: median_dx vs x_center
    # ---------------------------------
    fig1, ax1 = plt.subplots(figsize=(args.fig_w, args.fig_h))

    for label, d in drift_data.items():
        x = pd.to_numeric(d["x_center"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(d["median_dx"], errors="coerce").to_numpy(float)

        ok = np.isfinite(x) & np.isfinite(y)
        ax1.plot(
            x[ok],
            y[ok],
            lw=args.line_width,
            marker=args.marker,
            ms=args.marker_size,
            label=label,
        )

        if args.show_ci:
            add_ci_band(ax1, x, d, "median_lo", "median_hi", args.ci_alpha)

    ax1.axhline(0.0, linewidth=1.0)
    set_axis_labels(ax1, "Initial ln frequency", "median(Δx)", args)
    maybe_set_title(ax1, f"Median displacement by initial frequency (Δt={args.dt})", args)
    style_axes(ax1)
    place_legend(fig1, ax1, args.legend_out, args.legend_ncol, args.legend_x, args.legend_y)
    apply_fixed_layout(fig1, args)
    save_figure(fig1, out_dir / f"compare_median_dx_dt{args.dt}", args.fmt, args.dpi)
    plt.close(fig1)

    # ---------------------------------
    # Plot 2: MSD vs x_center
    # ---------------------------------
    fig2, ax2 = plt.subplots(figsize=(args.fig_w, args.fig_h))

    for label, d in msd_data.items():
        x = pd.to_numeric(d["x_center"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(d["msd"], errors="coerce").to_numpy(float)

        ok = np.isfinite(x) & np.isfinite(y)
        ax2.plot(
            x[ok],
            y[ok],
            lw=args.line_width,
            marker=args.marker,
            ms=args.marker_size,
            label=label,
        )

        if args.show_ci:
            add_ci_band(ax2, x, d, "msd_lo", "msd_hi", args.ci_alpha)

    set_axis_labels(ax2, "Initial ln frequency", "MSD = E[(Δx)^2]", args)
    maybe_set_title(ax2, f"MSD by initial frequency (Δt={args.dt})", args)
    style_axes(ax2)
    place_legend(fig2, ax2, args.legend_out, args.legend_ncol, args.legend_x, args.legend_y)
    apply_fixed_layout(fig2, args)
    save_figure(fig2, out_dir / f"compare_msd_dt{args.dt}", args.fmt, args.dpi)
    plt.close(fig2)

    # ---------------------------------
    # Plot 3: P(Δx > 0) vs x_center
    # ---------------------------------
    fig3, ax3 = plt.subplots(figsize=(args.fig_w, args.fig_h))

    for label, d in drift_data.items():
        x = pd.to_numeric(d["x_center"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(d["ppos"], errors="coerce").to_numpy(float)

        ok = np.isfinite(x) & np.isfinite(y)
        ax3.plot(
            x[ok],
            y[ok],
            lw=args.line_width,
            marker=args.marker,
            ms=args.marker_size,
            label=label,
        )

        if args.show_ci:
            add_ci_band(ax3, x, d, "ppos_lo", "ppos_hi", args.ci_alpha)

    ax3.axhline(0.5, linewidth=1.0)
    set_axis_labels(ax3, "Initial ln frequency", "P(Δx > 0)", args)
    maybe_set_title(ax3, f"Probability of positive displacement by initial frequency (Δt={args.dt})", args)
    style_axes(ax3)
    place_legend(fig3, ax3, args.legend_out, args.legend_ncol, args.legend_x, args.legend_y)
    apply_fixed_layout(fig3, args)
    save_figure(fig3, out_dir / f"compare_ppos_dt{args.dt}", args.fmt, args.dpi)
    plt.close(fig3)

    print(f"[SAVED] {out_dir / f'compare_median_dx_dt{args.dt}.{args.fmt}'}")
    print(f"[SAVED] {out_dir / f'compare_msd_dt{args.dt}.{args.fmt}'}")
    print(f"[SAVED] {out_dir / f'compare_ppos_dt{args.dt}.{args.fmt}'}")


if __name__ == "__main__":
    main()