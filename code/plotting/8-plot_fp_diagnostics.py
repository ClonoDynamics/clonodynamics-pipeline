#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
8-plot_fp_diagnostics.py
ROBUST v3.1 - plot FP diagnostic curves from fp_<tag>.grid.csv

PURPOSE
-------
Plot FP diagnostic quantities produced by 8-fp_diagnostics_build.py.

This version adds:
1) explicit global title control
   - if --title is provided, it is shown on every selected plot
   - if --title is omitted, no title is shown

2) explicit subtitle control
   - if --subtitle_mode is provided, each plot shows its predefined subtitle
   - if --subtitle_mode is omitted, no subtitle is shown

3) explicit legend control
   - by default, legends are not shown
   - --legend_all shows legends on all generated plots
   - --legend_plot shows legends only on selected plot types

4) per-plot figure size / dpi controls
   - global defaults can still be provided
   - each plot can override width, height, and dpi independently

SUPPORTED FIGURES
-----------------
- D(x)
- S(x)
- J(x)
- p_emp vs p_closed
- p_emp / p_closed

BACKWARD COMPATIBILITY
----------------------
Older grids are still supported:
    D      -> D_final
    J_emp  -> J_smooth
    S_emp  -> S_smooth

EXAMPLES
--------
Show all plots, no title, no subtitles, no legends:
python plot_fp_diagnostics.py \
  --grid_csv ./fp.grid.csv \
  --out_dir ./figures \
  --plot_all

Show a global title and built-in subtitles on all plots:
python plot_diagnostics.py \
  --grid_csv ./fp.grid.csv \
  --out_dir ./figures \
  --plot_all \
  --title "FP diagnostics — hinge, TT, dt=1, p=0.02" \
  --subtitle_mode

Show legends only on D and density plots:
python plot_diagnostics.py \
  --grid_csv ./fp.grid.csv \
  --out_dir ./figures \
  --plot_all \
  --legend_plot D density

Set different sizes for individual plots:
python plot_diagnostics.py \
  --grid_csv ./fp.grid.csv \
  --out_dir ./figures \
  --plot_all \
  --fig_w 5.5 --fig_h 3.2 --dpi 300 \
  --fig_w_D 5.0 --fig_h_D 3.0 --dpi_D 400 \
  --fig_w_density 6.0 --fig_h_density 3.5 --dpi_density 300

  python ./code/plotting/plot_fp_diagnostics_.py \
--grid_csv ./results/8-fp_stationary/TT_p02_w_clone/fp.grid.csv \
--out_dir ./figures/8-fp_stationary/TT_p02_w_clone \
--plot_all \
--title "TT, p=0.02" \
--fig_w 4 --fig_h 3 --dpi 300 \
--annotate_metrics \
--legend_all \
--summary_json ./results/8-fp_stationary/TT_p02_w_clone/fp.summary.json \
--no_grid
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# IO / style helpers
# ============================================================

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def apply_rcparams(font: str, font_size: float):
    plt.rcParams.update({
        "font.family": font,
        "font.size": font_size,
        "axes.titlesize": font_size,
        "axes.labelsize": font_size,
        "legend.fontsize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def parse_lim(s: str):
    if not s:
        return None
    a, b = s.split(",")
    return float(a), float(b)


def save_fig(fig, outpath: Path, dpi: int, also_pdf: bool):
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    if also_pdf:
        fig.savefig(outpath.with_suffix(".pdf"))


# ============================================================
# Grid compatibility layer
# ============================================================

def get_grid_columns(grid: pd.DataFrame) -> dict:
    cols = set(grid.columns)

    out = {
        "x": "x" if "x" in cols else None,
        "p_emp": "p_emp" if "p_emp" in cols else None,
        "p_closed": "p_closed" if "p_closed" in cols else None,
        "ratio": "ratio_pemp_over_pclosed" if "ratio_pemp_over_pclosed" in cols else None,
        "D_raw_interp": "D_raw_interp" if "D_raw_interp" in cols else None,
        "D_smooth": "D_smooth" if "D_smooth" in cols else None,
        "D_final": "D_final" if "D_final" in cols else ("D" if "D" in cols else None),
        "J_raw": "J_raw" if "J_raw" in cols else None,
        "J_smooth": "J_smooth" if "J_smooth" in cols else ("J_emp" if "J_emp" in cols else None),
        "S_raw": "S_raw" if "S_raw" in cols else None,
        "S_smooth": "S_smooth" if "S_smooth" in cols else ("S_emp" if "S_emp" in cols else None),
    }

    required = ["x", "p_emp", "p_closed", "ratio"]
    missing = [k for k in required if out[k] is None]
    if missing:
        raise ValueError(
            f"Missing required grid fields: {missing}. Available columns: {list(grid.columns)}"
        )

    return out


# ============================================================
# Plot helpers
# ============================================================

def should_show_legend(plot_key: str, args) -> bool:
    if args.legend_all:
        return True
    return plot_key in set(args.legend_plot or [])


def maybe_add_legend(ax, args, show: bool):
    if not show:
        return

    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return

    if args.legend_mode == "inside":
        ax.legend(loc=args.legend_loc, frameon=args.legend_frameon)
        return

    if args.legend_mode == "outside":
        ax.legend(
            loc=args.legend_loc,
            bbox_to_anchor=(args.legend_bbox_x, args.legend_bbox_y),
            frameon=args.legend_frameon,
            borderaxespad=args.legend_borderaxespad,
        )
        return

    raise ValueError(f"Unknown --legend_mode {args.legend_mode}")


def set_title_and_subtitle(ax, args, subtitle: str):
    pieces = []
    if args.title:
        pieces.append(args.title)
    if args.subtitle_mode:
        pieces.append(subtitle)

    if pieces:
        ax.set_title("\n".join(pieces))


def finalize_axes(ax, args, xlim=None, ylims=None, default_ylim=None):
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylims is not None:
        ax.set_ylim(*ylims)
    elif default_ylim is not None:
        ax.set_ylim(*default_ylim)

    if not args.no_grid:
        ax.grid(True, linewidth=0.5, alpha=0.4)


def get_plot_dims(args, key: str):
    fig_w = getattr(args, f"fig_w_{key}")
    fig_h = getattr(args, f"fig_h_{key}")
    dpi = getattr(args, f"dpi_{key}")

    if fig_w is None:
        fig_w = args.fig_w
    if fig_h is None:
        fig_h = args.fig_h
    if dpi is None:
        dpi = args.dpi

    return fig_w, fig_h, dpi


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(description="Plot FP diagnostic curves from fp_<tag>.grid.csv.")

    ap.add_argument("--grid_csv", required=True, help="fp_<tag>.grid.csv from build script")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--prefix", default="", help="Optional prefix for output filenames")

    # global aesthetics
    ap.add_argument("--font", default="Arial")
    ap.add_argument("--font_size", type=float, default=10.0)
    ap.add_argument("--fig_w", type=float, default=6.0, help="Default figure width in inches")
    ap.add_argument("--fig_h", type=float, default=4.0, help="Default figure height in inches")
    ap.add_argument("--dpi", type=int, default=300, help="Default figure dpi")
    ap.add_argument("--lw", type=float, default=1.8)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--also_pdf", action="store_true")

    # per-plot size / dpi overrides
    for key in ["D", "S", "J", "density", "ratio"]:
        ap.add_argument(f"--fig_w_{key}", type=float, default=None, help=f"Width in inches for {key} plot")
        ap.add_argument(f"--fig_h_{key}", type=float, default=None, help=f"Height in inches for {key} plot")
        ap.add_argument(f"--dpi_{key}", type=int, default=None, help=f"DPI for {key} plot")

    # legend controls
    ap.add_argument(
        "--legend_mode",
        default="inside",
        choices=["inside", "outside"],
        help="Placement mode for legends when legends are enabled",
    )
    ap.add_argument("--legend_loc", default="best", help="matplotlib legend loc")
    ap.add_argument("--legend_bbox_x", type=float, default=1.02)
    ap.add_argument("--legend_bbox_y", type=float, default=1.0)
    ap.add_argument("--legend_frameon", action="store_true")
    ap.add_argument("--legend_borderaxespad", type=float, default=0.0)
    ap.add_argument("--legend_all", action="store_true", help="Show legends on all generated plots")
    ap.add_argument(
        "--legend_plot",
        nargs="+",
        choices=["D", "S", "J", "density", "ratio"],
        default=[],
        help=(
            "Show legends only on the specified plots. "
            "Example: --legend_plot D density ratio"
        ),
    )

    # titles / axes
    ap.add_argument("--title", default="", help="Global title shown on all plots only if provided")
    ap.add_argument(
        "--subtitle_mode",
        action="store_true",
        help="Show predefined plot-specific subtitles on all generated plots",
    )
    ap.add_argument("--xlim", default="", help="x limits as xmin,xmax")
    ap.add_argument("--ylim", default="", help="generic y limits as ymin,ymax for linear plots")
    ap.add_argument("--ylim_D", default="", help="y limits for D plot")
    ap.add_argument("--ylim_J", default="", help="y limits for J plot")
    ap.add_argument("--ylim_S", default="", help="y limits for S plot")
    ap.add_argument("--ylim_ratio", default="", help="y limits for ratio plot")
    ap.add_argument("--ylog_min", type=float, default=1e-8, help="Minimum y for semilogy density plots")
    ap.add_argument("--no_grid", action="store_true")

    # plot selection
    ap.add_argument("--plot_D", action="store_true", help="Plot D_raw_interp, D_smooth, D_final")
    ap.add_argument("--plot_S", action="store_true", help="Plot source S(x)")
    ap.add_argument("--plot_J", action="store_true", help="Plot flux J(x)")
    ap.add_argument("--plot_density", action="store_true", help="Plot p_emp vs p_closed")
    ap.add_argument("--plot_ratio", action="store_true", help="Plot p_emp/p_closed ratio")
    ap.add_argument("--plot_all", action="store_true", help="Make all figures")

    # raw/smoothed controls
    ap.add_argument("--show_raw_D", action="store_true", help="Show D_raw_interp when available")
    ap.add_argument("--show_smooth_D", action="store_true", help="Show D_smooth when available")
    ap.add_argument("--show_final_D", action="store_true", help="Show D_final")

    ap.add_argument("--show_raw_J", action="store_true", help="Show J_raw when available")
    ap.add_argument("--show_smooth_J", action="store_true", help="Show J_smooth")

    ap.add_argument("--show_raw_S", action="store_true", help="Show S_raw when available")
    ap.add_argument("--show_smooth_S", action="store_true", help="Show S_smooth")

    # optional summary annotation
    ap.add_argument("--summary_json", default="", help="Optional fp_<tag>.summary.json for JS/KS annotation")
    ap.add_argument("--annotate_metrics", action="store_true", help="Annotate JS/KS on density plot")

    args = ap.parse_args()

    out = Path(args.out_dir)
    ensure_dir(out)
    apply_rcparams(args.font, args.font_size)

    grid = pd.read_csv(args.grid_csv)
    cmap = get_grid_columns(grid)
    x = pd.to_numeric(grid[cmap["x"]], errors="coerce").to_numpy(float)

    xlim = parse_lim(args.xlim)
    ylim = parse_lim(args.ylim)
    ylim_D = parse_lim(args.ylim_D)
    ylim_J = parse_lim(args.ylim_J)
    ylim_S = parse_lim(args.ylim_S)
    ylim_ratio = parse_lim(args.ylim_ratio)

    metrics = None
    if args.summary_json:
        with open(args.summary_json, "r", encoding="utf-8") as f:
            metrics = json.load(f)

    stem = Path(args.grid_csv).stem

    any_flag = args.plot_D or args.plot_S or args.plot_J or args.plot_density or args.plot_ratio or args.plot_all
    do_all = args.plot_all or (not any_flag)

    show_raw_D = args.show_raw_D or do_all
    show_smooth_D = args.show_smooth_D or do_all
    show_final_D = args.show_final_D or do_all

    show_raw_J = args.show_raw_J
    show_smooth_J = args.show_smooth_J or do_all
    if not args.show_raw_J and not args.show_smooth_J and not do_all:
        show_smooth_J = True

    show_raw_S = args.show_raw_S
    show_smooth_S = args.show_smooth_S or do_all
    if not args.show_raw_S and not args.show_smooth_S and not do_all:
        show_smooth_S = True

    subtitles = {
        "D": "Diffusion function",
        "S": "Source / sink diagnostic",
        "J": "Probability current",
        "density": "Empirical vs closed stationary density",
        "ratio": "Density ratio",
    }

    # --------------------------------------------------------
    # D(x)
    # --------------------------------------------------------
    if do_all or args.plot_D:
        fig_w, fig_h, dpi = get_plot_dims(args, "D")
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
        ax = fig.add_subplot(111)

        plotted = False

        if show_raw_D and cmap["D_raw_interp"] is not None:
            y = pd.to_numeric(grid[cmap["D_raw_interp"]], errors="coerce").to_numpy(float)
            ax.plot(x, y, linewidth=args.lw, alpha=args.alpha, label="D raw interp")
            plotted = True

        if show_smooth_D and cmap["D_smooth"] is not None:
            y = pd.to_numeric(grid[cmap["D_smooth"]], errors="coerce").to_numpy(float)
            ax.plot(x, y, linewidth=args.lw, alpha=args.alpha, label="D smooth")
            plotted = True

        if show_final_D and cmap["D_final"] is not None:
            y = pd.to_numeric(grid[cmap["D_final"]], errors="coerce").to_numpy(float)
            ax.plot(x, y, linewidth=args.lw, alpha=args.alpha, label="D final")
            plotted = True

        if not plotted:
            ax.text(0.5, 0.5, "No D columns available to plot", transform=ax.transAxes, ha="center", va="center")

        ax.set_xlabel("x = ln f")
        ax.set_ylabel("D(x)")
        set_title_and_subtitle(ax, args, subtitles["D"])
        maybe_add_legend(ax, args, should_show_legend("D", args))
        finalize_axes(ax, args, xlim=xlim, ylims=ylim_D, default_ylim=ylim)

        outpath = out / f"{args.prefix}{stem}.diffusion_D.png"
        save_fig(fig, outpath, dpi=dpi, also_pdf=args.also_pdf)
        plt.close(fig)

    # --------------------------------------------------------
    # S(x)
    # --------------------------------------------------------
    if do_all or args.plot_S:
        fig_w, fig_h, dpi = get_plot_dims(args, "S")
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
        ax = fig.add_subplot(111)

        plotted = False
        ax.axhline(0, linewidth=1)

        if show_raw_S and cmap["S_raw"] is not None:
            y = pd.to_numeric(grid[cmap["S_raw"]], errors="coerce").to_numpy(float)
            ax.plot(x, y, linewidth=args.lw, alpha=args.alpha, label="S raw")
            plotted = True

        if show_smooth_S and cmap["S_smooth"] is not None:
            y = pd.to_numeric(grid[cmap["S_smooth"]], errors="coerce").to_numpy(float)
            ax.plot(x, y, linewidth=args.lw, alpha=args.alpha, label="S smooth")
            plotted = True

        if not plotted:
            ax.text(0.5, 0.5, "No S columns available to plot", transform=ax.transAxes, ha="center", va="center")

        ax.set_xlabel("x = ln f")
        ax.set_ylabel("S(x) = dJ/dx")
        set_title_and_subtitle(ax, args, subtitles["S"])
        maybe_add_legend(ax, args, should_show_legend("S", args))
        finalize_axes(ax, args, xlim=xlim, ylims=ylim_S, default_ylim=ylim)

        outpath = out / f"{args.prefix}{stem}.source_S.png"
        save_fig(fig, outpath, dpi=dpi, also_pdf=args.also_pdf)
        plt.close(fig)

    # --------------------------------------------------------
    # J(x)
    # --------------------------------------------------------
    if do_all or args.plot_J:
        fig_w, fig_h, dpi = get_plot_dims(args, "J")
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
        ax = fig.add_subplot(111)

        plotted = False
        ax.axhline(0, linewidth=1)

        if show_raw_J and cmap["J_raw"] is not None:
            y = pd.to_numeric(grid[cmap["J_raw"]], errors="coerce").to_numpy(float)
            ax.plot(x, y, linewidth=args.lw, alpha=args.alpha, label="J raw")
            plotted = True

        if show_smooth_J and cmap["J_smooth"] is not None:
            y = pd.to_numeric(grid[cmap["J_smooth"]], errors="coerce").to_numpy(float)
            ax.plot(x, y, linewidth=args.lw, alpha=args.alpha, label="J smooth")
            plotted = True

        if not plotted:
            ax.text(0.5, 0.5, "No J columns available to plot", transform=ax.transAxes, ha="center", va="center")

        ax.set_xlabel("x = ln f")
        ax.set_ylabel("J(x)")
        set_title_and_subtitle(ax, args, subtitles["J"])
        maybe_add_legend(ax, args, should_show_legend("J", args))
        finalize_axes(ax, args, xlim=xlim, ylims=ylim_J, default_ylim=ylim)

        outpath = out / f"{args.prefix}{stem}.flux_J.png"
        save_fig(fig, outpath, dpi=dpi, also_pdf=args.also_pdf)
        plt.close(fig)

    # --------------------------------------------------------
    # p_emp vs p_closed
    # --------------------------------------------------------
    if do_all or args.plot_density:
        fig_w, fig_h, dpi = get_plot_dims(args, "density")
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
        ax = fig.add_subplot(111)

        p_emp = np.clip(pd.to_numeric(grid[cmap["p_emp"]], errors="coerce").to_numpy(float), args.ylog_min, np.inf)
        p_closed = np.clip(pd.to_numeric(grid[cmap["p_closed"]], errors="coerce").to_numpy(float), args.ylog_min, np.inf)

        ax.semilogy(x, p_emp, linewidth=args.lw, alpha=args.alpha, label="p_emp(x)")
        ax.semilogy(x, p_closed, linewidth=args.lw, alpha=args.alpha, label="p_closed(x)")
        ax.set_xlabel("x = ln f")
        ax.set_ylabel("density (log scale)")
        set_title_and_subtitle(ax, args, subtitles["density"])

        if args.annotate_metrics and metrics is not None:
            js = metrics.get("JS", None)
            ks = metrics.get("KS", None)
            txt = []
            if js is not None:
                txt.append(f"JS={js:.3f}")
            if ks is not None:
                txt.append(f"KS={ks:.3f}")
            if txt:
                ax.text(0.55, 0.9, "  ".join(txt), transform=ax.transAxes, ha="left", va="top")

        maybe_add_legend(ax, args, should_show_legend("density", args))
        finalize_axes(ax, args, xlim=xlim)

        outpath = out / f"{args.prefix}{stem}.p_emp_vs_p_closed.png"
        save_fig(fig, outpath, dpi=dpi, also_pdf=args.also_pdf)
        plt.close(fig)

    # --------------------------------------------------------
    # Ratio
    # --------------------------------------------------------
    if do_all or args.plot_ratio:
        fig_w, fig_h, dpi = get_plot_dims(args, "ratio")
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
        ax = fig.add_subplot(111)

        r = pd.to_numeric(grid[cmap["ratio"]], errors="coerce").to_numpy(float)
        ax.axhline(1, linewidth=1)
        ax.plot(x, r, linewidth=args.lw, alpha=args.alpha, label="p_emp / p_closed")
        ax.set_xlabel("x = ln f")
        ax.set_ylabel("p_emp / p_closed")
        set_title_and_subtitle(ax, args, subtitles["ratio"])
        maybe_add_legend(ax, args, should_show_legend("ratio", args))
        finalize_axes(ax, args, xlim=xlim, ylims=ylim_ratio, default_ylim=ylim)

        outpath = out / f"{args.prefix}{stem}.ratio.png"
        save_fig(fig, outpath, dpi=dpi, also_pdf=args.also_pdf)
        plt.close(fig)

    print(f"[DONE] Plots saved to: {out}")


if __name__ == "__main__":
    main()
