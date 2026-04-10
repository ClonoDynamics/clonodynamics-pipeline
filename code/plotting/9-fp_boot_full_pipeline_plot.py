#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
9-fp_boot_full_pipeline_plot_rewritten.py

Plot-only script for full-pipeline bootstrap outputs from:
  09_bootstrap_full_pipeline_flux.py

This rewritten version aligns the plotting interface with the updated
8-fp_diagnostics_plot workflow.

Main behavior
-------------
1) Global title
   - If --title is provided, the same global title is shown on every plot.
   - If --title is omitted, no title is shown.

2) Subtitles
   - If --subtitle_mode is provided, each plot shows a built-in subtitle:
       drift     -> Bootstrap drift estimate
       diffusion -> Bootstrap diffusion estimate
       J         -> Bootstrap flux estimate
       S         -> Bootstrap source estimate
       tau       -> Bootstrap distribution of stabilization time
   - If --subtitle_mode is omitted, no subtitle is shown.

3) Legends
   - By default, legends are not shown.
   - --legend_all shows legends on all generated plots.
   - --legend_plot allows legends only on selected plots.
     Example:
       --legend_plot drift diffusion tau
   - Allowed values:
       drift diffusion J S tau

4) Figure size / dpi
   - Global defaults:
       --fig_w --fig_h --dpi
   - Per-plot overrides:
       --fig_w_drift --fig_h_drift --dpi_drift
       --fig_w_diffusion --fig_h_diffusion --dpi_diffusion
       --fig_w_J --fig_h_J --dpi_J
       --fig_w_S --fig_h_S --dpi_S
       --fig_w_tau --fig_h_tau --dpi_tau

5) Plot selection
   - If no --plot_* flag is passed, all available inputs are plotted.
   - You can also request only specific plots with:
       --plot_drift --plot_diffusion --plot_J --plot_S --plot_tau
       or --plot_all

Typical usage
-------------
python ./code/plotting/9-fp_boot_full_pipeline_plot_rewritten.py \
  --drift_summary ./results/...drift_summary.csv \
  --diffusion_summary ./results/...diffusion_summary.csv \
  --J_summary ./results/...J_boot_summary.csv \
  --S_summary ./results/...S_boot_summary.csv \
  --drift_params ./results/...drift_boot_params.csv \
  --out_dir ./figures/9-fp_boot/hinge \
  --plot_all \
  --title "FP bootstrap diagnostics — hinge, TT, dt=1, p=0.02" \
  --subtitle_mode \
  --legend_plot drift diffusion tau \
  --fig_w 4.5 --fig_h 3.2 --dpi 300
"""

import argparse
from pathlib import Path
from typing import Optional, Tuple, Iterable, Set

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# utilities
# ---------------------------------------------------------------------

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def read_csv_optional(path: Optional[str], label: str) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{label} not found: {p}")
    return pd.read_csv(p)


def set_mpl_defaults(font: str, font_size: float) -> None:
    plt.rcParams.update({
        "font.family": font,
        "font.size": font_size,
        "axes.titlesize": font_size,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": font_size,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def parse_lim(s: Optional[str]) -> Optional[Tuple[float, float]]:
    if s is None or str(s).strip() == "":
        return None
    vals = [float(v.strip()) for v in s.split(",")]
    if len(vals) != 2:
        raise ValueError(f"Expected two comma-separated numbers, got: {s}")
    return (vals[0], vals[1])


def finite_tau_from_slope(slopes: np.ndarray) -> np.ndarray:
    slopes = np.asarray(slopes, dtype=float)
    tau = 1.0 / (-slopes)
    tau = tau[np.isfinite(tau) & (tau > 0)]
    return tau


def maybe_apply_limits(ax, xlim=None, ylim=None) -> None:
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)


def draw_titles(ax, global_title: str, subtitle: str) -> None:
    if global_title and subtitle:
        ax.set_title(f"{global_title}\n{subtitle}")
    elif global_title:
        ax.set_title(global_title)
    elif subtitle:
        ax.set_title(subtitle)


def save_fig(fig, outpath: Path, dpi: int, also_pdf: bool) -> None:
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    if also_pdf:
        fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def resolve_fig_params(
    args,
    plot_key: str,
) -> Tuple[float, float, int]:
    w = getattr(args, f"fig_w_{plot_key}")
    h = getattr(args, f"fig_h_{plot_key}")
    d = getattr(args, f"dpi_{plot_key}")
    return (
        args.fig_w if w is None else w,
        args.fig_h if h is None else h,
        args.dpi if d is None else d,
    )


def wanted_plots(args) -> Set[str]:
    any_specific = any([
        args.plot_drift,
        args.plot_diffusion,
        args.plot_J,
        args.plot_S,
        args.plot_tau,
        args.plot_all,
    ])
    if args.plot_all or not any_specific:
        return {"drift", "diffusion", "J", "S", "tau"}
    out = set()
    if args.plot_drift:
        out.add("drift")
    if args.plot_diffusion:
        out.add("diffusion")
    if args.plot_J:
        out.add("J")
    if args.plot_S:
        out.add("S")
    if args.plot_tau:
        out.add("tau")
    return out


def legend_targets(args) -> Set[str]:
    if args.legend_all:
        return {"drift", "diffusion", "J", "S", "tau"}
    return set(args.legend_plot or [])


def should_show_legend(plot_key: str, legend_set: Set[str]) -> bool:
    return plot_key in legend_set


def subtitle_for(plot_key: str, subtitle_mode: bool) -> str:
    if not subtitle_mode:
        return ""
    mapping = {
        "drift": "Bootstrap drift estimate",
        "diffusion": "Bootstrap diffusion estimate",
        "J": "Bootstrap flux estimate",
        "S": "Bootstrap source estimate",
        "tau": "Bootstrap distribution of stabilization time",
    }
    return mapping[plot_key]


def add_legend(ax, show: bool, frameon: bool, loc: str) -> None:
    if show:
        ax.legend(loc=loc, frameon=frameon)


# ---------------------------------------------------------------------
# plotting functions
# ---------------------------------------------------------------------

def plot_curve_with_ci(
    df: pd.DataFrame,
    xcol: str,
    ycol: str,
    lo_col: str,
    hi_col: str,
    outpath: Path,
    fig_w: float,
    fig_h: float,
    dpi: int,
    line_width: float,
    ci_alpha: float,
    line_label: str,
    ci_label: str,
    xlabel: str,
    ylabel: str,
    title: str,
    subtitle: str,
    xlim: Optional[Tuple[float, float]],
    ylim: Optional[Tuple[float, float]],
    show_hline0: bool,
    grid: bool,
    show_legend: bool,
    legend_frameon: bool,
    legend_loc: str,
    also_pdf: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    x = df[xcol].to_numpy(dtype=float)
    y = df[ycol].to_numpy(dtype=float)
    lo = df[lo_col].to_numpy(dtype=float)
    hi = df[hi_col].to_numpy(dtype=float)

    if show_hline0:
        ax.axhline(0.0, linewidth=0.8)

    ax.plot(x, y, linewidth=line_width, label=line_label)
    ax.fill_between(x, lo, hi, alpha=ci_alpha, label=ci_label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    draw_titles(ax, title, subtitle)
    maybe_apply_limits(ax, xlim=xlim, ylim=ylim)
    if grid:
        ax.grid(True, alpha=0.25, linewidth=0.5)
    add_legend(ax, show_legend, legend_frameon, legend_loc)
    save_fig(fig, outpath, dpi, also_pdf)


def plot_tau_histogram(
    params_df: pd.DataFrame,
    outpath: Path,
    fig_w: float,
    fig_h: float,
    dpi: int,
    bins: int,
    hist_alpha: float,
    line_width: float,
    show_mean: bool,
    show_median: bool,
    show_ci: bool,
    ci: float,
    xlabel: str,
    ylabel: str,
    title: str,
    subtitle: str,
    xlim: Optional[Tuple[float, float]],
    ylim: Optional[Tuple[float, float]],
    grid: bool,
    show_legend: bool,
    legend_frameon: bool,
    legend_loc: str,
    also_pdf: bool,
) -> None:
    if "slope" not in params_df.columns:
        raise ValueError("drift_boot_params.csv must contain column 'slope' to compute tau.")

    tau = finite_tau_from_slope(params_df["slope"].to_numpy(dtype=float))
    if len(tau) == 0:
        raise ValueError("No valid tau values could be computed from slope.")

    alpha = 0.5 * (1.0 - ci)
    tau_mean = float(np.mean(tau))
    tau_median = float(np.median(tau))
    tau_lo = float(np.quantile(tau, alpha))
    tau_hi = float(np.quantile(tau, 1.0 - alpha))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.hist(tau, bins=bins, alpha=hist_alpha, density=False, label="bootstrap")

    if show_mean:
        ax.axvline(tau_mean, linewidth=line_width, linestyle="-", label=f"mean = {tau_mean:.3f}")
    if show_median:
        ax.axvline(tau_median, linewidth=line_width, linestyle="--", label=f"median = {tau_median:.3f}")
    if show_ci:
        ax.axvline(tau_lo, linewidth=line_width, linestyle=":", label=f"CI low = {tau_lo:.3f}")
        ax.axvline(tau_hi, linewidth=line_width, linestyle=":", label=f"CI high = {tau_hi:.3f}")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    draw_titles(ax, title, subtitle)
    maybe_apply_limits(ax, xlim=xlim, ylim=ylim)
    if grid:
        ax.grid(True, alpha=0.25, linewidth=0.5)
    add_legend(ax, show_legend, legend_frameon, legend_loc)
    save_fig(fig, outpath, dpi, also_pdf)


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Plot full-pipeline bootstrap outputs with flexible titles, subtitles, legends, and per-plot figure sizes.")

    # input files
    ap.add_argument("--drift_summary", default=None, help="Path to drift_boot_summary.csv")
    ap.add_argument("--diffusion_summary", default=None, help="Path to diffusion_boot_summary.csv")
    ap.add_argument("--J_summary", default=None, help="Path to J_boot_summary.csv")
    ap.add_argument("--S_summary", default=None, help="Path to S_boot_summary.csv")
    ap.add_argument("--drift_params", default=None, help="Path to drift_boot_params.csv")

    # output
    ap.add_argument("--out_dir", required=True, help="Output directory for figures")
    ap.add_argument("--fmt", default="png", choices=["png", "pdf", "svg"])
    ap.add_argument("--also_pdf", action="store_true", help="Also save a PDF copy for each output")
    ap.add_argument("--prefix", default="", help="Optional prefix for output filenames")

    # plot selection
    ap.add_argument("--plot_drift", action="store_true")
    ap.add_argument("--plot_diffusion", action="store_true")
    ap.add_argument("--plot_J", action="store_true")
    ap.add_argument("--plot_S", action="store_true")
    ap.add_argument("--plot_tau", action="store_true")
    ap.add_argument("--plot_all", action="store_true")

    # global style
    ap.add_argument("--font", default="Arial")
    ap.add_argument("--font_size", type=float, default=10)
    ap.add_argument("--fig_w", type=float, default=4.0, help="Global figure width in inches")
    ap.add_argument("--fig_h", type=float, default=3.0, help="Global figure height in inches")
    ap.add_argument("--dpi", type=int, default=300, help="Global dpi")
    ap.add_argument("--line_width", type=float, default=1.4)
    ap.add_argument("--ci_alpha", type=float, default=0.25)
    ap.add_argument("--hist_alpha", type=float, default=0.5)

    # per-plot figure size / dpi overrides
    ap.add_argument("--fig_w_drift", type=float, default=None, help="Drift plot width in inches")
    ap.add_argument("--fig_h_drift", type=float, default=None, help="Drift plot height in inches")
    ap.add_argument("--dpi_drift", type=int, default=None, help="Drift plot dpi")

    ap.add_argument("--fig_w_diffusion", type=float, default=None, help="Diffusion plot width in inches")
    ap.add_argument("--fig_h_diffusion", type=float, default=None, help="Diffusion plot height in inches")
    ap.add_argument("--dpi_diffusion", type=int, default=None, help="Diffusion plot dpi")

    ap.add_argument("--fig_w_J", type=float, default=None, help="Flux J plot width in inches")
    ap.add_argument("--fig_h_J", type=float, default=None, help="Flux J plot height in inches")
    ap.add_argument("--dpi_J", type=int, default=None, help="Flux J plot dpi")

    ap.add_argument("--fig_w_S", type=float, default=None, help="Source S plot width in inches")
    ap.add_argument("--fig_h_S", type=float, default=None, help="Source S plot height in inches")
    ap.add_argument("--dpi_S", type=int, default=None, help="Source S plot dpi")

    ap.add_argument("--fig_w_tau", type=float, default=None, help="Tau histogram width in inches")
    ap.add_argument("--fig_h_tau", type=float, default=None, help="Tau histogram height in inches")
    ap.add_argument("--dpi_tau", type=int, default=None, help="Tau histogram dpi")

    # titles / subtitles
    ap.add_argument("--title", default="", help="If provided, shown on all generated plots; if omitted, no title is shown")
    ap.add_argument("--subtitle_mode", action="store_true", help="If provided, show built-in per-plot subtitles")

    # legends
    ap.add_argument("--legend_all", action="store_true", help="Show legends on all generated plots")
    ap.add_argument(
        "--legend_plot",
        nargs="+",
        choices=["drift", "diffusion", "J", "S", "tau"],
        default=None,
        help="Show legends only on selected plots. Example: --legend_plot drift diffusion tau",
    )
    ap.add_argument("--legend_loc", default="best", help="Legend location for plots where legend is enabled")
    ap.add_argument("--legend_frameon", action="store_true", help="Draw legend frame")

    # grid / reference lines
    ap.add_argument("--grid", action="store_true", default=False)
    ap.add_argument("--show_hline0_J", action="store_true", help="Show y=0 line on J plot")
    ap.add_argument("--show_hline0_S", action="store_true", help="Show y=0 line on S plot")
    ap.add_argument("--show_hline0_drift", action="store_true", help="Show y=0 line on drift plot")

    # labels
    ap.add_argument("--xlabel_x", default="x = ln f")
    ap.add_argument("--ylabel_b", default="b(x)")
    ap.add_argument("--ylabel_D", default="D(x)")
    ap.add_argument("--ylabel_J", default="J(x)")
    ap.add_argument("--ylabel_S", default="S(x) = dJ/dx")
    ap.add_argument("--xlabel_tau", default="tau = 1/(-slope) [weeks]")
    ap.add_argument("--ylabel_tau", default="count")

    # axis limits
    ap.add_argument("--xlim_drift", default=None)
    ap.add_argument("--ylim_drift", default=None)
    ap.add_argument("--xlim_diffusion", default=None)
    ap.add_argument("--ylim_diffusion", default=None)
    ap.add_argument("--xlim_J", default=None)
    ap.add_argument("--ylim_J", default=None)
    ap.add_argument("--xlim_S", default=None)
    ap.add_argument("--ylim_S", default=None)
    ap.add_argument("--xlim_tau", default=None)
    ap.add_argument("--ylim_tau", default=None)

    # tau histogram options
    ap.add_argument("--tau_bins", type=int, default=30)
    ap.add_argument("--tau_ci", type=float, default=0.95)
    ap.add_argument("--tau_show_mean", action="store_true", default=False)
    ap.add_argument("--tau_show_median", action="store_true", default=False)
    ap.add_argument("--tau_show_ci", action="store_true", default=False)

    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    ensure_dir(out_dir)
    set_mpl_defaults(args.font, args.font_size)

    requested = wanted_plots(args)
    legend_set = legend_targets(args)

    drift_df = read_csv_optional(args.drift_summary, "drift_summary")
    diff_df = read_csv_optional(args.diffusion_summary, "diffusion_summary")
    J_df = read_csv_optional(args.J_summary, "J_summary")
    S_df = read_csv_optional(args.S_summary, "S_summary")
    params_df = read_csv_optional(args.drift_params, "drift_params")

    produced = []

    # drift
    if "drift" in requested and drift_df is not None:
        required = {"x", "b_mean", "b_lo", "b_hi"}
        if not required.issubset(drift_df.columns):
            raise ValueError(f"drift_summary missing columns: {required - set(drift_df.columns)}")
        fig_w, fig_h, dpi = resolve_fig_params(args, "drift")
        outpath = out_dir / f"{args.prefix}fig_drift_bootstrap.{args.fmt}"
        plot_curve_with_ci(
            df=drift_df,
            xcol="x",
            ycol="b_mean",
            lo_col="b_lo",
            hi_col="b_hi",
            outpath=outpath,
            fig_w=fig_w,
            fig_h=fig_h,
            dpi=dpi,
            line_width=args.line_width,
            ci_alpha=args.ci_alpha,
            line_label="mean",
            ci_label="CI",
            xlabel=args.xlabel_x,
            ylabel=args.ylabel_b,
            title=args.title,
            subtitle=subtitle_for("drift", args.subtitle_mode),
            xlim=parse_lim(args.xlim_drift),
            ylim=parse_lim(args.ylim_drift),
            show_hline0=args.show_hline0_drift,
            grid=args.grid,
            show_legend=should_show_legend("drift", legend_set),
            legend_frameon=args.legend_frameon,
            legend_loc=args.legend_loc,
            also_pdf=args.also_pdf,
        )
        produced.append(outpath)

    # diffusion
    if "diffusion" in requested and diff_df is not None:
        required = {"x", "D_mean", "D_lo", "D_hi"}
        if not required.issubset(diff_df.columns):
            raise ValueError(f"diffusion_summary missing columns: {required - set(diff_df.columns)}")
        fig_w, fig_h, dpi = resolve_fig_params(args, "diffusion")
        outpath = out_dir / f"{args.prefix}fig_diffusion_bootstrap.{args.fmt}"
        plot_curve_with_ci(
            df=diff_df,
            xcol="x",
            ycol="D_mean",
            lo_col="D_lo",
            hi_col="D_hi",
            outpath=outpath,
            fig_w=fig_w,
            fig_h=fig_h,
            dpi=dpi,
            line_width=args.line_width,
            ci_alpha=args.ci_alpha,
            line_label="mean",
            ci_label="CI",
            xlabel=args.xlabel_x,
            ylabel=args.ylabel_D,
            title=args.title,
            subtitle=subtitle_for("diffusion", args.subtitle_mode),
            xlim=parse_lim(args.xlim_diffusion),
            ylim=parse_lim(args.ylim_diffusion),
            show_hline0=False,
            grid=args.grid,
            show_legend=should_show_legend("diffusion", legend_set),
            legend_frameon=args.legend_frameon,
            legend_loc=args.legend_loc,
            also_pdf=args.also_pdf,
        )
        produced.append(outpath)

    # J
    if "J" in requested and J_df is not None:
        required = {"x", "J_mean", "J_lo", "J_hi"}
        if not required.issubset(J_df.columns):
            raise ValueError(f"J_summary missing columns: {required - set(J_df.columns)}")
        fig_w, fig_h, dpi = resolve_fig_params(args, "J")
        outpath = out_dir / f"{args.prefix}fig_flux_J_bootstrap.{args.fmt}"
        plot_curve_with_ci(
            df=J_df,
            xcol="x",
            ycol="J_mean",
            lo_col="J_lo",
            hi_col="J_hi",
            outpath=outpath,
            fig_w=fig_w,
            fig_h=fig_h,
            dpi=dpi,
            line_width=args.line_width,
            ci_alpha=args.ci_alpha,
            line_label="mean",
            ci_label="CI",
            xlabel=args.xlabel_x,
            ylabel=args.ylabel_J,
            title=args.title,
            subtitle=subtitle_for("J", args.subtitle_mode),
            xlim=parse_lim(args.xlim_J),
            ylim=parse_lim(args.ylim_J),
            show_hline0=args.show_hline0_J,
            grid=args.grid,
            show_legend=should_show_legend("J", legend_set),
            legend_frameon=args.legend_frameon,
            legend_loc=args.legend_loc,
            also_pdf=args.also_pdf,
        )
        produced.append(outpath)

    # S
    if "S" in requested and S_df is not None:
        required = {"x", "S_mean", "S_lo", "S_hi"}
        if not required.issubset(S_df.columns):
            raise ValueError(f"S_summary missing columns: {required - set(S_df.columns)}")
        fig_w, fig_h, dpi = resolve_fig_params(args, "S")
        outpath = out_dir / f"{args.prefix}fig_source_S_bootstrap.{args.fmt}"
        plot_curve_with_ci(
            df=S_df,
            xcol="x",
            ycol="S_mean",
            lo_col="S_lo",
            hi_col="S_hi",
            outpath=outpath,
            fig_w=fig_w,
            fig_h=fig_h,
            dpi=dpi,
            line_width=args.line_width,
            ci_alpha=args.ci_alpha,
            line_label="mean",
            ci_label="CI",
            xlabel=args.xlabel_x,
            ylabel=args.ylabel_S,
            title=args.title,
            subtitle=subtitle_for("S", args.subtitle_mode),
            xlim=parse_lim(args.xlim_S),
            ylim=parse_lim(args.ylim_S),
            show_hline0=args.show_hline0_S,
            grid=args.grid,
            show_legend=should_show_legend("S", legend_set),
            legend_frameon=args.legend_frameon,
            legend_loc=args.legend_loc,
            also_pdf=args.also_pdf,
        )
        produced.append(outpath)

    # tau
    if "tau" in requested and params_df is not None:
        fig_w, fig_h, dpi = resolve_fig_params(args, "tau")
        outpath = out_dir / f"{args.prefix}fig_tau_distribution.{args.fmt}"
        plot_tau_histogram(
            params_df=params_df,
            outpath=outpath,
            fig_w=fig_w,
            fig_h=fig_h,
            dpi=dpi,
            bins=args.tau_bins,
            hist_alpha=args.hist_alpha,
            line_width=args.line_width,
            show_mean=args.tau_show_mean,
            show_median=args.tau_show_median,
            show_ci=args.tau_show_ci,
            ci=args.tau_ci,
            xlabel=args.xlabel_tau,
            ylabel=args.ylabel_tau,
            title=args.title,
            subtitle=subtitle_for("tau", args.subtitle_mode),
            xlim=parse_lim(args.xlim_tau),
            ylim=parse_lim(args.ylim_tau),
            grid=args.grid,
            show_legend=should_show_legend("tau", legend_set),
            legend_frameon=args.legend_frameon,
            legend_loc=args.legend_loc,
            also_pdf=args.also_pdf,
        )
        produced.append(outpath)

    if not produced:
        raise ValueError("No input files were provided, or requested plots had no matching input files. Nothing to plot.")

    print("[DONE] Produced figures:")
    for p in produced:
        print(f"  {p}")


if __name__ == "__main__":
    main()
