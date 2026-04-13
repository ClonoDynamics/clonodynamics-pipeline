#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_finite_time_diagnostics.py

Unified plotting for outputs produced by the unified analysis script
(local model-agnostic + finite-time summaries).

This version fixes a common issue of the original script:
changing figsize had little visible effect because figures were later
re-packed by tight_layout / constrained_layout and finally cropped by
bbox_inches="tight" at save time.

Main fixes
----------
- no tight_layout()
- no constrained_layout=True
- no bbox_inches="tight" in savefig
- explicit adaptive margins
- outside legends anchored to the figure
- axis labels spaced with labelpad
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =========================================================
# CLI
# =========================================================
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Unified plotting for local model-agnostic and finite-time outputs."
    )

    ap.add_argument("--in_dir", required=True, help="Directory produced by unified analysis script")
    ap.add_argument("--out_dir", default="", help="Output directory (default: in_dir)")
    ap.add_argument("--plot_mode", choices=["local", "finite_time", "all"], default="all")
    ap.add_argument("--csv_sep", default=",", help="CSV separator. Default ','")

    ap.add_argument("--fmt", choices=["png", "pdf"], default="png")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--transparent", action="store_true")

    ap.add_argument("--font", default="Arial")
    ap.add_argument("--font_size", type=float, default=10.0)

    ap.add_argument("--line_width", type=float, default=1.2)
    ap.add_argument("--marker_size", type=float, default=4.0)
    ap.add_argument("--ci_alpha", type=float, default=0.12)

    ap.add_argument("--dt_min_plot", type=int, default=0, help="For local panels")
    ap.add_argument("--dt_max_plot", type=int, default=10, help="For local panels")
    ap.add_argument("--dt_values", default="", help="For finite-time panels, comma-separated dt values")

    ap.add_argument("--no_ci", action="store_true", help="Do not draw CI bands even if present")
    ap.add_argument("--title_prefix", default="", help="Optional prefix added to titles")
    ap.add_argument("--title", action="store_true", help="Show titles on finite-time figures")
    ap.add_argument("--sqrt_diffusion", action="store_true", help="Finite-time: plot sqrt(var_mad_dx)")
    ap.add_argument("--zero_line", action="store_true", help="Add y=0 line where relevant")
    ap.add_argument("--xstar_lines", action="store_true", help="Finite-time: add x* vertical lines in drift overlay")

    ap.add_argument("--legend", action="store_true")
    ap.add_argument("--legend_out", action="store_true")
    ap.add_argument("--legend_loc", choices=["right", "bottom"], default="right")
    ap.add_argument("--legend_ncol", type=int, default=1)
    ap.add_argument("--legend_x", type=float, default=0.70, help="Figure-coordinate x anchor for outside legend")
    ap.add_argument("--legend_y", type=float, default=0.50, help="Figure-coordinate y anchor for outside legend")

    # Local panels
    ap.add_argument(
        "--panels",
        default="msd,msd_by_x,msd_by_freq,median,ppos",
        help="Local panels. Valid: msd,msd_by_x,msd_by_freq,median,ppos,counts"
    )
    ap.add_argument("--figsize", default="5,3.5", help="Default single size 'W,H'")
    ap.add_argument("--prefix", default="", help="Filename prefix for local single plots")
    ap.add_argument("--suffix", default="", help="Filename suffix for local single plots")
    ap.add_argument("--figsize_msd", default="")
    ap.add_argument("--figsize_msd_by_x", default="")
    ap.add_argument("--figsize_msd_by_freq", default="")
    ap.add_argument("--figsize_median", default="")
    ap.add_argument("--figsize_ppos", default="")
    ap.add_argument("--figsize_counts", default="")

    ap.add_argument("--make_multi", action="store_true", help="Also create local multipanel figure")
    ap.add_argument("--multi_layout", default="2x3")
    ap.add_argument("--multi_figsize", default="10,6.5")
    ap.add_argument("--multi_outname", default="model_agnostic_multipanel")
    ap.add_argument("--sharex", action="store_true")
    ap.add_argument("--sharey", action="store_true")

    ap.add_argument(
        "--counts_mode",
        choices=["total", "x", "f", "all"],
        default="all",
        help="What to show in counts panel"
    )
    ap.add_argument("--counts_logy", action="store_true")
    ap.add_argument("--counts_max_bins", type=int, default=12)

    # Finite-time panel sizes
    ap.add_argument("--figsize_default_w", type=float, default=5.0)
    ap.add_argument("--figsize_default_h", type=float, default=3.5)
    for panel in ["A", "B", "C", "D", "E", "F", "G"]:
        ap.add_argument(f"--fig{panel}_w", type=float, default=None)
        ap.add_argument(f"--fig{panel}_h", type=float, default=None)

    # Axis limits
    ap.add_argument("--xmin", type=float, default=None)
    ap.add_argument("--xmax", type=float, default=None)
    for panel in ["A", "B", "C", "D", "E", "F", "G"]:
        ap.add_argument(f"--ymin_{panel}", type=float, default=None)
        ap.add_argument(f"--ymax_{panel}", type=float, default=None)

    # Spacing controls
    ap.add_argument("--left_margin", type=float, default=0.16)
    ap.add_argument("--right_margin_in", type=float, default=0.97)
    ap.add_argument("--right_margin_out", type=float, default=0.68)
    ap.add_argument("--bottom_margin", type=float, default=0.17)
    ap.add_argument("--top_margin", type=float, default=0.90)
    ap.add_argument("--wspace", type=float, default=0.35)
    ap.add_argument("--hspace", type=float, default=0.35)

    ap.add_argument("--xlabel_pad", type=float, default=4.0)
    ap.add_argument("--ylabel_pad", type=float, default=9.0)
    ap.add_argument("--title_pad", type=float, default=8.0)

    ap.add_argument("--debug", action="store_true")
    return ap


# =========================================================
# Common utilities
# =========================================================
DT_ALIASES = ["Dt", "DT", "delta_t", "deltaT", "lag", "Lag", "Δt", "dT"]


def set_style(font: str, size: float, title_pad: float) -> None:
    plt.rcParams["font.family"] = font
    plt.rcParams["font.size"] = size
    plt.rcParams["axes.titlesize"] = size
    plt.rcParams["axes.labelsize"] = size
    plt.rcParams["legend.fontsize"] = size
    plt.rcParams["xtick.labelsize"] = size
    plt.rcParams["ytick.labelsize"] = size
    plt.rcParams["axes.titlepad"] = title_pad
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["figure.facecolor"] = "white"


def parse_figsize(s: str, default: Tuple[float, float]) -> Tuple[float, float]:
    s = (s or "").strip()
    if not s:
        return default
    w, h = [p.strip() for p in s.split(",")]
    return float(w), float(h)


def parse_layout(s: str) -> Tuple[int, int]:
    s = s.lower().strip()
    r, c = s.split("x", 1)
    return int(r), int(c)


def parse_panels(s: str) -> List[str]:
    panels = [p.strip() for p in (s or "").split(",") if p.strip()]
    valid = {"msd", "msd_by_x", "msd_by_freq", "median", "ppos", "counts"}
    bad = [p for p in panels if p not in valid]
    if bad:
        raise ValueError(f"Unknown panel(s): {bad}. Valid: {sorted(valid)}")
    return panels


def parse_dt(s: str) -> Optional[List[int]]:
    s = s.strip()
    if not s:
        return None
    return [int(x) for x in s.split(",")]


def read_csv_required(path: Path, sep: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path, sep=sep)


def read_csv_optional(path: Path, sep: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep=sep)


def ensure_dt(df: pd.DataFrame, name: str, path: Path) -> pd.DataFrame:
    if "dt" not in df.columns:
        for a in DT_ALIASES:
            if a in df.columns:
                df = df.copy()
                df["dt"] = df[a]
                break

    if "dt" not in df.columns:
        raise ValueError(
            f"[ERROR] {name} missing required column 'dt'.\n"
            f"  file: {path}\n"
            f"  columns: {list(df.columns)}\n"
            f"  head(5):\n{df.head(5).to_string(index=False)}\n"
        )

    out = df.copy()
    out["dt"] = pd.to_numeric(out["dt"], errors="coerce")
    out = out.dropna(subset=["dt"]).copy()
    out["dt"] = out["dt"].astype(int)
    return out


def require_columns(df: pd.DataFrame, cols: List[str], name: str, path: Path) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"[ERROR] {name} missing required columns {missing}\n"
            f"  file: {path}\n"
            f"  columns: {list(df.columns)}\n"
            f"  head(5):\n{df.head(5).to_string(index=False)}\n"
        )


def filter_dt_range(df: pd.DataFrame, dt_min: int, dt_max: int) -> pd.DataFrame:
    if df.empty or "dt" not in df.columns:
        return df
    d = df.copy()
    return d[(d["dt"] >= int(dt_min)) & (d["dt"] <= int(dt_max))].copy()


def save(fig: plt.Figure, path_base: Path, fmt: str, dpi: int, transparent: bool):
    if fmt == "png":
        fig.savefig(path_base.with_suffix(".png"), dpi=dpi, transparent=transparent)
    else:
        fig.savefig(path_base.with_suffix(".pdf"), transparent=transparent)


def style_axes(ax: plt.Axes) -> None:
    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(True)
    ax.tick_params(direction="out")


def set_axis_labels(ax: plt.Axes, xlabel: str, ylabel: str, args) -> None:
    ax.set_xlabel(xlabel, labelpad=args.xlabel_pad)
    ax.set_ylabel(ylabel, labelpad=args.ylabel_pad)


def place_legend_single(fig: plt.Figure, ax: plt.Axes, args):
    handles, labels = ax.get_legend_handles_labels()
    if not args.legend or not handles:
        return

    if args.legend_out:
        if args.legend_loc == "right":
            fig.legend(
                handles, labels,
                frameon=False,
                loc="center left",
                bbox_to_anchor=(args.legend_x, args.legend_y),
                bbox_transform=fig.transFigure,
                borderaxespad=0.0,
                ncol=max(1, int(args.legend_ncol)),
                handlelength=2.2,
            )
        else:
            fig.legend(
                handles, labels,
                frameon=False,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.02),
                bbox_transform=fig.transFigure,
                borderaxespad=0.0,
                ncol=max(1, int(args.legend_ncol)),
                handlelength=2.2,
            )
    else:
        ax.legend(frameon=False)


def place_legend_multi(fig: plt.Figure, axes: List[plt.Axes], args):
    if not args.legend:
        return

    handles = []
    labels = []
    seen = set()

    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        for hh, ll in zip(h, l):
            if ll not in seen:
                handles.append(hh)
                labels.append(ll)
                seen.add(ll)

    if not handles:
        return

    if args.legend_out:
        if args.legend_loc == "right":
            fig.legend(
                handles, labels,
                frameon=False,
                loc="center left",
                bbox_to_anchor=(args.legend_x, args.legend_y),
                bbox_transform=fig.transFigure,
                borderaxespad=0.0,
                ncol=max(1, int(args.legend_ncol)),
                handlelength=2.2,
            )
        else:
            fig.legend(
                handles, labels,
                frameon=False,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.02),
                bbox_transform=fig.transFigure,
                borderaxespad=0.0,
                ncol=max(1, int(args.legend_ncol)),
                handlelength=2.2,
            )


def _set_dt_ticks(ax, dt_values: np.ndarray):
    dts = sorted(np.unique(dt_values.astype(int)).tolist())
    ax.set_xticks(dts)
    ax.set_xticklabels([str(v) for v in dts])


def debug_df(tag: str, path: Path, df: pd.DataFrame) -> None:
    print(f"[DEBUG] {tag}")
    print(f"  path: {path}")
    print(f"  columns: {list(df.columns)}")
    print(f"  head:\n{df.head(3).to_string(index=False)}")


def get_figsize_ft(args, panel: str) -> Tuple[float, float]:
    w = getattr(args, f"fig{panel}_w")
    h = getattr(args, f"fig{panel}_h")
    if w is None:
        w = args.figsize_default_w
    if h is None:
        h = args.figsize_default_h
    return (w, h)


def apply_common_xlimits(ax, args):
    if args.xmin is not None or args.xmax is not None:
        ax.set_xlim(left=args.xmin, right=args.xmax)


def apply_panel_ylimits(ax, args, panel: str):
    ymin = getattr(args, f"ymin_{panel}")
    ymax = getattr(args, f"ymax_{panel}")
    if ymin is not None or ymax is not None:
        ax.set_ylim(bottom=ymin, top=ymax)


def adaptive_margins(width_in: float, height_in: float, args, has_title: bool, legend_out: bool, legend_loc: str):
    if legend_out and legend_loc == "right":
        if width_in < 5.0:
            right = min(args.right_margin_out, 0.64)
        elif width_in < 8.0:
            right = min(max(args.right_margin_out, 0.68), 0.72)
        else:
            right = min(max(args.right_margin_out, 0.72), 0.78)
    else:
        right = args.right_margin_in

    if width_in < 4.0:
        left = max(args.left_margin, 0.24)
    elif width_in < 5.0:
        left = max(args.left_margin, 0.20)
    else:
        left = max(args.left_margin, 0.16)

    if height_in < 3.0:
        bottom = max(args.bottom_margin, 0.24)
    elif height_in < 3.5:
        bottom = max(args.bottom_margin, 0.20)
    else:
        bottom = max(args.bottom_margin, 0.17)

    if legend_out and legend_loc == "bottom":
        bottom = max(bottom, 0.26)

    if height_in < 3.0:
        top = min(0.82, args.top_margin)
    elif height_in < 3.5:
        top = min(0.85, args.top_margin)
    else:
        top = min(0.88, args.top_margin)

    if has_title:
        top -= 0.04

    top = max(top, bottom + 0.18)
    return left, right, bottom, top


def apply_layout(fig: plt.Figure, width_in: float, height_in: float, args, has_title: bool):
    left, right, bottom, top = adaptive_margins(
        width_in=width_in,
        height_in=height_in,
        args=args,
        has_title=has_title,
        legend_out=args.legend_out,
        legend_loc=args.legend_loc,
    )
    fig.subplots_adjust(
        left=left,
        right=right,
        bottom=bottom,
        top=top,
        wspace=args.wspace,
        hspace=args.hspace,
    )


# =========================================================
# Local panels
# =========================================================
def panel_msd(ax, msd_df: pd.DataFrame, title_prefix: str, no_ci: bool, args):
    if msd_df.empty:
        ax.text(0.5, 0.5, "msd_by_dt.csv empty", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    d = msd_df.sort_values("dt")
    x = d["dt"].to_numpy(int)
    y = d["msd"].to_numpy(float)
    ax.plot(x, y, lw=args.line_width, marker="o", ms=args.marker_size, label="global")

    if (not no_ci) and {"msd_lo", "msd_hi"}.issubset(d.columns):
        lo = d["msd_lo"].to_numpy(float)
        hi = d["msd_hi"].to_numpy(float)
        if np.isfinite(lo).any() and np.isfinite(hi).any():
            ax.fill_between(x, lo, hi, alpha=args.ci_alpha)

    _set_dt_ticks(ax, x)
    set_axis_labels(ax, "Δt (weeks)", "MSD = E[(Δx)^2]", args)
    ttl = "MSD vs Δt"
    ax.set_title(f"{title_prefix}{ttl}" if title_prefix else ttl, pad=args.title_pad)
    style_axes(ax)


def panel_msd_by_x(ax, msd_x: pd.DataFrame, title_prefix: str, no_ci: bool, args):
    if msd_x.empty:
        ax.text(0.5, 0.5, "msd_by_dt_xbins.csv empty", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    for dt, g in msd_x.groupby("dt", sort=True):
        g = g.sort_values("x_center")
        x = g["x_center"].to_numpy(float)
        y = g["msd"].to_numpy(float)
        ax.plot(x, y, lw=args.line_width, marker="o", ms=args.marker_size, label=f"Δt={int(dt)}w")

        if (not no_ci) and {"msd_lo", "msd_hi"}.issubset(g.columns):
            lo = g["msd_lo"].to_numpy(float)
            hi = g["msd_hi"].to_numpy(float)
            if np.isfinite(lo).any() and np.isfinite(hi).any():
                ax.fill_between(x, lo, hi, alpha=args.ci_alpha)

    set_axis_labels(ax, "Initial ln frequency", "MSD = E[(Δx)^2]", args)
    ttl = "MSD vs initial frequency (by Δt)"
    ax.set_title(f"{title_prefix}{ttl}" if title_prefix else ttl, pad=args.title_pad)
    style_axes(ax)


def panel_msd_by_freq(ax, msd_f: pd.DataFrame, title_prefix: str, no_ci: bool, args):
    if msd_f.empty:
        ax.text(0.5, 0.5, "msd_by_dt_freqbins.csv empty", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    for f_id, g in msd_f.groupby("f_id", sort=True):
        g = g.sort_values("dt")
        x = g["dt"].to_numpy(int)
        y = g["msd"].to_numpy(float)
        ax.plot(x, y, lw=args.line_width, marker="o", ms=args.marker_size, label=f"freq bin {int(f_id)}")

        if (not no_ci) and {"msd_lo", "msd_hi"}.issubset(g.columns):
            lo = g["msd_lo"].to_numpy(float)
            hi = g["msd_hi"].to_numpy(float)
            if np.isfinite(lo).any() and np.isfinite(hi).any():
                ax.fill_between(x, lo, hi, alpha=args.ci_alpha)

    _set_dt_ticks(ax, msd_f["dt"].to_numpy(int))
    set_axis_labels(ax, "Δt (weeks)", "MSD = E[(Δx)^2]", args)
    ttl = "MSD vs Δt by frequency bin"
    ax.set_title(f"{title_prefix}{ttl}" if title_prefix else ttl, pad=args.title_pad)
    style_axes(ax)


def panel_median(ax, drift: pd.DataFrame, title_prefix: str, no_ci: bool, args):
    if drift.empty:
        ax.text(0.5, 0.5, "drift_by_dt_xbins.csv empty", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    for dt, g in drift.groupby("dt", sort=True):
        g = g.sort_values("x_center")
        x = g["x_center"].to_numpy(float)
        y = g["median_dx"].to_numpy(float)
        ax.plot(x, y, lw=args.line_width, label=f"Δt={int(dt)}w")

        if (not no_ci) and {"median_lo", "median_hi"}.issubset(g.columns):
            lo = g["median_lo"].to_numpy(float)
            hi = g["median_hi"].to_numpy(float)
            if np.isfinite(lo).any() and np.isfinite(hi).any():
                ax.fill_between(x, lo, hi, alpha=args.ci_alpha)

    ax.axhline(0.0, linewidth=1.0)
    set_axis_labels(ax, "Initial ln frequency", "median(Δx)", args)
    ttl = "Median(Δx) vs initial frequency (by Δt)"
    ax.set_title(f"{title_prefix}{ttl}" if title_prefix else ttl, pad=args.title_pad)
    style_axes(ax)


def panel_ppos(ax, drift: pd.DataFrame, title_prefix: str, no_ci: bool, args):
    if drift.empty:
        ax.text(0.5, 0.5, "drift_by_dt_xbins.csv empty", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    for dt, g in drift.groupby("dt", sort=True):
        g = g.sort_values("x_center")
        x = g["x_center"].to_numpy(float)
        y = g["ppos"].to_numpy(float)
        ax.plot(x, y, lw=args.line_width, label=f"Δt={int(dt)}w")

        if (not no_ci) and {"ppos_lo", "ppos_hi"}.issubset(g.columns):
            lo = g["ppos_lo"].to_numpy(float)
            hi = g["ppos_hi"].to_numpy(float)
            if np.isfinite(lo).any() and np.isfinite(hi).any():
                ax.fill_between(x, lo, hi, alpha=args.ci_alpha)

    ax.axhline(0.5, linewidth=1.0)
    set_axis_labels(ax, "Initial ln frequency", "P(Δx > 0)", args)
    ttl = "P(Δx>0) vs initial frequency (by Δt)"
    ax.set_title(f"{title_prefix}{ttl}" if title_prefix else ttl, pad=args.title_pad)
    style_axes(ax)


def panel_counts(ax, counts: pd.DataFrame, title_prefix: str, counts_mode: str, counts_logy: bool, counts_max_bins: int, args):
    if counts.empty:
        ax.text(0.5, 0.5, "dt_freqbin_counts.csv missing/empty", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    require_columns(counts, ["bin_type", "dt", "n"], "dt_freqbin_counts.csv", Path("dt_freqbin_counts.csv"))

    d = counts.copy()
    d["bin_type"] = d["bin_type"].astype(str)

    if counts_mode == "total":
        types = ["TOTAL"]
    elif counts_mode == "x":
        types = ["ln_f0"]
    elif counts_mode == "f":
        types = ["ln_fgeo"]
    else:
        types = ["TOTAL", "ln_f0", "ln_fgeo"]

    plotted = False

    for bt in types:
        dd = d[d["bin_type"] == bt].copy()
        if dd.empty:
            continue
        dd = dd.sort_values("dt")

        if bt == "TOTAL":
            ax.plot(dd["dt"].to_numpy(int), dd["n"].to_numpy(float),
                    marker="o", lw=args.line_width, ms=args.marker_size, label="TOTAL")
            plotted = True
            continue

        if "bin_id" not in dd.columns:
            continue

        dd["bin_id"] = pd.to_numeric(dd["bin_id"], errors="coerce")
        dd = dd.dropna(subset=["bin_id"]).copy()
        dd["bin_id"] = dd["bin_id"].astype(int)

        score = dd.groupby("bin_id")["n"].sum().sort_values(ascending=False)
        top_bins = score.head(int(counts_max_bins)).index.tolist()

        for bid in sorted(top_bins):
            g = dd[dd["bin_id"] == bid].sort_values("dt")
            ax.plot(g["dt"].to_numpy(int), g["n"].to_numpy(float),
                    marker="o", lw=args.line_width, ms=args.marker_size, label=f"{bt} bin {bid}")
            plotted = True

    if not plotted:
        ax.text(0.5, 0.5, "counts_mode matched no rows", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    _set_dt_ticks(ax, d["dt"].to_numpy(int))
    set_axis_labels(ax, "Δt (weeks)", "Number of transitions (n)", args)
    if counts_logy:
        ax.set_yscale("log")
    ttl = "Transition counts vs Δt"
    ax.set_title(f"{title_prefix}{ttl}" if title_prefix else ttl, pad=args.title_pad)
    style_axes(ax)


# =========================================================
# Finite-time panels
# =========================================================
def plot_drift_overlay(moments, summary, out_dir, args):
    w, h = get_figsize_ft(args, "A")
    fig, ax = plt.subplots(figsize=(w, h))

    for dt, g in moments.groupby("dt", sort=True):
        g = g.sort_values("x_center")
        ax.plot(g["x_center"], g["median_dx"], lw=args.line_width, label=f"Δt={dt}")

        if (not args.no_ci) and "median_lo" in g.columns and "median_hi" in g.columns:
            ax.fill_between(g["x_center"], g["median_lo"], g["median_hi"], alpha=args.ci_alpha)

    if args.zero_line:
        ax.axhline(0, lw=1.0)

    if args.xstar_lines:
        for _, r in summary.iterrows():
            if "x_star" in r and np.isfinite(r["x_star"]):
                ax.axvline(r["x_star"], ls="--", lw=1.0, alpha=0.6)

    set_axis_labels(ax, "Initial ln frequency", "median(Δx)", args)
    if args.title:
        ax.set_title("Finite-time drift curves", pad=args.title_pad)

    apply_common_xlimits(ax, args)
    apply_panel_ylimits(ax, args, "A")
    style_axes(ax)
    place_legend_single(fig, ax, args)
    apply_layout(fig, w, h, args, has_title=args.title)
    save(fig, out_dir / "figA_drift_overlay", args.fmt, args.dpi, args.transparent)
    plt.close(fig)


def plot_diffusion_overlay(moments, out_dir, args):
    w, h = get_figsize_ft(args, "B")
    fig, ax = plt.subplots(figsize=(w, h))

    for dt, g in moments.groupby("dt", sort=True):
        g = g.sort_values("x_center")

        if args.sqrt_diffusion:
            y = np.sqrt(np.clip(g["var_mad_dx"].to_numpy(), a_min=0, a_max=None))
            ylabel = r"$\sqrt{\mathrm{var}_{MAD}(\Delta x)}$"

            if (not args.no_ci) and "var_mad_lo" in g.columns and "var_mad_hi" in g.columns:
                lo = np.sqrt(np.clip(g["var_mad_lo"].to_numpy(), a_min=0, a_max=None))
                hi = np.sqrt(np.clip(g["var_mad_hi"].to_numpy(), a_min=0, a_max=None))
            else:
                lo = hi = None
        else:
            y = g["var_mad_dx"].to_numpy()
            ylabel = r"$\mathrm{var}_{MAD}(\Delta x)$"

            if (not args.no_ci) and "var_mad_lo" in g.columns and "var_mad_hi" in g.columns:
                lo = g["var_mad_lo"].to_numpy()
                hi = g["var_mad_hi"].to_numpy()
            else:
                lo = hi = None

        ax.plot(g["x_center"], y, lw=args.line_width, label=f"Δt={dt}")

        if lo is not None and hi is not None:
            ax.fill_between(g["x_center"], lo, hi, alpha=args.ci_alpha)

    set_axis_labels(ax, "Initial ln frequency", ylabel, args)
    if args.title:
        ax.set_title("Finite-time diffusion (robust MAD-based)", pad=args.title_pad)

    apply_common_xlimits(ax, args)
    apply_panel_ylimits(ax, args, "B")
    style_axes(ax)
    place_legend_single(fig, ax, args)
    apply_layout(fig, w, h, args, has_title=args.title)
    save(fig, out_dir / "figB_diffusion_overlay", args.fmt, args.dpi, args.transparent)
    plt.close(fig)


def plot_vs_dt(summary, col, ylabel, title, fname, panel, out_dir, args):
    w, h = get_figsize_ft(args, panel)
    fig, ax = plt.subplots(figsize=(w, h))

    summary = summary.sort_values("dt")
    ax.plot(summary["dt"], summary[col],
            marker="o", markersize=args.marker_size, lw=args.line_width)

    lo_col = col + "_lo"
    hi_col = col + "_hi"
    if (not args.no_ci) and lo_col in summary.columns and hi_col in summary.columns:
        ax.fill_between(summary["dt"], summary[lo_col], summary[hi_col], alpha=args.ci_alpha)

    set_axis_labels(ax, "Δt (weeks)", ylabel, args)
    ax.set_xticks(sorted(summary["dt"].unique()))

    vals = summary[col].to_numpy(float)
    if args.zero_line and np.isfinite(vals).any():
        if np.nanmin(vals) < 0 < np.nanmax(vals):
            ax.axhline(0, lw=1.0)

    if args.title:
        ax.set_title(title, pad=args.title_pad)

    apply_panel_ylimits(ax, args, panel)
    style_axes(ax)
    apply_layout(fig, w, h, args, has_title=args.title)
    save(fig, out_dir / fname, args.fmt, args.dpi, args.transparent)
    plt.close(fig)


# =========================================================
# Main
# =========================================================
def main():
    args = build_argparser().parse_args()

    in_dir = Path(args.in_dir).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    set_style(args.font, args.font_size, args.title_pad)

    if args.plot_mode in {"local", "all"}:
        panels = parse_panels(args.panels)

        drift_path = in_dir / "drift_by_dt_xbins.csv"
        msd_path = in_dir / "msd_by_dt.csv"
        msd_x_path = in_dir / "msd_by_dt_xbins.csv"
        msd_f_path = in_dir / "msd_by_dt_freqbins.csv"
        counts_path = in_dir / "dt_freqbin_counts.csv"

        drift = read_csv_required(drift_path, sep=args.csv_sep)
        msd = read_csv_required(msd_path, sep=args.csv_sep)
        msd_x = read_csv_required(msd_x_path, sep=args.csv_sep)
        msd_f = read_csv_required(msd_f_path, sep=args.csv_sep)
        counts = read_csv_optional(counts_path, sep=args.csv_sep)

        drift = ensure_dt(drift, "drift_by_dt_xbins.csv", drift_path)
        msd = ensure_dt(msd, "msd_by_dt.csv", msd_path)
        msd_x = ensure_dt(msd_x, "msd_by_dt_xbins.csv", msd_x_path)
        msd_f = ensure_dt(msd_f, "msd_by_dt_freqbins.csv", msd_f_path)
        if not counts.empty:
            counts = ensure_dt(counts, "dt_freqbin_counts.csv", counts_path)

        require_columns(msd, ["dt", "msd"], "msd_by_dt.csv", msd_path)
        require_columns(msd_x, ["dt", "x_id", "x_center", "msd"], "msd_by_dt_xbins.csv", msd_x_path)
        require_columns(msd_f, ["f_id", "dt", "msd"], "msd_by_dt_freqbins.csv", msd_f_path)
        require_columns(drift, ["dt", "x_center", "median_dx", "ppos"], "drift_by_dt_xbins.csv", drift_path)

        if args.debug:
            print("[DEBUG] LOCAL in_dir =", in_dir)
            debug_df("msd_by_dt", msd_path, msd)
            debug_df("msd_by_dt_xbins", msd_x_path, msd_x)
            debug_df("msd_by_dt_freqbins", msd_f_path, msd_f)
            debug_df("drift_by_dt_xbins", drift_path, drift)
            if not counts.empty:
                debug_df("dt_freqbin_counts", counts_path, counts)

        drift = filter_dt_range(drift, args.dt_min_plot, args.dt_max_plot)
        msd = filter_dt_range(msd, args.dt_min_plot, args.dt_max_plot)
        msd_x = filter_dt_range(msd_x, args.dt_min_plot, args.dt_max_plot)
        msd_f = filter_dt_range(msd_f, args.dt_min_plot, args.dt_max_plot)
        counts = filter_dt_range(counts, args.dt_min_plot, args.dt_max_plot) if not counts.empty else counts

        default_fs = parse_figsize(args.figsize, default=(5.0, 3.5))
        per_fs: Dict[str, Tuple[float, float]] = {
            "msd": parse_figsize(args.figsize_msd, default_fs),
            "msd_by_x": parse_figsize(args.figsize_msd_by_x, default_fs),
            "msd_by_freq": parse_figsize(args.figsize_msd_by_freq, default_fs),
            "median": parse_figsize(args.figsize_median, default_fs),
            "ppos": parse_figsize(args.figsize_ppos, default_fs),
            "counts": parse_figsize(args.figsize_counts, default_fs),
        }

        for p in panels:
            w, h = per_fs[p]
            fig, ax = plt.subplots(figsize=(w, h))

            if p == "msd":
                panel_msd(ax, msd, args.title_prefix, args.no_ci, args)
            elif p == "msd_by_x":
                panel_msd_by_x(ax, msd_x, args.title_prefix, args.no_ci, args)
            elif p == "msd_by_freq":
                panel_msd_by_freq(ax, msd_f, args.title_prefix, args.no_ci, args)
            elif p == "median":
                panel_median(ax, drift, args.title_prefix, args.no_ci, args)
            elif p == "ppos":
                panel_ppos(ax, drift, args.title_prefix, args.no_ci, args)
            elif p == "counts":
                panel_counts(ax, counts, args.title_prefix, args.counts_mode, args.counts_logy, args.counts_max_bins, args)

            place_legend_single(fig, ax, args)
            apply_layout(fig, w, h, args, has_title=True)

            fname = f"{args.prefix}{p}{args.suffix}"
            save(fig, out_dir / fname, args.fmt, args.dpi, args.transparent)
            plt.close(fig)

        print(f"[SAVED LOCAL SINGLE PLOTS] {out_dir.resolve()} ({args.fmt})")

        if args.make_multi:
            nrows, ncols = parse_layout(args.multi_layout)
            multi_fs = parse_figsize(args.multi_figsize, default=(10.0, 6.5))
            fig_w, fig_h = multi_fs

            fig, axes = plt.subplots(
                nrows=nrows,
                ncols=ncols,
                figsize=multi_fs,
                sharex=args.sharex,
                sharey=args.sharey,
                squeeze=False,
            )
            ax_list = axes.ravel().tolist()
            used_axes = []

            for i, p in enumerate(panels):
                if i >= len(ax_list):
                    break
                ax = ax_list[i]
                used_axes.append(ax)

                if p == "msd":
                    panel_msd(ax, msd, args.title_prefix, args.no_ci, args)
                elif p == "msd_by_x":
                    panel_msd_by_x(ax, msd_x, args.title_prefix, args.no_ci, args)
                elif p == "msd_by_freq":
                    panel_msd_by_freq(ax, msd_f, args.title_prefix, args.no_ci, args)
                elif p == "median":
                    panel_median(ax, drift, args.title_prefix, args.no_ci, args)
                elif p == "ppos":
                    panel_ppos(ax, drift, args.title_prefix, args.no_ci, args)
                elif p == "counts":
                    panel_counts(ax, counts, args.title_prefix, args.counts_mode, args.counts_logy, args.counts_max_bins, args)

            for j in range(len(panels), len(ax_list)):
                ax_list[j].axis("off")

            place_legend_multi(fig, used_axes, args)
            apply_layout(fig, fig_w, fig_h, args, has_title=True)

            multi_name = f"{args.prefix}{args.multi_outname}{args.suffix}"
            save(fig, out_dir / multi_name, args.fmt, args.dpi, args.transparent)
            plt.close(fig)

            print(f"[SAVED LOCAL MULTIPANEL] {(out_dir / multi_name).with_suffix('.' + args.fmt).resolve()}")

    if args.plot_mode in {"finite_time", "all"}:
        summary_path = in_dir / "finite_time_summary.csv"
        moments_path = in_dir / "finite_time_by_dt_xbins_boot.csv"

        summary = read_csv_required(summary_path, sep=args.csv_sep)
        moments = read_csv_required(moments_path, sep=args.csv_sep)

        require_columns(summary, ["dt"], "finite_time_summary.csv", summary_path)
        require_columns(moments, ["dt", "x_center", "median_dx", "var_mad_dx"], "finite_time_by_dt_xbins_boot.csv", moments_path)

        summary["dt"] = pd.to_numeric(summary["dt"], errors="coerce").astype("Int64")
        moments["dt"] = pd.to_numeric(moments["dt"], errors="coerce").astype("Int64")
        summary = summary.dropna(subset=["dt"]).copy()
        moments = moments.dropna(subset=["dt"]).copy()
        summary["dt"] = summary["dt"].astype(int)
        moments["dt"] = moments["dt"].astype(int)

        dt_keep = parse_dt(args.dt_values)
        if dt_keep:
            summary = summary[summary["dt"].isin(dt_keep)].copy()
            moments = moments[moments["dt"].isin(dt_keep)].copy()

        if args.debug:
            print("[DEBUG] FINITE_TIME in_dir =", in_dir)
            debug_df("finite_time_summary", summary_path, summary)
            debug_df("finite_time_by_dt_xbins_boot", moments_path, moments)

        plot_drift_overlay(moments, summary, out_dir, args)
        plot_diffusion_overlay(moments, out_dir, args)

        plot_vs_dt(summary, "x_star", r"$x^*(\Delta t)$", r"Zero-crossing $x^*(\Delta t)$", "figC_xstar_vs_dt", "C", out_dir, args)
        plot_vs_dt(summary, "slope", r"slope$(\Delta t)$", r"Local drift slope vs $\Delta t$", "figD_slope_vs_dt", "D", out_dir, args)
        plot_vs_dt(summary, "tau_drift", r"$\tau_{\mathrm{drift}}(\Delta t)$", r"Apparent drift timescale $\tau_{\mathrm{drift}}(\Delta t)$", "figE_tau_vs_dt", "E", out_dir, args)
        plot_vs_dt(summary, "var_mad_pooled", r"$\mathrm{var}_{MAD}^{pooled}(\Delta x)$", r"Pooled robust displacement variance", "figF_var_pooled_vs_dt", "F", out_dir, args)
        plot_vs_dt(summary, "D_mad_pooled", r"$D_{MAD}^{pooled}(\Delta t)$", r"Pooled finite-time diffusion proxy vs $\Delta t$", "figG_D_pooled_vs_dt", "G", out_dir, args)

        print(f"[SAVED FINITE-TIME FIGURES] {out_dir.resolve()} ({args.fmt})")


if __name__ == "__main__":
    main()