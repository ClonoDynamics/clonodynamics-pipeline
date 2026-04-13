#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_finite_time_fit.py

OU finite-time validation (PLOT ONLY):
- Reads ou_pred_vs_empirical.csv from finite_time_fit.py
- Produces:
    fig_ou_var_pooled_validation.(png/pdf)
    fig_ou_D_pooled_validation.(png/pdf)

Main feature:
- Optional legend outside on the right, with FIXED axes geometry across both plots.

Important design choice:
- NO tight_layout()
- NO bbox_inches="tight"
This avoids inconsistent panel sizes when legends/labels differ.

Example
-------
python ./...plot_finite_time_fit_plot.py \
  --in_csv ./results/4-finite_time_validation/unweighted/ou_pred_vs_empirical.csv \
  --out_dir ./figures/4-ou_finite_time_validation/unweighted \
  --fmt png --dpi 300 \
  --figsize_var 4,3 \
  --figsize_D 4,3 \
  --legend_outside_right
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -------------------------
# CLI
# -------------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Plot OU finite-time validation from ou_pred_vs_empirical.csv"
    )

    ap.add_argument("--in_csv", required=True, help="Path to ou_pred_vs_empirical.csv")
    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--csv_sep", default=",", help="CSV separator for input (default: ',')")

    ap.add_argument("--fmt", choices=["png", "pdf"], default="png")
    ap.add_argument("--dpi", type=int, default=300)

    # typography
    ap.add_argument("--font", default="Arial")
    ap.add_argument("--font_size", type=float, default=10.0)

    # figure sizes
    ap.add_argument("--figsize_var", default="4.8,3.2", help="Variance figure size 'w,h'")
    ap.add_argument("--figsize_D", default="4.8,3.2", help="Diffusion figure size 'w,h'")

    # styles
    ap.add_argument("--lw_emp", type=float, default=1.2, help="Line width empirical")
    ap.add_argument("--lw_pred", type=float, default=1.4, help="Line width prediction")
    ap.add_argument("--ms", type=float, default=4.5, help="Marker size")
    ap.add_argument("--alpha_emp_band", type=float, default=0.18, help="Alpha empirical band")
    ap.add_argument("--alpha_pred_band", type=float, default=0.12, help="Alpha prediction band")

    # legend
    ap.add_argument("--legend_loc", default="best",
                    help="Legend location when legend is inside")
    ap.add_argument("--legend_ncol", type=int, default=1)
    ap.add_argument("--legend_frameoff", action="store_true")
    ap.add_argument("--legend_outside_right", action="store_true",
                    help="Place legend outside, on the right, with fixed panel geometry")
    ap.add_argument("--legend_none", action="store_true",
                    help="Do not draw legend")

    # outside-right legend tuning
    ap.add_argument("--legend_out_x", type=float, default=1.02,
                    help="Outside-right legend bbox x")
    ap.add_argument("--legend_out_y", type=float, default=1.0,
                    help="Outside-right legend bbox y")
    ap.add_argument("--legend_borderaxespad", type=float, default=0.0)

    # manual figure margins
    ap.add_argument("--left", type=float, default=0.16)
    ap.add_argument("--right_inside", type=float, default=0.96)
    ap.add_argument("--right_outside", type=float, default=0.72)
    ap.add_argument("--bottom", type=float, default=0.18)
    ap.add_argument("--top", type=float, default=0.88)

    # axes and grid
    ap.add_argument("--grid", action="store_true")
    ap.add_argument("--xlim", default="", help="Optional xlim 'xmin,xmax'")
    ap.add_argument("--ylim_var", default="", help="Optional variance ylim 'ymin,ymax'")
    ap.add_argument("--ylim_D", default="", help="Optional diffusion ylim 'ymin,ymax'")
    ap.add_argument("--title_prefix", default="", help="Optional prefix for titles")

    return ap


# -------------------------
# Helpers
# -------------------------
def parse_pair(s: str) -> Optional[Tuple[float, float]]:
    s = (s or "").strip()
    if not s:
        return None
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 2:
        raise ValueError(f"Expected 'a,b', got: {s}")
    return float(parts[0]), float(parts[1])


def parse_figsize(s: str, default: Tuple[float, float]) -> Tuple[float, float]:
    s = (s or "").strip()
    if not s:
        return default
    parts = [p.strip() for p in s.split(",")]
    if len(parts) != 2:
        raise ValueError(f"Expected figsize 'w,h', got: {s}")
    return float(parts[0]), float(parts[1])


def set_matplotlib_defaults(font: str, font_size: float) -> None:
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


def has_any_finite(a: np.ndarray) -> bool:
    a = np.asarray(a, dtype=float)
    return np.isfinite(a).any()


def apply_common_axes_style(
    ax,
    xlim: Optional[Tuple[float, float]],
    ylim: Optional[Tuple[float, float]],
    grid: bool,
) -> None:
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if grid:
        ax.grid(True)


def add_legend(ax, args: argparse.Namespace) -> None:
    if args.legend_none:
        return

    kwargs = {"ncol": int(args.legend_ncol)}
    if args.legend_frameoff:
        kwargs["frameon"] = False

    if args.legend_outside_right:
        kwargs["loc"] = "upper left"
        kwargs["bbox_to_anchor"] = (float(args.legend_out_x), float(args.legend_out_y))
        kwargs["borderaxespad"] = float(args.legend_borderaxespad)
    else:
        kwargs["loc"] = args.legend_loc

    ax.legend(**kwargs)


def apply_fixed_layout(fig, args: argparse.Namespace) -> None:
    right = args.right_outside if args.legend_outside_right else args.right_inside
    fig.subplots_adjust(
        left=float(args.left),
        right=float(right),
        bottom=float(args.bottom),
        top=float(args.top),
    )


def save_fig(fig, path_base: Path, fmt: str, dpi: int) -> None:
    if fmt == "png":
        fig.savefig(path_base.with_suffix(".png"), dpi=dpi)
    else:
        fig.savefig(path_base.with_suffix(".pdf"))


def plot_validation_panel(
    *,
    x: np.ndarray,
    y_emp: np.ndarray,
    y_pred: np.ndarray,
    emp_lo: Optional[np.ndarray],
    emp_hi: Optional[np.ndarray],
    pred_lo: Optional[np.ndarray],
    pred_hi: Optional[np.ndarray],
    figsize: Tuple[float, float],
    title: str,
    xlabel: str,
    ylabel: str,
    xticks: list[int],
    xlim: Optional[Tuple[float, float]],
    ylim: Optional[Tuple[float, float]],
    args: argparse.Namespace,
    out_path_base: Path,
) -> None:
    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(
        x, y_emp,
        marker="o",
        ms=float(args.ms),
        lw=float(args.lw_emp),
        label="Empirical (MAD)",
    )

    if emp_lo is not None and emp_hi is not None:
        if has_any_finite(emp_lo) and has_any_finite(emp_hi):
            ax.fill_between(x, emp_lo, emp_hi, alpha=float(args.alpha_emp_band))

    ax.plot(
        x, y_pred,
        lw=float(args.lw_pred),
        label="OU prediction",
    )

    if pred_lo is not None and pred_hi is not None:
        if has_any_finite(pred_lo) and has_any_finite(pred_hi):
            ax.fill_between(x, pred_lo, pred_hi, alpha=float(args.alpha_pred_band))

    ax.set_title(f"{args.title_prefix}{title}" if args.title_prefix else title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(xticks)

    apply_common_axes_style(ax=ax, xlim=xlim, ylim=ylim, grid=args.grid)
    add_legend(ax, args)
    apply_fixed_layout(fig, args)
    save_fig(fig, out_path_base, args.fmt, args.dpi)
    plt.close(fig)


# -------------------------
# Main
# -------------------------
def main() -> None:
    args = build_argparser().parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    set_matplotlib_defaults(args.font, args.font_size)

    df = pd.read_csv(Path(args.in_csv).resolve(), sep=args.csv_sep)

    if "dt" not in df.columns:
        raise ValueError("Missing required column 'dt' in input CSV.")

    df["dt"] = pd.to_numeric(df["dt"], errors="coerce")
    df = df[np.isfinite(df["dt"])].copy()
    df["dt"] = df["dt"].astype(int)
    df = df.sort_values("dt").reset_index(drop=True)

    required_cols = ["var_mad_pooled", "D_mad_pooled", "var_ou_pred", "D_ou_pred"]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(
                f"Missing required column '{c}' in input CSV. "
                f"Did you run 4a_ou_finite_time_fit.py?"
            )

    x = df["dt"].to_numpy(int)
    xticks = sorted(np.unique(x).tolist())

    xlim = parse_pair(args.xlim)
    ylim_var = parse_pair(args.ylim_var)
    ylim_D = parse_pair(args.ylim_D)

    figsize_var = parse_figsize(args.figsize_var, (4.8, 3.2))
    figsize_D = parse_figsize(args.figsize_D, (4.8, 3.2))

    # Plot 1: finite-time variance
    plot_validation_panel(
        x=x,
        y_emp=df["var_mad_pooled"].to_numpy(float),
        y_pred=df["var_ou_pred"].to_numpy(float),
        emp_lo=df["var_mad_pooled_lo"].to_numpy(float) if "var_mad_pooled_lo" in df.columns else None,
        emp_hi=df["var_mad_pooled_hi"].to_numpy(float) if "var_mad_pooled_hi" in df.columns else None,
        pred_lo=df["var_ou_pred_lo"].to_numpy(float) if "var_ou_pred_lo" in df.columns else None,
        pred_hi=df["var_ou_pred_hi"].to_numpy(float) if "var_ou_pred_hi" in df.columns else None,
        figsize=figsize_var,
        title="Finite-time variance: empirical vs OU prediction",
        xlabel="Δt (weeks)",
        ylabel=r"$\mathrm{var}^{pooled}_{MAD}(\Delta x)$",
        xticks=xticks,
        xlim=xlim,
        ylim=ylim_var,
        args=args,
        out_path_base=out_dir / "fig_ou_var_pooled_validation",
    )

    # Plot 2: apparent diffusion
    plot_validation_panel(
        x=x,
        y_emp=df["D_mad_pooled"].to_numpy(float),
        y_pred=df["D_ou_pred"].to_numpy(float),
        emp_lo=df["D_mad_pooled_lo"].to_numpy(float) if "D_mad_pooled_lo" in df.columns else None,
        emp_hi=df["D_mad_pooled_hi"].to_numpy(float) if "D_mad_pooled_hi" in df.columns else None,
        pred_lo=df["D_ou_pred_lo"].to_numpy(float) if "D_ou_pred_lo" in df.columns else None,
        pred_hi=df["D_ou_pred_hi"].to_numpy(float) if "D_ou_pred_hi" in df.columns else None,
        figsize=figsize_D,
        title="Apparent diffusion: empirical vs OU prediction",
        xlabel="Δt (weeks)",
        ylabel=r"$D^{pooled}_{MAD}(\Delta t)=\mathrm{var}^{pooled}_{MAD}(\Delta x)/(2\Delta t)$",
        xticks=xticks,
        xlim=xlim,
        ylim=ylim_D,
        args=args,
        out_path_base=out_dir / "fig_ou_D_pooled_validation",
    )

    print("[SAVED]", out_dir / f"fig_ou_var_pooled_validation.{args.fmt}")
    print("[SAVED]", out_dir / f"fig_ou_D_pooled_validation.{args.fmt}")


if __name__ == "__main__":
    main()