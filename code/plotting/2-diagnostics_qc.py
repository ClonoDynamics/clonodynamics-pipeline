#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2-diagnostics_qc.py

Additional diagnostic plots for power-law fitting and gamma-based QC.

This script is complementary to:
    1-powerlaw_qc_plot.py

It intentionally generates diagnostic figures that are NOT already produced by
1-powerlaw_qc_plot.py, while matching its overall visual style, typography,
figure sizing logic, legend handling, and export behavior.

Reads
-----
From STEP A analysis:
  - A2_replicate_powerlaw_fits.csv
  - A2_pair_gamma_qc.csv
  - A2_gamma_qc_thresholds.csv        [optional; only for consistency / reporting]
  - A2_summary_by_subject.csv         [optional]
  - A2_summary_overall.csv            [optional]

Produces
--------
Single figures:
- figH_xmin_scatter.png
- figI_ks_scatter.png
- figJ_tail_fraction_scatter.png
- figK_depth_ratio_vs_delta_gamma.png
- figL_depth_ratio_vs_delta_xmin.png
- figM_gamma_mean_vs_depth_ratio.png
- figN_delta_xmin_by_pair.png
- figO_depth_ratio_by_pair.png
- figP_xmin_by_subject.png
- figQ_ks_by_subject.png
- figR_tail_fraction_by_subject.png
- figS_gamma_by_subject.png
- figT_delta_gamma_subject_time.png
- figU_delta_xmin_subject_time.png
- figV_depth_ratio_subject_time.png
- figW_qc_map.png

Optional multipanel:
- figS2_noise_qc_diagnostics.png

Example
-------
python stepA_noise_diagnostics_qc.py \
  --data_dir ./results_stepA \
  --fig_dir ./figures_stepA \
  --dpi 300 \
  --fig_w 4.8 \
  --fig_h 3.2 \
  --font_family Arial \
  --font_size 10 \
  --make_multipanel
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm


# ---------------------------------------------------------------------
# IO / utils
# ---------------------------------------------------------------------
def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def load_optional_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def parse_lim(s: str) -> Optional[Tuple[float, float]]:
    s = (s or "").strip()
    if not s:
        return None
    a, b = s.split(",")
    return float(a), float(b)


def maybe_apply_limits(ax: plt.Axes, xlim: Optional[Tuple[float, float]], ylim: Optional[Tuple[float, float]]) -> None:
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)


# ---------------------------------------------------------------------
# Style (matched to 1-powerlaw_qc_plot.py)
# ---------------------------------------------------------------------
def apply_rcparams(font_family: str, font_size: float, axes_w: float, tick_w: float, tick_len: float, dpi: int) -> None:
    plt.rcParams.update({
        "font.family": font_family,
        "font.size": font_size,
        "axes.labelsize": font_size,
        "axes.titlesize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": max(font_size - 1, 6),
        "axes.linewidth": axes_w,
        "xtick.major.width": tick_w,
        "ytick.major.width": tick_w,
        "xtick.major.size": tick_len,
        "ytick.major.size": tick_len,
        "figure.dpi": 150,
        "savefig.dpi": dpi,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.18, 1.04, label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontweight="bold",
        fontsize=plt.rcParams["font.size"] + 1,
    )


def apply_legend(ax: plt.Axes, args: argparse.Namespace) -> None:
    if args.legend_mode == "none":
        return
    if args.legend_mode == "outside":
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
    if args.legend_mode == "outside":
        fig.tight_layout(rect=[0.0, 0.0, float(args.legend_pad_right), 1.0])
    else:
        fig.tight_layout()


def save_png(fig: plt.Figure, path: Path, args: argparse.Namespace) -> None:
    if args.legend_mode == "outside":
        fig.savefig(path, dpi=int(args.dpi), bbox_inches="tight", pad_inches=float(args.legend_pad_inches))
    else:
        fig.savefig(path, dpi=int(args.dpi))


def save_pdf(fig: plt.Figure, path: Path, args: argparse.Namespace) -> None:
    if args.legend_mode == "outside":
        fig.savefig(path, bbox_inches="tight", pad_inches=float(args.legend_pad_inches))
    else:
        fig.savefig(path)


# ---------------------------------------------------------------------
# Colors (matched / compatible with 1-powerlaw_qc_plot.py)
# ---------------------------------------------------------------------
COLOR_PASS = "C2"
COLOR_FAIL = "C3"
COLOR_NEUTRAL = "0.35"
COLOR_XMIN = "C0"
COLOR_KS = "C1"
COLOR_TAIL = "C4"
COLOR_DEPTH = "C5"
COLOR_GAMMA = "C2"
COLOR_REF = "black"


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def qc_colors(series: pd.Series) -> list[str]:
    cmap = {"PASS": COLOR_PASS, "FAIL": COLOR_FAIL}
    return [cmap.get(str(v), COLOR_NEUTRAL) for v in series]


def finite_frame(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    d = df.copy()
    for c in cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    m = np.ones(len(d), dtype=bool)
    for c in cols:
        m &= np.isfinite(d[c].to_numpy(float))
    return d.loc[m].copy()


def pair_order(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    if "pair_id" not in d.columns:
        d["pair_id"] = d["subject"].astype(str) + "_" + d["time"].astype(str)
    return d.sort_values(["subject", "time"]).reset_index(drop=True)


def long_from_replicates(rep_df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    d = rep_df.copy()
    if value_col not in d.columns:
        return pd.DataFrame()
    out = d[["subject", "replica", value_col]].copy()
    out[value_col] = pd.to_numeric(out[value_col], errors="coerce")
    out = out[np.isfinite(out[value_col])].copy()
    return out


# ---------------------------------------------------------------------
# Panel functions
# ---------------------------------------------------------------------
def panel_identity_scatter(ax: plt.Axes, df: pd.DataFrame, xcol: str, ycol: str,
                           xlabel: str, ylabel: str, args: argparse.Namespace) -> None:
    d = finite_frame(df, [xcol, ycol])
    if d.empty:
        return
    colors = qc_colors(d.get("QC_traffic", pd.Series(["PASS"] * len(d))))
    x = d[xcol].to_numpy(float)
    y = d[ycol].to_numpy(float)
    ax.scatter(x, y, s=args.marker_s, alpha=args.marker_alpha, color=colors)
    lo = float(np.nanmin([np.nanmin(x), np.nanmin(y)]))
    hi = float(np.nanmax([np.nanmax(x), np.nanmax(y)]))
    pad = 0.05 * max(hi - lo, 1e-6)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linewidth=args.line_w,
            color=COLOR_REF, linestyle="--", label="identity")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def panel_scatter(ax: plt.Axes, df: pd.DataFrame, xcol: str, ycol: str,
                  xlabel: str, ylabel: str, color_line: Optional[str], args: argparse.Namespace) -> None:
    d = finite_frame(df, [xcol, ycol])
    if d.empty:
        return
    colors = qc_colors(d.get("QC_traffic", pd.Series(["PASS"] * len(d))))
    x = d[xcol].to_numpy(float)
    y = d[ycol].to_numpy(float)
    ax.scatter(x, y, s=args.marker_s, alpha=args.marker_alpha, color=colors)
    if color_line is not None and len(d) >= 2:
        order = np.argsort(x)
        ax.plot(x[order], y[order], linewidth=args.line_w, color=color_line, alpha=0.55)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def panel_pair_series(ax: plt.Axes, df: pd.DataFrame, value_col: str, ylabel: str,
                      color: str, args: argparse.Namespace) -> None:
    d = finite_frame(pair_order(df), [value_col])
    if d.empty:
        return
    vals = d[value_col].to_numpy(float)
    colors = qc_colors(d.get("QC_traffic", pd.Series(["PASS"] * len(d))))
    xpos = np.arange(len(d))
    ax.scatter(xpos, vals, s=args.marker_s, alpha=args.marker_alpha, color=colors)
    ax.plot(xpos, vals, linewidth=args.line_w, color=color, alpha=0.7)
    ax.set_xticks(xpos)
    ax.set_xticklabels(d["pair_id"].tolist(), rotation=90)
    ax.set_xlabel("Pair")
    ax.set_ylabel(ylabel)


def panel_stripbox(ax: plt.Axes, rep_df: pd.DataFrame, value_col: str, ylabel: str, color: str) -> None:
    d = long_from_replicates(rep_df, value_col)
    if d.empty:
        return
    subjects = sorted(d["subject"].astype(str).unique(), key=lambda s: int(s) if s.isdigit() else 10**9)
    data = [d.loc[d["subject"].astype(str) == sbj, value_col].dropna().to_numpy(float) for sbj in subjects]
    x = np.arange(1, len(subjects) + 1)
    ax.boxplot(data, positions=x, widths=0.6, patch_artist=False, showfliers=False)
    for i, vals in enumerate(data, start=1):
        if len(vals):
            rng = np.random.default_rng(12345 + i)
            jitter = rng.uniform(-0.12, 0.12, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals, s=12, alpha=0.8, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels(subjects)
    ax.set_xlabel("Subject")
    ax.set_ylabel(ylabel)


def panel_heatmap(ax: plt.Axes, df: pd.DataFrame, value_col: str, cbar_label: str) -> None:
    d = finite_frame(df, [value_col])
    if d.empty:
        return
    mat = (
        d.pivot_table(index="subject", columns="time", values=value_col, aggfunc="median")
         .sort_index(axis=0)
         .sort_index(axis=1)
    )
    im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto", interpolation="nearest", cmap="viridis")
    ax.set_xlabel("time")
    ax.set_ylabel("subject")
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels([str(c) for c in mat.columns])
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels([str(i) for i in mat.index])
    cbar = ax.figure.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(cbar_label)


def panel_qc_heatmap(ax: plt.Axes, df: pd.DataFrame) -> None:
    if "QC_traffic" not in df.columns:
        return
    d = pair_order(df)
    code_map = {"FAIL": 0, "PASS": 1}
    d["qc_code"] = d["QC_traffic"].map(code_map)
    d = d[d["qc_code"].notna()].copy()
    if d.empty:
        return
    mat = (
        d.pivot_table(index="subject", columns="time", values="qc_code", aggfunc="max")
         .sort_index(axis=0)
         .sort_index(axis=1)
    )
    cmap = ListedColormap([COLOR_FAIL, COLOR_PASS])
    norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)
    im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    ax.set_xlabel("time")
    ax.set_ylabel("subject")
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels([str(c) for c in mat.columns])
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels([str(i) for i in mat.index])
    cbar = ax.figure.colorbar(im, ax=ax, pad=0.02)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["FAIL", "PASS"])


# ---------------------------------------------------------------------
# Figure wrappers
# ---------------------------------------------------------------------
def make_single(fig_name: str, fig_dir: Path, args: argparse.Namespace, draw_fn) -> None:
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)
    draw_fn(ax)
    apply_layout(fig, args)
    save_png(fig, fig_dir / f"{fig_name}.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / f"{fig_name}.pdf", args)
    plt.close(fig)


def figH_xmin_scatter(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"xmin_rep1", "xmin_rep2"}.issubset(pair_df.columns):
        return
    def draw(ax):
        panel_identity_scatter(ax, pair_df, "xmin_rep1", "xmin_rep2", r"$x_{\min}$ (replicate 1)", r"$x_{\min}$ (replicate 2)", args)
        if not args.no_titles:
            ax.set_title(r"Replicate concordance of $x_{\min}$")
        maybe_apply_limits(ax, parse_lim(args.xlim_xmin), parse_lim(args.ylim_xmin))
        apply_legend(ax, args)
    make_single("figH_xmin_scatter", fig_dir, args, draw)


def figI_ks_scatter(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"ks_rep1", "ks_rep2"}.issubset(pair_df.columns):
        return
    def draw(ax):
        panel_identity_scatter(ax, pair_df, "ks_rep1", "ks_rep2", "KS distance (replicate 1)", "KS distance (replicate 2)", args)
        if not args.no_titles:
            ax.set_title("Replicate concordance of KS distance")
        maybe_apply_limits(ax, parse_lim(args.xlim_ks), parse_lim(args.ylim_ks))
        apply_legend(ax, args)
    make_single("figI_ks_scatter", fig_dir, args, draw)


def figJ_tail_fraction_scatter(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"tail_fraction_rep1", "tail_fraction_rep2"}.issubset(pair_df.columns):
        return
    def draw(ax):
        panel_identity_scatter(ax, pair_df, "tail_fraction_rep1", "tail_fraction_rep2", "Tail fraction (replicate 1)", "Tail fraction (replicate 2)", args)
        if not args.no_titles:
            ax.set_title("Replicate concordance of tail fraction")
        maybe_apply_limits(ax, parse_lim(args.xlim_tail_fraction), parse_lim(args.ylim_tail_fraction))
        apply_legend(ax, args)
    make_single("figJ_tail_fraction_scatter", fig_dir, args, draw)


def figK_depth_ratio_vs_delta_gamma(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"depth_ratio", "delta_gamma"}.issubset(pair_df.columns):
        return
    def draw(ax):
        panel_scatter(ax, pair_df, "depth_ratio", "delta_gamma", "Depth ratio", r"$|\Delta \gamma|$", COLOR_DEPTH, args)
        if not args.no_titles:
            ax.set_title(r"Depth imbalance versus $|\Delta \gamma|$")
        maybe_apply_limits(ax, parse_lim(args.xlim_depth_ratio), parse_lim(args.ylim_delta_gamma))
    make_single("figK_depth_ratio_vs_delta_gamma", fig_dir, args, draw)


def figL_depth_ratio_vs_delta_xmin(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"depth_ratio", "delta_xmin"}.issubset(pair_df.columns):
        return
    def draw(ax):
        panel_scatter(ax, pair_df, "depth_ratio", "delta_xmin", "Depth ratio", r"$|\Delta x_{\min}|$", COLOR_XMIN, args)
        if not args.no_titles:
            ax.set_title(r"Depth imbalance versus $|\Delta x_{\min}|$")
        maybe_apply_limits(ax, parse_lim(args.xlim_depth_ratio), parse_lim(args.ylim_delta_xmin))
    make_single("figL_depth_ratio_vs_delta_xmin", fig_dir, args, draw)


def figM_gamma_mean_vs_depth_ratio(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"depth_ratio", "gamma_mean"}.issubset(pair_df.columns):
        return
    def draw(ax):
        panel_scatter(ax, pair_df, "depth_ratio", "gamma_mean", "Depth ratio", r"Mean $\gamma$", COLOR_GAMMA, args)
        if not args.no_titles:
            ax.set_title(r"Depth imbalance versus mean $\gamma$")
        maybe_apply_limits(ax, parse_lim(args.xlim_depth_ratio), parse_lim(args.ylim_gamma_mean))
    make_single("figM_gamma_mean_vs_depth_ratio", fig_dir, args, draw)


def figN_delta_xmin_by_pair(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or "delta_xmin" not in pair_df.columns:
        return
    fig = plt.figure(figsize=(max(args.fig_w, args.fig_w_pairs), args.fig_h))
    ax = fig.add_subplot(111)
    panel_pair_series(ax, pair_df, "delta_xmin", r"$|\Delta x_{\min}|$", COLOR_XMIN, args)
    if not args.no_titles:
        ax.set_title(r"Pairwise discordance of $x_{\min}$")
    maybe_apply_limits(ax, None, parse_lim(args.ylim_delta_xmin))
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figN_delta_xmin_by_pair.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figN_delta_xmin_by_pair.pdf", args)
    plt.close(fig)


def figO_depth_ratio_by_pair(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or "depth_ratio" not in pair_df.columns:
        return
    fig = plt.figure(figsize=(max(args.fig_w, args.fig_w_pairs), args.fig_h))
    ax = fig.add_subplot(111)
    panel_pair_series(ax, pair_df, "depth_ratio", "Depth ratio", COLOR_DEPTH, args)
    if not args.no_titles:
        ax.set_title("Pairwise sequencing depth imbalance")
    maybe_apply_limits(ax, None, parse_lim(args.ylim_depth_ratio))
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figO_depth_ratio_by_pair.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figO_depth_ratio_by_pair.pdf", args)
    plt.close(fig)


def figP_xmin_by_subject(rep_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if rep_df.empty or "xmin_hat" not in rep_df.columns:
        return
    def draw(ax):
        panel_stripbox(ax, rep_df, "xmin_hat", r"$x_{\min}$", COLOR_XMIN)
        if not args.no_titles:
            ax.set_title(r"Distribution of $x_{\min}$ by subject")
    make_single("figP_xmin_by_subject", fig_dir, args, draw)


def figQ_ks_by_subject(rep_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if rep_df.empty or "ks_distance" not in rep_df.columns:
        return
    def draw(ax):
        panel_stripbox(ax, rep_df, "ks_distance", "KS distance", COLOR_KS)
        if not args.no_titles:
            ax.set_title("Distribution of KS distance by subject")
    make_single("figQ_ks_by_subject", fig_dir, args, draw)


def figR_tail_fraction_by_subject(rep_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if rep_df.empty or "tail_fraction" not in rep_df.columns:
        return
    def draw(ax):
        panel_stripbox(ax, rep_df, "tail_fraction", "Tail fraction", COLOR_TAIL)
        if not args.no_titles:
            ax.set_title("Distribution of tail fraction by subject")
    make_single("figR_tail_fraction_by_subject", fig_dir, args, draw)


def figS_gamma_by_subject(rep_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if rep_df.empty or "gamma_hat" not in rep_df.columns:
        return
    def draw(ax):
        panel_stripbox(ax, rep_df, "gamma_hat", r"$\gamma$", COLOR_GAMMA)
        if not args.no_titles:
            ax.set_title(r"Distribution of replicate-level $\gamma$ by subject")
    make_single("figS_gamma_by_subject", fig_dir, args, draw)


def figT_delta_gamma_subject_time(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"subject", "time", "delta_gamma"}.issubset(pair_df.columns):
        return
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)
    panel_heatmap(ax, pair_df, "delta_gamma", r"$|\Delta \gamma|$")
    if not args.no_titles:
        ax.set_title(r"Subject x time map of $|\Delta \gamma|$")
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figT_delta_gamma_subject_time.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figT_delta_gamma_subject_time.pdf", args)
    plt.close(fig)


def figU_delta_xmin_subject_time(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"subject", "time", "delta_xmin"}.issubset(pair_df.columns):
        return
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)
    panel_heatmap(ax, pair_df, "delta_xmin", r"$|\Delta x_{\min}|$")
    if not args.no_titles:
        ax.set_title(r"Subject x time map of $|\Delta x_{\min}|$")
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figU_delta_xmin_subject_time.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figU_delta_xmin_subject_time.pdf", args)
    plt.close(fig)


def figV_depth_ratio_subject_time(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or not {"subject", "time", "depth_ratio"}.issubset(pair_df.columns):
        return
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)
    panel_heatmap(ax, pair_df, "depth_ratio", "Depth ratio")
    if not args.no_titles:
        ax.set_title("Subject x time map of depth ratio")
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figV_depth_ratio_subject_time.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figV_depth_ratio_subject_time.pdf", args)
    plt.close(fig)


def figW_qc_map(pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_df.empty or "QC_traffic" not in pair_df.columns:
        return
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)
    panel_qc_heatmap(ax, pair_df)
    if not args.no_titles:
        ax.set_title("Subject x time QC map")
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figW_qc_map.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figW_qc_map.pdf", args)
    plt.close(fig)


def figS2_multipanel(rep_df: pd.DataFrame, pair_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    needed_pair = {"xmin_rep1", "xmin_rep2", "depth_ratio", "delta_gamma", "subject", "time", "QC_traffic"}
    needed_rep = {"subject", "xmin_hat"}
    if rep_df.empty or pair_df.empty or not needed_pair.issubset(pair_df.columns) or not needed_rep.issubset(rep_df.columns):
        return

    fig, axs = plt.subplots(2, 2, figsize=(9.6, 6.4))
    plt.subplots_adjust(left=0.08, right=0.98, bottom=0.10, top=0.95, wspace=0.35, hspace=0.40)

    panel_identity_scatter(axs[0, 0], pair_df, "xmin_rep1", "xmin_rep2",
                           r"$x_{\min}$ (replicate 1)", r"$x_{\min}$ (replicate 2)", args)
    add_panel_label(axs[0, 0], "A")
    axs[0, 0].legend(frameon=False, loc="best")

    panel_scatter(axs[0, 1], pair_df, "depth_ratio", "delta_gamma", "Depth ratio", r"$|\Delta \gamma|$", COLOR_DEPTH, args)
    add_panel_label(axs[0, 1], "B")

    panel_stripbox(axs[1, 0], rep_df, "xmin_hat", r"$x_{\min}$", COLOR_XMIN)
    add_panel_label(axs[1, 0], "C")

    panel_qc_heatmap(axs[1, 1], pair_df)
    add_panel_label(axs[1, 1], "D")

    save_png(fig, fig_dir / "figS2_noise_qc_diagnostics.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figS2_noise_qc_diagnostics.pdf", args)
    plt.close(fig)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="STEP A diagnostic plots complementary to 1-powerlaw_qc_plot.py.")
    ap.add_argument("--data_dir", required=True, help="Directory containing STEP A analysis CSV outputs.")
    ap.add_argument("--fig_dir", required=True, help="Directory for figures.")

    # Typography / style (matched to 1-powerlaw_qc_plot.py)
    ap.add_argument("--font_family", default="Arial")
    ap.add_argument("--font_size", type=float, default=10.0)
    ap.add_argument("--fig_w", type=float, default=4.8)
    ap.add_argument("--fig_h", type=float, default=3.2)
    ap.add_argument("--fig_w_pairs", type=float, default=8.0)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--pdf", action="store_true")

    ap.add_argument("--line_w", type=float, default=1.4)
    ap.add_argument("--axes_w", type=float, default=1.0)
    ap.add_argument("--tick_w", type=float, default=1.0)
    ap.add_argument("--tick_len", type=float, default=4.0)
    ap.add_argument("--marker_s", type=float, default=24.0)
    ap.add_argument("--marker_alpha", type=float, default=0.9)

    # Legends
    ap.add_argument("--legend_mode", choices=["inside", "outside", "none"], default="inside")
    ap.add_argument("--legend_loc", default="best")
    ap.add_argument("--legend_frame", type=int, default=0)
    ap.add_argument("--legend_fontsize", type=float, default=10.0)
    ap.add_argument("--legend_bbox_x", type=float, default=1.02)
    ap.add_argument("--legend_bbox_y", type=float, default=1.0)
    ap.add_argument("--legend_pad_right", type=float, default=0.78)
    ap.add_argument("--legend_pad_inches", type=float, default=0.05)

    # Titles
    ap.add_argument("--no_titles", action="store_true")

    # Axis limits
    ap.add_argument("--xlim_xmin", default="")
    ap.add_argument("--ylim_xmin", default="")
    ap.add_argument("--xlim_ks", default="")
    ap.add_argument("--ylim_ks", default="")
    ap.add_argument("--xlim_tail_fraction", default="")
    ap.add_argument("--ylim_tail_fraction", default="")
    ap.add_argument("--xlim_depth_ratio", default="")
    ap.add_argument("--ylim_depth_ratio", default="")
    ap.add_argument("--ylim_delta_gamma", default="")
    ap.add_argument("--ylim_delta_xmin", default="")
    ap.add_argument("--ylim_gamma_mean", default="")

    # Multipanel
    ap.add_argument("--make_multipanel", action="store_true")

    return ap.parse_args()


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    fig_dir = ensure_dir(Path(args.fig_dir))

    apply_rcparams(args.font_family, args.font_size, args.axes_w, args.tick_w, args.tick_len, args.dpi)

    rep_df = load_optional_csv(data_dir / "A2_replicate_powerlaw_fits.csv")
    pair_df = load_optional_csv(data_dir / "A2_pair_gamma_qc.csv")
    _thr_df = load_optional_csv(data_dir / "A2_gamma_qc_thresholds.csv")
    _summary_subj = load_optional_csv(data_dir / "A2_summary_by_subject.csv")
    _summary_overall = load_optional_csv(data_dir / "A2_summary_overall.csv")

    figH_xmin_scatter(pair_df, fig_dir, args)
    figI_ks_scatter(pair_df, fig_dir, args)
    figJ_tail_fraction_scatter(pair_df, fig_dir, args)
    figK_depth_ratio_vs_delta_gamma(pair_df, fig_dir, args)
    figL_depth_ratio_vs_delta_xmin(pair_df, fig_dir, args)
    figM_gamma_mean_vs_depth_ratio(pair_df, fig_dir, args)
    figN_delta_xmin_by_pair(pair_df, fig_dir, args)
    figO_depth_ratio_by_pair(pair_df, fig_dir, args)
    figP_xmin_by_subject(rep_df, fig_dir, args)
    figQ_ks_by_subject(rep_df, fig_dir, args)
    figR_tail_fraction_by_subject(rep_df, fig_dir, args)
    figS_gamma_by_subject(rep_df, fig_dir, args)
    figT_delta_gamma_subject_time(pair_df, fig_dir, args)
    figU_delta_xmin_subject_time(pair_df, fig_dir, args)
    figV_depth_ratio_subject_time(pair_df, fig_dir, args)
    figW_qc_map(pair_df, fig_dir, args)

    if args.make_multipanel:
        figS2_multipanel(rep_df, pair_df, fig_dir, args)

    print("✅ Completed STEP A additional diagnostics plotting.")
    print(f"[PATH] fig_dir: {fig_dir.resolve()}")
    for name in [
        "figH_xmin_scatter.png",
        "figI_ks_scatter.png",
        "figJ_tail_fraction_scatter.png",
        "figK_depth_ratio_vs_delta_gamma.png",
        "figL_depth_ratio_vs_delta_xmin.png",
        "figM_gamma_mean_vs_depth_ratio.png",
        "figN_delta_xmin_by_pair.png",
        "figO_depth_ratio_by_pair.png",
        "figP_xmin_by_subject.png",
        "figQ_ks_by_subject.png",
        "figR_tail_fraction_by_subject.png",
        "figS_gamma_by_subject.png",
        "figT_delta_gamma_subject_time.png",
        "figU_delta_xmin_subject_time.png",
        "figV_depth_ratio_subject_time.png",
        "figW_qc_map.png",
    ]:
        print(f"[OUT] {name}")
    if args.make_multipanel:
        print("[OUT] figS2_noise_qc_diagnostics.png")


if __name__ == "__main__":
    main()
