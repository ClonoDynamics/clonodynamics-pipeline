#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_powerlaw_qc.py

STEP A — plot-only script for:
(A1) heavy-tail validation and tail-model comparison
(A2) gamma-based replicate QC

Reads
-----
From STEP A analysis:
  - A1_repertoire_tail_fits.csv
  - A1_repertoire_threshold_scan.csv
  - A1_repertoire_model_compare.csv
  - A1_pair_tail_preference_summary.csv
  - A2_pair_gamma_qc.csv
  - A2_gamma_qc_thresholds.csv

Optionally, for panel A CCDF:
  - original repertoire files from --repertoires_dir

Produces
--------
Single figures:
- figA_ccdf_pooled.png
- figB_model_compare_hist.png
- figC_threshold_stability.png
- figD_gamma_scatter.png
- figE_delta_gamma_by_pair.png
- figF_delta_gamma_hist.png
- figG_tail_preference_by_pair.png

Optional multipanel:
- figS1_powerlaw_prior_like.png

Example
-------
python plot_powerlaw_qc.py \
  --data_dir ./results_stepA \
  --fig_dir ./figures_stepA \
  --repertoires_dir ./longitudinali \
  --dpi 300 \
  --fig_w 4.8 \
  --fig_h 3.2 \
  --font_family Arial \
  --font_size 10 \
  --make_multipanel
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PATTERN_DEFAULT = r"^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])(?:\.(?:csv|tsv|txt))?$"


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


def finite_xy(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)
    y = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(float)
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


# ---------------------------------------------------------------------
# Style
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
# Colors
# ---------------------------------------------------------------------
COLOR_POOLED = "black"
COLOR_SINGLE = "0.80"
COLOR_DELTA_LL = "C0"
COLOR_THRESHOLD = "C1"
COLOR_GAMMA = "C2"
COLOR_DELTA_GAMMA = "C3"
COLOR_PAIR_PREF = "C4"
COLOR_PASS = "C2"
COLOR_FAIL = "C3"
COLOR_REF = "black"


# ---------------------------------------------------------------------
# Repertoire reading for panel A
# ---------------------------------------------------------------------
def read_table_auto(path: Path) -> pd.DataFrame:
    last_err = None
    for enc in ["utf-8", "latin1"]:
        for sep in [",", "\t", ";", "|"]:
            try:
                df = pd.read_csv(path, sep=sep, encoding=enc)
                if df.shape[1] > 1:
                    return df
            except Exception as e:
                last_err = e
        try:
            df = pd.read_csv(path, sep=None, engine="python", encoding=enc)
            if df.shape[1] > 1:
                return df
        except Exception as e:
            last_err = e
    raise ValueError(f"Could not parse repertoire file: {path}") from last_err


def read_freqs(path: Path) -> np.ndarray:
    df = read_table_auto(path)
    cols = {c.lower(): c for c in df.columns}

    if "readfraction" in cols:
        x = pd.to_numeric(df[cols["readfraction"]], errors="coerce").to_numpy(dtype=float)
    elif "readcount" in cols:
        counts = pd.to_numeric(df[cols["readcount"]], errors="coerce").to_numpy(dtype=float)
        counts = counts[np.isfinite(counts) & (counts > 0)]
        if counts.size == 0:
            raise ValueError(f"{path.name}: no positive readCount values")
        x = counts / counts.sum()
    else:
        raise ValueError(f"{path.name}: missing readFraction/readCount column")

    x = x[np.isfinite(x) & (x > 0)]
    if x.size == 0:
        raise ValueError(f"{path.name}: no positive frequencies found")
    return np.sort(x)


def list_valid_repertoire_files(repertoires_dir: Path, pattern_str: str) -> List[Path]:
    pat = re.compile(pattern_str)
    valid: List[Path] = []

    for p in sorted(repertoires_dir.iterdir()):
        if not p.is_file():
            continue
        if p.name.startswith("."):
            continue
        if not pat.match(p.name):
            continue
        try:
            _ = read_freqs(p)
            valid.append(p)
        except Exception:
            continue

    return valid


def empirical_ccdf(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(x, dtype=float))
    x = x[np.isfinite(x) & (x > 0)]
    n = len(x)
    if n == 0:
        return np.array([]), np.array([])
    ccdf = (n - np.arange(n)) / n
    return x, ccdf


# ---------------------------------------------------------------------
# QC helpers
# ---------------------------------------------------------------------
def robust_mad(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    med = np.median(x)
    return float(np.median(np.abs(x - med)))


def robust_center_and_mad(x: np.ndarray) -> Tuple[float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan
    med = float(np.median(x))
    mad = robust_mad(x)
    return med, mad


def add_gamma_fail_line(ax: plt.Axes, vals: np.ndarray, cutoff: Optional[float]) -> None:
    if cutoff is None or not np.isfinite(cutoff):
        return
    med, mad = robust_center_and_mad(vals)
    if not np.isfinite(med) or not np.isfinite(mad) or mad < 1e-12:
        return
    ax.axhline(med + cutoff * mad, linewidth=1.0, linestyle="--", color=COLOR_REF)


# ---------------------------------------------------------------------
# Panel functions
# ---------------------------------------------------------------------
def panel_A_ccdf(ax: plt.Axes, files: List[Path]) -> None:
    pooled_list: List[np.ndarray] = []

    for f in files:
        x = read_freqs(f)
        pooled_list.append(x)
        xx, yy = empirical_ccdf(x)
        ax.plot(xx, yy, color=COLOR_SINGLE, alpha=0.20, linewidth=0.7)

    pooled = np.sort(np.concatenate(pooled_list))
    xx, yy = empirical_ccdf(pooled)
    ax.plot(xx, yy, color=COLOR_POOLED, linewidth=1.8, label="Pooled")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Clone frequency")
    ax.set_ylabel("CCDF")


def panel_B_model_compare(ax: plt.Axes, rep_compare: pd.DataFrame) -> None:
    vals = pd.to_numeric(rep_compare["mean_delta_loglik"], errors="coerce").dropna().to_numpy(float)
    if vals.size == 0:
        return

    ax.hist(vals, bins=40, density=True, color=COLOR_DELTA_LL, alpha=0.85)
    ax.axvline(0, linestyle="--", linewidth=1.0, color=COLOR_REF)
    ax.set_xlabel("Mean Δ log-likelihood per clonotype\n(PL − truncated log-normal)")
    ax.set_ylabel("Density")


def panel_C_threshold(ax: plt.Axes, scan_df: pd.DataFrame) -> None:
    df = scan_df.copy()
    df["threshold_quantile"] = pd.to_numeric(df["threshold_quantile"], errors="coerce")
    df["delta"] = pd.to_numeric(df["mean_delta_loglik_pl_minus_tln"], errors="coerce")
    df = df[np.isfinite(df["threshold_quantile"]) & np.isfinite(df["delta"])].copy()
    if df.empty:
        return

    grouped = df.groupby("threshold_quantile")["delta"]
    x = np.array(sorted(grouped.groups.keys()), dtype=float)
    med = grouped.median().reindex(x).to_numpy(dtype=float)
    q1 = grouped.quantile(0.25).reindex(x).to_numpy(dtype=float)
    q3 = grouped.quantile(0.75).reindex(x).to_numpy(dtype=float)

    ax.fill_between(x, q1, q3, alpha=0.22, color=COLOR_THRESHOLD, label="IQR")
    ax.plot(x, med, marker="o", markersize=5.0, linewidth=1.5, color=COLOR_THRESHOLD, label="Median")
    ax.axhline(0, linestyle="--", linewidth=1.0, color=COLOR_REF)

    ax.set_xlabel("Threshold quantile")
    ax.set_ylabel("Mean Δ log-likelihood per clonotype\n(PL − truncated log-normal)")


# ---------------------------------------------------------------------
# Single figure wrappers
# ---------------------------------------------------------------------
def figA_ccdf_pooled(files: List[Path], fig_dir: Path, args: argparse.Namespace) -> None:
    if not files:
        return
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)
    panel_A_ccdf(ax, files)
    if not args.no_titles:
        ax.set_title("Heavy-tailed clone-size structure")
    apply_legend(ax, args)
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figA_ccdf_pooled.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figA_ccdf_pooled.pdf", args)
    plt.close(fig)


def figB_model_compare_hist(rep_compare: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if rep_compare.empty:
        return
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)
    panel_B_model_compare(ax, rep_compare)
    if not args.no_titles:
        ax.set_title("Model comparison in the high-frequency tail")
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figB_model_compare_hist.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figB_model_compare_hist.pdf", args)
    plt.close(fig)


def figC_threshold_stability(scan_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if scan_df.empty:
        return
    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)
    panel_C_threshold(ax, scan_df)
    if not args.no_titles:
        ax.set_title("Stability of model preference across thresholds")
    apply_legend(ax, args)
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figC_threshold_stability.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figC_threshold_stability.pdf", args)
    plt.close(fig)


def figD_gamma_scatter(qc_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if qc_df.empty or "gamma_rep1" not in qc_df.columns or "gamma_rep2" not in qc_df.columns:
        return

    d = qc_df.copy()
    d2 = d[np.isfinite(pd.to_numeric(d["gamma_rep1"], errors="coerce")) & np.isfinite(pd.to_numeric(d["gamma_rep2"], errors="coerce"))].copy()
    if d2.empty:
        return

    x = pd.to_numeric(d2["gamma_rep1"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(d2["gamma_rep2"], errors="coerce").to_numpy(float)

    color_map = {"PASS": COLOR_PASS, "FAIL": COLOR_FAIL}
    colors = [color_map.get(v, COLOR_GAMMA) for v in d2.get("QC_traffic", pd.Series(["PASS"] * len(d2)))]

    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)

    ax.scatter(x, y, s=args.marker_s, alpha=args.marker_alpha, color=colors)
    lo = float(np.nanmin([np.nanmin(x), np.nanmin(y)]))
    hi = float(np.nanmax([np.nanmax(x), np.nanmax(y)]))
    pad = 0.05 * max(hi - lo, 1e-6)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linewidth=args.line_w, color=COLOR_REF, linestyle="--", label="identity")

    ax.set_xlabel(r"$\gamma$ (replicate 1)")
    ax.set_ylabel(r"$\gamma$ (replicate 2)")
    if not args.no_titles:
        ax.set_title("Replicate concordance of the power-law exponent")

    maybe_apply_limits(ax, parse_lim(args.xlim_gamma), parse_lim(args.ylim_gamma))
    apply_legend(ax, args)
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figD_gamma_scatter.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figD_gamma_scatter.pdf", args)
    plt.close(fig)


def figE_delta_gamma_by_pair(qc_df: pd.DataFrame, thr_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if qc_df.empty or "delta_gamma" not in qc_df.columns:
        return

    d = qc_df.copy()
    d = d[np.isfinite(pd.to_numeric(d["delta_gamma"], errors="coerce"))].copy()
    if d.empty:
        return

    d["pair_id"] = d["subject"].astype(str) + "_" + d["time"].astype(str)
    d = d.sort_values(["subject", "time"]).reset_index(drop=True)

    vals = pd.to_numeric(d["delta_gamma"], errors="coerce").to_numpy(float)
    qc = d.get("QC_traffic", pd.Series(["PASS"] * len(d)))
    color_map = {"PASS": COLOR_PASS, "FAIL": COLOR_FAIL}
    colors = [color_map.get(v, COLOR_DELTA_GAMMA) for v in qc]

    fig = plt.figure(figsize=(max(args.fig_w, args.fig_w_pairs), args.fig_h))
    ax = fig.add_subplot(111)

    ax.scatter(np.arange(len(d)), vals, s=args.marker_s, alpha=args.marker_alpha, color=colors)
    ax.plot(np.arange(len(d)), vals, linewidth=args.line_w, color=COLOR_DELTA_GAMMA, alpha=0.7)

    cutoff = None
    if not thr_df.empty and "gamma_mad_cutoff" in thr_df.columns:
        x = pd.to_numeric(thr_df["gamma_mad_cutoff"], errors="coerce").dropna()
        if not x.empty:
            cutoff = float(x.iloc[0])
    add_gamma_fail_line(ax, vals, cutoff)

    ax.set_xticks(np.arange(len(d)))
    ax.set_xticklabels(d["pair_id"].tolist(), rotation=90)
    ax.set_xlabel("Pair")
    ax.set_ylabel(r"$|\Delta \gamma|$")
    if not args.no_titles:
        ax.set_title("Pairwise discordance of the power-law exponent")

    maybe_apply_limits(ax, None, parse_lim(args.ylim_delta_gamma))
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figE_delta_gamma_by_pair.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figE_delta_gamma_by_pair.pdf", args)
    plt.close(fig)


def figF_delta_gamma_hist(qc_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if qc_df.empty or "delta_gamma" not in qc_df.columns:
        return

    vals = pd.to_numeric(qc_df["delta_gamma"], errors="coerce").dropna().to_numpy(float)
    if vals.size == 0:
        return

    fig = plt.figure(figsize=(args.fig_w, args.fig_h))
    ax = fig.add_subplot(111)

    ax.hist(vals, bins=int(args.hist_bins), density=False, color=COLOR_DELTA_GAMMA, alpha=0.85, label="pairs")
    ax.set_xlabel(r"$|\Delta \gamma|$")
    ax.set_ylabel("Count")
    if not args.no_titles:
        ax.set_title("Distribution of pairwise gamma discordance")

    maybe_apply_limits(ax, parse_lim(args.xlim_delta_gamma_hist), parse_lim(args.ylim_delta_gamma_hist))
    apply_legend(ax, args)
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figF_delta_gamma_hist.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figF_delta_gamma_hist.pdf", args)
    plt.close(fig)


def figG_tail_preference_by_pair(pair_pref: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if pair_pref.empty or "mean_frac_thresholds_pl_better" not in pair_pref.columns:
        return

    d = pair_pref.copy()
    d["pair_id"] = d["subject"].astype(str) + "_" + d["time"].astype(str)
    d = d.sort_values(["subject", "time"]).reset_index(drop=True)

    vals = pd.to_numeric(d["mean_frac_thresholds_pl_better"], errors="coerce").to_numpy(float)
    m = np.isfinite(vals)
    if not np.any(m):
        return
    d = d.loc[m].reset_index(drop=True)
    vals = vals[m]

    fig = plt.figure(figsize=(max(args.fig_w, args.fig_w_pairs), args.fig_h))
    ax = fig.add_subplot(111)

    ax.scatter(np.arange(len(d)), vals, s=args.marker_s, alpha=args.marker_alpha, color=COLOR_PAIR_PREF)
    ax.plot(np.arange(len(d)), vals, linewidth=args.line_w, color=COLOR_PAIR_PREF, alpha=0.7)
    ax.axhline(0.5, linewidth=args.refline_w, color=COLOR_REF, linestyle="--")

    ax.set_xticks(np.arange(len(d)))
    ax.set_xticklabels(d["pair_id"].tolist(), rotation=90)
    ax.set_xlabel("Pair")
    ax.set_ylabel("Mean fraction of thresholds with PL better")
    if not args.no_titles:
        ax.set_title("Stability of power-law preference across replicate pairs")

    maybe_apply_limits(ax, None, parse_lim(args.ylim_pair_preference))
    apply_layout(fig, args)
    save_png(fig, fig_dir / "figG_tail_preference_by_pair.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figG_tail_preference_by_pair.pdf", args)
    plt.close(fig)


def figS1_multipanel(files: List[Path], rep_compare: pd.DataFrame, scan_df: pd.DataFrame, fig_dir: Path, args: argparse.Namespace) -> None:
    if not files or rep_compare.empty or scan_df.empty:
        return

    fig, axs = plt.subplots(1, 3, figsize=(12.0, 3.2))
    plt.subplots_adjust(left=0.07, right=0.985, bottom=0.18, top=0.92, wspace=0.38)

    panel_A_ccdf(axs[0], files)
    add_panel_label(axs[0], "A")
    axs[0].legend(frameon=False, loc="lower left")

    panel_B_model_compare(axs[1], rep_compare)
    add_panel_label(axs[1], "B")

    panel_C_threshold(axs[2], scan_df)
    add_panel_label(axs[2], "C")
    axs[2].legend(frameon=False, loc="best")

    save_png(fig, fig_dir / "figS1_powerlaw_prior_like.png", args)
    if args.pdf:
        save_pdf(fig, fig_dir / "figS1_powerlaw_prior_like.pdf", args)
    plt.close(fig)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="STEP A plot-only script for heavy-tail validation and gamma QC.")
    ap.add_argument("--data_dir", required=True, help="Directory containing STEP A analysis CSV outputs.")
    ap.add_argument("--fig_dir", required=True, help="Directory for figures.")

    ap.add_argument("--repertoires_dir", default="", help="Optional directory with original repertoire files for panel A CCDF.")
    ap.add_argument("--pattern", default=PATTERN_DEFAULT, help=f"Regex for repertoire filenames (default: {PATTERN_DEFAULT})")

    # Typography / style
    ap.add_argument("--font_family", default="Arial")
    ap.add_argument("--font_size", type=float, default=10.0)
    ap.add_argument("--fig_w", type=float, default=4.8)
    ap.add_argument("--fig_h", type=float, default=3.2)
    ap.add_argument("--fig_w_pairs", type=float, default=8.0)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--pdf", action="store_true")

    ap.add_argument("--line_w", type=float, default=1.4)
    ap.add_argument("--refline_w", type=float, default=1.0)
    ap.add_argument("--axes_w", type=float, default=1.0)
    ap.add_argument("--tick_w", type=float, default=1.0)
    ap.add_argument("--tick_len", type=float, default=4.0)
    ap.add_argument("--marker_s", type=float, default=24.0)
    ap.add_argument("--marker_alpha", type=float, default=0.9)
    ap.add_argument("--hist_bins", type=int, default=20)

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
    ap.add_argument("--xlim_gamma", default="")
    ap.add_argument("--ylim_gamma", default="")
    ap.add_argument("--ylim_delta_gamma", default="")
    ap.add_argument("--xlim_delta_gamma_hist", default="")
    ap.add_argument("--ylim_delta_gamma_hist", default="")
    ap.add_argument("--ylim_pair_preference", default="")

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

    rep_tail = load_optional_csv(data_dir / "A1_repertoire_tail_fits.csv")
    scan_df = load_optional_csv(data_dir / "A1_repertoire_threshold_scan.csv")
    rep_compare = load_optional_csv(data_dir / "A1_repertoire_model_compare.csv")
    pair_pref = load_optional_csv(data_dir / "A1_pair_tail_preference_summary.csv")
    qc_df = load_optional_csv(data_dir / "A2_pair_gamma_qc.csv")
    thr_df = load_optional_csv(data_dir / "A2_gamma_qc_thresholds.csv")

    files: List[Path] = []
    if args.repertoires_dir:
        rep_dir = Path(args.repertoires_dir)
        if rep_dir.exists():
            files = list_valid_repertoire_files(rep_dir, args.pattern)

    if files:
        figA_ccdf_pooled(files, fig_dir, args)
    figB_model_compare_hist(rep_compare, fig_dir, args)
    figC_threshold_stability(scan_df, fig_dir, args)
    figD_gamma_scatter(qc_df, fig_dir, args)
    figE_delta_gamma_by_pair(qc_df, thr_df, fig_dir, args)
    figF_delta_gamma_hist(qc_df, fig_dir, args)
    figG_tail_preference_by_pair(pair_pref, fig_dir, args)

    if args.make_multipanel and files and (not rep_compare.empty) and (not scan_df.empty):
        figS1_multipanel(files, rep_compare, scan_df, fig_dir, args)

    print("✅ Completed STEP A plotting.")
    print(f"[PATH] fig_dir: {fig_dir.resolve()}")
    if files:
        print("[OUT] figA_ccdf_pooled.png")
    print("[OUT] figB_model_compare_hist.png")
    print("[OUT] figC_threshold_stability.png")
    print("[OUT] figD_gamma_scatter.png")
    print("[OUT] figE_delta_gamma_by_pair.png")
    print("[OUT] figF_delta_gamma_hist.png")
    print("[OUT] figG_tail_preference_by_pair.png")
    if args.make_multipanel and files and (not rep_compare.empty) and (not scan_df.empty):
        print("[OUT] figS1_powerlaw_prior_like.png")


if __name__ == "__main__":
    main()