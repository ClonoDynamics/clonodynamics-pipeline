#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_transitions_alpha_summary.py

Publication-ready rewrite of transitions_alpha_summary.py.

What this script does
---------------------
1) Reads one or more transition CSV files produced at different alpha thresholds.
2) Infers alpha from:
   - explicit CLI arguments, or
   - parent folder names like p_02, or
   - filename tokens like p_02.
3) Infers transition design from filename:
   - adjacent
   - lag_sampled
   - all_transitions
4) Counts transition classes (TT, TF, FT, FF) globally.
5) Stratifies transitions by x0 bins based on the FIRST time-point frequency.
6) Produces publication-ready summary plots with configurable aesthetics.

Key improvements over the original script
-----------------------------------------
- Full CLI control of figure size in inches and dpi.
- Full CLI control of typography (Arial 10 by default).
- Configurable line widths, marker size, marker alpha, axes width, tick width/length.
- Configurable legend mode (inside / outside / none).
- Configurable panel titles.
- Configurable x/y limits for all major plot families.
- Optional PDF export.
- Consistent styling across all figures.
- Safe handling of x0 scale (raw, ln, log10, or auto inference).

Expected input columns
----------------------
Required:
    obs_class

Optional:
    x0
    freq_t0
    f0
    dt
    subject

Main outputs
------------
Data tables:
    transition_class_summary_by_alpha.csv
    transition_total_by_alpha.csv
    transition_class_summary_by_alpha_x0bin.csv
    x0_fixed_bins.csv

Figures:
    obsclass_counts_vs_alpha_<design>.png
    obsclass_props_vs_alpha_<design>.png
    obsclass_props_stacked_<design>.png
    TT_count_vs_alpha_<design>.png
    TT_prop_vs_alpha_<design>.png
    TT_retention_vs_alpha_<design>.png
    TT_prop_by_x0_vs_alpha_<design>.png
    TT_count_by_x0_vs_alpha_<design>.png
    TT_count_by_x0_vs_alpha_<design>_logy.png

Example
-------
python plot_transitions_alpha_summary.py \
  --glob "results/p_*/transitions_*.csv" \
  --out-data-dir ./tables_alpha \
  --out-fig-dir ./figures_alpha \
  --font Arial \
  --fontsize 10 \
  --fig-w 5.8 \
  --fig-h 3.8 \
  --fig-w-wide 7.0 \
  --fig-h-wide 4.2 \
  --dpi 300 \
  --pdf \
  --legend-mode outside
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator


CLASS_ORDER = ["TT", "TF", "FT", "FF"]
DEFAULT_CLASS_COLORS = {
    "TT": "C0",
    "TF": "C1",
    "FT": "C2",
    "FF": "C3",
}


# =========================================================
# Helpers
# =========================================================

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def set_mpl_defaults(args: argparse.Namespace) -> None:
    plt.rcParams.update({
        "font.family": args.font,
        "font.size": args.fontsize,
        "axes.titlesize": args.fontsize,
        "axes.labelsize": args.fontsize,
        "xtick.labelsize": args.fontsize,
        "ytick.labelsize": args.fontsize,
        "legend.fontsize": args.legend_fontsize,
        "axes.linewidth": args.axes_w,
        "xtick.major.width": args.tick_w,
        "ytick.major.width": args.tick_w,
        "xtick.major.size": args.tick_len,
        "ytick.major.size": args.tick_len,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": args.dpi,
    })


def format_alpha_label(v: float) -> str:
    if not np.isfinite(v):
        return ""
    return f"{v:.6f}".rstrip("0").rstrip(".")


def alpha_from_token(token: str) -> float:
    token = token.strip()
    if "." in token:
        return float(token)
    if not re.fullmatch(r"\d+", token):
        raise ValueError(f"Invalid alpha token: {token}")
    if token == "1":
        return 1.0
    return int(token) / (10 ** len(token))


def infer_alpha_from_path(
    fp: Path,
    cli_alpha_label: Optional[str] = None,
    cli_alpha_value: Optional[float] = None,
) -> tuple[str, float]:
    if cli_alpha_label is not None and cli_alpha_value is not None:
        return str(cli_alpha_label), float(cli_alpha_value)

    parent = fp.parent.name
    m_parent = re.search(r"(p_(?P<token>\d+|\d*\.\d+))", parent)
    if m_parent:
        alpha_label = m_parent.group(1)
        alpha_float = alpha_from_token(m_parent.group("token"))
        return alpha_label, float(alpha_float)

    name = fp.name
    m_name = re.search(r"(p_(?P<token>\d+|\d*\.\d+))", name)
    if m_name:
        alpha_label = m_name.group(1)
        alpha_float = alpha_from_token(m_name.group("token"))
        return alpha_label, float(alpha_float)

    raise ValueError(
        f"Cannot infer alpha for file: {fp}\n"
        f"Expected one of:\n"
        f"  - parent folder like p_02\n"
        f"  - filename containing p_02\n"
        f"  - explicit CLI arguments --alpha-label and --alpha-value"
    )


def infer_design_from_filename(fp: Path) -> str:
    name = fp.name.lower()
    if "adj" in name:
        return "adjacent"
    if "lag" in name:
        return "lag_sampled"
    if "all" in name:
        return "all_transitions"
    return "unknown"


def interval_labels(edges: List[float]) -> List[str]:
    return [f"[{edges[i]}, {edges[i+1]})" for i in range(len(edges) - 1)]


def interval_labels_with_tails(edges: List[float]) -> List[str]:
    labels = [f"< {edges[0]}"]
    labels.extend(interval_labels(edges))
    labels.append(f">= {edges[-1]}")
    return labels


def apply_alpha_axis(ax: plt.Axes, x: np.ndarray) -> None:
    ax.set_xscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([format_alpha_label(v) for v in x], rotation=45, ha="right")
    ax.xaxis.set_minor_locator(NullLocator())


def parse_lim(s: str) -> Optional[Tuple[float, float]]:
    s = (s or "").strip()
    if not s:
        return None
    a, b = s.split(",")
    return float(a), float(b)


def maybe_apply_limits(
    ax: plt.Axes,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
) -> None:
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)


def apply_legend(ax: plt.Axes, args: argparse.Namespace, title: Optional[str] = None) -> None:
    if args.legend_mode == "none":
        return
    kwargs = {
        "frameon": bool(args.legend_frame),
        "fontsize": args.legend_fontsize,
        "title": title,
        "title_fontsize": args.legend_fontsize,
    }
    if args.legend_mode == "outside":
        ax.legend(
            loc=args.legend_loc,
            bbox_to_anchor=(args.legend_bbox_x, args.legend_bbox_y),
            borderaxespad=0.0,
            **kwargs,
        )
    else:
        ax.legend(loc=args.legend_loc, **kwargs)


def apply_layout(fig: plt.Figure, args: argparse.Namespace) -> None:
    if args.legend_mode == "outside":
        fig.tight_layout(rect=[0.0, 0.0, args.legend_pad_right, 1.0])
    else:
        fig.tight_layout()


def save_figure(fig: plt.Figure, path_png: Path, args: argparse.Namespace) -> None:
    path_png.parent.mkdir(parents=True, exist_ok=True)
    if args.legend_mode == "outside":
        fig.savefig(path_png, dpi=args.dpi, bbox_inches="tight", pad_inches=args.legend_pad_inches)
    else:
        fig.savefig(path_png, dpi=args.dpi, bbox_inches="tight")
    print(f"[SAVED FIGURE] {path_png}")

    if args.pdf:
        pdf_path = path_png.with_suffix(".pdf")
        if args.legend_mode == "outside":
            fig.savefig(pdf_path, bbox_inches="tight", pad_inches=args.legend_pad_inches)
        else:
            fig.savefig(pdf_path, bbox_inches="tight")
        print(f"[SAVED FIGURE] {pdf_path}")

    plt.close(fig)


# =========================================================
# x0 handling
# =========================================================

def _safe_log10(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    s = s.where(s > 0)
    return np.log10(s)


def _infer_x0_scale_from_series(x: pd.Series) -> str:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if x.empty:
        return "unknown"

    q01, q50, q99 = x.quantile([0.01, 0.50, 0.99]).tolist()

    if q01 >= 0 and q99 <= 1.0:
        return "raw"
    if q99 <= 0 and q01 >= -15 and q50 >= -8.5:
        return "log10"
    if q99 <= 0 and q01 < -10:
        return "ln"
    return "log10"


def compute_x0_log10(df: pd.DataFrame, x0_mode: str = "auto") -> Tuple[pd.Series, str]:
    if "freq_t0" in df.columns:
        return _safe_log10(df["freq_t0"]), "freq_t0(raw)->log10"

    if "f0" in df.columns:
        return _safe_log10(df["f0"]), "f0(raw)->log10"

    if "x0" not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float), "missing"

    x = pd.to_numeric(df["x0"], errors="coerce")

    mode = _infer_x0_scale_from_series(x) if x0_mode == "auto" else x0_mode

    if mode == "raw":
        return _safe_log10(x), "x0(raw)->log10"
    if mode == "ln":
        return x / math.log(10.0), "x0(ln)->log10"
    if mode == "log10":
        return x, "x0(log10)"

    return x, f"x0({mode})"


# =========================================================
# Reading
# =========================================================

def read_transition_csv(
    fp: Path,
    x0_mode: str = "auto",
    cli_alpha_label: Optional[str] = None,
    cli_alpha_value: Optional[float] = None,
) -> pd.DataFrame:
    df = pd.read_csv(fp)

    if "obs_class" not in df.columns:
        raise ValueError(f"{fp}: required column 'obs_class' not found.")

    alpha_label, alpha_float = infer_alpha_from_path(
        fp,
        cli_alpha_label=cli_alpha_label,
        cli_alpha_value=cli_alpha_value,
    )
    design = infer_design_from_filename(fp)

    out = df.copy()
    out["source_file"] = str(fp)
    out["source_name"] = fp.name
    out["alpha_label"] = alpha_label
    out["alpha_float"] = alpha_float
    out["design"] = design
    out["obs_class"] = out["obs_class"].astype(str).str.strip()

    x0_log10, x0_source = compute_x0_log10(out, x0_mode=x0_mode)
    out["x0_log10"] = x0_log10
    out["x0_source"] = x0_source

    if "dt" in out.columns:
        out["dt"] = pd.to_numeric(out["dt"], errors="coerce")
    if "subject" in out.columns:
        out["subject"] = out["subject"].astype(str).str.strip()

    return out


def load_inputs(
    paths: List[Path],
    x0_mode: str = "auto",
    cli_alpha_label: Optional[str] = None,
    cli_alpha_value: Optional[float] = None,
) -> pd.DataFrame:
    dfs = [
        read_transition_csv(
            fp,
            x0_mode=x0_mode,
            cli_alpha_label=cli_alpha_label,
            cli_alpha_value=cli_alpha_value,
        )
        for fp in paths
    ]
    df = pd.concat(dfs, ignore_index=True)

    print(f"[INFO] Loaded {len(df)} rows from {len(paths)} files")
    print(f"[INFO] Alpha levels: {sorted(df['alpha_label'].dropna().unique().tolist())}")
    print(f"[INFO] Designs: {sorted(df['design'].dropna().unique().tolist())}")
    print(f"[INFO] Obs classes found: {sorted(df['obs_class'].dropna().unique().tolist())}")
    print(f"[INFO] x0 sources found: {sorted(df['x0_source'].dropna().unique().tolist())}")
    return df


# =========================================================
# Summaries
# =========================================================

def build_global_summary(df: pd.DataFrame) -> pd.DataFrame:
    alphas = (
        df[["alpha_label", "alpha_float"]]
        .drop_duplicates()
        .sort_values(["alpha_float", "alpha_label"])
        .reset_index(drop=True)
    )
    designs = sorted(df["design"].dropna().unique().tolist())

    rows = []
    for design in designs:
        gd = df[df["design"] == design].copy()
        for _, ar in alphas.iterrows():
            a_lab = ar["alpha_label"]
            a_val = ar["alpha_float"]
            ga = gd[gd["alpha_label"] == a_lab].copy()

            n_total = len(ga)
            counts = ga["obs_class"].value_counts(dropna=False).to_dict()

            for c in CLASS_ORDER:
                n = int(counts.get(c, 0))
                prop = n / n_total if n_total > 0 else np.nan
                rows.append({
                    "alpha_label": a_lab,
                    "alpha_float": a_val,
                    "design": design,
                    "obs_class": c,
                    "n": n,
                    "n_total": n_total,
                    "prop": prop,
                })

    return pd.DataFrame(rows).sort_values(
        ["design", "alpha_float", "obs_class"]
    ).reset_index(drop=True)


def assign_x0_bins(df: pd.DataFrame, edges: List[float], include_tails: bool = True) -> pd.DataFrame:
    if "x0_log10" not in df.columns:
        return pd.DataFrame()

    out = df.dropna(subset=["x0_log10"]).copy()
    if out.empty:
        return out

    if include_tails:
        cut_edges = [-np.inf] + list(edges) + [np.inf]
        labels = interval_labels_with_tails(edges)
    else:
        cut_edges = list(edges)
        labels = interval_labels(edges)

    out["x0_bin"] = pd.cut(
        out["x0_log10"],
        bins=cut_edges,
        right=False,
        include_lowest=True,
        labels=labels,
    )
    out = out.dropna(subset=["x0_bin"]).copy()
    out["x0_bin_label"] = out["x0_bin"].astype(str)

    label_to_order = {lab: i for i, lab in enumerate(labels)}
    out["x0_bin_order"] = out["x0_bin_label"].map(label_to_order)
    return out


def build_x0_summary(df_binned: pd.DataFrame, all_bin_labels: List[str]) -> pd.DataFrame:
    if df_binned.empty:
        return pd.DataFrame()

    alphas = (
        df_binned[["alpha_label", "alpha_float"]]
        .drop_duplicates()
        .sort_values(["alpha_float", "alpha_label"])
        .reset_index(drop=True)
    )
    designs = sorted(df_binned["design"].dropna().unique().tolist())

    rows = []
    for design in designs:
        gd = df_binned[df_binned["design"] == design].copy()
        for _, ar in alphas.iterrows():
            a_lab = ar["alpha_label"]
            a_val = ar["alpha_float"]
            ga = gd[gd["alpha_label"] == a_lab].copy()

            for bin_label in all_bin_labels:
                gb = ga[ga["x0_bin_label"] == bin_label].copy()
                n_total = len(gb)
                counts = gb["obs_class"].value_counts(dropna=False).to_dict()

                for c in CLASS_ORDER:
                    n = int(counts.get(c, 0))
                    prop = n / n_total if n_total > 0 else np.nan
                    rows.append({
                        "alpha_label": a_lab,
                        "alpha_float": a_val,
                        "design": design,
                        "x0_bin_label": bin_label,
                        "x0_bin_order": all_bin_labels.index(bin_label),
                        "obs_class": c,
                        "n": n,
                        "n_total": n_total,
                        "prop": prop,
                    })

    return pd.DataFrame(rows).sort_values(
        ["design", "alpha_float", "x0_bin_order", "obs_class"]
    ).reset_index(drop=True)


# =========================================================
# Plot builders
# =========================================================

def make_single_ax_figure(figsize: Tuple[float, float]) -> tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


def plot_class_counts_vs_alpha(summary: pd.DataFrame, design: str, out_dir: Path, args: argparse.Namespace):
    g = summary[summary["design"] == design].copy()
    if g.empty:
        return

    x = np.array(sorted(g["alpha_float"].dropna().unique().tolist()), dtype=float)
    fig, ax = make_single_ax_figure((args.fig_w, args.fig_h))

    for c in CLASS_ORDER:
        gc = g[g["obs_class"] == c].sort_values("alpha_float")
        ax.plot(
            gc["alpha_float"],
            gc["n"],
            marker="o",
            linewidth=args.line_w,
            markersize=args.marker_s,
            alpha=args.marker_alpha,
            label=c,
            color=DEFAULT_CLASS_COLORS[c],
        )

    ax.set_xlabel("alpha")
    ax.set_ylabel("n transitions")
    if not args.no_titles:
        ax.set_title(f"{design} | transition counts")
    apply_alpha_axis(ax, x)
    maybe_apply_limits(ax, parse_lim(args.xlim_alpha), parse_lim(args.ylim_counts))
    apply_legend(ax, args)
    apply_layout(fig, args)
    save_figure(fig, out_dir / f"obsclass_counts_vs_alpha_{design}.png", args)


def plot_class_props_vs_alpha(summary: pd.DataFrame, design: str, out_dir: Path, args: argparse.Namespace):
    g = summary[summary["design"] == design].copy()
    if g.empty:
        return

    x = np.array(sorted(g["alpha_float"].dropna().unique().tolist()), dtype=float)
    fig, ax = make_single_ax_figure((args.fig_w, args.fig_h))

    for c in CLASS_ORDER:
        gc = g[g["obs_class"] == c].sort_values("alpha_float")
        ax.plot(
            gc["alpha_float"],
            gc["prop"],
            marker="o",
            linewidth=args.line_w,
            markersize=args.marker_s,
            alpha=args.marker_alpha,
            label=c,
            color=DEFAULT_CLASS_COLORS[c],
        )

    ax.set_xlabel("alpha")
    ax.set_ylabel("proportion")
    if not args.no_titles:
        ax.set_title(f"{design} | transition proportions")
    apply_alpha_axis(ax, x)
    maybe_apply_limits(ax, parse_lim(args.xlim_alpha), parse_lim(args.ylim_props))
    apply_legend(ax, args)
    apply_layout(fig, args)
    save_figure(fig, out_dir / f"obsclass_props_vs_alpha_{design}.png", args)


def plot_class_props_stacked(summary: pd.DataFrame, design: str, out_dir: Path, args: argparse.Namespace):
    g = summary[summary["design"] == design].copy()
    if g.empty:
        return

    alphas = (
        g[["alpha_label", "alpha_float"]]
        .drop_duplicates()
        .sort_values(["alpha_float", "alpha_label"])
        .reset_index(drop=True)
    )

    x = np.arange(len(alphas))
    labels = [format_alpha_label(v) for v in alphas["alpha_float"].to_numpy(dtype=float)]
    fig, ax = make_single_ax_figure((args.fig_w, args.fig_h))
    bottom = np.zeros(len(alphas), dtype=float)

    for c in CLASS_ORDER:
        vals = []
        for _, ar in alphas.iterrows():
            vv = g[(g["alpha_label"] == ar["alpha_label"]) & (g["obs_class"] == c)]["prop"]
            vals.append(float(vv.iloc[0]) if len(vv) else 0.0)
        vals = np.array(vals, dtype=float)
        ax.bar(x, vals, bottom=bottom, label=c, alpha=args.bar_alpha, color=DEFAULT_CLASS_COLORS[c])
        bottom += vals

    ax.set_xlabel("alpha")
    ax.set_ylabel("proportion")
    if not args.no_titles:
        ax.set_title(f"{design} | stacked transition composition")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    maybe_apply_limits(ax, None, parse_lim(args.ylim_props))
    apply_legend(ax, args)
    apply_layout(fig, args)
    save_figure(fig, out_dir / f"obsclass_props_stacked_{design}.png", args)


def plot_tt_summary(summary: pd.DataFrame, design: str, out_dir: Path, args: argparse.Namespace):
    g = summary[(summary["design"] == design) & (summary["obs_class"] == "TT")].copy()
    if g.empty:
        return

    g = g.sort_values("alpha_float").reset_index(drop=True)
    x = g["alpha_float"].to_numpy(dtype=float)

    fig, ax = make_single_ax_figure((args.fig_w, args.fig_h))
    ax.plot(x, g["n"], marker="o", linewidth=args.line_w, markersize=args.marker_s, alpha=args.marker_alpha, color=DEFAULT_CLASS_COLORS["TT"])
    ax.set_xlabel("alpha")
    ax.set_ylabel("TT count")
    if not args.no_titles:
        ax.set_title(f"{design} | TT absolute abundance")
    apply_alpha_axis(ax, x)
    maybe_apply_limits(ax, parse_lim(args.xlim_alpha), parse_lim(args.ylim_tt_count))
    apply_layout(fig, args)
    save_figure(fig, out_dir / f"TT_count_vs_alpha_{design}.png", args)

    fig, ax = make_single_ax_figure((args.fig_w, args.fig_h))
    ax.plot(x, g["prop"], marker="o", linewidth=args.line_w, markersize=args.marker_s, alpha=args.marker_alpha, color=DEFAULT_CLASS_COLORS["TT"])
    ax.set_xlabel("alpha")
    ax.set_ylabel("TT proportion")
    if not args.no_titles:
        ax.set_title(f"{design} | TT relative abundance")
    apply_alpha_axis(ax, x)
    maybe_apply_limits(ax, parse_lim(args.xlim_alpha), parse_lim(args.ylim_tt_prop))
    apply_layout(fig, args)
    save_figure(fig, out_dir / f"TT_prop_vs_alpha_{design}.png", args)

    ref = float(g["n"].iloc[-1]) if len(g) else np.nan
    retention = g["n"] / ref if np.isfinite(ref) and ref > 0 else np.nan

    fig, ax = make_single_ax_figure((args.fig_w, args.fig_h))
    ax.plot(x, retention, marker="o", linewidth=args.line_w, markersize=args.marker_s, alpha=args.marker_alpha, color=DEFAULT_CLASS_COLORS["TT"])
    ax.set_xlabel("alpha")
    ax.set_ylabel("TT retention")
    if not args.no_titles:
        ax.set_title(f"{design} | TT retention (vs max alpha)")
    apply_alpha_axis(ax, x)
    maybe_apply_limits(ax, parse_lim(args.xlim_alpha), parse_lim(args.ylim_tt_retention))
    apply_layout(fig, args)
    save_figure(fig, out_dir / f"TT_retention_vs_alpha_{design}.png", args)


def plot_tt_props_by_x0(summary_x0: pd.DataFrame, design: str, out_dir: Path, bin_labels: List[str], args: argparse.Namespace):
    g = summary_x0[(summary_x0["design"] == design) & (summary_x0["obs_class"] == "TT")].copy()
    if g.empty:
        return

    alphas = (
        g[["alpha_label", "alpha_float"]]
        .drop_duplicates()
        .sort_values(["alpha_float", "alpha_label"])
        .reset_index(drop=True)
    )
    x = np.array(sorted(g["alpha_float"].dropna().unique().tolist()), dtype=float)
    fig, ax = make_single_ax_figure((args.fig_w_wide, args.fig_h_wide))

    for i, bin_label in enumerate(bin_labels):
        gb = g[g["x0_bin_label"] == bin_label].sort_values("alpha_float")
        y = []
        for _, ar in alphas.iterrows():
            vv = gb[gb["alpha_label"] == ar["alpha_label"]]["prop"]
            y.append(float(vv.iloc[0]) if len(vv) else np.nan)
        ax.plot(
            x,
            np.array(y, dtype=float),
            marker="o",
            linewidth=args.line_w,
            markersize=args.marker_s,
            alpha=args.marker_alpha,
            label=bin_label,
            color=f"C{i % 10}",
        )

    ax.set_xlabel("alpha")
    ax.set_ylabel("TT proportion")
    if not args.no_titles:
        ax.set_title(f"{design} | TT proportion by x0 bin")
    apply_alpha_axis(ax, x)
    maybe_apply_limits(ax, parse_lim(args.xlim_alpha), parse_lim(args.ylim_tt_prop_x0))
    apply_legend(ax, args)
    apply_layout(fig, args)
    save_figure(fig, out_dir / f"TT_prop_by_x0_vs_alpha_{design}.png", args)


def plot_tt_counts_across_bins(summary_x0: pd.DataFrame, design: str, out_dir: Path, bin_labels: List[str], args: argparse.Namespace):
    g = summary_x0[(summary_x0["design"] == design) & (summary_x0["obs_class"] == "TT")].copy()
    if g.empty:
        return

    alphas = (
        g[["alpha_label", "alpha_float"]]
        .drop_duplicates()
        .sort_values(["alpha_float", "alpha_label"])
        .reset_index(drop=True)
    )

    x = np.arange(len(bin_labels))
    series = []
    for _, ar in alphas.iterrows():
        y = []
        for bin_label in bin_labels:
            vv = g[(g["alpha_label"] == ar["alpha_label"]) & (g["x0_bin_label"] == bin_label)]["n"]
            y.append(float(vv.iloc[0]) if len(vv) else 0.0)
        series.append((ar["alpha_label"], np.array(y, dtype=float)))

    fig, ax = make_single_ax_figure((args.fig_w_wide, args.fig_h_wide))
    for i, (alpha_label, y) in enumerate(series):
        ax.plot(
            x,
            y,
            marker="o",
            linewidth=args.line_w,
            markersize=args.marker_s,
            alpha=args.marker_alpha,
            label=alpha_label,
            color=f"C{i % 10}",
        )

    ax.set_xlabel("x0 bin (log10 frequency at t0)")
    ax.set_ylabel("TT absolute count")
    if not args.no_titles:
        ax.set_title(f"{design} | TT absolute abundance across x0 bins")
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=45, ha="right")
    maybe_apply_limits(ax, None, parse_lim(args.ylim_tt_count_x0))
    apply_legend(ax, args, title="alpha")
    apply_layout(fig, args)
    save_figure(fig, out_dir / f"TT_count_by_x0_vs_alpha_{design}.png", args)

    if any(np.any(y > 0) for _, y in series):
        fig, ax = make_single_ax_figure((args.fig_w_wide, args.fig_h_wide))
        for i, (alpha_label, y) in enumerate(series):
            y_plot = np.where(y > 0, y, np.nan)
            ax.plot(
                x,
                y_plot,
                marker="o",
                linewidth=args.line_w,
                markersize=args.marker_s,
                alpha=args.marker_alpha,
                label=alpha_label,
                color=f"C{i % 10}",
            )

        ax.set_yscale("log")
        ax.set_xlabel("x0 bin (log10 frequency at t0)")
        ax.set_ylabel("TT absolute count (log scale)")
        if not args.no_titles:
            ax.set_title(f"{design} | TT absolute abundance across x0 bins (log y)")
        ax.set_xticks(x)
        ax.set_xticklabels(bin_labels, rotation=45, ha="right")
        maybe_apply_limits(ax, None, parse_lim(args.ylim_tt_count_x0_log))
        apply_legend(ax, args, title="alpha")
        apply_layout(fig, args)
        save_figure(fig, out_dir / f"TT_count_by_x0_vs_alpha_{design}_logy.png", args)


# =========================================================
# Diagnostics
# =========================================================

def print_consistency_checks(df: pd.DataFrame, summary_x0: Optional[pd.DataFrame]) -> None:
    print("\n[CHECK] Consistency between total TT and binned TT")
    designs = sorted(df["design"].dropna().unique().tolist())

    for design in designs:
        gd = df[df["design"] == design].copy()
        alphas = sorted(gd["alpha_label"].dropna().unique().tolist())
        for alpha_label in alphas:
            gda = gd[gd["alpha_label"] == alpha_label].copy()
            tt_total = int((gda["obs_class"] == "TT").sum())
            tt_with_x0 = int(((gda["obs_class"] == "TT") & gda["x0_log10"].notna()).sum())

            if summary_x0 is not None and not summary_x0.empty:
                tt_binned = summary_x0[
                    (summary_x0["design"] == design) &
                    (summary_x0["alpha_label"] == alpha_label) &
                    (summary_x0["obs_class"] == "TT")
                ]["n"].sum()
            else:
                tt_binned = np.nan

            print(
                f"  design={design:15s} alpha={alpha_label:6s} "
                f"TT_total={tt_total:8d} TT_with_x0={tt_with_x0:8d} "
                f"TT_binned={int(tt_binned) if pd.notna(tt_binned) else 'NA':>8}"
            )


# =========================================================
# CLI
# =========================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read transitions across alpha thresholds, count transition classes, stratify by x0, and generate publication-ready plots."
    )

    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--inputs", nargs="+", type=Path, help="List of transition CSV files.")
    src.add_argument("--glob", type=str, help="Glob pattern for transition CSV files.")

    parser.add_argument("--out-data-dir", type=Path, required=True, help="Directory where CSV tables will be saved.")
    parser.add_argument("--out-fig-dir", type=Path, required=True, help="Directory where figures will be saved.")

    parser.add_argument("--alpha-label", type=str, default=None, help="Explicit alpha label to assign to all input files, e.g. p_02.")
    parser.add_argument("--alpha-value", type=float, default=None, help="Explicit alpha numeric value to assign to all input files, e.g. 0.02.")

    parser.add_argument("--font", default="Arial")
    parser.add_argument("--fontsize", type=int, default=10)

    parser.add_argument("--fig-w", type=float, default=5.8, help="Default figure width in inches.")
    parser.add_argument("--fig-h", type=float, default=3.8, help="Default figure height in inches.")
    parser.add_argument("--fig-w-wide", type=float, default=7.0, help="Wide figure width in inches.")
    parser.add_argument("--fig-h-wide", type=float, default=4.2, help="Wide figure height in inches.")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--pdf", action="store_true")

    parser.add_argument("--line-w", type=float, default=1.8)
    parser.add_argument("--marker-s", type=float, default=4.5, help="Marker size for line plots.")
    parser.add_argument("--marker-alpha", type=float, default=0.95)
    parser.add_argument("--bar-alpha", type=float, default=0.95)
    parser.add_argument("--axes-w", type=float, default=1.0)
    parser.add_argument("--tick-w", type=float, default=1.0)
    parser.add_argument("--tick-len", type=float, default=4.0)

    parser.add_argument("--legend-mode", choices=["inside", "outside", "none"], default="outside")
    parser.add_argument("--legend-loc", default="upper left")
    parser.add_argument("--legend-frame", type=int, default=0)
    parser.add_argument("--legend-fontsize", type=float, default=9.0)
    parser.add_argument("--legend-bbox-x", type=float, default=1.02)
    parser.add_argument("--legend-bbox-y", type=float, default=1.0)
    parser.add_argument("--legend-pad-right", type=float, default=0.78)
    parser.add_argument("--legend-pad-inches", type=float, default=0.05)

    parser.add_argument("--no-titles", action="store_true")

    parser.add_argument("--xlim-alpha", default="")
    parser.add_argument("--ylim-counts", default="")
    parser.add_argument("--ylim-props", default="")
    parser.add_argument("--ylim-tt-count", default="")
    parser.add_argument("--ylim-tt-prop", default="")
    parser.add_argument("--ylim-tt-retention", default="")
    parser.add_argument("--ylim-tt-prop-x0", default="")
    parser.add_argument("--ylim-tt-count-x0", default="")
    parser.add_argument("--ylim-tt-count-x0-log", default="")

    parser.add_argument(
        "--x0-bin-edges",
        nargs="+",
        type=float,
        default=[-6, -5, -4, -3, -2, -1],
        help="Internal fixed bin edges on log10 frequency at t0. Underflow and overflow bins are added if tail bins are enabled.",
    )
    parser.add_argument(
        "--include-tail-bins",
        action="store_true",
        default=True,
        help="Include underflow and overflow bins so that the sum across bins matches all transitions with valid x0.",
    )
    parser.add_argument(
        "--x0-mode",
        choices=["auto", "raw", "ln", "log10"],
        default="auto",
        help="How to interpret x0 when freq_t0/f0 are absent.",
    )

    return parser.parse_args()


# =========================================================
# Main
# =========================================================

def main() -> None:
    args = parse_args()

    if (args.alpha_label is None) ^ (args.alpha_value is None):
        raise ValueError("You must provide both --alpha-label and --alpha-value together, or neither.")

    data_dir = ensure_dir(args.out_data_dir)
    fig_dir = ensure_dir(args.out_fig_dir)

    set_mpl_defaults(args)

    if args.glob:
        paths = sorted(Path(".").glob(args.glob))
        if not paths:
            raise FileNotFoundError(f"No files matched glob: {args.glob}")
    else:
        paths = [Path(p) for p in args.inputs]
        for fp in paths:
            if not fp.exists():
                raise FileNotFoundError(str(fp))

    edges = list(args.x0_bin_edges)
    if len(edges) < 2:
        raise ValueError("Need at least two x0 bin edges.")
    if sorted(edges) != edges:
        raise ValueError("x0 bin edges must be increasing.")

    df = load_inputs(
        paths,
        x0_mode=args.x0_mode,
        cli_alpha_label=args.alpha_label,
        cli_alpha_value=args.alpha_value,
    )

    summary = build_global_summary(df)
    summary_path = data_dir / "transition_class_summary_by_alpha.csv"
    summary.to_csv(summary_path, index=False)
    print(f"[SAVED DATA] {summary_path}")

    totals = (
        summary[["alpha_label", "alpha_float", "design", "n_total"]]
        .drop_duplicates()
        .sort_values(["design", "alpha_float"])
        .reset_index(drop=True)
    )
    totals_path = data_dir / "transition_total_by_alpha.csv"
    totals.to_csv(totals_path, index=False)
    print(f"[SAVED DATA] {totals_path}")

    all_bin_labels = interval_labels_with_tails(edges) if args.include_tail_bins else interval_labels(edges)
    bins_df = pd.DataFrame({
        "x0_bin_order": np.arange(len(all_bin_labels), dtype=int),
        "x0_bin_label": all_bin_labels,
    })
    bins_path = data_dir / "x0_fixed_bins.csv"
    bins_df.to_csv(bins_path, index=False)
    print(f"[SAVED DATA] {bins_path}")

    designs = sorted(summary["design"].dropna().unique().tolist())
    for design in designs:
        plot_class_counts_vs_alpha(summary, design, fig_dir, args)
        plot_class_props_vs_alpha(summary, design, fig_dir, args)
        plot_class_props_stacked(summary, design, fig_dir, args)
        plot_tt_summary(summary, design, fig_dir, args)

    summary_x0 = None
    if "x0_log10" in df.columns:
        df_x0 = assign_x0_bins(df, edges, include_tails=args.include_tail_bins)
        if not df_x0.empty:
            summary_x0 = build_x0_summary(df_x0, all_bin_labels)
            summary_x0_path = data_dir / "transition_class_summary_by_alpha_x0bin.csv"
            summary_x0.to_csv(summary_x0_path, index=False)
            print(f"[SAVED DATA] {summary_x0_path}")

            for design in designs:
                plot_tt_props_by_x0(summary_x0, design, fig_dir, all_bin_labels, args)
                plot_tt_counts_across_bins(summary_x0, design, fig_dir, all_bin_labels, args)

    print(f"[DONE] Data saved to: {data_dir}")
    print(f"[DONE] Figures saved to: {fig_dir}")
    print(f"[INFO] TOTAL transitions: {len(df)}")
    print(f"[INFO] TOTAL TT: {len(df[df['obs_class'] == 'TT'])}")
    print(f"[INFO] Rows with valid x0_log10: {df['x0_log10'].notna().sum()}")
    print(f"[INFO] TT with valid x0_log10: {df.loc[df['obs_class'] == 'TT', 'x0_log10'].notna().sum()}")

    print_consistency_checks(df, summary_x0)


if __name__ == "__main__":
    main()
