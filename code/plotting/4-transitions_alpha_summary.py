#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transitions_alpha_summary.py

Unified script to:
1) read transition CSVs across multiple alpha thresholds (p_* folders or filenames)
2) count transition classes (TT, TF, FT, FF)
3) stratify transitions by bins of the FIRST time-point frequency of the transition
4) generate summary plots

Key features
------------
- x0 binning is explicitly based on the first time-point frequency whenever
  available (freq_t0 or f0), converted to log10.
- If only x0 is present, the script tries to infer whether x0 is raw frequency,
  ln(frequency), or log10(frequency), and converts it to log10 consistently.
- Underflow / overflow bins are included by default, so the sum across bins
  matches the total number of transitions with valid x0 information.
- The script prints per-design consistency checks between total TT and binned TT.
- DATA and FIGURES are saved in separate output directories.
- alpha can be inferred in one of three ways:
    1) from parent folder name, e.g. p_005, p_01, p_02, p_05
    2) from filename, if it contains a token like p_005 / p_01 / p_02 / p_05
    3) from explicit CLI arguments (--alpha-label and --alpha-value) when
       all input files belong to the same alpha level

Expected input
--------------
- one or more transition CSV files
- alpha is identified using the following priority:
    1) explicit CLI arguments (--alpha-label, --alpha-value), if provided
    2) parent folder name, e.g. p_005, p_01, p_02, p_05
    3) filename token containing p_005 / p_01 / p_02 / p_05
- design inferred from filename: adj / lag / all

Examples of valid file layouts
------------------------------
A) alpha in parent folder:
    results/p_02/transitions_adj.csv

B) alpha in filename:
    results/4-transitions_global/transitions_adj_p_02.csv

C) alpha provided explicitly at launch:
    python transitions_alpha_summary.py \
      --inputs results/4-transitions_global/transitions_adj.csv \
      --alpha-label p_02 \
      --alpha-value 0.02 \
      --out-data-dir ./tables \
      --out-fig-dir ./figures

Required column:
    obs_class

Optional columns:
    x0
    freq_t0
    f0
    dt
    subject

Main outputs
------------
data dir:
    transition_class_summary_by_alpha.csv
    transition_total_by_alpha.csv
    transition_class_summary_by_alpha_x0bin.csv
    x0_fixed_bins.csv

figure dir:
    obsclass_counts_vs_alpha_<design>.png
    obsclass_props_vs_alpha_<design>.png
    obsclass_props_stacked_<design>.png
    TT_count_vs_alpha_<design>.png
    TT_prop_vs_alpha_<design>.png
    TT_retention_vs_alpha_<design>.png
    TT_prop_by_x0_vs_alpha_<design>.png
    TT_count_by_x0_vs_alpha_<design>.png
    TT_count_by_x0_vs_alpha_<design>_logy.png
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


# =========================================================
# Helpers
# =========================================================

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def set_mpl_defaults(font_family: str = "Arial", font_size: int = 10) -> None:
    plt.rcParams.update({
        "font.family": font_family,
        "font.size": font_size,
        "axes.titlesize": font_size,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": font_size,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
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
    """
    Infer alpha with priority:
    1) explicit CLI alpha
    2) parent folder name like p_02
    3) filename containing token like p_02
    """
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


def apply_alpha_axis(ax, x: np.ndarray):
    ax.set_xscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([format_alpha_label(v) for v in x], rotation=45, ha="right")
    ax.xaxis.set_minor_locator(NullLocator())


def savefig(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVED FIGURE] {path}")


# =========================================================
# x0 handling
# =========================================================

def _safe_log10(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    s = s.where(s > 0)
    return np.log10(s)


def _infer_x0_scale_from_series(x: pd.Series) -> str:
    """
    Infer whether x0 is raw frequency, ln(f), or log10(f).
    Heuristic, grounded in typical repertoire frequency ranges.
    """
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
    """
    Return x0 on a common log10 scale and a source description.
    Priority:
    1) freq_t0 or f0 (explicit first-time frequency)
    2) x0 (auto-inferred or user-declared scale)
    """
    if "freq_t0" in df.columns:
        return _safe_log10(df["freq_t0"]), "freq_t0(raw)->log10"

    if "f0" in df.columns:
        return _safe_log10(df["f0"]), "f0(raw)->log10"

    if "x0" not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float), "missing"

    x = pd.to_numeric(df["x0"], errors="coerce")

    if x0_mode == "auto":
        mode = _infer_x0_scale_from_series(x)
    else:
        mode = x0_mode

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
# Deterministic counts
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

    out = pd.DataFrame(rows).sort_values(
        ["design", "alpha_float", "obs_class"]
    ).reset_index(drop=True)

    return out


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

    out = pd.DataFrame(rows).sort_values(
        ["design", "alpha_float", "x0_bin_order", "obs_class"]
    ).reset_index(drop=True)

    return out


# =========================================================
# Plots
# =========================================================

def plot_class_counts_vs_alpha(summary: pd.DataFrame, design: str, out_dir: Path):
    g = summary[summary["design"] == design].copy()
    if g.empty:
        return

    x = np.array(sorted(g["alpha_float"].dropna().unique().tolist()), dtype=float)

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    for c in CLASS_ORDER:
        gc = g[g["obs_class"] == c].sort_values("alpha_float")
        ax.plot(gc["alpha_float"], gc["n"], marker="o", linewidth=2, label=c)

    ax.set_xlabel("alpha")
    ax.set_ylabel("n transitions")
    ax.set_title(f"{design} | transition counts")
    apply_alpha_axis(ax, x)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, out_dir / f"obsclass_counts_vs_alpha_{design}.png")


def plot_class_props_vs_alpha(summary: pd.DataFrame, design: str, out_dir: Path):
    g = summary[summary["design"] == design].copy()
    if g.empty:
        return

    x = np.array(sorted(g["alpha_float"].dropna().unique().tolist()), dtype=float)

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    for c in CLASS_ORDER:
        gc = g[g["obs_class"] == c].sort_values("alpha_float")
        ax.plot(gc["alpha_float"], gc["prop"], marker="o", linewidth=2, label=c)

    ax.set_xlabel("alpha")
    ax.set_ylabel("proportion")
    ax.set_title(f"{design} | transition proportions")
    apply_alpha_axis(ax, x)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, out_dir / f"obsclass_props_vs_alpha_{design}.png")


def plot_class_props_stacked(summary: pd.DataFrame, design: str, out_dir: Path):
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

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    bottom = np.zeros(len(alphas), dtype=float)

    for c in CLASS_ORDER:
        vals = []
        for _, ar in alphas.iterrows():
            vv = g[(g["alpha_label"] == ar["alpha_label"]) & (g["obs_class"] == c)]["prop"]
            vals.append(float(vv.iloc[0]) if len(vv) else 0.0)
        vals = np.array(vals, dtype=float)
        ax.bar(x, vals, bottom=bottom, label=c)
        bottom += vals

    ax.set_xlabel("alpha")
    ax.set_ylabel("proportion")
    ax.set_title(f"{design} | stacked transition composition")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, out_dir / f"obsclass_props_stacked_{design}.png")


def plot_tt_summary(summary: pd.DataFrame, design: str, out_dir: Path):
    g = summary[(summary["design"] == design) & (summary["obs_class"] == "TT")].copy()
    if g.empty:
        return

    g = g.sort_values("alpha_float").reset_index(drop=True)
    x = g["alpha_float"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(5.2, 3.5))
    ax.plot(x, g["n"], marker="o", linewidth=2)
    ax.set_xlabel("alpha")
    ax.set_ylabel("TT count")
    ax.set_title(f"{design} | TT absolute abundance")
    apply_alpha_axis(ax, x)
    savefig(fig, out_dir / f"TT_count_vs_alpha_{design}.png")

    fig, ax = plt.subplots(figsize=(5.2, 3.5))
    ax.plot(x, g["prop"], marker="o", linewidth=2)
    ax.set_xlabel("alpha")
    ax.set_ylabel("TT proportion")
    ax.set_title(f"{design} | TT relative abundance")
    apply_alpha_axis(ax, x)
    savefig(fig, out_dir / f"TT_prop_vs_alpha_{design}.png")

    ref = float(g["n"].iloc[-1]) if len(g) else np.nan
    retention = g["n"] / ref if np.isfinite(ref) and ref > 0 else np.nan

    fig, ax = plt.subplots(figsize=(5.2, 3.5))
    ax.plot(x, retention, marker="o", linewidth=2)
    ax.set_xlabel("alpha")
    ax.set_ylabel("TT retention")
    ax.set_title(f"{design} | TT retention (vs max alpha)")
    apply_alpha_axis(ax, x)
    savefig(fig, out_dir / f"TT_retention_vs_alpha_{design}.png")


def plot_tt_props_by_x0(summary_x0: pd.DataFrame, design: str, out_dir: Path, bin_labels: List[str]):
    g = summary_x0[(summary_x0["design"] == design) & (summary_x0["obs_class"] == "TT")].copy()
    if g.empty:
        return

    alphas = (
        g[["alpha_label", "alpha_float"]]
        .drop_duplicates()
        .sort_values(["alpha_float", "alpha_label"])
    )

    x = np.array(sorted(g["alpha_float"].dropna().unique().tolist()), dtype=float)

    fig, ax = plt.subplots(figsize=(6.6, 4.0))

    for bin_label in bin_labels:
        gb = g[g["x0_bin_label"] == bin_label].sort_values("alpha_float")
        y = []
        for _, ar in alphas.iterrows():
            vv = gb[gb["alpha_label"] == ar["alpha_label"]]["prop"]
            y.append(float(vv.iloc[0]) if len(vv) else np.nan)

        ax.plot(x, np.array(y, dtype=float), marker="o", linewidth=2, label=bin_label)

    ax.set_xlabel("alpha")
    ax.set_ylabel("TT proportion")
    ax.set_title(f"{design} | TT proportion by x0 bin")
    apply_alpha_axis(ax, x)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, out_dir / f"TT_prop_by_x0_vs_alpha_{design}.png")


def plot_tt_counts_across_bins(summary_x0: pd.DataFrame, design: str, out_dir: Path, bin_labels: List[str]):
    g = summary_x0[
        (summary_x0["design"] == design) &
        (summary_x0["obs_class"] == "TT")
    ].copy()
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
            vv = g[
                (g["alpha_label"] == ar["alpha_label"]) &
                (g["x0_bin_label"] == bin_label)
            ]["n"]
            y.append(float(vv.iloc[0]) if len(vv) else 0.0)
        series.append((ar["alpha_label"], np.array(y, dtype=float)))

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for alpha_label, y in series:
        ax.plot(x, y, marker="o", linewidth=2, label=alpha_label)

    ax.set_xlabel("x0 bin (log10 frequency at t0)")
    ax.set_ylabel("TT absolute count")
    ax.set_title(f"{design} | TT absolute abundance across x0 bins")
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=45, ha="right")
    ax.legend(title="alpha", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, out_dir / f"TT_count_by_x0_vs_alpha_{design}.png")

    if any(np.any(y > 0) for _, y in series):
        fig, ax = plt.subplots(figsize=(7.0, 4.2))
        for alpha_label, y in series:
            y_plot = np.where(y > 0, y, np.nan)
            ax.plot(x, y_plot, marker="o", linewidth=2, label=alpha_label)

        ax.set_yscale("log")
        ax.set_xlabel("x0 bin (log10 frequency at t0)")
        ax.set_ylabel("TT absolute count (log scale)")
        ax.set_title(f"{design} | TT absolute abundance across x0 bins (log y)")
        ax.set_xticks(x)
        ax.set_xticklabels(bin_labels, rotation=45, ha="right")
        ax.legend(title="alpha", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        savefig(fig, out_dir / f"TT_count_by_x0_vs_alpha_{design}_logy.png")


# =========================================================
# Diagnostics
# =========================================================

def print_consistency_checks(df: pd.DataFrame, summary_x0: Optional[pd.DataFrame]):
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
# Main
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="Read transitions across alpha thresholds, count transition classes correctly, and plot summaries."
    )

    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--inputs", nargs="+", type=Path, help="List of transition CSV files.")
    src.add_argument("--glob", type=str, help="Glob pattern for transition CSV files.")

    parser.add_argument(
        "--out-data-dir",
        type=Path,
        required=True,
        help="Directory where CSV tables will be saved."
    )
    parser.add_argument(
        "--out-fig-dir",
        type=Path,
        required=True,
        help="Directory where figures will be saved."
    )

    parser.add_argument(
        "--alpha-label",
        type=str,
        default=None,
        help="Explicit alpha label to assign to all input files, e.g. p_02."
    )
    parser.add_argument(
        "--alpha-value",
        type=float,
        default=None,
        help="Explicit alpha numeric value to assign to all input files, e.g. 0.02."
    )

    parser.add_argument("--font", default="Arial")
    parser.add_argument("--fontsize", type=int, default=10)

    parser.add_argument(
        "--x0-bin-edges",
        nargs="+",
        type=float,
        default=[-6, -5, -4, -3, -2, -1],
        help=(
            "Internal fixed bin edges on log10 frequency at t0. "
            "With default settings, underflow and overflow bins are also added."
        )
    )
    parser.add_argument(
        "--include-tail-bins",
        action="store_true",
        default=True,
        help=(
            "Include underflow and overflow bins so that the sum across bins "
            "matches all transitions with valid x0."
        )
    )
    parser.add_argument(
        "--x0-mode",
        choices=["auto", "raw", "ln", "log10"],
        default="auto",
        help="How to interpret x0 when freq_t0/f0 are absent."
    )

    args = parser.parse_args()

    if (args.alpha_label is None) ^ (args.alpha_value is None):
        raise ValueError(
            "You must provide both --alpha-label and --alpha-value together, or neither of them."
        )

    data_dir = ensure_dir(args.out_data_dir)
    fig_dir = ensure_dir(args.out_fig_dir)

    set_mpl_defaults(args.font, args.fontsize)

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

    if args.include_tail_bins:
        all_bin_labels = interval_labels_with_tails(edges)
    else:
        all_bin_labels = interval_labels(edges)

    bins_df = pd.DataFrame({
        "x0_bin_order": np.arange(len(all_bin_labels), dtype=int),
        "x0_bin_label": all_bin_labels,
    })
    bins_path = data_dir / "x0_fixed_bins.csv"
    bins_df.to_csv(bins_path, index=False)
    print(f"[SAVED DATA] {bins_path}")

    designs = sorted(summary["design"].dropna().unique().tolist())

    for design in designs:
        plot_class_counts_vs_alpha(summary, design, fig_dir)
        plot_class_props_vs_alpha(summary, design, fig_dir)
        plot_class_props_stacked(summary, design, fig_dir)
        plot_tt_summary(summary, design, fig_dir)

    summary_x0 = None
    if "x0_log10" in df.columns:
        df_x0 = assign_x0_bins(df, edges, include_tails=args.include_tail_bins)
        if not df_x0.empty:
            summary_x0 = build_x0_summary(df_x0, all_bin_labels)
            summary_x0_path = data_dir / "transition_class_summary_by_alpha_x0bin.csv"
            summary_x0.to_csv(summary_x0_path, index=False)
            print(f"[SAVED DATA] {summary_x0_path}")

            for design in designs:
                plot_tt_props_by_x0(summary_x0, design, fig_dir, all_bin_labels)
                plot_tt_counts_across_bins(summary_x0, design, fig_dir, all_bin_labels)

    print(f"[DONE] Data saved to: {data_dir}")
    print(f"[DONE] Figures saved to: {fig_dir}")
    print(f"[INFO] TOTAL transitions: {len(df)}")
    print(f"[INFO] TOTAL TT: {len(df[df['obs_class'] == 'TT'])}")
    print(f"[INFO] Rows with valid x0_log10: {df['x0_log10'].notna().sum()}")
    print(f"[INFO] TT with valid x0_log10: {df.loc[df['obs_class'] == 'TT', 'x0_log10'].notna().sum()}")

    print_consistency_checks(df, summary_x0)


if __name__ == "__main__":
    main()