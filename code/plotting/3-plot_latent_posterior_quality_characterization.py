#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3-plot_latent_posterior_quality_characterization.py
===================================================

Plot-only companion to ClonoDynamics Step 3:

    3-latent_posterior_quality_characterization.py

Role in the pipeline
--------------------
Step 3 characterizes the statistical quality of the latent clonotype states
inferred by Step 2. This plotter visualizes those Step-3 outputs.

The script is deliberately plot-only:

    - it does not refit the Step-2 latent model;
    - it does not redefine detectability;
    - it does not score or filter clonotypes;
    - it does not modify the Step-3 QC table;
    - it does not generate inputs for Step 4.

For large clonotype-level tables, a deterministic systematic sample is used
only for rendering histograms, hexbin panels, and smooth curves. Numerical
annotations are read from full-dataset Step-3 summary tables whenever
available.

Analysis inputs
---------------
Required --analysis-dir: output directory produced by

    3-latent_posterior_quality_characterization.py

Required files:

    clonotype_latent_posterior_qc.parquet
        Clonotype-timepoint posterior-quality table.

    clonotype_latent_posterior_qc_metric_summary.csv
        Full-dataset univariate QC summaries used for annotations.

    clonotype_latent_posterior_qc_key_relationships.csv
        Full-dataset selected Spearman relationships.

    qc_metric_spearman_correlations.csv
        Full-dataset QC correlation matrix.

    clonotype_latent_posterior_qc_by_abundance.csv
        Abundance-binned posterior-quality summaries.

    clonotype_latent_posterior_qc_by_replicate_difference.csv
        Replicate-difference-binned posterior-quality summaries.

    clonotype_latent_posterior_qc_abundance_reference.csv
        Abundance reference values used by the abundance panels.

Generated figures
-----------------
Main-manuscript panels (Fig. 2D-F):

    fig2D_abundance_vs_posterior_width
    fig2E_abundance_vs_replicate_overlap
    fig2F_state_vs_union_detectability
    fig2_DEF_latent_state_inference

Supplementary Fig. S1 panels:

    suppS1A_posterior_width_distribution
    suppS1B_replicate_overlap_distribution
    suppS1C_detectability_correction_distribution
    suppS1D_detectability_correction_vs_posterior_width
    suppS1E_detectability_correction_vs_state_detectability
    suppS1F_qc_metric_correlation

The former entropy/support and replicate-detection diagnostic plots are no
longer rendered by default because they are redundant with the selected
main/supplementary evidence. Their numerical summaries remain in Step-3 outputs.

Output directories
------------------
By default:

    <analysis-dir>/figures_latent_posterior_quality/main/
    <analysis-dir>/figures_latent_posterior_quality/supplementary/

A different root can be supplied with:

    --figures-dir

Each generated panel is exported as PDF and PNG.

Figure customization
--------------------
Two dictionaries near the top of this script are the primary publication
layout controls:

    FIGURE_SIZES_INCHES
        Individual (width, height) values in inches.

    AXIS_CONFIG
        Optional per-panel axis controls:
            xlim
            ylim
            xscale
            yscale
            xticks
            yticks

Edit these dictionaries to modify one figure without changing the analytical
logic or data selection.

Example
-------
python3 ./code/plotting/3-plot_latent_posterior_quality_characterization.py \
    --analysis-dir \
    ./dataset_longitudinal_clonodynamics_results/longitudinal/3-latent_posterior_quality_characterization \
    --figures-dir ./figures/3-latent-posterior-quality
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

# Consistent manuscript-scale typography and drawing order.
MAIN_BLUE = "#1f77b4"
SECONDARY_GRAY = "0.75"
STATS_FONT = 7.0
LABEL_FONT = 8.0
TICK_FONT = 7.0
LEGEND_FONT = 7.0
ANNOT_FONT = 8.0
ZORDER_SECONDARY = 0
ZORDER_BAND = 2
ZORDER_MAIN = 5
ZORDER_REFERENCE = 6
ZORDER_TEXT = 10

# -----------------------------------------------------------------------------
# Individual figure sizes in inches: (width, height).
# Edit any entry to resize only that panel. Defaults reproduce the current script.
# -----------------------------------------------------------------------------
FIGURE_SIZES_INCHES = {
    # Main manuscript panels (Fig. 2D-F)
    "fig2D_abundance_vs_posterior_width": (3.2, 2.4),
    "fig2E_abundance_vs_replicate_overlap": (3.2, 2.4),
    "fig2F_state_vs_union_detectability": (3.4, 2.4),
    "fig2_DEF_latent_state_inference": (8.4, 2.7),

    # Supplementary Fig. S1 panels
    "suppS1A_posterior_width_distribution": (3.2, 2.4),
    "suppS1B_replicate_overlap_distribution": (3.2, 2.4),
    "suppS1C_detectability_correction_distribution": (3.2, 2.4),
    "suppS1D_detectability_correction_vs_posterior_width": (3.2, 2.4),
    "suppS1E_detectability_correction_vs_state_detectability": (3.2, 2.4),
    "suppS1F_qc_metric_correlation": (4.2, 3.0),
}



# Optional manual axis settings for each output. Values set to None preserve
# the data-driven defaults established inside the plotting function.
AXIS_CONFIG = {
    "fig2D_abundance_vs_posterior_width": {
        "xlim": None, "ylim": None, "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
    "fig2E_abundance_vs_replicate_overlap": {
        "xlim": None, "ylim": (0.0, 1.0), "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
    "fig2F_state_vs_union_detectability": {
        "xlim": None, "ylim": None, "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
    "suppS1A_posterior_width_distribution": {
        "xlim": None, "ylim": None, "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
    "suppS1B_replicate_overlap_distribution": {
        "xlim": (0.0, 1.0), "ylim": None, "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
    "suppS1C_detectability_correction_distribution": {
        "xlim": None, "ylim": None, "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
    "suppS1D_detectability_correction_vs_posterior_width": {
        "xlim": None, "ylim": None, "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
    "suppS1E_detectability_correction_vs_state_detectability": {
        "xlim": None, "ylim": None, "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
    "suppS1F_qc_metric_correlation": {
        "xlim": None, "ylim": None, "xscale": None, "yscale": None,
        "xticks": None, "yticks": None,
    },
}



def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path



def clear_plot_outputs(output_dir: Path) -> None:
    """Remove stale rendered panels from a previous plotter version."""
    for pattern in ("*.pdf", "*.png", "*.svg"):
        for path in output_dir.glob(pattern):
            path.unlink()


def set_rcparams(font: str, font_size: float) -> None:
    plt.rcParams.update({
        "font.family": font,
        "font.size": font_size,
        "axes.labelsize": LABEL_FONT,
        "xtick.labelsize": TICK_FONT,
        "ytick.labelsize": TICK_FONT,
        "legend.fontsize": LEGEND_FONT,
        "axes.linewidth": 0.8,
        "axes.prop_cycle": plt.cycler(color=[MAIN_BLUE]),
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def clean_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def apply_axis_config(ax: plt.Axes, stem: str) -> None:
    """Apply optional final axis overrides from AXIS_CONFIG."""
    cfg = AXIS_CONFIG.get(stem, {})

    if cfg.get("xscale"):
        ax.set_xscale(str(cfg["xscale"]))
    if cfg.get("yscale"):
        ax.set_yscale(str(cfg["yscale"]))
    if cfg.get("xlim") is not None:
        ax.set_xlim(*cfg["xlim"])
    if cfg.get("ylim") is not None:
        ax.set_ylim(*cfg["ylim"])
    if cfg.get("xticks") is not None:
        ax.set_xticks(cfg["xticks"])
    if cfg.get("yticks") is not None:
        ax.set_yticks(cfg["yticks"])


def save(fig: plt.Figure, output_dir: Path, stem: str, dpi: int) -> None:
    # Keep a small, uniform padding around every exported panel.
    fig.savefig(output_dir / (stem + ".pdf"), bbox_inches="tight", pad_inches=0.04, facecolor="white")
    fig.savefig(output_dir / (stem + ".png"), dpi=dpi, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    plt.close(fig)


def systematic_sample(parquet_path: Path, columns, sample_n: int, seed: int) -> pd.DataFrame:
    base = pl.scan_parquet(parquet_path).select(columns)
    n_rows = int(base.select(pl.len()).collect(engine="streaming").item())
    if n_rows <= sample_n:
        return base.collect(engine="streaming").to_pandas()
    stride = max(n_rows // sample_n, 1)
    offset = int(seed) % stride
    return (
        pl.scan_parquet(parquet_path)
        .select(columns)
        .with_row_index("_row_index")
        .filter((pl.col("_row_index") % stride) == offset)
        .select(columns)
        .head(sample_n)
        .collect(engine="streaming")
        .to_pandas()
    )


def finite_series(values) -> pd.Series:
    return pd.Series(values).replace([np.inf, -np.inf], np.nan).dropna()


def add_summary_margin(fig: plt.Figure, text: str) -> None:
    """Place a compact summary box in the reserved right margin."""
    fig.text(
        0.975, 0.94, text, ha="left", va="top",
        fontsize=STATS_FONT, fontfamily="monospace", linespacing=1.10,
        zorder=ZORDER_TEXT,
        bbox=dict(
            facecolor="white", edgecolor="0.82", linewidth=0.45,
            alpha=0.80, boxstyle="round,pad=0.22"
        ),
    )


def hist_plot(values, xlabel: str, stem: str, output_dir: Path, dpi: int,
              bins: int = 70, xlim=None, annotate: bool = True,
              full_summary=None) -> None:
    values = finite_series(values)
    # The histogram may use a plotting sample, but manuscript annotations come
    # from the full-dataset summary generated by 3-latent_posterior_quality_characterization.py.
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    fig.subplots_adjust(left=0.15, right=0.76, bottom=0.20, top=0.96)
    ax.hist(values, bins=bins, density=True, color=MAIN_BLUE, edgecolor="none", linewidth=0, zorder=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    if xlim is not None:
        ax.set_xlim(*xlim)
    if annotate and len(values):
        if full_summary is not None:
            med = float(full_summary["median"])
            q25 = float(full_summary["q25"])
            q75 = float(full_summary["q75"])
            q95 = float(full_summary["q95"])
        else:
            q25, med, q75, q95 = np.quantile(values, [0.25, 0.5, 0.75, 0.95])
        text = ("Median  %7.3g\n"
                "IQR     %7.3g–%.3g\n"
                "P95     %7.3g") % (med, q25, q75, q95)
        add_summary_margin(fig, text)
    apply_axis_config(ax, stem)
    clean_axes(ax)
    save(fig, output_dir, stem, dpi)


def interval_curve(df: pd.DataFrame, left: str, right: str, median: str,
                   q25: str, q75: str, xlabel: str, ylabel: str, stem: str,
                   output_dir: Path, dpi: int, ylim=None, show_n: bool = False,
                   vertical_ref=None, horizontal_ref=None,
                   annotation_text: Optional[str] = None) -> None:
    d = df.sort_values(left).copy()
    x = ((d[left] + d[right]) / 2.0).to_numpy(float)
    med = d[median].to_numpy(float)
    lo = d[q25].to_numpy(float)
    hi = d[q75].to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(med) & np.isfinite(lo) & np.isfinite(hi)
    x, med, lo, hi = x[valid], med[valid], lo[valid], hi[valid]
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    ax.fill_between(x, lo, hi, color=MAIN_BLUE, alpha=0.25, linewidth=0, zorder=ZORDER_BAND)
    ax.plot(x, med, color=MAIN_BLUE, marker="o", markersize=3, linewidth=0.9, zorder=ZORDER_MAIN)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if vertical_ref is not None and np.isfinite(vertical_ref):
        ax.axvline(vertical_ref, linestyle="--", linewidth=0.7, color="0.35",
                   zorder=ZORDER_REFERENCE)
    if horizontal_ref is not None and np.isfinite(horizontal_ref):
        ax.axhline(horizontal_ref, linestyle="--", linewidth=0.7, color="0.35",
                   zorder=ZORDER_REFERENCE)
    if annotation_text:
        ax.text(
            0.5, 0.8, annotation_text, transform=ax.transAxes,
            ha="left", va="top", fontsize=STATS_FONT,
            bbox=dict(facecolor="white", edgecolor="0.82", linewidth=0.45,
                      alpha=0.82, boxstyle="round,pad=0.22"),
            zorder=ZORDER_TEXT,
        )
    if show_n and "n_clonotype_timepoints" in d.columns:
        n = d.loc[valid, "n_clonotype_timepoints"].to_numpy(float)
        if len(n):
            ax2 = ax.twinx()
            ax2.plot(x, n, linewidth=1.0, alpha=0.60, color=SECONDARY_GRAY,
                     linestyle=":", zorder=ZORDER_SECONDARY)
            ax2.patch.set_visible(False)
            ax.set_zorder(2)
            ax.patch.set_visible(False)
            ax2.set_yscale("log")
            ax2.set_ylabel("Clonotype-timepoints")
            ax2.spines["top"].set_visible(False)
    apply_axis_config(ax, stem)
    clean_axes(ax)
    fig.tight_layout(pad=1.0)
    save(fig, output_dir, stem, dpi)


def lowess_curve(df: pd.DataFrame, x: str, y: str, xlabel: str, ylabel: str,
                 stem: str, output_dir: Path, dpi: int, frac: float = 0.25,
                 max_points: int = 50000, grid_n: int = 120,
                 full_spearman: Optional[float] = None) -> None:
    """Plot a deterministic locally weighted linear smoother without extra dependencies.

    The fit uses tricube weights over the nearest ``frac`` fraction of observations.
    A systematic subset is used for computational stability when the sampled QC table
    is very large. The shaded band is a descriptive local IQR, not a confidence interval.
    """
    d = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna().sort_values(x)
    if len(d) < 20:
        raise ValueError("LOWESS panel requires at least 20 finite observations.")
    if len(d) > max_points:
        idx = np.linspace(0, len(d) - 1, max_points, dtype=int)
        d = d.iloc[idx]

    xv = d[x].to_numpy(float)
    yv = d[y].to_numpy(float)
    x_grid = np.linspace(np.quantile(xv, 0.01), np.quantile(xv, 0.99), grid_n)
    k = max(int(np.ceil(frac * len(xv))), 20)

    smooth = np.empty_like(x_grid)
    q25 = np.empty_like(x_grid)
    q75 = np.empty_like(x_grid)
    for i, x0 in enumerate(x_grid):
        dist = np.abs(xv - x0)
        bandwidth = np.partition(dist, k - 1)[k - 1]
        if not np.isfinite(bandwidth) or bandwidth <= 0:
            bandwidth = np.max(dist)
        if bandwidth <= 0:
            smooth[i] = np.median(yv)
            q25[i], q75[i] = np.quantile(yv, [0.25, 0.75])
            continue
        mask = dist <= bandwidth
        xx = xv[mask]
        yy = yv[mask]
        u = dist[mask] / bandwidth
        w = (1.0 - u ** 3) ** 3
        # Weighted local linear regression centered at x0.
        z = xx - x0
        sw = w.sum()
        swz = np.sum(w * z)
        swz2 = np.sum(w * z * z)
        swy = np.sum(w * yy)
        swzy = np.sum(w * z * yy)
        denom = sw * swz2 - swz * swz
        if sw <= 0 or abs(denom) < 1e-14:
            smooth[i] = np.average(yy, weights=w) if sw > 0 else np.median(yy)
        else:
            smooth[i] = (swy * swz2 - swzy * swz) / denom
        q25[i], q75[i] = np.quantile(yy, [0.25, 0.75])

    rho = (float(full_spearman) if full_spearman is not None
           else pd.Series(xv).corr(pd.Series(yv), method="spearman"))
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    ax.fill_between(x_grid, q25, q75, color=MAIN_BLUE, alpha=0.18,
                    linewidth=0, zorder=ZORDER_BAND)
    ax.plot(x_grid, smooth, color=MAIN_BLUE, linewidth=1.1, zorder=ZORDER_MAIN)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.text(0.03, 0.96, "Spearman ρ = %.3f" % rho,
            transform=ax.transAxes, ha="left", va="top", fontsize=STATS_FONT,
            bbox=dict(facecolor="white", edgecolor="0.82", linewidth=0.45,
                      alpha=0.80, boxstyle="round,pad=0.22"))
    apply_axis_config(ax, stem)
    clean_axes(ax)
    fig.tight_layout(pad=1.0)
    save(fig, output_dir, stem, dpi)


def state_vs_union_plot(df: pd.DataFrame, output_dir: Path, dpi: int,
                        stem: str = "fig2F_state_vs_union_detectability",
                        annotation_text: Optional[str] = None) -> None:
    """Compare posterior-integrated state detectability with the simple union approximation."""
    cols = ["cqc_union_detectability", "cqc_state_detectability"]
    d = df[cols].replace([np.inf, -np.inf], np.nan).dropna()
    if d.empty:
        raise ValueError("State-vs-union panel requires finite detectability values.")
    x = d[cols[0]].to_numpy(float)
    y = d[cols[1]].to_numpy(float)
    lo = float(min(np.quantile(x, 0.001), np.quantile(y, 0.001)))
    hi = float(max(np.quantile(x, 0.999), np.quantile(y, 0.999)))
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    hb = ax.hexbin(x, y, gridsize=60, mincnt=1, bins="log")
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=0.8, color="0.35",
            zorder=ZORDER_REFERENCE)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Replicate-union detectability")
    ax.set_ylabel("Latent-state detectability")
    if annotation_text:
        ax.text(
            0.03, 0.96, annotation_text, transform=ax.transAxes,
            ha="left", va="top", fontsize=STATS_FONT,
            bbox=dict(facecolor="white", edgecolor="0.82", linewidth=0.45,
                      alpha=0.82, boxstyle="round,pad=0.22"),
            zorder=ZORDER_TEXT,
        )
    apply_axis_config(ax, stem)
    clean_axes(ax)
    cbar = fig.colorbar(hb, ax=ax)
    cbar.set_label("Count (log scale)")
    fig.tight_layout(pad=1.0)
    save(fig, output_dir, stem, dpi)



def make_main_strip(
    sample: pd.DataFrame,
    abundance: pd.DataFrame,
    x_median: float,
    median_width: float,
    median_bc: float,
    median_abs_delta: float,
    output_dir: Path,
    dpi: int,
    bc_reference: float,
) -> None:
    """Create the vector Fig. 2D-F strip used for final manuscript assembly."""
    fig, axes = plt.subplots(
        1, 3, figsize=FIGURE_SIZES_INCHES["fig2_DEF_latent_state_inference"]
    )

    # D: abundance-dependent posterior width.
    d = abundance.sort_values("x_bin_left").copy()
    x = ((d["x_bin_left"] + d["x_bin_right"]) / 2.0).to_numpy(float)
    med = d["posterior_width95_median"].to_numpy(float)
    lo = d["posterior_width95_q25"].to_numpy(float)
    hi = d["posterior_width95_q75"].to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(med) & np.isfinite(lo) & np.isfinite(hi)
    ax = axes[0]
    ax.fill_between(x[valid], lo[valid], hi[valid], color=MAIN_BLUE, alpha=0.25,
                    linewidth=0, zorder=ZORDER_BAND)
    ax.plot(x[valid], med[valid], color=MAIN_BLUE, marker="o", markersize=3,
            linewidth=0.9, zorder=ZORDER_MAIN)
    ax.set_xlabel("Latent log-frequency bin midpoint")
    ax.set_ylabel("Posterior 95% width (x)")
    ax.text(0.1, 0.96, f"overall median = {median_width:.2f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=STATS_FONT,
            bbox=dict(facecolor="white", edgecolor="0.82", linewidth=0.45,
                      alpha=0.82, boxstyle="round,pad=0.22"))
    apply_axis_config(ax, "fig2D_abundance_vs_posterior_width")
    clean_axes(ax)

    # E: abundance-dependent replicate posterior overlap.
    ax = axes[1]
    med = d["replicate_bc_median"].to_numpy(float)
    lo = d["replicate_bc_q25"].to_numpy(float)
    hi = d["replicate_bc_q75"].to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(med) & np.isfinite(lo) & np.isfinite(hi)
    ax.fill_between(x[valid], lo[valid], hi[valid], color=MAIN_BLUE, alpha=0.25,
                    linewidth=0, zorder=ZORDER_BAND)
    ax.plot(x[valid], med[valid], color=MAIN_BLUE, marker="o", markersize=3,
            linewidth=0.9, zorder=ZORDER_MAIN)
    ax.set_xlabel("Latent log-frequency bin midpoint")
    ax.set_ylabel("Replicate posterior overlap (BC)")
    ax.text(0.03, 0.96, f"overall median BC = {median_bc:.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=STATS_FONT,
            bbox=dict(facecolor="white", edgecolor="0.82", linewidth=0.45,
                      alpha=0.82, boxstyle="round,pad=0.22"))
    apply_axis_config(ax, "fig2E_abundance_vs_replicate_overlap")
    clean_axes(ax)

    # F: state detectability versus replicate union.
    ax = axes[2]
    cols = ["cqc_union_detectability", "cqc_state_detectability"]
    dd = sample[cols].replace([np.inf, -np.inf], np.nan).dropna()
    xx = dd[cols[0]].to_numpy(float)
    yy = dd[cols[1]].to_numpy(float)
    lo_xy = float(min(np.quantile(xx, 0.001), np.quantile(yy, 0.001)))
    hi_xy = float(max(np.quantile(xx, 0.999), np.quantile(yy, 0.999)))
    hb = ax.hexbin(xx, yy, gridsize=55, mincnt=1, bins="log")
    ax.plot([lo_xy, hi_xy], [lo_xy, hi_xy], linestyle="--", linewidth=0.8,
            color="0.35", zorder=ZORDER_REFERENCE)
    ax.set_xlim(lo_xy, hi_xy)
    ax.set_ylim(lo_xy, hi_xy)
    ax.set_xlabel("Replicate-union detectability")
    ax.set_ylabel("Latent-state detectability")
    ax.text(0.03, 0.96, f"median |Δ| = {median_abs_delta:.5f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=STATS_FONT,
            bbox=dict(facecolor="white", edgecolor="0.82", linewidth=0.45,
                      alpha=0.82, boxstyle="round,pad=0.22"))
    apply_axis_config(ax, "fig2F_state_vs_union_detectability")
    clean_axes(ax)
    cbar = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Count (log scale)")

    for ax, label in zip(axes, ["D", "E", "F"]):
        ax.text(-0.16, 1.05, label, transform=ax.transAxes, ha="left", va="top",
                fontsize=9.0, fontweight="bold")

    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.22, top=0.96, wspace=0.42)
    save(fig, output_dir, "fig2_DEF_latent_state_inference", dpi)

def load_metric_summaries(path: Path) -> dict:
    table = pd.read_csv(path)
    if "metric" not in table.columns:
        raise ValueError("Metric summary is missing the 'metric' column: %s" % path)
    required = {"q25", "median", "q75", "q95"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(
            "Metric summary lacks full-dataset annotation fields %s. "
            "Rerun 3-latent_posterior_quality_characterization.py." % sorted(missing)
        )
    return {str(row["metric"]): row.to_dict() for _, row in table.iterrows()}


def load_key_relationships(path: Path) -> dict:
    table = pd.read_csv(path)
    required = {"relationship", "spearman_rho"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError("Key-relationship table is missing columns: %s" % sorted(missing))
    return {
        str(row["relationship"]): float(row["spearman_rho"])
        for _, row in table.iterrows()
    }


def load_spearman_matrix(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    if "metric" not in table.columns:
        raise ValueError("Spearman matrix is missing the 'metric' column: %s" % path)
    corr = table.set_index("metric")
    corr = corr.apply(pd.to_numeric, errors="coerce")
    corr = corr.loc[corr.index, corr.index]
    return corr


def correlation_heatmap(corr: pd.DataFrame, output_dir: Path, dpi: int,
                         stem: str = "suppS1F_qc_metric_correlation") -> None:
    labels = {
        "cqc_posterior_width95_x": "Width",
        "cqc_posterior_entropy": "Entropy",
        "cqc_state_detectability": "State detectability",
        "cqc_union_detectability": "Union detectability",
        "cqc_replicate_bhattacharyya": "Overlap",
        "cqc_effective_posterior_support_simpson": "Inv. Simpson",
    }
    keep = [c for c in labels if c in corr.index and c in corr.columns]
    corr = corr.loc[keep, keep]
    values = corr.to_numpy(float)
    # Keep only the lower triangle and omit the diagonal.
    mask = np.triu(np.ones_like(values, dtype=bool), k=0)
    masked = np.ma.array(values, mask=mask)

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("white")
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    im = ax.imshow(masked, vmin=-1, vmax=1, aspect="equal", cmap=cmap)
    ticks = np.arange(len(keep))
    names = [labels[c] for c in keep]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=TICK_FONT)
    ax.set_yticklabels(names, fontsize=TICK_FONT)
    for i in range(len(corr)):
        for j in range(len(corr)):
            if i > j:
                value = corr.iloc[i, j]
                text_color = "white" if abs(value) >= 0.55 else "black"
                ax.text(j, i, "%.2f" % value, ha="center", va="center",
                        fontsize=ANNOT_FONT, color=text_color, zorder=ZORDER_TEXT)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Spearman correlation")
    apply_axis_config(ax, stem)
    fig.tight_layout(pad=1.0)
    save(fig, output_dir, stem, dpi)

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Plot ClonoDynamics Step-3 latent posterior-quality characterization outputs.")
    p.add_argument("--analysis-dir", type=Path, required=True)
    p.add_argument("--figures-dir", type=Path, default=None,
                   help="Default: <analysis-dir>/figures_latent_posterior_quality")
    p.add_argument("--sample-n", type=int, default=250000)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--font", default="Arial")
    p.add_argument("--font-size", type=float, default=8.0)
    p.add_argument("--bc-reference", type=float, default=0.95,
                   help="Horizontal visual reference in abundance-vs-overlap panel.")
    p.add_argument("--replicate-delta-reference", type=float, default=2.0,
                   help="Unlabelled vertical reference in replicate-difference panel.")
    p.add_argument("--show-bin-counts", action="store_true",
                   help="Add a secondary log-scale count curve to abundance panels.")
    p.add_argument("--include-qc-panels", action="store_true",
                   help="Retained for compatibility; subject/pair plots remain diagnostic only.")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    analysis_dir = args.analysis_dir
    root = ensure_dir(args.figures_dir or analysis_dir / "figures_latent_posterior_quality")
    main_dir = ensure_dir(root / "main")
    supplementary_dir = ensure_dir(root / "supplementary")
    # Prevent removed or renumbered panels from surviving from earlier runs.
    clear_plot_outputs(main_dir)
    clear_plot_outputs(supplementary_dir)
    set_rcparams(args.font, args.font_size)

    paths = {
        "parquet": analysis_dir / "clonotype_latent_posterior_qc.parquet",
        "metric_summary": analysis_dir / "clonotype_latent_posterior_qc_metric_summary.csv",
        "key_relationships": analysis_dir / "clonotype_latent_posterior_qc_key_relationships.csv",
        "spearman_matrix": analysis_dir / "qc_metric_spearman_correlations.csv",
        "abundance": analysis_dir / "clonotype_latent_posterior_qc_by_abundance.csv",
        "abundance_reference": analysis_dir / "clonotype_latent_posterior_qc_abundance_reference.csv",
    }
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing Step-3 outputs. Rerun 3-latent_posterior_quality_characterization.py first:\n" +
            "\n".join(missing)
        )

    schema = pl.scan_parquet(paths["parquet"]).collect_schema().names()
    sample_cols = [
        "cqc_posterior_width95_x",
        "cqc_state_detectability",
        "cqc_union_detectability",
        "cqc_abs_delta_detectability",
        "cqc_replicate_bhattacharyya",
    ]
    sample_cols = [c for c in sample_cols if c in schema]
    sample = systematic_sample(paths["parquet"], sample_cols, args.sample_n, args.seed)
    metric_summaries = load_metric_summaries(paths["metric_summary"])
    key_relationships = load_key_relationships(paths["key_relationships"])
    full_spearman = load_spearman_matrix(paths["spearman_matrix"])

    def full_summary(metric: str):
        return metric_summaries.get(metric)

    abundance = pd.read_csv(paths["abundance"])
    abundance_ref = pd.read_csv(paths["abundance_reference"])
    x_median = float(abundance_ref.loc[0, "x_median"])

    width_summary = full_summary("cqc_posterior_width95_x")
    bc_summary = full_summary("cqc_replicate_bhattacharyya")
    delta_summary = full_summary("cqc_abs_delta_detectability")
    if width_summary is None or bc_summary is None or delta_summary is None:
        raise ValueError(
            "Required full-dataset summaries for width, overlap, or detectability "
            "discrepancy are missing from the Step-3 metric summary."
        )

    median_width = float(width_summary["median"])
    median_bc = float(bc_summary["median"])
    median_abs_delta = float(delta_summary["median"])

    # ------------------------------------------------------------------
    # Main Figure 2D-F
    # ------------------------------------------------------------------
    interval_curve(
        abundance, "x_bin_left", "x_bin_right",
        "posterior_width95_median", "posterior_width95_q25", "posterior_width95_q75",
        "Latent log-frequency bin midpoint", "Posterior 95% width (x)",
        "fig2D_abundance_vs_posterior_width", main_dir, args.dpi,
        show_n=args.show_bin_counts,
        annotation_text=f"overall median = {median_width:.2f}",
    )
    interval_curve(
        abundance, "x_bin_left", "x_bin_right",
        "replicate_bc_median", "replicate_bc_q25", "replicate_bc_q75",
        "Latent log-frequency bin midpoint", "Replicate posterior overlap (BC)",
        "fig2E_abundance_vs_replicate_overlap", main_dir, args.dpi, ylim=(0, 1),
        show_n=args.show_bin_counts,
        annotation_text=f"overall median BC = {median_bc:.3f}",
    )
    state_vs_union_plot(
        sample, main_dir, args.dpi,
        "fig2F_state_vs_union_detectability",
        annotation_text=f"median |Δ| = {median_abs_delta:.5f}",
    )
    make_main_strip(
        sample=sample,
        abundance=abundance,
        x_median=x_median,
        median_width=median_width,
        median_bc=median_bc,
        median_abs_delta=median_abs_delta,
        output_dir=main_dir,
        dpi=args.dpi,
        bc_reference=args.bc_reference,
    )

    # ------------------------------------------------------------------
    # Supplementary Figure S1A-F: extended latent-state diagnostics
    # ------------------------------------------------------------------
    hist_plot(
        sample["cqc_posterior_width95_x"], "Posterior 95% width (x)",
        "suppS1A_posterior_width_distribution", supplementary_dir, args.dpi,
        full_summary=width_summary,
    )
    hist_plot(
        sample["cqc_replicate_bhattacharyya"],
        "Replicate posterior overlap (Bhattacharyya coefficient)",
        "suppS1B_replicate_overlap_distribution", supplementary_dir, args.dpi,
        xlim=(0, 1), full_summary=bc_summary,
    )
    hist_plot(
        sample["cqc_abs_delta_detectability"],
        "|State detectability − union approximation|",
        "suppS1C_detectability_correction_distribution", supplementary_dir, args.dpi,
        full_summary=delta_summary,
    )
    lowess_curve(
        sample, "cqc_posterior_width95_x", "cqc_abs_delta_detectability",
        "Posterior 95% width (x)",
        "|State detectability −\nunion approximation|",
        "suppS1D_detectability_correction_vs_posterior_width",
        supplementary_dir, args.dpi,
        full_spearman=key_relationships["posterior_width95_vs_abs_delta_detectability"],
    )
    lowess_curve(
        sample, "cqc_state_detectability", "cqc_abs_delta_detectability",
        "Latent-state detectability",
        "|State detectability −\nunion approximation|",
        "suppS1E_detectability_correction_vs_state_detectability",
        supplementary_dir, args.dpi,
        full_spearman=key_relationships["state_detectability_vs_abs_delta_detectability"],
    )
    correlation_heatmap(
        full_spearman, supplementary_dir, args.dpi,
        "suppS1F_qc_metric_correlation",
    )

    print("[Step-3 plot] Main Fig. 2D-F:", main_dir)
    print("[Step-3 plot] Supplementary Fig. S1:", supplementary_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())