#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_step03_latent_state_diagnostics.py
========================================================

Plot-only companion to ClonoDynamics Step 3:

    3-latent_state_and_observation_model_diagnostics.py

Scientific role
---------------
Step 3 is a diagnostic layer. It characterizes the model-based latent consensus
representation, technical-replicate posterior agreement, realized positive-
replicate support, and model-based observation diagnostics without filtering
or redefining clonotype states.

This plotter is aligned with the replicate-resolved ClonoDynamics architecture.
It therefore does NOT present the latent consensus as a universal dynamical
state and does NOT treat posterior-predictive detectability as biological
presence. In the production pipeline:

    - x_latent is a model-based consensus representation used primarily for
      conditioning geometry and uncertainty characterization;
    - primary forward dynamics are estimated downstream from replicate-
      decoupled observed AB/BA measurements;
    - primary fluctuations are estimated downstream from cross-replicate
      observed displacement covariance on common4;
    - realized technical detection, model-based posterior-predictive
      detectability, and p_value-based operational observability remain distinct
      information layers.

The plotter never refits Step 2, never changes Step-3 tables, never scores or
filters clonotypes, and never generates Step-4 inputs.

Production figure emphasis
--------------------------
The figure set treats the state-versus-factorized-union comparison as a
model-structure diagnostic rather than an external validation or correction: p_detect_state is compared with a
factorized-union quantity generated from the same fitted observation model,
so the comparison is a model-structure diagnostic rather than an independent
validation or correction.

More importantly, the ONE_POSITIVE analysis identifies a strong dependence of
the latent-versus-positive representation shift on relative replicate depth.
That result is central to the methodological pivot away from treating the
paired latent consensus as the universally preferred abundance representation.
The production plotter emphasizes the paired ONE_POSITIVE depth
coupling to the main figure and moves the factorized-union diagnostic to the
Supplement.

Production outputs
------------------
Only single-panel figures are generated. Figure assembly is intentionally left
to a separate manuscript-assembly step.

Main manuscript candidates (continuing Figure 2 after Step-1 panels A-C):

    fig2D_abundance_vs_posterior_width
        Abundance dependence of the paired latent-posterior 95% width.

    fig2E_abundance_vs_replicate_overlap
        Abundance dependence of the exact discrete Bhattacharyya overlap
        between the two single-replicate latent posteriors.

    fig2F_one_positive_paired_depth_coupling
        Across the 57 subject-timepoints, log(depth_rep1/depth_rep2) versus
        the difference between the median latent-minus-positive shift when
        replicate 1 is the zero side and when replicate 2 is the zero side.
        Full-data Spearman, Pearson and R-squared annotations are read from the
        Step-3 paired diagnostic summary.

Supplementary diagnostic panels:

    suppS1A_posterior_width_distribution
        Global distribution of latent-posterior 95% width.

    suppS1B_replicate_overlap_distribution
        Global distribution of exact replicate-posterior Bhattacharyya overlap.

    suppS1C_latent_vs_positive_by_detection_support
        ECDF of x_latent - x_positive for ONE_POSITIVE versus TWO_POSITIVE
        states. Plotting curves use a deterministic sample; sample-independent
        counts and medians are read from the full-data positive-state summary.

    suppS1D_one_positive_relative_depth_coupling
        Pair-aware ONE_POSITIVE diagnostic across 114 subject x time x zero-side
        units: log(depth_zero/depth_positive) versus the unit median
        x_latent - x_positive. The full-data Spearman statistic comes from the
        aggregated Step-3 diagnostic, so clonotype rows are not treated as
        independent depth observations.

    suppS1E_state_vs_factorized_union_detectability
        Model-structure diagnostic comparing posterior-integrated state
        detectability with the factorized union of replicate-specific marginal
        predictive detectabilities. This is not an external calibration.

    suppS1F_qc_metric_correlation
        Full-data Spearman correlation matrix for the selected Step-3
        diagnostic metrics, using canonical factorized-union terminology.

Canonical Step-3 inputs
-----------------------
The analysis directory must contain:

    clonotype_latent_posterior_qc.parquet
    clonotype_latent_posterior_qc_metric_summary.csv
    clonotype_latent_posterior_qc_by_abundance.csv
    clonotype_latent_posterior_qc_by_positive_state.csv
    one_positive_zero_depth_by_subject_time_replica.csv
    one_positive_zero_depth_diagnostic.csv
    one_positive_paired_replica_depth_by_subject_time.csv
    one_positive_paired_replica_depth_diagnostic.csv
    state_detectability_factorized_union_diagnostic.csv
    qc_metric_spearman_correlations.csv

The plotter deliberately requires canonical Step-3 v11 field names rather than
legacy aliases such as cqc_union_detectability or cqc_abs_delta_detectability.
A schema failure should therefore be interpreted as evidence that an older
Step-3 output directory is being used.

Sampling policy
---------------
The 16-million-state QC parquet is sampled deterministically only for graphical
rendering of distributions and the factorized-union hexbin. Manuscript-facing
summary statistics are read from the complete Step-3 CSV summaries. The
ONE_POSITIVE depth-coupling panels use already aggregated subject-time analysis
units and therefore do not require clonotype-level sampling.

Output directory
----------------
By default:

    <analysis-dir>/figures_latent_state_and_observation_model_diagnostics/
        main/
        supplementary/

A different root can be supplied with --figures-dir.

Typical usage
-------------
python3 ./code/plotting/plot_step03_latent_state_diagnostics.py \
    --analysis-dir \
    ./dataset_longitudinal_clonodynamics_results/longitudinal/3-latent_state_and_observation_model_diagnostics \
    --figures-dir ./figures/3-latent-state-and-observation-model-diagnostics \
    --png

Interpretive boundary
---------------------
These figures document properties and limitations of the Step-2/Step-3 state
representation layer. They do not estimate longitudinal forward drift,
fluctuation covariance or temporal scaling and should not be used to infer
biological presence/absence from technical observation states.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl


SCRIPT_VERSION = "v13-layout-box-position-2026-09-09"

# -----------------------------------------------------------------------------
# Publication geometry and styling
# -----------------------------------------------------------------------------

MAIN_BLUE = "#1f77b4"
SECONDARY_ORANGE = "#ff7f0e"
SECONDARY_GRAY = "0.70"
STATS_FONT = 7.0
LABEL_FONT = 8.0
TICK_FONT = 7.0
LEGEND_FONT = 7.0
ANNOT_FONT = 8.0
ZORDER_BAND = 2
ZORDER_MAIN = 5
ZORDER_REFERENCE = 6
ZORDER_TEXT = 10

FIGURE_SIZES_INCHES = {
    # Main manuscript candidates.
    "fig2D_abundance_vs_posterior_width": (3.2, 2.4),
    "fig2E_abundance_vs_replicate_overlap": (3.2, 2.4),
    "fig2F_one_positive_paired_depth_coupling": (3.2, 2.4),

    # Supplementary panels.
    "suppS1A_posterior_width_distribution": (3.2, 2.4),
    "suppS1B_replicate_overlap_distribution": (3.2, 2.4),
    "suppS1C_latent_vs_positive_by_detection_support": (3.2, 2.4),
    "suppS1D_one_positive_relative_depth_coupling": (3.2, 2.4),
    "suppS1E_state_vs_factorized_union_detectability": (3.2, 2.4),
    "suppS1F_qc_metric_correlation": (3.2, 2.4),
}

AXIS_CONFIG = {
    stem: {
        "xlim": None,
        "ylim": None,
        "xscale": None,
        "yscale": None,
        "xticks": None,
        "yticks": None,
    }
    for stem in FIGURE_SIZES_INCHES
}
AXIS_CONFIG["fig2E_abundance_vs_replicate_overlap"]["ylim"] = (0.0, 1.0)
AXIS_CONFIG["suppS1B_replicate_overlap_distribution"]["xlim"] = (0.0, 1.0)


# Per-panel positions for internal statistics/annotation boxes.
# Coordinates are in axes units: x=0 is left, x=1 is right;
# y=0 is bottom, y=1 is top.  Keeping these settings here makes it possible
# to move a box without touching the plotting logic.
BOX_CONFIG = {
    # Main panels: keep boxes away from the principal curves/fit line.
    "fig2D_abundance_vs_posterior_width": {
        "x": 0.62, "y": 0.08, "ha": "left", "va": "bottom",
    },
    "fig2E_abundance_vs_replicate_overlap": {
        "x": 0.04, "y": 0.08, "ha": "left", "va": "bottom",
    },
    "fig2F_one_positive_paired_depth_coupling": {
        "x": 0.98, "y": 0.96, "ha": "right", "va": "top",
    },

    # Histogram summaries: the upper-left region is visually cleaner for the
    # current distributions than the upper-right tail.
    "suppS1A_posterior_width_distribution": {
        "x": 0.03, "y": 0.96, "ha": "left", "va": "top",
    },
    "suppS1B_replicate_overlap_distribution": {
        "x": 0.03, "y": 0.96, "ha": "left", "va": "top",
    },

    # ONE_POSITIVE and detectability diagnostics.
    "suppS1D_one_positive_relative_depth_coupling": {
        "x": 1, "y": 0.98, "ha": "right", "va": "top",
    },
    "suppS1E_state_vs_factorized_union_detectability": {
        "x": 0.03, "y": 0.96, "ha": "left", "va": "top",
    },
}


def get_box_config(stem: str, default=None) -> dict:
    """Return per-panel annotation-box coordinates and alignment."""
    if default is None:
        default = {"x": 0.03, "y": 0.96, "ha": "left", "va": "top"}
    out = dict(default)
    out.update(BOX_CONFIG.get(stem, {}))
    return out


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def clear_plot_outputs(output_dir: Path) -> None:
    """Remove stale rendered panels from earlier plotter versions."""
    for pattern in ("*.pdf", "*.png", "*.svg"):
        for path in output_dir.glob(pattern):
            path.unlink()


def set_rcparams(font: str, font_size: float, dpi: int) -> None:
    plt.rcParams.update(
        {
            "font.family": font,
            "font.size": font_size,
            "axes.labelsize": LABEL_FONT,
            "xtick.labelsize": TICK_FONT,
            "ytick.labelsize": TICK_FONT,
            "legend.fontsize": LEGEND_FONT,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.dpi": 150,
            "savefig.dpi": dpi,
            "savefig.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
        }
    )


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def apply_axis_config(ax: plt.Axes, stem: str) -> None:
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


def save_figure(
    fig: plt.Figure,
    output_dir: Path,
    stem: str,
    args: argparse.Namespace,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    kwargs = dict(bbox_inches="tight", pad_inches=0.04, facecolor="white")
    if args.pdf:
        path = output_dir / f"{stem}.pdf"
        fig.savefig(path, **kwargs)
        print(f"[SAVED] {path}")
    if args.png:
        path = output_dir / f"{stem}.png"
        fig.savefig(path, dpi=int(args.dpi), **kwargs)
        print(f"[SAVED] {path}")
    if args.svg:
        path = output_dir / f"{stem}.svg"
        fig.savefig(path, **kwargs)
        print(f"[SAVED] {path}")
    plt.close(fig)


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Required Step-3 output not found: {path}")
    return path


def require_columns(frame: pd.DataFrame, required: Sequence[str], label: str) -> None:
    missing = [c for c in required if c not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def systematic_sample(
    parquet_path: Path,
    columns: Sequence[str],
    sample_n: int,
    seed: int,
) -> pd.DataFrame:
    """Deterministic systematic sample used only for graphical rendering."""
    scan = pl.scan_parquet(parquet_path)
    schema = set(scan.collect_schema().names())
    missing = [c for c in columns if c not in schema]
    if missing:
        raise ValueError(
            "Step-3 QC parquet is missing canonical production fields: "
            + ", ".join(missing)
        )

    base = scan.select(list(columns))
    n_rows = int(base.select(pl.len()).collect(engine="streaming").item())
    if n_rows <= int(sample_n):
        return base.collect(engine="streaming").to_pandas()

    stride = max(n_rows // int(sample_n), 1)
    offset = int(seed) % stride
    return (
        scan.select(list(columns))
        .with_row_index("_row_index")
        .filter((pl.col("_row_index") % stride) == offset)
        .select(list(columns))
        .head(int(sample_n))
        .collect(engine="streaming")
        .to_pandas()
    )


def finite_series(values) -> pd.Series:
    return pd.Series(values).replace([np.inf, -np.inf], np.nan).dropna()


def load_metric_summaries(path: Path) -> dict:
    table = pd.read_csv(path)
    require_columns(table, ["metric", "q25", "median", "q75", "q95"], path.name)
    return {str(row["metric"]): row.to_dict() for _, row in table.iterrows()}


def load_spearman_matrix(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    require_columns(table, ["metric"], path.name)
    corr = table.set_index("metric")
    corr = corr.apply(pd.to_numeric, errors="coerce")
    common = [x for x in corr.index if x in corr.columns]
    return corr.loc[common, common]


# -----------------------------------------------------------------------------
# Core plotting functions
# -----------------------------------------------------------------------------


def interval_curve(
    df: pd.DataFrame,
    left: str,
    right: str,
    median: str,
    q25: str,
    q75: str,
    xlabel: str,
    ylabel: str,
    stem: str,
    output_dir: Path,
    args: argparse.Namespace,
    annotation_text: Optional[str] = None,
) -> None:
    require_columns(df, [left, right, median, q25, q75], "abundance summary")
    d = df.sort_values(left).copy()
    x = ((d[left] + d[right]) / 2.0).to_numpy(float)
    med = pd.to_numeric(d[median], errors="coerce").to_numpy(float)
    lo = pd.to_numeric(d[q25], errors="coerce").to_numpy(float)
    hi = pd.to_numeric(d[q75], errors="coerce").to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(med) & np.isfinite(lo) & np.isfinite(hi)

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    ax.fill_between(
        x[valid], lo[valid], hi[valid],
        color=MAIN_BLUE, alpha=0.25, linewidth=0, zorder=ZORDER_BAND,
    )
    ax.plot(
        x[valid], med[valid],
        color=MAIN_BLUE, marker="o", markersize=3,
        linewidth=0.9, zorder=ZORDER_MAIN,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if annotation_text:
        box = get_box_config(stem)
        ax.text(
            box["x"], box["y"], annotation_text,
            transform=ax.transAxes,
            ha=box["ha"], va=box["va"], fontsize=STATS_FONT,
            bbox=dict(
                facecolor="white", edgecolor="0.82", linewidth=0.45,
                alpha=0.85, boxstyle="round,pad=0.20",
            ),
            zorder=ZORDER_TEXT,
        )

    if args.show_bin_counts and "n_clonotype_timepoints" in d.columns:
        n = pd.to_numeric(d.loc[valid, "n_clonotype_timepoints"], errors="coerce").to_numpy(float)
        ax2 = ax.twinx()
        ax2.plot(
            x[valid], n,
            linewidth=0.9, alpha=0.55, color=SECONDARY_GRAY,
            linestyle=":", zorder=1,
        )
        ax2.set_yscale("log")
        ax2.set_ylabel("Clonotype-timepoints")
        ax2.spines["top"].set_visible(False)
        ax2.patch.set_visible(False)

    apply_axis_config(ax, stem)
    clean_axes(ax)
    fig.tight_layout(pad=1.0)
    save_figure(fig, output_dir, stem, args)


def histogram_panel(
    values,
    xlabel: str,
    stem: str,
    output_dir: Path,
    args: argparse.Namespace,
    full_summary: dict,
    bins: int = 70,
) -> None:
    x = finite_series(values)
    if x.empty:
        raise ValueError(f"No finite values available for {stem}")

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    ax.hist(x, bins=bins, density=True, color=MAIN_BLUE, edgecolor="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")

    med = float(full_summary["median"])
    q25 = float(full_summary["q25"])
    q75 = float(full_summary["q75"])
    q95 = float(full_summary["q95"])
    ax.axvline(med, color="0.25", linestyle="--", linewidth=0.7)
    box = get_box_config(stem, {"x": 0.97, "y": 0.96, "ha": "right", "va": "top"})
    ax.text(
        box["x"], box["y"],
        f"median = {med:.3g}\nIQR = {q25:.3g}–{q75:.3g}\nP95 = {q95:.3g}",
        transform=ax.transAxes,
        ha=box["ha"], va=box["va"], fontsize=STATS_FONT,
        bbox=dict(
            facecolor="white", edgecolor="0.82", linewidth=0.45,
            alpha=0.85, boxstyle="round,pad=0.20",
        ),
    )

    apply_axis_config(ax, stem)
    clean_axes(ax)
    fig.tight_layout(pad=1.0)
    save_figure(fig, output_dir, stem, args)


def ecdf_by_positive_state(
    sample: pd.DataFrame,
    positive_summary: pd.DataFrame,
    output_dir: Path,
    args: argparse.Namespace,
    stem: str = "suppS1C_latent_vs_positive_by_detection_support",
) -> None:
    require_columns(
        sample,
        ["positive_state_class", "cqc_latent_minus_positive_x"],
        "QC plotting sample",
    )
    require_columns(
        positive_summary,
        [
            "positive_state_class",
            "n_clonotype_timepoints",
            "latent_minus_positive_x_median",
            "fraction_of_all_states",
        ],
        "positive-state summary",
    )

    styles = {
        "ONE_POSITIVE": (MAIN_BLUE, "ONE_POSITIVE"),
        "TWO_POSITIVE": (SECONDARY_ORANGE, "TWO_POSITIVE"),
    }

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    for state in ["ONE_POSITIVE", "TWO_POSITIVE"]:
        d = sample.loc[
            sample["positive_state_class"].astype(str) == state,
            "cqc_latent_minus_positive_x",
        ]
        x = np.sort(finite_series(d).to_numpy(float))
        if x.size == 0:
            raise ValueError(f"No sampled values for {state}")
        y = np.arange(1, x.size + 1, dtype=float) / x.size
        color, label = styles[state]

        full = positive_summary.loc[
            positive_summary["positive_state_class"].astype(str) == state
        ]
        if len(full) != 1:
            raise ValueError(f"Expected one full-data summary row for {state}")
        n_full = int(full.iloc[0]["n_clonotype_timepoints"])
        frac_full = float(full.iloc[0]["fraction_of_all_states"])
        med_full = float(full.iloc[0]["latent_minus_positive_x_median"])

        ax.plot(
            x, y,
            color=color, linewidth=1.1,
            label=label,
        )
        ax.axvline(med_full, color=color, linestyle="--", linewidth=0.7)

    ax.axvline(0.0, color="0.35", linestyle=":", linewidth=0.7)
    ax.set_xlabel(r"$x_{latent} - x_{positive}$")
    ax.set_ylabel("ECDF")
    ax.legend(frameon=False, loc="lower right", fontsize=LEGEND_FONT)
    apply_axis_config(ax, stem)
    clean_axes(ax)
    fig.tight_layout(pad=1.0)
    save_figure(fig, output_dir, stem, args)


def one_positive_relative_depth_plot(
    units: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: Path,
    args: argparse.Namespace,
    stem: str = "suppS1D_one_positive_relative_depth_coupling",
) -> None:
    require_columns(
        units,
        [
            "delta_log_depth_zero_minus_positive",
            "median_latent_minus_positive_x",
            "zero_replicate",
        ],
        "ONE_POSITIVE zero-side unit table",
    )
    require_columns(
        summary,
        [
            "scope",
            "n_analysis_units",
            "spearman_relative_log_depth_vs_unit_median_latent_minus_positive_x",
        ],
        "ONE_POSITIVE zero-side diagnostic",
    )

    d = units.copy()
    x = pd.to_numeric(d["delta_log_depth_zero_minus_positive"], errors="coerce")
    y = pd.to_numeric(d["median_latent_minus_positive_x"], errors="coerce")
    keep = np.isfinite(x) & np.isfinite(y)
    d = d.loc[keep].copy()
    x = x.loc[keep].to_numpy(float)
    y = y.loc[keep].to_numpy(float)

    all_row = summary.loc[summary["scope"].astype(str) == "ALL_ONE_POSITIVE_UNITS"]
    if len(all_row) != 1:
        raise ValueError("Expected ALL_ONE_POSITIVE_UNITS row in zero-depth diagnostic")
    rho = float(all_row.iloc[0]["spearman_relative_log_depth_vs_unit_median_latent_minus_positive_x"])
    n_units = int(all_row.iloc[0]["n_analysis_units"])

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    for side, color, label in [
        ("REP1", MAIN_BLUE, "replicate 1 is zero"),
        ("REP2", SECONDARY_ORANGE, "replicate 2 is zero"),
    ]:
        mask = d["zero_replicate"].astype(str).to_numpy() == side
        ax.scatter(
            x[mask], y[mask],
            s=15, alpha=0.80, color=color, edgecolors="none", label=label,
        )

    ax.axvline(0.0, color="0.35", linestyle="--", linewidth=0.7)
    ax.axhline(0.0, color="0.35", linestyle=":", linewidth=0.7)
    ax.set_xlabel(r"$\log(N_{zero}) - \log(N_{positive})$")
    ax.set_ylabel(r"median $(x_{latent} - x_{positive})$")
    box = get_box_config(stem)
    ax.text(
        box["x"], box["y"],
        f"n = {n_units} orientation units\nSpearman ρ = {rho:.3f}",
        transform=ax.transAxes,
        ha=box["ha"], va=box["va"], fontsize=STATS_FONT,
        bbox=dict(
            facecolor="white", edgecolor="0.82", linewidth=0.45,
            alpha=0.85, boxstyle="round,pad=0.20",
        ),
    )
    ax.legend(frameon=False, loc="lower left", fontsize=LEGEND_FONT)
    apply_axis_config(ax, stem)
    clean_axes(ax)
    fig.tight_layout(pad=1.0)
    save_figure(fig, output_dir, stem, args)


def one_positive_paired_depth_plot(
    paired: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: Path,
    args: argparse.Namespace,
    stem: str = "fig2F_one_positive_paired_depth_coupling",
) -> None:
    require_columns(
        paired,
        [
            "log_depth_ratio_rep1_over_rep2",
            "delta_shift_zero_rep1_minus_zero_rep2",
        ],
        "ONE_POSITIVE paired orientation table",
    )
    required_summary = [
        "n_complete_for_primary_relation",
        "spearman_log_depth_ratio_rep1_over_rep2_vs_delta_shift_zero_rep1_minus_zero_rep2",
        "pearson_log_depth_ratio_rep1_over_rep2_vs_delta_shift_zero_rep1_minus_zero_rep2",
        "linear_slope_log_depth_ratio_rep1_over_rep2_vs_delta_shift_zero_rep1_minus_zero_rep2",
        "linear_intercept_log_depth_ratio_rep1_over_rep2_vs_delta_shift_zero_rep1_minus_zero_rep2",
        "linear_r2_log_depth_ratio_rep1_over_rep2_vs_delta_shift_zero_rep1_minus_zero_rep2",
    ]
    require_columns(summary, required_summary, "ONE_POSITIVE paired diagnostic")
    if len(summary) != 1:
        raise ValueError("Expected exactly one row in paired diagnostic summary")

    x = pd.to_numeric(paired["log_depth_ratio_rep1_over_rep2"], errors="coerce")
    y = pd.to_numeric(paired["delta_shift_zero_rep1_minus_zero_rep2"], errors="coerce")
    keep = np.isfinite(x) & np.isfinite(y)
    x = x.loc[keep].to_numpy(float)
    y = y.loc[keep].to_numpy(float)

    row = summary.iloc[0]
    n = int(row["n_complete_for_primary_relation"])
    rho = float(row[required_summary[1]])
    pearson = float(row[required_summary[2]])
    slope = float(row[required_summary[3]])
    intercept = float(row[required_summary[4]])
    r2 = float(row[required_summary[5]])

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    ax.scatter(x, y, s=18, alpha=0.85, color=MAIN_BLUE, edgecolors="none")

    if x.size:
        grid = np.linspace(float(np.min(x)), float(np.max(x)), 100)
        ax.plot(grid, intercept + slope * grid, color="0.20", linewidth=0.9)

    ax.axvline(0.0, color="0.45", linestyle="--", linewidth=0.7)
    ax.axhline(0.0, color="0.45", linestyle=":", linewidth=0.7)
    ax.set_xlabel(r"$\log(N_1/N_2)$")
    ax.set_ylabel("Median shift difference\n(zero rep. 1 − zero rep. 2)")
    box = get_box_config(stem)
    ax.text(
        box["x"], box["y"],
        f"n = {n}\nSpearman ρ = {rho:.3f}\nPearson r = {pearson:.3f}\n$R^2$ = {r2:.3f}",
        transform=ax.transAxes,
        ha=box["ha"], va=box["va"], fontsize=STATS_FONT,
        bbox=dict(
            facecolor="white", edgecolor="0.82", linewidth=0.45,
            alpha=0.85, boxstyle="round,pad=0.20",
        ),
        zorder=ZORDER_TEXT,
    )
    apply_axis_config(ax, stem)
    clean_axes(ax)
    fig.tight_layout(pad=1.0)
    save_figure(fig, output_dir, stem, args)


def state_vs_factorized_union_plot(
    sample: pd.DataFrame,
    diagnostic: pd.DataFrame,
    output_dir: Path,
    args: argparse.Namespace,
    stem: str = "suppS1E_state_vs_factorized_union_detectability",
) -> None:
    require_columns(
        sample,
        ["cqc_factorized_union_detectability", "cqc_state_detectability"],
        "QC plotting sample",
    )
    require_columns(
        diagnostic,
        [
            "median_absolute_delta_state_vs_factorized_union",
            "pearson_state_vs_factorized_union",
            "spearman_state_vs_factorized_union",
        ],
        "factorized-union diagnostic",
    )
    if len(diagnostic) != 1:
        raise ValueError("Expected exactly one factorized-union diagnostic row")

    d = sample[
        ["cqc_factorized_union_detectability", "cqc_state_detectability"]
    ].replace([np.inf, -np.inf], np.nan).dropna()
    if d.empty:
        raise ValueError("No finite state/factorized-union detectability pairs")

    x = d["cqc_factorized_union_detectability"].to_numpy(float)
    y = d["cqc_state_detectability"].to_numpy(float)
    lo = float(min(np.quantile(x, 0.001), np.quantile(y, 0.001)))
    hi = float(max(np.quantile(x, 0.999), np.quantile(y, 0.999)))

    row = diagnostic.iloc[0]
    med_abs = float(row["median_absolute_delta_state_vs_factorized_union"])
    pearson = float(row["pearson_state_vs_factorized_union"])
    rho = float(row["spearman_state_vs_factorized_union"])

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
    hb = ax.hexbin(x, y, gridsize=60, mincnt=1, bins="log", cmap="viridis")
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=0.8, color="0.35")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Factorized-union detectability")
    ax.set_ylabel("State detectability")
    box = get_box_config(stem)
    ax.text(
        box["x"], box["y"],
        f"median |Δ| = {med_abs:.4f}\nSpearman ρ = {rho:.4f}\nPearson r = {pearson:.4f}",
        transform=ax.transAxes,
        ha=box["ha"], va=box["va"], fontsize=STATS_FONT,
        bbox=dict(
            facecolor="white", edgecolor="0.82", linewidth=0.45,
            alpha=0.85, boxstyle="round,pad=0.20",
        ),
    )
    cbar = fig.colorbar(hb, ax=ax)
    cbar.set_label("Count (log scale)")
    apply_axis_config(ax, stem)
    clean_axes(ax)
    fig.tight_layout(pad=1.0)
    save_figure(fig, output_dir, stem, args)


def correlation_heatmap(
    corr: pd.DataFrame,
    output_dir: Path,
    args: argparse.Namespace,
    stem: str = "suppS1F_qc_metric_correlation",
) -> None:
    labels = {
        "cqc_posterior_width95_x": "Width",
        "cqc_posterior_entropy": "Entropy",
        "cqc_state_detectability": "State\ndetectability",
        "cqc_factorized_union_detectability": "Factorized\nunion",
        "cqc_replicate_bhattacharyya": "Replicate\noverlap",
        "cqc_effective_posterior_support_simpson": "Inv. Simpson",
    }
    keep = [c for c in labels if c in corr.index and c in corr.columns]
    if len(keep) < 2:
        raise ValueError("Too few canonical metrics available for correlation heatmap")

    d = corr.loc[keep, keep]
    values = d.to_numpy(float)
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

    for i in range(len(d)):
        for j in range(len(d)):
            if i > j and np.isfinite(d.iloc[i, j]):
                value = float(d.iloc[i, j])
                text_color = "white" if abs(value) >= 0.55 else "black"
                ax.text(
                    j, i, f"{value:.2f}",
                    ha="center", va="center",
                    fontsize=ANNOT_FONT, color=text_color,
                )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Spearman correlation")
    apply_axis_config(ax, stem)
    fig.tight_layout(pad=1.0)
    save_figure(fig, output_dir, stem, args)


# -----------------------------------------------------------------------------
# CLI and main
# -----------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Plot ClonoDynamics Step-3 latent-state and observation-model "
            "diagnostics using canonical replicate-resolved outputs."
        )
    )
    p.add_argument("--analysis-dir", type=Path, required=True)
    p.add_argument(
        "--figures-dir",
        type=Path,
        default=None,
        help=(
            "Default: <analysis-dir>/"
            "figures_latent_state_and_observation_model_diagnostics"
        ),
    )
    p.add_argument("--sample-n", type=int, default=250000)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--font", default="Arial")
    p.add_argument("--font-size", type=float, default=8.0)
    p.add_argument(
        "--show-bin-counts",
        action="store_true",
        help="Add a secondary log-scale count curve to abundance panels.",
    )
    p.add_argument(
        "--pdf", action=argparse.BooleanOptionalAction, default=True
    )
    p.add_argument(
        "--png", action=argparse.BooleanOptionalAction, default=False
    )
    p.add_argument(
        "--svg", action=argparse.BooleanOptionalAction, default=False
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    analysis_dir = args.analysis_dir.expanduser().resolve()
    if not analysis_dir.exists():
        raise FileNotFoundError(f"Analysis directory not found: {analysis_dir}")

    root = ensure_dir(
        args.figures_dir.expanduser().resolve()
        if args.figures_dir is not None
        else analysis_dir / "figures_latent_state_and_observation_model_diagnostics"
    )
    main_dir = ensure_dir(root / "main")
    supplementary_dir = ensure_dir(root / "supplementary")
    clear_plot_outputs(main_dir)
    clear_plot_outputs(supplementary_dir)
    set_rcparams(args.font, args.font_size, args.dpi)

    paths = {
        "parquet": require_file(analysis_dir / "clonotype_latent_posterior_qc.parquet"),
        "metric_summary": require_file(
            analysis_dir / "clonotype_latent_posterior_qc_metric_summary.csv"
        ),
        "abundance": require_file(
            analysis_dir / "clonotype_latent_posterior_qc_by_abundance.csv"
        ),
        "positive_summary": require_file(
            analysis_dir / "clonotype_latent_posterior_qc_by_positive_state.csv"
        ),
        "zero_units": require_file(
            analysis_dir / "one_positive_zero_depth_by_subject_time_replica.csv"
        ),
        "zero_diagnostic": require_file(
            analysis_dir / "one_positive_zero_depth_diagnostic.csv"
        ),
        "paired_units": require_file(
            analysis_dir / "one_positive_paired_replica_depth_by_subject_time.csv"
        ),
        "paired_diagnostic": require_file(
            analysis_dir / "one_positive_paired_replica_depth_diagnostic.csv"
        ),
        "factorized_diagnostic": require_file(
            analysis_dir / "state_detectability_factorized_union_diagnostic.csv"
        ),
        "spearman_matrix": require_file(
            analysis_dir / "qc_metric_spearman_correlations.csv"
        ),
    }

    metric_summaries = load_metric_summaries(paths["metric_summary"])
    abundance = pd.read_csv(paths["abundance"])
    positive_summary = pd.read_csv(paths["positive_summary"])
    zero_units = pd.read_csv(paths["zero_units"])
    zero_diagnostic = pd.read_csv(paths["zero_diagnostic"])
    paired_units = pd.read_csv(paths["paired_units"])
    paired_diagnostic = pd.read_csv(paths["paired_diagnostic"])
    factorized_diagnostic = pd.read_csv(paths["factorized_diagnostic"])
    full_spearman = load_spearman_matrix(paths["spearman_matrix"])

    canonical_sample_cols = [
        "cqc_posterior_width95_x",
        "cqc_replicate_bhattacharyya",
        "positive_state_class",
        "cqc_latent_minus_positive_x",
        "cqc_factorized_union_detectability",
        "cqc_state_detectability",
    ]
    sample = systematic_sample(
        paths["parquet"],
        canonical_sample_cols,
        sample_n=int(args.sample_n),
        seed=int(args.seed),
    )

    width_summary = metric_summaries.get("cqc_posterior_width95_x")
    bc_summary = metric_summaries.get("cqc_replicate_bhattacharyya")
    if width_summary is None or bc_summary is None:
        raise ValueError(
            "Metric summary is missing posterior-width or exact-overlap summaries."
        )

    median_width = float(width_summary["median"])
    median_bc = float(bc_summary["median"])

    # ------------------------------------------------------------------
    # Main manuscript candidates: D-F
    # ------------------------------------------------------------------
    interval_curve(
        abundance,
        "x_bin_left", "x_bin_right",
        "posterior_width95_median",
        "posterior_width95_q25",
        "posterior_width95_q75",
        "Latent log-frequency bin midpoint",
        "Posterior 95% width (x)",
        "fig2D_abundance_vs_posterior_width",
        main_dir,
        args,
        annotation_text=f"overall median = {median_width:.2f}",
    )

    interval_curve(
        abundance,
        "x_bin_left", "x_bin_right",
        "replicate_bc_median",
        "replicate_bc_q25",
        "replicate_bc_q75",
        "Latent log-frequency bin midpoint",
        "Replicate posterior overlap (BC)",
        "fig2E_abundance_vs_replicate_overlap",
        main_dir,
        args,
        annotation_text=f"overall median BC = {median_bc:.3f}",
    )

    one_positive_paired_depth_plot(
        paired_units,
        paired_diagnostic,
        main_dir,
        args,
        "fig2F_one_positive_paired_depth_coupling",
    )

    # ------------------------------------------------------------------
    # Supplementary diagnostics
    # ------------------------------------------------------------------
    histogram_panel(
        sample["cqc_posterior_width95_x"],
        "Posterior 95% width (x)",
        "suppS1A_posterior_width_distribution",
        supplementary_dir,
        args,
        width_summary,
    )

    histogram_panel(
        sample["cqc_replicate_bhattacharyya"],
        "Replicate posterior overlap (Bhattacharyya coefficient)",
        "suppS1B_replicate_overlap_distribution",
        supplementary_dir,
        args,
        bc_summary,
    )

    ecdf_by_positive_state(
        sample,
        positive_summary,
        supplementary_dir,
        args,
        "suppS1C_latent_vs_positive_by_detection_support",
    )

    one_positive_relative_depth_plot(
        zero_units,
        zero_diagnostic,
        supplementary_dir,
        args,
        "suppS1D_one_positive_relative_depth_coupling",
    )

    state_vs_factorized_union_plot(
        sample,
        factorized_diagnostic,
        supplementary_dir,
        args,
        "suppS1E_state_vs_factorized_union_detectability",
    )

    correlation_heatmap(
        full_spearman,
        supplementary_dir,
        args,
        "suppS1F_qc_metric_correlation",
    )

    # Key run summary from full-data tables.
    one = positive_summary.loc[
        positive_summary["positive_state_class"].astype(str) == "ONE_POSITIVE"
    ].iloc[0]
    paired = paired_diagnostic.iloc[0]
    zero_all = zero_diagnostic.loc[
        zero_diagnostic["scope"].astype(str) == "ALL_ONE_POSITIVE_UNITS"
    ].iloc[0]

    print("[STEP3-PLOT] Completed canonical single-panel output set.")
    print(f"[STEP3-PLOT] version={SCRIPT_VERSION}")
    print(f"[STEP3-PLOT] analysis={analysis_dir}")
    print(f"[STEP3-PLOT] main={main_dir}")
    print(f"[STEP3-PLOT] supplementary={supplementary_dir}")
    print(
        "[STEP3-PLOT] ONE_POSITIVE: n={:,}; fraction={:.4f}; "
        "median latent-positive shift={:.4f}".format(
            int(one["n_clonotype_timepoints"]),
            float(one["fraction_of_all_states"]),
            float(one["latent_minus_positive_x_median"]),
        )
    )
    print(
        "[STEP3-PLOT] relative-depth diagnostic: n_units={}; Spearman={:.6f}".format(
            int(zero_all["n_analysis_units"]),
            float(
                zero_all[
                    "spearman_relative_log_depth_vs_unit_median_latent_minus_positive_x"
                ]
            ),
        )
    )
    print(
        "[STEP3-PLOT] paired depth diagnostic: n={}; Spearman={:.6f}; "
        "Pearson={:.6f}; R2={:.6f}".format(
            int(paired["n_complete_for_primary_relation"]),
            float(
                paired[
                    "spearman_log_depth_ratio_rep1_over_rep2_vs_delta_shift_zero_rep1_minus_zero_rep2"
                ]
            ),
            float(
                paired[
                    "pearson_log_depth_ratio_rep1_over_rep2_vs_delta_shift_zero_rep1_minus_zero_rep2"
                ]
            ),
            float(
                paired[
                    "linear_r2_log_depth_ratio_rep1_over_rep2_vs_delta_shift_zero_rep1_minus_zero_rep2"
                ]
            ),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
