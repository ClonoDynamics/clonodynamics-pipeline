#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_step13_temporal_scaling.py
===============================

Publication plotter for final ClonoDynamics Step 13 temporal-scaling analysis.

Final contract:
- Step 13 v3 fixed primary `cross_cov` cohort/core;
- complete-case cohort defined by `cross_cov` only;
- technical comparators cannot redefine subject composition;
- positive-incremental model constrains D >= 0 only, with K unconstrained.

All publication aesthetics remain centralized in user-editable dictionaries;
Arial is the default font and all text defaults to 8 pt.

The script writes SINGLE-PANEL PDF and PNG files only.

Recommended main Figure 6
-------------------------
Figure6A_primary_crosscov_vs_lag
    Primary equal-bin / equal-subject cross-replicate covariance across
    temporal lags 1–5, with biological subject-bootstrap intervals.

Figure6B_temporal_decomposition
    Primary temporal curves for cross_cov, same_var_mean and
    replicate_specific_excess.

Figure6C_weighting_sensitivity
    Primary equal-bin estimator versus df-weighted-within-subject sensitivity
    for cross_cov.

Figure6D_subject_slopes
    Subject-specific cross_cov temporal slopes, with the cohort slope and its
    subject-bootstrap interval.

Figure6E_binwise_slopes
    Abundance-resolved temporal slopes across the fixed Step-13 abundance core.

Recommended Supplementary Figure 7
----------------------------------
SuppFigS7A_core_selection
    Lag-specific subject coverage across Step-11 abundance bins, with the
    selected fixed core highlighted.

SuppFigS7B_bootstrap_slope_distribution
    Joint subject-bootstrap distribution of the primary cross_cov temporal
    slope.

SuppFigS7C_model_comparison
    Observed AICc values for constant, signed-linear and positive-incremental
    models for the primary cross_cov temporal curve.

SuppFigS7D_subject_temporal_curves
    Subject-level primary cross_cov temporal trajectories within the fixed
    abundance core.

Typical run
-----------
python3 plot_step13_temporal_scaling.py \
    --step13-dir ./dataset_longitudinal_results/13-temporal_fluctuation_scaling \
    --outdir ./figures/step13_temporal_scaling

All dimensions, line widths, colors, limits, legend placement and annotation
positions are editable in the dictionaries below.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt


# =============================================================================
# USER-EDITABLE PUBLICATION AESTHETICS
# =============================================================================
#
# All dimensions below are intended to be edited manually from this block.
# Figure sizes are in inches. Font sizes and line widths are in points.
#

FONT = {
    "family": "Arial",
    "default_pt": 8.0,
    "axis_label_pt": 8.0,
    "tick_label_pt": 7.0,
    "legend_pt": 7.0,
    "annotation_pt": 7.0,
    "title_pt": 7.0,
}

EXPORT = {
    "pdf_font_type": 42,
    "ps_font_type": 42,
    "png_dpi": 600,
    "save_formats": ("pdf", "png"),
    "bbox_inches": "tight",
    "pad_inches": 0.06,
}

# Individual panel dimensions in inches.
FIGSIZE_INCHES: Dict[str, Tuple[float, float]] = {
    "Figure6A_primary_crosscov_vs_lag": (3.2, 2.4),
    "Figure6B_temporal_decomposition": (3.2, 2.4),
    "Figure6C_weighting_sensitivity": (3.2, 2.4),
    "Figure6D_subject_slopes": (3.2, 2.4),
    "Figure6E_binwise_slopes": (3.2, 2.4),
    "SuppFigS7A_core_selection": (3.2, 2.4),
    "SuppFigS7B_bootstrap_slope_distribution": (3.2, 2.4),
    "SuppFigS7C_model_comparison": (3.2, 2.4),
    "SuppFigS7D_subject_temporal_curves": (3.2, 2.4),
}

LINE_WIDTHS = {
    "standard": 1.20,
    "primary": 1.50,
    "secondary": 1.00,
    "errorbar": 1.00,
    "zero": 0.80,
    "axis": 0.80,
    "tick": 0.75,
    "ci_boundary": 0.80,
}

MARKER_SIZES = {
    "standard": 3.5,
    "point": 4.5,
    "subject": 3.0,
    "scatter_area": 22.0,   # matplotlib scatter uses points^2
    "bar_width": 0.62,
    "errorbar_capsize": 2.5,
}

ALPHAS = {
    "ribbon": 0.15,
    "cohort_ci": 0.12,
    "subject_line": 0.42,
    "histogram": 0.78,
    "zero_line": 0.65,
    "core_span": 0.08,
    "bootstrap_ci_span": 0.07,
    "bar": 0.88,
}

AXIS_STYLE = {
    "tick_direction": "out",
    "major_tick_length": 3.5,
    "title_pad": 6.0,
    "grid": False,
}

COLORS = {
    "cross": "#0072B2",
    "same": "#D55E00",
    "excess": "#CC79A7",
    "primary": "#0072B2",
    "sensitivity": "#009E73",
    "subject": "#4D4D4D",
    "cohort": "#222222",
    "bootstrap": "#7A7A7A",
    "core": "#0072B2",
    "noncore": "#BDBDBD",
    "zero": "#666666",
    "constant": "#7A7A7A",
    "signed_linear": "#0072B2",
    "positive_incremental": "#D55E00",
}

LINESTYLES = {
    "cross": "-",
    "same": "-",
    "excess": "--",
    "primary": "-",
    "sensitivity": "--",
    "subject": "-",
    "zero": "--",
}

MARKERS = {
    "cross": "o",
    "same": "s",
    "excess": "^",
    "primary": "o",
    "sensitivity": "s",
    "subject": "o",
}

LEGEND_SPECS = {
    "Figure6B_temporal_decomposition": {
        "loc": "lower center",
        "bbox_to_anchor": (0.50, 1.02),
        "ncol": 2,
        "fontsize": 8.0,
    },
    "Figure6C_weighting_sensitivity": {
        "loc": "lower center",
        "bbox_to_anchor": (0.50, 1.02),
        "ncol": 2,
        "fontsize": 8.0,
    },
    "Figure6D_subject_slopes": {
        "loc": "lower center",
        "bbox_to_anchor": (0.50, 1.02),
        "ncol": 2,
        "fontsize": 8.0,
    },
    "SuppFigS7A_core_selection": {
        "loc": "lower center",
        "bbox_to_anchor": (0.50, 1.02),
        "ncol": 3,
        "fontsize": 8.0,
    },
    "SuppFigS7D_subject_temporal_curves": {
        "loc": "lower center",
        "bbox_to_anchor": (0.50, 1.02),
        "ncol": 2,
        "fontsize": 8.0,
    },
}

SHOW_TITLES = False

PLOT_TITLES = {
    "Figure6A_primary_crosscov_vs_lag":
        "Temporal scaling of shared fluctuation covariance",
    "Figure6B_temporal_decomposition":
        "Temporal fluctuation decomposition",
    "Figure6C_weighting_sensitivity":
        "Weighting sensitivity of temporal scaling",
    "Figure6D_subject_slopes":
        "Subject-specific temporal slopes",
    "Figure6E_binwise_slopes":
        "Abundance-resolved temporal slopes",
    "SuppFigS7A_core_selection":
        "Fixed abundance-core selection",
    "SuppFigS7B_bootstrap_slope_distribution":
        "Bootstrap slope distribution",
    "SuppFigS7C_model_comparison":
        "Finite-lag model comparison",
    "SuppFigS7D_subject_temporal_curves":
        "Subject-level temporal trajectories",
}

XLIMS = {
    "Figure6A_primary_crosscov_vs_lag": (0.7, 5.3),
    "Figure6B_temporal_decomposition": (0.7, 5.3),
    "Figure6C_weighting_sensitivity": (0.7, 5.3),
    "Figure6D_subject_slopes": None,
    "Figure6E_binwise_slopes": None,
    "SuppFigS7A_core_selection": None,
    "SuppFigS7B_bootstrap_slope_distribution": (-0.1, 0.1),
    "SuppFigS7C_model_comparison": None,
    "SuppFigS7D_subject_temporal_curves": (0.7, 5.3),
}

YLIMS = {
    "Figure6A_primary_crosscov_vs_lag": (0.3, 0.7),
    "Figure6B_temporal_decomposition": None,
    "Figure6C_weighting_sensitivity": None,
    "Figure6D_subject_slopes": (-0.15, 0.1),
    "Figure6E_binwise_slopes": (-0.25, 0.1),
    "SuppFigS7A_core_selection": (0.0, 1.05),
    "SuppFigS7B_bootstrap_slope_distribution": None,
    "SuppFigS7C_model_comparison": None,
    "SuppFigS7D_subject_temporal_curves": (0,1),
}

SHOW_ZERO_LINE = True
SHOW_TEMPORAL_LINEAR_DESCRIPTOR = True

HISTOGRAM = {
    "bins": 35,
}

ANNOTATIONS = {
    "Figure6A_primary_crosscov_vs_lag": {
        "enabled": True,
        "xy_axes": (0.03, 0.035),
        "ha": "left",
        "va": "bottom",
        "fontsize": 8.0,
        "format": (
            "slope = {slope:.4f} week$^{{-1}}$\n"
            "95% bootstrap CI [{lo:.4f}, {hi:.4f}]\n"
            "bootstrap P(slope > 0) = {ppos:.4f}"
        ),
    },
    "Figure6D_subject_slopes": {
        "enabled": True,
        "xy_axes": (0.97, 0.95),
        "ha": "right",
        "va": "top",
        "fontsize": 8.0,
        "format": (
            "{nneg}/{n} subject slopes < 0\n"
            "exact sign-flip P = {p:.4f}"
        ),
    },
    "Figure6E_binwise_slopes": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "ha": "left",
        "va": "bottom",
        "fontsize": 8.0,
        "format": (
            "{nneg}/{n} observed bin slopes < 0\n"
            "{nci}/{n} bootstrap CIs entirely < 0"
        ),
    },
    "SuppFigS7B_bootstrap_slope_distribution": {
        "enabled": True,
        "xy_axes": (0.53, 0.95),
        "ha": "left",
        "va": "top",
        "fontsize": 8.0,
        "format": (
            "slope = {slope:.4f}\n"
            "95% CI [{lo:.4f}, {hi:.4f}]"
        ),
    },
}


# =============================================================================
# Matplotlib helpers
# =============================================================================

def configure_matplotlib() -> None:
    mpl.rcParams.update({
        "font.family": FONT["family"],
        "font.size": FONT["default_pt"],
        "axes.labelsize": FONT["axis_label_pt"],
        "axes.titlesize": FONT["title_pt"],
        "xtick.labelsize": FONT["tick_label_pt"],
        "ytick.labelsize": FONT["tick_label_pt"],
        "legend.fontsize": FONT["legend_pt"],
        "pdf.fonttype": EXPORT["pdf_font_type"],
        "ps.fonttype": EXPORT["ps_font_type"],
        "axes.linewidth": LINE_WIDTHS["axis"],
        "xtick.major.width": LINE_WIDTHS["tick"],
        "ytick.major.width": LINE_WIDTHS["tick"],
        "xtick.major.size": AXIS_STYLE["major_tick_length"],
        "ytick.major.size": AXIS_STYLE["major_tick_length"],
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def make_figure(panel: str) -> Tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(
        figsize=FIGSIZE_INCHES[panel],
        constrained_layout=False,
    )
    style_axis(ax)
    return fig, ax


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction=AXIS_STYLE["tick_direction"])
    ax.grid(AXIS_STYLE["grid"])


def maybe_zero_line(ax: plt.Axes, orientation: str = "horizontal") -> None:
    if not SHOW_ZERO_LINE:
        return

    kwargs = dict(
        color=COLORS["zero"],
        lw=LINE_WIDTHS["zero"],
        ls=LINESTYLES["zero"],
        alpha=ALPHAS["zero_line"],
        zorder=0,
    )

    if orientation == "horizontal":
        ax.axhline(0.0, **kwargs)
    else:
        ax.axvline(0.0, **kwargs)


def maybe_title(ax: plt.Axes, panel: str) -> None:
    if SHOW_TITLES:
        ax.set_title(
            PLOT_TITLES.get(panel, panel),
            pad=AXIS_STYLE["title_pad"],
        )


def apply_limits(ax: plt.Axes, panel: str) -> None:
    xlim = XLIMS.get(panel)
    ylim = YLIMS.get(panel)

    if xlim is not None:
        current = ax.get_xlim()
        lo, hi = xlim
        ax.set_xlim(
            current[0] if lo is None else lo,
            current[1] if hi is None else hi,
        )

    if ylim is not None:
        current = ax.get_ylim()
        lo, hi = ylim
        ax.set_ylim(
            current[0] if lo is None else lo,
            current[1] if hi is None else hi,
        )


def add_annotation(ax: plt.Axes, panel: str, text: str) -> None:
    spec = ANNOTATIONS.get(panel, {})
    if not spec.get("enabled", False):
        return

    x, y = spec.get("xy_axes", (0.03, 0.04))
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha=spec.get("ha", "left"),
        va=spec.get("va", "bottom"),
        fontsize=spec.get("fontsize", FONT["annotation_pt"]),
    )


def add_panel_legend(ax: plt.Axes, panel: str) -> None:
    spec = LEGEND_SPECS.get(panel)
    if spec is None:
        ax.legend(frameon=False, fontsize=FONT["legend_pt"])
        return

    ax.legend(
        frameon=False,
        loc=spec.get("loc", "best"),
        bbox_to_anchor=spec.get("bbox_to_anchor"),
        ncol=int(spec.get("ncol", 1)),
        fontsize=spec.get("fontsize", FONT["legend_pt"]),
        columnspacing=1.2,
        handletextpad=0.6,
        borderaxespad=0.0,
    )


def save_figure(fig: plt.Figure, outdir: Path, panel: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    for fmt in EXPORT["save_formats"]:
        path = outdir / f"{panel}.{fmt}"
        kwargs = dict(
            bbox_inches=EXPORT["bbox_inches"],
            pad_inches=EXPORT["pad_inches"],
        )
        if fmt.lower() == "png":
            kwargs["dpi"] = EXPORT["png_dpi"]
        fig.savefig(path, **kwargs)

    plt.close(fig)


# =============================================================================
# IO
# =============================================================================

def read_json_required(
    path: Path,
    label: str,
) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def bool_series(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "t", "yes", "y"])
    )


def validate_final_contract(
    step13_dir: Path,
) -> Dict[str, object]:
    config = read_json_required(
        step13_dir / "00_run_config.json",
        "Step 13 run config",
    )
    signature = read_json_required(
        step13_dir / "00_run_signature.json",
        "Step 13 run signature",
    )

    if config.get("primary_metric") != "cross_cov":
        raise ValueError(
            "Final plotter requires Step 13 primary_metric='cross_cov'."
        )
    if config.get("primary_estimator") != "equal_bin_equal_subject":
        raise ValueError(
            "Final plotter requires equal_bin_equal_subject primary estimator."
        )

    complete = config.get("complete_case", {})
    if complete.get("defined_by_metric") != "cross_cov":
        raise ValueError(
            "Final plotter requires complete-case cohort defined by cross_cov."
        )
    if complete.get("technical_comparators_cannot_redefine_cohort") is not True:
        raise ValueError(
            "Final plotter requires frozen primary cohort for comparators."
        )

    models = config.get("model_comparison", {})
    if models.get("positive_incremental_constraint") != (
        "slope D >= 0 only; intercept K unconstrained"
    ):
        raise ValueError(
            "Final plotter requires corrected positive-incremental constraint."
        )

    payload = signature.get("signature_payload", {})
    if payload.get("primary_complete_case_defined_by") != "cross_cov":
        raise ValueError(
            "Step 13 signature does not certify cross_cov-defined cohort."
        )

    selected_bins = [
        int(x)
        for x in config.get("core_selection", {}).get("selected_bins", [])
    ]
    subjects = [
        int(x)
        for x in complete.get("subjects", [])
    ]

    if not selected_bins or not subjects:
        raise ValueError("Step 13 config has empty primary core or cohort.")

    return {
        "step13_script_version": config.get("script_version"),
        "step13_analysis_signature": config.get("analysis_signature"),
        "core_min_subject_fraction": float(
            config.get("core_selection", {}).get(
                "core_min_subject_fraction", np.nan
            )
        ),
        "selected_core_bins": ";".join(map(str, selected_bins)),
        "n_core_bins": int(len(selected_bins)),
        "complete_case_subjects": ";".join(map(str, subjects)),
        "n_complete_case_subjects": int(len(subjects)),
        "primary_metric": config.get("primary_metric"),
        "primary_estimator": config.get("primary_estimator"),
        "sensitivity_estimator": config.get("sensitivity_estimator"),
    }


def read_csv_required(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")

    d = pd.read_csv(path)
    if d.empty:
        raise ValueError(f"{label} is empty: {path}")
    return d


def require_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
    label: str,
) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(
            f"{label}: missing required columns {missing}; "
            f"available={list(frame.columns)}"
        )


def write_source_data(
    frame: pd.DataFrame,
    source_dir: Path,
    panel: str,
) -> None:
    source_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    frame.to_csv(
        source_dir / f"{panel}_source_data.csv",
        index=False,
    )


def load_inputs(step13_dir: Path) -> Dict[str, object]:
    contract = validate_final_contract(step13_dir)
    return {
        "contract": contract,
        "summary": read_csv_required(
            step13_dir / "00_analysis_summary.csv",
            "Step 13 analysis summary",
        ),
        "core": read_csv_required(
            step13_dir / "01_core_bin_selection.csv",
            "Step 13 core-bin selection",
        ),
        "complete": read_csv_required(
            step13_dir / "02_complete_case_subjects.csv",
            "Step 13 complete-case subjects",
        ),
        "subject_core": read_csv_required(
            step13_dir / "03_subject_core_metric_by_dt.csv",
            "Step 13 subject core metrics",
        ),
        "cohort": read_csv_required(
            step13_dir / "04_cohort_core_metric_by_dt.csv",
            "Step 13 cohort core metrics",
        ),
        "slope": read_csv_required(
            step13_dir / "05_temporal_slope_summary.csv",
            "Step 13 temporal slope summary",
        ),
        "subject_slopes": read_csv_required(
            step13_dir / "06_subject_slopes.csv",
            "Step 13 subject slopes",
        ),
        "signflip": read_csv_required(
            step13_dir / "07_exact_signflip_tests.csv",
            "Step 13 sign-flip tests",
        ),
        "binwise": read_csv_required(
            step13_dir / "08_binwise_temporal_slopes.csv",
            "Step 13 binwise temporal slopes",
        ),
        "models": read_csv_required(
            step13_dir / "09_model_comparison.csv",
            "Step 13 model comparison",
        ),
        "bootstrap": read_csv_required(
            step13_dir / "10_joint_subject_bootstrap.csv",
            "Step 13 joint subject bootstrap",
        ),
        "step14_contract": read_csv_required(
            step13_dir / "11_step14_contract.csv",
            "Step 13 Step-14 contract",
        ),
        "support_audit": read_csv_required(
            step13_dir / "12_core_cell_support_audit.csv",
            "Step 13 core-cell support audit",
        ),
    }


# =============================================================================
# Data selectors
# =============================================================================

def primary_cohort(cohort: pd.DataFrame, metric: str) -> pd.DataFrame:
    d = cohort[
        cohort["estimator"].astype(str).eq("equal_bin_equal_subject")
        & cohort["metric"].astype(str).eq(metric)
    ].copy().sort_values("dt")

    if d.empty:
        raise ValueError(f"No primary cohort rows for metric={metric}")

    return d


def sensitivity_cohort(cohort: pd.DataFrame, metric: str) -> pd.DataFrame:
    d = cohort[
        cohort["estimator"].astype(str).eq(
            "df_weighted_within_subject_equal_subject"
        )
        & cohort["metric"].astype(str).eq(metric)
    ].copy().sort_values("dt")

    if d.empty:
        raise ValueError(f"No sensitivity cohort rows for metric={metric}")

    return d


def primary_slope_row(slope: pd.DataFrame, metric: str) -> pd.Series:
    d = slope[
        slope["estimator"].astype(str).eq("equal_bin_equal_subject")
        & slope["metric"].astype(str).eq(metric)
    ]

    if d.empty:
        raise ValueError(f"No primary slope row for metric={metric}")

    return d.iloc[0]


# =============================================================================
# Figure 6A
# =============================================================================

def plot_figure6A(
    cohort: pd.DataFrame,
    slope: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure6A_primary_crosscov_vs_lag"

    d = primary_cohort(
        cohort,
        "cross_cov",
    )
    write_source_data(d, source_dir, panel)

    row = primary_slope_row(
        slope,
        "cross_cov",
    )

    fig, ax = make_figure(panel)

    x = pd.to_numeric(
        d["dt"], errors="coerce"
    ).to_numpy(float)
    y = pd.to_numeric(
        d["value"], errors="coerce"
    ).to_numpy(float)
    lo = pd.to_numeric(
        d["bootstrap_q025"], errors="coerce"
    ).to_numpy(float)
    hi = pd.to_numeric(
        d["bootstrap_q975"], errors="coerce"
    ).to_numpy(float)

    ax.errorbar(
        x,
        y,
        yerr=np.vstack([y - lo, hi - y]),
        fmt=MARKERS["primary"],
        color=COLORS["primary"],
        lw=LINE_WIDTHS["errorbar"],
        ms=MARKER_SIZES["point"],
        capsize=MARKER_SIZES["errorbar_capsize"],
        label="Cross-replicate covariance",
        zorder=3,
    )

    if SHOW_TEMPORAL_LINEAR_DESCRIPTOR:
        intercept = float(
            row["observed_intercept_at_dt1"]
        )
        slope_value = float(
            row["observed_slope_per_week"]
        )
        xline = np.linspace(
            np.min(x),
            np.max(x),
            200,
        )
        yline = (
            intercept
            + slope_value * (xline - 1.0)
        )
        ax.plot(
            xline,
            yline,
            color=COLORS["primary"],
            lw=LINE_WIDTHS["secondary"],
            ls="--",
            zorder=2,
        )

    # Covariance is far above zero here; forcing the zero baseline into
    # the axis wastes vertical resolution and visually suppresses the lag effect.
    ax.set_xlabel("Temporal lag (weeks)")
    ax.set_ylabel(
        "Core-averaged cross covariance"
    )
    ax.set_xticks([1, 2, 3, 4, 5])

    spec = ANNOTATIONS[panel]
    text = spec["format"].format(
        slope=float(
            row["observed_slope_per_week"]
        ),
        lo=float(
            row["bootstrap_q025"]
        ),
        hi=float(
            row["bootstrap_q975"]
        ),
        ppos=float(
            row[
                "bootstrap_positive_slope_fraction"
            ]
        ),
    )
    add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 6B
# =============================================================================

def plot_figure6B(
    cohort: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure6B_temporal_decomposition"

    rows = []
    for metric in [
        "cross_cov",
        "same_var_mean",
        "replicate_specific_excess",
    ]:
        d = primary_cohort(cohort, metric)
        rows.append(d)

    source = pd.concat(
        rows,
        ignore_index=True,
    )
    write_source_data(source, source_dir, panel)

    fig, ax = make_figure(panel)

    specs = [
        (
            "cross_cov",
            "Cross-replicate covariance",
            COLORS["cross"],
            LINESTYLES["cross"],
            MARKERS["cross"],
        ),
        (
            "same_var_mean",
            "Within-replicate variance",
            COLORS["same"],
            LINESTYLES["same"],
            MARKERS["same"],
        ),
        (
            "replicate_specific_excess",
            "Replicate-specific excess",
            COLORS["excess"],
            LINESTYLES["excess"],
            MARKERS["excess"],
        ),
    ]

    for metric, label, color, ls, marker in specs:
        d = primary_cohort(cohort, metric)

        x = d["dt"].to_numpy(float)
        y = d["value"].to_numpy(float)
        lo = d["bootstrap_q025"].to_numpy(float)
        hi = d["bootstrap_q975"].to_numpy(float)

        ax.errorbar(
            x,
            y,
            yerr=np.vstack([y - lo, hi - y]),
            fmt=marker,
            color=color,
            lw=LINE_WIDTHS["errorbar"],
            ms=MARKER_SIZES["standard"],
            capsize=MARKER_SIZES["errorbar_capsize"],
            label=label,
            zorder=3,
        )

        ax.plot(
            x,
            y,
            color=color,
            lw=LINE_WIDTHS["standard"],
            ls=ls,
            zorder=2,
        )

    ax.set_xlabel("Temporal lag (weeks)")
    ax.set_ylabel(
        "Core-averaged variance / covariance"
    )
    ax.set_xticks([1, 2, 3, 4, 5])

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    add_panel_legend(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 6C
# =============================================================================

def plot_figure6C(
    cohort: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure6C_weighting_sensitivity"

    primary = primary_cohort(
        cohort,
        "cross_cov",
    )
    sensitivity = sensitivity_cohort(
        cohort,
        "cross_cov",
    )

    source = pd.concat(
        [primary, sensitivity],
        ignore_index=True,
    )
    write_source_data(source, source_dir, panel)

    fig, ax = make_figure(panel)

    specs = [
        (
            primary,
            "Equal-bin / equal-subject",
            COLORS["primary"],
            LINESTYLES["primary"],
            MARKERS["primary"],
        ),
        (
            sensitivity,
            "DF-weighted within subject",
            COLORS["sensitivity"],
            LINESTYLES["sensitivity"],
            MARKERS["sensitivity"],
        ),
    ]

    for d, label, color, ls, marker in specs:
        x = d["dt"].to_numpy(float)
        y = d["value"].to_numpy(float)
        lo = d["bootstrap_q025"].to_numpy(float)
        hi = d["bootstrap_q975"].to_numpy(float)

        ax.errorbar(
            x,
            y,
            yerr=np.vstack([y - lo, hi - y]),
            fmt=marker,
            color=color,
            lw=LINE_WIDTHS["errorbar"],
            ms=MARKER_SIZES["standard"],
            capsize=MARKER_SIZES["errorbar_capsize"],
            label=label,
        )
        ax.plot(
            x,
            y,
            color=color,
            lw=LINE_WIDTHS["standard"],
            ls=ls,
        )

    ax.set_xlabel("Temporal lag (weeks)")
    ax.set_ylabel(
        "Core-averaged cross-replicate covariance"
    )
    ax.set_xticks([1, 2, 3, 4, 5])

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    add_panel_legend(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 6D
# =============================================================================

def plot_figure6D(
    subject_slopes: pd.DataFrame,
    slope: pd.DataFrame,
    signflip: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure6D_subject_slopes"

    d = subject_slopes[
        subject_slopes["estimator"].astype(str).eq(
            "equal_bin_equal_subject"
        )
        & subject_slopes["metric"].astype(str).eq(
            "cross_cov"
        )
    ].copy().sort_values("subject")

    write_source_data(d, source_dir, panel)

    slope_row = primary_slope_row(
        slope,
        "cross_cov",
    )

    test = signflip[
        signflip["estimator"].astype(str).eq(
            "equal_bin_equal_subject"
        )
        & signflip["metric"].astype(str).eq(
            "cross_cov"
        )
    ]
    if test.empty:
        raise ValueError(
            "Primary cross_cov sign-flip row missing."
        )
    test = test.iloc[0]

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    subjects = [f"P{int(x):02d}" for x in d["subject"].tolist()]
    x = np.arange(len(subjects))
    y = d["slope_per_week"].to_numpy(float)

    ax.scatter(
        x,
        y,
        s=MARKER_SIZES["scatter_area"],
        color=COLORS["subject"],
        zorder=3,
    )

    cohort_slope = float(
        slope_row["observed_slope_per_week"]
    )
    lo = float(
        slope_row["bootstrap_q025"]
    )
    hi = float(
        slope_row["bootstrap_q975"]
    )

    ax.axhspan(
        lo,
        hi,
        color=COLORS["cohort"],
        alpha=ALPHAS["cohort_ci"],
        linewidth=0,
        label="Cohort 95% bootstrap CI",
        zorder=1,
    )
    ax.axhline(
        cohort_slope,
        color=COLORS["cohort"],
        lw=LINE_WIDTHS["primary"],
        label="Cohort slope",
        zorder=2,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(subjects)
    ax.set_xlabel("Subject")
    ax.set_ylabel(
        r"Subject-specific slope (week$^{-1}$)"
    )

    spec = ANNOTATIONS[panel]
    text = spec["format"].format(
        nneg=int(
            test["n_negative_subject_slopes"]
        ),
        n=int(test["n_subjects"]),
        p=float(test["two_sided_p"]),
    )
    add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    add_panel_legend(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 6E
# =============================================================================

def plot_figure6E(
    binwise: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure6E_binwise_slopes"

    require_columns(
        binwise,
        [
            "x_center",
            "observed_slope_per_week",
            "bootstrap_q025",
            "bootstrap_q975",
        ],
        panel,
    )

    d = binwise.copy().sort_values("x_center")
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    x = d["x_center"].to_numpy(float)
    y = d["observed_slope_per_week"].to_numpy(float)
    lo = d["bootstrap_q025"].to_numpy(float)
    hi = d["bootstrap_q975"].to_numpy(float)

    ax.errorbar(
        x,
        y,
        yerr=np.vstack([y - lo, hi - y]),
        fmt="o",
        color=COLORS["primary"],
        lw=LINE_WIDTHS["errorbar"],
        ms=MARKER_SIZES["standard"],
        capsize=MARKER_SIZES["errorbar_capsize"],
    )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        r"Cross-covariance slope (week$^{-1}$)"
    )

    n_total = int(len(d))
    n_negative = int(
        np.sum(
            pd.to_numeric(
                d["observed_slope_per_week"],
                errors="coerce",
            ).to_numpy(float) < 0
        )
    )
    n_ci_negative = int(
        np.sum(
            pd.to_numeric(
                d["bootstrap_q975"],
                errors="coerce",
            ).to_numpy(float) < 0
        )
    )

    spec = ANNOTATIONS[panel]
    text = spec["format"].format(
        nneg=n_negative,
        nci=n_ci_negative,
        n=n_total,
    )
    add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S7A
# =============================================================================

def plot_supp_s7A(
    core: pd.DataFrame,
    core_min_subject_fraction: float,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS7A_core_selection"

    require_columns(
        core,
        [
            "x_center",
            "selected_for_primary_core",
            "coverage_fraction_dt1",
            "coverage_fraction_dt2",
            "coverage_fraction_dt3",
            "coverage_fraction_dt4",
            "coverage_fraction_dt5",
        ],
        panel,
    )

    d = core.copy().sort_values("x_center")
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)

    x = d["x_center"].to_numpy(float)

    for dt, color in zip(
        [1, 2, 3, 4, 5],
        [
            "#0072B2",
            "#009E73",
            "#E69F00",
            "#CC79A7",
            "#D55E00",
        ],
    ):
        ax.plot(
            x,
            d[f"coverage_fraction_dt{dt}"].to_numpy(float),
            lw=LINE_WIDTHS["standard"],
            marker="o",
            ms=3.5,
            color=color,
            label=f"{dt} week" if dt == 1 else f"{dt} weeks",
        )

    selected = d[
        d["selected_for_primary_core"].astype(str)
        .str.lower()
        .isin(["true", "1"])
    ]

    if not selected.empty:
        ax.axvspan(
            float(selected["x_left"].min()),
            float(selected["x_right"].max()),
            color=COLORS["core"],
            alpha=ALPHAS["core_span"],
            linewidth=0,
            label="Primary abundance core",
        )

    ax.axhline(
        float(core_min_subject_fraction),
        color=COLORS["zero"],
        lw=LINE_WIDTHS["zero"],
        ls="--",
        alpha=ALPHAS["zero_line"],
    )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel("Subject coverage fraction")

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    add_panel_legend(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S7B
# =============================================================================

def plot_supp_s7B(
    bootstrap: pd.DataFrame,
    slope: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS7B_bootstrap_slope_distribution"

    d = bootstrap[
        bootstrap["estimator"].astype(str).eq(
            "equal_bin_equal_subject"
        )
        & bootstrap["metric"].astype(str).eq(
            "cross_cov"
        )
    ].copy()

    vals = pd.to_numeric(
        d["slope_per_week"],
        errors="coerce",
    ).to_numpy(float)
    vals = vals[np.isfinite(vals)]

    if vals.size == 0:
        raise ValueError("No finite primary bootstrap slopes.")

    source = pd.DataFrame(
        {"slope_per_week": vals}
    )
    write_source_data(source, source_dir, panel)

    row = primary_slope_row(
        slope,
        "cross_cov",
    )

    fig, ax = make_figure(panel)
    maybe_zero_line(ax, orientation="vertical")

    ax.hist(
        vals,
        bins=HISTOGRAM["bins"],
        density=True,
        color=COLORS["bootstrap"],
        alpha=ALPHAS["histogram"],
        edgecolor="none",
    )

    obs = float(
        row["observed_slope_per_week"]
    )
    lo = float(row["bootstrap_q025"])
    hi = float(row["bootstrap_q975"])

    ax.axvspan(
        lo,
        hi,
        color=COLORS["primary"],
        alpha=ALPHAS["bootstrap_ci_span"],
        linewidth=0,
    )
    ax.axvline(
        lo,
        color=COLORS["primary"],
        lw=LINE_WIDTHS["ci_boundary"],
        ls=":",
        alpha=0.75,
    )
    ax.axvline(
        hi,
        color=COLORS["primary"],
        lw=LINE_WIDTHS["ci_boundary"],
        ls=":",
        alpha=0.75,
    )
    ax.axvline(
        obs,
        color=COLORS["primary"],
        lw=LINE_WIDTHS["primary"],
    )

    ax.set_xlabel(
        r"Bootstrap cross-covariance slope (week$^{-1}$)"
    )
    ax.set_ylabel("Bootstrap density")

    spec = ANNOTATIONS[panel]
    text = spec["format"].format(
        slope=obs,
        lo=lo,
        hi=hi,
    )
    add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S7C
# =============================================================================

def plot_supp_s7C(
    models: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS7C_model_comparison"

    d = models[
        models["estimator"].astype(str).eq(
            "equal_bin_equal_subject"
        )
        & models["metric"].astype(str).eq(
            "cross_cov"
        )
    ].copy()

    if d.empty:
        raise ValueError(
            "No primary cross_cov model-comparison rows."
        )

    order = [
        "constant",
        "signed_linear",
        "positive_incremental",
    ]
    d["model_order"] = d["model"].map(
        {m: i for i, m in enumerate(order)}
    )
    d = d.sort_values("model_order")

    if "delta_aicc" not in d.columns:
        d["delta_aicc"] = (
            pd.to_numeric(
                d["observed_aicc"],
                errors="coerce",
            )
            - pd.to_numeric(
                d["observed_aicc"],
                errors="coerce",
            ).min()
        )
    else:
        d["delta_aicc"] = pd.to_numeric(
            d["delta_aicc"],
            errors="coerce",
        )

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)

    x = np.arange(len(d))
    y = d["delta_aicc"].to_numpy(float)

    colors = [
        COLORS.get(str(model), "#777777")
        for model in d["model"]
    ]

    ax.bar(
        x,
        y,
        color=colors,
        width=MARKER_SIZES["bar_width"],
        alpha=ALPHAS["bar"],
    )
    ax.scatter(
        x,
        y,
        color=colors,
        s=26,
        zorder=3,
    )

    for xi, yi in zip(x, y):
        ax.text(
            xi,
            yi + max(0.35, 0.025 * max(y.max(), 1.0)),
            f"{yi:.1f}",
            ha="center",
            va="bottom",
            fontsize=FONT["annotation_pt"],
        )

    labels = [
        "Constant",
        "Signed\nlinear",
        "Positive\nincremental",
    ]

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r"$\Delta$AICc relative to best model")
    ax.set_ylim(
        0.0,
        max(1.0, float(y.max()) * 1.16),
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S7D
# =============================================================================

def plot_supp_s7D(
    subject_core: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS7D_subject_temporal_curves"

    d = subject_core[
        subject_core["estimator"].astype(str).eq(
            "equal_bin_equal_subject"
        )
        & subject_core["metric"].astype(str).eq(
            "cross_cov"
        )
    ].copy().sort_values(["subject", "dt"])

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)

    for subject, g in d.groupby("subject", sort=True):
        ax.plot(
            g["dt"].to_numpy(float),
            g["value"].to_numpy(float),
            color=COLORS["subject"],
            alpha=ALPHAS["subject_line"],
            lw=LINE_WIDTHS["secondary"],
            marker="o",
            ms=MARKER_SIZES["subject"],
        )

    cohort = (
        d.groupby("dt", as_index=False)
        .agg(value=("value", "mean"))
        .sort_values("dt")
    )
    ax.plot(
        cohort["dt"].to_numpy(float),
        cohort["value"].to_numpy(float),
        color=COLORS["primary"],
        lw=LINE_WIDTHS["primary"],
        marker="o",
        ms=MARKER_SIZES["standard"],
        label="Cohort mean",
        zorder=4,
    )

    # Dummy handle describing the thin gray subject curves.
    ax.plot(
        [],
        [],
        color=COLORS["subject"],
        alpha=ALPHAS["subject_line"],
        lw=LINE_WIDTHS["secondary"],
        marker="o",
        ms=MARKER_SIZES["subject"],
        label="Individual subjects",
    )

    ax.set_xlabel("Temporal lag (weeks)")
    ax.set_ylabel(
        "Subject core-averaged cross covariance"
    )
    ax.set_xticks([1, 2, 3, 4, 5])

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    add_panel_legend(ax, panel)
    save_figure(fig, outdir, panel)


def write_plot_contract_audit(
    outdir: Path,
    contract: Mapping[str, object],
    core: pd.DataFrame,
    complete: pd.DataFrame,
    support_audit: pd.DataFrame,
    binwise: pd.DataFrame,
    subject_slopes: pd.DataFrame,
) -> None:
    selected = core[
        bool_series(
            core["selected_for_primary_core"]
        )
    ].copy()

    selected_subjects = complete[
        bool_series(
            complete["selected_for_primary_complete_case"]
        )
    ].copy()

    primary_subject_slopes = subject_slopes[
        subject_slopes["estimator"].astype(str).eq(
            "equal_bin_equal_subject"
        )
        & subject_slopes["metric"].astype(str).eq("cross_cov")
    ].copy()

    row = {
        **contract,
        "core_rows_marked_selected": int(len(selected)),
        "complete_case_rows_marked_selected": int(len(selected_subjects)),
        "support_audit_min_subject_bin_n": float(
            pd.to_numeric(
                support_audit["min_subject_bin_n"],
                errors="coerce",
            ).min()
        ),
        "n_binwise_slopes": int(len(binwise)),
        "n_binwise_observed_negative": int(
            (
                pd.to_numeric(
                    binwise["observed_slope_per_week"],
                    errors="coerce",
                ) < 0
            ).sum()
        ),
        "n_binwise_ci_entirely_negative": int(
            (
                pd.to_numeric(
                    binwise["bootstrap_q975"],
                    errors="coerce",
                ) < 0
            ).sum()
        ),
        "n_subject_slopes": int(len(primary_subject_slopes)),
        "n_subject_slopes_negative": int(
            (
                pd.to_numeric(
                    primary_subject_slopes["slope_per_week"],
                    errors="coerce",
                ) < 0
            ).sum()
        ),
    }

    if int(row["core_rows_marked_selected"]) != int(contract["n_core_bins"]):
        raise ValueError(
            "Plotter audit: selected core-bin count differs from Step 13 config."
        )
    if int(row["complete_case_rows_marked_selected"]) != int(
        contract["n_complete_case_subjects"]
    ):
        raise ValueError(
            "Plotter audit: complete-case subject count differs from Step 13 config."
        )

    pd.DataFrame([row]).to_csv(
        outdir / "00_plot_contract_audit.csv",
        index=False,
    )


# =============================================================================
# CLI / main
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Create single-panel publication plots for final Step-13 "
            "temporal-scaling analysis."
        ),
    )
    p.add_argument(
        "--step13-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    configure_matplotlib()

    step13_dir = (
        args.step13_dir
        .expanduser()
        .resolve(strict=True)
    )
    outdir = (
        args.outdir
        .expanduser()
        .resolve()
    )

    main_dir = outdir / "main_figure"
    supplementary_dir = outdir / "supplementary_figure"
    source_dir = outdir / "figure_source_data"

    data = load_inputs(step13_dir)

    plot_figure6A(
        data["cohort"],
        data["slope"],
        main_dir,
        source_dir,
    )

    plot_figure6B(
        data["cohort"],
        main_dir,
        source_dir,
    )

    plot_figure6C(
        data["cohort"],
        main_dir,
        source_dir,
    )

    plot_figure6D(
        data["subject_slopes"],
        data["slope"],
        data["signflip"],
        main_dir,
        source_dir,
    )

    plot_figure6E(
        data["binwise"],
        main_dir,
        source_dir,
    )

    plot_supp_s7A(
        data["core"],
        float(
            data["contract"]["core_min_subject_fraction"]
        ),
        supplementary_dir,
        source_dir,
    )

    plot_supp_s7B(
        data["bootstrap"],
        data["slope"],
        supplementary_dir,
        source_dir,
    )

    plot_supp_s7C(
        data["models"],
        supplementary_dir,
        source_dir,
    )

    plot_supp_s7D(
        data["subject_core"],
        supplementary_dir,
        source_dir,
    )

    write_plot_contract_audit(
        outdir,
        data["contract"],
        data["core"],
        data["complete"],
        data["support_audit"],
        data["binwise"],
        data["subject_slopes"],
    )

    print("\n[DONE] Step 13 publication plots")
    print("Main Figure 6 panels    :", main_dir)
    print("Supplementary Figure 7 :", supplementary_dir)
    print("Exact plotted data      :", source_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
