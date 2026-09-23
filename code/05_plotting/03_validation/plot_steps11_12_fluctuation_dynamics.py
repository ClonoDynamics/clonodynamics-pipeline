#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_steps11_12_fluctuation_dynamics.py
==========================================

Publication plotter for final ClonoDynamics Step 11–12 results.

Final contract:
- Step 11 v2 fixed-support subject bootstrap;
- Step 12 v5 direct matched native Step-8 support;
- pseudo missing bins are never interpolated or imputed.

The script writes SINGLE-PANEL PDF and PNG files only. Panels are intended for
final assembly outside Python.

Recommended main figure
-----------------------
Figure5A_real_crosscov_by_abundance
    Genuine longitudinal cross-replicate displacement covariance versus
    xmid_latent for dt = 1–5.

Figure5B_real_fluctuation_decomposition_dt1
    Genuine longitudinal cross-replicate covariance, mean within-replicate
    variance, and replicate-specific excess at dt = 1.

Figure5C_real_vs_pseudo_dt1
    Matched-support longitudinal versus pseudo cross-replicate covariance at
    dt = 1.

Figure5D_matched_crosscov_by_lag
    Abundance-averaged matched-support longitudinal and pseudo covariance at
    dt = 1–5, with biological bootstrap and pseudo randomization intervals.

Figure5E_joint_null_distribution
    Distribution of the equal-lag pseudo technical-null statistic across
    configurations eligible at all five lags, with the longitudinal estimate
    and subject-bootstrap interval.

Recommended supplementary figure
--------------------------------
SuppFigS6A_shared_fraction_by_abundance
    Longitudinal shared fraction versus xmid_latent for dt = 1–5.

SuppFigS6B_common4_support_by_lag
    Fraction of the full transition universe in common4 versus dt.

SuppFigS6C_grid_exceedance_by_lag
    Fraction of the matched xmid grid for which longitudinal covariance exceeds
    the pseudo 97.5th randomization percentile.

SuppFigS6D_pseudo_eligibility_by_lag
    Fraction of the 2,000 pseudo configurations eligible for the matched-support
    empirical null at each lag.

Inputs
------
--step11-dir
    Final Step-11 output directory.

--step12-dir
    Final Step-12 output directory.

Typical run
-----------
python3 plot_steps11_12_fluctuation_dynamics.py \
    --step11-dir ./dataset_longitudinal_results/11-cross_replicate_fluctuation_dynamics \
    --step12-dir ./dataset_longitudinal_results/12-longitudinal_vs_pseudo_cross_replicate_fluctuations \
    --outdir ./figures/step11_step12_fluctuations

All figure dimensions, colors, line widths, axis limits, annotation positions,
and legend positions are editable in the dictionaries below.
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

FONT_FAMILY = "Arial"
FONT_SIZE = 8
AXIS_LABEL_SIZE = 8
TICK_LABEL_SIZE = 7
LEGEND_FONT_SIZE = 7
ANNOTATION_FONT_SIZE = 7

PDF_FONT_TYPE = 42
PS_FONT_TYPE = 42
PNG_DPI = 600

SAVE_FORMATS = ("pdf", "png")

FIGSIZE_INCHES: Dict[str, Tuple[float, float]] = {
    "Figure5A_real_crosscov_by_abundance": (3.2, 2.4),
    "Figure5B_real_fluctuation_decomposition_dt1": (3.2, 2.4),
    "Figure5C_real_vs_pseudo_dt1": (3.2, 2.4),
    "Figure5D_matched_crosscov_by_lag": (3.2, 2.4),
    "Figure5E_joint_null_distribution": (3.2, 2.4),
    "SuppFigS6A_shared_fraction_by_abundance": (3.2, 2.4),
    "SuppFigS6B_common4_support_by_lag": (3.2, 2.4),
    "SuppFigS6C_grid_exceedance_by_lag": (3.2, 2.4),
    "SuppFigS6D_pseudo_eligibility_by_lag": (3.2, 2.4),
}

COLORS = {
    "dt1": "#0072B2",
    "dt2": "#009E73",
    "dt3": "#E69F00",
    "dt4": "#CC79A7",
    "dt5": "#D55E00",
    "cross": "#0072B2",
    "same": "#D55E00",
    "excess": "#CC79A7",
    "zero": "#666666",
    "support": "#333333",
    "bootstrap_hist": "#E69F00",
    "longitudinal": "#0072B2",
    "pseudo": "#D55E00",
    "pseudo_fill": "#F4C29A",
    "longitudinal_fill": "#B9D7F0",
}

LINESTYLES = {
    "dt1": "-",
    "dt2": "--",
    "dt3": "-.",
    "dt4": ":",
    "dt5": (0, (5, 2)),
    "cross": "-",
    "same": "-",
    "excess": "--",
    "longitudinal": "-",
    "pseudo": "--",
    "zero": "--",
}

MARKERS = {
    "dt1": "o",
    "dt2": "s",
    "dt3": "^",
    "dt4": "D",
    "dt5": "v",
    "longitudinal": "o",
    "pseudo": "s",
}

LINE_WIDTH = 1
PRIMARY_LINE_WIDTH = 1.5
SECONDARY_LINE_WIDTH = 0.5
ZERO_LINE_WIDTH = 0.5
ERRORBAR_LINE_WIDTH = 0.5
IDENTITY_LINE_WIDTH = 1.0
CAP_SIZE = 2.0
MARKER_SIZE = 2.5
POINT_SIZE = 3.0
SCATTER_SIZE = 28

RIBBON_ALPHA = 0.28
PSEUDO_RIBBON_ALPHA = 0.28
COHORT_CI_ALPHA = 0.10
HIST_ALPHA = 0.75

ZERO_LINE_ALPHA = 0.65
IDENTITY_LINE_ALPHA = 0.70


LEGEND_LOC = "upper left"
LEGEND_BBOX = (1.02, 1.0)


SHOW_TITLES = False

PLOT_TITLES = {
    "Figure5A_real_crosscov_by_abundance":
        "Longitudinal shared fluctuation covariance",
    "Figure5B_real_fluctuation_decomposition_dt1":
        "Longitudinal fluctuation decomposition",
    "Figure5C_real_vs_pseudo_dt1":
        "Longitudinal versus pseudo technical null at one week",
    "Figure5D_matched_crosscov_by_lag":
        "Matched-support cross-replicate covariance across lags",
    "Figure5E_joint_null_distribution":
        "Equal-lag technical-null comparison",
    "SuppFigS6A_shared_fraction_by_abundance":
        "Shared fraction across abundance and lag",
    "SuppFigS6B_common4_support_by_lag":
        "common4 support across temporal lags",
    "SuppFigS6C_grid_exceedance_by_lag":
        "Matched native-bin exceedance of pseudo 97.5th percentile",
    "SuppFigS6D_pseudo_eligibility_by_lag":
        "Pseudo matched-support eligibility",
}

XLIMS = {
    "Figure5A_real_crosscov_by_abundance": None,
    "Figure5B_real_fluctuation_decomposition_dt1": None,
    "Figure5C_real_vs_pseudo_dt1": None,
    "Figure5D_matched_crosscov_by_lag": (0.7, 5.3),
    "Figure5E_joint_null_distribution": None,
    "SuppFigS6A_shared_fraction_by_abundance": None,
    "SuppFigS6B_common4_support_by_lag": (0.7, 5.3),
    "SuppFigS6C_grid_exceedance_by_lag": (0.7, 5.3),
    "SuppFigS6D_pseudo_eligibility_by_lag": (0.7, 5.3),
}

YLIMS = {
    "Figure5A_real_crosscov_by_abundance": None,
    "Figure5B_real_fluctuation_decomposition_dt1": None,
    "Figure5C_real_vs_pseudo_dt1": None,
    "Figure5D_matched_crosscov_by_lag": None,
    "Figure5E_joint_null_distribution": None,
    "SuppFigS6A_shared_fraction_by_abundance": None,
    "SuppFigS6B_common4_support_by_lag": (0.0, None),
    "SuppFigS6C_grid_exceedance_by_lag": (0.0, 1.05),
    "SuppFigS6D_pseudo_eligibility_by_lag": (0.0, 1.05),
}

# Display switches
SHOW_FIG5A_RIBBONS = False
SHOW_FIG5B_RIBBONS = True
SHOW_FIG5C_LONGITUDINAL_CI = True
SHOW_FIG5C_PSEUDO_ENVELOPE = True

SHOW_ZERO_LINE = True

HISTOGRAM_BINS = 35

ANNOTATIONS = {
    "Figure5C_real_vs_pseudo_dt1": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "ha": "left",
        "va": "bottom",
        "format": (
            r"$P_{{emp}}$ = {p:.5f}" "\n"
            "{nout:d}/{ntot:d} matched bins above pseudo 97.5%"
        ),
    },
    "Figure5D_matched_crosscov_by_lag": {
        "enabled": False,
        "xy_axes": (0.03, 0.04),
        "ha": "left",
        "va": "bottom",
        "format": "",
    },
    "Figure5E_joint_null_distribution": {
        "enabled": True,
        "xy_axes": (0.77, 0.95),
        "ha": "right",
        "va": "top",
        "format": (
            "longitudinal = {real:.3f}\n"
            "pseudo median = {pseudo:.3f}\n"
            r"$P_{{emp}}$ = {p:.4f}"
        ),
    },
}


# =============================================================================
# Matplotlib helpers
# =============================================================================

def configure_matplotlib() -> None:
    mpl.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.size": FONT_SIZE,
        "axes.labelsize": AXIS_LABEL_SIZE,
        "xtick.labelsize": TICK_LABEL_SIZE,
        "ytick.labelsize": TICK_LABEL_SIZE,
        "legend.fontsize": LEGEND_FONT_SIZE,
        "pdf.fonttype": PDF_FONT_TYPE,
        "ps.fonttype": PS_FONT_TYPE,
        "axes.linewidth": 0.9,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 4.0,
        "ytick.major.size": 4.0,
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
    ax.tick_params(direction="out")
    ax.grid(False)


def maybe_zero_line(
    ax: plt.Axes,
    orientation: str = "horizontal",
) -> None:
    if not SHOW_ZERO_LINE:
        return

    kwargs = dict(
        color=COLORS["zero"],
        lw=ZERO_LINE_WIDTH,
        ls=LINESTYLES["zero"],
        alpha=ZERO_LINE_ALPHA,
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
            pad=8,
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


def add_annotation(
    ax: plt.Axes,
    panel: str,
    text: str,
) -> None:
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
        fontsize=ANNOTATION_FONT_SIZE,
    )


def save_figure(
    fig: plt.Figure,
    outdir: Path,
    panel: str,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    for fmt in SAVE_FORMATS:
        path = outdir / f"{panel}.{fmt}"
        kwargs = dict(
            bbox_inches="tight",
            pad_inches=0.06,
        )
        if fmt.lower() == "png":
            kwargs["dpi"] = PNG_DPI

        fig.savefig(path, **kwargs)

    plt.close(fig)


# =============================================================================
# IO helpers
# =============================================================================

def read_json_required(
    path: Path,
    label: str,
) -> dict:
    if not path.is_file():
        raise FileNotFoundError(
            f"{label} not found: {path}"
        )
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_final_contracts(
    step11_dir: Path,
    step12_dir: Path,
) -> Dict[str, object]:
    step11 = read_json_required(
        step11_dir / "00_run_config.json",
        "Step 11 run config",
    )
    step12 = read_json_required(
        step12_dir / "00_run_config.json",
        "Step 12 run config",
    )

    if step11.get("primary_support") != "common4":
        raise ValueError("Final plotter requires Step 11 common4 support.")
    if step11.get("primary_conditioning") != "xmid_latent":
        raise ValueError("Final plotter requires Step 11 xmid_latent conditioning.")
    est11 = step11.get("primary_estimands", {})
    if est11.get("cross_cov_truncated_at_zero") is not False:
        raise ValueError("Final plotter requires signed/untruncated Step 11 cross_cov.")

    if step12.get("primary_support") != "common4":
        raise ValueError("Final plotter requires Step 12 common4 support.")
    if step12.get("primary_conditioning") != "xmid_latent":
        raise ValueError("Final plotter requires Step 12 xmid_latent conditioning.")

    matching = step12.get("matching", {})
    if matching.get("pseudo_interpolation") is not False:
        raise ValueError(
            "Final plotter requires Step 12 v5 pseudo_interpolation=false."
        )
    if matching.get("pseudo_imputation") is not False:
        raise ValueError(
            "Final plotter requires Step 12 v5 pseudo_imputation=false."
        )
    if matching.get("pseudo_stable_bin_fraction_denominator") != (
        "full Step-8 randomization ensemble"
    ):
        raise ValueError(
            "Final plotter requires full Step-8 ensemble denominator "
            "for stable native-bin support."
        )

    pseudo_null = step12.get("pseudo_null", {})
    if pseudo_null.get("production_requires_complete_matched_native_bins") is not True:
        raise ValueError(
            "Final plotter requires complete matched native-bin pseudo coverage."
        )
    if float(
        pseudo_null.get("min_matched_native_bin_coverage", np.nan)
    ) != 1.0:
        raise ValueError(
            "Final plotter requires min matched native-bin coverage = 1.0."
        )

    return {
        "step11_script_version": step11.get("script_version"),
        "step11_analysis_signature": step11.get("analysis_signature"),
        "step12_script_version": step12.get("script_version"),
        "step12_analysis_signature": step12.get("analysis_signature"),
        "full_pseudo_configuration_count": int(
            pseudo_null.get("full_configuration_count", -1)
        ),
        "stable_bin_min_valid_fraction": matching.get(
            "pseudo_stable_bin_min_valid_fraction"
        ),
        "min_shared_native_bins": matching.get(
            "min_shared_native_bins"
        ),
        "min_matched_native_bin_coverage": pseudo_null.get(
            "min_matched_native_bin_coverage"
        ),
    }


def read_csv_required(
    path: Path,
    label: str,
) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(
            f"{label} not found: {path}"
        )

    d = pd.read_csv(path)

    if d.empty:
        raise ValueError(
            f"{label} is empty: {path}"
        )

    return d


def require_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
    label: str,
) -> None:
    missing = [
        c for c in columns
        if c not in frame.columns
    ]
    if missing:
        raise ValueError(
            f"{label}: missing required columns {missing}; "
            f"available={list(frame.columns)}"
        )


def bool_series(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(
            ["true", "1", "t", "yes", "y"]
        )
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
        source_dir
        / f"{panel}_source_data.csv",
        index=False,
    )


def load_inputs(
    step11_dir: Path,
    step12_dir: Path,
) -> Dict[str, object]:
    contract = validate_final_contracts(
        step11_dir,
        step12_dir,
    )
    return {
        "contract": contract,
        "step11_by_bin_dt": read_csv_required(
            step11_dir
            / "02_fluctuation_by_bin_dt.csv",
            "Step 11 fluctuation by bin/dt",
        ),
        "step11_support": read_csv_required(
            step11_dir
            / "06_support_by_dt.csv",
            "Step 11 support by dt",
        ),
        "step12_summary": read_csv_required(
            step12_dir
            / "00_analysis_summary.csv",
            "Step 12 analysis summary",
        ),
        "step12_support": read_csv_required(
            step12_dir
            / "01_shared_support_by_dt.csv",
            "Step 12 shared support",
        ),
        "step12_pointwise": read_csv_required(
            step12_dir
            / "02_pointwise_cross_covariance_comparison.csv",
            "Step 12 pointwise comparison",
        ),
        "step12_pseudo_metrics": read_csv_required(
            step12_dir
            / "03_pseudo_configuration_matched_support_metrics.csv",
            "Step 12 pseudo configuration metrics",
        ),
        "step12_tests": read_csv_required(
            step12_dir
            / "04_empirical_tests_by_dt.csv",
            "Step 12 empirical tests",
        ),
        "step12_joint": read_csv_required(
            step12_dir
            / "06_joint_comparison.csv",
            "Step 12 joint comparison",
        ),
        "step12_audit": read_csv_required(
            step12_dir
            / "07_support_audit.csv",
            "Step 12 support audit",
        ),
        "step12_stable_bins": read_csv_required(
            step12_dir
            / "07b_stable_native_bins_by_dt.csv",
            "Step 12 stable native-bin audit",
        ),
        "step12_native_global": read_csv_required(
            step12_dir
            / "08_native_global_by_dt.csv",
            "Step 12 native global by dt",
        ),
    }


# =============================================================================
# Figure 5A
# =============================================================================

def plot_figure5A(
    step11: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure5A_real_crosscov_by_abundance"

    require_columns(
        step11,
        [
            "dt",
            "x_center",
            "cross_cov",
            "cross_cov_ci025",
            "cross_cov_ci975",
            "meets_min_n",
        ],
        panel,
    )

    d = step11[
        bool_series(step11["meets_min_n"])
    ].copy()

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    for dt in [1, 2, 3, 4, 5]:
        g = d[
            pd.to_numeric(
                d["dt"],
                errors="coerce",
            ).eq(dt)
        ].copy().sort_values(
            "x_center"
        )

        if g.empty:
            continue

        key = f"dt{dt}"
        color = COLORS[key]

        x = pd.to_numeric(
            g["x_center"],
            errors="coerce",
        ).to_numpy(float)
        y = pd.to_numeric(
            g["cross_cov"],
            errors="coerce",
        ).to_numpy(float)

        if SHOW_FIG5A_RIBBONS:
            lo = pd.to_numeric(
                g["cross_cov_ci025"],
                errors="coerce",
            ).to_numpy(float)
            hi = pd.to_numeric(
                g["cross_cov_ci975"],
                errors="coerce",
            ).to_numpy(float)

            ax.fill_between(
                x,
                lo,
                hi,
                color=color,
                alpha=RIBBON_ALPHA,
                linewidth=0,
            )

        ax.plot(
            x,
            y,
            color=color,
            lw=LINE_WIDTH,
            ls=LINESTYLES[key],
            marker=MARKERS[key],
            ms=MARKER_SIZE,
            label=f"{dt} week"
            if dt == 1
            else f"{dt} weeks",
        )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        "Cross-replicate displacement\ncovariance"
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )

    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Figure 5B
# =============================================================================

def plot_figure5B(
    step11: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure5B_real_fluctuation_decomposition_dt1"

    require_columns(
        step11,
        [
            "dt",
            "x_center",
            "cross_cov",
            "same_var_mean",
            "replicate_specific_excess",
            "cross_cov_ci025",
            "cross_cov_ci975",
            "same_var_mean_ci025",
            "same_var_mean_ci975",
            "replicate_specific_excess_ci025",
            "replicate_specific_excess_ci975",
            "meets_min_n",
        ],
        panel,
    )

    d = step11[
        pd.to_numeric(
            step11["dt"],
            errors="coerce",
        ).eq(1)
        & bool_series(
            step11["meets_min_n"]
        )
    ].copy().sort_values(
        "x_center"
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    x = pd.to_numeric(
        d["x_center"],
        errors="coerce",
    ).to_numpy(float)

    specs = [
        (
            "cross_cov",
            "cross_cov_ci025",
            "cross_cov_ci975",
            "Cross-replicate covariance",
            COLORS["cross"],
            LINESTYLES["cross"],
            PRIMARY_LINE_WIDTH,
        ),
        (
            "same_var_mean",
            "same_var_mean_ci025",
            "same_var_mean_ci975",
            "Mean within-replicate variance",
            COLORS["same"],
            LINESTYLES["same"],
            LINE_WIDTH,
        ),
        (
            "replicate_specific_excess",
            "replicate_specific_excess_ci025",
            "replicate_specific_excess_ci975",
            "Replicate-specific excess",
            COLORS["excess"],
            LINESTYLES["excess"],
            LINE_WIDTH,
        ),
    ]

    for (
        metric,
        lo_col,
        hi_col,
        label,
        color,
        ls,
        lw,
    ) in specs:

        y = pd.to_numeric(
            d[metric],
            errors="coerce",
        ).to_numpy(float)

        if SHOW_FIG5B_RIBBONS:
            lo = pd.to_numeric(
                d[lo_col],
                errors="coerce",
            ).to_numpy(float)
            hi = pd.to_numeric(
                d[hi_col],
                errors="coerce",
            ).to_numpy(float)

            ax.fill_between(
                x,
                lo,
                hi,
                color=color,
                alpha=RIBBON_ALPHA,
                linewidth=0,
            )

        ax.plot(
            x,
            y,
            color=color,
            lw=lw,
            ls=ls,
            label=label,
        )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        "Displacement variance / covariance"
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )

    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Figure 5C
# =============================================================================

def plot_figure5C(
    pointwise: pd.DataFrame,
    tests: pd.DataFrame,
    audit: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure5C_real_vs_pseudo_dt1"

    require_columns(
        pointwise,
        [
            "dt",
            "x",
            "longitudinal_cross_cov",
            "longitudinal_cross_cov_ci025",
            "longitudinal_cross_cov_ci975",
            "pseudo_cross_cov_median",
            "pseudo_cross_cov_q025",
            "pseudo_cross_cov_q975",
        ],
        panel,
    )

    d = pointwise[
        pd.to_numeric(
            pointwise["dt"],
            errors="coerce",
        ).eq(1)
    ].copy().sort_values(
        "x"
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    x = pd.to_numeric(
        d["x"],
        errors="coerce",
    ).to_numpy(float)

    real = pd.to_numeric(
        d["longitudinal_cross_cov"],
        errors="coerce",
    ).to_numpy(float)

    pseudo = pd.to_numeric(
        d["pseudo_cross_cov_median"],
        errors="coerce",
    ).to_numpy(float)

    if SHOW_FIG5C_PSEUDO_ENVELOPE:
        pseudo_lo = pd.to_numeric(
            d["pseudo_cross_cov_q025"],
            errors="coerce",
        ).to_numpy(float)
        pseudo_hi = pd.to_numeric(
            d["pseudo_cross_cov_q975"],
            errors="coerce",
        ).to_numpy(float)

        ax.fill_between(
            x,
            pseudo_lo,
            pseudo_hi,
            color=COLORS["pseudo_fill"],
            alpha=0.38,
            linewidth=0,
            label="Pseudo 2.5–97.5% randomization interval",
        )

    if SHOW_FIG5C_LONGITUDINAL_CI:
        real_lo = pd.to_numeric(
            d["longitudinal_cross_cov_ci025"],
            errors="coerce",
        ).to_numpy(float)
        real_hi = pd.to_numeric(
            d["longitudinal_cross_cov_ci975"],
            errors="coerce",
        ).to_numpy(float)

        ax.fill_between(
            x,
            real_lo,
            real_hi,
            color=COLORS["longitudinal_fill"],
            alpha=0.38,
            linewidth=0,
            label="Longitudinal 95% subject-bootstrap CI",
        )

    ax.plot(
        x,
        pseudo,
        color=COLORS["pseudo"],
        lw=LINE_WIDTH,
        ls=LINESTYLES["pseudo"],
        label="Pseudo median",
        zorder=3,
    )

    ax.plot(
        x,
        real,
        color=COLORS["longitudinal"],
        lw=PRIMARY_LINE_WIDTH,
        ls=LINESTYLES["longitudinal"],
        label="Longitudinal",
        zorder=4,
    )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        "Cross-replicate displacement\ncovariance"
    )

    p_row = tests[
        pd.to_numeric(
            tests["dt"],
            errors="coerce",
        ).eq(1)
    ]
    a_row = audit[
        pd.to_numeric(
            audit["dt"],
            errors="coerce",
        ).eq(1)
    ]

    if not p_row.empty and not a_row.empty:
        spec = ANNOTATIONS[panel]
        text = spec["format"].format(
            p=float(
                p_row.iloc[0]["empirical_p"]
            ),
            nout=int(
                a_row.iloc[0][
                    "n_longitudinal_points_above_pseudo_q975"
                ]
            ),
            ntot=int(
                a_row.iloc[0][
                    "n_matched_native_bins"
                ]
            ),
        )
        add_annotation(
            ax,
            panel,
            text,
        )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )

    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Figure 5D
# =============================================================================

def plot_figure5D(
    tests: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure5D_matched_crosscov_by_lag"

    require_columns(
        tests,
        [
            "dt",
            "longitudinal_value",
            "bootstrap_q025",
            "bootstrap_q975",
            "pseudo_median",
            "pseudo_q025",
            "pseudo_q975",
            "empirical_p",
            "pseudo_n",
        ],
        panel,
    )

    d = tests.copy().sort_values(
        "dt"
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    x = pd.to_numeric(
        d["dt"],
        errors="coerce",
    ).to_numpy(float)

    real = pd.to_numeric(
        d["longitudinal_value"],
        errors="coerce",
    ).to_numpy(float)
    real_lo = pd.to_numeric(
        d["bootstrap_q025"],
        errors="coerce",
    ).to_numpy(float)
    real_hi = pd.to_numeric(
        d["bootstrap_q975"],
        errors="coerce",
    ).to_numpy(float)

    pseudo = pd.to_numeric(
        d["pseudo_median"],
        errors="coerce",
    ).to_numpy(float)
    pseudo_lo = pd.to_numeric(
        d["pseudo_q025"],
        errors="coerce",
    ).to_numpy(float)
    pseudo_hi = pd.to_numeric(
        d["pseudo_q975"],
        errors="coerce",
    ).to_numpy(float)

    offset = 0.07

    ax.errorbar(
        x - offset,
        real,
        yerr=np.vstack(
            [
                real - real_lo,
                real_hi - real,
            ]
        ),
        fmt=MARKERS["longitudinal"],
        color=COLORS["longitudinal"],
        lw=ERRORBAR_LINE_WIDTH,
        ms=POINT_SIZE,
        capsize=CAP_SIZE,
        label="Longitudinal",
        zorder=4,
    )

    ax.errorbar(
        x + offset,
        pseudo,
        yerr=np.vstack(
            [
                pseudo - pseudo_lo,
                pseudo_hi - pseudo,
            ]
        ),
        fmt=MARKERS["pseudo"],
        color=COLORS["pseudo"],
        lw=ERRORBAR_LINE_WIDTH,
        ms=POINT_SIZE,
        capsize=CAP_SIZE,
        label="Pseudo technical null",
        zorder=3,
    )

    ax.set_xlabel(
        "Temporal lag (weeks)"
    )
    ax.set_ylabel(
        "Matched-bin mean cross-replicate\ncovariance"
    )
    ax.set_xticks(
        sorted(
            pd.to_numeric(
                d["dt"],
                errors="coerce",
            ).dropna().astype(int).unique()
        )
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )

    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Figure 5E
# =============================================================================

def derive_joint_pseudo_distribution(
    pseudo_metrics: pd.DataFrame,
) -> pd.DataFrame:
    require_columns(
        pseudo_metrics,
        [
            "configuration_id",
            "dt",
            "eligible_for_null",
            "abundance_averaged_cross_cov",
        ],
        "Step 12 pseudo configuration metrics",
    )

    d = pseudo_metrics[
        bool_series(
            pseudo_metrics[
                "eligible_for_null"
            ]
        )
    ].copy()

    pivot = d.pivot_table(
        index="configuration_id",
        columns="dt",
        values="abundance_averaged_cross_cov",
        aggfunc="first",
    )

    needed = [1, 2, 3, 4, 5]

    for dt in needed:
        if dt not in pivot.columns:
            raise ValueError(
                f"Joint pseudo distribution missing dt={dt}"
            )

    complete = (
        pivot[needed]
        .dropna()
        .copy()
    )

    complete[
        "equal_lag_mean_cross_cov"
    ] = complete[needed].mean(
        axis=1
    )

    return (
        complete[
            ["equal_lag_mean_cross_cov"]
        ]
        .reset_index()
    )


def plot_figure5E(
    pseudo_metrics: pd.DataFrame,
    joint: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure5E_joint_null_distribution"

    require_columns(
        joint,
        [
            "longitudinal_value",
            "longitudinal_bootstrap_q025",
            "longitudinal_bootstrap_q975",
            "pseudo_median",
            "empirical_p",
        ],
        panel,
    )

    pseudo_joint = (
        derive_joint_pseudo_distribution(
            pseudo_metrics
        )
    )

    row = joint.iloc[0]

    real = float(
        row["longitudinal_value"]
    )
    real_lo = float(
        row["longitudinal_bootstrap_q025"]
    )
    real_hi = float(
        row["longitudinal_bootstrap_q975"]
    )
    pseudo_median = float(
        row["pseudo_median"]
    )
    p_emp = float(
        row["empirical_p"]
    )

    source = pseudo_joint.copy()
    source["longitudinal_value"] = real
    source["longitudinal_bootstrap_q025"] = real_lo
    source["longitudinal_bootstrap_q975"] = real_hi
    source["pseudo_median_reported"] = pseudo_median
    source["empirical_p"] = p_emp

    write_source_data(
        source,
        source_dir,
        panel,
    )

    vals = pd.to_numeric(
        pseudo_joint[
            "equal_lag_mean_cross_cov"
        ],
        errors="coerce",
    ).to_numpy(float)
    vals = vals[
        np.isfinite(vals)
    ]

    fig, ax = make_figure(panel)

    ax.hist(
        vals,
        bins=HISTOGRAM_BINS,
        density=True,
        color=COLORS["bootstrap_hist"],
        alpha=HIST_ALPHA,
        edgecolor="none",
        label="Pseudo complete-configuration null",
    )

    ax.axvspan(
        real_lo,
        real_hi,
        color=COLORS["longitudinal"],
        alpha=COHORT_CI_ALPHA,
        linewidth=0,
        label="Longitudinal 95% subject-bootstrap CI",
    )

    ax.axvline(
        real,
        color=COLORS["longitudinal"],
        lw=PRIMARY_LINE_WIDTH,
        label="Longitudinal",
    )

    ax.axvline(
        pseudo_median,
        color=COLORS["pseudo"],
        lw=LINE_WIDTH,
        ls=LINESTYLES["pseudo"],
        label="Pseudo median",
    )

    ax.set_xlabel(
        "Equal-lag mean abundance-averaged\ncross-replicate covariance"
    )
    ax.set_ylabel(
        "Pseudo randomization density"
    )

    spec = ANNOTATIONS[panel]
    text = spec["format"].format(
        real=real,
        pseudo=pseudo_median,
        p=p_emp,
    )

    add_annotation(
        ax,
        panel,
        text,
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    ax.legend(
        frameon=False,
        bbox_to_anchor=(1.02, 0.9),
        loc=LEGEND_LOC,
    )

    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S6A
# =============================================================================

def plot_supp_s6A(
    step11: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS6A_shared_fraction_by_abundance"

    require_columns(
        step11,
        [
            "dt",
            "x_center",
            "shared_fraction",
            "meets_min_n",
        ],
        panel,
    )

    d = step11[
        bool_series(
            step11["meets_min_n"]
        )
    ].copy()

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    for dt in [1, 2, 3, 4, 5]:
        g = d[
            pd.to_numeric(
                d["dt"],
                errors="coerce",
            ).eq(dt)
        ].copy().sort_values(
            "x_center"
        )

        if g.empty:
            continue

        key = f"dt{dt}"

        ax.plot(
            pd.to_numeric(
                g["x_center"],
                errors="coerce",
            ),
            pd.to_numeric(
                g["shared_fraction"],
                errors="coerce",
            ),
            color=COLORS[key],
            lw=LINE_WIDTH,
            ls=LINESTYLES[key],
            marker=MARKERS[key],
            ms=MARKER_SIZE,
            label=f"{dt} week"
            if dt == 1
            else f"{dt} weeks",
        )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        r"Shared fraction, $V_{\mathrm{cross}}/V_{\mathrm{same}}$"
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )

    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S6B
# =============================================================================

def plot_supp_s6B(
    support: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS6B_common4_support_by_lag"

    require_columns(
        support,
        [
            "dt",
            "fraction_common4",
            "n_common4",
            "n_subjects",
            "n_subject_intervals",
        ],
        panel,
    )

    d = support.copy().sort_values(
        "dt"
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)

    x = pd.to_numeric(
        d["dt"],
        errors="coerce",
    ).to_numpy(float)

    y = pd.to_numeric(
        d["fraction_common4"],
        errors="coerce",
    ).to_numpy(float)

    ax.plot(
        x,
        y,
        color=COLORS["support"],
        lw=LINE_WIDTH,
        marker="o",
        ms=MARKER_SIZE,
    )

    ax.set_xlabel(
        "Temporal lag (weeks)"
    )
    ax.set_ylabel(
        "Fraction of transitions in common4"
    )
    ax.set_xticks(
        sorted(
            d["dt"].astype(int).unique()
        )
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S6C
# =============================================================================

def plot_supp_s6C(
    audit: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS6C_grid_exceedance_by_lag"

    require_columns(
        audit,
        [
            "dt",
            "fraction_longitudinal_points_above_pseudo_q975",
        ],
        panel,
    )

    d = audit.copy().sort_values(
        "dt"
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)

    x = pd.to_numeric(
        d["dt"],
        errors="coerce",
    ).to_numpy(float)

    y = pd.to_numeric(
        d[
            "fraction_longitudinal_points_above_pseudo_q975"
        ],
        errors="coerce",
    ).to_numpy(float)

    ax.plot(
        x,
        y,
        color=COLORS["cross"],
        lw=LINE_WIDTH,
        marker="o",
        ms=MARKER_SIZE,
    )

    ax.set_xlabel(
        "Temporal lag (weeks)"
    )
    ax.set_ylabel(
        "Fraction of matched native bins above\npseudo 97.5th percentile"
    )
    ax.set_xticks(
        sorted(
            d["dt"].astype(int).unique()
        )
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S6D
# =============================================================================

def plot_supp_s6D(
    audit: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS6D_pseudo_eligibility_by_lag"

    require_columns(
        audit,
        [
            "dt",
            "fraction_eligible",
            "n_eligible_pseudo_configurations",
            "n_selected_pseudo_configurations",
        ],
        panel,
    )

    d = audit.copy().sort_values(
        "dt"
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)

    x = pd.to_numeric(
        d["dt"],
        errors="coerce",
    ).to_numpy(float)

    y = pd.to_numeric(
        d["fraction_eligible"],
        errors="coerce",
    ).to_numpy(float)

    ax.plot(
        x,
        y,
        color=COLORS["pseudo"],
        lw=LINE_WIDTH,
        marker="o",
        ms=MARKER_SIZE,
    )

    ax.set_xlabel(
        "Temporal lag (weeks)"
    )
    ax.set_ylabel(
        "Fraction of pseudo\nconfigurations eligible"
    )
    ax.set_xticks(
        sorted(
            d["dt"].astype(int).unique()
        )
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)

    save_figure(
        fig,
        outdir,
        panel,
    )


def write_plot_contract_audit(
    outdir: Path,
    contract: Mapping[str, object],
    stable_bins: pd.DataFrame,
    support_audit: pd.DataFrame,
    joint: pd.DataFrame,
    pseudo_metrics: pd.DataFrame,
) -> None:
    require_columns(
        stable_bins,
        [
            "dt",
            "bin_id",
            "valid_fraction_full_ensemble",
            "used_in_matched_analysis",
        ],
        "Step 12 stable native-bin audit",
    )

    used = stable_bins[
        bool_series(
            stable_bins["used_in_matched_analysis"]
        )
    ].copy()

    # Reconstruct the complete-across-lag pseudo null count independently.
    eligible = pseudo_metrics[
        bool_series(
            pseudo_metrics["eligible_for_null"]
        )
    ].copy()

    pivot = eligible.pivot_table(
        index="configuration_id",
        columns="dt",
        values="abundance_averaged_cross_cov",
        aggfunc="first",
    )
    needed = [1, 2, 3, 4, 5]
    if all(dt in pivot.columns for dt in needed):
        complete_n = int(
            pivot[needed].dropna().shape[0]
        )
    else:
        complete_n = 0

    reported_joint_n = np.nan
    if not joint.empty and "pseudo_n_complete_configurations" in joint.columns:
        reported_joint_n = int(
            joint.iloc[0]["pseudo_n_complete_configurations"]
        )

    if (
        np.isfinite(reported_joint_n)
        and int(reported_joint_n) != int(complete_n)
    ):
        raise ValueError(
            "Joint pseudo complete-configuration count does not match "
            "independent plotter reconstruction."
        )

    rows = []
    for dt, g in used.groupby("dt", sort=True):
        a = support_audit[
            pd.to_numeric(
                support_audit["dt"], errors="coerce"
            ).eq(int(dt))
        ]
        rows.append(
            {
                **contract,
                "dt": int(dt),
                "n_matched_native_bins": int(len(g)),
                "min_valid_fraction_full_ensemble": float(
                    pd.to_numeric(
                        g["valid_fraction_full_ensemble"],
                        errors="coerce",
                    ).min()
                ),
                "pseudo_eligible_n": (
                    int(
                        a.iloc[0][
                            "n_eligible_pseudo_configurations"
                        ]
                    )
                    if not a.empty
                    else np.nan
                ),
                "fraction_longitudinal_bins_above_pseudo_q975": (
                    float(
                        a.iloc[0][
                            "fraction_longitudinal_points_above_pseudo_q975"
                        ]
                    )
                    if not a.empty
                    else np.nan
                ),
                "joint_complete_pseudo_n": int(complete_n),
            }
        )

    pd.DataFrame(rows).to_csv(
        outdir / "00_plot_contract_audit.csv",
        index=False,
    )


# =============================================================================
# Main
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Create single-panel publication plots for final Step 11–12 "
            "cross-replicate fluctuation analyses."
        ),
    )

    p.add_argument(
        "--step11-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--step12-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    return p


def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(argv)

    configure_matplotlib()

    step11_dir = (
        args.step11_dir
        .expanduser()
        .resolve(strict=True)
    )

    step12_dir = (
        args.step12_dir
        .expanduser()
        .resolve(strict=True)
    )

    outdir = (
        args.outdir
        .expanduser()
        .resolve()
    )

    main_dir = (
        outdir / "main_figure"
    )

    supplementary_dir = (
        outdir / "supplementary_figure"
    )

    source_dir = (
        outdir / "figure_source_data"
    )

    data = load_inputs(
        step11_dir,
        step12_dir,
    )

    # Main Figure 5.
    plot_figure5A(
        data["step11_by_bin_dt"],
        main_dir,
        source_dir,
    )

    plot_figure5B(
        data["step11_by_bin_dt"],
        main_dir,
        source_dir,
    )

    plot_figure5C(
        data["step12_pointwise"],
        data["step12_tests"],
        data["step12_audit"],
        main_dir,
        source_dir,
    )

    plot_figure5D(
        data["step12_tests"],
        main_dir,
        source_dir,
    )

    plot_figure5E(
        data["step12_pseudo_metrics"],
        data["step12_joint"],
        main_dir,
        source_dir,
    )

    # Supplementary Figure 6.
    plot_supp_s6A(
        data["step11_by_bin_dt"],
        supplementary_dir,
        source_dir,
    )

    plot_supp_s6B(
        data["step11_support"],
        supplementary_dir,
        source_dir,
    )

    plot_supp_s6C(
        data["step12_audit"],
        supplementary_dir,
        source_dir,
    )

    plot_supp_s6D(
        data["step12_audit"],
        supplementary_dir,
        source_dir,
    )

    write_plot_contract_audit(
        outdir,
        data["contract"],
        data["step12_stable_bins"],
        data["step12_audit"],
        data["step12_joint"],
        data["step12_pseudo_metrics"],
    )

    print(
        "\n[DONE] Step 11–12 publication plots"
    )
    print(
        "Main Figure 5 panels    :",
        main_dir,
    )
    print(
        "Supplementary Figure 6 :",
        supplementary_dir,
    )
    print(
        "Exact plotted data      :",
        source_dir,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
