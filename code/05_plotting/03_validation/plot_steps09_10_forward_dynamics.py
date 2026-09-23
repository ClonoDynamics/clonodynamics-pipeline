#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_steps09_10_forward_dynamics.py
=====================================

Publication plotter for ClonoDynamics Step 9–10.

Final contract:
- Step 9 v3 exact equal-weight AB/BA;
- Step 10 v5 compact Step-7 comparison with no pseudo missing-bin
  interpolation in the absolute branch.

The script writes SINGLE-PANEL publication figures only, so the panels can be
assembled later in Illustrator / Inkscape / Affinity / PowerPoint without
changing the analysis.

Main Figure 4
-------------
Figure4A_real_forward
    Genuine longitudinal observed replicate-decoupled forward profile:
    AB, BA and equal-weight cross_combined.

Figure4B_matched_same_vs_cross
    Matched same-measure versus replicate-decoupled forward profile on common4.

Figure4C_absolute_real_vs_pseudo
    Genuine longitudinal versus pseudo technical-null comparison on the shared
    absolute observed-log-frequency support.

Figure4D_percentile_real_vs_pseudo
    Genuine longitudinal versus pseudo technical-null comparison in the
    within-unit transition-level abundance-percentile geometry.

Figure4E_percentile_slope_comparison
    Compact comparison of longitudinal percentile slope (subject-bootstrap
    interval) versus the pseudo configuration-level slope distribution
    (randomization interval).

Supplementary Figure 5
----------------------
SuppFigS5A_subject_slopes
    Subject-specific primary forward slopes.

SuppFigS5B_LOSO_slopes
    Leave-one-subject-out primary forward slopes.

SuppFigS5C_AB_BA_concordance
    AB versus BA bin-wise mean-displacement concordance.

SuppFigS5D_crossover_bootstrap
    Subject-bootstrap distribution of the longitudinal percentile crossover.

Inputs
------
Step 9 directory:
    02_primary_forward_by_bin.csv
    03_matched_same_vs_cross_by_bin.csv
    04_forward_global_summary.csv
    05_subject_forward_by_bin.csv
    06_loso_forward_by_bin.csv
    07_AB_BA_concordance.csv
    09_bootstrap_linear_descriptors.csv

Step 10 directory:
    01_absolute_abundance/
        02_shared_grid.csv
        03_pseudo_null_metrics.csv
        04_empirical_tests.csv
        06_comparison_summary.csv

    02_percentile_geometry/
        03_pseudo_percentile_envelope.csv
        04_pseudo_null_metrics.csv
        05_empirical_tests.csv
        06_pointwise_comparison.csv
        07_longitudinal_bootstrap_linear_descriptors.csv
        09_comparison_summary.csv

Typical run
-----------
python3 plot_steps09_10_forward_dynamics.py \
    --step9-dir ./dataset_longitudinal_results/9-observed_replicate_decoupled_forward_drift \
    --step10-dir ./dataset_longitudinal_results/10-longitudinal_vs_pseudo_forward_null \
    --outdir ./figures/step9_step10_forward

Outputs
-------
<outdir>/
    main_figure/
    supplementary_figure/
    figure_source_data/

All sizes, line widths, colors, axis limits, legend placement and annotation
positions are editable in the dictionaries below.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# =============================================================================
# USER-EDITABLE AESTHETICS
# =============================================================================

FONT_FAMILY = "Arial"
FONT_SIZE = 8
AXIS_LABEL_SIZE = 8
TICK_LABEL_SIZE = 8
LEGEND_FONT_SIZE = 8
ANNOTATION_FONT_SIZE = 8

PDF_FONT_TYPE = 42
PS_FONT_TYPE = 42
PNG_DPI = 600

SAVE_FORMATS = ("pdf", "png")

# ---------------------------------------------------------------------
# Figure dimensions in inches
# ---------------------------------------------------------------------

FIGSIZE_INCHES: Dict[str, Tuple[float, float]] = {
    "Figure4A_real_forward": (3.2, 2.4),
    "Figure4B_matched_same_vs_cross": (3.2, 2.4),
    "Figure4C_absolute_real_vs_pseudo": (3.2, 2.4),
    "Figure4D_percentile_real_vs_pseudo": (3.2, 2.4),
    "Figure4E_percentile_slope_comparison": (3.2, 2.4),
    "SuppFigS5A_subject_slopes": (3.2, 2.4),
    "SuppFigS5B_LOSO_slopes": (3.2, 2.4),
    "SuppFigS5C_AB_BA_concordance": (3.2, 2.4),
    "SuppFigS5D_crossover_bootstrap": (3.2, 2.4),
}

# ---------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------

COLORS = {
    # Step 9 folds
    "AB": "#0072B2",
    "BA": "#D55E00",
    "combined": "#222222",

    # same-versus-cross
    "same": "#D55E00",
    "cross": "#0072B2",

    # real-versus-pseudo
    "longitudinal": "#D55E00",
    "pseudo": "#0072B2",
    "pseudo_fill": "#86A3B4",
    "longitudinal_fill":"#D55E00",

    # robustness panels
    "subject": "#0072B2",
    "loso": "#009E73",
    "bootstrap": "#666666",
    "zero": "#666666",
    "identity": "#777777",
}

# ---------------------------------------------------------------------
# Line styles / markers
# ---------------------------------------------------------------------

LINESTYLES = {
    "AB": "-",
    "BA": "--",
    "combined": "-",
    "same": "-",
    "cross": "-",
    "longitudinal": "-",
    "pseudo": "--",
    "zero": "--",
    "identity": "--",
}

MARKERS = {
    "AB": "o",
    "BA": "s",
    "combined": None,
    "same": None,
    "cross": None,
    "subject": "o",
    "loso": "o",
    "slope_real": "o",
    "slope_pseudo": "s",
}

# ---------------------------------------------------------------------
# Line / marker widths
# ---------------------------------------------------------------------

LINE_WIDTH = 1
PRIMARY_LINE_WIDTH = 1
SECONDARY_LINE_WIDTH = 0.5
ZERO_LINE_WIDTH = 0.5
ERRORBAR_LINE_WIDTH = 0.5
IDENTITY_LINE_WIDTH = 1.0
CAP_SIZE = 2.0
MARKER_SIZE = 2.5
POINT_SIZE = 3.0
SCATTER_SIZE = 28

RIBBON_ALPHA = 0.18
PSEUDO_RIBBON_ALPHA = 0.28
COHORT_CI_ALPHA = 0.20
HIST_ALPHA = 0.75

ZERO_LINE_ALPHA = 0.65
IDENTITY_LINE_ALPHA = 0.70

# ---------------------------------------------------------------------
# Legend placement
# ---------------------------------------------------------------------

LEGEND_LOC = "upper left"
LEGEND_BBOX = (1.02, 1.01)

# ---------------------------------------------------------------------
# Optional titles. Usually off for manuscript assembly.
# ---------------------------------------------------------------------

SHOW_TITLES = False

PLOT_TITLES = {
    "Figure4A_real_forward":
        "Observed replicate-decoupled forward dynamics",
    "Figure4B_matched_same_vs_cross":
        "Matched same- versus cross-replicate forward dynamics",
    "Figure4C_absolute_real_vs_pseudo":
        "Longitudinal versus pseudo technical null: absolute abundance",
    "Figure4D_percentile_real_vs_pseudo":
        "Longitudinal versus pseudo technical null: abundance percentile",
    "Figure4E_percentile_slope_comparison":
        "Percentile-space slope comparison",
    "SuppFigS5A_subject_slopes":
        "Subject-specific forward slopes",
    "SuppFigS5B_LOSO_slopes":
        "Leave-one-subject-out forward slopes",
    "SuppFigS5C_AB_BA_concordance":
        "AB/BA fold concordance",
    "SuppFigS5D_crossover_bootstrap":
        "Bootstrap crossover distribution",
}

# ---------------------------------------------------------------------
# Axis limits
# Leave as None for automatic limits.
# Tuple entries can use None for only one side.
# ---------------------------------------------------------------------

XLIMS = {
    "Figure4A_real_forward": None,
    "Figure4B_matched_same_vs_cross": None,
    "Figure4C_absolute_real_vs_pseudo": None,
    "Figure4D_percentile_real_vs_pseudo": (0.0, 1.0),
    "Figure4E_percentile_slope_comparison": None,
    "SuppFigS5A_subject_slopes": None,
    "SuppFigS5B_LOSO_slopes": None,
    "SuppFigS5C_AB_BA_concordance": None,
    "SuppFigS5D_crossover_bootstrap": (0.0, 1.0),
}

YLIMS = {
    "Figure4A_real_forward": None,
    "Figure4B_matched_same_vs_cross": None,
    "Figure4C_absolute_real_vs_pseudo": None,
    "Figure4D_percentile_real_vs_pseudo": None,
    "Figure4E_percentile_slope_comparison": None,
    "SuppFigS5A_subject_slopes": None,
    "SuppFigS5B_LOSO_slopes": None,
    "SuppFigS5C_AB_BA_concordance": None,
    "SuppFigS5D_crossover_bootstrap": None,
}

# ---------------------------------------------------------------------
# Display switches
# ---------------------------------------------------------------------

SHOW_AB_BA_CI = False
SHOW_COMBINED_CI = True

SHOW_MATCHED_CI = True

SHOW_ABSOLUTE_LONGITUDINAL_CI = True
SHOW_ABSOLUTE_PSEUDO_ENVELOPE = True

SHOW_PERCENTILE_LONGITUDINAL_CI = True
SHOW_PERCENTILE_PSEUDO_ENVELOPE = True

SHOW_ZERO_LINE = True

# ---------------------------------------------------------------------
# Annotation dictionary
# All annotations are optional. Set enabled=False to remove them.
# Coordinates are in axes fraction coordinates.
# ---------------------------------------------------------------------

ANNOTATIONS = {
    "Figure4A_real_forward": {
        "enabled": True,
        "xy_axes": (0.3, 0.74),
        "format": (
            "slope = {slope:.3f}\n"
            "95% bootstrap CI [{lo:.3f}, {hi:.3f}]\n"
            r"$x_{{0}}$ = {xzero:.2f}"
        ),
    },
    "Figure4B_matched_same_vs_cross": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "format": (
            "same slope = {same_slope:.3f}\n"
            "cross slope = {cross_slope:.3f}"
        ),
    },
    "Figure4C_absolute_real_vs_pseudo": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "format": (
            "slope = {slope:.3f}\n"
            r"$P_{{emp}}$(slope) = {p:.5f}" "\n"
            "{nout:d}/{ntot:d} bins outside pseudo 95%"
        ),
    },
    "Figure4D_percentile_real_vs_pseudo": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "format": (
            r"$p_{{cross}}$ = {cross:.3f}" "\n"
            "95% bootstrap CI [{lo:.3f}, {hi:.3f}]\n"
            r"$P_{{emp}}$(slope) = {p:.5f}"
        ),
    },
    "Figure4E_percentile_slope_comparison": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "format": (
            r"$P_{{emp}}$ = {p:.5f}"
        ),
    },
    "SuppFigS5C_AB_BA_concordance": {
        "enabled": True,
        "xy_axes": (0.04, 0.96),
        "format": (
            "Pearson r = {r:.3f}\n"
            "n = {n:d} bins"
        ),
        "va": "top",
    },
    "SuppFigS5D_crossover_bootstrap": {
        "enabled": True,
        "xy_axes": (0.5, 0.95),
        "format": (
            "crossover = {cross:.3f}\n"
            "95% bootstrap CI [{lo:.3f}, {hi:.3f}]"
        ),
        "va": "top",
    },
}

# Number of histogram bins for bootstrap distribution.
HISTOGRAM_BINS = 35


# =============================================================================
# Matplotlib configuration
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


def maybe_zero_line(ax: plt.Axes, orientation: str = "horizontal") -> None:
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
    va = spec.get("va", "bottom")

    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="left",
        va=va,
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
        raise FileNotFoundError(f"{label} not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_final_contracts(
    step9_dir: Path,
    step10_dir: Path,
) -> Dict[str, object]:
    step9 = read_json_required(
        step9_dir / "00_run_config.json",
        "Step 9 run config",
    )
    step10 = read_json_required(
        step10_dir / "00_run_config.json",
        "Step 10 run config",
    )

    s9_est = step9.get("primary_estimand", {})
    if (
        s9_est.get("AB_BA_combination_policy")
        != "exact_equal_fold_weight_after_binning"
        or s9_est.get("AB_BA_row_count_weighting") is not False
    ):
        raise ValueError(
            "Plotter requires final Step 9 exact equal-weight AB/BA output."
        )

    s10_est = step10.get("forward_estimand", {})
    if (
        s10_est.get("AB_BA_combination_policy")
        != "exact_equal_fold_weight_after_binning"
        or s10_est.get("AB_BA_row_count_weighting") is not False
    ):
        raise ValueError(
            "Plotter requires final Step 10 exact equal-weight AB/BA output."
        )

    absolute = step10.get("absolute_abundance", {})
    if absolute.get("pseudo_missing_bin_interpolation") is not False:
        raise ValueError(
            "Plotter requires Step 10 v5 absolute branch with "
            "pseudo_missing_bin_interpolation=false."
        )
    if absolute.get("pseudo_missing_bin_imputation") is not False:
        raise ValueError(
            "Plotter requires Step 10 v5 absolute branch with "
            "pseudo_missing_bin_imputation=false."
        )
    if absolute.get("production_requires_complete_shared_grid") is not True:
        raise ValueError(
            "Plotter requires complete shared-grid production inference."
        )

    pseudo = step10.get("pseudo_null", {})
    if int(pseudo.get("n_configurations_selected", -1)) <= 0:
        raise ValueError("Invalid Step 10 pseudo ensemble size.")

    return {
        "step9_script_version": step9.get("script_version"),
        "step9_analysis_signature": step9.get("analysis_signature"),
        "step10_script_version": step10.get("script_version"),
        "step10_analysis_signature": step10.get("analysis_signature"),
        "n_pseudo_configurations": int(
            pseudo.get("n_configurations_selected")
        ),
        "absolute_min_pseudo_native_bin_coverage": absolute.get(
            "min_pseudo_native_bin_coverage"
        ),
        "absolute_min_configuration_grid_coverage": absolute.get(
            "min_configuration_grid_coverage"
        ),
    }


def read_csv_required(
    path: Path,
    label: str,
) -> pd.DataFrame:
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


def bool_series(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "t", "yes", "y"])
    )


def write_source_data(
    frame: pd.DataFrame,
    source_dir: Path,
    panel: str,
) -> None:
    source_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        source_dir / f"{panel}_source_data.csv",
        index=False,
    )


def linear_fit_xy(
    x: np.ndarray,
    y: np.ndarray,
) -> Tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]

    if len(x) < 2 or np.ptp(x) <= 0:
        return np.nan, np.nan, np.nan

    slope, intercept = np.polyfit(x, y, 1)
    slope = float(slope)
    intercept = float(intercept)
    xzero = (
        float(-intercept / slope)
        if np.isfinite(slope) and abs(slope) > 1e-12
        else np.nan
    )
    return slope, intercept, xzero


# =============================================================================
# Load all inputs
# =============================================================================

def load_inputs(
    step9_dir: Path,
    step10_dir: Path,
) -> Dict[str, object]:

    absolute_dir = step10_dir / "01_absolute_abundance"
    percentile_dir = step10_dir / "02_percentile_geometry"
    contract = validate_final_contracts(
        step9_dir,
        step10_dir,
    )

    return {
        "contract": contract,
        # Step 9
        "step9_primary": read_csv_required(
            step9_dir / "02_primary_forward_by_bin.csv",
            "Step 9 primary forward",
        ),
        "step9_matched": read_csv_required(
            step9_dir / "03_matched_same_vs_cross_by_bin.csv",
            "Step 9 matched same-vs-cross",
        ),
        "step9_global": read_csv_required(
            step9_dir / "04_forward_global_summary.csv",
            "Step 9 global summary",
        ),
        "step9_subject": read_csv_required(
            step9_dir / "05_subject_forward_by_bin.csv",
            "Step 9 subject forward",
        ),
        "step9_loso": read_csv_required(
            step9_dir / "06_loso_forward_by_bin.csv",
            "Step 9 LOSO forward",
        ),
        "step9_abba": read_csv_required(
            step9_dir / "07_AB_BA_concordance.csv",
            "Step 9 AB/BA concordance",
        ),
        "step9_boot": read_csv_required(
            step9_dir / "09_bootstrap_linear_descriptors.csv",
            "Step 9 bootstrap linear descriptors",
        ),

        # Step 10 absolute branch
        "absolute_shared": read_csv_required(
            absolute_dir / "02_shared_grid.csv",
            "Step 10 absolute shared grid",
        ),
        "absolute_null": read_csv_required(
            absolute_dir / "03_pseudo_null_metrics.csv",
            "Step 10 absolute pseudo-null metrics",
        ),
        "absolute_tests": read_csv_required(
            absolute_dir / "04_empirical_tests.csv",
            "Step 10 absolute empirical tests",
        ),
        "absolute_summary": read_csv_required(
            absolute_dir / "06_comparison_summary.csv",
            "Step 10 absolute comparison summary",
        ),
        "absolute_support_audit": read_csv_required(
            absolute_dir / "07_pseudo_native_bin_support.csv",
            "Step 10 absolute native-bin support audit",
        ),

        # Step 10 percentile branch
        "percentile_envelope": read_csv_required(
            percentile_dir / "03_pseudo_percentile_envelope.csv",
            "Step 10 pseudo percentile envelope",
        ),
        "percentile_null": read_csv_required(
            percentile_dir / "04_pseudo_null_metrics.csv",
            "Step 10 percentile pseudo-null metrics",
        ),
        "percentile_tests": read_csv_required(
            percentile_dir / "05_empirical_tests.csv",
            "Step 10 percentile empirical tests",
        ),
        "percentile_pointwise": read_csv_required(
            percentile_dir / "06_pointwise_comparison.csv",
            "Step 10 percentile pointwise comparison",
        ),
        "percentile_boot": read_csv_required(
            percentile_dir / "07_longitudinal_bootstrap_linear_descriptors.csv",
            "Step 10 longitudinal percentile bootstrap descriptors",
        ),
        "percentile_summary": read_csv_required(
            percentile_dir / "09_comparison_summary.csv",
            "Step 10 percentile comparison summary",
        ),
    }


# =============================================================================
# Figure 4A
# =============================================================================

def plot_figure4A(
    primary: pd.DataFrame,
    step9_global: pd.DataFrame,
    step9_boot: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure4A_real_forward"

    require_columns(
        primary,
        [
            "representation",
            "x_center",
            "mean_dx",
            "mean_dx_ci025",
            "mean_dx_ci975",
            "meets_min_n",
        ],
        panel,
    )

    reps = ["cross_AB", "cross_BA", "cross_combined"]
    d = primary[
        primary["representation"].astype(str).isin(reps)
    ].copy()
    d = d[bool_series(d["meets_min_n"])].copy()

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    specs = {
        "cross_AB": (
            "AB",
            COLORS["AB"],
            LINESTYLES["AB"],
            MARKERS["AB"],
            SECONDARY_LINE_WIDTH,
        ),
        "cross_BA": (
            "BA",
            COLORS["BA"],
            LINESTYLES["BA"],
            MARKERS["BA"],
            SECONDARY_LINE_WIDTH,
        ),
        "cross_combined": (
            "AB/BA combined",
            COLORS["combined"],
            LINESTYLES["combined"],
            MARKERS["combined"],
            PRIMARY_LINE_WIDTH,
        ),
    }

    for rep in reps:
        g = d[
            d["representation"].astype(str).eq(rep)
        ].sort_values("x_center")

        if g.empty:
            continue

        label, color, ls, marker, lw = specs[rep]

        x = pd.to_numeric(
            g["x_center"], errors="coerce"
        ).to_numpy(float)
        y = pd.to_numeric(
            g["mean_dx"], errors="coerce"
        ).to_numpy(float)
        lo = pd.to_numeric(
            g["mean_dx_ci025"], errors="coerce"
        ).to_numpy(float)
        hi = pd.to_numeric(
            g["mean_dx_ci975"], errors="coerce"
        ).to_numpy(float)

        show_ci = (
            SHOW_COMBINED_CI
            if rep == "cross_combined"
            else SHOW_AB_BA_CI
        )

        if show_ci:
            ax.fill_between(
                x,
                lo,
                hi,
                color=color,
                alpha=(
                    RIBBON_ALPHA
                    if rep == "cross_combined"
                    else RIBBON_ALPHA * 0.5
                ),
                linewidth=0,
                zorder=1,
            )

        ax.plot(
            x,
            y,
            color=color,
            lw=lw,
            ls=ls,
            marker=marker,
            ms=MARKER_SIZE if marker else 0,
            label=label,
            zorder=3,
        )

    ax.set_xlabel(
        r"Observed initial log-frequency, $x_0$"
    )
    ax.set_ylabel(
        r"Mean observed displacement, $\langle \Delta x \rangle$"
    )

    # Annotation from the point estimate and bootstrap.
    global_row = step9_global[
        step9_global["representation"]
        .astype(str)
        .eq("cross_combined")
    ]
    boot_row = step9_boot[
        step9_boot["representation"]
        .astype(str)
        .eq("cross_combined")
        & step9_boot["metric"].astype(str).eq("slope")
    ]

    xzero_row = step9_global[
        step9_global["representation"]
        .astype(str)
        .eq("cross_combined")
    ]

    if (
        not global_row.empty
        and not boot_row.empty
        and not xzero_row.empty
    ):
        spec = ANNOTATIONS[panel]
        text = spec["format"].format(
            slope=float(
                global_row.iloc[0][
                    "slope_mean_dx_vs_x"
                ]
            ),
            lo=float(
                boot_row.iloc[0][
                    "bootstrap_q025"
                ]
            ),
            hi=float(
                boot_row.iloc[0][
                    "bootstrap_q975"
                ]
            ),
            xzero=float(
                xzero_row.iloc[0]["x_zero"]
            ),
        )
        add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 4B
# =============================================================================

def plot_figure4B(
    matched: pd.DataFrame,
    step9_global: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure4B_matched_same_vs_cross"

    require_columns(
        matched,
        [
            "x_center",
            "mean_dx_same",
            "mean_dx_cross",
            "same_mean_dx_ci025",
            "same_mean_dx_ci975",
            "cross_mean_dx_ci025",
            "cross_mean_dx_ci975",
            "valid_same",
            "valid_cross",
        ],
        panel,
    )

    d = matched[
        bool_series(matched["valid_same"])
        & bool_series(matched["valid_cross"])
    ].copy().sort_values("x_center")

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    x = pd.to_numeric(
        d["x_center"], errors="coerce"
    ).to_numpy(float)

    y_same = pd.to_numeric(
        d["mean_dx_same"], errors="coerce"
    ).to_numpy(float)
    y_cross = pd.to_numeric(
        d["mean_dx_cross"], errors="coerce"
    ).to_numpy(float)

    if SHOW_MATCHED_CI:
        same_lo = pd.to_numeric(
            d["same_mean_dx_ci025"],
            errors="coerce",
        ).to_numpy(float)
        same_hi = pd.to_numeric(
            d["same_mean_dx_ci975"],
            errors="coerce",
        ).to_numpy(float)
        cross_lo = pd.to_numeric(
            d["cross_mean_dx_ci025"],
            errors="coerce",
        ).to_numpy(float)
        cross_hi = pd.to_numeric(
            d["cross_mean_dx_ci975"],
            errors="coerce",
        ).to_numpy(float)

        ax.fill_between(
            x,
            same_lo,
            same_hi,
            color=COLORS["same"],
            alpha=RIBBON_ALPHA,
            linewidth=0,
        )
        ax.fill_between(
            x,
            cross_lo,
            cross_hi,
            color=COLORS["cross"],
            alpha=RIBBON_ALPHA,
            linewidth=0,
        )

    ax.plot(
        x,
        y_same,
        color=COLORS["same"],
        lw=PRIMARY_LINE_WIDTH,
        ls=LINESTYLES["same"],
        label="Same-measure",
    )
    ax.plot(
        x,
        y_cross,
        color=COLORS["cross"],
        lw=PRIMARY_LINE_WIDTH,
        ls=LINESTYLES["cross"],
        label="Replicate-decoupled",
    )

    ax.set_xlabel(
        r"Observed initial log-frequency, $x_0$"
    )
    ax.set_ylabel(
        r"Mean observed displacement, $\langle \Delta x \rangle$"
    )

    same_row = step9_global[
        step9_global["representation"]
        .astype(str)
        .eq("same_combined")
    ]
    cross_row = step9_global[
        step9_global["representation"]
        .astype(str)
        .eq("cross_combined_matched")
    ]

    if not same_row.empty and not cross_row.empty:
        spec = ANNOTATIONS[panel]
        text = spec["format"].format(
            same_slope=float(
                same_row.iloc[0][
                    "slope_mean_dx_vs_x"
                ]
            ),
            cross_slope=float(
                cross_row.iloc[0][
                    "slope_mean_dx_vs_x"
                ]
            ),
        )
        add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 4C
# =============================================================================

def plot_figure4C(
    shared: pd.DataFrame,
    tests: pd.DataFrame,
    summary: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure4C_absolute_real_vs_pseudo"

    require_columns(
        shared,
        [
            "x",
            "longitudinal_mean_dx",
            "longitudinal_ci025",
            "longitudinal_ci975",
            "pseudo_median_dx",
            "pseudo_q025",
            "pseudo_q975",
        ],
        panel,
    )

    d = shared.copy().sort_values("x")
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    x = pd.to_numeric(
        d["x"], errors="coerce"
    ).to_numpy(float)

    real = pd.to_numeric(
        d["longitudinal_mean_dx"],
        errors="coerce",
    ).to_numpy(float)
    real_lo = pd.to_numeric(
        d["longitudinal_ci025"],
        errors="coerce",
    ).to_numpy(float)
    real_hi = pd.to_numeric(
        d["longitudinal_ci975"],
        errors="coerce",
    ).to_numpy(float)

    pseudo = pd.to_numeric(
        d["pseudo_median_dx"],
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

    if SHOW_ABSOLUTE_PSEUDO_ENVELOPE:
        ax.fill_between(
            x,
            pseudo_lo,
            pseudo_hi,
            color=COLORS["pseudo_fill"],
            alpha=PSEUDO_RIBBON_ALPHA,
            linewidth=0,
            label="Pseudo 2.5–97.5% randomization interval",
        )

    if SHOW_ABSOLUTE_LONGITUDINAL_CI:
        ax.fill_between(
            x,
            real_lo,
            real_hi,
            color=COLORS["longitudinal_fill"],
            alpha=COHORT_CI_ALPHA,
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
        r"Observed initial log-frequency, $x_0$"
    )
    ax.set_ylabel(
        r"Mean observed displacement, $\langle \Delta x \rangle$"
    )

    p_row = tests[
        tests["metric"].astype(str).eq("slope")
    ]
    if not p_row.empty and not summary.empty:
        row = summary.iloc[0]
        spec = ANNOTATIONS[panel]
        text = spec["format"].format(
            slope=float(row["longitudinal_slope"]),
            p=float(p_row.iloc[0]["empirical_p"]),
            nout=int(row["n_outside_pseudo_95"]),
            ntot=int(row["n_shared_grid_points"]),
        )
        add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 4D
# =============================================================================

def plot_figure4D(
    pointwise: pd.DataFrame,
    tests: pd.DataFrame,
    summary: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure4D_percentile_real_vs_pseudo"

    require_columns(
        pointwise,
        [
            "percentile_center",
            "mean_dx",
            "mean_dx_ci025",
            "mean_dx_ci975",
            "pseudo_median_dx",
            "pseudo_q025",
            "pseudo_q975",
        ],
        panel,
    )

    d = pointwise.copy().sort_values("percentile_center")
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    x = pd.to_numeric(
        d["percentile_center"],
        errors="coerce",
    ).to_numpy(float)

    real = pd.to_numeric(
        d["mean_dx"],
        errors="coerce",
    ).to_numpy(float)
    real_lo = pd.to_numeric(
        d["mean_dx_ci025"],
        errors="coerce",
    ).to_numpy(float)
    real_hi = pd.to_numeric(
        d["mean_dx_ci975"],
        errors="coerce",
    ).to_numpy(float)

    pseudo = pd.to_numeric(
        d["pseudo_median_dx"],
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

    if SHOW_PERCENTILE_PSEUDO_ENVELOPE:
        ax.fill_between(
            x,
            pseudo_lo,
            pseudo_hi,
            color=COLORS["pseudo_fill"],
            alpha=PSEUDO_RIBBON_ALPHA,
            linewidth=0,
            label="Pseudo 2.5–97.5% randomization interval",
        )

    if SHOW_PERCENTILE_LONGITUDINAL_CI:
        ax.fill_between(
            x,
            real_lo,
            real_hi,
            color=COLORS["longitudinal_fill"],
            alpha=COHORT_CI_ALPHA,
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
        marker="o",
        ms=MARKER_SIZE,
        label="Longitudinal",
        zorder=4,
    )

    ax.set_xlabel(
        "Within-unit initial-abundance percentile"
    )
    ax.set_ylabel(
        r"Mean observed displacement, $\langle \Delta x \rangle$"
    )

    p_row = tests[
        tests["metric"].astype(str).eq("slope")
    ]
    if not summary.empty and not p_row.empty:
        row = summary.iloc[0]
        spec = ANNOTATIONS[panel]
        text = spec["format"].format(
            cross=float(
                row[
                    "longitudinal_crossover_percentile"
                ]
            ),
            lo=float(
                row[
                    "longitudinal_crossover_bootstrap_q025"
                ]
            ),
            hi=float(
                row[
                    "longitudinal_crossover_bootstrap_q975"
                ]
            ),
            p=float(p_row.iloc[0]["empirical_p"]),
        )
        add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 4E
# =============================================================================

def plot_figure4E(
    pseudo_null: pd.DataFrame,
    tests: pd.DataFrame,
    summary: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure4E_percentile_slope_comparison"

    require_columns(
        pseudo_null,
        ["slope"],
        panel,
    )

    slopes = pd.to_numeric(
        pseudo_null["slope"],
        errors="coerce",
    ).to_numpy(float)
    slopes = slopes[np.isfinite(slopes)]

    if slopes.size == 0:
        raise ValueError(
            f"{panel}: no finite pseudo slopes"
        )

    if summary.empty:
        raise ValueError(
            f"{panel}: percentile summary is empty"
        )

    row = summary.iloc[0]

    long_slope = float(
        row["longitudinal_slope"]
    )
    long_lo = float(
        row["longitudinal_slope_bootstrap_q025"]
    )
    long_hi = float(
        row["longitudinal_slope_bootstrap_q975"]
    )

    pseudo_med = float(np.median(slopes))
    pseudo_lo = float(
        np.quantile(slopes, 0.025)
    )
    pseudo_hi = float(
        np.quantile(slopes, 0.975)
    )

    source = pd.DataFrame(
        [
            {
                "group": "Longitudinal",
                "center": long_slope,
                "lower": long_lo,
                "upper": long_hi,
                "interval_type": "95% subject-bootstrap CI",
                "n": np.nan,
            },
            {
                "group": "Pseudo technical null",
                "center": pseudo_med,
                "lower": pseudo_lo,
                "upper": pseudo_hi,
                "interval_type": "2.5–97.5% randomization interval",
                "n": int(slopes.size),
            },
        ]
    )
    write_source_data(source, source_dir, panel)

    fig, ax = make_figure(panel)
    maybe_zero_line(ax, orientation="vertical")

    y_positions = [1, 0]

    specs = [
        (
            "Longitudinal",
            long_slope,
            long_lo,
            long_hi,
            COLORS["longitudinal"],
            MARKERS["slope_real"],
            y_positions[0],
        ),
        (
            "Pseudo technical null",
            pseudo_med,
            pseudo_lo,
            pseudo_hi,
            COLORS["pseudo"],
            MARKERS["slope_pseudo"],
            y_positions[1],
        ),
    ]

    for (
        label,
        center,
        lo,
        hi,
        color,
        marker,
        y,
    ) in specs:
        ax.errorbar(
            center,
            y,
            xerr=np.asarray(
                [[center - lo], [hi - center]]
            ),
            fmt=marker,
            color=color,
            ms=POINT_SIZE,
            lw=ERRORBAR_LINE_WIDTH,
            capsize=CAP_SIZE,
            zorder=3,
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [
            "Longitudinal",
            "Pseudo technical null",
        ]
    )
    ax.set_xlabel(
        "Slope of mean displacement vs abundance percentile"
    )
    ax.set_ylabel("")

    p_row = tests[
        tests["metric"].astype(str).eq("slope")
    ]
    if not p_row.empty:
        spec = ANNOTATIONS[panel]
        text = spec["format"].format(
            p=float(
                p_row.iloc[0]["empirical_p"]
            )
        )
        add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S5A
# =============================================================================

def slope_by_group(
    frame: pd.DataFrame,
    group_col: str,
    x_col: str = "x_center",
    y_col: str = "mean_dx",
) -> pd.DataFrame:
    rows = []

    for group, g in frame.groupby(
        group_col,
        sort=True,
    ):
        x = pd.to_numeric(
            g[x_col],
            errors="coerce",
        ).to_numpy(float)
        y = pd.to_numeric(
            g[y_col],
            errors="coerce",
        ).to_numpy(float)

        slope, intercept, xzero = linear_fit_xy(
            x,
            y,
        )

        rows.append(
            {
                group_col: group,
                "n_bins": int(
                    np.sum(
                        np.isfinite(x)
                        & np.isfinite(y)
                    )
                ),
                "slope": slope,
                "intercept": intercept,
                "x_zero": xzero,
            }
        )

    return pd.DataFrame(rows)


def plot_supp_s5A(
    subject_frame: pd.DataFrame,
    step9_global: pd.DataFrame,
    step9_boot: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS5A_subject_slopes"

    require_columns(
        subject_frame,
        ["subject", "x_center", "mean_dx"],
        panel,
    )

    slopes = slope_by_group(
        subject_frame,
        "subject",
    )
    write_source_data(slopes, source_dir, panel)

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    subjects = slopes["subject"].astype(str).tolist()
    x = np.arange(len(subjects))
    y = slopes["slope"].to_numpy(float)

    ax.scatter(
        x,
        y,
        s=SCATTER_SIZE,
        color=COLORS["subject"],
        zorder=3,
    )

    global_row = step9_global[
        step9_global["representation"]
        .astype(str)
        .eq("cross_combined")
    ]
    boot_row = step9_boot[
        step9_boot["representation"]
        .astype(str)
        .eq("cross_combined")
        & step9_boot["metric"].astype(str).eq("slope")
    ]

    if not global_row.empty:
        cohort = float(
            global_row.iloc[0]["slope_mean_dx_vs_x"]
        )
        ax.axhline(
            cohort,
            color=COLORS["combined"],
            lw=PRIMARY_LINE_WIDTH,
            ls="-",
            label="Cohort",
            zorder=2,
        )

    if not boot_row.empty:
        lo = float(
            boot_row.iloc[0]["bootstrap_q025"]
        )
        hi = float(
            boot_row.iloc[0]["bootstrap_q975"]
        )
        ax.axhspan(
            lo,
            hi,
            color=COLORS["combined"],
            alpha=COHORT_CI_ALPHA,
            linewidth=0,
            label="Cohort 95% bootstrap CI",
            zorder=1,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(
        subjects,
        rotation=0,
    )
    ax.set_xlabel("Subject")
    ax.set_ylabel(
        "Slope of mean displacement vs\nobserved log-frequency"
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S5B
# =============================================================================

def plot_supp_s5B(
    loso_frame: pd.DataFrame,
    step9_global: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS5B_LOSO_slopes"

    require_columns(
        loso_frame,
        [
            "omitted_subject",
            "x_center",
            "mean_dx",
        ],
        panel,
    )

    slopes = slope_by_group(
        loso_frame,
        "omitted_subject",
    )
    write_source_data(slopes, source_dir, panel)

    fig, ax = make_figure(panel)
    maybe_zero_line(ax)

    labels = [
        str(x)
        for x in slopes[
            "omitted_subject"
        ].tolist()
    ]
    x = np.arange(len(labels))
    y = slopes["slope"].to_numpy(float)

    ax.scatter(
        x,
        y,
        s=SCATTER_SIZE,
        color=COLORS["loso"],
        zorder=3,
    )

    global_row = step9_global[
        step9_global["representation"]
        .astype(str)
        .eq("cross_combined")
    ]
    if not global_row.empty:
        cohort = float(
            global_row.iloc[0]["slope_mean_dx_vs_x"]
        )
        ax.axhline(
            cohort,
            color=COLORS["combined"],
            lw=PRIMARY_LINE_WIDTH,
            ls="-",
            label="Full cohort",
            zorder=2,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Omitted subject")
    ax.set_ylabel(
        "LOSO slope of mean displacement vs\nobserved log-frequency"
    )

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S5C
# =============================================================================

def plot_supp_s5C(
    primary: pd.DataFrame,
    concordance: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS5C_AB_BA_concordance"

    ab = primary[
        primary["representation"]
        .astype(str)
        .eq("cross_AB")
    ][
        ["bin", "mean_dx", "meets_min_n"]
    ].rename(
        columns={
            "mean_dx": "mean_dx_AB",
            "meets_min_n": "valid_AB",
        }
    )

    ba = primary[
        primary["representation"]
        .astype(str)
        .eq("cross_BA")
    ][
        ["bin", "mean_dx", "meets_min_n"]
    ].rename(
        columns={
            "mean_dx": "mean_dx_BA",
            "meets_min_n": "valid_BA",
        }
    )

    d = ab.merge(
        ba,
        on="bin",
        how="inner",
    )
    d = d[
        bool_series(d["valid_AB"])
        & bool_series(d["valid_BA"])
    ].copy()

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)

    x = pd.to_numeric(
        d["mean_dx_AB"],
        errors="coerce",
    ).to_numpy(float)
    y = pd.to_numeric(
        d["mean_dx_BA"],
        errors="coerce",
    ).to_numpy(float)

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]

    if len(x) == 0:
        raise ValueError(
            f"{panel}: no matched AB/BA bins"
        )

    lim_min = float(
        min(np.min(x), np.min(y))
    )
    lim_max = float(
        max(np.max(x), np.max(y))
    )
    pad = 0.05 * max(
        lim_max - lim_min,
        1e-6,
    )
    lim_min -= pad
    lim_max += pad

    ax.plot(
        [lim_min, lim_max],
        [lim_min, lim_max],
        color=COLORS["identity"],
        lw=IDENTITY_LINE_WIDTH,
        ls=LINESTYLES["identity"],
        alpha=IDENTITY_LINE_ALPHA,
        zorder=1,
    )

    ax.scatter(
        x,
        y,
        s=SCATTER_SIZE,
        color=COLORS["combined"],
        zorder=3,
    )

    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)

    ax.set_xlabel(
        r"AB mean displacement, $\langle \Delta x \rangle$"
    )
    ax.set_ylabel(
        r"BA mean displacement, $\langle \Delta x \rangle$"
    )

    if not concordance.empty:
        row = concordance.iloc[0]
        spec = ANNOTATIONS[panel]
        text = spec["format"].format(
            r=float(row["pearson_r_AB_BA"]),
            n=int(row["n_matched_bins"]),
        )
        add_annotation(ax, panel, text)

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S5D
# =============================================================================

def plot_supp_s5D(
    percentile_boot: pd.DataFrame,
    percentile_summary: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS5D_crossover_bootstrap"

    require_columns(
        percentile_boot,
        [
            "bootstrap",
            "crossover_percentile",
        ],
        panel,
    )

    vals = pd.to_numeric(
        percentile_boot["crossover_percentile"],
        errors="coerce",
    ).to_numpy(float)
    vals = vals[
        np.isfinite(vals)
        & (vals >= 0)
        & (vals <= 1)
    ]

    if vals.size == 0:
        raise ValueError(
            f"{panel}: no finite crossover values in [0,1]"
        )

    source = pd.DataFrame(
        {
            "crossover_percentile": vals
        }
    )
    write_source_data(source, source_dir, panel)

    fig, ax = make_figure(panel)

    ax.hist(
        vals,
        bins=HISTOGRAM_BINS,
        density=True,
        color=COLORS["bootstrap"],
        alpha=HIST_ALPHA,
        edgecolor="none",
    )

    if not percentile_summary.empty:
        row = percentile_summary.iloc[0]
        cross = float(
            row["longitudinal_crossover_percentile"]
        )
        lo = float(
            row[
                "longitudinal_crossover_bootstrap_q025"
            ]
        )
        hi = float(
            row[
                "longitudinal_crossover_bootstrap_q975"
            ]
        )

        ax.axvspan(
            lo,
            hi,
            color=COLORS["combined"],
            alpha=COHORT_CI_ALPHA,
            linewidth=0,
        )
        ax.axvline(
            cross,
            color=COLORS["combined"],
            lw=PRIMARY_LINE_WIDTH,
        )

        spec = ANNOTATIONS[panel]
        text = spec["format"].format(
            cross=cross,
            lo=lo,
            hi=hi,
        )
        add_annotation(ax, panel, text)

    ax.set_xlabel(
        "Longitudinal crossover abundance percentile"
    )
    ax.set_ylabel("Bootstrap density")

    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


def write_plot_contract_audit(
    outdir: Path,
    contract: Mapping[str, object],
    absolute_support: pd.DataFrame,
    absolute_shared: pd.DataFrame,
) -> None:
    require_columns(
        absolute_support,
        [
            "bin",
            "x_center",
            "valid_coverage_fraction",
            "used_in_shared_absolute_analysis",
        ],
        "absolute native-bin support audit",
    )

    used = absolute_support[
        bool_series(
            absolute_support["used_in_shared_absolute_analysis"]
        )
    ].copy()

    audit = {
        **contract,
        "absolute_shared_bins": int(len(absolute_shared)),
        "absolute_support_bins_marked_used": int(len(used)),
        "absolute_used_min_valid_coverage": (
            float(
                pd.to_numeric(
                    used["valid_coverage_fraction"],
                    errors="coerce",
                ).min()
            )
            if not used.empty else np.nan
        ),
        "absolute_shared_pseudo_n_min": (
            int(
                pd.to_numeric(
                    absolute_shared["pseudo_n_configurations"],
                    errors="coerce",
                ).min()
            )
            if "pseudo_n_configurations" in absolute_shared.columns
            else np.nan
        ),
        "absolute_shared_pseudo_n_max": (
            int(
                pd.to_numeric(
                    absolute_shared["pseudo_n_configurations"],
                    errors="coerce",
                ).max()
            )
            if "pseudo_n_configurations" in absolute_shared.columns
            else np.nan
        ),
    }

    pd.DataFrame([audit]).to_csv(
        outdir / "00_plot_contract_audit.csv",
        index=False,
    )


# =============================================================================
# CLI / main
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Create single-panel publication figures for ClonoDynamics "
            "Step 9–10 forward dynamics."
        ),
    )

    parser.add_argument(
        "--step9-dir",
        required=True,
        type=Path,
        help="Final Step-9 output directory.",
    )
    parser.add_argument(
        "--step10-dir",
        required=True,
        type=Path,
        help="Final Step-10 output directory.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for publication plots.",
    )

    return parser


def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(argv)

    configure_matplotlib()

    step9_dir = (
        args.step9_dir
        .expanduser()
        .resolve(strict=True)
    )
    step10_dir = (
        args.step10_dir
        .expanduser()
        .resolve(strict=True)
    )
    outdir = (
        args.outdir
        .expanduser()
        .resolve()
    )

    main_dir = outdir / "main_figure"
    supplementary_dir = (
        outdir / "supplementary_figure"
    )
    source_dir = (
        outdir / "figure_source_data"
    )

    data = load_inputs(
        step9_dir,
        step10_dir,
    )

    # Main Figure 4
    plot_figure4A(
        data["step9_primary"],
        data["step9_global"],
        data["step9_boot"],
        main_dir,
        source_dir,
    )

    plot_figure4B(
        data["step9_matched"],
        data["step9_global"],
        main_dir,
        source_dir,
    )

    plot_figure4C(
        data["absolute_shared"],
        data["absolute_tests"],
        data["absolute_summary"],
        main_dir,
        source_dir,
    )

    plot_figure4D(
        data["percentile_pointwise"],
        data["percentile_tests"],
        data["percentile_summary"],
        main_dir,
        source_dir,
    )

    plot_figure4E(
        data["percentile_null"],
        data["percentile_tests"],
        data["percentile_summary"],
        main_dir,
        source_dir,
    )

    # Supplementary Figure 5
    plot_supp_s5A(
        data["step9_subject"],
        data["step9_global"],
        data["step9_boot"],
        supplementary_dir,
        source_dir,
    )

    plot_supp_s5B(
        data["step9_loso"],
        data["step9_global"],
        supplementary_dir,
        source_dir,
    )

    plot_supp_s5C(
        data["step9_primary"],
        data["step9_abba"],
        supplementary_dir,
        source_dir,
    )

    plot_supp_s5D(
        data["percentile_boot"],
        data["percentile_summary"],
        supplementary_dir,
        source_dir,
    )

    write_plot_contract_audit(
        outdir,
        data["contract"],
        data["absolute_support_audit"],
        data["absolute_shared"],
    )

    print("\n[DONE] Step 9–10 publication plots")
    print("Main Figure 4 panels    :", main_dir)
    print("Supplementary Figure 5 :", supplementary_dir)
    print("Exact plotted data      :", source_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
