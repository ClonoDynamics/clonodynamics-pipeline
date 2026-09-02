#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
8-plot_pseudo_longitudinal_representation_assessment.py
=======================================================

Plot-only companion to ClonoDynamics Step 8:

    8-pseudo_longitudinal_representation_assessment.py

Scientific role
---------------
Step 8 evaluates alternative finite-time displacement representations under the
pseudo-longitudinal technical reference while keeping transition conditioning
fixed at `xstar_latent`.

The assessment deliberately separates:

    1. directional-null adequacy;
    2. residual sequencing-depth dependence;
    3. closure / count behavior;
    4. descriptive fluctuation magnitude.

These properties are not collapsed into a single score and no automatic
representation winner is assigned.

This script is plotting-only. It visualizes Step-8 outputs and does not read the
Step-5 transition table, redefine xstar bins, recompute correlations or
bootstrap effects, rerun sensitivity scenarios, or rank representations.

Position in the pipeline
------------------------
Pseudo-longitudinal branch:

    Step 5  latent transitions
        |
        +--> Step 7  conditioning / replicate-decoupling / ordering benchmark
        |
        +--> Step 8  xstar-conditioned displacement-representation assessment
                         |
                         +--> Step 12 longitudinal-vs-pseudo fluctuation baseline

Step 8 does not computationally consume Step-7 outputs. Both Steps use the
pseudo-longitudinal data for different methodological questions.

Analysis input
--------------
Required:

    --analysis-dir

pointing to the output directory produced by:

    8-pseudo_longitudinal_representation_assessment.py

Required Step-8 tables
----------------------
    03_pooled_representation_curves.csv
        Per-representation, per-xstar-bin displacement summaries.

    08_depth_dependence_summary.csv
        Precomputed residual sequencing-depth dependence by representation.

    09_closure_diagnostic_curves.csv
        Precomputed closure-component summaries across xstar bins.

    10_closure_depth_summary.csv
        Precomputed whole-dataset closure-versus-depth association.

Optional Step-8 table
---------------------
    12_sensitivity_representation_assessment.csv
        Representation descriptors across primary and alternative abundance
        ranges. Sensitivity panels are skipped when this file is absent.

Other Step-8 outputs remain analysis-layer products and are not recomputed by
this plotter, including:

    04_pooled_representation_assessment.csv
    05_subject_representation_assessment.csv
    06_subject_bootstrap_pairwise_effects.csv
    07_pairwise_displacement_agreement.csv
    11_representation_assessment.csv
    13_sensitivity_bootstrap_pairwise_effects.csv
    representation_assessment_report.md

Generated figures
-----------------
Main composite:

    00_main_representation_benchmark
        A  P(dx > 0) - 0.5 versus xstar
        B  mean squared displacement versus xstar
        C  |Spearman(dx, delta log depth)| by representation

Supplementary composite:

    Sx_representation_diagnostics
        A  mean displacement versus xstar
        B  MAD(dx) versus xstar
        C  mean closure component versus xstar, overlaid with
           mean(observed-frequency dx) - mean(log-count dx)
        D  precomputed closure-depth association

Standalone figures:

    01_primary_mean_dx
    02_primary_median_dx
    03_primary_ppos_bias
    04_depth_change_dependence
    05_primary_msd
    06_primary_mad_dx
    07_closure_median
    08_closure_mad
    09_sensitivity_depth_dependence
    10_sensitivity_directional_mean_bias
    11_sensitivity_directional_median_bias
    12_sensitivity_directional_ppos_bias
    13_closure_mean_identity
    14_closure_depth_association

Non-figure outputs
------------------
Only graphical provenance is written:

    figure_manifest.csv
    plot_config.json

No new scientific summary table is created by the plotting layer.

Display-only transformations
----------------------------
The plotter may calculate only quantities needed to render an existing Step-8
comparison, including:

    P(dx > 0) - 0.5

and:

    mean(observed-frequency dx) - mean(log-count dx)

for the closure identity overlay. These are not exported as new analysis
results.

Figure customization
--------------------
Edit:

    FIGURE_SIZES_INCHES
        Per-figure (width, height) in inches.

    AXIS_CONFIG
        Per-panel:
            xlim
            ylim
            xscale
            yscale
            xtick_step
            ytick_step

Composite keys such as `MainA` and `SxC` may be modified independently from
standalone panels.

Example
-------
python3 ./code/plotting/8-plot_pseudo_longitudinal_representation_assessment.py \
    --analysis-dir \
    ./dataset_pseudo_clonodynamics_results/pseudo-longitudinal/8-pseudo_longitudinal_representation_assessment \
    --figures-dir ./figures/8-pseudo-longitudinal-representation-assessment \
    --formats pdf,png \
    --figure-set all
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


SCRIPT_VERSION = "v1-step8-public-plotter-2026-08-25"


# =============================================================================
# Figure dimensions and axis configuration
# =============================================================================

FIGURE_SIZES_INCHES = {
    "00_main_representation_benchmark": (9.0, 2.5),
    "Sx_representation_diagnostics": (7.2, 6.0),
    "01_primary_mean_dx": (3.2, 2.4),
    "02_primary_median_dx": (3.2, 2.4),
    "03_primary_ppos_bias": (3.2, 2.4),
    "04_depth_change_dependence": (3.2, 2.4),
    "05_primary_msd": (3.2, 2.4),
    "06_primary_mad_dx": (3.2, 2.4),
    "07_closure_median": (3.2, 2.4),
    "08_closure_mad": (3.2, 2.4),
    "09_sensitivity_depth_dependence": (3.2, 2.4),
    "10_sensitivity_directional_mean_bias": (3.2, 2.4),
    "11_sensitivity_directional_median_bias": (3.2, 2.4),
    "12_sensitivity_directional_ppos_bias": (3.2, 2.4),
    "13_closure_mean_identity": (3.2, 2.4),
    "14_closure_depth_association": (3.2, 2.4),
}


def _axis_cfg(
    xlim=None,
    ylim=None,
    xscale=None,
    yscale=None,
    xtick_step=None,
    ytick_step=None,
) -> Dict[str, object]:
    return {
        "xlim": xlim,
        "ylim": ylim,
        "xscale": xscale,
        "yscale": yscale,
        "xtick_step": xtick_step,
        "ytick_step": ytick_step,
    }


# Edit this dictionary to tune axes without touching plotting functions.
# Composite-panel keys (MainA-C, SxA-D) override the corresponding standalone
# settings only inside the composite figure.
AXIS_CONFIG = {
    "MainA": _axis_cfg(ylim=(-0.20, 0.20), ytick_step=0.05),
    "MainB": _axis_cfg(ylim=(0.0, 0.40), ytick_step=0.05),
    "MainC": _axis_cfg(ylim=(0.0, 0.25), ytick_step=0.05),
    "SxA": _axis_cfg(ylim=(-0.20, 0.20), ytick_step=0.05),
    "SxB": _axis_cfg(ylim=(0.0, 0.50)),
    "SxC": _axis_cfg(ylim=(-0.20, 0.20), ytick_step=0.05),
    "SxD": _axis_cfg(xlim=(-0.60, 0.05), xtick_step=0.10),
    "00_main_representation_benchmark": _axis_cfg(),
    "Sx_representation_diagnostics": _axis_cfg(),
    "01_primary_mean_dx": _axis_cfg(ylim=(-0.20, 0.20), ytick_step=0.05),
    "02_primary_median_dx": _axis_cfg(ylim=(-0.20, 0.20), ytick_step=0.05),
    "03_primary_ppos_bias": _axis_cfg(ylim=(-0.20, 0.20), ytick_step=0.05),
    "04_depth_change_dependence": _axis_cfg(ylim=(0.0, 0.25), ytick_step=0.05),
    "05_primary_msd": _axis_cfg(ylim=(0.0, 0.40), ytick_step=0.05),
    "06_primary_mad_dx": _axis_cfg(ylim=(0.0, 0.50)),
    "07_closure_median": _axis_cfg(ylim=(-0.20, 0.20)),
    "08_closure_mad": _axis_cfg(ylim=(0.0, 0.24)),
    "09_sensitivity_depth_dependence": _axis_cfg(ylim=(0.0, 0.25), ytick_step=0.05),
    "10_sensitivity_directional_mean_bias": _axis_cfg(),
    "11_sensitivity_directional_median_bias": _axis_cfg(),
    "12_sensitivity_directional_ppos_bias": _axis_cfg(),
    "13_closure_mean_identity": _axis_cfg(ylim=(-0.20, 0.20), ytick_step=0.05),
    "14_closure_depth_association": _axis_cfg(xlim=(-0.60, 0.05), xtick_step=0.10),
}


REPRESENTATION_ORDER = ["latent_frequency", "observed_frequency", "log_count"]
REPRESENTATION_LABELS = {
    "latent_frequency": "Latent frequency",
    "observed_frequency": "Observed frequency",
    "log_count": "Log-count",
}

SCENARIO_ORDER = ["primary", "low_boundary_trim", "central_range", "reliable_range"]
SCENARIO_LABELS = {
    "primary": "Primary",
    "low_boundary_trim": "Low-boundary\ntrim",
    "central_range": "Central\nrange",
    "reliable_range": "Reliable\nrange",
}

# Files required to render the Step-8 figure set.
REQUIRED_FILES = {
    "curves": "03_pooled_representation_curves.csv",
    "depth": "08_depth_dependence_summary.csv",
    "closure": "09_closure_diagnostic_curves.csv",
    "closure_depth": "10_closure_depth_summary.csv",
}
# Optional Step-8 sensitivity input; corresponding panels are skipped when absent.
OPTIONAL_FILES = {
    "sensitivity": "12_sensitivity_representation_assessment.csv",
}


# =============================================================================
# Style helpers
# =============================================================================


def setup_style(font_size: float) -> None:
    try:
        from matplotlib import font_manager
        font_manager.findfont("Arial", fallback_to_default=False)
        font_family = "Arial"
    except Exception:
        font_family = "DejaVu Sans"

    plt.rcParams.update(
        {
            "font.family": font_family,
            "font.size": font_size,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size,
            "xtick.labelsize": max(font_size - 1.0, 5.0),
            "ytick.labelsize": max(font_size - 1.0, 5.0),
            "legend.fontsize": max(font_size - 1.0, 5.0),
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")


def candidate_color_map() -> Dict[str, str]:
    colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if len(colors) < len(REPRESENTATION_ORDER):
        colors = ["C{}".format(i) for i in range(len(REPRESENTATION_ORDER))]
    return {
        representation: colors[index]
        for index, representation in enumerate(REPRESENTATION_ORDER)
    }


def axis_legend(ax, title: str) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles,
            labels,
            title=title,
            frameon=False,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
        )


def parse_limits(text: str) -> Optional[Tuple[Optional[float], Optional[float]]]:
    if text is None or not str(text).strip():
        return None
    parts = [part.strip() for part in str(text).split(",")]
    if len(parts) != 2:
        raise ValueError("Axis limits must be specified as lo,hi")

    def parse_value(value: str) -> Optional[float]:
        if value.lower() in {"", "none", "nan"}:
            return None
        return float(value)

    return parse_value(parts[0]), parse_value(parts[1])


def apply_axis_settings(
    ax,
    key: str,
    fallback_xlim: str = "",
    fallback_ylim: str = "",
) -> None:
    settings = AXIS_CONFIG.get(key, {})
    xlim = settings.get("xlim")
    ylim = settings.get("ylim")

    if settings.get("xscale"):
        ax.set_xscale(str(settings["xscale"]))
    if settings.get("yscale"):
        ax.set_yscale(str(settings["yscale"]))

    if xlim is None:
        xlim = parse_limits(fallback_xlim)
    if ylim is None:
        ylim = parse_limits(fallback_ylim)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if settings.get("xtick_step") is not None:
        ax.xaxis.set_major_locator(MultipleLocator(float(settings["xtick_step"])))
    if settings.get("ytick_step") is not None:
        ax.yaxis.set_major_locator(MultipleLocator(float(settings["ytick_step"])))


def add_panel_label(ax, label: str, fontsize: float = 10.0) -> None:
    ax.text(
        -0.16,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight="bold",
        va="top",
        ha="left",
    )


def save_figure(fig, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def parse_formats(
    formats_text: str,
    legacy_fmt: Optional[str],
    legacy_also_png: bool,
) -> List[str]:
    """Resolve canonical --formats while retaining inexpensive CLI compatibility."""
    if legacy_fmt is not None:
        values = [legacy_fmt]
        if legacy_also_png and legacy_fmt != "png":
            values.append("png")
    else:
        values = [
            value.strip().lower()
            for value in str(formats_text).split(",")
            if value.strip()
        ]

    allowed = {"pdf", "png", "svg"}
    invalid = [value for value in values if value not in allowed]
    if invalid:
        raise ValueError(
            "Unsupported figure format(s): {}. Use pdf,png,svg.".format(invalid)
        )

    out: List[str] = []
    for value in values:
        if value not in out:
            out.append(value)
    if not out:
        raise ValueError("At least one figure format is required.")
    return out


# =============================================================================
# I/O
# =============================================================================


def read_csv_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            "Required representation-assessment table not found: {}".format(path)
        )
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError("Representation-assessment table is empty: {}".format(path))
    return df


def read_csv_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def load_tables(analysis_dir: Path) -> Dict[str, pd.DataFrame]:
    tables: Dict[str, pd.DataFrame] = {}
    for key, filename in REQUIRED_FILES.items():
        tables[key] = read_csv_required(analysis_dir / filename)
    for key, filename in OPTIONAL_FILES.items():
        tables[key] = read_csv_optional(analysis_dir / filename)
    return tables


def representation_order(values: Sequence[str]) -> List[str]:
    observed = [str(v) for v in values]
    ordered = [v for v in REPRESENTATION_ORDER if v in observed]
    ordered.extend(sorted(v for v in observed if v not in ordered))
    return ordered


def representation_label(value: str) -> str:
    return REPRESENTATION_LABELS.get(str(value), str(value).replace("_", " ").title())


# =============================================================================
# Plotting helpers
# =============================================================================


def draw_representation_curves(
    ax,
    curves: pd.DataFrame,
    metric: str,
    line_width: float,
    marker_size: float,
) -> None:
    colors = candidate_color_map()
    order = representation_order(curves["representation"].dropna().unique())
    for representation in order:
        d = curves[curves["representation"].astype(str) == representation].copy()
        d = d.sort_values("x_center")
        x = pd.to_numeric(d["x_center"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(d[metric], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        if not np.any(valid):
            continue
        ax.plot(
            x[valid],
            y[valid],
            marker="o",
            linewidth=line_width,
            markersize=marker_size,
            label=representation_label(representation),
            color=colors.get(representation),
        )


def prepare_depth_bars(depth: pd.DataFrame) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, List[object]]:
    d = depth.copy()
    order = representation_order(d["representation"].dropna().unique())
    d["_order"] = d["representation"].astype(str).map({v: i for i, v in enumerate(order)})
    d = d.sort_values("_order")
    values = pd.to_numeric(
        d["abs_spearman_dx_vs_depth_change"], errors="coerce"
    ).to_numpy(dtype=float)
    positions = np.arange(len(d), dtype=float)
    colors = candidate_color_map()
    bar_colors = [colors.get(str(rep)) for rep in d["representation"]]
    return d, positions, values, bar_colors


def closure_identity_table(curves: pd.DataFrame, closure: pd.DataFrame) -> pd.DataFrame:
    required_curve_cols = {"representation", "bin_id", "x_center", "mean_dx"}
    required_closure_cols = {"bin_id", "x_center", "mean_closure"}
    if not required_curve_cols.issubset(curves.columns):
        raise ValueError("Curves table lacks columns needed for closure identity panel")
    if not required_closure_cols.issubset(closure.columns):
        raise ValueError("Closure table lacks columns needed for closure identity panel")

    obs = curves.loc[
        curves["representation"].astype(str) == "observed_frequency",
        ["bin_id", "x_center", "mean_dx"],
    ].copy()
    obs = obs.rename(columns={"x_center": "x_center_obs", "mean_dx": "mean_dx_observed"})
    cnt = curves.loc[
        curves["representation"].astype(str) == "log_count",
        ["bin_id", "mean_dx"],
    ].copy()
    cnt = cnt.rename(columns={"mean_dx": "mean_dx_count"})
    clo = closure[["bin_id", "x_center", "mean_closure"]].copy()

    out = clo.merge(obs, on="bin_id", how="inner").merge(cnt, on="bin_id", how="inner")
    out["observed_minus_count"] = (
        pd.to_numeric(out["mean_dx_observed"], errors="coerce")
        - pd.to_numeric(out["mean_dx_count"], errors="coerce")
    )
    out["identity_abs_error"] = np.abs(
        pd.to_numeric(out["mean_closure"], errors="coerce")
        - pd.to_numeric(out["observed_minus_count"], errors="coerce")
    )
    return out.sort_values("x_center")


def draw_closure_identity(
    ax,
    curves: pd.DataFrame,
    closure: pd.DataFrame,
    line_width: float,
    marker_size: float,
) -> pd.DataFrame:
    d = closure_identity_table(curves, closure)
    ax.plot(
        d["x_center"],
        d["mean_closure"],
        marker="o",
        linewidth=line_width,
        markersize=marker_size,
        label="Mean closure component",
    )
    ax.plot(
        d["x_center"],
        d["observed_minus_count"],
        linewidth=max(line_width * 0.85, 0.5),
        linestyle="--",
        label="Observed − log-count mean displacement",
    )
    ax.axhline(0.0, linewidth=0.8, linestyle="--", color="0.35")
    ax.set_xlabel(r"Latent abundance, $x^\star$")
    ax.set_ylabel("Mean closure component")
    ax.legend(frameon=False, fontsize=max(float(plt.rcParams["font.size"]) - 1.0, 6.0))
    return d


def closure_depth_values(closure_depth: pd.DataFrame) -> Tuple[float, float, int]:
    if closure_depth.empty:
        raise ValueError("Closure-depth summary is empty")
    row = closure_depth.iloc[0]
    rho = float(pd.to_numeric(pd.Series([row.get("spearman_closure_vs_depth_change")]), errors="coerce").iloc[0])
    slope = float(pd.to_numeric(pd.Series([row.get("ols_slope_closure_vs_depth_change")]), errors="coerce").iloc[0])
    n_raw = pd.to_numeric(pd.Series([row.get("n_depth_valid")]), errors="coerce").iloc[0]
    n = int(n_raw) if np.isfinite(n_raw) else 0
    return rho, slope, n


def draw_closure_depth_association(ax, closure_depth: pd.DataFrame) -> Tuple[float, float, int]:
    """Plot the precomputed closure-depth association as an effect-size summary.

    The analysis table contains aggregate association statistics only, not raw or
    binned closure-versus-depth observations.  The panel therefore displays the
    Spearman effect size directly, without a stem/error-bar-like connector that
    could be mistaken for an uncertainty interval.
    """
    rho, slope, n = closure_depth_values(closure_depth)

    # Horizontal effect-size display: x is the association estimate, y is only a
    # categorical row.  The dashed vertical line denotes the null association.
    ax.plot([rho], [0.0], marker="o", markersize=6.0, linestyle="none")
    ax.axvline(0.0, linewidth=0.8, linestyle="--", color="0.35")
    # No categorical y tick is shown: with a single effect-size estimate the
    # quantity is already defined by the x-axis label, and suppressing the
    # category label prevents visual spill-over into the neighbouring panel.
    ax.set_yticks([])
    ax.set_ylim(-0.65, 0.65)
    ax.set_xlabel(r"Spearman $\rho$(closure, $\Delta\log$ depth)")
    ax.set_ylabel("")

    annotation = "OLS slope = {:.3f}\nn = {:,}".format(slope, n)
    ax.text(0.97, 0.96, annotation, transform=ax.transAxes, ha="right", va="top")
    return rho, slope, n


# =============================================================================
# Main and supplementary composites
# =============================================================================


def plot_main_composite(
    curves: pd.DataFrame,
    depth: pd.DataFrame,
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> str:
    key = "00_main_representation_benchmark"
    fig, axes = plt.subplots(1, 3, figsize=FIGURE_SIZES_INCHES[key])
    axA, axB, axC = axes

    # A: signed directional-null bias.
    d = curves.copy()
    d["ppos_bias_signed"] = pd.to_numeric(d["ppos"], errors="coerce") - 0.5
    draw_representation_curves(axA, d, "ppos_bias_signed", args.line_width, args.marker_size)
    axA.axhline(0.0, linewidth=args.reference_line_width, linestyle="--")
    axA.set_xlabel(r"Latent abundance, $x^\star$")
    axA.set_ylabel(r"P($\Delta x>0$) - 0.5")
    apply_axis_settings(axA, "MainA", fallback_xlim=args.abundance_xlim)
    style_axes(axA)
    add_panel_label(axA, "A")

    # B: finite-time fluctuation magnitude.
    draw_representation_curves(axB, curves, "msd", args.line_width, args.marker_size)
    axB.set_xlabel(r"Latent abundance, $x^\star$")
    axB.set_ylabel("Mean squared displacement")
    apply_axis_settings(axB, "MainB", fallback_xlim=args.abundance_xlim)
    style_axes(axB)
    add_panel_label(axB, "B")

    # C: residual sequencing-depth dependence.
    dd, positions, vals, bar_colors = prepare_depth_bars(depth)
    axC.bar(positions, vals, width=0.68, color=bar_colors)
    axC.set_xticks(positions)
    axC.set_xticklabels(
        [representation_label(v) for v in dd["representation"]], rotation=28, ha="right"
    )
    axC.set_ylabel(r"|Spearman($\Delta x$, $\Delta\log$ depth)|")
    apply_axis_settings(axC, "MainC")
    style_axes(axC)
    add_panel_label(axC, "C")

    handles, labels = axA.get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.50, 1.08),
            ncol=len(handles),
        )

    fig.subplots_adjust(wspace=0.48, top=0.82)
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename


def plot_supplementary_composite(
    curves: pd.DataFrame,
    closure: pd.DataFrame,
    closure_depth: pd.DataFrame,
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> Tuple[str, pd.DataFrame]:
    key = "Sx_representation_diagnostics"
    fig, axes = plt.subplots(2, 2, figsize=FIGURE_SIZES_INCHES[key])
    axA, axB, axC, axD = axes.ravel()

    # A: conditional mean displacement.
    draw_representation_curves(axA, curves, "mean_dx", args.line_width, args.marker_size)
    axA.axhline(0.0, linewidth=args.reference_line_width, linestyle="--")
    axA.set_xlabel(r"Latent abundance, $x^\star$")
    axA.set_ylabel(r"Mean $\Delta x$")
    apply_axis_settings(axA, "SxA", fallback_xlim=args.abundance_xlim)
    style_axes(axA)
    add_panel_label(axA, "A")
    axA.legend(frameon=False, title="Representation")

    # B: robust fluctuation magnitude.
    draw_representation_curves(axB, curves, "mad_dx", args.line_width, args.marker_size)
    axB.set_xlabel(r"Latent abundance, $x^\star$")
    axB.set_ylabel(r"MAD($\Delta x$)")
    apply_axis_settings(axB, "SxB", fallback_xlim=args.abundance_xlim)
    style_axes(axB)
    add_panel_label(axB, "B")
    axB.legend(frameon=False, title="Representation")

    # C: closure identity.
    identity = draw_closure_identity(
        axC, curves, closure, args.line_width, args.marker_size
    )
    apply_axis_settings(axC, "SxC", fallback_xlim=args.abundance_xlim)
    style_axes(axC)
    add_panel_label(axC, "C")

    # D: precomputed closure-depth association effect-size summary.
    draw_closure_depth_association(axD, closure_depth)
    apply_axis_settings(axD, "SxD")
    style_axes(axD)
    add_panel_label(axD, "D")

    fig.subplots_adjust(wspace=0.40, hspace=0.42)
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename, identity


# =============================================================================
# Standalone figures
# =============================================================================


def plot_primary_curve(
    curves: pd.DataFrame,
    metric: str,
    ylabel: str,
    baseline: Optional[float],
    key: str,
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> str:
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[key])
    draw_representation_curves(
        ax,
        curves,
        metric=metric,
        line_width=args.line_width,
        marker_size=args.marker_size,
    )
    if baseline is not None:
        ax.axhline(baseline, linewidth=args.reference_line_width, linestyle="--")
    ax.set_xlabel(r"Latent abundance, $x^\star$")
    ax.set_ylabel(ylabel)
    apply_axis_settings(ax, key, fallback_xlim=args.abundance_xlim)
    style_axes(ax)
    axis_legend(ax, "Representation")
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename


def plot_ppos_bias_curve(
    curves: pd.DataFrame,
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> str:
    key = "03_primary_ppos_bias"
    d = curves.copy()
    d["ppos_bias_signed"] = pd.to_numeric(d["ppos"], errors="coerce") - 0.5
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[key])
    draw_representation_curves(
        ax,
        d,
        metric="ppos_bias_signed",
        line_width=args.line_width,
        marker_size=args.marker_size,
    )
    ax.axhline(0.0, linewidth=args.reference_line_width, linestyle="--")
    ax.set_xlabel(r"Latent abundance, $x^\star$")
    ax.set_ylabel(r"P($\Delta x>0$) - 0.5")
    apply_axis_settings(ax, key, fallback_xlim=args.abundance_xlim)
    style_axes(ax)
    axis_legend(ax, "Representation")
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename


def plot_depth_dependence(
    depth: pd.DataFrame,
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> str:
    key = "04_depth_change_dependence"
    d, positions, values, bar_colors = prepare_depth_bars(depth)
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[key])
    ax.bar(positions, values, width=0.68, color=bar_colors)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [representation_label(v) for v in d["representation"]], rotation=25, ha="right"
    )
    ax.set_ylabel(r"|Spearman($\Delta x$, $\Delta\log$ depth)|")
    apply_axis_settings(ax, key)
    style_axes(ax)
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename


def plot_closure_curve(
    closure: pd.DataFrame,
    metric: str,
    ylabel: str,
    key: str,
    baseline: Optional[float],
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> Optional[str]:
    if closure.empty or metric not in closure.columns:
        return None
    d = closure.sort_values("x_center")
    x = pd.to_numeric(d["x_center"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(d[metric], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if not np.any(valid):
        return None
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[key])
    ax.plot(
        x[valid],
        y[valid],
        marker="o",
        linewidth=args.line_width,
        markersize=args.marker_size,
    )
    if baseline is not None:
        ax.axhline(baseline, linewidth=args.reference_line_width, linestyle="--")
    ax.set_xlabel(r"Latent abundance, $x^\star$")
    ax.set_ylabel(ylabel)
    apply_axis_settings(ax, key, fallback_xlim=args.abundance_xlim)
    style_axes(ax)
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename


def plot_closure_mean_identity(
    curves: pd.DataFrame,
    closure: pd.DataFrame,
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> Tuple[str, pd.DataFrame]:
    key = "13_closure_mean_identity"
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[key])
    identity = draw_closure_identity(ax, curves, closure, args.line_width, args.marker_size)
    apply_axis_settings(ax, key, fallback_xlim=args.abundance_xlim)
    style_axes(ax)
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename, identity


def plot_closure_depth_association(
    closure_depth: pd.DataFrame,
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> str:
    key = "14_closure_depth_association"
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[key])
    draw_closure_depth_association(ax, closure_depth)
    apply_axis_settings(ax, key)
    style_axes(ax)
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename


def ordered_scenarios(d: pd.DataFrame) -> List[str]:
    observed = set(d["scenario"].dropna().astype(str))
    scenarios = [s for s in SCENARIO_ORDER if s in observed]
    scenarios.extend(sorted(s for s in observed if s not in scenarios))
    return scenarios


def plot_sensitivity_metric(
    sensitivity: pd.DataFrame,
    metric: str,
    ylabel: str,
    key: str,
    figdir: Path,
    fmt: str,
    args: argparse.Namespace,
) -> Optional[str]:
    if sensitivity.empty or metric not in sensitivity.columns:
        return None
    d = sensitivity.copy()
    scenarios = ordered_scenarios(d)
    reps = representation_order(d["representation"].dropna().unique())
    if not scenarios or not reps:
        return None

    colors = candidate_color_map()
    x = np.arange(len(scenarios), dtype=float)
    width = 0.8 / max(len(reps), 1)
    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[key])

    for index, rep in enumerate(reps):
        values = []
        for scenario in scenarios:
            row = d[
                (d["scenario"].astype(str) == scenario)
                & (d["representation"].astype(str) == rep)
            ]
            if row.empty:
                values.append(np.nan)
            else:
                values.append(float(pd.to_numeric(row[metric], errors="coerce").iloc[0]))
        offsets = x - 0.4 + width / 2.0 + index * width
        ax.bar(
            offsets,
            values,
            width=width,
            label=representation_label(rep),
            color=colors.get(rep),
        )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [SCENARIO_LABELS.get(s, s.replace("_", " ").title()) for s in scenarios]
    )
    ax.set_ylabel(ylabel)
    apply_axis_settings(ax, key)
    style_axes(ax)
    axis_legend(ax, "Representation")
    filename = "{}.{}".format(key, fmt)
    save_figure(fig, figdir / filename, args.dpi)
    return filename


# =============================================================================
# CLI and main
# =============================================================================


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Plot ClonoDynamics Step-8 pseudo-longitudinal "
            "representation-assessment outputs."
        ),
    )
    parser.add_argument(
        "--analysis-dir",
        required=True,
        type=Path,
        help=(
            "Output directory produced by "
            "8-pseudo_longitudinal_representation_assessment.py."
        ),
    )
    parser.add_argument(
        "--figures-dir",
        "--figdir",
        dest="figures_dir",
        type=Path,
        default=None,
        help=(
            "Figure output directory. Default: "
            "<analysis-dir>/figures_pseudo_longitudinal_representation_assessment"
        ),
    )
    parser.add_argument(
        "--formats",
        default="pdf,png",
        help="Comma-separated formats from: pdf,png,svg.",
    )

    # Legacy aliases retained only for local compatibility.
    parser.add_argument("--fmt", default=None, choices=["pdf", "png", "svg"], help=argparse.SUPPRESS)
    parser.add_argument("--also-png", action="store_true", help=argparse.SUPPRESS)

    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-size", type=float, default=8.0)
    parser.add_argument("--line-width", type=float, default=1.0)
    parser.add_argument("--reference-line-width", type=float, default=0.8)
    parser.add_argument("--marker-size", type=float, default=3.0)
    parser.add_argument(
        "--abundance-xlim",
        default="",
        help=(
            "Optional display-only fallback range lo,hi. "
            "AXIS_CONFIG takes precedence."
        ),
    )
    parser.add_argument(
        "--figure-set",
        default="all",
        choices=["all", "main", "supplementary", "standalone"],
        help="Choose manuscript composites, standalones, or all figures.",
    )
    parser.add_argument(
        "--only",
        default="",
        help="Optional comma-separated figure keys; overrides --figure-set.",
    )
    return parser


def default_selected_keys(figure_set: str) -> List[str]:
    if figure_set == "main":
        return ["00_main_representation_benchmark"]
    if figure_set == "supplementary":
        return ["Sx_representation_diagnostics"]
    if figure_set == "standalone":
        return [
            "01_primary_mean_dx",
            "02_primary_median_dx",
            "03_primary_ppos_bias",
            "04_depth_change_dependence",
            "05_primary_msd",
            "06_primary_mad_dx",
            "07_closure_median",
            "08_closure_mad",
            "09_sensitivity_depth_dependence",
            "10_sensitivity_directional_mean_bias",
            "11_sensitivity_directional_median_bias",
            "12_sensitivity_directional_ppos_bias",
            "13_closure_mean_identity",
            "14_closure_depth_association",
        ]
    return []


def should_plot(key: str, selected: Sequence[str], figure_set: str) -> bool:
    if selected:
        return key in set(selected)
    return figure_set == "all" or key in set(default_selected_keys(figure_set))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    analysis_dir = args.analysis_dir.expanduser().resolve()
    if not analysis_dir.is_dir():
        raise FileNotFoundError(
            "Step-8 analysis directory not found: {}".format(analysis_dir)
        )

    figures_dir = (
        args.figures_dir.expanduser().resolve()
        if args.figures_dir is not None
        else analysis_dir / "figures_pseudo_longitudinal_representation_assessment"
    )
    figures_dir.mkdir(parents=True, exist_ok=True)

    formats = parse_formats(args.formats, args.fmt, bool(args.also_png))
    setup_style(args.font_size)

    tables = load_tables(analysis_dir)
    curves = tables["curves"]
    depth = tables["depth"]
    closure = tables["closure"]
    closure_depth = tables["closure_depth"]
    sensitivity = tables["sensitivity"]
    selected = [v.strip() for v in args.only.split(",") if v.strip()]

    if sensitivity.empty:
        print(
            "[NOTE] optional 12_sensitivity_representation_assessment.csv "
            "not available; sensitivity panels will be skipped."
        )

    records: List[Dict[str, object]] = []

    for fmt in formats:
        if should_plot("00_main_representation_benchmark", selected, args.figure_set):
            filename = plot_main_composite(curves, depth, figures_dir, fmt, args)
            records.append({
                "figure": filename,
                "family": "main",
                "format": fmt,
                "source_table": "03_pooled_representation_curves.csv;08_depth_dependence_summary.csv",
            })

        if should_plot("Sx_representation_diagnostics", selected, args.figure_set):
            filename, _ = plot_supplementary_composite(
                curves, closure, closure_depth, figures_dir, fmt, args
            )
            records.append({
                "figure": filename,
                "family": "supplementary",
                "format": fmt,
                "source_table": (
                    "03_pooled_representation_curves.csv;"
                    "09_closure_diagnostic_curves.csv;"
                    "10_closure_depth_summary.csv"
                ),
            })

        for key, metric, ylabel, baseline in [
            ("01_primary_mean_dx", "mean_dx", r"Mean $\Delta x$", 0.0),
            ("02_primary_median_dx", "median_dx", r"Median $\Delta x$", 0.0),
        ]:
            if should_plot(key, selected, args.figure_set):
                filename = plot_primary_curve(
                    curves, metric, ylabel, baseline, key, figures_dir, fmt, args
                )
                records.append({
                    "figure": filename, "family": key, "format": fmt,
                    "source_table": "03_pooled_representation_curves.csv",
                })

        if should_plot("03_primary_ppos_bias", selected, args.figure_set):
            filename = plot_ppos_bias_curve(curves, figures_dir, fmt, args)
            records.append({
                "figure": filename, "family": "03_primary_ppos_bias", "format": fmt,
                "source_table": "03_pooled_representation_curves.csv",
            })

        if should_plot("04_depth_change_dependence", selected, args.figure_set):
            filename = plot_depth_dependence(depth, figures_dir, fmt, args)
            records.append({
                "figure": filename, "family": "04_depth_change_dependence", "format": fmt,
                "source_table": "08_depth_dependence_summary.csv",
            })

        for key, metric, ylabel, baseline in [
            ("05_primary_msd", "msd", "Mean squared displacement", None),
            ("06_primary_mad_dx", "mad_dx", r"MAD($\Delta x$)", None),
        ]:
            if should_plot(key, selected, args.figure_set):
                filename = plot_primary_curve(
                    curves, metric, ylabel, baseline, key, figures_dir, fmt, args
                )
                records.append({
                    "figure": filename, "family": key, "format": fmt,
                    "source_table": "03_pooled_representation_curves.csv",
                })

        if should_plot("07_closure_median", selected, args.figure_set):
            filename = plot_closure_curve(
                closure, "median_closure", "Median closure component",
                "07_closure_median", 0.0, figures_dir, fmt, args,
            )
            if filename is not None:
                records.append({
                    "figure": filename, "family": "07_closure_median", "format": fmt,
                    "source_table": "09_closure_diagnostic_curves.csv",
                })

        if should_plot("08_closure_mad", selected, args.figure_set):
            filename = plot_closure_curve(
                closure, "mad_closure", "MAD(closure component)",
                "08_closure_mad", None, figures_dir, fmt, args,
            )
            if filename is not None:
                records.append({
                    "figure": filename, "family": "08_closure_mad", "format": fmt,
                    "source_table": "09_closure_diagnostic_curves.csv",
                })

        if should_plot("13_closure_mean_identity", selected, args.figure_set):
            filename, _ = plot_closure_mean_identity(
                curves, closure, figures_dir, fmt, args
            )
            records.append({
                "figure": filename, "family": "13_closure_mean_identity", "format": fmt,
                "source_table": "03_pooled_representation_curves.csv;09_closure_diagnostic_curves.csv",
            })

        if should_plot("14_closure_depth_association", selected, args.figure_set):
            filename = plot_closure_depth_association(
                closure_depth, figures_dir, fmt, args
            )
            records.append({
                "figure": filename, "family": "14_closure_depth_association", "format": fmt,
                "source_table": "10_closure_depth_summary.csv",
            })

        for key, metric, ylabel in [
            (
                "09_sensitivity_depth_dependence",
                "abs_spearman_dx_vs_depth_change",
                r"|Spearman($\Delta x$, $\Delta\log$ depth)|",
            ),
            (
                "10_sensitivity_directional_mean_bias",
                "mean_abs_mean_dx",
                r"Mean |mean $\Delta x$|",
            ),
            (
                "11_sensitivity_directional_median_bias",
                "mean_abs_median_dx",
                r"Mean |median $\Delta x$|",
            ),
            (
                "12_sensitivity_directional_ppos_bias",
                "mean_abs_ppos_bias",
                r"Mean |P($\Delta x>0$)-0.5|",
            ),
        ]:
            if should_plot(key, selected, args.figure_set):
                filename = plot_sensitivity_metric(
                    sensitivity, metric, ylabel, key, figures_dir, fmt, args
                )
                if filename is not None:
                    records.append({
                        "figure": filename, "family": key, "format": fmt,
                        "source_table": "12_sensitivity_representation_assessment.csv",
                    })

    pd.DataFrame(records).to_csv(
        figures_dir / "figure_manifest.csv",
        index=False,
    )

    plot_config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "analysis_source": "8-pseudo_longitudinal_representation_assessment.py",
        "analysis_dir": str(analysis_dir),
        "figures_dir": str(figures_dir),
        "formats": formats,
        "font_size": args.font_size,
        "line_width": args.line_width,
        "reference_line_width": args.reference_line_width,
        "marker_size": args.marker_size,
        "dpi": args.dpi,
        "figure_set": args.figure_set,
        "only": selected,
        "figure_sizes_inches": FIGURE_SIZES_INCHES,
        "axis_config": AXIS_CONFIG,
        "conditioning_coordinate": "xstar_latent",
        "automatic_representation_winner": False,
        "dispersion_is_selection_criterion": False,
        "closure_depth_statistics_are_precomputed": True,
        "closure_identity_is_display_only": True,
    }
    (figures_dir / "plot_config.json").write_text(
        json.dumps(plot_config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("[INFO] {}".format(SCRIPT_VERSION))
    print("[INFO] Step-8 representation-assessment plotting complete")
    print("[INFO] Analysis directory: {}".format(analysis_dir))
    print("[INFO] Figures written to: {}".format(figures_dir))
    print("[INFO] No representation ranking was recomputed by the plotter")
    return 0


if __name__ == "__main__":
    sys.exit(main())
