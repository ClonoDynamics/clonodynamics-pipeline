#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
13-plot_temporal_fluctuation_scaling.py
=======================================

Plot-only companion to ClonoDynamics Step 13:

    13-temporal_fluctuation_scaling.py

Scientific role
---------------
The plotter visualizes finite-lag model-comparison outputs already computed by
Step 13 for the genuine longitudinal Step-11 fluctuation data.

Primary architecture:

    representation = latent
    conditioning   = xstar
    mode           = TT
    primary metric = Var(dx | xstar, dt)
    complementary  = MSD(dx | xstar, dt)

Primary estimator:

    transition_weighted

Sensitivity estimator:

    equal_subject_weighted

The plotter does not refit temporal models, rerun the joint subject bootstrap,
redefine the temporal abundance core, recompute AIC/AICc or model-preference
fractions, or subtract the Step-12 pseudo technical baseline.

Required Step-13 inputs
-----------------------
    03_observed_six_model_fits_long.csv
    05_model_comparison_by_bin.csv
    06_model_comparison_summary.csv

Optional provenance inputs
--------------------------
    02_tested_abundance_bins.csv
    07_interpretation_aid.csv
    08_subject_availability_and_requirements_by_dt.csv

Figures
-------
Main:

    00_main_temporal_fluctuation_scaling_ABC
        A  observed conditional-variance profiles across dt and xstar
        B  constant-model bootstrap preference fraction for Var
        C  signed-linear Var slope with joint-bootstrap interval

Supplementary:

    Sx_temporal_fluctuation_scaling_diagnostics_ABCD_ABCD
        A  observed MSD profiles
        B  constant-model bootstrap preference fraction for MSD
        C  signed-linear MSD slope with joint-bootstrap interval
        D  transition-weighted vs equal-subject-weighted Var slope

Standalone:

    01_variance_temporal_heatmap
    02_variance_constant_preference
    03_variance_signed_slope
    04_msd_temporal_heatmap
    05_msd_constant_preference
    06_msd_signed_slope
    07_variance_slope_weighting_sensitivity
    08_variance_temporal_profiles
    09_msd_temporal_profiles

Non-figure outputs
------------------
    figure_manifest.csv
    plot_config.json

Rendering boundary
------------------
The temporal-profile matrix is reconstructed only by parsing the `dt_values`
and `metric_values` strings already stored in the observed Step-13 fit table.
The constant-model row is used only to avoid duplicate copies of the same
upstream observed profile. No smoothing or fitting is performed.

Figure customization
--------------------
Use:

    FIGURE_SIZES_INCHES
    AXIS_CONFIG

with per-panel controls:

    xlim
    ylim
    xscale
    yscale
    xtick_step
    ytick_step

Example
-------
python3 ./code/plotting/13-plot_temporal_fluctuation_scaling.py \
    --analysis-dir ./results/13-temporal-fluctuation-scaling \
    --figures-dir ./figures/13-temporal-fluctuation-scaling \
    --formats pdf,png
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


SCRIPT_VERSION = "v2-step13-public-plotter-2026-08-25"
PRIMARY_ESTIMATOR = "transition_weighted"
SENSITIVITY_ESTIMATOR = "equal_subject_weighted"
PRIMARY_METRIC = "var_dx"
COMPLEMENTARY_METRIC = "msd"


# =============================================================================
# Figure configuration
# =============================================================================

FIGURE_SIZES_INCHES = {
    "00_main_temporal_fluctuation_scaling_ABC": (9, 2.4),
    "Sx_temporal_fluctuation_scaling_diagnostics_ABCD": (6.8, 5.4),

    "01_variance_temporal_heatmap": (3.2, 2.4),
    "02_variance_constant_preference": (3.2, 2.4),
    "03_variance_signed_slope": (3.2, 2.4),
    "04_msd_temporal_heatmap": (3.2, 2.4),
    "05_msd_constant_preference": (3.2, 2.4),
    "06_msd_signed_slope": (3.2, 2.4),
    "07_variance_slope_weighting_sensitivity": (3.2, 2.4),
    "08_variance_temporal_profiles": (3.2, 2.4),
    "09_msd_temporal_profiles": (3.2, 2.4),
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


AXIS_CONFIG = {
    "MainA": _axis_cfg(),
    "MainB": _axis_cfg(ylim=(0.0, 1.02)),
    "MainC": _axis_cfg(),

    "SxA": _axis_cfg(),
    "SxB": _axis_cfg(ylim=(0.0, 1.02)),
    "SxC": _axis_cfg(),
    "SxD": _axis_cfg(),

    "01_variance_temporal_heatmap": _axis_cfg(),
    "02_variance_constant_preference": _axis_cfg(ylim=(0.0, 1.02)),
    "03_variance_signed_slope": _axis_cfg(),
    "04_msd_temporal_heatmap": _axis_cfg(),
    "05_msd_constant_preference": _axis_cfg(ylim=(0.0, 1.02)),
    "06_msd_signed_slope": _axis_cfg(),
    "07_variance_slope_weighting_sensitivity": _axis_cfg(),
    "08_variance_temporal_profiles": _axis_cfg(xtick_step=1.0),
    "09_msd_temporal_profiles": _axis_cfg(xtick_step=1.0),
}


REQUIRED_FILE_PREFIXES = {
    "fits": "03_observed_six_model_fits_long",
    "comparison": "05_model_comparison_by_bin",
    "summary": "06_model_comparison_summary",
}

OPTIONAL_FILE_PREFIXES = {
    "tested_bins": "02_tested_abundance_bins",
    "interpretation": "07_interpretation_aid",
    "availability": "08_subject_availability_and_requirements_by_dt",
}


# =============================================================================
# Style / I/O helpers
# =============================================================================


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def setup_style(font_size: float) -> str:
    try:
        from matplotlib import font_manager

        font_manager.findfont("Arial", fallback_to_default=False)
        family = "Arial"
    except Exception:
        family = "DejaVu Sans"
        eprint("[WARN] Arial not found; using DejaVu Sans fallback.")

    tick_font = max(float(font_size) - 1.0, 5.0)
    plt.rcParams.update(
        {
            "font.family": family,
            "font.size": float(font_size),
            "axes.labelsize": float(font_size),
            "axes.titlesize": float(font_size),
            "xtick.labelsize": tick_font,
            "ytick.labelsize": tick_font,
            "legend.fontsize": tick_font,
            "legend.title_fontsize": tick_font,
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
    return family


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")


def add_panel_label(ax, label: str, font_size: float) -> None:
    ax.text(
        -0.16,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=max(float(font_size) + 1.0, 8.0),
        fontweight="bold",
        va="top",
        ha="left",
    )


def apply_axis_config(ax, key: str) -> None:
    cfg = AXIS_CONFIG.get(key, {})

    if cfg.get("xscale"):
        ax.set_xscale(str(cfg["xscale"]))
    if cfg.get("yscale"):
        ax.set_yscale(str(cfg["yscale"]))

    if cfg.get("xlim") is not None:
        lo, hi = cfg["xlim"]
        old_lo, old_hi = ax.get_xlim()
        ax.set_xlim(old_lo if lo is None else lo, old_hi if hi is None else hi)

    if cfg.get("ylim") is not None:
        lo, hi = cfg["ylim"]
        old_lo, old_hi = ax.get_ylim()
        ax.set_ylim(old_lo if lo is None else lo, old_hi if hi is None else hi)

    if cfg.get("xtick_step") is not None:
        ax.xaxis.set_major_locator(MultipleLocator(float(cfg["xtick_step"])))
    if cfg.get("ytick_step") is not None:
        ax.yaxis.set_major_locator(MultipleLocator(float(cfg["ytick_step"])))


def resolve_csv(analysis_dir: Path, prefix: str) -> Path:
    exact = analysis_dir / (prefix + ".csv")
    if exact.exists():
        return exact

    matches = sorted(analysis_dir.glob(prefix + "*.csv"))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(
            "Required Step-13 output not found for prefix {!r} in {}".format(
                prefix, analysis_dir
            )
        )
    raise RuntimeError(
        "Multiple Step-13 files match prefix {!r}: {}. "
        "Use a directory containing canonical outputs.".format(
            prefix, [str(p) for p in matches]
        )
    )


def read_csv_required(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError("Step-13 output is empty: {}".format(path))
    return df


def load_tables(
    analysis_dir: Path,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, str]]:
    tables: Dict[str, pd.DataFrame] = {}
    paths: Dict[str, str] = {}

    for key, prefix in REQUIRED_FILE_PREFIXES.items():
        path = resolve_csv(analysis_dir, prefix)
        tables[key] = read_csv_required(path)
        paths[key] = str(path)

    for key, prefix in OPTIONAL_FILE_PREFIXES.items():
        path = analysis_dir / (prefix + ".csv")
        if path.exists():
            tables[key] = pd.read_csv(path)
            paths[key] = str(path)
        else:
            tables[key] = pd.DataFrame()

    return tables, paths


def require_columns(df: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            "{} is missing columns {}. Available: {}".format(
                label, missing, list(df.columns)
            )
        )


def parse_formats(text: str) -> List[str]:
    out: List[str] = []
    for value in str(text).split(","):
        value = value.strip().lower().lstrip(".")
        if not value:
            continue
        if value not in {"pdf", "png", "svg"}:
            raise ValueError("Unsupported format: {}".format(value))
        if value not in out:
            out.append(value)
    return out or ["pdf", "png"]


def save_figure(
    fig,
    outdir: Path,
    basename: str,
    formats: Sequence[str],
    dpi: int,
) -> List[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    for fmt in formats:
        path = outdir / "{}.{}".format(basename, fmt)
        kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            kwargs["dpi"] = int(dpi)
        fig.savefig(path, **kwargs)
        paths.append(path)
        print("[WRITE] {}".format(path))
    plt.close(fig)
    return paths


def parse_float_list(value: object) -> np.ndarray:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.array([], dtype=float)
    text = str(value).strip()
    if not text:
        return np.array([], dtype=float)
    return np.asarray([float(x) for x in text.split(",") if str(x).strip()], dtype=float)


def subset_table(
    df: pd.DataFrame,
    estimator: str,
    metric: str,
) -> pd.DataFrame:
    require_columns(df, ["estimator", "metric"], "Step-13 table")
    out = df.loc[
        (df["estimator"].astype(str) == str(estimator))
        & (df["metric"].astype(str) == str(metric))
    ].copy()
    if out.empty:
        raise ValueError(
            "No rows for estimator={!r}, metric={!r}".format(estimator, metric)
        )
    return out


def summary_row(summary: pd.DataFrame, estimator: str, metric: str) -> pd.Series:
    d = subset_table(summary, estimator, metric)
    if len(d) != 1:
        raise ValueError(
            "Expected one summary row for estimator={}, metric={}; found {}".format(
                estimator, metric, len(d)
            )
        )
    return d.iloc[0]


# =============================================================================
# Data preparation for plotting only
# =============================================================================


def observed_profile_matrix(
    fits: pd.DataFrame,
    estimator: str,
    metric: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return dt values, sorted xstar bin centers, and observed metric matrix.

    metric_values and dt_values are written upstream and are identical across
    model rows for a given estimator x metric x bin. We use the constant-model
    row only to avoid duplicates; no model fit is recomputed here.
    """
    d = subset_table(fits, estimator, metric)
    require_columns(
        d,
        ["model", "bin_id", "bin_mid", "dt_values", "metric_values"],
        "03_observed_six_model_fits_long.csv",
    )
    d = d.loc[d["model"].astype(str) == "constant"].copy()
    if d.empty:
        raise ValueError("Constant-model rows absent from observed fit table")
    d = d.sort_values("bin_mid")

    dt_ref: Optional[np.ndarray] = None
    rows: List[np.ndarray] = []
    xstar: List[float] = []

    for _, row in d.iterrows():
        dt = parse_float_list(row["dt_values"])
        values = parse_float_list(row["metric_values"])
        if dt.size != values.size or dt.size == 0:
            raise ValueError(
                "Malformed dt_values/metric_values for bin {}".format(row["bin_id"])
            )
        if dt_ref is None:
            dt_ref = dt
        elif dt.size != dt_ref.size or not np.allclose(dt, dt_ref):
            raise ValueError("Inconsistent dt grid across abundance bins")
        rows.append(values)
        xstar.append(float(row["bin_mid"]))

    assert dt_ref is not None
    return dt_ref, np.asarray(xstar, dtype=float), np.vstack(rows)


def comparison_subset(
    comparison: pd.DataFrame,
    estimator: str,
    metric: str,
) -> pd.DataFrame:
    d = subset_table(comparison, estimator, metric)
    require_columns(
        d,
        [
            "bin_id",
            "bin_mid",
            "observed_preferred_model",
            "constant_preference_fraction",
            "signed_linear_preference_fraction",
            "signed_linear_slope_median",
            "signed_linear_slope_q_low",
            "signed_linear_slope_q_high",
            "signed_linear_slope_ci_below_zero",
            "signed_linear_slope_ci_above_zero",
            "signed_trend_classification",
        ],
        "05_model_comparison_by_bin.csv",
    )
    return d.sort_values("bin_mid").reset_index(drop=True)


# =============================================================================
# Plotting functions
# =============================================================================


def draw_temporal_heatmap(
    ax,
    fits: pd.DataFrame,
    estimator: str,
    metric: str,
    axis_key: str,
    cbar_title: str,
) -> None:
    dt, xstar, matrix = observed_profile_matrix(fits, estimator, metric)

    # imshow uses the exact upstream matrix; no smoothing/interpolation.
    y_half = (
        0.5 * float(np.median(np.diff(xstar)))
        if xstar.size > 1
        else 0.5
    )
    image = ax.imshow(
        matrix,
        aspect="auto",
        origin="lower",
        extent=(
            float(dt.min()) - 0.5,
            float(dt.max()) + 0.5,
            float(xstar.min()) - y_half,
            float(xstar.max()) + y_half,
        ),
    )
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.ax.set_title(cbar_title, fontsize=plt.rcParams["font.size"], pad=3)
    cbar.ax.tick_params(width=0.7, length=2.5)

    ax.set_xlabel(r"Temporal lag, $\Delta t$ (weeks)")
    ax.set_ylabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_xticks(dt)
    style_axes(ax)
    apply_axis_config(ax, axis_key)


def draw_constant_preference(
    ax,
    comparison: pd.DataFrame,
    summary: pd.DataFrame,
    estimator: str,
    metric: str,
    axis_key: str,
) -> None:
    d = comparison_subset(comparison, estimator, metric)
    s = summary_row(summary, estimator, metric)

    x = pd.to_numeric(d["bin_mid"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(
        d["constant_preference_fraction"], errors="coerce"
    ).to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)

    line, = ax.plot(x[valid], y[valid], marker="o", linewidth=1.15, markersize=3.0)

    threshold = float(s.get("model_preference_threshold", np.nan))
    if np.isfinite(threshold):
        ax.axhline(threshold, linewidth=0.7, linestyle=":", label="Preference threshold")

    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel("Constant-model bootstrap\npreference fraction")
    style_axes(ax)
    apply_axis_config(ax, axis_key)

    n_const = int(s["n_observed_constant_best"])
    n_bins = int(s["n_bins_tested"])
    median_pref = float(s["median_constant_preference_fraction"])
    ax.text(
        0.03,
        0.05,
        "Observed constant best: {}/{}\nMedian bootstrap preference = {:.3f}".format(
            n_const, n_bins, median_pref
        ),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=max(float(plt.rcParams["font.size"]) - 0.5, 5.0),
    )


def draw_signed_slope(
    ax,
    comparison: pd.DataFrame,
    summary: pd.DataFrame,
    estimator: str,
    metric: str,
    axis_key: str,
    show_annotation: bool = True,
) -> None:
    d = comparison_subset(comparison, estimator, metric)
    s = summary_row(summary, estimator, metric)

    x = pd.to_numeric(d["bin_mid"], errors="coerce").to_numpy(dtype=float)
    median = pd.to_numeric(d["signed_linear_slope_median"], errors="coerce").to_numpy(dtype=float)
    low = pd.to_numeric(d["signed_linear_slope_q_low"], errors="coerce").to_numpy(dtype=float)
    high = pd.to_numeric(d["signed_linear_slope_q_high"], errors="coerce").to_numpy(dtype=float)
    supported_neg = d["signed_linear_slope_ci_below_zero"].astype(bool).to_numpy()
    supported_pos = d["signed_linear_slope_ci_above_zero"].astype(bool).to_numpy()

    valid = np.isfinite(x) & np.isfinite(median) & np.isfinite(low) & np.isfinite(high)
    yerr = np.vstack([median - low, high - median])

    ax.errorbar(
        x[valid],
        median[valid],
        yerr=yerr[:, valid],
        fmt="o",
        markersize=3.0,
        linewidth=0.9,
        capsize=2.0,
        label="Bootstrap median and 95% interval",
    )

    neg = valid & supported_neg
    pos = valid & supported_pos
    if np.any(neg):
        ax.plot(
            x[neg],
            median[neg],
            linestyle="none",
            marker="D",
            markersize=4.2,
            label="95% interval below 0",
        )
    if np.any(pos):
        ax.plot(
            x[pos],
            median[pos],
            linestyle="none",
            marker="s",
            markersize=4.2,
            label="95% interval above 0",
        )

    ax.axhline(0.0, linewidth=0.7, linestyle=":")
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel(r"Signed-linear slope, $D$")
    style_axes(ax)
    apply_axis_config(ax, axis_key)

    if show_annotation:
        n_below = int(s["n_bins_signed_linear_slope_ci_below_zero"])
        n_above = int(s["n_bins_signed_linear_slope_ci_above_zero"])
        n_bins = int(s["n_bins_tested"])
        ax.text(
            0.03,
            0.05,
            "95% CI < 0: {}/{}\n95% CI > 0: {}/{}".format(
                n_below, n_bins, n_above, n_bins
            ),
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=max(float(plt.rcParams["font.size"]) - 0.5, 5.0),
        )


def draw_weighting_sensitivity(
    ax,
    comparison: pd.DataFrame,
    metric: str,
    axis_key: str,
) -> None:
    d_primary = comparison_subset(comparison, PRIMARY_ESTIMATOR, metric)
    d_sens = comparison_subset(comparison, SENSITIVITY_ESTIMATOR, metric)

    # Inner merge only aligns already-defined abundance bins; no statistic is recomputed.
    cols = [
        "bin_id",
        "bin_mid",
        "signed_linear_slope_median",
        "signed_linear_slope_q_low",
        "signed_linear_slope_q_high",
    ]
    a = d_primary[cols].rename(
        columns={
            "bin_mid": "bin_mid_primary",
            "signed_linear_slope_median": "median_primary",
            "signed_linear_slope_q_low": "low_primary",
            "signed_linear_slope_q_high": "high_primary",
        }
    )
    b = d_sens[cols].rename(
        columns={
            "bin_mid": "bin_mid_sensitivity",
            "signed_linear_slope_median": "median_sensitivity",
            "signed_linear_slope_q_low": "low_sensitivity",
            "signed_linear_slope_q_high": "high_sensitivity",
        }
    )
    m = a.merge(b, on="bin_id", how="inner").sort_values("bin_mid_primary")

    x = pd.to_numeric(m["bin_mid_primary"], errors="coerce").to_numpy(dtype=float)
    offset = 0.025

    for label, median_col, low_col, high_col, shift in [
        (
            "Transition-weighted",
            "median_primary",
            "low_primary",
            "high_primary",
            -offset,
        ),
        (
            "Equal-subject weighted",
            "median_sensitivity",
            "low_sensitivity",
            "high_sensitivity",
            offset,
        ),
    ]:
        med = pd.to_numeric(m[median_col], errors="coerce").to_numpy(dtype=float)
        lo = pd.to_numeric(m[low_col], errors="coerce").to_numpy(dtype=float)
        hi = pd.to_numeric(m[high_col], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(med) & np.isfinite(lo) & np.isfinite(hi)
        ax.errorbar(
            x[valid] + shift,
            med[valid],
            yerr=np.vstack([med - lo, hi - med])[:, valid],
            fmt="o",
            markersize=2.8,
            linewidth=0.8,
            capsize=1.8,
            label=label,
        )

    ax.axhline(0.0, linewidth=0.7, linestyle=":")
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel(r"Signed-linear slope, $D$")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False, loc="lower left")


def draw_temporal_profiles(
    ax,
    fits: pd.DataFrame,
    estimator: str,
    metric: str,
    axis_key: str,
) -> None:
    dt, xstar, matrix = observed_profile_matrix(fits, estimator, metric)

    for i in range(matrix.shape[0]):
        ax.plot(
            dt,
            matrix[i, :],
            marker="o",
            markersize=2.4,
            linewidth=0.9,
            label=r"$x^\star$={:.2f}".format(xstar[i]),
        )

    ax.set_xlabel(r"Temporal lag, $\Delta t$ (weeks)")
    ylabel = r"Var$(\Delta x\mid x^\star,\Delta t)$" if metric == "var_dx" else r"MSD$(\Delta x\mid x^\star,\Delta t)$"
    ax.set_ylabel(ylabel)
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        title=r"Abundance bin",
    )


# =============================================================================
# Composite figures
# =============================================================================


def plot_main(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    dpi: int,
    font_size: float,
) -> List[Path]:
    fig, axes = plt.subplots(
        1,
        3,
        figsize=FIGURE_SIZES_INCHES["00_main_temporal_fluctuation_scaling_ABC"],
    )

    draw_temporal_heatmap(
        axes[0],
        tables["fits"],
        PRIMARY_ESTIMATOR,
        PRIMARY_METRIC,
        "MainA",
        "Var",
    )
    draw_constant_preference(
        axes[1],
        tables["comparison"],
        tables["summary"],
        PRIMARY_ESTIMATOR,
        PRIMARY_METRIC,
        "MainB",
    )
    draw_signed_slope(
        axes[2],
        tables["comparison"],
        tables["summary"],
        PRIMARY_ESTIMATOR,
        PRIMARY_METRIC,
        "MainC",
    )

    for ax, label in zip(axes, ["A", "B", "C"]):
        add_panel_label(ax, label, font_size)

    fig.subplots_adjust(wspace=0.48)
    return save_figure(
        fig,
        outdir,
        "00_main_temporal_fluctuation_scaling_ABC",
        formats,
        dpi,
    )


def plot_supplementary(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    dpi: int,
    font_size: float,
) -> List[Path]:
    fig, axes = plt.subplots(
        2,
        2,
        figsize=FIGURE_SIZES_INCHES["Sx_temporal_fluctuation_scaling_diagnostics_ABCD"],
    )
    axA, axB, axC, axD = axes.ravel()

    draw_temporal_heatmap(
        axA,
        tables["fits"],
        PRIMARY_ESTIMATOR,
        COMPLEMENTARY_METRIC,
        "SxA",
        "MSD",
    )
    draw_constant_preference(
        axB,
        tables["comparison"],
        tables["summary"],
        PRIMARY_ESTIMATOR,
        COMPLEMENTARY_METRIC,
        "SxB",
    )
    draw_signed_slope(
        axC,
        tables["comparison"],
        tables["summary"],
        PRIMARY_ESTIMATOR,
        COMPLEMENTARY_METRIC,
        "SxC",
    )
    draw_weighting_sensitivity(
        axD,
        tables["comparison"],
        PRIMARY_METRIC,
        "SxD",
    )

    for ax, label in zip([axA, axB, axC, axD], ["A", "B", "C", "D"]):
        add_panel_label(ax, label, font_size)

    fig.subplots_adjust(wspace=0.42, hspace=0.50)
    return save_figure(
        fig,
        outdir,
        "Sx_temporal_fluctuation_scaling_diagnostics_ABCD",
        formats,
        dpi,
    )


# =============================================================================
# Standalone exports
# =============================================================================


def plot_standalones(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    dpi: int,
) -> List[Path]:
    written: List[Path] = []

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["01_variance_temporal_heatmap"])
    draw_temporal_heatmap(
        ax,
        tables["fits"],
        PRIMARY_ESTIMATOR,
        PRIMARY_METRIC,
        "01_variance_temporal_heatmap",
        "Var",
    )
    written.extend(save_figure(fig, outdir, "01_variance_temporal_heatmap", formats, dpi))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["02_variance_constant_preference"])
    draw_constant_preference(
        ax,
        tables["comparison"],
        tables["summary"],
        PRIMARY_ESTIMATOR,
        PRIMARY_METRIC,
        "02_variance_constant_preference",
    )
    written.extend(save_figure(fig, outdir, "02_variance_constant_preference", formats, dpi))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["03_variance_signed_slope"])
    draw_signed_slope(
        ax,
        tables["comparison"],
        tables["summary"],
        PRIMARY_ESTIMATOR,
        PRIMARY_METRIC,
        "03_variance_signed_slope",
    )
    written.extend(save_figure(fig, outdir, "03_variance_signed_slope", formats, dpi))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["04_msd_temporal_heatmap"])
    draw_temporal_heatmap(
        ax,
        tables["fits"],
        PRIMARY_ESTIMATOR,
        COMPLEMENTARY_METRIC,
        "04_msd_temporal_heatmap",
        "MSD",
    )
    written.extend(save_figure(fig, outdir, "04_msd_temporal_heatmap", formats, dpi))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["05_msd_constant_preference"])
    draw_constant_preference(
        ax,
        tables["comparison"],
        tables["summary"],
        PRIMARY_ESTIMATOR,
        COMPLEMENTARY_METRIC,
        "05_msd_constant_preference",
    )
    written.extend(save_figure(fig, outdir, "05_msd_constant_preference", formats, dpi))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["06_msd_signed_slope"])
    draw_signed_slope(
        ax,
        tables["comparison"],
        tables["summary"],
        PRIMARY_ESTIMATOR,
        COMPLEMENTARY_METRIC,
        "06_msd_signed_slope",
    )
    written.extend(save_figure(fig, outdir, "06_msd_signed_slope", formats, dpi))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["07_variance_slope_weighting_sensitivity"])
    draw_weighting_sensitivity(
        ax,
        tables["comparison"],
        PRIMARY_METRIC,
        "07_variance_slope_weighting_sensitivity",
    )
    written.extend(
        save_figure(fig, outdir, "07_variance_slope_weighting_sensitivity", formats, dpi)
    )

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["08_variance_temporal_profiles"])
    draw_temporal_profiles(
        ax,
        tables["fits"],
        PRIMARY_ESTIMATOR,
        PRIMARY_METRIC,
        "08_variance_temporal_profiles",
    )
    written.extend(save_figure(fig, outdir, "08_variance_temporal_profiles", formats, dpi))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES["09_msd_temporal_profiles"])
    draw_temporal_profiles(
        ax,
        tables["fits"],
        PRIMARY_ESTIMATOR,
        COMPLEMENTARY_METRIC,
        "09_msd_temporal_profiles",
    )
    written.extend(save_figure(fig, outdir, "09_msd_temporal_profiles", formats, dpi))

    return written


# =============================================================================
# CLI
# =============================================================================


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Plot ClonoDynamics Step-13 temporal fluctuation-scaling outputs.",
    )
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument(
        "--figures-dir",
        "--outdir",
        dest="figures_dir",
        type=Path,
        default=None,
        help="Default: <analysis-dir>/figures_temporal_fluctuation_scaling",
    )
    parser.add_argument(
        "--figure-set",
        choices=["all", "main", "supplementary", "standalone"],
        default="all",
    )
    parser.add_argument("--formats", default="pdf,png")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-size", type=float, default=8.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    analysis_dir = args.analysis_dir.expanduser().resolve()
    if not analysis_dir.is_dir():
        raise NotADirectoryError(
            "Step-13 analysis directory not found: {}".format(analysis_dir)
        )

    figures_dir = (
        args.figures_dir.expanduser().resolve()
        if args.figures_dir is not None
        else analysis_dir / "figures_temporal_fluctuation_scaling"
    )
    figures_dir.mkdir(parents=True, exist_ok=True)

    font_family = setup_style(args.font_size)
    formats = parse_formats(args.formats)
    tables, input_paths = load_tables(analysis_dir)

    written: List[Path] = []

    if args.figure_set in {"all", "main"}:
        written.extend(
            plot_main(tables, figures_dir, formats, args.dpi, args.font_size)
        )

    if args.figure_set in {"all", "supplementary"}:
        written.extend(
            plot_supplementary(
                tables, figures_dir, formats, args.dpi, args.font_size
            )
        )

    if args.figure_set in {"all", "standalone"}:
        written.extend(
            plot_standalones(tables, figures_dir, formats, args.dpi)
        )

    pd.DataFrame(
        {
            "file": [path.name for path in written],
            "figure": [path.stem for path in written],
            "format": [path.suffix.lstrip(".") for path in written],
        }
    ).to_csv(figures_dir / "figure_manifest.csv", index=False)

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "analysis_source": "13-temporal_fluctuation_scaling.py",
        "analysis_dir": str(analysis_dir),
        "figures_dir": str(figures_dir),
        "input_paths": input_paths,
        "figure_set": args.figure_set,
        "formats": formats,
        "dpi": int(args.dpi),
        "font_size": float(args.font_size),
        "font_family_used": font_family,
        "figure_sizes_inches": {
            key: list(value)
            for key, value in FIGURE_SIZES_INCHES.items()
        },
        "axis_config": AXIS_CONFIG,
        "analysis_free": True,
        "primary_estimator": PRIMARY_ESTIMATOR,
        "sensitivity_estimator": SENSITIVITY_ESTIMATOR,
        "primary_metric": PRIMARY_METRIC,
        "complementary_metric": COMPLEMENTARY_METRIC,
        "temporal_models_refitted": False,
        "joint_subject_bootstrap_rerun": False,
        "temporal_core_redefined": False,
        "step12_pseudo_baseline_subtracted": False,
    }
    (figures_dir / "plot_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("[DONE] Step-13 plotting complete")
    print("[INFO] Analysis directory: {}".format(analysis_dir))
    print("[INFO] Figures directory: {}".format(figures_dir))
    print("[INFO] {} rendered files".format(len(written)))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("[ERROR] {}".format(exc), file=sys.stderr)
        raise
