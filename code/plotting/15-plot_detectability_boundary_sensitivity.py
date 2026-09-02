#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
15-plot_detectability_boundary_sensitivity.py
=============================================

Plot-only companion to ClonoDynamics Step 15:

    15-detectability_boundary_sensitivity.py

The plotter does not recompute or reclassify any scientific quantity.
Standalone panels use a common 3.2 x 2.4 inch canvas. Typography and axis
geometry are harmonized across the Step-14 and Step-15 plotters:

    font family       Arial (with DejaVu Sans fallback)
    font size         8 pt
    axis linewidth    0.8 pt
    tick linewidth    0.8 pt
    data linewidth    1.2 pt

Multipanel composites retain larger canvases derived from the same panel unit.

Python >= 3.9
Required: numpy, pandas, matplotlib
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd


SCRIPT_VERSION = "v4-step15-harmonized-style-2026-09-01"
PRIMARY_ESTIMATOR = "transition_weighted"
PRIMARY_METRIC = "var_dx"
COMPLEMENTARY_METRIC = "msd"

# =============================================================================
# Shared visual specification
# =============================================================================

DEFAULT_PANEL_SIZE_INCHES: Tuple[float, float] = (3.2, 2.4)
DEFAULT_FONT_FAMILY = "Arial"
DEFAULT_FONT_SIZE = 8.0
AXES_LINEWIDTH = 0.8
TICK_LINEWIDTH = 0.8
TICK_LENGTH = 2.5
DATA_LINEWIDTH = 1.2
REFERENCE_LINEWIDTH = 0.8
ERRORBAR_LINEWIDTH = 0.8
MARKER_SIZE = 2.8
ERRORBAR_MARKER_SIZE = 3.0
ERRORBAR_CAPSIZE = 2.0

FIGURE_SIZES_INCHES: Dict[str, Tuple[float, float]] = {
    "00_main_detectability_boundary_sensitivity_ABCD": (7.2, 5.6),
    "Sx_detectability_boundary_decomposition_ABCD": (7.2, 5.6),
    "Sy_detectability_boundary_robustness_ABCD": (7.2, 5.6),
    "01_crossing_burden_by_lag": DEFAULT_PANEL_SIZE_INCHES,
    "02_crossing_imbalance_by_lag": DEFAULT_PANEL_SIZE_INCHES,
    "03_forward_TT_vs_T0PLUS": DEFAULT_PANEL_SIZE_INCHES,
    "04_forward_T0PLUS_minus_TT": DEFAULT_PANEL_SIZE_INCHES,
    "05_forward_TT_TF_T0PLUS_decomposition": DEFAULT_PANEL_SIZE_INCHES,
    "06_forward_detectability_mixture_weights": DEFAULT_PANEL_SIZE_INCHES,
    "07_forward_weighted_mixture_contributions": DEFAULT_PANEL_SIZE_INCHES,
    "08_variance_NON_FF_minus_TT": DEFAULT_PANEL_SIZE_INCHES,
    "09_msd_NON_FF_minus_TT": DEFAULT_PANEL_SIZE_INCHES,
    "10_variance_temporal_slopes_TT_vs_NON_FF": DEFAULT_PANEL_SIZE_INCHES,
    "11_msd_temporal_slopes_TT_vs_NON_FF": DEFAULT_PANEL_SIZE_INCHES,
    "12_variance_constant_preference_TT_vs_NON_FF": DEFAULT_PANEL_SIZE_INCHES,
    "13_observation_class_composition_by_dt": DEFAULT_PANEL_SIZE_INCHES,
    "14_fluctuation_mode_support": DEFAULT_PANEL_SIZE_INCHES,
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


AXIS_CONFIG: Dict[str, Dict[str, object]] = {
    "MainA": _axis_cfg(ylim=(0.0, 1.02)),
    "MainB": _axis_cfg(),
    "MainC": _axis_cfg(),
    "MainD": _axis_cfg(),
    "SxA": _axis_cfg(ylim=(-1.02, 1.02)),
    "SxB": _axis_cfg(),
    "SxC": _axis_cfg(ylim=(0.0, 1.02)),
    "SxD": _axis_cfg(),
    "SyA": _axis_cfg(),
    "SyB": _axis_cfg(),
    "SyC": _axis_cfg(),
    "SyD": _axis_cfg(ylim=(0.0, 1.02)),
    "01_crossing_burden_by_lag": _axis_cfg(ylim=(0.0, 1.02)),
    "02_crossing_imbalance_by_lag": _axis_cfg(ylim=(-1.02, 1.02)),
    "03_forward_TT_vs_T0PLUS": _axis_cfg(),
    "04_forward_T0PLUS_minus_TT": _axis_cfg(),
    "05_forward_TT_TF_T0PLUS_decomposition": _axis_cfg(),
    "06_forward_detectability_mixture_weights": _axis_cfg(ylim=(0.0, 1.02)),
    "07_forward_weighted_mixture_contributions": _axis_cfg(),
    "08_variance_NON_FF_minus_TT": _axis_cfg(),
    "09_msd_NON_FF_minus_TT": _axis_cfg(),
    "10_variance_temporal_slopes_TT_vs_NON_FF": _axis_cfg(),
    "11_msd_temporal_slopes_TT_vs_NON_FF": _axis_cfg(),
    "12_variance_constant_preference_TT_vs_NON_FF": _axis_cfg(
        ylim=(0.0, 1.02)
    ),
    "13_observation_class_composition_by_dt": _axis_cfg(
        ylim=(0.0, 1.02), xtick_step=1.0
    ),
    "14_fluctuation_mode_support": _axis_cfg(xtick_step=1.0),
}


REQUIRED_TABLES = {
    "crossing": "05_crossing_metrics_by_dt_bin",
    "forward_curves": "08_forward_TT_vs_T0plus_curves",
    "forward_difference": "09_forward_T0plus_minus_TT_by_bin",
    "forward_decomposition": "10b_forward_TT_TF_decomposition_by_bin",
    "forward_decomposition_summary": "10c_forward_decomposition_summary",
}

OPTIONAL_TABLES = {
    "class_composition": "02_class_composition_by_dt",
    "fluctuation_support": "03_fluctuation_mode_support_by_dt",
    "fluctuation_curves": "13_combined_step11_fluctuation_curves",
    "fluctuation_difference": "14_fluctuation_mode_difference_vs_TT",
    "fluctuation_summary": "15_fluctuation_mode_robustness_summary",
    "temporal_core": "18_temporal_scaling_mode_comparison_core",
    "temporal_summary": "19_temporal_scaling_sensitivity_summary",
    "step15_summary": "20_step15_summary",
}


# =============================================================================
# I/O and generic helpers
# =============================================================================


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def find_table(directory: Path, stem: str, required: bool) -> Optional[Path]:
    for suffix in (".csv", ".parquet", ".pq"):
        path = directory / (stem + suffix)
        if path.exists():
            return path
    if required:
        raise FileNotFoundError("Required Step-15 table not found: {}".format(stem))
    return None


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() in {".parquet", ".pq"}:
        try:
            return pd.read_parquet(path)
        except Exception:
            fallback = path.with_suffix(".csv")
            if fallback.exists():
                return pd.read_csv(fallback)
            raise
    raise ValueError("Unsupported table format: {}".format(path))


def load_tables(analysis_dir: Path) -> Tuple[Dict[str, pd.DataFrame], Dict[str, str]]:
    tables: Dict[str, pd.DataFrame] = {}
    paths: Dict[str, str] = {}

    for key, stem in REQUIRED_TABLES.items():
        path = find_table(analysis_dir, stem, True)
        assert path is not None
        frame = read_table(path)
        if frame.empty:
            raise ValueError("Required Step-15 table is empty: {}".format(path))
        tables[key] = frame
        paths[key] = str(path)

    for key, stem in OPTIONAL_TABLES.items():
        path = find_table(analysis_dir, stem, False)
        if path is None:
            tables[key] = pd.DataFrame()
        else:
            tables[key] = read_table(path)
            paths[key] = str(path)

    return tables, paths


def require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(
            "{} missing required columns {}. Available: {}".format(
                label, missing, list(frame.columns)
            )
        )


def numeric(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    out = frame.copy()
    for column in columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out.replace([np.inf, -np.inf], np.nan)


def resolve_font_family(requested_family: str) -> str:
    requested = str(requested_family or DEFAULT_FONT_FAMILY).strip()
    try:
        from matplotlib import font_manager

        font_manager.findfont(requested, fallback_to_default=False)
        return requested
    except Exception:
        fallback = "DejaVu Sans"
        eprint(
            "[WARN] Font {!r} not found; using {}.".format(
                requested, fallback
            )
        )
        return fallback


def setup_style(font_family: str, font_size: float) -> str:
    family = resolve_font_family(font_family)
    size = float(font_size)

    plt.rcParams.update(
        {
            "font.family": family,
            "font.size": size,
            "axes.labelsize": size,
            "axes.titlesize": size,
            "xtick.labelsize": size,
            "ytick.labelsize": size,
            "legend.fontsize": size,
            "legend.title_fontsize": size,
            "axes.linewidth": AXES_LINEWIDTH,
            "xtick.major.width": TICK_LINEWIDTH,
            "ytick.major.width": TICK_LINEWIDTH,
            "xtick.minor.width": TICK_LINEWIDTH,
            "ytick.minor.width": TICK_LINEWIDTH,
            "xtick.major.size": TICK_LENGTH,
            "ytick.major.size": TICK_LENGTH,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    return family


def style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(AXES_LINEWIDTH)
    ax.tick_params(
        direction="out",
        width=TICK_LINEWIDTH,
        length=TICK_LENGTH,
    )


def apply_axis_config(ax: plt.Axes, key: str) -> None:
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


def add_panel_label(ax: plt.Axes, label: str, font_size: float) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=max(float(font_size) + 1.0, 8.0),
        fontweight="bold",
        va="top",
        ha="left",
    )


def parse_formats(text: str) -> List[str]:
    values = [
        value.strip().lower().lstrip(".")
        for value in str(text).split(",")
        if value.strip()
    ]
    allowed = {"pdf", "png", "svg"}
    invalid = [value for value in values if value not in allowed]
    if invalid:
        raise ValueError("Unsupported figure format(s): {}".format(invalid))

    out: List[str] = []
    for value in values:
        if value not in out:
            out.append(value)
    return out or ["pdf", "png"]


def save_figure(
    fig: plt.Figure,
    outdir: Path,
    stem: str,
    formats: Sequence[str],
    dpi: int,
) -> List[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    # Fit labels, ticks and legends inside one-panel canvases without changing
    # their physical size. bbox_inches="tight" is deliberately avoided because
    # it would make nominally identical panels acquire different export sizes.
    if fig.get_figheight() <= DEFAULT_PANEL_SIZE_INCHES[1] + 1e-9:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            fig.tight_layout(pad=0.55)

    for fmt in formats:
        path = outdir / "{}.{}".format(stem, fmt)
        kwargs = {}
        if fmt == "png":
            kwargs["dpi"] = int(dpi)
        fig.savefig(path, **kwargs)
        written.append(path)
        print("[SAVED] {}".format(path))

    plt.close(fig)
    return written


def reference_horizontal(ax: plt.Axes, value: float) -> None:
    ax.axhline(
        value,
        linewidth=REFERENCE_LINEWIDTH,
        linestyle="--",
        color=ax.spines["bottom"].get_edgecolor(),
    )


def show_unavailable(ax: plt.Axes, text: str) -> None:
    ax.axis("off")
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center")


def dt_label(value: int) -> str:
    return r"$\Delta t={}$".format(int(value))


# =============================================================================
# Validation / normalized subsets
# =============================================================================


def validate_required_tables(tables: Dict[str, pd.DataFrame]) -> None:
    require_columns(
        tables["crossing"],
        [
            "dt",
            "bin_mid",
            "n_subjects",
            "n_NON_FF",
            "n_crossing",
            "crossing_burden_non_ff",
            "crossing_imbalance",
        ],
        "05_crossing_metrics_by_dt_bin",
    )
    require_columns(
        tables["forward_curves"],
        ["transition_mode", "bin", "x_center", "mean_dx"],
        "08_forward_TT_vs_T0plus_curves",
    )
    require_columns(
        tables["forward_difference"],
        [
            "bin",
            "x_center",
            "delta_mean_dx_T0PLUS_minus_TT",
            "delta_mean_dx_ci_low",
            "delta_mean_dx_ci_high",
        ],
        "09_forward_T0plus_minus_TT_by_bin",
    )
    require_columns(
        tables["forward_decomposition"],
        [
            "bin",
            "x_center",
            "fold",
            "p_TT_given_T0_detectable",
            "p_TF_given_T0_detectable",
            "mean_dx_TT",
            "mean_dx_TF_derived",
            "mean_dx_T0PLUS_direct",
            "TT_contribution_to_T0PLUS_mean",
            "TF_contribution_to_T0PLUS_mean",
            "mean_dx_T0PLUS_reconstructed",
            "reconstruction_residual",
        ],
        "10b_forward_TT_TF_decomposition_by_bin",
    )
    require_columns(
        tables["forward_decomposition_summary"],
        [
            "mean_p_TF_given_T0_detectable",
            "mean_dx_TF_across_bins",
            "mean_abs_TT_contribution",
            "mean_abs_TF_contribution",
            "max_abs_reconstruction_error",
            "pearson_direct_vs_reconstructed",
        ],
        "10c_forward_decomposition_summary",
    )


def combined_decomposition(frame: pd.DataFrame) -> pd.DataFrame:
    data = frame.loc[frame["fold"].astype(str) == "combined"].copy()
    data = numeric(
        data,
        [
            "bin",
            "x_center",
            "p_TT_given_T0_detectable",
            "p_TF_given_T0_detectable",
            "mean_dx_TT",
            "mean_dx_TF_derived",
            "mean_dx_T0PLUS_direct",
            "TT_contribution_to_T0PLUS_mean",
            "TF_contribution_to_T0PLUS_mean",
            "mean_dx_T0PLUS_reconstructed",
            "reconstruction_residual",
        ],
    ).sort_values("bin")
    if data.empty:
        raise ValueError("Forward decomposition has no fold='combined' rows.")
    return data


def temporal_subset(
    frame: pd.DataFrame,
    metric: str,
    estimator: str,
) -> pd.DataFrame:
    require_columns(
        frame,
        [
            "transition_mode",
            "metric",
            "estimator",
            "bin_id",
            "bin_mid",
            "signed_linear_slope_median",
            "signed_linear_slope_q_low",
            "signed_linear_slope_q_high",
        ],
        "18_temporal_scaling_mode_comparison_core",
    )
    data = frame.loc[
        (frame["metric"].astype(str) == metric)
        & (frame["estimator"].astype(str) == estimator)
        & (frame["transition_mode"].astype(str).isin(["TT", "NON_FF"]))
    ].copy()
    data = numeric(
        data,
        [
            "bin_id",
            "bin_mid",
            "signed_linear_slope_median",
            "signed_linear_slope_q_low",
            "signed_linear_slope_q_high",
            "constant_preference_fraction",
        ],
    )
    if data.empty:
        raise ValueError(
            "No temporal-core rows for metric={} estimator={}.".format(
                metric, estimator
            )
        )
    return data.sort_values(["transition_mode", "bin_id"])


# =============================================================================
# Crossing descriptors
# =============================================================================


def draw_crossing_burden(
    ax: plt.Axes,
    crossing: pd.DataFrame,
    axis_key: str,
    min_nonff: int,
    min_subjects: int,
    show_legend: bool = True,
) -> None:
    data = numeric(
        crossing,
        ["dt", "bin_mid", "n_subjects", "n_NON_FF", "crossing_burden_non_ff"],
    )
    data = data.loc[
        (data["n_NON_FF"] >= int(min_nonff))
        & (data["n_subjects"] >= int(min_subjects))
        & np.isfinite(data["bin_mid"])
        & np.isfinite(data["crossing_burden_non_ff"])
    ].copy()
    if data.empty:
        raise ValueError("No crossing-burden rows remain after display filtering.")

    markers = ["o", "s", "^", "D", "v", "P", "X"]
    for idx, dt in enumerate(
        sorted(data["dt"].dropna().astype(int).unique().tolist())
    ):
        subset = data.loc[data["dt"].astype(int) == dt].sort_values("bin_mid")
        ax.plot(
            subset["bin_mid"],
            subset["crossing_burden_non_ff"],
            marker=markers[idx % len(markers)],
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=dt_label(dt),
        )
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel("Crossing burden among NON_FF")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    if show_legend:
        ax.legend(frameon=False, title="Lag")


def draw_crossing_imbalance(
    ax: plt.Axes,
    crossing: pd.DataFrame,
    axis_key: str,
    min_crossing_events: int,
    min_subjects: int,
    show_legend: bool = True,
) -> None:
    data = numeric(
        crossing,
        ["dt", "bin_mid", "n_subjects", "n_crossing", "crossing_imbalance"],
    )
    data = data.loc[
        (data["n_crossing"] >= int(min_crossing_events))
        & (data["n_subjects"] >= int(min_subjects))
        & np.isfinite(data["bin_mid"])
        & np.isfinite(data["crossing_imbalance"])
    ].copy()
    if data.empty:
        raise ValueError("No crossing-imbalance rows remain after display filtering.")

    markers = ["o", "s", "^", "D", "v", "P", "X"]
    for idx, dt in enumerate(
        sorted(data["dt"].dropna().astype(int).unique().tolist())
    ):
        subset = data.loc[data["dt"].astype(int) == dt].sort_values("bin_mid")
        ax.plot(
            subset["bin_mid"],
            subset["crossing_imbalance"],
            marker=markers[idx % len(markers)],
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=dt_label(dt),
        )
    reference_horizontal(ax, 0.0)
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel(r"Crossing imbalance, $J$")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    if show_legend:
        ax.legend(frameon=False, title="Lag")


# =============================================================================
# Forward sensitivity / decomposition
# =============================================================================


def draw_forward_curves(
    ax: plt.Axes,
    curves: pd.DataFrame,
    axis_key: str,
) -> None:
    data = numeric(
        curves,
        ["bin", "x_center", "mean_dx", "mean_dx_ci_low", "mean_dx_ci_high"],
    )
    for mode, label, marker in (
        ("TT", "TT", "o"),
        ("T0PLUS", "T0PLUS = TT+TF", "s"),
    ):
        subset = data.loc[data["transition_mode"].astype(str) == mode].sort_values(
            "bin"
        )
        if subset.empty:
            continue
        x = subset["x_center"].to_numpy(float)
        y = subset["mean_dx"].to_numpy(float)
        valid = np.isfinite(x) & np.isfinite(y)
        line = ax.plot(
            x[valid],
            y[valid],
            marker=marker,
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=label,
        )[0]
        if {"mean_dx_ci_low", "mean_dx_ci_high"}.issubset(subset.columns):
            low = subset["mean_dx_ci_low"].to_numpy(float)
            high = subset["mean_dx_ci_high"].to_numpy(float)
            band = np.isfinite(x) & np.isfinite(low) & np.isfinite(high)
            if band.any():
                ax.fill_between(
                    x[band],
                    low[band],
                    high[band],
                    color=line.get_color(),
                    alpha=0.12,
                    linewidth=0,
                )
    reference_horizontal(ax, 0.0)
    ax.set_xlabel(r"Initial latent abundance, $x_0$")
    ax.set_ylabel("Mean latent displacement")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)


def draw_forward_difference(
    ax: plt.Axes,
    difference: pd.DataFrame,
    axis_key: str,
) -> None:
    data = numeric(
        difference,
        [
            "bin",
            "x_center",
            "delta_mean_dx_T0PLUS_minus_TT",
            "delta_mean_dx_boot_median",
            "delta_mean_dx_ci_low",
            "delta_mean_dx_ci_high",
        ],
    ).sort_values("bin")
    x = data["x_center"].to_numpy(float)
    point_col = (
        "delta_mean_dx_boot_median"
        if "delta_mean_dx_boot_median" in data.columns
        else "delta_mean_dx_T0PLUS_minus_TT"
    )
    y = data[point_col].to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(y)
    line = ax.plot(
        x[valid],
        y[valid],
        marker="o",
        markersize=MARKER_SIZE,
        linewidth=DATA_LINEWIDTH,
    )[0]
    if {"delta_mean_dx_ci_low", "delta_mean_dx_ci_high"}.issubset(data.columns):
        low = data["delta_mean_dx_ci_low"].to_numpy(float)
        high = data["delta_mean_dx_ci_high"].to_numpy(float)
        band = np.isfinite(x) & np.isfinite(low) & np.isfinite(high)
        if band.any():
            ax.fill_between(
                x[band],
                low[band],
                high[band],
                color=line.get_color(),
                alpha=0.14,
                linewidth=0,
            )
    reference_horizontal(ax, 0.0)
    ax.set_xlabel(r"Initial latent abundance, $x_0$")
    ax.set_ylabel("T0PLUS − TT\nmean displacement")
    style_axes(ax)
    apply_axis_config(ax, axis_key)


def draw_forward_decomposition(
    ax: plt.Axes,
    decomposition: pd.DataFrame,
    summary: pd.DataFrame,
    axis_key: str,
    show_reconstruction_error: bool = True,
) -> None:
    data = combined_decomposition(decomposition)
    x = data["x_center"].to_numpy(float)
    for column, label, marker, linestyle in (
        ("mean_dx_TT", "TT", "o", "-"),
        ("mean_dx_TF_derived", "Derived TF", "s", ":"),
        ("mean_dx_T0PLUS_direct", "T0PLUS", "^", "--"),
    ):
        y = data[column].to_numpy(float)
        valid = np.isfinite(x) & np.isfinite(y)
        ax.plot(
            x[valid],
            y[valid],
            marker=marker,
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            linestyle=linestyle,
            label=label,
        )
    reference_horizontal(ax, 0.0)
    ax.set_xlabel(r"Initial latent abundance, $x_0$")
    ax.set_ylabel("Mean latent displacement")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)

    if not summary.empty:
        row = summary.iloc[0]
        p_tf = pd.to_numeric(
            pd.Series([row.get("mean_p_TF_given_T0_detectable")]),
            errors="coerce",
        ).iloc[0]
        error = pd.to_numeric(
            pd.Series([row.get("max_abs_reconstruction_error")]),
            errors="coerce",
        ).iloc[0]
        parts: List[str] = []
        if np.isfinite(p_tf):
            parts.append(r"mean $P(TF\mid T_0)$={:.3f}".format(float(p_tf)))
        if show_reconstruction_error and np.isfinite(error):
            parts.append("max recon. err={:.2g}".format(float(error)))
        if parts:
            ax.set_title(
                ";  ".join(parts),
                loc="left",
                pad=6.0,
                fontsize=float(plt.rcParams["font.size"]),
                fontweight="normal",
            )


def draw_forward_weights(
    ax: plt.Axes,
    decomposition: pd.DataFrame,
    axis_key: str,
) -> None:
    data = combined_decomposition(decomposition)
    for column, label, marker in (
        ("p_TT_given_T0_detectable", "TT fraction", "o"),
        ("p_TF_given_T0_detectable", "TF fraction", "s"),
    ):
        ax.plot(
            data["x_center"],
            data[column],
            marker=marker,
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=label,
        )
    ax.set_xlabel(r"Initial latent abundance, $x_0$")
    ax.set_ylabel("Conditional transition fraction")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)


def draw_forward_contributions(
    ax: plt.Axes,
    decomposition: pd.DataFrame,
    axis_key: str,
) -> None:
    data = combined_decomposition(decomposition)
    for column, label, marker, linestyle in (
        ("TT_contribution_to_T0PLUS_mean", "TT contribution", "o", "-"),
        ("TF_contribution_to_T0PLUS_mean", "TF contribution", "s", ":"),
        ("mean_dx_T0PLUS_direct", "Direct T0PLUS", "^", "--"),
    ):
        ax.plot(
            data["x_center"],
            data[column],
            marker=marker,
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            linestyle=linestyle,
            label=label,
        )
    reference_horizontal(ax, 0.0)
    ax.set_xlabel(r"Initial latent abundance, $x_0$")
    ax.set_ylabel(r"Contribution to $E[\Delta x\mid T_0]$")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)


# =============================================================================
# Fluctuation sensitivity
# =============================================================================


def draw_fluctuation_difference(
    ax: plt.Axes,
    difference: pd.DataFrame,
    metric: str,
    axis_key: str,
    show_legend: bool = True,
) -> None:
    require_columns(
        difference,
        ["metric", "dt", "bin_mid", "difference_NON_FF_minus_TT"],
        "14_fluctuation_mode_difference_vs_TT",
    )
    data = difference.loc[difference["metric"].astype(str) == metric].copy()
    data = numeric(data, ["dt", "bin_mid", "difference_NON_FF_minus_TT"])
    if data.empty:
        show_unavailable(ax, "No {} difference output".format(metric))
        return

    markers = ["o", "s", "^", "D", "v", "P", "X"]
    for idx, dt in enumerate(
        sorted(data["dt"].dropna().astype(int).unique().tolist())
    ):
        subset = data.loc[data["dt"].astype(int) == dt].sort_values("bin_mid")
        ax.plot(
            subset["bin_mid"],
            subset["difference_NON_FF_minus_TT"],
            marker=markers[idx % len(markers)],
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=dt_label(dt),
        )
    reference_horizontal(ax, 0.0)
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    if metric == "var_dx":
        ax.set_ylabel("Excess conditional variance,\nNON_FF − TT")
    else:
        ax.set_ylabel("Excess MSD,\nNON_FF − TT")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    if show_legend:
        ax.legend(frameon=False, title="Lag")


def draw_fluctuation_curves_at_dt(
    ax: plt.Axes,
    curves: pd.DataFrame,
    metric: str,
    dt: int,
    axis_key: str,
) -> None:
    require_columns(
        curves,
        ["transition_mode", "dt", "bin_mid", metric],
        "13_combined_step11_fluctuation_curves",
    )
    data = numeric(curves, ["dt", "bin_mid", metric])
    data = data.loc[data["dt"] == int(dt)].copy()
    if data.empty:
        show_unavailable(ax, "No Step-11 sensitivity curves at dt={}".format(dt))
        return

    for mode, label, marker in (
        ("TT", "TT", "o"),
        ("NON_FF", "NON_FF = TT+TF+FT", "s"),
    ):
        subset = data.loc[data["transition_mode"].astype(str) == mode].sort_values(
            "bin_mid"
        )
        if subset.empty:
            continue
        ax.plot(
            subset["bin_mid"],
            subset[metric],
            marker=marker,
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=label,
        )
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel(r"Var($\Delta x$)" if metric == "var_dx" else "MSD")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)


# =============================================================================
# Temporal-scaling sensitivity
# =============================================================================


def draw_temporal_slopes(
    ax: plt.Axes,
    temporal_core: pd.DataFrame,
    metric: str,
    estimator: str,
    axis_key: str,
) -> None:
    data = temporal_subset(temporal_core, metric, estimator)
    specs = {"TT": (-0.025, "o"), "NON_FF": (0.025, "s")}

    for mode in ("TT", "NON_FF"):
        subset = data.loc[data["transition_mode"].astype(str) == mode].sort_values(
            "bin_id"
        )
        if subset.empty:
            continue
        offset, marker = specs[mode]
        x = subset["bin_mid"].to_numpy(float) + offset
        med = subset["signed_linear_slope_median"].to_numpy(float)
        low = subset["signed_linear_slope_q_low"].to_numpy(float)
        high = subset["signed_linear_slope_q_high"].to_numpy(float)
        valid = np.isfinite(x) & np.isfinite(med) & np.isfinite(low) & np.isfinite(high)
        yerr = np.vstack(
            [
                np.maximum(med - low, 0.0),
                np.maximum(high - med, 0.0),
            ]
        )
        ax.errorbar(
            x[valid],
            med[valid],
            yerr=yerr[:, valid],
            fmt=marker,
            markersize=ERRORBAR_MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            elinewidth=ERRORBAR_LINEWIDTH,
            capsize=ERRORBAR_CAPSIZE,
            label=mode,
        )
    reference_horizontal(ax, 0.0)
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel(r"Signed-linear slope, $D$")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)


def draw_constant_preference(
    ax: plt.Axes,
    temporal_core: pd.DataFrame,
    metric: str,
    estimator: str,
    axis_key: str,
) -> None:
    data = temporal_subset(temporal_core, metric, estimator)
    if "constant_preference_fraction" not in data.columns:
        show_unavailable(ax, "Constant-preference output unavailable")
        return

    for mode, marker in (("TT", "o"), ("NON_FF", "s")):
        subset = data.loc[data["transition_mode"].astype(str) == mode].sort_values(
            "bin_id"
        )
        ax.plot(
            subset["bin_mid"],
            subset["constant_preference_fraction"],
            marker=marker,
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=mode,
        )
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel("Constant-model preference fraction")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)


# =============================================================================
# Structural diagnostics
# =============================================================================


def draw_class_composition(
    ax: plt.Axes,
    table: pd.DataFrame,
    axis_key: str,
) -> None:
    require_columns(
        table,
        ["dt", "obs_class", "fraction"],
        "02_class_composition_by_dt",
    )
    data = numeric(table, ["dt", "fraction"])
    classes = [
        value
        for value in ("TT", "TF", "FT", "FF")
        if value in set(data["obs_class"].astype(str))
    ]
    for obs_class in classes:
        subset = data.loc[data["obs_class"].astype(str) == obs_class].sort_values(
            "dt"
        )
        ax.plot(
            subset["dt"],
            subset["fraction"],
            marker="o",
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=obs_class,
        )
    ax.set_xlabel(r"Temporal lag, $\Delta t$")
    ax.set_ylabel("Transition-class fraction")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)


def draw_mode_support(
    ax: plt.Axes,
    table: pd.DataFrame,
    axis_key: str,
) -> None:
    require_columns(
        table,
        ["transition_mode", "dt", "n_transitions"],
        "03_fluctuation_mode_support_by_dt",
    )
    data = numeric(table, ["dt", "n_transitions"])
    for mode, marker in (("TT", "o"), ("NON_FF", "s")):
        subset = data.loc[data["transition_mode"].astype(str) == mode].sort_values(
            "dt"
        )
        ax.plot(
            subset["dt"],
            subset["n_transitions"],
            marker=marker,
            markersize=MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            label=mode,
        )
    ax.set_xlabel(r"Temporal lag, $\Delta t$")
    ax.set_ylabel("Transitions")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False)


# =============================================================================
# Figure assembly
# =============================================================================


def downstream_available(tables: Dict[str, pd.DataFrame]) -> bool:
    return (
        not tables["fluctuation_difference"].empty
        and not tables["temporal_core"].empty
    )


def plot_main(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    args: argparse.Namespace,
) -> List[Path]:
    stem = "00_main_detectability_boundary_sensitivity_ABCD"
    fig, axes = plt.subplots(2, 2, figsize=FIGURE_SIZES_INCHES[stem])
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    draw_crossing_burden(
        ax_a,
        tables["crossing"],
        "MainA",
        args.display_min_nonff,
        args.display_min_subjects,
    )
    draw_forward_decomposition(
        ax_b,
        tables["forward_decomposition"],
        tables["forward_decomposition_summary"],
        "MainB",
        show_reconstruction_error=False,
    )

    if not tables["fluctuation_difference"].empty:
        draw_fluctuation_difference(
            ax_c,
            tables["fluctuation_difference"],
            PRIMARY_METRIC,
            "MainC",
        )
    else:
        show_unavailable(
            ax_c,
            "Run Step 15 with --run-step11\nfor fluctuation sensitivity",
        )

    if not tables["temporal_core"].empty:
        draw_temporal_slopes(
            ax_d,
            tables["temporal_core"],
            PRIMARY_METRIC,
            PRIMARY_ESTIMATOR,
            "MainD",
        )
    else:
        show_unavailable(
            ax_d,
            "Run Step 15 with --run-step13\nfor temporal-scaling sensitivity",
        )

    for ax, label in zip(axes.ravel(), "ABCD"):
        add_panel_label(ax, label, args.font_size)

    fig.subplots_adjust(wspace=0.46, hspace=0.46)
    return save_figure(fig, outdir, stem, formats, args.dpi)


def plot_decomposition(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    args: argparse.Namespace,
) -> List[Path]:
    stem = "Sx_detectability_boundary_decomposition_ABCD"
    fig, axes = plt.subplots(2, 2, figsize=FIGURE_SIZES_INCHES[stem])
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    draw_crossing_imbalance(
        ax_a,
        tables["crossing"],
        "SxA",
        args.display_min_crossing,
        args.display_min_subjects,
    )
    draw_forward_difference(ax_b, tables["forward_difference"], "SxB")
    draw_forward_weights(ax_c, tables["forward_decomposition"], "SxC")
    draw_forward_contributions(ax_d, tables["forward_decomposition"], "SxD")

    for ax, label in zip(axes.ravel(), "ABCD"):
        add_panel_label(ax, label, args.font_size)

    fig.subplots_adjust(wspace=0.46, hspace=0.46)
    return save_figure(fig, outdir, stem, formats, args.dpi)


def plot_robustness(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    args: argparse.Namespace,
) -> List[Path]:
    if not downstream_available(tables):
        warnings.warn(
            "Step-11/13 sensitivity outputs are incomplete; "
            "Sy robustness composite skipped."
        )
        return []

    stem = "Sy_detectability_boundary_robustness_ABCD"
    fig, axes = plt.subplots(2, 2, figsize=FIGURE_SIZES_INCHES[stem])
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    if not tables["fluctuation_curves"].empty:
        available_dt = pd.to_numeric(
            tables["fluctuation_curves"]["dt"], errors="coerce"
        ).dropna().astype(int)
        display_dt = (
            args.display_dt
            if args.display_dt is not None
            else int(available_dt.min())
        )
        draw_fluctuation_curves_at_dt(
            ax_a,
            tables["fluctuation_curves"],
            PRIMARY_METRIC,
            display_dt,
            "SyA",
        )
    else:
        show_unavailable(ax_a, "Step-11 combined curves unavailable")

    draw_fluctuation_difference(
        ax_b,
        tables["fluctuation_difference"],
        COMPLEMENTARY_METRIC,
        "SyB",
    )
    draw_temporal_slopes(
        ax_c,
        tables["temporal_core"],
        COMPLEMENTARY_METRIC,
        PRIMARY_ESTIMATOR,
        "SyC",
    )
    draw_constant_preference(
        ax_d,
        tables["temporal_core"],
        PRIMARY_METRIC,
        PRIMARY_ESTIMATOR,
        "SyD",
    )

    for ax, label in zip(axes.ravel(), "ABCD"):
        add_panel_label(ax, label, args.font_size)

    fig.subplots_adjust(wspace=0.46, hspace=0.46)
    return save_figure(fig, outdir, stem, formats, args.dpi)


def plot_standalones(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    args: argparse.Namespace,
) -> List[Path]:
    written: List[Path] = []

    specs = [
        (
            "01_crossing_burden_by_lag",
            lambda ax: draw_crossing_burden(
                ax,
                tables["crossing"],
                "01_crossing_burden_by_lag",
                args.display_min_nonff,
                args.display_min_subjects,
            ),
        ),
        (
            "02_crossing_imbalance_by_lag",
            lambda ax: draw_crossing_imbalance(
                ax,
                tables["crossing"],
                "02_crossing_imbalance_by_lag",
                args.display_min_crossing,
                args.display_min_subjects,
            ),
        ),
        (
            "03_forward_TT_vs_T0PLUS",
            lambda ax: draw_forward_curves(
                ax,
                tables["forward_curves"],
                "03_forward_TT_vs_T0PLUS",
            ),
        ),
        (
            "04_forward_T0PLUS_minus_TT",
            lambda ax: draw_forward_difference(
                ax,
                tables["forward_difference"],
                "04_forward_T0PLUS_minus_TT",
            ),
        ),
        (
            "05_forward_TT_TF_T0PLUS_decomposition",
            lambda ax: draw_forward_decomposition(
                ax,
                tables["forward_decomposition"],
                tables["forward_decomposition_summary"],
                "05_forward_TT_TF_T0PLUS_decomposition",
            ),
        ),
        (
            "06_forward_detectability_mixture_weights",
            lambda ax: draw_forward_weights(
                ax,
                tables["forward_decomposition"],
                "06_forward_detectability_mixture_weights",
            ),
        ),
        (
            "07_forward_weighted_mixture_contributions",
            lambda ax: draw_forward_contributions(
                ax,
                tables["forward_decomposition"],
                "07_forward_weighted_mixture_contributions",
            ),
        ),
    ]

    for stem, function in specs:
        fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
        function(ax)
        written.extend(save_figure(fig, outdir, stem, formats, args.dpi))

    if not tables["fluctuation_difference"].empty:
        for stem, metric in (
            ("08_variance_NON_FF_minus_TT", PRIMARY_METRIC),
            ("09_msd_NON_FF_minus_TT", COMPLEMENTARY_METRIC),
        ):
            fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
            draw_fluctuation_difference(
                ax,
                tables["fluctuation_difference"],
                metric,
                stem,
            )
            written.extend(save_figure(fig, outdir, stem, formats, args.dpi))

    if not tables["temporal_core"].empty:
        for stem, metric in (
            ("10_variance_temporal_slopes_TT_vs_NON_FF", PRIMARY_METRIC),
            ("11_msd_temporal_slopes_TT_vs_NON_FF", COMPLEMENTARY_METRIC),
        ):
            fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
            draw_temporal_slopes(
                ax,
                tables["temporal_core"],
                metric,
                PRIMARY_ESTIMATOR,
                stem,
            )
            written.extend(save_figure(fig, outdir, stem, formats, args.dpi))

        stem = "12_variance_constant_preference_TT_vs_NON_FF"
        fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
        draw_constant_preference(
            ax,
            tables["temporal_core"],
            PRIMARY_METRIC,
            PRIMARY_ESTIMATOR,
            stem,
        )
        written.extend(save_figure(fig, outdir, stem, formats, args.dpi))

    if not tables["class_composition"].empty:
        stem = "13_observation_class_composition_by_dt"
        fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
        draw_class_composition(ax, tables["class_composition"], stem)
        written.extend(save_figure(fig, outdir, stem, formats, args.dpi))

    if not tables["fluctuation_support"].empty:
        stem = "14_fluctuation_mode_support"
        fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
        draw_mode_support(ax, tables["fluctuation_support"], stem)
        written.extend(save_figure(fig, outdir, stem, formats, args.dpi))

    return written


# =============================================================================
# CLI / main
# =============================================================================


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Plot ClonoDynamics Step-15 detectability-boundary sensitivity outputs."
        ),
    )
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument(
        "--figures-dir",
        "--outdir",
        dest="figures_dir",
        type=Path,
        default=None,
        help=(
            "Default: <analysis-dir>/figures_detectability_boundary_sensitivity"
        ),
    )
    parser.add_argument(
        "--figure-set",
        choices=["all", "main", "decomposition", "robustness", "standalone"],
        default="all",
    )
    parser.add_argument("--formats", default="pdf,png")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-family", default=DEFAULT_FONT_FAMILY)
    parser.add_argument("--font-size", type=float, default=DEFAULT_FONT_SIZE)
    parser.add_argument(
        "--display-min-nonff",
        type=int,
        default=30,
        help="Display-only minimum NON_FF events for crossing-burden curves.",
    )
    parser.add_argument(
        "--display-min-crossing",
        type=int,
        default=10,
        help="Display-only minimum TF+FT events for crossing-imbalance curves.",
    )
    parser.add_argument(
        "--display-min-subjects",
        type=int,
        default=2,
        help="Display-only minimum subjects for crossing panels.",
    )
    parser.add_argument(
        "--display-dt",
        type=int,
        default=None,
        help=(
            "Lag used in the TT-vs-NON_FF Var curve panel; default uses the "
            "smallest available lag."
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    analysis_dir = args.analysis_dir.expanduser().resolve()
    if not analysis_dir.is_dir():
        raise NotADirectoryError(analysis_dir)

    figures_dir = (
        args.figures_dir.expanduser().resolve()
        if args.figures_dir is not None
        else analysis_dir / "figures_detectability_boundary_sensitivity"
    )
    figures_dir.mkdir(parents=True, exist_ok=True)

    formats = parse_formats(args.formats)
    font_used = setup_style(args.font_family, args.font_size)
    tables, input_paths = load_tables(analysis_dir)
    validate_required_tables(tables)

    written: List[Path] = []
    if args.figure_set in {"all", "main"}:
        written.extend(plot_main(tables, figures_dir, formats, args))
    if args.figure_set in {"all", "decomposition"}:
        written.extend(plot_decomposition(tables, figures_dir, formats, args))
    if args.figure_set in {"all", "robustness"}:
        written.extend(plot_robustness(tables, figures_dir, formats, args))
    if args.figure_set in {"all", "standalone"}:
        written.extend(plot_standalones(tables, figures_dir, formats, args))

    pd.DataFrame(
        {
            "file": [path.name for path in written],
            "figure": [path.stem for path in written],
            "format": [path.suffix.lstrip(".") for path in written],
        }
    ).to_csv(figures_dir / "figure_manifest.csv", index=False)

    analysis_config = {}
    config_path = analysis_dir / "00_run_config.json"
    if config_path.exists():
        try:
            analysis_config = json.loads(config_path.read_text(encoding="utf-8"))
        except Exception:
            analysis_config = {}

    plot_config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "analysis_source": "15-detectability_boundary_sensitivity.py",
        "analysis_dir": str(analysis_dir),
        "figures_dir": str(figures_dir),
        "input_paths": input_paths,
        "figure_set": args.figure_set,
        "formats": formats,
        "dpi": int(args.dpi),
        "font_family_requested": str(args.font_family),
        "font_family_used": font_used,
        "font_size": float(args.font_size),
        "default_panel_size_inches": list(DEFAULT_PANEL_SIZE_INCHES),
        "axes_linewidth": AXES_LINEWIDTH,
        "tick_linewidth": TICK_LINEWIDTH,
        "data_linewidth": DATA_LINEWIDTH,
        "display_min_nonff": int(args.display_min_nonff),
        "display_min_crossing": int(args.display_min_crossing),
        "display_min_subjects": int(args.display_min_subjects),
        "display_dt": args.display_dt,
        "figure_sizes_inches": {
            key: list(value) for key, value in FIGURE_SIZES_INCHES.items()
        },
        "axis_config": AXIS_CONFIG,
        "analysis_free": True,
        "detectability_labels_reclassified": False,
        "alpha_sweep_performed": False,
        "forward_estimand": "TT vs T0PLUS=TT+TF",
        "fluctuation_estimand": "TT vs NON_FF=TT+TF+FT",
        "primary_metric": PRIMARY_METRIC,
        "primary_estimator": PRIMARY_ESTIMATOR,
        "analysis_script_version": analysis_config.get("script_version"),
    }
    (figures_dir / "plot_config.json").write_text(
        json.dumps(plot_config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("[DONE] Step-15 plotting complete")
    print("[INFO] Analysis directory: {}".format(analysis_dir))
    print("[INFO] Figures directory: {}".format(figures_dir))
    print("[INFO] {} rendered files".format(len(written)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
