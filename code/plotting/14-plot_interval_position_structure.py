#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
14-plot_interval_position_structure.py
======================================

Plot-only companion to ClonoDynamics Step 14:

    14-interval_position_structure.py

The plotter reads finalized Step-14 outputs and never recomputes scientific
quantities. Standalone panels use a common 3.2 x 2.4 inch canvas. Typography
and axis geometry are harmonized across the Step-14 and Step-15 plotters:

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


SCRIPT_VERSION = "v6-step14-harmonized-style-2026-09-01"
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

# Every one-panel output has exactly the same requested canvas size. Composite
# outputs remain larger because they contain two or four panels.
FIGURE_SIZES_INCHES: Dict[str, Tuple[float, float]] = {
    "00_main_interval_position_structure_ABCD": (7.2, 5.6),
    "Sx_interval_composition_diagnostics_ABCD": (7.2, 5.6),
    "01_var_interval_position_q_heatmap": DEFAULT_PANEL_SIZE_INCHES,
    "02_var_anchor_available_slopes": DEFAULT_PANEL_SIZE_INCHES,
    "03_var_anchor_matched_slopes": DEFAULT_PANEL_SIZE_INCHES,
    "04_msd_interval_position_q_heatmap": DEFAULT_PANEL_SIZE_INCHES,
    "05_msd_anchor_available_slopes": DEFAULT_PANEL_SIZE_INCHES,
    "06_msd_anchor_matched_slopes": DEFAULT_PANEL_SIZE_INCHES,
    "07_var_secondary_interval_q_heatmaps": (6.4, 2.4),
    "08_interval_transition_counts": DEFAULT_PANEL_SIZE_INCHES,
    "09_interval_subject_counts": DEFAULT_PANEL_SIZE_INCHES,
    "10_anchor_complete_subjects": DEFAULT_PANEL_SIZE_INCHES,
    "11_matched_var_subject_support": DEFAULT_PANEL_SIZE_INCHES,
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
    "MainA": _axis_cfg(),
    "MainB": _axis_cfg(),
    "MainC": _axis_cfg(),
    "MainD": _axis_cfg(),
    "SxA": _axis_cfg(xtick_step=1.0),
    "SxB": _axis_cfg(xtick_step=1.0),
    "SxC": _axis_cfg(),
    "SxD": _axis_cfg(),
    "01_var_interval_position_q_heatmap": _axis_cfg(),
    "02_var_anchor_available_slopes": _axis_cfg(),
    "03_var_anchor_matched_slopes": _axis_cfg(),
    "04_msd_interval_position_q_heatmap": _axis_cfg(),
    "05_msd_anchor_available_slopes": _axis_cfg(),
    "06_msd_anchor_matched_slopes": _axis_cfg(),
    "07a_var_calendar_slope_q": _axis_cfg(),
    "07b_var_alternation_q": _axis_cfg(),
    "08_interval_transition_counts": _axis_cfg(xtick_step=1.0),
    "09_interval_subject_counts": _axis_cfg(xtick_step=1.0),
    "10_anchor_complete_subjects": _axis_cfg(),
    "11_matched_var_subject_support": _axis_cfg(),
}


REQUIRED_TABLES = {
    "core_bins": "02_core_bin_definitions.csv",
    "permutation_tests": "06_interval_position_permutation_tests.csv",
    "anchored_slopes": "09_anchored_signed_slope_bootstrap.csv",
    "analysis_summary": "10_analysis_summary.csv",
}

OPTIONAL_TABLES = {
    "composition": "03_interval_composition_by_pair.csv",
    "permutation_summary": "07_interval_position_test_summary.csv",
    "anchored_profiles": "08_anchored_temporal_profiles.csv",
    "matching": "11_anchor_subject_matching_diagnostic.csv",
}


# =============================================================================
# I/O and style
# =============================================================================


def eprint(*args, **kwargs) -> None:
    print(*args, file=sys.stderr, **kwargs)


def read_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError("Required Step-14 table not found: {}".format(path))
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError("Required Step-14 table is empty: {}".format(path))
    return frame


def read_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def load_tables(analysis_dir: Path) -> Tuple[Dict[str, pd.DataFrame], Dict[str, str]]:
    tables: Dict[str, pd.DataFrame] = {}
    paths: Dict[str, str] = {}

    for key, filename in REQUIRED_TABLES.items():
        path = analysis_dir / filename
        tables[key] = read_required(path)
        paths[key] = str(path)

    for key, filename in OPTIONAL_TABLES.items():
        path = analysis_dir / filename
        tables[key] = read_optional(path)
        if path.exists():
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


def bool_mask(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).to_numpy(dtype=bool)
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "yes", "y", "t"])
        .to_numpy(dtype=bool)
    )


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
        1.20,
        label,
        transform=ax.transAxes,
        fontsize=max(float(font_size) + 1.0, 8.0),
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
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
        raise ValueError("Unsupported format(s): {}".format(invalid))

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


# =============================================================================
# Primary interval-position test
# =============================================================================


def summary_row(summary: pd.DataFrame, metric: str) -> pd.Series:
    require_columns(summary, ["metric"], "10_analysis_summary.csv")
    subset = summary.loc[summary["metric"].astype(str) == str(metric)]
    if len(subset) != 1:
        raise ValueError(
            "Expected one Step-14 summary row for metric {!r}; found {}".format(
                metric, len(subset)
            )
        )
    return subset.iloc[0]


def prepare_q_matrix(
    tests: pd.DataFrame,
    core_bins: pd.DataFrame,
    metric: str,
    q_col: str,
) -> Tuple[np.ndarray, List[int], np.ndarray]:
    require_columns(
        tests,
        ["metric", "dt", "bin_id", "test_status", q_col],
        "06_interval_position_permutation_tests.csv",
    )
    require_columns(
        core_bins,
        ["bin_id", "bin_mid"],
        "02_core_bin_definitions.csv",
    )

    data = tests.loc[
        (tests["metric"].astype(str) == str(metric))
        & (tests["test_status"].astype(str) == "ok")
    ].copy()

    data[q_col] = pd.to_numeric(data[q_col], errors="coerce")
    data["dt"] = pd.to_numeric(data["dt"], errors="coerce")
    data["bin_id"] = pd.to_numeric(data["bin_id"], errors="coerce")
    data = data.dropna(subset=[q_col, "dt", "bin_id"])

    if data.empty:
        raise ValueError(
            "No valid interval-position tests for metric={!r}, q={!r}".format(
                metric, q_col
            )
        )

    dt_values = sorted(data["dt"].astype(int).unique().tolist())
    bins = core_bins.copy()
    bins["bin_id"] = pd.to_numeric(bins["bin_id"], errors="coerce")
    bins["bin_mid"] = pd.to_numeric(bins["bin_mid"], errors="coerce")
    bins = bins.dropna(subset=["bin_id", "bin_mid"]).sort_values("bin_id")

    bin_ids = bins["bin_id"].astype(int).tolist()
    mids = bins["bin_mid"].to_numpy(dtype=float)
    matrix = np.full((len(bin_ids), len(dt_values)), np.nan, dtype=float)
    by_bin = {value: i for i, value in enumerate(bin_ids)}
    by_dt = {value: i for i, value in enumerate(dt_values)}

    for _, row in data.iterrows():
        bid = int(row["bin_id"])
        dt = int(row["dt"])
        if bid in by_bin and dt in by_dt:
            matrix[by_bin[bid], by_dt[dt]] = float(row[q_col])

    return matrix, dt_values, mids


def q_summary_values(row: pd.Series, q_col: str) -> Tuple[int, float]:
    mapping = {
        "q_bh_rms_position_structure": (
            "permutation_n_rms_position_structure_q_lt_0_05",
            "permutation_min_rms_position_structure_q",
        ),
        "q_bh_abs_calendar_slope": (
            "permutation_n_calendar_slope_q_lt_0_05",
            "permutation_min_calendar_slope_q",
        ),
        "q_bh_alternation": (
            "permutation_n_alternation_q_lt_0_05",
            "permutation_min_alternation_q",
        ),
    }

    count_col, min_col = mapping[q_col]
    n_sig = pd.to_numeric(pd.Series([row.get(count_col)]), errors="coerce").iloc[0]
    min_q = pd.to_numeric(pd.Series([row.get(min_col)]), errors="coerce").iloc[0]

    return (
        int(n_sig) if np.isfinite(n_sig) else 0,
        float(min_q) if np.isfinite(min_q) else np.nan,
    )


def draw_q_heatmap(
    ax: plt.Axes,
    tests: pd.DataFrame,
    core_bins: pd.DataFrame,
    summary: pd.DataFrame,
    metric: str,
    q_col: str,
    axis_key: str,
    annotate_summary: bool = True,
) -> None:
    matrix, dt_values, mids = prepare_q_matrix(tests, core_bins, metric, q_col)

    image = ax.imshow(
        matrix,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        vmin=0.0,
        vmax=1.0,
    )

    ax.set_xticks(np.arange(len(dt_values)))
    ax.set_xticklabels([str(value) for value in dt_values])

    n = len(mids)
    y_idx = np.unique(np.linspace(0, n - 1, min(6, n)).round().astype(int))
    ax.set_yticks(y_idx)
    ax.set_yticklabels(["{:.2f}".format(mids[i]) for i in y_idx])

    ax.set_xlabel(r"Temporal lag, $\Delta t$")
    ax.set_ylabel(r"Latent abundance, $x^\star$")
    style_axes(ax)
    apply_axis_config(ax, axis_key)

    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("BH-adjusted q")
    cbar.ax.tick_params(
        direction="out",
        width=TICK_LINEWIDTH,
        length=TICK_LENGTH,
    )
    cbar.outline.set_linewidth(AXES_LINEWIDTH)

    if annotate_summary:
        row = summary_row(summary, metric)
        n_sig, min_q = q_summary_values(row, q_col)
        n_tests = int(np.isfinite(matrix).sum())

        subtitle = r"$q<0.05$: {}/{}".format(n_sig, n_tests)
        if np.isfinite(min_q):
            subtitle += r";  $q_{{\min}}={:.3f}$".format(min_q)

        ax.set_title(
            subtitle,
            loc="left",
            pad=7.0,
            fontsize=float(plt.rcParams["font.size"]),
            fontweight="normal",
        )


# =============================================================================
# Anchored signed slopes
# =============================================================================


def get_anchor_rows(
    slopes: pd.DataFrame,
    metric: str,
    subject_matching: str,
    anchor_type: str,
) -> pd.DataFrame:
    require_columns(
        slopes,
        [
            "anchor_type",
            "subject_matching",
            "metric",
            "bin_id",
            "bin_mid",
            "test_status",
            "signed_slope_median",
            "signed_slope_q_low",
            "signed_slope_q_high",
        ],
        "09_anchored_signed_slope_bootstrap.csv",
    )

    data = slopes.loc[
        (slopes["metric"].astype(str) == str(metric))
        & (slopes["subject_matching"].astype(str) == str(subject_matching))
        & (slopes["anchor_type"].astype(str) == str(anchor_type))
        & (slopes["test_status"].astype(str) == "ok")
    ].copy()

    for column in (
        "bin_mid",
        "signed_slope_median",
        "signed_slope_q_low",
        "signed_slope_q_high",
    ):
        data[column] = pd.to_numeric(data[column], errors="coerce")

    data = data.dropna(
        subset=[
            "bin_mid",
            "signed_slope_median",
            "signed_slope_q_low",
            "signed_slope_q_high",
        ]
    ).sort_values("bin_mid")

    if data.empty:
        raise ValueError(
            "No anchored slopes for metric={}, matching={}, anchor={}".format(
                metric, subject_matching, anchor_type
            )
        )
    return data


def anchor_summary_text(
    summary: pd.DataFrame,
    metric: str,
    subject_matching: str,
) -> str:
    row = summary_row(summary, metric)
    prefix = "" if subject_matching == "available" else "matched_"

    start_total = int(row.get(prefix + "common_start_n_bins_valid", 0))
    end_total = int(row.get(prefix + "common_end_n_bins_valid", 0))
    start_neg = int(row.get(prefix + "common_start_n_negative_ci", 0))
    start_pos = int(row.get(prefix + "common_start_n_positive_ci", 0))
    end_neg = int(row.get(prefix + "common_end_n_negative_ci", 0))
    end_pos = int(row.get(prefix + "common_end_n_positive_ci", 0))
    both_neg = int(
        row.get(
            "n_bins_negative_in_both_anchors"
            if subject_matching == "available"
            else "matched_n_bins_negative_in_both_anchors",
            0,
        )
    )

    return (
        "Start CI<0/CI>0: {}/{} / {}/{};  End: {}/{} / {}/{}\n"
        "Both anchors CI<0: {}/{}"
    ).format(
        start_neg,
        start_total,
        start_pos,
        start_total,
        end_neg,
        end_total,
        end_pos,
        end_total,
        both_neg,
        min(start_total, end_total),
    )


def draw_anchor_slopes(
    ax: plt.Axes,
    slopes: pd.DataFrame,
    summary: pd.DataFrame,
    metric: str,
    subject_matching: str,
    axis_key: str,
    show_legend: bool = True,
) -> None:
    series = [
        (
            get_anchor_rows(slopes, metric, subject_matching, "common_start"),
            "Common start",
            "o",
            -0.025,
        ),
        (
            get_anchor_rows(slopes, metric, subject_matching, "common_end"),
            "Common end",
            "s",
            0.025,
        ),
    ]

    for data, label, marker, offset in series:
        x = data["bin_mid"].to_numpy(dtype=float) + offset
        y = data["signed_slope_median"].to_numpy(dtype=float)
        low = data["signed_slope_q_low"].to_numpy(dtype=float)
        high = data["signed_slope_q_high"].to_numpy(dtype=float)
        yerr = np.vstack(
            [
                np.maximum(y - low, 0.0),
                np.maximum(high - y, 0.0),
            ]
        )

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=marker,
            markersize=ERRORBAR_MARKER_SIZE,
            linewidth=DATA_LINEWIDTH,
            elinewidth=ERRORBAR_LINEWIDTH,
            capsize=ERRORBAR_CAPSIZE,
            label=label,
        )

    ax.axhline(0.0, linewidth=REFERENCE_LINEWIDTH, linestyle="--")
    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel(r"Signed lag slope, $D$")
    style_axes(ax)
    apply_axis_config(ax, axis_key)

    ax.set_title(
        anchor_summary_text(summary, metric, subject_matching),
        loc="left",
        pad=7.0,
        fontsize=float(plt.rcParams["font.size"]),
        fontweight="normal",
        linespacing=1.15,
    )

    if show_legend:
        ax.legend(
            frameon=False,
            loc="best",
            borderaxespad=0.25,
            handlelength=1.4,
            labelspacing=0.25,
        )


# =============================================================================
# Composition diagnostics
# =============================================================================


def draw_interval_composition(
    ax: plt.Axes,
    composition: pd.DataFrame,
    metric: str,
    ylabel: str,
    axis_key: str,
) -> None:
    require_columns(
        composition,
        ["dt", "t0", metric],
        "03_interval_composition_by_pair.csv",
    )

    data = composition.copy()
    data["dt"] = pd.to_numeric(data["dt"], errors="coerce")
    data["t0"] = pd.to_numeric(data["t0"], errors="coerce")
    data[metric] = pd.to_numeric(data[metric], errors="coerce")
    data = data.dropna(subset=["dt", "t0", metric])

    for dt in sorted(data["dt"].astype(int).unique().tolist()):
        subset = data.loc[data["dt"].astype(int) == dt].sort_values("t0")
        ax.plot(
            subset["t0"],
            subset[metric],
            marker="o",
            linewidth=DATA_LINEWIDTH,
            markersize=MARKER_SIZE,
            label=r"$\Delta t={}$".format(dt),
        )

    ax.set_xlabel("Calendar interval start, $t_0$")
    ax.set_ylabel(ylabel)
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False, loc="best")


def draw_anchor_complete_subjects(
    ax: plt.Axes,
    matching: pd.DataFrame,
    axis_key: str,
) -> None:
    require_columns(
        matching,
        ["anchor_type", "subject", "complete_all_requested_dt"],
        "11_anchor_subject_matching_diagnostic.csv",
    )

    rows = []
    for anchor_type, group in matching.groupby("anchor_type", sort=True):
        complete = bool_mask(group["complete_all_requested_dt"])
        rows.append(
            {
                "anchor_type": str(anchor_type),
                "n_total": int(len(group)),
                "n_complete": int(np.sum(complete)),
            }
        )

    data = pd.DataFrame(rows)
    x = np.arange(len(data))

    ax.plot(
        x,
        data["n_complete"],
        marker="o",
        linewidth=DATA_LINEWIDTH,
        markersize=MARKER_SIZE,
        label="Complete at all requested lags",
    )
    ax.plot(
        x,
        data["n_total"],
        marker="s",
        linewidth=DATA_LINEWIDTH,
        markersize=MARKER_SIZE,
        linestyle="--",
        label="Structurally available subjects",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(
        [value.replace("_", " ") for value in data["anchor_type"].astype(str)],
        rotation=0,
    )
    ax.set_ylabel("Subjects")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False, loc="best")


def draw_matched_subject_support(
    ax: plt.Axes,
    slopes: pd.DataFrame,
    metric: str,
    axis_key: str,
) -> None:
    require_columns(
        slopes,
        [
            "anchor_type",
            "subject_matching",
            "metric",
            "bin_mid",
            "n_subjects_analysis",
        ],
        "09_anchored_signed_slope_bootstrap.csv",
    )

    data = slopes.loc[
        (slopes["metric"].astype(str) == metric)
        & (slopes["subject_matching"].astype(str) == "matched_all_dt")
    ].copy()
    data["bin_mid"] = pd.to_numeric(data["bin_mid"], errors="coerce")
    data["n_subjects_analysis"] = pd.to_numeric(
        data["n_subjects_analysis"], errors="coerce"
    )

    for anchor_type, label in (
        ("common_start", "Common start"),
        ("common_end", "Common end"),
    ):
        subset = data.loc[data["anchor_type"].astype(str) == anchor_type].dropna(
            subset=["bin_mid", "n_subjects_analysis"]
        ).sort_values("bin_mid")
        ax.plot(
            subset["bin_mid"],
            subset["n_subjects_analysis"],
            marker="o",
            linewidth=DATA_LINEWIDTH,
            markersize=MARKER_SIZE,
            label=label,
        )

    ax.set_xlabel(r"Transition-centred latent abundance, $x^\star$")
    ax.set_ylabel("Matched subjects")
    style_axes(ax)
    apply_axis_config(ax, axis_key)
    ax.legend(frameon=False, loc="best")


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
    stem = "00_main_interval_position_structure_ABCD"
    fig, axes = plt.subplots(2, 2, figsize=FIGURE_SIZES_INCHES[stem])
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    draw_q_heatmap(
        ax_a,
        tables["permutation_tests"],
        tables["core_bins"],
        tables["analysis_summary"],
        PRIMARY_METRIC,
        "q_bh_rms_position_structure",
        "MainA",
        True,
    )
    draw_anchor_slopes(
        ax_b,
        tables["anchored_slopes"],
        tables["analysis_summary"],
        PRIMARY_METRIC,
        "available",
        "MainB",
        False,
    )
    draw_anchor_slopes(
        ax_c,
        tables["anchored_slopes"],
        tables["analysis_summary"],
        PRIMARY_METRIC,
        "matched_all_dt",
        "MainC",
        False,
    )
    draw_anchor_slopes(
        ax_d,
        tables["anchored_slopes"],
        tables["analysis_summary"],
        COMPLEMENTARY_METRIC,
        "matched_all_dt",
        "MainD",
        False,
    )

    for ax, label in zip(axes.ravel(), "ABCD"):
        add_panel_label(ax, label, font_size)

    handles, labels = ax_b.get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.755, 0.985),
            ncol=2,
            frameon=False,
            handlelength=1.6,
            columnspacing=1.8,
        )

    fig.subplots_adjust(
        left=0.095,
        right=0.965,
        bottom=0.095,
        top=0.805,
        wspace=0.50,
        hspace=0.66,
    )
    return save_figure(fig, outdir, stem, formats, dpi)


def composition_available(tables: Dict[str, pd.DataFrame]) -> bool:
    return not tables["composition"].empty and not tables["matching"].empty


def plot_composition(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    dpi: int,
    font_size: float,
) -> List[Path]:
    if not composition_available(tables):
        warnings.warn(
            "Optional Step-14 composition/matching tables are unavailable; "
            "composition composite skipped."
        )
        return []

    stem = "Sx_interval_composition_diagnostics_ABCD"
    fig, axes = plt.subplots(2, 2, figsize=FIGURE_SIZES_INCHES[stem])
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    draw_interval_composition(
        ax_a,
        tables["composition"],
        "n_transitions",
        "Retained transitions",
        "SxA",
    )
    draw_interval_composition(
        ax_b,
        tables["composition"],
        "n_subjects",
        "Subjects",
        "SxB",
    )
    draw_anchor_complete_subjects(ax_c, tables["matching"], "SxC")
    draw_matched_subject_support(
        ax_d,
        tables["anchored_slopes"],
        PRIMARY_METRIC,
        "SxD",
    )

    for ax, label in zip(axes.ravel(), "ABCD"):
        add_panel_label(ax, label, font_size)

    handles, labels = ax_b.get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.755, 0.985),
            ncol=2,
            frameon=False,
            handlelength=1.6,
            columnspacing=1.8,
        )

    fig.subplots_adjust(
        left=0.095,
        right=0.965,
        bottom=0.095,
        top=0.805,
        wspace=0.50,
        hspace=0.66,
    )
    return save_figure(fig, outdir, stem, formats, dpi)


# =============================================================================
# Standalone figures
# =============================================================================


def plot_standalones(
    tables: Dict[str, pd.DataFrame],
    outdir: Path,
    formats: Sequence[str],
    dpi: int,
) -> List[Path]:
    written: List[Path] = []

    for stem, metric in (
        ("01_var_interval_position_q_heatmap", PRIMARY_METRIC),
        ("04_msd_interval_position_q_heatmap", COMPLEMENTARY_METRIC),
    ):
        fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
        draw_q_heatmap(
            ax,
            tables["permutation_tests"],
            tables["core_bins"],
            tables["analysis_summary"],
            metric,
            "q_bh_rms_position_structure",
            stem,
            True,
        )
        written.extend(save_figure(fig, outdir, stem, formats, dpi))

    for stem, metric, matching in (
        ("02_var_anchor_available_slopes", PRIMARY_METRIC, "available"),
        ("03_var_anchor_matched_slopes", PRIMARY_METRIC, "matched_all_dt"),
        ("05_msd_anchor_available_slopes", COMPLEMENTARY_METRIC, "available"),
        ("06_msd_anchor_matched_slopes", COMPLEMENTARY_METRIC, "matched_all_dt"),
    ):
        fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
        draw_anchor_slopes(
            ax,
            tables["anchored_slopes"],
            tables["analysis_summary"],
            metric,
            matching,
            stem,
            True,
        )
        written.extend(save_figure(fig, outdir, stem, formats, dpi))

    stem = "07_var_secondary_interval_q_heatmaps"
    fig, axes = plt.subplots(1, 2, figsize=FIGURE_SIZES_INCHES[stem])
    draw_q_heatmap(
        axes[0],
        tables["permutation_tests"],
        tables["core_bins"],
        tables["analysis_summary"],
        PRIMARY_METRIC,
        "q_bh_abs_calendar_slope",
        "07a_var_calendar_slope_q",
        True,
    )
    axes[0].set_title("Calendar-slope test")
    draw_q_heatmap(
        axes[1],
        tables["permutation_tests"],
        tables["core_bins"],
        tables["analysis_summary"],
        PRIMARY_METRIC,
        "q_bh_alternation",
        "07b_var_alternation_q",
        True,
    )
    axes[1].set_title("Exploratory alternation")
    fig.subplots_adjust(wspace=0.45)
    written.extend(save_figure(fig, outdir, stem, formats, dpi))

    if composition_available(tables):
        for stem, metric, ylabel in (
            ("08_interval_transition_counts", "n_transitions", "Retained transitions"),
            ("09_interval_subject_counts", "n_subjects", "Subjects"),
        ):
            fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
            draw_interval_composition(
                ax,
                tables["composition"],
                metric,
                ylabel,
                stem,
            )
            written.extend(save_figure(fig, outdir, stem, formats, dpi))

        stem = "10_anchor_complete_subjects"
        fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
        draw_anchor_complete_subjects(ax, tables["matching"], stem)
        written.extend(save_figure(fig, outdir, stem, formats, dpi))

        stem = "11_matched_var_subject_support"
        fig, ax = plt.subplots(figsize=FIGURE_SIZES_INCHES[stem])
        draw_matched_subject_support(
            ax,
            tables["anchored_slopes"],
            PRIMARY_METRIC,
            stem,
        )
        written.extend(save_figure(fig, outdir, stem, formats, dpi))

    return written


# =============================================================================
# CLI
# =============================================================================


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Plot ClonoDynamics Step-14 interval-position and anchored "
            "composition-control outputs."
        ),
    )
    parser.add_argument(
        "--analysis-dir",
        required=True,
        type=Path,
        help="Output directory produced by 14-interval_position_structure.py.",
    )
    parser.add_argument(
        "--figures-dir",
        "--outdir",
        dest="figures_dir",
        type=Path,
        default=None,
        help=(
            "Figure output directory. Default: "
            "<analysis-dir>/figures_interval_position_structure"
        ),
    )
    parser.add_argument(
        "--figure-set",
        choices=["all", "main", "composition", "standalone"],
        default="all",
    )
    parser.add_argument("--formats", default="pdf,png")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-family", default=DEFAULT_FONT_FAMILY)
    parser.add_argument("--font-size", type=float, default=DEFAULT_FONT_SIZE)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    analysis_dir = args.analysis_dir.expanduser().resolve()
    if not analysis_dir.is_dir():
        raise NotADirectoryError(
            "Step-14 analysis directory not found: {}".format(analysis_dir)
        )

    figures_dir = (
        args.figures_dir.expanduser().resolve()
        if args.figures_dir is not None
        else analysis_dir / "figures_interval_position_structure"
    )
    figures_dir.mkdir(parents=True, exist_ok=True)

    font_family = setup_style(args.font_family, args.font_size)
    formats = parse_formats(args.formats)
    tables, input_paths = load_tables(analysis_dir)

    written: List[Path] = []
    if args.figure_set in {"all", "main"}:
        written.extend(
            plot_main(
                tables,
                figures_dir,
                formats,
                args.dpi,
                args.font_size,
            )
        )
    if args.figure_set in {"all", "composition"}:
        written.extend(
            plot_composition(
                tables,
                figures_dir,
                formats,
                args.dpi,
                args.font_size,
            )
        )
    if args.figure_set in {"all", "standalone"}:
        written.extend(
            plot_standalones(
                tables,
                figures_dir,
                formats,
                args.dpi,
            )
        )

    pd.DataFrame(
        {
            "file": [path.name for path in written],
            "figure": [path.stem for path in written],
            "format": [path.suffix.lstrip(".") for path in written],
        }
    ).to_csv(figures_dir / "figure_manifest.csv", index=False)

    run_config_path = analysis_dir / "00_run_config.json"
    run_config = {}
    if run_config_path.exists():
        try:
            run_config = json.loads(run_config_path.read_text(encoding="utf-8"))
        except Exception:
            run_config = {}

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "analysis_source": "14-interval_position_structure.py",
        "analysis_dir": str(analysis_dir),
        "figures_dir": str(figures_dir),
        "input_paths": input_paths,
        "figure_set": args.figure_set,
        "formats": formats,
        "dpi": int(args.dpi),
        "font_family_requested": str(args.font_family),
        "font_family_used": font_family,
        "font_size": float(args.font_size),
        "default_panel_size_inches": list(DEFAULT_PANEL_SIZE_INCHES),
        "axes_linewidth": AXES_LINEWIDTH,
        "tick_linewidth": TICK_LINEWIDTH,
        "data_linewidth": DATA_LINEWIDTH,
        "figure_sizes_inches": {
            key: list(value) for key, value in FIGURE_SIZES_INCHES.items()
        },
        "axis_config": AXIS_CONFIG,
        "analysis_free": True,
        "primary_metric": PRIMARY_METRIC,
        "complementary_metric": COMPLEMENTARY_METRIC,
        "primary_interval_position_test": (
            "BH-adjusted q for synchronized RMS calendar-position structure"
        ),
        "anchor_slope_source": "09_anchored_signed_slope_bootstrap.csv",
        "permutations_rerun": False,
        "bh_correction_recomputed": False,
        "anchor_bootstrap_rerun": False,
        "step13_core_redefined": False,
        "analysis_script_version": run_config.get("script_version"),
    }

    (figures_dir / "plot_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("[DONE] Step-14 plotting complete")
    print("[INFO] Analysis directory: {}".format(analysis_dir))
    print("[INFO] Figures directory: {}".format(figures_dir))
    print("[INFO] {} rendered files".format(len(written)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
