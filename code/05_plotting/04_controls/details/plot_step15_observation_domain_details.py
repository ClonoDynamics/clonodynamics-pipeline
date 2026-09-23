#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_step15_observation_domain_details.py
=================================================

Publication plotter for the finalized ClonoDynamics Step 15 v6 fixed-threshold
operational observation-domain sensitivity analysis.

Final contract:
- exact Step-9 equal-weight AB/BA with both folds required;
- Step-11 common4 / xmid_latent signed cross_cov;
- Step-13/14 frozen temporal subject/core contract;
- bootstrap cell support rechecked within every subject resample;
- fixed common operational temporal core required at all lags.

The plotter is analysis-free. It reads finalized Step-15 CSV outputs and never
reclassifies transitions, reruns bootstrap analyses, or changes the temporal
common core.

Recommended main Figure 8
-------------------------
Figure8A_crossing_burden
    Operational observation-boundary crossing fraction among NON_FF transitions
    as a function of xmid_latent, shown for temporal lags 1–5 weeks.

Figure8B_forward_domains
    Replicate-decoupled observed forward drift at dt=1 for PRIMARY, T0PLUS, TT,
    and TF operational domains using the fixed Step-9 abundance grid.

Figure8C_crosscov_domains_dt1
    Cross-replicate fluctuation covariance at dt=1 for PRIMARY, NON_FF, and TT
    using the fixed Step-11 xmid_latent grid.

Figure8D_temporal_slope_forest
    Matched-core temporal slopes for PRIMARY, NON_FF, and TT on the common
    operational temporal core.

Recommended Supplementary Figure 9
----------------------------------
SuppFigS9A_class_composition
    TT/TF/FT/FF transition-class composition across lag.

SuppFigS9B_crossing_imbalance
    Directional TF/FT crossing imbalance across abundance and lag.

SuppFigS9C_forward_differences
    Paired-bootstrap forward differences for T0PLUS-PRIMARY, TT-PRIMARY, and
    T0PLUS-TT.

SuppFigS9D_crosscov_difference_dt1
    Paired-bootstrap cross_cov differences at dt=1 for NON_FF-PRIMARY and
    TT-PRIMARY.

SuppFigS9E_samevar_domains_dt1
    same_var_mean at dt=1 across PRIMARY, NON_FF, and TT.

SuppFigS9F_excess_domains_dt1
    replicate_specific_excess at dt=1 across PRIMARY, NON_FF, and TT.

SuppFigS9G_temporal_profiles
    Matched-core cross_cov profiles across dt=1–5 for PRIMARY, NON_FF, and TT.

Every panel is exported as an independent PDF/PNG. Exact plotted source tables
are also written to <outdir>/figure_source_data.

Typical run
-----------
python3 ./code/05_plotting/04_controls/plot_step15_observation_domain_details.py \
    --step15-dir \
    ./dataset_longitudinal_results/15-detectability_boundary_sensitivity \
    --outdir \
    ./figures/step15_detectability_boundary_sensitivity
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

FONT = {
    "family": "Arial",
    "default_pt": 8.0,
    "axis_label_pt": 8.0,
    "tick_label_pt": 8.0,
    "legend_pt": 7.0,
    "annotation_pt": 7,
    "title_pt": 8.0,
}

EXPORT = {
    "pdf_font_type": 42,
    "ps_font_type": 42,
    "png_dpi": 600,
    "save_formats": ("pdf", "png"),
}

FIGSIZE_INCHES: Dict[str, Tuple[float, float]] = {
    "Figure8A_crossing_burden": (3.40, 2.60),
    "Figure8B_forward_domains": (3.40, 2.60),
    "Figure8C_crosscov_domains_dt1": (3.40, 2.60),
    "Figure8D_temporal_slope_forest": (3.65, 2.60),
    "SuppFigS9A_class_composition": (3.40, 2.60),
    "SuppFigS9B_crossing_imbalance": (3.40, 2.60),
    "SuppFigS9C_forward_differences": (3.40, 2.60),
    "SuppFigS9D_crosscov_difference_dt1": (3.40, 2.60),
    "SuppFigS9E_samevar_domains_dt1": (3.40, 2.60),
    "SuppFigS9F_excess_domains_dt1": (3.40, 2.60),
    "SuppFigS9G_temporal_profiles": (3.40, 2.60),
}

LINE_WIDTHS = {
    "standard": 0.6,
    "primary": 1,
    "errorbar": 0.80,
    "reference": 0.80,
    "axis": 0.80,
    "tick": 0.75,
}

MARKER_SIZES = {
    "standard": 3.6,
    "point": 4.8,
    "errorbar_capsize": 2.2,
}

ALPHAS = {
    "ci_band": 0.18,
    "zero_line": 0.70,
    "reference_line": 0.70,
}

AXIS_STYLE = {
    "tick_direction": "out",
    "major_tick_length": 3.5,
    "grid": False,
}

COLORS = {
    # Temporal lags
    "dt1": "#0072B2",
    "dt2": "#E69F00",
    "dt3": "#009E73",
    "dt4": "#CC79A7",
    "dt5": "#D55E00",

    # Operational domains
    "PRIMARY": "#4D4D4D",
    "T0PLUS": "#E69F00",
    "TT": "#0072B2",
    "TF": "#D55E00",
    "NON_FF": "#009E73",

    # Operational classes
    "FF": "#BDBDBD",
    "FT": "#56B4E9",

    "zero": "#777777",
}

LINESTYLES = {
    "PRIMARY": "-",
    "T0PLUS": "-",
    "TT": "-",
    "TF": ":",
    "NON_FF": "-",
    "zero": "--",
}

MARKERS = {
    "PRIMARY": "o",
    "T0PLUS": "s",
    "TT": "^",
    "TF": "D",
    "NON_FF": "s",
    "dt1": "o",
    "dt2": "s",
    "dt3": "^",
    "dt4": "D",
    "dt5": "v",
}

XLIMS = {
    "Figure8A_crossing_burden": None,
    "Figure8B_forward_domains": None,
    "Figure8C_crosscov_domains_dt1": None,
    "Figure8D_temporal_slope_forest": None,
    "SuppFigS9A_class_composition": (0.8, 5.2),
    "SuppFigS9B_crossing_imbalance": None,
    "SuppFigS9C_forward_differences": None,
    "SuppFigS9D_crosscov_difference_dt1": None,
    "SuppFigS9E_samevar_domains_dt1": None,
    "SuppFigS9F_excess_domains_dt1": None,
    "SuppFigS9G_temporal_profiles": (0.8, 5.2),
}

YLIMS = {
    "Figure8A_crossing_burden": (0.0, 1.02),
    "Figure8B_forward_domains": None,
    "Figure8C_crosscov_domains_dt1": None,
    "Figure8D_temporal_slope_forest": None,
    "SuppFigS9A_class_composition": (0.0, 1.1),
    "SuppFigS9B_crossing_imbalance": (-1.02, 1.02),
    "SuppFigS9C_forward_differences": (-1.45, 0.80),
    "SuppFigS9D_crosscov_difference_dt1": None,
    "SuppFigS9E_samevar_domains_dt1": (0.25, 2.45),
    "SuppFigS9F_excess_domains_dt1": None,
    "SuppFigS9G_temporal_profiles": (0.20, 0.82),
}

LEGEND_SPECS = {
    "Figure8A_crossing_burden": {"loc": "upper right", "ncol": 1},
    "Figure8B_forward_domains": {"loc": "upper right", "ncol": 1},
    "Figure8C_crosscov_domains_dt1": {"loc": "upper right", "ncol": 1},
    "SuppFigS9A_class_composition": {"loc": "upper right", "ncol": 1},
    "SuppFigS9B_crossing_imbalance": {"loc": "lower right", "ncol": 1},
    "SuppFigS9C_forward_differences": {"loc": "lower right", "ncol": 1},
    "SuppFigS9D_crosscov_difference_dt1": {"loc": "best", "ncol": 1},
    "SuppFigS9E_samevar_domains_dt1": {"loc": "upper right", "ncol": 1},
    "SuppFigS9F_excess_domains_dt1": {"loc": "upper right", "ncol": 1},
    "SuppFigS9G_temporal_profiles": {"loc": "upper right", "ncol": 1},
}

ANNOTATIONS = {
    "Figure8A_crossing_burden": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "text": r"$\alpha=0.05$",
    },
    "Figure8B_forward_domains": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "text": r"$\Delta t=1$ week",
    },
    "Figure8C_crosscov_domains_dt1": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "text": r"$\Delta t=1$ week",
    },
    "Figure8D_temporal_slope_forest": {
        "enabled": True,
        "xy_axes": (0.97, 0.04),
        "text": "",
        "ha": "right",
    },
}

SHOW_TITLES = False

PLOT_TITLES = {
    "Figure8A_crossing_burden":
        "Operational boundary crossings",
    "Figure8B_forward_domains":
        "Forward-drift observation-domain sensitivity",
    "Figure8C_crosscov_domains_dt1":
        "Cross-replicate covariance observation-domain sensitivity",
    "Figure8D_temporal_slope_forest":
        "Temporal-scaling observation-domain sensitivity",
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


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(
        direction=AXIS_STYLE["tick_direction"],
    )
    ax.grid(AXIS_STYLE["grid"])


def make_figure(panel: str) -> Tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(
        figsize=FIGSIZE_INCHES[panel],
        constrained_layout=False,
    )
    style_axis(ax)
    return fig, ax


def maybe_title(ax: plt.Axes, panel: str) -> None:
    if SHOW_TITLES:
        ax.set_title(
            PLOT_TITLES.get(panel, panel),
            pad=6,
        )


def apply_limits(ax: plt.Axes, panel: str) -> None:
    xlim = XLIMS.get(panel)
    ylim = YLIMS.get(panel)

    if xlim is not None:
        ax.set_xlim(*xlim)

    if ylim is not None:
        ax.set_ylim(*ylim)


def add_legend(ax: plt.Axes, panel: str) -> None:
    spec = LEGEND_SPECS.get(
        panel,
        {},
    )
    ax.legend(
        frameon=False,
        loc=spec.get(
            "loc",
            "best",
        ),
        ncol=int(
            spec.get(
                "ncol",
                1,
            )
        ),
        borderaxespad=0.25,
        handlelength=1.5,
        labelspacing=0.25,
    )


def zero_line(
    ax: plt.Axes,
    orientation: str = "horizontal",
) -> None:
    kwargs = dict(
        color=COLORS["zero"],
        lw=LINE_WIDTHS["reference"],
        ls=LINESTYLES["zero"],
        alpha=ALPHAS["zero_line"],
        zorder=0,
    )

    if orientation == "horizontal":
        ax.axhline(
            0.0,
            **kwargs,
        )
    else:
        ax.axvline(
            0.0,
            **kwargs,
        )


def add_annotation(
    ax: plt.Axes,
    panel: str,
    text_override: Optional[str] = None,
) -> None:
    spec = ANNOTATIONS.get(
        panel,
        {},
    )

    if not spec.get(
        "enabled",
        False,
    ):
        return

    x, y = spec.get(
        "xy_axes",
        (0.03, 0.04),
    )

    ax.text(
        x,
        y,
        (
            spec.get("text", "")
            if text_override is None
            else str(text_override)
        ),
        transform=ax.transAxes,
        ha=spec.get(
            "ha",
            "left",
        ),
        va=spec.get(
            "va",
            "bottom",
        ),
        fontsize=FONT[
            "annotation_pt"
        ],
    )


def save_figure(
    fig: plt.Figure,
    outdir: Path,
    panel: str,
) -> None:
    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    try:
        fig.tight_layout(
            pad=0.55
        )
    except Exception:
        pass

    for fmt in EXPORT[
        "save_formats"
    ]:
        path = (
            outdir
            / f"{panel}.{fmt}"
        )
        kwargs = {}
        if fmt.lower() == "png":
            kwargs[
                "dpi"
            ] = EXPORT[
                "png_dpi"
            ]

        fig.savefig(
            path,
            **kwargs,
        )

    plt.close(fig)


# =============================================================================
# I/O helpers
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
        obj = json.load(handle)
    if not isinstance(obj, dict):
        raise ValueError(
            f"{label} must contain a JSON object: {path}"
        )
    return obj


def bool_series(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "t", "yes", "y"])
    )


def validate_final_contract(
    step15_dir: Path,
) -> Dict[str, object]:
    config = read_json_required(
        step15_dir / "00_run_config.json",
        "Step 15 run config",
    )
    signature = read_json_required(
        step15_dir / "00_run_signature.json",
        "Step 15 run signature",
    )

    if config.get("primary_estimands_unchanged") is not True:
        raise ValueError(
            "Final plotter requires primary_estimands_unchanged=true."
        )
    if float(config.get("reference_operational_alpha", np.nan)) != 0.05:
        raise ValueError(
            "This Step-15 plotter expects the reference alpha=0.05 run."
        )

    forward = config.get("forward", {})
    if forward.get("single_fold_fallback") is not False:
        raise ValueError(
            "Final plotter requires Step-15 v6 with no single-fold fallback."
        )
    if forward.get("fold_combination") != (
        "exact 0.5 AB + 0.5 BA after binning; both folds required"
    ):
        raise ValueError(
            "Final plotter requires exact equal-weight AB/BA fold contract."
        )

    boot = config.get("bootstrap", {})
    if boot.get("cell_support_rechecked_within_each_draw") is not True:
        raise ValueError(
            "Final plotter requires per-bootstrap support rechecking."
        )
    if boot.get("temporal_slope_requires_complete_fixed_core_all_dt") is not True:
        raise ValueError(
            "Final plotter requires fixed-core temporal bootstrap support."
        )

    if config.get("analysis_signature") != signature.get("analysis_signature"):
        raise ValueError(
            "Step-15 run config/signature analysis signatures disagree."
        )

    payload = signature.get("signature_payload", {})
    if payload.get("forward_fold_combination") != (
        "exact equal AB/BA weight; both folds required"
    ):
        raise ValueError(
            "Step-15 signature does not certify final forward fold contract."
        )
    if payload.get("bootstrap_cell_support_rechecked_each_draw") is not True:
        raise ValueError(
            "Step-15 signature does not certify bootstrap support contract."
        )
    if payload.get("temporal_bootstrap_requires_complete_fixed_core_all_dt") is not True:
        raise ValueError(
            "Step-15 signature does not certify fixed temporal-core bootstrap."
        )

    temporal = config.get("temporal_sensitivity", {})
    common_core = [
        int(x)
        for x in temporal.get("common_operational_core_bins", [])
    ]
    subjects = [
        int(x)
        for x in temporal.get("subject_universe", [])
    ]

    if not common_core or not subjects:
        raise ValueError(
            "Step-15 final run has empty temporal core or subject universe."
        )

    return {
        "step15_script_version": config.get("script_version"),
        "step15_analysis_signature": config.get("analysis_signature"),
        "reference_alpha": float(config.get("reference_operational_alpha")),
        "n_common_operational_core_bins": int(len(common_core)),
        "common_operational_core_bins": ";".join(map(str, common_core)),
        "n_temporal_subjects": int(len(subjects)),
        "temporal_subjects": ";".join(map(str, subjects)),
        "n_bootstrap": int(boot.get("n_bootstrap", 0)),
        "step9_analysis_signature": payload.get("step9_analysis_signature"),
        "step11_analysis_signature": payload.get("step11_analysis_signature"),
        "step13_analysis_signature": payload.get("step13_analysis_signature"),
        "step14_analysis_signature": payload.get("step14_analysis_signature"),
    }


def read_required(
    step15_dir: Path,
    filename: str,
) -> pd.DataFrame:
    path = (
        step15_dir
        / filename
    )

    if not path.is_file():
        raise FileNotFoundError(
            path
        )

    frame = pd.read_csv(
        path
    )

    if frame.empty:
        raise ValueError(
            f"Empty Step-15 table: {path}"
        )

    return frame


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
    step15_dir: Path,
) -> Dict[str, object]:
    contract = validate_final_contract(
        step15_dir
    )
    return {
        "contract": contract,
        "class_composition":
            read_required(
                step15_dir,
                "01_operational_class_composition_by_dt.csv",
            ),
        "crossing":
            read_required(
                step15_dir,
                "02_boundary_crossing_by_dt_bin.csv",
            ),
        "forward":
            read_required(
                step15_dir,
                "03_forward_domain_by_bin.csv",
            ),
        "forward_summary":
            read_required(
                step15_dir,
                "04_forward_domain_summary.csv",
            ),
        "forward_difference":
            read_required(
                step15_dir,
                "05_forward_difference_vs_primary.csv",
            ),
        "forward_audit":
            read_required(
                step15_dir,
                "06_forward_primary_audit.csv",
            ),
        "fluctuation":
            read_required(
                step15_dir,
                "07_fluctuation_domain_by_bin_dt.csv",
            ),
        "fluctuation_difference":
            read_required(
                step15_dir,
                "08_fluctuation_difference_vs_primary.csv",
            ),
        "fluctuation_audit":
            read_required(
                step15_dir,
                "10_fluctuation_primary_audit.csv",
            ),
        "temporal_core":
            read_required(
                step15_dir,
                "11_temporal_common_core.csv",
            ),
        "temporal_profiles":
            read_required(
                step15_dir,
                "12_temporal_domain_by_dt.csv",
            ),
        "temporal_slopes":
            read_required(
                step15_dir,
                "13_temporal_domain_slope_summary.csv",
            ),
        "temporal_support_audit":
            read_required(
                step15_dir,
                "14_temporal_support_audit.csv",
            ),
    }


# =============================================================================
# Main Figure 8A
# =============================================================================

def plot_figure8A(
    crossing: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure8A_crossing_burden"

    d = crossing[
        [
            "dt",
            "bin",
            "x_center",
            "n_NON_FF",
            "n_crossing",
            "crossing_fraction_non_ff",
        ]
    ].copy().sort_values(
        [
            "dt",
            "x_center",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )

    for dt in sorted(
        d[
            "dt"
        ].astype(int).unique()
    ):
        g = d[
            d[
                "dt"
            ].astype(int).eq(
                int(dt)
            )
        ].sort_values(
            "x_center"
        )

        key = f"dt{int(dt)}"

        ax.plot(
            g[
                "x_center"
            ],
            g[
                "crossing_fraction_non_ff"
            ],
            color=COLORS[
                key
            ],
            marker=MARKERS[
                key
            ],
            ms=MARKER_SIZES[
                "standard"
            ],
            lw=LINE_WIDTHS[
                "standard"
            ],
            label=(
                f"{dt} week"
                if int(dt) == 1
                else f"{dt} weeks"
            ),
        )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        "Boundary-crossing fraction among NON_FF"
    )

    add_legend(
        ax,
        panel,
    )
    add_annotation(
        ax,
        panel,
    )
    maybe_title(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Main Figure 8B
# =============================================================================

def plot_figure8B(
    forward: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure8B_forward_domains"

    d = forward.copy().sort_values(
        [
            "mode",
            "x_center",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )
    zero_line(
        ax
    )

    mode_order = [
        "PRIMARY",
        "T0PLUS",
        "TT",
        "TF",
    ]

    for mode in mode_order:
        g = d[
            d[
                "mode"
            ].astype(str).eq(
                mode
            )
        ].sort_values(
            "x_center"
        )

        if g.empty:
            continue

        lw = (
            LINE_WIDTHS[
                "primary"
            ]
            if mode == "PRIMARY"
            else LINE_WIDTHS[
                "standard"
            ]
        )

        ax.plot(
            g[
                "x_center"
            ],
            g[
                "mean_dx"
            ],
            color=COLORS[
                mode
            ],
            ls=LINESTYLES[
                mode
            ],
            marker=MARKERS[
                mode
            ],
            ms=MARKER_SIZES[
                "standard"
            ],
            lw=lw,
            label=mode,
        )

        # PRIMARY band only, to avoid four overlapping uncertainty ribbons.
        if mode == "PRIMARY":
            ax.fill_between(
                g[
                    "x_center"
                ].to_numpy(float),
                g[
                    "bootstrap_q025"
                ].to_numpy(float),
                g[
                    "bootstrap_q975"
                ].to_numpy(float),
                color=COLORS[
                    mode
                ],
                alpha=ALPHAS[
                    "ci_band"
                ],
                linewidth=0,
            )

    ax.set_xlabel(
        r"Initial observed log-frequency, $x_0$"
    )
    ax.set_ylabel(
        r"Mean replicate-decoupled displacement\n$E[\Delta x]$"
    )

    add_legend(
        ax,
        panel,
    )
    add_annotation(
        ax,
        panel,
    )
    maybe_title(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Domain metric curve helper
# =============================================================================

def plot_domain_metric_curve(
    fluctuation: pd.DataFrame,
    metric: str,
    panel: str,
    ylabel: str,
    outdir: Path,
    source_dir: Path,
    dt: int = 1,
) -> None:
    d = fluctuation[
        fluctuation[
            "dt"
        ].astype(int).eq(
            int(dt)
        )
    ].copy().sort_values(
        [
            "mode",
            "x_center",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )

    if metric == "cross_cov":
        zero_line(
            ax
        )

    for mode in [
        "PRIMARY",
        "NON_FF",
        "TT",
    ]:
        g = d[
            d[
                "mode"
            ].astype(str).eq(
                mode
            )
        ].sort_values(
            "x_center"
        )

        if g.empty:
            continue

        lw = (
            LINE_WIDTHS[
                "primary"
            ]
            if mode == "PRIMARY"
            else LINE_WIDTHS[
                "standard"
            ]
        )

        ax.plot(
            g[
                "x_center"
            ],
            g[
                metric
            ],
            color=COLORS[
                mode
            ],
            ls=LINESTYLES[
                mode
            ],
            marker=MARKERS[
                mode
            ],
            ms=MARKER_SIZES[
                "standard"
            ],
            lw=lw,
            label=mode,
        )

        if mode == "PRIMARY":
            low_col = (
                f"{metric}_ci025"
            )
            high_col = (
                f"{metric}_ci975"
            )

            if (
                low_col in g.columns
                and high_col in g.columns
            ):
                ax.fill_between(
                    g[
                        "x_center"
                    ].to_numpy(float),
                    g[
                        low_col
                    ].to_numpy(float),
                    g[
                        high_col
                    ].to_numpy(float),
                    color=COLORS[
                        mode
                    ],
                    alpha=ALPHAS[
                        "ci_band"
                    ],
                    linewidth=0,
                )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        ylabel
    )

    add_legend(
        ax,
        panel,
    )
    add_annotation(
        ax,
        panel,
    )
    maybe_title(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Main Figure 8D
# =============================================================================

def plot_figure8D(
    temporal_slopes: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure8D_temporal_slope_forest"

    mode_order = [
        "PRIMARY",
        "NON_FF",
        "TT",
    ]

    d = temporal_slopes[
        temporal_slopes[
            "mode"
        ].astype(str).isin(
            mode_order
        )
    ].copy()

    order = {
        mode: i
        for i, mode in enumerate(
            mode_order
        )
    }

    d[
        "_order"
    ] = d[
        "mode"
    ].map(
        order
    )
    d = d.sort_values(
        "_order"
    ).drop(
        columns=[
            "_order",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )
    zero_line(
        ax,
        orientation="vertical",
    )

    y = np.arange(
        len(
            mode_order
        )
    )[::-1]

    for yy, mode in zip(
        y,
        mode_order,
    ):
        row = d[
            d[
                "mode"
            ].astype(str).eq(
                mode
            )
        ].iloc[0]

        estimate = float(
            row[
                "observed_slope_per_week"
            ]
        )
        lo = float(
            row[
                "bootstrap_q025"
            ]
        )
        hi = float(
            row[
                "bootstrap_q975"
            ]
        )

        ax.plot(
            [
                lo,
                hi,
            ],
            [
                yy,
                yy,
            ],
            color=COLORS[
                mode
            ],
            lw=LINE_WIDTHS[
                "errorbar"
            ],
        )

        ax.plot(
            estimate,
            yy,
            marker=MARKERS[
                mode
            ],
            color=COLORS[
                mode
            ],
            ms=MARKER_SIZES[
                "point"
            ],
            linestyle="None",
        )

    ax.set_yticks(
        y
    )
    ax.set_yticklabels(
        mode_order
    )
    ax.set_xlabel(
        r"Temporal slope of cross-replicate covariance (week$^{-1}$)"
    )

    core_bins = int(
        pd.to_numeric(
            d["n_common_core_bins"],
            errors="coerce",
        ).dropna().iloc[0]
    )
    n_subjects = int(
        pd.to_numeric(
            d["n_subjects_universe"],
            errors="coerce",
        ).dropna().iloc[0]
    )

    add_annotation(
        ax,
        panel,
        text_override=(
            f"{core_bins}-bin common operational core\n"
            f"{n_subjects} subjects"
        ),
    )
    maybe_title(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S9A
# =============================================================================

def plot_supp_s9A(
    class_composition: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS9A_class_composition"

    d = class_composition.copy().sort_values(
        [
            "obs_class_step15",
            "dt",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )

    class_colors = {
        "TT": COLORS[
            "TT"
        ],
        "TF": COLORS[
            "TF"
        ],
        "FT": COLORS[
            "FT"
        ],
        "FF": COLORS[
            "FF"
        ],
    }

    class_markers = {
        "TT": "^",
        "TF": "D",
        "FT": "s",
        "FF": "o",
    }

    for cls in [
        "TT",
        "TF",
        "FT",
        "FF",
    ]:
        g = d[
            d[
                "obs_class_step15"
            ].astype(str).eq(
                cls
            )
        ].sort_values(
            "dt"
        )

        ax.plot(
            g[
                "dt"
            ],
            g[
                "fraction_within_dt"
            ],
            color=class_colors[
                cls
            ],
            marker=class_markers[
                cls
            ],
            ms=MARKER_SIZES[
                "standard"
            ],
            lw=LINE_WIDTHS[
                "standard"
            ],
            label=cls,
        )

    ax.set_xlabel(
        "Temporal lag (weeks)"
    )
    ax.set_ylabel(
        "Fraction of transitions"
    )
    ax.set_xticks(
        [
            1,
            2,
            3,
            4,
            5,
        ]
    )

    add_legend(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S9B
# =============================================================================

def plot_supp_s9B(
    crossing: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS9B_crossing_imbalance"

    d = crossing[
        [
            "dt",
            "bin",
            "x_center",
            "n_crossing",
            "crossing_imbalance",
        ]
    ].copy().sort_values(
        [
            "dt",
            "x_center",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )
    zero_line(
        ax
    )

    for dt in sorted(
        d[
            "dt"
        ].astype(int).unique()
    ):
        g = d[
            d[
                "dt"
            ].astype(int).eq(
                int(dt)
            )
        ].sort_values(
            "x_center"
        )

        key = f"dt{int(dt)}"

        ax.plot(
            g[
                "x_center"
            ],
            g[
                "crossing_imbalance"
            ],
            color=COLORS[
                key
            ],
            marker=MARKERS[
                key
            ],
            ms=MARKER_SIZES[
                "standard"
            ],
            lw=LINE_WIDTHS[
                "standard"
            ],
            label=(
                f"{dt} week"
                if int(dt) == 1
                else f"{dt} weeks"
            ),
        )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        r"Crossing imbalance, $(FT-TF)/(FT+TF)$"
    )

    add_legend(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S9C
# =============================================================================

def plot_supp_s9C(
    forward_difference: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS9C_forward_differences"

    d = forward_difference.copy().sort_values(
        [
            "contrast",
            "x_center",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )
    zero_line(
        ax
    )

    specs = [
        (
            "T0PLUS-PRIMARY",
            COLORS[
                "T0PLUS"
            ],
            "s",
        ),
        (
            "TT-PRIMARY",
            COLORS[
                "TT"
            ],
            "^",
        ),
        (
            "T0PLUS-TT",
            COLORS[
                "TF"
            ],
            "D",
        ),
    ]

    for contrast, color, marker in specs:
        g = d[
            d[
                "contrast"
            ].astype(str).eq(
                contrast
            )
        ].sort_values(
            "x_center"
        )

        ax.plot(
            g[
                "x_center"
            ],
            g[
                "difference"
            ],
            color=color,
            marker=marker,
            ms=MARKER_SIZES[
                "standard"
            ],
            lw=LINE_WIDTHS[
                "standard"
            ],
            label=contrast,
        )

        ax.fill_between(
            g["x_center"].to_numpy(float),
            g["bootstrap_difference_q025"].to_numpy(float),
            g["bootstrap_difference_q975"].to_numpy(float),
            color=color,
            alpha=0.10,
            linewidth=0,
        )

    ax.set_xlabel(
        r"Initial observed log-frequency, $x_0$"
    )
    ax.set_ylabel(
        r"Difference in mean displacement"
    )

    add_legend(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S9D
# =============================================================================

def plot_supp_s9D(
    fluctuation_difference: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS9D_crosscov_difference_dt1"

    d = fluctuation_difference[
        fluctuation_difference[
            "metric"
        ].astype(str).eq(
            "cross_cov"
        )
        & fluctuation_difference[
            "dt"
        ].astype(int).eq(
            1
        )
    ].copy().sort_values(
        [
            "mode",
            "x_center",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )
    zero_line(
        ax
    )

    for mode in [
        "NON_FF",
        "TT",
    ]:
        g = d[
            d[
                "mode"
            ].astype(str).eq(
                mode
            )
        ].sort_values(
            "x_center"
        )

        ax.plot(
            g[
                "x_center"
            ],
            g[
                "difference"
            ],
            color=COLORS[
                mode
            ],
            marker=MARKERS[
                mode
            ],
            ms=MARKER_SIZES[
                "standard"
            ],
            lw=LINE_WIDTHS[
                "standard"
            ],
            label=(
                f"{mode} - PRIMARY"
            ),
        )

        ax.fill_between(
            g[
                "x_center"
            ].to_numpy(float),
            g[
                "bootstrap_difference_q025"
            ].to_numpy(float),
            g[
                "bootstrap_difference_q975"
            ].to_numpy(float),
            color=COLORS[
                mode
            ],
            alpha=ALPHAS[
                "ci_band"
            ],
            linewidth=0,
        )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        r"Difference in cross-replicate covariance"
    )

    add_legend(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


# =============================================================================
# Supplementary Figure S9G
# =============================================================================

def plot_supp_s9G(
    temporal_profiles: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS9G_temporal_profiles"

    d = temporal_profiles.copy().sort_values(
        [
            "mode",
            "dt",
        ]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(
        panel
    )

    for mode in [
        "PRIMARY",
        "NON_FF",
        "TT",
    ]:
        g = d[
            d[
                "mode"
            ].astype(str).eq(
                mode
            )
        ].sort_values(
            "dt"
        )

        ax.plot(
            g[
                "dt"
            ],
            g[
                "value"
            ],
            color=COLORS[
                mode
            ],
            marker=MARKERS[
                mode
            ],
            ms=MARKER_SIZES[
                "standard"
            ],
            lw=(
                LINE_WIDTHS[
                    "primary"
                ]
                if mode == "PRIMARY"
                else LINE_WIDTHS[
                    "standard"
                ]
            ),
            label=mode,
        )

    ax.set_xlabel(
        "Temporal lag (weeks)"
    )
    ax.set_ylabel(
        "Common-core cross-replicate covariance"
    )
    ax.set_xticks(
        [
            1,
            2,
            3,
            4,
            5,
        ]
    )

    add_legend(
        ax,
        panel,
    )
    apply_limits(
        ax,
        panel,
    )
    save_figure(
        fig,
        outdir,
        panel,
    )


def write_plot_contract_audit(
    outdir: Path,
    contract: Mapping[str, object],
    forward_audit: pd.DataFrame,
    fluctuation_audit: pd.DataFrame,
    temporal_core: pd.DataFrame,
    temporal_slopes: pd.DataFrame,
    temporal_support_audit: pd.DataFrame,
) -> None:
    fwd = pd.to_numeric(
        forward_audit["difference"],
        errors="coerce",
    ).to_numpy(float)
    fwd = fwd[np.isfinite(fwd)]

    fluct_cols = [
        c for c in fluctuation_audit.columns
        if c.endswith("_difference")
    ]
    fluct_vals = []
    for col in fluct_cols:
        x = pd.to_numeric(
            fluctuation_audit[col],
            errors="coerce",
        ).to_numpy(float)
        fluct_vals.extend(
            x[np.isfinite(x)].tolist()
        )

    selected = temporal_core[
        bool_series(
            temporal_core[
                "selected_common_operational_temporal_core"
            ]
        )
    ].copy()

    row = {
        **contract,
        "forward_primary_audit_max_abs_diff": (
            float(np.max(np.abs(fwd)))
            if fwd.size
            else np.nan
        ),
        "fluctuation_primary_audit_max_abs_diff": (
            float(np.max(np.abs(fluct_vals)))
            if fluct_vals
            else np.nan
        ),
        "temporal_core_rows_marked_selected": int(len(selected)),
        "temporal_support_rows": int(len(temporal_support_audit)),
    }

    for mode in ["PRIMARY", "NON_FF", "TT"]:
        g = temporal_slopes[
            temporal_slopes["mode"].astype(str).eq(mode)
        ]
        if len(g) == 1:
            row[f"{mode}_temporal_slope"] = float(
                g.iloc[0]["observed_slope_per_week"]
            )
            row[f"{mode}_bootstrap_q025"] = float(
                g.iloc[0]["bootstrap_q025"]
            )
            row[f"{mode}_bootstrap_q975"] = float(
                g.iloc[0]["bootstrap_q975"]
            )
            row[f"{mode}_n_bootstrap_valid"] = int(
                g.iloc[0]["n_bootstrap_valid"]
            )

    if int(row["temporal_core_rows_marked_selected"]) != int(
        contract["n_common_operational_core_bins"]
    ):
        raise ValueError(
            "Plotter audit: temporal core count differs from Step-15 config."
        )

    if (
        np.isfinite(row["forward_primary_audit_max_abs_diff"])
        and row["forward_primary_audit_max_abs_diff"] > 1e-10
    ):
        raise ValueError(
            "Plotter audit: forward PRIMARY no longer matches Step 9."
        )

    if (
        np.isfinite(row["fluctuation_primary_audit_max_abs_diff"])
        and row["fluctuation_primary_audit_max_abs_diff"] > 1e-10
    ):
        raise ValueError(
            "Plotter audit: fluctuation PRIMARY no longer matches Step 11."
        )

    pd.DataFrame([row]).to_csv(
        outdir / "00_plot_contract_audit.csv",
        index=False,
    )


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Plot finalized Step-15 operational observation-domain sensitivity."
        ),
    )

    p.add_argument(
        "--step15-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--figure-set",
        choices=[
            "all",
            "main",
            "supplementary",
        ],
        default="all",
    )

    return p


def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(
        argv
    )

    configure_matplotlib()

    step15_dir = (
        args.step15_dir
        .expanduser()
        .resolve(strict=True)
    )

    outdir = (
        args.outdir
        .expanduser()
        .resolve()
    )

    main_dir = (
        outdir
        / "main_figure"
    )
    supp_dir = (
        outdir
        / "supplementary_figure"
    )
    source_dir = (
        outdir
        / "figure_source_data"
    )

    data = load_inputs(
        step15_dir
    )

    if args.figure_set in {
        "all",
        "main",
    }:
        plot_figure8A(
            data[
                "crossing"
            ],
            main_dir,
            source_dir,
        )

        plot_figure8B(
            data[
                "forward"
            ],
            main_dir,
            source_dir,
        )

        plot_domain_metric_curve(
            data[
                "fluctuation"
            ],
            "cross_cov",
            "Figure8C_crosscov_domains_dt1",
            "Cross-replicate covariance",
            main_dir,
            source_dir,
            dt=1,
        )

        plot_figure8D(
            data[
                "temporal_slopes"
            ],
            main_dir,
            source_dir,
        )

    if args.figure_set in {
        "all",
        "supplementary",
    }:
        plot_supp_s9A(
            data[
                "class_composition"
            ],
            supp_dir,
            source_dir,
        )

        plot_supp_s9B(
            data[
                "crossing"
            ],
            supp_dir,
            source_dir,
        )

        plot_supp_s9C(
            data[
                "forward_difference"
            ],
            supp_dir,
            source_dir,
        )

        plot_supp_s9D(
            data[
                "fluctuation_difference"
            ],
            supp_dir,
            source_dir,
        )

        plot_domain_metric_curve(
            data[
                "fluctuation"
            ],
            "same_var_mean",
            "SuppFigS9E_samevar_domains_dt1",
            "Mean within-replicate variance",
            supp_dir,
            source_dir,
            dt=1,
        )

        plot_domain_metric_curve(
            data[
                "fluctuation"
            ],
            "replicate_specific_excess",
            "SuppFigS9F_excess_domains_dt1",
            "Replicate-specific excess",
            supp_dir,
            source_dir,
            dt=1,
        )

        plot_supp_s9G(
            data[
                "temporal_profiles"
            ],
            supp_dir,
            source_dir,
        )

    write_plot_contract_audit(
        outdir,
        data["contract"],
        data["forward_audit"],
        data["fluctuation_audit"],
        data["temporal_core"],
        data["temporal_slopes"],
        data["temporal_support_audit"],
    )

    print(
        "\n[DONE] Step 15 publication plotting complete"
    )
    print(
        "Main Figure 8:",
        main_dir,
    )
    print(
        "Supplementary Figure 9:",
        supp_dir,
    )
    print(
        "Exact plotted source data:",
        source_dir,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
