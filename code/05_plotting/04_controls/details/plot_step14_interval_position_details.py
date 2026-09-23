#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_step14_interval_position_details.py
==========================================

Publication plotter for the finalized ClonoDynamics Step 14
interval-position / anchored-composition control.

Final contract:
- Step 14 v2 inherits Step-13 fixed core and complete-case subjects;
- primary interval-bin support is inherited from Step 13;
- Step 11 / Step 13 signatures are propagated;
- calendar-position permutation uses the same within-subject exchangeability
  design, vectorized only for runtime.

The plotter is analysis-free: it reads finalized Step-14 outputs and does not
recompute scientific quantities.

Recommended main Figure 7
-------------------------
Figure7A_calendar_position_profiles
    Fixed-core cross-replicate covariance across calendar interval start
    positions for the lag values with >=3 valid interval positions.

Figure7B_anchor_slope_forest
    Step-13 pooled temporal slope together with common-start/common-end anchored
    controls, including available-subject and matched-all-dt analyses.

Figure7C_common_start_subject_slopes
    Subject-specific temporal slopes in the strict matched common-start control.

Figure7D_common_end_subject_slopes
    Subject-specific temporal slopes in the strict matched common-end control.

Recommended Supplementary Figure 8
----------------------------------
SuppFigS8A_core_position_qvalues
    Core-level calendar-position BH-adjusted q values across fluctuation
    metrics and temporal lags.

SuppFigS8B_binwise_position_qvalues
    Abundance-resolved cross_cov calendar-position BH-adjusted q values.

SuppFigS8C_common_start_binwise_slopes
    Abundance-resolved matched common-start temporal slopes.

SuppFigS8D_common_end_binwise_slopes
    Abundance-resolved matched common-end temporal slopes.

Typical run
-----------
python3 plot_step14_interval_position_details.py \
    --step14-dir ./dataset_longitudinal_results/14-interval_position_structure \
    --outdir ./figures/step14_interval_position_structure

All visual settings are centralized in the dictionaries below.
Default typography is Arial 8 pt.
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
    "legend_pt": 8.0,
    "annotation_pt": 8.0,
    "title_pt": 8.0,
}

EXPORT = {
    "pdf_font_type": 42,
    "ps_font_type": 42,
    "png_dpi": 600,
    "save_formats": ("pdf", "png"),
    # Keep fixed physical canvas sizes. Do not use bbox_inches="tight".
    "pad_inches": 0.05,
}

FIGSIZE_INCHES: Dict[str, Tuple[float, float]] = {
    "Figure7A_calendar_position_profiles": (3.40, 2.60),
    "Figure7B_anchor_slope_forest": (3.70, 2.60),
    "Figure7C_common_start_subject_slopes": (3.40, 2.60),
    "Figure7D_common_end_subject_slopes": (3.40, 2.60),
    "SuppFigS8A_core_position_qvalues": (3.40, 2.60),
    "SuppFigS8B_binwise_position_qvalues": (3.40, 2.60),
    "SuppFigS8C_common_start_binwise_slopes": (3.40, 2.60),
    "SuppFigS8D_common_end_binwise_slopes": (3.40, 2.60),
}

LINE_WIDTHS = {
    "standard": 1,
    "primary": 1.20,
    "reference": 0.8,
    "errorbar": 0.8,
    "zero": 0.80,
    "axis": 0.80,
    "tick": 0.75,
}

MARKER_SIZES = {
    "standard": 3.8,
    "point": 4.8,
    "subject": 4.5,
    "scatter_area": 24.0,
    "errorbar_capsize": 2.5,
}

ALPHAS = {
    "cohort_ci": 0.15,
    "subject_point": 0.95,
    "zero_line": 0.65,
    "q_threshold": 0.75,
}

AXIS_STYLE = {
    "tick_direction": "out",
    "major_tick_length": 3.5,
    "grid": False,
}

COLORS = {
    "dt1": "#0072B2",
    "dt2": "#E69F00",
    "dt3": "#009E73",
    "reference": "#4D4D4D",
    "common_start": "#0072B2",
    "common_end": "#D55E00",
    "same_var_mean": "#E69F00",
    "excess": "#CC79A7",
    "zero": "#777777",
    "subject": "#4D4D4D",
}

LINESTYLES = {
    "dt1": "-",
    "dt2": "-",
    "dt3": "-",
    "available": "-",
    "matched": "-",
    "zero": "--",
    "q_threshold": "--",
}

MARKERS = {
    "dt1": "o",
    "dt2": "s",
    "dt3": "^",
    "available": "o",
    "matched": "s",
    "subject": "o",
}

XLIMS = {
    "Figure7A_calendar_position_profiles": (0.8, 5.2),
    "Figure7B_anchor_slope_forest": None,
    "Figure7C_common_start_subject_slopes": None,
    "Figure7D_common_end_subject_slopes": None,
    "SuppFigS8A_core_position_qvalues": (0.7, 3.3),
    "SuppFigS8B_binwise_position_qvalues": None,
    "SuppFigS8C_common_start_binwise_slopes": None,
    "SuppFigS8D_common_end_binwise_slopes": None,
}

YLIMS = {
    "Figure7A_calendar_position_profiles": None,
    "Figure7B_anchor_slope_forest": None,
    "Figure7C_common_start_subject_slopes": None,
    "Figure7D_common_end_subject_slopes": None,
    "SuppFigS8A_core_position_qvalues": (0.0, 1.02),
    "SuppFigS8B_binwise_position_qvalues": (0.0, 1.5),
    "SuppFigS8C_common_start_binwise_slopes": None,
    "SuppFigS8D_common_end_binwise_slopes": None,
}

LEGEND_SPECS = {
    "Figure7A_calendar_position_profiles": {
        "loc": "upper left",
        "ncol": 1,
        "fontsize": 8.0,
    },
    "SuppFigS8A_core_position_qvalues": {
        "loc": "upper left",
        "ncol": 1,
        "fontsize": 8.0,
    },
    "SuppFigS8B_binwise_position_qvalues": {
        "loc": "upper left",
        "ncol": 1,
        "fontsize": 8.0,
    },
}

ANNOTATIONS = {
    "Figure7A_calendar_position_profiles": {
        "enabled": True,
        "xy_axes": (0.03, 0.035),
        "ha": "left",
        "va": "bottom",
        "fontsize": 7.5,
    },
    "Figure7C_common_start_subject_slopes": {
        "enabled": True,
        "xy_axes": (0.5, 0.02),
        "ha": "right",
        "va": "bottom",
        "fontsize": 7.5,
    },
    "Figure7D_common_end_subject_slopes": {
        "enabled": True,
        "xy_axes": (0.97, 0.035),
        "ha": "right",
        "va": "bottom",
        "fontsize": 7.5,
    },
}

SHOW_TITLES = False

PLOT_TITLES = {
    "Figure7A_calendar_position_profiles":
        "Calendar-position profiles",
    "Figure7B_anchor_slope_forest":
        "Anchored temporal-slope controls",
    "Figure7C_common_start_subject_slopes":
        "Matched common-start subject slopes",
    "Figure7D_common_end_subject_slopes":
        "Matched common-end subject slopes",
    "SuppFigS8A_core_position_qvalues":
        "Core-level calendar-position tests",
    "SuppFigS8B_binwise_position_qvalues":
        "Abundance-resolved calendar-position tests",
    "SuppFigS8C_common_start_binwise_slopes":
        "Matched common-start abundance-resolved slopes",
    "SuppFigS8D_common_end_binwise_slopes":
        "Matched common-end abundance-resolved slopes",
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
        cur = ax.get_xlim()
        lo, hi = xlim
        ax.set_xlim(
            cur[0] if lo is None else lo,
            cur[1] if hi is None else hi,
        )

    if ylim is not None:
        cur = ax.get_ylim()
        lo, hi = ylim
        ax.set_ylim(
            cur[0] if lo is None else lo,
            cur[1] if hi is None else hi,
        )


def zero_line(
    ax: plt.Axes,
    orientation: str = "horizontal",
) -> None:
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


def add_annotation(
    ax: plt.Axes,
    panel: str,
    text: str,
) -> None:
    spec = ANNOTATIONS.get(panel, {})
    if not spec.get("enabled", False):
        return

    x, y = spec.get(
        "xy_axes",
        (0.03, 0.04),
    )

    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha=spec.get("ha", "left"),
        va=spec.get("va", "bottom"),
        fontsize=spec.get(
            "fontsize",
            FONT["annotation_pt"],
        ),
    )


def add_panel_legend(
    ax: plt.Axes,
    panel: str,
) -> None:
    spec = LEGEND_SPECS.get(panel, {})
    ax.legend(
        frameon=False,
        loc=spec.get("loc", "best"),
        ncol=int(spec.get("ncol", 1)),
        fontsize=spec.get(
            "fontsize",
            FONT["legend_pt"],
        ),
        borderaxespad=0.25,
        handlelength=1.5,
        labelspacing=0.25,
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

    # Preserve the exact physical canvas dimensions in FIGSIZE_INCHES.
    try:
        fig.tight_layout(
            pad=0.55
        )
    except Exception:
        pass

    for fmt in EXPORT["save_formats"]:
        path = outdir / f"{panel}.{fmt}"
        kwargs = {}
        if fmt.lower() == "png":
            kwargs["dpi"] = EXPORT["png_dpi"]
        fig.savefig(
            path,
            **kwargs,
        )

    plt.close(fig)


# =============================================================================
# I/O
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
    step14_dir: Path,
) -> Dict[str, object]:
    config = read_json_required(
        step14_dir / "00_run_config.json",
        "Step 14 run config",
    )
    signature = read_json_required(
        step14_dir / "00_run_signature.json",
        "Step 14 run signature",
    )

    if config.get("primary_metric") != "cross_cov":
        raise ValueError(
            "Final plotter requires Step 14 primary_metric='cross_cov'."
        )
    if config.get("primary_support") != "common4":
        raise ValueError(
            "Final plotter requires Step 14 primary_support='common4'."
        )
    if config.get("primary_conditioning") != "xmid_latent":
        raise ValueError(
            "Final plotter requires Step 14 primary_conditioning='xmid_latent'."
        )
    if config.get("operational_class_filter") is not None:
        raise ValueError(
            "Final plotter requires no TT/TF/FT/FF primary filter."
        )
    if config.get("step13_core_redefined") is not False:
        raise ValueError(
            "Final plotter requires inherited Step-13 abundance core."
        )
    if (
        config.get(
            "step13_complete_case_subject_universe_redefined"
        )
        is not False
    ):
        raise ValueError(
            "Final plotter requires inherited Step-13 subject universe."
        )

    payload = signature.get("signature_payload", {})
    params = payload.get("parameters", {})

    if params.get("interval_min_n_source") != (
        "inherited_step13_min_subject_cell_n"
    ):
        raise ValueError(
            "Final plotter requires Step-14 primary interval support "
            "inherited from Step 13."
        )

    if int(
        params.get("min_subject_interval_bin_n", -1)
    ) != int(
        config.get("min_subject_interval_bin_n", -2)
    ):
        raise ValueError(
            "Step-14 run config and signature disagree on interval-bin min n."
        )

    if (
        config.get("analysis_signature")
        != signature.get("analysis_signature")
    ):
        raise ValueError(
            "Step-14 run config/signature analysis signatures disagree."
        )

    cal = config.get("calendar_position_test", {})
    if cal.get("permutation_implementation") != (
        "vectorized; identical within-subject exchangeability design"
    ):
        raise ValueError(
            "Final plotter requires finalized vectorized Step-14 permutation contract."
        )

    subjects = [
        int(x)
        for x in config.get("complete_case_subjects", [])
    ]
    core_bins = [
        int(x)
        for x in config.get("core_bins", [])
    ]

    if not subjects or not core_bins:
        raise ValueError(
            "Step-14 config has empty inherited subject/core definitions."
        )

    return {
        "step14_script_version": config.get("script_version"),
        "step14_analysis_signature": config.get("analysis_signature"),
        "step11_analysis_signature": payload.get(
            "step11_analysis_signature"
        ),
        "step13_analysis_signature": payload.get(
            "step13_analysis_signature"
        ),
        "n_complete_case_subjects": int(len(subjects)),
        "complete_case_subjects": ";".join(map(str, subjects)),
        "n_core_bins": int(len(core_bins)),
        "core_bins": ";".join(map(str, core_bins)),
        "min_subject_interval_bin_n": int(
            config.get("min_subject_interval_bin_n")
        ),
        "n_permutations": int(
            cal.get("n_permutations", 0)
        ),
        "n_bootstrap": int(
            config.get("anchored_control", {}).get(
                "n_bootstrap", 0
            )
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
    step14_dir: Path,
) -> Dict[str, object]:
    contract = validate_final_contract(
        step14_dir
    )
    return {
        "contract": contract,
        "summary": read_csv_required(
            step14_dir
            / "00_analysis_summary.csv",
            "Step-14 analysis summary",
        ),
        "cohort_profiles": read_csv_required(
            step14_dir
            / "05_cohort_interval_profiles.csv",
            "Step-14 cohort interval profiles",
        ),
        "core_tests": read_csv_required(
            step14_dir
            / "06_core_interval_position_tests.csv",
            "Step-14 core interval-position tests",
        ),
        "binwise_tests": read_csv_required(
            step14_dir
            / "07_binwise_interval_position_tests.csv",
            "Step-14 binwise interval-position tests",
        ),
        "anchor_profiles": read_csv_required(
            step14_dir
            / "08_anchored_core_profiles.csv",
            "Step-14 anchored core profiles",
        ),
        "anchor_slopes": read_csv_required(
            step14_dir
            / "09_anchored_core_slope_bootstrap.csv",
            "Step-14 anchored core slope bootstrap",
        ),
        "subject_slopes": read_csv_required(
            step14_dir
            / "10_anchored_subject_slopes.csv",
            "Step-14 anchored subject slopes",
        ),
        "matching": read_csv_required(
            step14_dir
            / "11_anchor_subject_matching_diagnostic.csv",
            "Step-14 anchor subject-matching diagnostic",
        ),
        "binwise_anchor": read_csv_required(
            step14_dir
            / "12_binwise_anchor_slope_bootstrap.csv",
            "Step-14 binwise anchor slopes",
        ),
        "support_audit": read_csv_required(
            step14_dir
            / "14_interval_support_audit.csv",
            "Step-14 interval support audit",
        ),
    }


# =============================================================================
# Main Figure 7A
# =============================================================================

def plot_figure7A(
    cohort_profiles: pd.DataFrame,
    core_tests: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure7A_calendar_position_profiles"

    require_columns(
        cohort_profiles,
        [
            "metric",
            "dt",
            "t0",
            "equal_subject_value",
        ],
        panel,
    )

    require_columns(
        core_tests,
        [
            "metric",
            "dt",
            "test_status",
            "q_bh_rms_position_structure",
        ],
        panel,
    )

    valid_dt = (
        core_tests[
            core_tests["metric"].astype(str).eq(
                "cross_cov"
            )
            & core_tests["test_status"].astype(str).eq(
                "ok"
            )
        ]["dt"]
        .astype(int)
        .tolist()
    )

    d = cohort_profiles[
        cohort_profiles["metric"].astype(str).eq(
            "cross_cov"
        )
        & cohort_profiles["dt"].astype(int).isin(
            valid_dt
        )
    ].copy().sort_values(
        ["dt", "t0"]
    )

    qrows = core_tests[
        core_tests["metric"].astype(str).eq(
            "cross_cov"
        )
        & core_tests["test_status"].astype(str).eq(
            "ok"
        )
    ][
        [
            "dt",
            "q_bh_rms_position_structure",
            "q_bh_abs_calendar_slope",
        ]
    ].copy()

    d = d.merge(
        qrows,
        on="dt",
        how="left",
        validate="many_to_one",
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)

    for dt in sorted(
        d["dt"].astype(int).unique()
    ):
        g = d[
            d["dt"].astype(int).eq(
                int(dt)
            )
        ].copy().sort_values(
            "t0"
        )

        key = f"dt{int(dt)}"
        if key not in COLORS:
            continue

        ax.plot(
            pd.to_numeric(
                g["t0"],
                errors="coerce",
            ),
            pd.to_numeric(
                g["equal_subject_value"],
                errors="coerce",
            ),
            color=COLORS[key],
            lw=LINE_WIDTHS["standard"],
            ls=LINESTYLES[key],
            marker=MARKERS[key],
            ms=MARKER_SIZES["standard"],
            label=(
                f"{dt} week"
                if int(dt) == 1
                else f"{dt} weeks"
            ),
        )

    qrows = qrows.sort_values(
        "dt"
    )

    min_q = float(
        pd.to_numeric(
            qrows["q_bh_rms_position_structure"],
            errors="coerce",
        ).min()
    )
    tested = ",".join(
        str(int(x))
        for x in qrows["dt"].tolist()
    )

    annotation = (
        f"tested lags: {tested} weeks\n"
        rf"min BH $q_{{RMS}}$ = {min_q:.3f}"
    )

    add_annotation(
        ax,
        panel,
        annotation,
    )

    ax.set_xlabel(
        r"Interval start position, $t_0$"
    )
    ax.set_ylabel(
        "Fixed-core cross-replicate covariance"
    )
    ax.set_xticks(
        [1, 2, 3, 4, 5]
    )

    add_panel_legend(
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
# Main Figure 7B
# =============================================================================

def build_forest_source(
    anchor_slopes: pd.DataFrame,
) -> pd.DataFrame:
    d = anchor_slopes[
        anchor_slopes["metric"].astype(str).eq(
            "cross_cov"
        )
    ].copy()

    if d.empty:
        raise ValueError(
            "No cross_cov rows in anchored slope table."
        )

    ref = d.iloc[0]

    rows = [
        {
            "display_order": 0,
            "label": "Step 13 pooled reference",
            "anchor_type": "step13_reference",
            "subject_matching": "reference",
            "n_subjects": np.nan,
            "estimate": float(
                ref[
                    "step13_observed_slope_per_week"
                ]
            ),
            "ci_low": float(
                ref[
                    "step13_bootstrap_q025"
                ]
            ),
            "ci_high": float(
                ref[
                    "step13_bootstrap_q975"
                ]
            ),
        }
    ]

    specs = [
        (
            "common_start",
            "available",
            "Common start, available",
            1,
        ),
        (
            "common_end",
            "available",
            "Common end, available",
            2,
        ),
        (
            "common_start",
            "matched_all_dt",
            "Common start, matched",
            3,
        ),
        (
            "common_end",
            "matched_all_dt",
            "Common end, matched",
            4,
        ),
    ]

    for anchor, matching, label, order in specs:
        g = d[
            d["anchor_type"].astype(str).eq(
                anchor
            )
            & d["subject_matching"].astype(str).eq(
                matching
            )
        ]

        if len(g) != 1:
            raise ValueError(
                f"Expected one row for {anchor}/{matching}; found {len(g)}"
            )

        r = g.iloc[0]

        rows.append(
            {
                "display_order": int(order),
                "label": label,
                "anchor_type": anchor,
                "subject_matching": matching,
                "n_subjects": int(
                    r[
                        "n_subjects_analysis"
                    ]
                ),
                "estimate": float(
                    r[
                        "observed_slope_per_week"
                    ]
                ),
                "ci_low": float(
                    r[
                        "bootstrap_q025"
                    ]
                ),
                "ci_high": float(
                    r[
                        "bootstrap_q975"
                    ]
                ),
            }
        )

    return pd.DataFrame(rows).sort_values(
        "display_order"
    )


def plot_figure7B(
    anchor_slopes: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure7B_anchor_slope_forest"

    d = build_forest_source(
        anchor_slopes
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)
    zero_line(
        ax,
        orientation="vertical",
    )

    y_positions = np.arange(
        len(d)
    )[::-1]

    for y, row in zip(
        y_positions,
        d.itertuples(
            index=False
        ),
    ):
        if row.anchor_type == "step13_reference":
            color = COLORS["reference"]
            marker = "o"
        else:
            color = COLORS[
                row.anchor_type
            ]
            marker = (
                MARKERS["available"]
                if row.subject_matching == "available"
                else MARKERS["matched"]
            )

        ax.plot(
            [
                float(
                    row.ci_low
                ),
                float(
                    row.ci_high
                ),
            ],
            [y, y],
            color=color,
            lw=LINE_WIDTHS["errorbar"],
        )

        ax.plot(
            float(
                row.estimate
            ),
            y,
            marker=marker,
            color=color,
            ms=MARKER_SIZES["point"],
            linestyle="None",
        )

    ax.set_yticks(
        y_positions
    )
    labels = []
    for row in d.itertuples(index=False):
        if np.isfinite(row.n_subjects):
            labels.append(
                f"{row.label} (n={int(row.n_subjects)})"
            )
        else:
            labels.append(
                row.label
            )

    ax.set_yticklabels(
        labels
    )
    ax.set_xlabel(
        r"Temporal slope (week$^{-1}$)"
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
# Main Figure 7C/D
# =============================================================================

def matched_anchor_row(
    anchor_slopes: pd.DataFrame,
    anchor: str,
) -> pd.Series:
    g = anchor_slopes[
        anchor_slopes["metric"].astype(str).eq(
            "cross_cov"
        )
        & anchor_slopes["anchor_type"].astype(str).eq(
            anchor
        )
        & anchor_slopes["subject_matching"].astype(str).eq(
            "matched_all_dt"
        )
    ]

    if len(g) != 1:
        raise ValueError(
            f"Expected one matched cross_cov row for {anchor}; found {len(g)}"
        )

    return g.iloc[0]


def plot_subject_anchor(
    subject_slopes: pd.DataFrame,
    anchor_slopes: pd.DataFrame,
    anchor: str,
    panel: str,
    outdir: Path,
    source_dir: Path,
) -> None:
    d = subject_slopes[
        subject_slopes["metric"].astype(str).eq(
            "cross_cov"
        )
        & subject_slopes["anchor_type"].astype(str).eq(
            anchor
        )
    ].copy().sort_values(
        "subject"
    )

    row = matched_anchor_row(
        anchor_slopes,
        anchor,
    )

    ids = [
        int(x)
        for x in str(
            row[
                "analysis_subject_ids"
            ]
        ).split(";")
        if str(x).strip()
    ]

    d = d[
        d["subject"].astype(int).isin(
            ids
        )
    ].copy().sort_values(
        "subject"
    )

    if d.empty:
        raise ValueError(
            f"No matched subject slopes for {anchor}."
        )

    source = d.copy()
    source[
        "cohort_observed_slope"
    ] = float(
        row[
            "observed_slope_per_week"
        ]
    )
    source[
        "cohort_bootstrap_q025"
    ] = float(
        row[
            "bootstrap_q025"
        ]
    )
    source[
        "cohort_bootstrap_q975"
    ] = float(
        row[
            "bootstrap_q975"
        ]
    )
    source[
        "exact_signflip_two_sided_p"
    ] = float(
        row[
            "signflip_two_sided_p"
        ]
    )
    source[
        "n_negative_subject_slopes"
    ] = int(
        row[
            "signflip_n_negative_subject_slopes"
        ]
    )

    write_source_data(
        source,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)
    zero_line(ax)

    color = COLORS[anchor]

    ax.axhspan(
        float(
            row[
                "bootstrap_q025"
            ]
        ),
        float(
            row[
                "bootstrap_q975"
            ]
        ),
        color=color,
        alpha=ALPHAS["cohort_ci"],
        linewidth=0,
        zorder=0,
    )

    ax.axhline(
        float(
            row[
                "observed_slope_per_week"
            ]
        ),
        color=color,
        lw=LINE_WIDTHS["primary"],
        zorder=1,
    )

    x = np.arange(
        len(d)
    )
    y = pd.to_numeric(
        d[
            "slope_per_week"
        ],
        errors="coerce",
    ).to_numpy(float)

    ax.scatter(
        x,
        y,
        s=MARKER_SIZES["scatter_area"],
        color=COLORS["subject"],
        alpha=ALPHAS["subject_point"],
        zorder=3,
    )

    ax.set_xticks(
        x
    )
    ax.set_xticklabels(
        [
            f"P{int(s):02d}"
            for s in d[
                "subject"
            ].tolist()
        ]
    )
    ax.set_xlabel(
        "Subject"
    )
    ax.set_ylabel(
        r"Subject-specific slope (week$^{-1}$)"
    )

    annotation = (
        f"{int(row['signflip_n_negative_subject_slopes'])}/"
        f"{int(row['signflip_n_subjects'])} slopes < 0\n"
        f"exact sign-flip $P$ = "
        f"{float(row['signflip_two_sided_p']):.4f}"
    )

    add_annotation(
        ax,
        panel,
        annotation,
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
# Supplementary Figure S8A
# =============================================================================

def plot_supp_s8A(
    core_tests: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS8A_core_position_qvalues"

    d = core_tests[
        core_tests["test_status"].astype(str).eq(
            "ok"
        )
    ].copy().sort_values(
        ["metric", "dt"]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)

    metric_specs = [
        (
            "cross_cov",
            "Cross-replicate covariance",
            COLORS["common_start"],
            "o",
        ),
        (
            "same_var_mean",
            "Within-replicate variance",
            COLORS["same_var_mean"],
            "s",
        ),
        (
            "replicate_specific_excess",
            "Replicate-specific excess",
            COLORS["excess"],
            "^",
        ),
    ]

    for metric, label, color, marker in metric_specs:
        g = d[
            d["metric"].astype(str).eq(
                metric
            )
        ].sort_values(
            "dt"
        )

        if g.empty:
            continue

        ax.plot(
            g["dt"],
            g[
                "q_bh_rms_position_structure"
            ],
            color=color,
            lw=LINE_WIDTHS["standard"],
            marker=marker,
            ms=MARKER_SIZES["standard"],
            label=label,
        )

    ax.axhline(
        0.05,
        color=COLORS["zero"],
        lw=LINE_WIDTHS["zero"],
        ls=LINESTYLES["q_threshold"],
        alpha=ALPHAS["q_threshold"],
    )

    ax.set_xlabel(
        "Temporal lag (weeks)"
    )
    ax.set_ylabel(
        r"BH-adjusted $q_{\mathrm{RMS}}$"
    )
    ax.set_xticks(
        [1, 2, 3]
    )

    add_panel_legend(
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
# Supplementary Figure S8B
# =============================================================================

def plot_supp_s8B(
    binwise_tests: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS8B_binwise_position_qvalues"

    d = binwise_tests[
        binwise_tests["metric"].astype(str).eq(
            "cross_cov"
        )
        & binwise_tests["test_status"].astype(str).eq(
            "ok"
        )
    ].copy().sort_values(
        ["dt", "x_center"]
    )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)

    for dt in sorted(
        d["dt"].astype(int).unique()
    ):
        key = f"dt{int(dt)}"
        if key not in COLORS:
            continue

        g = d[
            d["dt"].astype(int).eq(
                int(dt)
            )
        ].sort_values(
            "x_center"
        )

        ax.plot(
            g[
                "x_center"
            ],
            g[
                "q_bh_rms_position_structure"
            ],
            color=COLORS[key],
            lw=LINE_WIDTHS["standard"],
            marker=MARKERS[key],
            ms=MARKER_SIZES["standard"],
            label=(
                f"{dt} week"
                if int(dt) == 1
                else f"{dt} weeks"
            ),
        )

    ax.axhline(
        0.05,
        color=COLORS["zero"],
        lw=LINE_WIDTHS["zero"],
        ls=LINESTYLES["q_threshold"],
        alpha=ALPHAS["q_threshold"],
    )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        r"BH-adjusted $q_{\mathrm{RMS}}$"
    )

    add_panel_legend(
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
# Supplementary Figure S8C/D
# =============================================================================

def plot_binwise_anchor(
    binwise_anchor: pd.DataFrame,
    anchor: str,
    panel: str,
    outdir: Path,
    source_dir: Path,
) -> None:
    d = binwise_anchor[
        binwise_anchor["metric"].astype(str).eq(
            "cross_cov"
        )
        & binwise_anchor["anchor_type"].astype(str).eq(
            anchor
        )
        & binwise_anchor["subject_matching"].astype(str).eq(
            "matched_all_dt"
        )
        & binwise_anchor["test_status"].astype(str).eq(
            "ok"
        )
    ].copy().sort_values(
        "x_center"
    )

    if d.empty:
        raise ValueError(
            f"No matched binwise anchor slopes for {anchor}."
        )

    write_source_data(
        d,
        source_dir,
        panel,
    )

    fig, ax = make_figure(panel)
    zero_line(ax)

    x = pd.to_numeric(
        d[
            "x_center"
        ],
        errors="coerce",
    ).to_numpy(float)
    y = pd.to_numeric(
        d[
            "observed_slope_per_week"
        ],
        errors="coerce",
    ).to_numpy(float)
    lo = pd.to_numeric(
        d[
            "bootstrap_q025"
        ],
        errors="coerce",
    ).to_numpy(float)
    hi = pd.to_numeric(
        d[
            "bootstrap_q975"
        ],
        errors="coerce",
    ).to_numpy(float)

    color = COLORS[
        anchor
    ]

    ax.errorbar(
        x,
        y,
        yerr=np.vstack(
            [
                y - lo,
                hi - y,
            ]
        ),
        fmt="o",
        color=color,
        lw=LINE_WIDTHS["errorbar"],
        ms=MARKER_SIZES["standard"],
        capsize=MARKER_SIZES[
            "errorbar_capsize"
        ],
    )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(
        r"Anchored cross-covariance slope (week$^{-1}$)"
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


def write_plot_contract_audit(
    outdir: Path,
    contract: Mapping[str, object],
    core_tests: pd.DataFrame,
    anchor_slopes: pd.DataFrame,
    support_audit: pd.DataFrame,
    binwise_tests: pd.DataFrame,
) -> None:
    primary_core = core_tests[
        core_tests["metric"].astype(str).eq("cross_cov")
        & core_tests["test_status"].astype(str).eq("ok")
    ].copy()

    matched = anchor_slopes[
        anchor_slopes["metric"].astype(str).eq("cross_cov")
        & anchor_slopes["subject_matching"].astype(str).eq(
            "matched_all_dt"
        )
    ].copy()

    b = binwise_tests[
        binwise_tests["metric"].astype(str).eq("cross_cov")
        & binwise_tests["test_status"].astype(str).eq("ok")
    ].copy()

    row = {
        **contract,
        "n_primary_calendar_tests_valid": int(len(primary_core)),
        "min_primary_calendar_q_rms": (
            float(
                pd.to_numeric(
                    primary_core[
                        "q_bh_rms_position_structure"
                    ],
                    errors="coerce",
                ).min()
            )
            if not primary_core.empty
            else np.nan
        ),
        "n_primary_calendar_q_rms_lt_0_05": int(
            (
                pd.to_numeric(
                    primary_core[
                        "q_bh_rms_position_structure"
                    ],
                    errors="coerce",
                ) < 0.05
            ).sum()
        ) if not primary_core.empty else 0,
        "n_binwise_crosscov_tests_valid": int(len(b)),
        "n_binwise_crosscov_q_rms_lt_0_05": int(
            (
                pd.to_numeric(
                    b["q_bh_rms_position_structure"],
                    errors="coerce",
                ) < 0.05
            ).sum()
        ) if not b.empty else 0,
        "common_start_matched_n": np.nan,
        "common_start_matched_slope": np.nan,
        "common_end_matched_n": np.nan,
        "common_end_matched_slope": np.nan,
        "support_audit_min_primary_interval_bin_n": float(
            pd.to_numeric(
                support_audit[
                    "min_primary_interval_bin_n"
                ],
                errors="coerce",
            ).min()
        ),
    }

    for anchor in ["common_start", "common_end"]:
        g = matched[
            matched["anchor_type"].astype(str).eq(anchor)
        ]
        if len(g) == 1:
            row[f"{anchor}_matched_n"] = int(
                g.iloc[0]["n_subjects_analysis"]
            )
            row[f"{anchor}_matched_slope"] = float(
                g.iloc[0]["observed_slope_per_week"]
            )

    if int(row["support_audit_min_primary_interval_bin_n"]) < int(
        contract["min_subject_interval_bin_n"]
    ):
        raise ValueError(
            "Plotter audit: Step-14 support audit falls below configured "
            "interval-bin minimum."
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
            "Create publication plots for finalized Step-14 "
            "interval-position and anchored-composition controls."
        ),
    )

    p.add_argument(
        "--step14-dir",
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

    step14_dir = (
        args.step14_dir
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
    supplementary_dir = (
        outdir
        / "supplementary_figure"
    )
    source_dir = (
        outdir
        / "figure_source_data"
    )

    data = load_inputs(
        step14_dir
    )

    if args.figure_set in {
        "all",
        "main",
    }:
        plot_figure7A(
            data[
                "cohort_profiles"
            ],
            data[
                "core_tests"
            ],
            main_dir,
            source_dir,
        )

        plot_figure7B(
            data[
                "anchor_slopes"
            ],
            main_dir,
            source_dir,
        )

        plot_subject_anchor(
            data[
                "subject_slopes"
            ],
            data[
                "anchor_slopes"
            ],
            "common_start",
            "Figure7C_common_start_subject_slopes",
            main_dir,
            source_dir,
        )

        plot_subject_anchor(
            data[
                "subject_slopes"
            ],
            data[
                "anchor_slopes"
            ],
            "common_end",
            "Figure7D_common_end_subject_slopes",
            main_dir,
            source_dir,
        )

    if args.figure_set in {
        "all",
        "supplementary",
    }:
        plot_supp_s8A(
            data[
                "core_tests"
            ],
            supplementary_dir,
            source_dir,
        )

        plot_supp_s8B(
            data[
                "binwise_tests"
            ],
            supplementary_dir,
            source_dir,
        )

        plot_binwise_anchor(
            data[
                "binwise_anchor"
            ],
            "common_start",
            "SuppFigS8C_common_start_binwise_slopes",
            supplementary_dir,
            source_dir,
        )

        plot_binwise_anchor(
            data[
                "binwise_anchor"
            ],
            "common_end",
            "SuppFigS8D_common_end_binwise_slopes",
            supplementary_dir,
            source_dir,
        )

    write_plot_contract_audit(
        outdir,
        data["contract"],
        data["core_tests"],
        data["anchor_slopes"],
        data["support_audit"],
        data["binwise_tests"],
    )

    print(
        "\n[DONE] Step 14 publication plots"
    )
    print(
        "Main Figure 7          :",
        main_dir,
    )
    print(
        "Supplementary Figure 8:",
        supplementary_dir,
    )
    print(
        "Exact plotted data     :",
        source_dir,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
