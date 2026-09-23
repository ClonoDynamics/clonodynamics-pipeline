#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_step16_threshold_robustness_details.py
================================================

Publication plotter for finalized ClonoDynamics Step 16 v3.

Final contract:
- finalized Step-15 v6-forward-safe implementation;
- one 7-bin cross-alpha forward common support;
- one 11-bin cross-alpha temporal common core;
- complete fixed support required within every forward/temporal bootstrap draw;
- alpha=0.05 reference outputs externally audited against finalized Step 15.

The plotter is analysis-free: it reads finalized Step-16 tables, restricts
displayed rows only according to already-finalized support flags, and never
reruns transition classification, bootstrap inference, or model fitting.

Recommended main Figure 9
-------------------------
Figure9A_class_composition
    Pooled TT/TF/FT/FF fractions across the operational threshold sweep.

Figure9B_forward_slopes
    Common-support abundance slopes for T0PLUS, TT and TF across alpha, with
    95% biological subject-bootstrap intervals.

Figure9C_TT_crosscov_dt1
    Abundance-resolved TT cross-replicate covariance at dt=1 across alpha.

Figure9D_temporal_slopes
    Cross-alpha common-core temporal slopes for PRIMARY, NON_FF and TT across
    alpha, with 95% joint biological subject-bootstrap intervals.

Recommended Supplementary Figure 10
-----------------------------------
S10A TT retention relative to alpha=0.05.
S10B T0PLUS forward curves across alpha on the fixed common forward support.
S10C TT forward curves across alpha on the fixed common forward support.
S10D TF forward curves across alpha on the fixed common forward support.
S10E paired forward-slope differences versus alpha=0.05.
S10F NON_FF cross_cov abundance curves at dt=1 across alpha.
S10G TT cross_cov paired differences versus alpha=0.05 at dt=1.
S10H TT temporal covariance profiles across alpha on the common temporal core.

Every panel is exported individually as PDF and PNG. The exact plotted rows are
also written as CSV source data.

Typical run
-----------
python3 ./code/05_plotting/04_controls/plot_step16_threshold_robustness_details.py \
    --step16-dir \
    ./dataset_longitudinal_results/16-observation_threshold_robustness \
    --outdir \
    ./figures/step16_observation_threshold_robustness
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
    "tick_label_pt": 7.0,
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
    "Figure9A_class_composition": (3.40, 2.60),
    "Figure9B_forward_slopes": (3.40, 2.60),
    "Figure9C_TT_crosscov_dt1": (3.40, 2.60),
    "Figure9D_temporal_slopes": (3.40, 2.60),
    "SuppFigS10A_TT_retention": (3.40, 2.60),
    "SuppFigS10B_T0PLUS_forward_curves": (3.40, 2.60),
    "SuppFigS10C_TT_forward_curves": (3.40, 2.60),
    "SuppFigS10D_TF_forward_curves": (3.40, 2.60),
    "SuppFigS10E_forward_slope_differences": (3.40, 2.60),
    "SuppFigS10F_NONFF_crosscov_dt1": (3.40, 2.60),
    "SuppFigS10G_TT_crosscov_difference_dt1": (3.40, 2.60),
    "SuppFigS10H_TT_temporal_profiles": (3.40, 2.60),
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
    "ci_band": 0.15,
    "reference_line": 0.75,
    "zero_line": 0.70,
}

AXIS_STYLE = {
    "tick_direction": "out",
    "major_tick_length": 3.5,
    "grid": False,
}

COLORS = {
    # Thresholds
    "a0p1": "#0072B2",
    "a0p05": "#232222",
    "a0p025": "#009E73",
    "a0p01": "#D55E00",

    # Operational classes/domains
    "TT": "#04527F",
    "TF": "#07A82F",
    "FT": "#D00EC9",
    "FF": "#F24949",
    "T0PLUS": "#E62A00",
    "PRIMARY": "#171717",
    "NON_FF": "#009E73",

    "zero": "#777777",
}

LINESTYLES = {
    "reference_alpha": "--",
    "zero": "--",
    "PRIMARY": "-",
    "NON_FF": "-",
    "TT": "-",
    "TF": "-",
    "T0PLUS": "-",
}

MARKERS = {
    "TT": "^",
    "TF": "D",
    "FT": "s",
    "FF": "o",
    "T0PLUS": "s",
    "PRIMARY": "o",
    "NON_FF": "s",
    "a0p1": "o",
    "a0p05": "s",
    "a0p025": "^",
    "a0p01": "D",
}

ALPHA_ORDER = [0.01, 0.025, 0.05, 0.10]
REFERENCE_ALPHA = 0.05

XLIMS = {
    "Figure9A_class_composition": (0.005, 0.105),
    "Figure9B_forward_slopes": (0.005, 0.105),
    "Figure9C_TT_crosscov_dt1": None,
    "Figure9D_temporal_slopes": (0.005, 0.105),
    "SuppFigS10A_TT_retention": (0.005, 0.105),
    "SuppFigS10B_T0PLUS_forward_curves": None,
    "SuppFigS10C_TT_forward_curves": None,
    "SuppFigS10D_TF_forward_curves": None,
    "SuppFigS10E_forward_slope_differences": (0.005, 0.105),
    "SuppFigS10F_NONFF_crosscov_dt1": None,
    "SuppFigS10G_TT_crosscov_difference_dt1": None,
    "SuppFigS10H_TT_temporal_profiles": (0.8, 5.2),
}

YLIMS = {
    "Figure9A_class_composition": (0.0, 1.02),
    "Figure9B_forward_slopes": None,
    "Figure9C_TT_crosscov_dt1": (-0.2, 1.25),
    "Figure9D_temporal_slopes": None,
    "SuppFigS10A_TT_retention": (0.0, 1.05),
    "SuppFigS10B_T0PLUS_forward_curves": None,
    "SuppFigS10C_TT_forward_curves": None,
    "SuppFigS10D_TF_forward_curves": None,
    "SuppFigS10E_forward_slope_differences": None,
    "SuppFigS10F_NONFF_crosscov_dt1": None,
    "SuppFigS10G_TT_crosscov_difference_dt1": None,
    "SuppFigS10H_TT_temporal_profiles": (0.10, 0.55),
}

LEGEND_SPECS = {
    "Figure9A_class_composition": {"loc": "upper right", "ncol": 1},
    "Figure9B_forward_slopes": {"loc": "lower left", "ncol": 1},
    "Figure9C_TT_crosscov_dt1": {"loc": "upper right", "ncol": 1},
    "Figure9D_temporal_slopes": {"loc": "lower right", "ncol": 1},
    "SuppFigS10B_T0PLUS_forward_curves": {"loc": "best", "ncol": 1},
    "SuppFigS10C_TT_forward_curves": {"loc": "best", "ncol": 1},
    "SuppFigS10D_TF_forward_curves": {"loc": "best", "ncol": 1},
    "SuppFigS10E_forward_slope_differences": {"loc": "best", "ncol": 1},
    "SuppFigS10F_NONFF_crosscov_dt1": {"loc": "upper right", "ncol": 1},
    "SuppFigS10G_TT_crosscov_difference_dt1": {"loc": "best", "ncol": 1},
    "SuppFigS10H_TT_temporal_profiles": {"loc": "upper right", "ncol": 1},
}

ANNOTATIONS = {
    "Figure9B_forward_slopes": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "text": "",
        "ha": "left",
    },
    "Figure9C_TT_crosscov_dt1": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "text": r"$\Delta t=1$ week",
        "ha": "left",
    },
    "Figure9D_temporal_slopes": {
        "enabled": True,
        "xy_axes": (0.03, 0.04),
        "text": "",
        "ha": "left",
    },
}

SHOW_TITLES = False


# =============================================================================
# Generic helpers
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
    ax.tick_params(direction=AXIS_STYLE["tick_direction"])
    ax.grid(AXIS_STYLE["grid"])


def make_figure(panel: str) -> Tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(
        figsize=FIGSIZE_INCHES[panel],
        constrained_layout=False,
    )
    style_axis(ax)
    return fig, ax


def apply_limits(ax: plt.Axes, panel: str) -> None:
    if XLIMS.get(panel) is not None:
        ax.set_xlim(*XLIMS[panel])
    if YLIMS.get(panel) is not None:
        ax.set_ylim(*YLIMS[panel])


def zero_line(
    ax: plt.Axes,
    orientation: str = "horizontal",
) -> None:
    kwargs = {
        "color": COLORS["zero"],
        "lw": LINE_WIDTHS["reference"],
        "ls": LINESTYLES["zero"],
        "alpha": ALPHAS["zero_line"],
        "zorder": 0,
    }
    if orientation == "horizontal":
        ax.axhline(0.0, **kwargs)
    else:
        ax.axvline(0.0, **kwargs)


def reference_alpha_line(ax: plt.Axes) -> None:
    ax.axvline(
        REFERENCE_ALPHA,
        color=COLORS["a0p05"],
        lw=LINE_WIDTHS["reference"],
        ls=LINESTYLES["reference_alpha"],
        alpha=ALPHAS["reference_line"],
        zorder=0,
    )


def set_alpha_ticks(ax: plt.Axes) -> None:
    ax.set_xticks(ALPHA_ORDER)
    ax.set_xticklabels(["0.01", "0.025", "0.05", "0.10"])


def add_legend(ax: plt.Axes, panel: str) -> None:
    spec = LEGEND_SPECS.get(panel, {})
    ax.legend(
        frameon=False,
        loc=spec.get("loc", "best"),
        ncol=int(spec.get("ncol", 1)),
        borderaxespad=0.25,
        handlelength=1.5,
        labelspacing=0.25,
    )


def add_annotation(
    ax: plt.Axes,
    panel: str,
    text_override: Optional[str] = None,
) -> None:
    spec = ANNOTATIONS.get(panel, {})
    if not spec.get("enabled", False):
        return
    x, y = spec.get("xy_axes", (0.03, 0.04))
    ax.text(
        x,
        y,
        (
            spec.get("text", "")
            if text_override is None
            else str(text_override)
        ),
        transform=ax.transAxes,
        ha=spec.get("ha", "left"),
        va=spec.get("va", "bottom"),
        fontsize=FONT["annotation_pt"],
    )


def save_figure(
    fig: plt.Figure,
    outdir: Path,
    panel: str,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        fig.tight_layout(pad=0.55)
    except Exception:
        pass

    for fmt in EXPORT["save_formats"]:
        path = outdir / f"{panel}.{fmt}"
        kwargs = {}
        if fmt.lower() == "png":
            kwargs["dpi"] = EXPORT["png_dpi"]
        fig.savefig(path, **kwargs)

    plt.close(fig)


def alpha_key(alpha: float) -> str:
    if np.isclose(alpha, 0.10):
        return "a0p1"
    if np.isclose(alpha, 0.05):
        return "a0p05"
    if np.isclose(alpha, 0.025):
        return "a0p025"
    if np.isclose(alpha, 0.01):
        return "a0p01"
    raise ValueError(f"Unsupported alpha for color mapping: {alpha}")


def alpha_label(alpha: float) -> str:
    if np.isclose(alpha, REFERENCE_ALPHA):
        return r"$\alpha=0.05$ (reference)"
    return rf"$\alpha={alpha:g}$"


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
    step16_dir: Path,
) -> Dict[str, object]:
    config = read_json_required(
        step16_dir / "00_run_config.json",
        "Step 16 run config",
    )
    signature = read_json_required(
        step16_dir / "00_run_signature.json",
        "Step 16 run signature",
    )

    if config.get("analysis_signature") != signature.get("analysis_signature"):
        raise ValueError(
            "Step-16 run config/signature analysis signatures disagree."
        )

    if config.get("scientific_scope") != (
        "operational_observation_threshold_robustness"
    ):
        raise ValueError(
            "Unexpected Step-16 scientific scope."
        )

    if config.get("primary_estimands_changed") is not False:
        raise ValueError(
            "Final Step-16 plotter requires unchanged primary estimands."
        )

    if not str(
        config.get("step15_script_version", "")
    ).startswith("v6-"):
        raise ValueError(
            "Final Step-16 plotter requires finalized Step-15 v6 implementation."
        )

    expected_alphas = [0.1, 0.05, 0.025, 0.01]
    observed_alphas = [
        float(x)
        for x in config.get("alphas", [])
    ]
    if len(observed_alphas) != len(expected_alphas) or not np.allclose(
        sorted(observed_alphas),
        sorted(expected_alphas),
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError(
            "Final plotter expects alpha sweep 0.10, 0.05, 0.025, 0.01."
        )

    if not np.isclose(
        float(config.get("reference_alpha", np.nan)),
        REFERENCE_ALPHA,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError(
            "Final plotter expects reference alpha=0.05."
        )

    boot = config.get("bootstrap", {})
    if boot.get(
        "forward_complete_common_support_required_per_draw"
    ) is not True:
        raise ValueError(
            "Final Step-16 forward bootstrap support contract is missing."
        )
    if boot.get(
        "temporal_complete_common_core_all_dt_required_per_draw"
    ) is not True:
        raise ValueError(
            "Final Step-16 temporal bootstrap support contract is missing."
        )
    if boot.get("same_draws_across_alpha") is not True:
        raise ValueError(
            "Final Step-16 paired cross-alpha bootstrap contract is missing."
        )

    forward = config.get("forward", {})
    forward_bins = [
        int(x)
        for x in forward.get(
            "cross_alpha_common_support_bins",
            [],
        )
    ]

    temporal = config.get("temporal", {})
    temporal_bins = [
        int(x)
        for x in temporal.get(
            "cross_alpha_common_core_bins",
            [],
        )
    ]
    subjects = [
        int(x)
        for x in temporal.get(
            "subject_universe",
            [],
        )
    ]

    if not forward_bins or not temporal_bins or not subjects:
        raise ValueError(
            "Step-16 final run has empty common support/core/subject universe."
        )

    payload = signature.get("signature_payload", {})
    if payload.get(
        "forward_bootstrap_requires_complete_common_support"
    ) is not True:
        raise ValueError(
            "Step-16 signature does not certify fixed forward bootstrap support."
        )
    if payload.get(
        "temporal_bootstrap_requires_complete_common_core_all_dt"
    ) is not True:
        raise ValueError(
            "Step-16 signature does not certify fixed temporal bootstrap support."
        )

    return {
        "step16_script_version": config.get("script_version"),
        "step16_analysis_signature": config.get("analysis_signature"),
        "step15_script_version": config.get("step15_script_version"),
        "reference_alpha": float(config.get("reference_alpha")),
        "n_forward_common_bins": int(len(forward_bins)),
        "forward_common_bins": ";".join(map(str, forward_bins)),
        "n_temporal_common_bins": int(len(temporal_bins)),
        "temporal_common_bins": ";".join(map(str, temporal_bins)),
        "n_temporal_subjects": int(len(subjects)),
        "temporal_subjects": ";".join(map(str, subjects)),
        "n_bootstrap": int(boot.get("n_bootstrap", 0)),
        "step9_analysis_signature": payload.get("step9_analysis_signature"),
        "step11_analysis_signature": payload.get("step11_analysis_signature"),
        "step13_analysis_signature": payload.get("step13_analysis_signature"),
        "step14_analysis_signature": payload.get("step14_analysis_signature"),
        "step15_reference_analysis_signature": payload.get(
            "step15_reference_analysis_signature"
        ),
    }


def read_required(step16_dir: Path, filename: str) -> pd.DataFrame:
    path = step16_dir / filename
    if not path.is_file():
        raise FileNotFoundError(path)
    d = pd.read_csv(path)
    if d.empty:
        raise ValueError(f"Empty Step-16 table: {path}")
    return d


def load_inputs(step16_dir: Path) -> Dict[str, object]:
    contract = validate_final_contract(
        step16_dir
    )
    return {
        "contract": contract,
        "composition": read_required(
            step16_dir,
            "02_class_composition_pooled_by_alpha.csv",
        ),
        "retention": read_required(
            step16_dir,
            "03_tt_retention_by_alpha.csv",
        ),
        "forward_curves": read_required(
            step16_dir,
            "04_forward_domain_cross_alpha_by_bin.csv",
        ),
        "forward_slopes": read_required(
            step16_dir,
            "05_forward_slope_cross_alpha.csv",
        ),
        "forward_support": read_required(
            step16_dir,
            "05b_forward_cross_alpha_common_support.csv",
        ),
        "forward_slope_diff": read_required(
            step16_dir,
            "07_forward_slope_difference_vs_reference_alpha.csv",
        ),
        "fluctuation": read_required(
            step16_dir,
            "08_fluctuation_domain_cross_alpha_by_bin_dt.csv",
        ),
        "fluctuation_diff": read_required(
            step16_dir,
            "09_fluctuation_difference_vs_reference_alpha.csv",
        ),
        "temporal_core": read_required(
            step16_dir,
            "11_cross_alpha_temporal_common_core.csv",
        ),
        "temporal_profiles": read_required(
            step16_dir,
            "12_temporal_profiles_cross_alpha.csv",
        ),
        "temporal_slopes": read_required(
            step16_dir,
            "13_temporal_slopes_cross_alpha.csv",
        ),
        "temporal_slope_diff": read_required(
            step16_dir,
            "14_temporal_slope_difference_vs_reference_alpha.csv",
        ),
        "internal_audit": read_required(
            step16_dir,
            "15_reference_alpha_internal_audit.csv",
        ),
        "external_audit": read_required(
            step16_dir,
            "16_reference_step15_output_audit.csv",
        ),
        "support_summary": read_required(
            step16_dir,
            "17_support_summary_by_alpha.csv",
        ),
    }


# =============================================================================
# Figure 9A
# =============================================================================

def plot_figure9A(
    composition: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure9A_class_composition"
    d = composition.copy().sort_values(["obs_class_step15", "alpha"])

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    reference_alpha_line(ax)

    for cls in ["TT", "TF", "FT", "FF"]:
        g = d[
            d["obs_class_step15"].astype(str).eq(cls)
        ].sort_values("alpha")

        ax.plot(
            g["alpha"],
            g["fraction"],
            color=COLORS[cls],
            marker=MARKERS[cls],
            ms=MARKER_SIZES["standard"],
            lw=LINE_WIDTHS["standard"],
            label=cls,
        )

    ax.set_xlabel(r"Operational observation threshold, $\alpha$")
    ax.set_ylabel("Fraction of transitions")
    set_alpha_ticks(ax)
    add_legend(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 9B
# =============================================================================

def plot_figure9B(
    forward_slopes: pd.DataFrame,
    contract: Mapping[str, object],
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure9B_forward_slopes"

    d = forward_slopes.copy().sort_values(["mode", "alpha"])
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)
    reference_alpha_line(ax)

    for mode in ["T0PLUS", "TT", "TF"]:
        g = d[
            d["mode"].astype(str).eq(mode)
        ].sort_values("alpha")

        x = g["alpha"].to_numpy(float)
        y = g["observed_slope"].to_numpy(float)
        lo = g["bootstrap_q025"].to_numpy(float)
        hi = g["bootstrap_q975"].to_numpy(float)

        ax.errorbar(
            x,
            y,
            yerr=np.vstack([y - lo, hi - y]),
            color=COLORS[mode],
            marker=MARKERS[mode],
            ms=MARKER_SIZES["standard"],
            lw=LINE_WIDTHS["standard"],
            elinewidth=LINE_WIDTHS["errorbar"],
            capsize=MARKER_SIZES["errorbar_capsize"],
            label=mode,
        )

    ax.set_xlabel(r"Operational observation threshold, $\alpha$")
    ax.set_ylabel(r"Forward abundance slope")
    set_alpha_ticks(ax)
    add_legend(ax, panel)
    add_annotation(
        ax,
        panel,
        text_override=(
            f"{int(contract['n_forward_common_bins'])}-bin "
            "cross-alpha common support"
        ),
    )
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 9C / cross-covariance curves
# =============================================================================

def plot_crosscov_alpha_curves(
    fluctuation: pd.DataFrame,
    mode: str,
    panel: str,
    outdir: Path,
    source_dir: Path,
    dt: int = 1,
) -> None:
    mask = (
        fluctuation["mode"].astype(str).eq(mode)
        & fluctuation["dt"].astype(int).eq(int(dt))
    )

    if "meets_min_support" in fluctuation.columns:
        mask &= bool_series(
            fluctuation["meets_min_support"]
        )

    d = fluctuation[
        mask
    ].copy().sort_values(["alpha", "x_center"])

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)

    for alpha in ALPHA_ORDER[::-1]:
        g = d[
            np.isclose(
                d["alpha"].to_numpy(float),
                float(alpha),
                rtol=0.0,
                atol=1e-12,
            )
        ].sort_values("x_center")

        if g.empty:
            continue

        key = alpha_key(alpha)

        ax.plot(
            g["x_center"],
            g["cross_cov"],
            color=COLORS[key],
            marker=MARKERS[key],
            ms=MARKER_SIZES["standard"],
            lw=(
                LINE_WIDTHS["primary"]
                if np.isclose(alpha, REFERENCE_ALPHA)
                else LINE_WIDTHS["standard"]
            ),
            label=alpha_label(alpha),
        )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel("Cross-replicate covariance")
    add_legend(ax, panel)
    add_annotation(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 9D
# =============================================================================

def plot_figure9D(
    temporal_slopes: pd.DataFrame,
    contract: Mapping[str, object],
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure9D_temporal_slopes"

    d = temporal_slopes.copy().sort_values(["mode", "alpha"])
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)
    reference_alpha_line(ax)

    for mode in ["PRIMARY", "NON_FF", "TT"]:
        g = d[
            d["mode"].astype(str).eq(mode)
        ].sort_values("alpha")

        x = g["alpha"].to_numpy(float)
        y = g["observed_slope_per_week"].to_numpy(float)
        lo = g["bootstrap_q025"].to_numpy(float)
        hi = g["bootstrap_q975"].to_numpy(float)

        ax.errorbar(
            x,
            y,
            yerr=np.vstack([y - lo, hi - y]),
            color=COLORS[mode],
            marker=MARKERS[mode],
            ms=MARKER_SIZES["standard"],
            lw=(
                LINE_WIDTHS["primary"]
                if mode == "PRIMARY"
                else LINE_WIDTHS["standard"]
            ),
            elinewidth=LINE_WIDTHS["errorbar"],
            capsize=MARKER_SIZES["errorbar_capsize"],
            label=mode,
        )

    ax.set_xlabel(r"Operational observation threshold, $\alpha$")
    ax.set_ylabel(r"Temporal covariance slope (week$^{-1}$)")
    set_alpha_ticks(ax)
    add_legend(ax, panel)
    add_annotation(
        ax,
        panel,
        text_override=(
            f"{int(contract['n_temporal_common_bins'])}-bin "
            "cross-alpha temporal core\n"
            f"{int(contract['n_temporal_subjects'])} subjects"
        ),
    )
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary S10A
# =============================================================================

def plot_s10A(
    retention: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS10A_TT_retention"

    d = retention.copy().sort_values("alpha")
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    reference_alpha_line(ax)

    ax.plot(
        d["alpha"],
        d["fraction_reference_tt_retained"],
        color=COLORS["TT"],
        marker=MARKERS["TT"],
        ms=MARKER_SIZES["standard"],
        lw=LINE_WIDTHS["standard"],
    )

    ax.set_xlabel(r"Operational observation threshold, $\alpha$")
    ax.set_ylabel("Fraction of reference TT retained")
    set_alpha_ticks(ax)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary forward curves
# =============================================================================

def plot_forward_mode_curves(
    forward_curves: pd.DataFrame,
    forward_support: pd.DataFrame,
    mode: str,
    panel: str,
    outdir: Path,
    source_dir: Path,
) -> None:
    common_bins = (
        forward_support.loc[
            forward_support[
                "selected_cross_alpha_forward_common_support"
            ].astype(bool),
            "bin",
        ]
        .astype(int)
        .tolist()
    )

    d = forward_curves[
        forward_curves["mode"].astype(str).eq(mode)
        & forward_curves["bin"].astype(int).isin(common_bins)
    ].copy().sort_values(["alpha", "bin"])

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)

    for alpha in ALPHA_ORDER[::-1]:
        g = d[
            np.isclose(
                d["alpha"].to_numpy(float),
                float(alpha),
                rtol=0.0,
                atol=1e-12,
            )
        ].sort_values("x_center")

        if g.empty:
            continue

        key = alpha_key(alpha)

        ax.plot(
            g["x_center"],
            g["mean_dx"],
            color=COLORS[key],
            marker=MARKERS[key],
            ms=MARKER_SIZES["standard"],
            lw=(
                LINE_WIDTHS["primary"]
                if np.isclose(alpha, REFERENCE_ALPHA)
                else LINE_WIDTHS["standard"]
            ),
            label=alpha_label(alpha),
        )

    ax.set_xlabel(r"Initial observed log-frequency, $x_0$")
    ax.set_ylabel(rf"{mode} mean replicate-decoupled displacement")
    add_legend(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary S10E
# =============================================================================

def plot_s10E(
    forward_slope_diff: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS10E_forward_slope_differences"

    d = forward_slope_diff.copy().sort_values(["mode", "alpha"])
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)
    reference_alpha_line(ax)

    for mode in ["T0PLUS", "TT", "TF"]:
        g = d[
            d["mode"].astype(str).eq(mode)
        ].sort_values("alpha")

        x = g["alpha"].to_numpy(float)
        y = g["bootstrap_slope_difference_median"].to_numpy(float)
        lo = g["bootstrap_slope_difference_q025"].to_numpy(float)
        hi = g["bootstrap_slope_difference_q975"].to_numpy(float)

        ax.errorbar(
            x,
            y,
            yerr=np.vstack([y - lo, hi - y]),
            color=COLORS[mode],
            marker=MARKERS[mode],
            ms=MARKER_SIZES["standard"],
            lw=LINE_WIDTHS["standard"],
            elinewidth=LINE_WIDTHS["errorbar"],
            capsize=MARKER_SIZES["errorbar_capsize"],
            label=mode,
        )

    ax.set_xlabel(r"Operational observation threshold, $\alpha$")
    ax.set_ylabel(r"Slope difference versus $\alpha=0.05$")
    set_alpha_ticks(ax)
    add_legend(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary S10G
# =============================================================================

def plot_s10G(
    fluctuation_diff: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS10G_TT_crosscov_difference_dt1"

    d = fluctuation_diff[
        fluctuation_diff["mode"].astype(str).eq("TT")
        & fluctuation_diff["metric"].astype(str).eq("cross_cov")
        & fluctuation_diff["dt"].astype(int).eq(1)
    ].copy().sort_values(["alpha", "x_center"])

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)

    for alpha in [0.10, 0.025, 0.01]:
        g = d[
            np.isclose(
                d["alpha"].to_numpy(float),
                float(alpha),
                rtol=0.0,
                atol=1e-12,
            )
        ].sort_values("x_center")

        if g.empty:
            continue

        key = alpha_key(alpha)

        ax.plot(
            g["x_center"],
            g["difference"],
            color=COLORS[key],
            marker=MARKERS[key],
            ms=MARKER_SIZES["standard"],
            lw=LINE_WIDTHS["standard"],
            label=alpha_label(alpha),
        )

        ax.fill_between(
            g["x_center"].to_numpy(float),
            g["bootstrap_difference_q025"].to_numpy(float),
            g["bootstrap_difference_q975"].to_numpy(float),
            color=COLORS[key],
            alpha=ALPHAS["ci_band"],
            linewidth=0,
        )

    ax.set_xlabel(
        r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$"
    )
    ax.set_ylabel(r"TT cross-covariance difference vs $\alpha=0.05$")
    add_legend(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary S10H
# =============================================================================

def plot_s10H(
    temporal_profiles: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "SuppFigS10H_TT_temporal_profiles"

    d = temporal_profiles[
        temporal_profiles["mode"].astype(str).eq("TT")
    ].copy().sort_values(["alpha", "dt"])

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)

    for alpha in ALPHA_ORDER[::-1]:
        g = d[
            np.isclose(
                d["alpha"].to_numpy(float),
                float(alpha),
                rtol=0.0,
                atol=1e-12,
            )
        ].sort_values("dt")

        key = alpha_key(alpha)

        ax.plot(
            g["dt"],
            g["cross_cov"],
            color=COLORS[key],
            marker=MARKERS[key],
            ms=MARKER_SIZES["standard"],
            lw=(
                LINE_WIDTHS["primary"]
                if np.isclose(alpha, REFERENCE_ALPHA)
                else LINE_WIDTHS["standard"]
            ),
            label=alpha_label(alpha),
        )

    ax.set_xlabel("Temporal lag (weeks)")
    ax.set_ylabel("TT common-core cross-replicate covariance")
    ax.set_xticks([1, 2, 3, 4, 5])
    add_legend(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


def write_plot_contract_audit(
    outdir: Path,
    contract: Mapping[str, object],
    forward_support: pd.DataFrame,
    temporal_core: pd.DataFrame,
    temporal_slopes: pd.DataFrame,
    internal_audit: pd.DataFrame,
    external_audit: pd.DataFrame,
    support_summary: pd.DataFrame,
) -> None:
    forward_selected = forward_support[
        bool_series(
            forward_support[
                "selected_cross_alpha_forward_common_support"
            ]
        )
    ].copy()

    temporal_selected = temporal_core[
        bool_series(
            temporal_core[
                "selected_cross_alpha_common_core"
            ]
        )
    ].copy()

    identity = external_audit[
        external_audit["component"].astype(str).isin(
            [
                "forward_domain_by_bin",
                "fluctuation_domain_by_bin_dt",
            ]
        )
    ].copy()

    max_external_identity_diff = (
        float(
            pd.to_numeric(
                identity["max_abs_difference"],
                errors="coerce",
            ).max()
        )
        if not identity.empty
        else np.nan
    )

    internal_max_abs_diff = float(
        pd.to_numeric(
            internal_audit["absolute_difference"],
            errors="coerce",
        ).max()
    )

    primary = temporal_slopes[
        temporal_slopes["mode"].astype(str).eq("PRIMARY")
    ].copy()

    primary_slope_range = (
        float(
            pd.to_numeric(
                primary["observed_slope_per_week"],
                errors="coerce",
            ).max()
            - pd.to_numeric(
                primary["observed_slope_per_week"],
                errors="coerce",
            ).min()
        )
        if not primary.empty
        else np.nan
    )

    row = {
        **contract,
        "forward_support_rows_selected": int(len(forward_selected)),
        "temporal_core_rows_selected": int(len(temporal_selected)),
        "internal_reference_audit_max_abs_difference": internal_max_abs_diff,
        "external_step15_identity_audit_max_abs_difference":
            max_external_identity_diff,
        "primary_temporal_slope_range_across_alpha": primary_slope_range,
        "support_summary_min_frozen_core_bins_valid_all_modes_all_dt": int(
            pd.to_numeric(
                support_summary[
                    "n_frozen_core_bins_valid_all_modes_all_dt"
                ],
                errors="coerce",
            ).min()
        ),
    }

    if int(row["forward_support_rows_selected"]) != int(
        contract["n_forward_common_bins"]
    ):
        raise ValueError(
            "Plotter audit: forward common-support count mismatch."
        )

    if int(row["temporal_core_rows_selected"]) != int(
        contract["n_temporal_common_bins"]
    ):
        raise ValueError(
            "Plotter audit: temporal common-core count mismatch."
        )

    if np.isfinite(max_external_identity_diff) and (
        max_external_identity_diff > 1e-10
    ):
        raise ValueError(
            "Plotter audit: alpha=0.05 no longer reproduces finalized Step 15."
        )

    if internal_max_abs_diff != 0:
        raise ValueError(
            "Plotter audit: Step-16 internal reference-alpha audit failed."
        )

    if np.isfinite(primary_slope_range) and primary_slope_range > 1e-12:
        raise ValueError(
            "Plotter audit: PRIMARY temporal slope is not alpha invariant."
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
            "Plot finalized Step-16 operational observation-threshold robustness."
        ),
    )

    p.add_argument(
        "--step16-dir",
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
        choices=["all", "main", "supplementary"],
        default="all",
    )

    return p


def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(argv)
    configure_matplotlib()

    step16_dir = args.step16_dir.expanduser().resolve(strict=True)
    outdir = args.outdir.expanduser().resolve()

    main_dir = outdir / "main_figure"
    supp_dir = outdir / "supplementary_figure"
    source_dir = outdir / "figure_source_data"

    data = load_inputs(step16_dir)

    if args.figure_set in {"all", "main"}:
        plot_figure9A(
            data["composition"],
            main_dir,
            source_dir,
        )
        plot_figure9B(
            data["forward_slopes"],
            data["contract"],
            main_dir,
            source_dir,
        )
        plot_crosscov_alpha_curves(
            data["fluctuation"],
            "TT",
            "Figure9C_TT_crosscov_dt1",
            main_dir,
            source_dir,
            dt=1,
        )
        plot_figure9D(
            data["temporal_slopes"],
            data["contract"],
            main_dir,
            source_dir,
        )

    if args.figure_set in {"all", "supplementary"}:
        plot_s10A(
            data["retention"],
            supp_dir,
            source_dir,
        )
        plot_forward_mode_curves(
            data["forward_curves"],
            data["forward_support"],
            "T0PLUS",
            "SuppFigS10B_T0PLUS_forward_curves",
            supp_dir,
            source_dir,
        )
        plot_forward_mode_curves(
            data["forward_curves"],
            data["forward_support"],
            "TT",
            "SuppFigS10C_TT_forward_curves",
            supp_dir,
            source_dir,
        )
        plot_forward_mode_curves(
            data["forward_curves"],
            data["forward_support"],
            "TF",
            "SuppFigS10D_TF_forward_curves",
            supp_dir,
            source_dir,
        )
        plot_s10E(
            data["forward_slope_diff"],
            supp_dir,
            source_dir,
        )
        plot_crosscov_alpha_curves(
            data["fluctuation"],
            "NON_FF",
            "SuppFigS10F_NONFF_crosscov_dt1",
            supp_dir,
            source_dir,
            dt=1,
        )
        plot_s10G(
            data["fluctuation_diff"],
            supp_dir,
            source_dir,
        )
        plot_s10H(
            data["temporal_profiles"],
            supp_dir,
            source_dir,
        )

    write_plot_contract_audit(
        outdir,
        data["contract"],
        data["forward_support"],
        data["temporal_core"],
        data["temporal_slopes"],
        data["internal_audit"],
        data["external_audit"],
        data["support_summary"],
    )

    print("\n[DONE] Step 16 publication plotting complete")
    print("Main Figure 9:", main_dir)
    print("Supplementary Figure 10:", supp_dir)
    print("Exact plotted source data:", source_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
