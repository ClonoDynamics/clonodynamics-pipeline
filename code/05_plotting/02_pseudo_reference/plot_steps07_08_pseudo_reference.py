#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_steps07_08_pseudo_reference.py
========================================

Publication plotter for the autonomous pseudo-longitudinal technical-null
chapter of ClonoDynamics.

This plotter is compatible with the compact-native equal-weight AB/BA Step 7
and the fast pair-transition-bank Step 8.

The script consumes the FINAL Step-7 and Step-8 summary outputs and writes
SINGLE-PANEL figures only. Panels can then be assembled externally into the
main and supplementary figures.

Recommended main figure
-----------------------
Figure 3A  Primary observed replicate-decoupled forward drift (equal-weight AB/BA)
Figure 3B  Matched same-measure versus cross-replicate forward drift on common4
Figure 3C  Cross-replicate displacement covariance versus abundance across pseudo-lags
Figure 3D  Shared covariance versus within-replicate variance at pseudo-lag 1
Figure 3E  Pooled versus equal-interval cross-covariance pseudo-lag slopes

Recommended supplementary figure
--------------------------------
Figure S4A Conditioning-geometry diagnostic (x0 / xmid / xstar)
Figure S4B common4 support across pseudo-lags
Figure S4C pooled cross-covariance across pseudo-lags

The plotter also writes the exact source data used for each panel to
`<outdir>/figure_source_data/`.

No biological confidence intervals are shown. All 2.5th–97.5th intervals are
randomization intervals across pseudo configurations.

Typical usage
-------------
python3 plot_steps07_08_pseudo_reference.py \
    --step7-dir ./dataset_longitudinal_results/7-pseudo_forward_technical_null_compact \
    --step8-dir ./dataset_longitudinal_results/8-pseudo_fluctuation_conditioning_validation \
    --outdir ./figures/pseudo_technical_null

All dimensions, line widths, marker sizes, fonts, legend placement, y-limits,
and colors are editable near the top of the script.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# =============================================================================
# USER-EDITABLE PUBLICATION AESTHETICS
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

# Single-panel dimensions in inches.
FIGSIZE_INCHES: Dict[str, Tuple[float, float]] = {
    "Figure3A_forward_primary": (3.2, 2.4),
    "Figure3B_forward_matched": (3.2, 2.4),
    "Figure3C_cross_cov_by_lag": (3.2, 2.4),
    "Figure3D_fluctuation_decomposition_dt1": (3.2, 2.4),
    "Figure3E_pseudolag_slopes": (3.2, 2.4),
    "FigureS4A_conditioning_geometry": (3.2, 2.4),
    "FigureS4B_common4_support": (3.2, 2.4),
    "FigureS4C_pooled_cross_cov_by_lag": (3.2, 2.4),
}

# Line / marker constants.
LINE_WIDTH = 1
PRIMARY_LINE_WIDTH = 0.7
SECONDARY_LINE_WIDTH = 0.5
ZERO_LINE_WIDTH = 0.5
ERRORBAR_LINE_WIDTH = 0.5
CAP_SIZE = 2.0
MARKER_SIZE = 2.5
POINT_SIZE = 3.0

RIBBON_ALPHA = 0.15
SECONDARY_ALPHA = 0.45
ZERO_LINE_ALPHA = 0.65

# Publication colors. Color-blind-safe / high-contrast choices.
COLORS = {
    "cross": "#0072B2",
    "same": "#D55E00",
    "combined": "#222222",
    "excess": "#CC79A7",
    "x0": "#D55E00",
    "xmid": "#0072B2",
    "xstar": "#009E73",
    "support": "#333333",
    "lag1": "#0072B2",
    "lag2": "#009E73",
    "lag3": "#E69F00",
    "lag4": "#CC79A7",
    "lag5": "#D55E00",
}

LINESTYLES = {
    "AB": "-",
    "BA": "--",
    "combined": "-",
    "cross": "-",
    "same": "-",
    "excess": "--",
    "x0": "-",
    "xmid": "-",
    "xstar": "--",
    "lag1": "-",
    "lag2": "--",
    "lag3": "-.",
    "lag4": ":",
    "lag5": (0, (5, 2)),
}

MARKERS = {
    "AB": "o",
    "BA": "s",
    "combined": None,
    "cross": "o",
    "same": "s",
    "excess": "^",
    "lag1": "o",
    "lag2": "s",
    "lag3": "^",
    "lag4": "D",
    "lag5": "v",
}

# Plot titles are off by default because panel letters / titles are normally
# added only during final figure assembly.
SHOW_PLOT_TITLES = False

PLOT_TITLES = {
    "Figure3A_forward_primary": "Replicate-decoupled forward technical null",
    "Figure3B_forward_matched": "Matched same- versus cross-replicate conditioning",
    "Figure3C_cross_cov_by_lag": "Shared fluctuation covariance across pseudo-lags",
    "Figure3D_fluctuation_decomposition_dt1": "Shared versus within-replicate fluctuation",
    "Figure3E_pseudolag_slopes": "Shared-covariance pseudo-lag robustness",
    "FigureS4A_conditioning_geometry": "Unrestricted conditioning-geometry diagnostic",
    "FigureS4B_common4_support": "common4 support across pseudo-lags",
    "FigureS4C_pooled_cross_cov_by_lag": "Pooled cross-covariance across pseudo-lags",
}

# Axis limits are deliberately centralized and editable. Leave as None for auto.
YLIMS = {
    "Figure3A_forward_primary": None,
    "Figure3B_forward_matched": None,
    "Figure3C_cross_cov_by_lag": (-0.1, 0.1),
    "Figure3D_fluctuation_decomposition_dt1": None,
    "Figure3E_pseudolag_slopes": None,
    "FigureS4A_conditioning_geometry": (-0.6, 0.4),
    "FigureS4B_common4_support": (-0.05, 0.2),
    "FigureS4C_pooled_cross_cov_by_lag": (-0.5, 0.8),
}

XLIMS = {
    "Figure3A_forward_primary": None,
    "Figure3B_forward_matched": None,
    "Figure3C_cross_cov_by_lag": None,
    "Figure3D_fluctuation_decomposition_dt1": None,
    "Figure3E_pseudolag_slopes": None,
    "FigureS4A_conditioning_geometry": None,
    "FigureS4B_common4_support": (0.7, 5.3),
    "FigureS4C_pooled_cross_cov_by_lag": (0.7, 5.3),
}

# Step-7 publication guard: a displayed abundance bin must be supported by
# at least this fraction of the full randomization ensemble for the combined
# estimator. With N=2000 and 0.80, the current data exclude only the poorly
# supported sixth forward bin (767/2000 configurations).
MIN_STEP7_CONFIG_FRACTION_PER_BIN = 0.80

# Show the canonical AB/BA folds as muted diagnostics in Figure 3A. The
# scientific estimator is the equal-weight combined curve.
SHOW_AB_BA_FOLDS_MAIN = False

# For abundance-resolved Step-8 curves, retain a bin only when at least this
# fraction of the randomization ensemble provides a valid estimate. Unsupported
# bins are kept as explicit gaps in the plotted curve rather than connected over.
MIN_STEP8_CONFIG_FRACTION_PER_BIN = 0.50

# Ribbons can be toggled independently.
SHOW_PRIMARY_FORWARD_RIBBON = True
SHOW_MATCHED_FORWARD_RIBBONS = True
SHOW_CROSS_COV_LAG_RIBBONS = False
SHOW_DECOMPOSITION_RIBBONS = False

# Save both vector PDF and high-resolution PNG.
SAVE_FORMATS = ("pdf", "png")

# Legends outside the axes.
LEGEND_BBOX = (1.02, 1.0)
LEGEND_LOC = "upper left"


# =============================================================================
# Styling helpers
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


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")
    ax.grid(False)


def zero_line(ax: plt.Axes, orientation: str = "horizontal") -> None:
    kwargs = dict(
        color="#666666",
        lw=ZERO_LINE_WIDTH,
        alpha=ZERO_LINE_ALPHA,
        ls="--",
        zorder=0,
    )
    if orientation == "horizontal":
        ax.axhline(0.0, **kwargs)
    else:
        ax.axvline(0.0, **kwargs)


def apply_limits(ax: plt.Axes, panel_name: str) -> None:
    xlim = XLIMS.get(panel_name)
    ylim = YLIMS.get(panel_name)

    if xlim is not None:
        lo, hi = xlim
        current = ax.get_xlim()
        ax.set_xlim(
            current[0] if lo is None else lo,
            current[1] if hi is None else hi,
        )

    if ylim is not None:
        lo, hi = ylim
        current = ax.get_ylim()
        ax.set_ylim(
            current[0] if lo is None else lo,
            current[1] if hi is None else hi,
        )


def maybe_title(ax: plt.Axes, panel_name: str) -> None:
    if SHOW_PLOT_TITLES:
        ax.set_title(PLOT_TITLES.get(panel_name, panel_name), pad=8)


def save_figure(
    fig: plt.Figure,
    outdir: Path,
    panel_name: str,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for fmt in SAVE_FORMATS:
        path = outdir / f"{panel_name}.{fmt}"
        kwargs = dict(
            bbox_inches="tight",
            pad_inches=0.06,
        )
        if fmt.lower() == "png":
            kwargs["dpi"] = PNG_DPI
        fig.savefig(path, **kwargs)
    plt.close(fig)


def make_figure(panel_name: str) -> Tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(
        figsize=FIGSIZE_INCHES[panel_name],
        constrained_layout=False,
    )
    style_axis(ax)
    return fig, ax


# =============================================================================
# IO / validation
# =============================================================================

def read_nonempty_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size == 0:
        raise ValueError(f"CSV is empty: {path}")
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError(f"CSV has no rows: {path}")
    return frame


def require_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
    source_name: str,
) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(
            f"{source_name}: missing required columns {missing}"
        )


def write_source_data(
    frame: pd.DataFrame,
    source_dir: Path,
    panel_name: str,
) -> None:
    source_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        source_dir / f"{panel_name}_source_data.csv",
        index=False,
    )



def validate_step7_contract(step7_dir: Path) -> None:
    """
    Validate the compact-native Step-7 publication contract.

    The run signature is used when available, but plotting remains possible
    from the final CSV outputs alone.
    """
    signature = step7_dir / "00_run_signature.json"
    if signature.is_file():
        import json
        obj = json.loads(signature.read_text(encoding="utf-8"))
        payload = obj.get("signature_payload", {})
        policy = payload.get("AB_BA_combination_policy")
        row_weighting = payload.get("AB_BA_row_count_weighting")

        if policy != "exact_equal_fold_weight_after_binning":
            raise ValueError(
                "Step 7 run does not declare exact equal-weight AB/BA "
                f"combination after binning: {policy!r}"
            )
        if row_weighting is not False:
            raise ValueError(
                "Step 7 run does not explicitly disable AB/BA row-count "
                f"weighting: {row_weighting!r}"
            )


def step7_supported_bins(
    envelope: pd.DataFrame,
    section: str,
    combined_representation: str,
) -> Tuple[set, int, int]:
    """
    Return bins retained for publication from the combined estimator.

    The threshold is a fraction of the full Step-7 ensemble, inferred from the
    maximum `n_configurations` in the envelope. The same retained-bin set is
    applied to all representations in a panel so AB/BA/same/cross are compared
    on identical abundance support.
    """
    require_columns(
        envelope,
        ["section", "representation", "bin", "n_configurations"],
        "Step 7 curve envelope",
    )

    section_rows = envelope[
        envelope["section"].astype(str).eq(section)
    ].copy()
    if section_rows.empty:
        raise ValueError(f"No Step-7 rows found for section {section!r}")

    n_total = int(
        pd.to_numeric(
            section_rows["n_configurations"], errors="coerce"
        ).max()
    )
    if n_total <= 0:
        raise ValueError(f"Could not infer Step-7 ensemble size for {section}")

    threshold = int(np.ceil(MIN_STEP7_CONFIG_FRACTION_PER_BIN * n_total))

    combined = section_rows[
        section_rows["representation"].astype(str).eq(combined_representation)
    ].copy()
    if combined.empty:
        raise ValueError(
            f"No {combined_representation!r} rows found in Step-7 section "
            f"{section!r}"
        )

    nconf = pd.to_numeric(combined["n_configurations"], errors="coerce")
    bins = set(
        pd.to_numeric(
            combined.loc[nconf >= threshold, "bin"], errors="coerce"
        ).dropna().astype(int)
    )
    if not bins:
        raise ValueError(
            f"No Step-7 bins meet publication support threshold "
            f"{threshold}/{n_total}"
        )

    return bins, n_total, threshold


def mark_step7_support(
    frame: pd.DataFrame,
    supported_bins: set,
    n_total: int,
    threshold: int,
) -> pd.DataFrame:
    out = frame.copy()
    out["plot_supported"] = (
        pd.to_numeric(out["bin"], errors="coerce")
        .fillna(-1)
        .astype(int)
        .isin(supported_bins)
    )
    out["step7_total_configurations"] = int(n_total)
    out["step7_min_configurations_for_plot"] = int(threshold)
    return out


def validate_step8_contract(step8_dir: Path) -> None:
    """Validate the fast pair-transition-bank Step-8 publication contract."""
    config_path = step8_dir / "00_run_config.json"
    if not config_path.is_file():
        return

    import json
    config = json.loads(config_path.read_text(encoding="utf-8"))

    if config.get("primary_support") != "common4":
        raise ValueError(
            "Step 8 run is not common4-based: "
            f"{config.get('primary_support')!r}"
        )
    if config.get("primary_conditioning") != "xmid_latent":
        raise ValueError(
            "Step 8 primary conditioning is not xmid_latent: "
            f"{config.get('primary_conditioning')!r}"
        )

    estimands = config.get("primary_estimands", {})
    if estimands.get("cross_cov_truncated_at_zero") is not False:
        raise ValueError(
            "Step 8 must retain signed cross_cov without clipping."
        )

    for key in (
        "step2_aggregate_only_rerun",
        "step4_rerun",
        "step5_rerun",
    ):
        if key in config and config.get(key) is not False:
            raise ValueError(
                f"Step 8 unexpectedly reports {key}=True"
            )


def load_inputs(
    step7_dir: Path,
    step8_dir: Path,
) -> Dict[str, pd.DataFrame]:
    """
    Load only files actually consumed by this plotter.

    Step 7 is read from the compact-native ensemble envelope. The old
    `06_matched_same_vs_cross_by_configuration.csv` dependency is intentionally
    removed because it is not used for plotting and no longer exists under that
    name in the current Step-7 contract.
    """
    validate_step7_contract(step7_dir)
    validate_step8_contract(step8_dir)

    inputs = {
        "s7_envelope": read_nonempty_csv(
            step7_dir / "04_curve_randomization_envelope.csv"
        ),
        "s8_dt_summary": read_nonempty_csv(
            step8_dir / "02_ensemble_dt_summary.csv"
        ),
        "s8_binned": read_nonempty_csv(
            step8_dir / "04_binned_randomization_envelope.csv"
        ),
        "s8_slopes": read_nonempty_csv(
            step8_dir / "06_temporal_slope_randomization_summary.csv"
        ),
        "s8_support": read_nonempty_csv(
            step8_dir / "09_support_by_configuration_dt.csv"
        ),
    }
    return inputs


# =============================================================================
# Figure 3A — primary observed replicate-decoupled forward
# =============================================================================

def plot_figure3A_forward_primary(
    s7_envelope: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure3A_forward_primary"

    require_columns(
        s7_envelope,
        [
            "section",
            "representation",
            "bin",
            "n_configurations",
            "x_center_median",
            "mean_dx_median",
            "mean_dx_randomization_q025",
            "mean_dx_randomization_q975",
        ],
        "Step 7 curve envelope",
    )

    reps = ["cross_AB", "cross_BA", "cross_combined"]
    d = s7_envelope[
        s7_envelope["section"].astype(str).eq("primary_forward")
        & s7_envelope["representation"].astype(str).isin(reps)
    ].copy()

    if d.empty:
        raise ValueError("No primary_forward rows found in Step 7 envelope")

    supported_bins, n_total, threshold = step7_supported_bins(
        s7_envelope,
        section="primary_forward",
        combined_representation="cross_combined",
    )
    d = mark_step7_support(d, supported_bins, n_total, threshold)
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)

    specs = {
        "cross_AB": (
            "AB fold",
            COLORS["cross"],
            LINESTYLES["AB"],
            MARKERS["AB"],
            SECONDARY_LINE_WIDTH,
            SECONDARY_ALPHA,
        ),
        "cross_BA": (
            "BA fold",
            COLORS["same"],
            LINESTYLES["BA"],
            MARKERS["BA"],
            SECONDARY_LINE_WIDTH,
            SECONDARY_ALPHA,
        ),
        "cross_combined": (
            "Equal-weight AB/BA",
            COLORS["combined"],
            LINESTYLES["combined"],
            MARKERS["combined"],
            LINE_WIDTH,
            1.0,
        ),
    }

    for rep in reps:
        if rep != "cross_combined" and not SHOW_AB_BA_FOLDS_MAIN:
            continue

        g = (
            d[d["representation"].astype(str).eq(rep)]
            .sort_values("bin")
            .copy()
        )
        label, color, ls, marker, lw, alpha = specs[rep]

        x = pd.to_numeric(g["x_center_median"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(g["mean_dx_median"], errors="coerce").to_numpy(float)
        supported = g["plot_supported"].to_numpy(dtype=bool)

        # Preserve the x grid but create a visible break at unsupported bins.
        y = np.where(supported, y, np.nan)

        if rep == "cross_combined" and SHOW_PRIMARY_FORWARD_RIBBON:
            lo = pd.to_numeric(
                g["mean_dx_randomization_q025"], errors="coerce"
            ).to_numpy(float)
            hi = pd.to_numeric(
                g["mean_dx_randomization_q975"], errors="coerce"
            ).to_numpy(float)
            lo = np.where(supported, lo, np.nan)
            hi = np.where(supported, hi, np.nan)
            ax.fill_between(
                x, lo, hi,
                color=color,
                alpha=RIBBON_ALPHA,
                linewidth=0,
                zorder=1,
            )

        ax.plot(
            x,
            y,
            label=label,
            color=color,
            lw=lw,
            ls=ls,
            marker=marker,
            ms=MARKER_SIZE if marker else 0,
            alpha=alpha,
            zorder=3 if rep == "cross_combined" else 2,
        )

    ax.set_xlabel(r"Observed initial log-frequency, $x_0$")
    ax.set_ylabel(r"Mean observed displacement, $\langle \Delta x \rangle$")
    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 3B — matched same vs cross
# =============================================================================

def plot_figure3B_forward_matched(
    s7_envelope: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure3B_forward_matched"

    require_columns(
        s7_envelope,
        [
            "section",
            "representation",
            "bin",
            "n_configurations",
            "x_center_median",
            "mean_dx_median",
            "mean_dx_randomization_q025",
            "mean_dx_randomization_q975",
        ],
        "Step 7 curve envelope",
    )

    reps = ["same_combined", "cross_combined_matched"]
    d = s7_envelope[
        s7_envelope["section"].astype(str).eq("matched_same_vs_cross")
        & s7_envelope["representation"].astype(str).isin(reps)
    ].copy()

    if d.empty:
        raise ValueError("No matched same-vs-cross rows found in Step 7 envelope")

    supported_bins, n_total, threshold = step7_supported_bins(
        s7_envelope,
        section="matched_same_vs_cross",
        combined_representation="cross_combined_matched",
    )
    d = mark_step7_support(d, supported_bins, n_total, threshold)
    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)

    specs = {
        "same_combined": (
            "Same-measure",
            COLORS["same"],
            LINESTYLES["same"],
            LINE_WIDTH,
        ),
        "cross_combined_matched": (
            "Replicate-decoupled",
            COLORS["cross"],
            LINESTYLES["cross"],
            LINE_WIDTH,
        ),
    }

    for rep in reps:
        g = (
            d[d["representation"].astype(str).eq(rep)]
            .sort_values("bin")
            .copy()
        )
        label, color, ls, lw = specs[rep]

        x = pd.to_numeric(g["x_center_median"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(g["mean_dx_median"], errors="coerce").to_numpy(float)
        lo = pd.to_numeric(
            g["mean_dx_randomization_q025"], errors="coerce"
        ).to_numpy(float)
        hi = pd.to_numeric(
            g["mean_dx_randomization_q975"], errors="coerce"
        ).to_numpy(float)
        supported = g["plot_supported"].to_numpy(dtype=bool)

        y = np.where(supported, y, np.nan)
        lo = np.where(supported, lo, np.nan)
        hi = np.where(supported, hi, np.nan)

        if SHOW_MATCHED_FORWARD_RIBBONS:
            ax.fill_between(
                x, lo, hi,
                color=color,
                alpha=RIBBON_ALPHA,
                linewidth=0,
                zorder=1,
            )

        ax.plot(
            x,
            y,
            color=color,
            lw=lw,
            ls=ls,
            label=label,
            zorder=3,
        )

    ax.set_xlabel(r"Observed initial log-frequency, $x_0$")
    ax.set_ylabel(r"Mean observed displacement, $\langle \Delta x \rangle$")
    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 3C — cross covariance by abundance and pseudo-lag
# =============================================================================

def _filtered_step8_binned(
    s8_binned: pd.DataFrame,
    conditioning: str,
) -> pd.DataFrame:
    require_columns(
        s8_binned,
        [
            "conditioning",
            "dt",
            "bin_id",
            "x_center",
            "n_configurations",
            "n_configurations_valid",
            "cross_cov",
            "cross_cov_randomization_q025",
            "cross_cov_randomization_q975",
            "same_var_mean",
            "replicate_specific_excess",
        ],
        "Step 8 binned randomization envelope",
    )

    d = s8_binned[
        s8_binned["conditioning"].astype(str).eq(conditioning)
    ].copy()
    if d.empty:
        return d

    n_total = int(
        pd.to_numeric(
            d["n_configurations"], errors="coerce"
        ).max()
    )
    if n_total <= 0:
        raise ValueError("Could not infer Step-8 ensemble size")

    threshold = int(
        np.ceil(MIN_STEP8_CONFIG_FRACTION_PER_BIN * n_total)
    )
    d["n_configurations_valid"] = pd.to_numeric(
        d["n_configurations_valid"], errors="coerce"
    )
    d["plot_supported"] = (
        d["n_configurations_valid"] >= threshold
    )
    d["step8_total_configurations"] = n_total
    d["step8_min_configurations_for_plot"] = threshold
    return d


def plot_figure3C_cross_cov_by_lag(
    s8_binned: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure3C_cross_cov_by_lag"
    d = _filtered_step8_binned(s8_binned, "xmid_latent")

    if d.empty:
        raise ValueError("No supported xmid_latent Step-8 binned rows")

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)

    for dt in sorted(pd.to_numeric(d["dt"], errors="coerce").dropna().astype(int).unique()):
        g = d[pd.to_numeric(d["dt"], errors="coerce").eq(dt)].sort_values("x_center")
        key = f"lag{dt}"
        color = COLORS.get(key, "#333333")
        ls = LINESTYLES.get(key, "-")
        marker = MARKERS.get(key, None)

        x = pd.to_numeric(g["x_center"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(g["cross_cov"], errors="coerce").to_numpy(float)
        supported = g["plot_supported"].to_numpy(dtype=bool)
        y = np.where(supported, y, np.nan)

        if SHOW_CROSS_COV_LAG_RIBBONS:
            lo = pd.to_numeric(
                g["cross_cov_randomization_q025"], errors="coerce"
            ).to_numpy(float)
            hi = pd.to_numeric(
                g["cross_cov_randomization_q975"], errors="coerce"
            ).to_numpy(float)
            lo = np.where(supported, lo, np.nan)
            hi = np.where(supported, hi, np.nan)
            ax.fill_between(
                x, lo, hi,
                color=color,
                alpha=RIBBON_ALPHA * 0.7,
                linewidth=0,
            )

        ax.plot(
            x,
            y,
            color=color,
            lw=LINE_WIDTH,
            ls=ls,
            marker=marker,
            ms=MARKER_SIZE,
            label=f"Pseudo-lag {dt}",
        )

    ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
    ax.set_ylabel("Cross-replicate displacement\ncovariance")
    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 3D — fluctuation decomposition at dt=1
# =============================================================================

def plot_figure3D_fluctuation_decomposition_dt1(
    s8_binned: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure3D_fluctuation_decomposition_dt1"
    d = _filtered_step8_binned(s8_binned, "xmid_latent")
    d = d[pd.to_numeric(d["dt"], errors="coerce").eq(1)].copy()

    if d.empty:
        raise ValueError("No supported xmid_latent dt=1 rows for Step 8")

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)

    metric_specs = [
        (
            "cross_cov",
            "cross_cov_randomization_q025",
            "cross_cov_randomization_q975",
            "Cross-replicate covariance",
            COLORS["cross"],
            LINESTYLES["cross"],
            PRIMARY_LINE_WIDTH,
        ),
        (
            "same_var_mean",
            "same_var_mean_randomization_q025",
            "same_var_mean_randomization_q975",
            "Mean within-replicate variance",
            COLORS["same"],
            LINESTYLES["same"],
            LINE_WIDTH,
        ),
    ]

    x = pd.to_numeric(d["x_center"], errors="coerce").to_numpy(float)

    for metric, lo_col, hi_col, label, color, ls, lw in metric_specs:
        require_columns(d, [metric, lo_col, hi_col], panel)

        y = pd.to_numeric(d[metric], errors="coerce").to_numpy(float)
        supported = d["plot_supported"].to_numpy(dtype=bool)
        y = np.where(supported, y, np.nan)

        if SHOW_DECOMPOSITION_RIBBONS:
            lo = pd.to_numeric(d[lo_col], errors="coerce").to_numpy(float)
            hi = pd.to_numeric(d[hi_col], errors="coerce").to_numpy(float)
            lo = np.where(supported, lo, np.nan)
            hi = np.where(supported, hi, np.nan)
            ax.fill_between(
                x, lo, hi,
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

    ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
    ax.set_ylabel("Displacement variance / covariance")
    maybe_title(ax, panel)
    apply_limits(ax, panel)
    ax.legend(
        frameon=False,
        bbox_to_anchor=LEGEND_BBOX,
        loc=LEGEND_LOC,
    )
    save_figure(fig, outdir, panel)


# =============================================================================
# Figure 3E — pseudo-lag slope summary
# =============================================================================

def plot_figure3E_pseudolag_slopes(
    s8_slopes: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "Figure3E_pseudolag_slopes"

    require_columns(
        s8_slopes,
        [
            "metric",
            "n_configurations",
            "median",
            "randomization_q025",
            "randomization_q975",
        ],
        "Step 8 temporal slope summary",
    )

    selected = [
        (
            "cross_cov_slope_per_pseudolag",
            "Pooled cross-replicate covariance",
            COLORS["cross"],
            "o",
        ),
        (
            "equal_interval_cross_cov_slope_per_pseudolag",
            "Equal-interval cross-replicate covariance",
            COLORS["xstar"],
            "s",
        ),
    ]

    records = []
    for metric, label, color, marker in selected:
        g = s8_slopes[s8_slopes["metric"].astype(str).eq(metric)]
        if len(g) != 1:
            raise ValueError(
                f"Expected one row for {metric}, found {len(g)}"
            )
        row = g.iloc[0]
        records.append(
            {
                "metric": metric,
                "label": label,
                "color": color,
                "marker": marker,
                "n_configurations": int(row["n_configurations"]),
                "median": float(row["median"]),
                "q025": float(row["randomization_q025"]),
                "q975": float(row["randomization_q975"]),
            }
        )

    source = pd.DataFrame(records)
    write_source_data(source, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax, orientation="vertical")

    y_positions = np.arange(len(records))[::-1]

    for y, record in zip(y_positions, records):
        lower = record["median"] - record["q025"]
        upper = record["q975"] - record["median"]

        ax.errorbar(
            record["median"],
            y,
            xerr=np.array([[lower], [upper]]),
            fmt=record["marker"],
            color=record["color"],
            ms=POINT_SIZE,
            lw=ERRORBAR_LINE_WIDTH,
            capsize=CAP_SIZE,
            zorder=3,
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels([r["label"] for r in records])
    ax.set_xlabel("Slope per pseudo-lag")
    ax.set_ylabel("")
    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S4A — conditioning geometry
# =============================================================================

def plot_figureS4A_conditioning_geometry(
    s7_envelope: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "FigureS4A_conditioning_geometry"

    require_columns(
        s7_envelope,
        [
            "section",
            "representation",
            "bin",
            "n_configurations",
            "x_center_median",
            "mean_dx_median",
        ],
        "Step 7 curve envelope",
    )

    reps = ["x0_latent", "xmid_latent", "xstar_latent"]
    d = s7_envelope[
        s7_envelope["section"].astype(str).eq(
            "conditioning_geometry_unrestricted"
        )
        & s7_envelope["representation"].astype(str).isin(reps)
    ].copy()

    if d.empty:
        raise ValueError(
            "No conditioning_geometry_unrestricted rows found in Step 7 envelope"
        )

    # The current conditioning block has >90% configuration support in all bins,
    # but retain the same transparent support flag in the exported source data.
    n_total = int(
        pd.to_numeric(d["n_configurations"], errors="coerce").max()
    )
    threshold = int(np.ceil(MIN_STEP7_CONFIG_FRACTION_PER_BIN * n_total))
    d["plot_supported"] = (
        pd.to_numeric(d["n_configurations"], errors="coerce") >= threshold
    )
    d["step7_total_configurations"] = n_total
    d["step7_min_configurations_for_plot"] = threshold

    write_source_data(d, source_dir, panel)

    specs = {
        "x0_latent": (
            r"$x_0$",
            COLORS["x0"],
            LINESTYLES["x0"],
        ),
        "xmid_latent": (
            r"$x_{\mathrm{mid}}$",
            COLORS["xmid"],
            LINESTYLES["xmid"],
        ),
        "xstar_latent": (
            r"$x_{\star}$",
            COLORS["xstar"],
            LINESTYLES["xstar"],
        ),
    }

    fig, ax = make_figure(panel)
    zero_line(ax)

    for rep in reps:
        g = (
            d[d["representation"].astype(str).eq(rep)]
            .sort_values("bin")
            .copy()
        )
        label, color, ls = specs[rep]

        x = pd.to_numeric(g["x_center_median"], errors="coerce").to_numpy(float)
        y = pd.to_numeric(g["mean_dx_median"], errors="coerce").to_numpy(float)
        supported = g["plot_supported"].to_numpy(dtype=bool)
        y = np.where(supported, y, np.nan)

        ax.plot(
            x,
            y,
            color=color,
            lw=LINE_WIDTH,
            ls=ls,
            label=label,
        )

    ax.set_xlabel("Conditioning coordinate")
    ax.set_ylabel(
        r"Mean latent displacement, $\langle \Delta x_{\mathrm{latent}} \rangle$"
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
# Supplementary Figure S4B — common4 support
# =============================================================================

def plot_figureS4B_common4_support(
    s8_support: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "FigureS4B_common4_support"

    require_columns(
        s8_support,
        ["configuration_id", "dt", "fraction_common4"],
        "Step 8 support table",
    )

    rows = []
    for dt, g in s8_support.groupby("dt", sort=True):
        x = pd.to_numeric(
            g["fraction_common4"], errors="coerce"
        ).dropna().to_numpy(float)
        if x.size == 0:
            continue
        rows.append(
            {
                "dt": int(dt),
                "n_configurations": int(len(x)),
                "median": float(np.median(x)),
                "q025": float(np.quantile(x, 0.025)),
                "q975": float(np.quantile(x, 0.975)),
            }
        )

    source = pd.DataFrame(rows)
    write_source_data(source, source_dir, panel)

    fig, ax = make_figure(panel)

    x = source["dt"].to_numpy(float)
    y = source["median"].to_numpy(float)
    lower = y - source["q025"].to_numpy(float)
    upper = source["q975"].to_numpy(float) - y

    ax.errorbar(
        x,
        y,
        yerr=np.vstack([lower, upper]),
        fmt="o-",
        color=COLORS["support"],
        lw=LINE_WIDTH,
        ms=MARKER_SIZE,
        capsize=CAP_SIZE,
    )

    ax.set_xlabel("Pseudo-lag")
    ax.set_ylabel("Fraction of transitions in common4")
    ax.set_xticks(sorted(source["dt"].astype(int).unique()))
    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Supplementary Figure S4C — pooled cross covariance by pseudo-lag
# =============================================================================

def plot_figureS4C_pooled_cross_cov_by_lag(
    s8_dt_summary: pd.DataFrame,
    outdir: Path,
    source_dir: Path,
) -> None:
    panel = "FigureS4C_pooled_cross_cov_by_lag"

    require_columns(
        s8_dt_summary,
        [
            "dt",
            "metric",
            "n_configurations",
            "median",
            "randomization_q025",
            "randomization_q975",
        ],
        "Step 8 dt summary",
    )

    d = s8_dt_summary[
        s8_dt_summary["metric"].astype(str).eq("cross_cov")
    ].copy().sort_values("dt")

    if d.empty:
        raise ValueError("No cross_cov rows in Step 8 dt summary")

    write_source_data(d, source_dir, panel)

    fig, ax = make_figure(panel)
    zero_line(ax)

    x = pd.to_numeric(d["dt"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(d["median"], errors="coerce").to_numpy(float)
    lo = pd.to_numeric(
        d["randomization_q025"], errors="coerce"
    ).to_numpy(float)
    hi = pd.to_numeric(
        d["randomization_q975"], errors="coerce"
    ).to_numpy(float)

    lower = y - lo
    upper = hi - y

    ax.errorbar(
        x,
        y,
        yerr=np.vstack([lower, upper]),
        fmt="o-",
        color=COLORS["cross"],
        lw=LINE_WIDTH,
        ms=MARKER_SIZE,
        capsize=CAP_SIZE,
    )

    ax.set_xlabel("Pseudo-lag")
    ax.set_ylabel("Cross-replicate displacement\ncovariance")
    ax.set_xticks(sorted(pd.to_numeric(d["dt"], errors="coerce").astype(int).unique()))
    maybe_title(ax, panel)
    apply_limits(ax, panel)
    save_figure(fig, outdir, panel)


# =============================================================================
# Main
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Create single-panel publication plots for final pseudo Step 7/8 "
            "technical-null analyses."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--step7-dir", required=True, type=Path)
    p.add_argument("--step8-dir", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    configure_matplotlib()

    step7_dir = args.step7_dir.expanduser().resolve(strict=True)
    step8_dir = args.step8_dir.expanduser().resolve(strict=True)
    outdir = args.outdir.expanduser().resolve()

    main_dir = outdir / "main_figure"
    supplementary_dir = outdir / "supplementary_figure"
    source_dir = outdir / "figure_source_data"

    inputs = load_inputs(step7_dir, step8_dir)

    print("[INFO] Step 7 contract: compact-native equal-weight AB/BA")
    print("[INFO] Step 7 minimum bin support fraction:", MIN_STEP7_CONFIG_FRACTION_PER_BIN)
    print("[INFO] Step 8 contract: fast pair-transition-bank")
    print("[INFO] Step 8 minimum bin support fraction:", MIN_STEP8_CONFIG_FRACTION_PER_BIN)

    plot_figure3A_forward_primary(
        inputs["s7_envelope"],
        main_dir,
        source_dir,
    )
    plot_figure3B_forward_matched(
        inputs["s7_envelope"],
        main_dir,
        source_dir,
    )
    plot_figure3C_cross_cov_by_lag(
        inputs["s8_binned"],
        main_dir,
        source_dir,
    )
    plot_figure3D_fluctuation_decomposition_dt1(
        inputs["s8_binned"],
        main_dir,
        source_dir,
    )
    plot_figure3E_pseudolag_slopes(
        inputs["s8_slopes"],
        main_dir,
        source_dir,
    )

    plot_figureS4A_conditioning_geometry(
        inputs["s7_envelope"],
        supplementary_dir,
        source_dir,
    )
    plot_figureS4B_common4_support(
        inputs["s8_support"],
        supplementary_dir,
        source_dir,
    )
    plot_figureS4C_pooled_cross_cov_by_lag(
        inputs["s8_dt_summary"],
        supplementary_dir,
        source_dir,
    )

    # Explicit plot-support audit for Step 8.
    s8_audit = _filtered_step8_binned(
        inputs["s8_binned"], "xmid_latent"
    )[
        [
            "dt", "bin_id", "x_center",
            "n_configurations", "n_configurations_valid",
            "plot_supported", "step8_total_configurations",
            "step8_min_configurations_for_plot",
        ]
    ].copy()
    s8_audit.to_csv(
        outdir / "00_step8_plot_support_audit.csv",
        index=False,
    )

    print("\n[DONE] pseudo technical-null publication plots")
    print("Main panels          :", main_dir)
    print("Supplementary panels :", supplementary_dir)
    print("Figure source data   :", source_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
