#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ClonoDynamics: plot the integrated Figure 8 from finalized Steps 14–16.

Purpose
-------
Create the finalized main Figure 8 panels and the three subsequent supplementary
figures without rerunning any scientific analysis.

Main Figure 8:
A. Step-13 reference and common-start/common-end matched-subject slopes.
B. Step-15 PRIMARY, T0PLUS and TT observed forward profiles at one week.
C. Step-16 one-week common-core cross-covariance versus operational threshold.
D. Step-16 common-core temporal slopes versus operational threshold.

Supplementary Figure 9 (Step 14):
A. Core-level calendar-position q values.
B. Abundance-resolved calendar-position q values for cross-replicate covariance.
C. Common-start abundance-resolved temporal slopes.
D. Common-end abundance-resolved temporal slopes.

Supplementary Figure 10 (Step 15, reference alpha):
A. Operational-class composition across temporal lags.
B. Boundary-crossing imbalance across abundance and lag.
C. Forward-profile differences between operational domains.
D. One-week cross-replicate covariance differences versus PRIMARY.
E. One-week within-replicate displacement variance by domain.
F. One-week replicate-specific excess by domain.
G. Common-core temporal covariance profiles by domain.

Supplementary Figure 11 (Step 16, cross-threshold sensitivity):
A. Retention of the reference TT population across operational thresholds.
B. T0PLUS forward profiles on the shared cross-alpha forward support.
C. TT forward profiles on the same shared support.
D. TF forward profiles on the same shared support.
E. Forward-slope differences relative to the reference alpha.
F. NON_FF one-week cross-replicate covariance profiles across thresholds.
G. TT one-week covariance differences relative to the reference alpha.
H. TT common-core temporal covariance profiles across thresholds.

Inputs (analysis directories, not figure directories)
---------------------------------------------------
--step14-dir:
    00_run_config.json
    02_core_bin_definitions.csv
    09_anchored_core_slope_bootstrap.csv
--step15-dir:
    00_run_config.json
    03_forward_domain_by_bin.csv
--step16-dir:
    00_run_config.json
    03_tt_retention_by_alpha.csv
    04_forward_domain_cross_alpha_by_bin.csv
    05b_forward_cross_alpha_common_support.csv
    07_forward_slope_difference_vs_reference_alpha.csv
    08_fluctuation_domain_cross_alpha_by_bin_dt.csv
    09_fluctuation_difference_vs_reference_alpha.csv
    11_cross_alpha_temporal_common_core.csv
    12_temporal_profiles_cross_alpha.csv
    13_temporal_slopes_cross_alpha.csv
--step13-dir [optional, recommended]:
    00_run_config.json
    05_temporal_slope_summary.csv

The schemas are those of the replicate-resolved cross_cov branch. Legacy
var_dx/xstar/TT-primary tables are not substituted. Input files are read-only.
The optional Step-13 input independently checks the numerical reference,
empirical core and subject universe inherited by Step 14.

Scientific boundaries
---------------------
No refitting, bootstrap, new covariance, smoothing, interpolation or rebinning
is performed. Unsupported forward points remain gaps. No numerical result is
hard-coded. Panel A retains the three prespecified comparisons regardless of
sign or whether an interval crosses zero. Panel B displays the PRIMARY
bootstrap band only; TF is omitted from this main panel, not from the analysis.
Panels C/D retain the Step-16 common-core pooled-covariance sensitivity estimator;
it is not replaced by the Step-13 equal-bin/equal-subject reference.

IMPORTANT: the inspected Step-16 profile table exports POINT ESTIMATES ONLY.
Panel C therefore has NO confidence bars. The program does not derive amplitude
intervals from slope intervals or from abundance-bin intervals. Panel D uses
its existing slope intervals. A PRIMARY horizontal reference in C/D is drawn
only after checking that the supplied PRIMARY values are threshold-invariant.

Outputs
-------
Four individual main Figure 8 panels; Supplementary Figures 9, 10 and 11 as both
assembled multipanel figures and separate panels; exact selected-source CSVs for
every displayed supplementary panel; input hashes, plotting checks, configuration,
manifest and an English figure note. No scientific statistic or Excel workbook is
generated. --validate-only writes source selections/checks but no figures. Existing
output folders are refused unless --overwrite is explicit and the previous folder
was created by this script; unrelated files are kept.

Customization
-------------
Main-panel aesthetics are in FONT, FIGURE_SIZES, LINE_WIDTHS, SERIES,
AXIS_CONFIG, LEGEND_CONFIG and BOX_CONFIG. Supplementary axis and legend
placement are controlled independently by SUPP_AXIS_CONFIG and
SUPP_LEGEND_CONFIG; layouts/styles use SUPP_FIGURE_SIZES, SUPP9_METRICS,
SUPP10_CLASSES, SUPP10_CONTRASTS and SUPP11_STEMS.
BOX_CONFIG applies to main Figure 8 only. Colours are taken from the active
Matplotlib colour cycle with stable series/domain mappings.

Dependencies: Python >=3.9, NumPy, pandas, Matplotlib. No original plotter imports,
LaTeX installation, large transition files, spreadsheets or ground truth needed.

Example (run from the project root)
----------------------------------
python code/05_plotting/04_controls/plot_steps14_16_robustness_summary.py \\
  --step14-dir dataset_longitudinal_results/14-interval_position_structure \\
  --step15-dir dataset_longitudinal_results/15-detectability_boundary_sensitivity \\
  --step16-dir dataset_longitudinal_results/16-observation_threshold_robustness \\
  --step13-dir dataset_longitudinal_results/13-temporal_fluctuation_scaling \\
  --outdir figures/figure8_robustness
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import shlex
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator, ScalarFormatter

VERSION = "1.4.0-figure8-supp9-supp10-supp11-axis-legend-config-2026-09-22"
OWNER = "clonodynamics_figure8_robustness_plotter"
RTOL, ATOL = 1e-8, 1e-11

# ======================== USER-EDITABLE AESTHETICS ===========================
FONT = {"family": "Arial", "size": 8, "axis": 8, "tick": 7,
        "legend": 7, "annotation": 7, "panel_label": 10}
FIGURE_SIZES = {p: (3.40, 2.6) for p in "ABCD"}  # inches; fixed canvas
LINE_WIDTHS = {"curve": 1.15, "primary": 1.45, "ci": 0.9,
               "reference": 0.75, "axis": 0.8}
SERIES = {
    "PRIMARY": {"cycle_index": 0, "marker": "o", "linestyle": "-"},
    "T0PLUS": {"cycle_index": 1, "marker": "s", "linestyle": "-"},
    "TT": {"cycle_index": 2, "marker": "^", "linestyle": "-"},
    "NON_FF": {"cycle_index": 3, "marker": "D", "linestyle": "-"},
    "common_start": {"cycle_index": 1, "marker": "s", "linestyle": "None"},
    "common_end": {"cycle_index": 2, "marker": "^", "linestyle": "None"},
    "TF": {"cycle_index": 4, "marker": "D", "linestyle": "-"},
}
MARKER_SIZE = 4.0
CI_ALPHA = 0.15
SHOW_PANEL_LABELS = True
SHOW_TITLES = False
# Axis configuration for the four main Figure 8 panels.
#
# xlim / ylim:
#     None                       -> automatic limits
#     (minimum, maximum)         -> fixed displayed interval
# xtick_step / ytick_step:
#     None                       -> automatic tick spacing
#     positive number            -> fixed regular tick spacing
# xticks / yticks:
#     None                       -> use automatic/step-based ticks
#     [v1, v2, ...]              -> explicit tick positions (takes precedence over *_step)
#
# Examples:
#     "A": {"xlim": (-0.10, 0.03), "ylim": None,
#           "xtick_step": 0.025, "ytick_step": None,
#           "xticks": None, "yticks": None}
#
#     "C": {"xlim": (0.005, 0.105), "ylim": (0.0, 0.60),
#           "xtick_step": None, "ytick_step": 0.10,
#           "xticks": [0.01, 0.025, 0.05, 0.10], "yticks": None}
AXIS_CONFIG = {
    "A": {"xlim": None, "ylim": None, "xtick_step": None, "ytick_step": None,
          "xticks": None, "yticks": None},
    "B": {"xlim": None, "ylim": None, "xtick_step": None, "ytick_step": None,
          "xticks": None, "yticks": None},
    "C": {"xlim": None, "ylim": None, "xtick_step": None, "ytick_step": None,
          "xticks": None, "yticks": None},
    "D": {"xlim": None, "ylim": None, "xtick_step": None, "ytick_step": None,
          "xticks": None, "yticks": None},
}

# Independent axis configuration for every Supplementary Figure 9–11 panel.
# Keys are S9A...S9D, S10A...S10G and S11A...S11H.
# The same xlim/ylim/tick rules described above apply here.
def _axis_defaults():
    return {"xlim": None, "ylim": None,
            "xtick_step": None, "ytick_step": None,
            "xticks": None, "yticks": None}

SUPP_AXIS_CONFIG = {
    **{f"S9{p}": _axis_defaults() for p in "ABCD"},
    **{f"S10{p}": _axis_defaults() for p in "ABCDEFG"},
    **{f"S11{p}": _axis_defaults() for p in "ABCDEFGH"},
}
SUPP_AXIS_CONFIG["S9A"] = {
    "xlim": (0.5, 3.5),
    "ylim": (0.0, 0.8),
    "xtick_step": 0.5,
    "ytick_step": 0.2,
    "xticks": None,
    "yticks": None,
}
SUPP_AXIS_CONFIG["S9B"] = {
    "xlim": None,
    "ylim": (0.0, 1.3),
    "xtick_step": 0.5,
    "ytick_step": 0.2,
    "xticks": None,
    "yticks": None,
}
SUPP_AXIS_CONFIG["S10A"] = {
    "xlim": None,
    "ylim": (0.0, 1.2),
    "xtick_step": 0.5,
    "ytick_step": 0.2,
    "xticks": None,
    "yticks": None,
}

SUPP_AXIS_CONFIG["S10B"] = {
    "xlim": None,
    "ylim": (-1.5, 1.1),
    "xtick_step": 0.5,
    "ytick_step": 0.2,
    "xticks": None,
    "yticks": None,
}
SUPP_AXIS_CONFIG["S11H"] = {
    "xlim": None,
    "ylim": (0.0, 0.5),
    "xtick_step": 0.5,
    "ytick_step": 0.2,
    "xticks": None,
    "yticks": None,
}
LEGEND_CONFIG = {
    "A": {"enabled": False, "loc": "best", "bbox_to_anchor": None, "ncol": 1},
    "B": {"enabled": True, "loc": "upper right", "bbox_to_anchor": None, "ncol": 1},
    "C": {"enabled": True, "loc": "best", "bbox_to_anchor": None, "ncol": 1},
    "D": {"enabled": True, "loc": "best", "bbox_to_anchor": None, "ncol": 1},
}

# Independent legend configuration for every Supplementary Figure 9–11 panel.
# `loc` accepts standard Matplotlib positions. For finer placement, set
# `bbox_to_anchor=(x, y)` in axes coordinates. Panels without a legend are
# disabled by default. These settings apply to both individual and composite plots.
SUPP_LEGEND_CONFIG = {
    "S9A":  {"enabled": True,  "loc": "upper left",  "bbox_to_anchor": None, "ncol": 1},
    "S9B":  {"enabled": True,  "loc": "upper left",  "bbox_to_anchor": None, "ncol": 1},
    "S9C":  {"enabled": False, "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
    "S9D":  {"enabled": False, "loc": "best",        "bbox_to_anchor": None, "ncol": 1},

    "S10A": {"enabled": True,  "loc": "upper right", "bbox_to_anchor": None, "ncol": 1},
    "S10B": {"enabled": True,  "loc": "lower left",  "bbox_to_anchor": None, "ncol": 1},
    "S10C": {"enabled": True,  "loc": "lower right", "bbox_to_anchor": None, "ncol": 1},
    "S10D": {"enabled": True,  "loc": "upper right", "bbox_to_anchor": None, "ncol": 1},
    "S10E": {"enabled": True,  "loc": "upper right", "bbox_to_anchor": None, "ncol": 1},
    "S10F": {"enabled": True,  "loc": "upper right", "bbox_to_anchor": None, "ncol": 1},
    "S10G": {"enabled": True,  "loc": "best",        "bbox_to_anchor": None, "ncol": 1},

    "S11A": {"enabled": False, "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
    "S11B": {"enabled": True,  "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
    "S11C": {"enabled": True,  "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
    "S11D": {"enabled": True,  "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
    "S11E": {"enabled": True,  "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
    "S11F": {"enabled": True,  "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
    "S11G": {"enabled": True,  "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
    "S11H": {"enabled": True,  "loc": "best",        "bbox_to_anchor": None, "ncol": 1},
}
# Change x/y here to move each annotation; use text="..." to override its text.
BOX_CONFIG= {
    "A": {"enabled": False, "x": 0.98, "y": 0.03, "ha": "right", "va": "bottom", "text": None},
    "B": {"enabled": True, "x": 0.98, "y": 0.03, "ha": "right", "va": "bottom", "text": None},
    "C": {"enabled": True, "x": 0.03, "y": 0.03, "ha": "left", "va": "bottom", "text": None},
    "D": {"enabled": True, "x": 0.98, "y": 0.98, "ha": "right", "va": "top", "text": None},
}

STEMS = {
    "A": "Figure8A_interval_controls",
    "B": "Figure8B_forward_domains",
    "C": "Figure8C_covariance_amplitude_by_threshold",
    "D": "Figure8D_temporal_slopes_by_threshold",
}
TITLES = {"A": "Interval-composition controls", "B": "Observation-domain selection",
          "C": "Fluctuation amplitude", "D": "Temporal accumulation"}

SUPP_FIGURE_SIZES = {
    "S9_panel": (3.40, 2.60),
    "S9_composite": (7.20, 5.80),
    "S10_panel": (3.40, 2.60),
    "S10_composite": (7.20, 10.40),
    "S11_panel": (3.40, 2.60),
    "S11_composite": (7.20, 10.40),
}
SUPP_FIGURE_SIZES["10G"] = {
    "S10_panel": (3.40, 2.60),}

SUPP9_METRICS = {
    "cross_cov": {"label": "Cross-replicate covariance", "cycle_index": 0, "marker": "o"},
    "same_var_mean": {"label": "Within-replicate variance", "cycle_index": 1, "marker": "s"},
    "replicate_specific_excess": {"label": "Replicate-specific excess", "cycle_index": 2, "marker": "^"},
}
SUPP10_CLASSES = {
    "TT": {"cycle_index": 2, "marker": "^"},
    "TF": {"cycle_index": 1, "marker": "D"},
    "FT": {"cycle_index": 4, "marker": "s"},
    "FF": {"cycle_index": 7, "marker": "o"},
}
SUPP10_CONTRASTS = {
    "T0PLUS-PRIMARY": {"cycle_index": 1, "marker": "s"},
    "TT-PRIMARY": {"cycle_index": 2, "marker": "^"},
    "T0PLUS-TT": {"cycle_index": 4, "marker": "D"},
}
SUPP9_STEMS = {
    "A": "SupplementaryFigure9A_core_position_qvalues",
    "B": "SupplementaryFigure9B_binwise_position_qvalues",
    "C": "SupplementaryFigure9C_common_start_binwise_slopes",
    "D": "SupplementaryFigure9D_common_end_binwise_slopes",
}
SUPP10_STEMS = {
    "A": "SupplementaryFigure10A_class_composition",
    "B": "SupplementaryFigure10B_crossing_imbalance",
    "C": "SupplementaryFigure10C_forward_differences",
    "D": "SupplementaryFigure10D_crosscov_difference_dt1",
    "E": "SupplementaryFigure10E_samevar_domains_dt1",
    "F": "SupplementaryFigure10F_excess_domains_dt1",
    "G": "SupplementaryFigure10G_temporal_profiles",
}
SUPP11_STEMS = {
    "A": "SupplementaryFigure11A_TT_retention",
    "B": "SupplementaryFigure11B_T0PLUS_forward_curves",
    "C": "SupplementaryFigure11C_TT_forward_curves",
    "D": "SupplementaryFigure11D_TF_forward_curves",
    "E": "SupplementaryFigure11E_forward_slope_differences",
    "F": "SupplementaryFigure11F_NONFF_crosscov_dt1",
    "G": "SupplementaryFigure11G_TT_crosscov_difference_dt1",
    "H": "SupplementaryFigure11H_TT_temporal_profiles",
}
SHOW_SUPP_INDIVIDUAL_PANEL_LABELS = True
SHOW_SUPP_COMPOSITE_PANEL_LABELS = False


class InputError(ValueError):
    """An absent, incompatible or internally inconsistent source input."""


class Audit:
    def __init__(self) -> None:
        self.rows: List[Dict[str, Any]] = []
        self.inputs: List[Dict[str, Any]] = []

    def check(self, name: str, condition: Any, detail: str = "") -> None:
        ok = bool(condition)
        self.rows.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": detail})
        if not ok:
            raise InputError(f"{name}: {detail}")

    def note(self, name: str, detail: str) -> None:
        self.rows.append({"check": name, "status": "NOTE", "detail": detail})

    def register(self, path: Path, rows: Optional[int] = None) -> None:
        self.inputs.append({"path": str(path), "filename": path.name,
                            "sha256": digest(path), "size_bytes": path.stat().st_size,
                            "data_rows": rows})


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_csv(folder: Path, name: str, audit: Audit) -> pd.DataFrame:
    path = folder / name
    if not path.is_file():
        raise InputError(f"Required file missing: {path}\nSupply the finalized analysis folder, not its figure folder.")
    frame = pd.read_csv(path)
    audit.check(f"Nonempty {name}", not frame.empty, str(path))
    audit.register(path, len(frame))
    # One-based row number in the ORIGINAL CSV, including its header.
    if "source_csv_row" in frame:
        raise InputError(f"Reserved column source_csv_row is already present in {path}")
    frame.insert(0, "source_csv_row", np.arange(2, len(frame) + 2))
    return frame


def load_json(folder: Path, name: str, audit: Audit) -> Dict[str, Any]:
    path = folder / name
    if not path.is_file():
        raise InputError(f"Required metadata missing: {path}")
    with path.open(encoding="utf-8") as f:
        obj = json.load(f)
    audit.check(f"JSON object {folder.name}/{name}", isinstance(obj, dict))
    audit.register(path)
    return obj


def require_columns(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise InputError(f"{label}: missing columns {missing}. Available: {list(frame.columns)}. "
                         "Legacy var_dx/xstar tables are not accepted as cross_cov tables.")


def get(obj: Dict[str, Any], dotted: str) -> Any:
    value: Any = obj
    for key in dotted.split("."):
        if not isinstance(value, dict) or key not in value:
            raise InputError(f"Required configuration field is absent: {dotted}")
        value = value[key]
    return value


def integers(values: Any, label: str) -> List[int]:
    if isinstance(values, str):
        values = [s for s in re.split(r"[;,\s]+", values.strip()) if s]
    try:
        a = np.asarray(values, dtype=float)
    except (ValueError, TypeError) as exc:
        raise InputError(f"{label}: invalid integer list") from exc
    if a.ndim != 1 or len(a) == 0 or not np.isfinite(a).all() or not np.equal(a, np.floor(a)).all():
        raise InputError(f"{label}: expected a nonempty list of finite integers")
    if len(set(a.tolist())) != len(a):
        raise InputError(f"{label}: duplicate identifiers")
    return sorted(a.astype(int).tolist())


def numeric(frame: pd.DataFrame, cols: Sequence[str], label: str,
            allow_nan: bool = False, integer: bool = False) -> None:
    require_columns(frame, cols, label)
    for col in cols:
        try:
            frame[col] = pd.to_numeric(frame[col], errors="raise")
        except (ValueError, TypeError) as exc:
            raise InputError(f"{label}.{col}: not numeric") from exc
        a = frame[col].to_numpy(float)
        if np.isinf(a).any() or (not allow_nan and not np.isfinite(a).all()):
            raise InputError(f"{label}.{col}: unexpected non-finite value")
        if integer:
            if not np.equal(a, np.floor(a)).all():
                raise InputError(f"{label}.{col}: integer values required")
            frame[col] = a.astype(int)


def same(a: Any, b: Any) -> bool:
    aa, bb = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return aa.shape == bb.shape and bool(np.allclose(aa, bb, rtol=RTOL, atol=ATOL, equal_nan=True))


def unique(frame: pd.DataFrame, keys: Sequence[str], label: str, audit: Audit) -> None:
    require_columns(frame, keys, label)
    audit.check(f"Unique keys: {label}", not frame.duplicated(list(keys)).any(), str(list(keys)))


def validate_ci(frame: pd.DataFrame, lo: str, hi: str, label: str, allow_nan: bool = False) -> None:
    numeric(frame, [lo, hi], label, allow_nan=allow_nan)
    left, right = frame[lo].to_numpy(float), frame[hi].to_numpy(float)
    if not np.array_equal(np.isnan(left), np.isnan(right)):
        raise InputError(f"{label}: only one endpoint of a confidence interval is missing")
    valid = np.isfinite(left) & np.isfinite(right)
    if np.any(left[valid] > right[valid]):
        raise InputError(f"{label}: lower confidence limit exceeds upper limit")
    # A percentile interval need not contain the point estimate: do not clamp it.


def constant(frame: pd.DataFrame, col: str, expected: Any, label: str, audit: Audit) -> None:
    numeric(frame, [col], label)
    audit.check(label, same(frame[col].to_numpy(float), np.full(len(frame), expected)),
                f"Expected {col}={expected}")


def source(frame: pd.DataFrame, path: str) -> pd.DataFrame:
    out = frame.copy()
    out.insert(0, "source_table", path)
    return out


def load_panels(step14: Path, step15: Path, step16: Path,
                step13: Optional[Path], audit: Audit) -> Tuple[Dict[str, pd.DataFrame], Dict[str, Any]]:
    c14 = load_json(step14, "00_run_config.json", audit)
    c15 = load_json(step15, "00_run_config.json", audit)
    c16 = load_json(step16, "00_run_config.json", audit)
    audit.check("Step14 metric", get(c14, "primary_metric") == "cross_cov")
    audit.check("Step14 observation support", get(c14, "primary_support") == "common4")
    audit.check("Step14 conditioning", get(c14, "primary_conditioning") == "xmid_latent")
    audit.check("Step14 has no operational filter", get(c14, "operational_class_filter") is None)
    for label, cfg in [("Step15", c15), ("Step16", c16)]:
        audit.check(f"{label} common4 support", get(cfg, "fluctuations.base_support") == "common4")
        audit.check(f"{label} xmid conditioning", get(cfg, "fluctuations.conditioning") == "xmid_latent")
        audit.check(f"{label} no pseudo subtraction", get(cfg, "step12_used_computationally") is False)
        audit.check(f"{label} one-week forward", get(cfg, "forward.dt") == 1)
    alpha = float(get(c15, "reference_operational_alpha"))
    audit.check("Same reference threshold", same(alpha, get(c16, "reference_alpha")))
    alphas = sorted(float(v) for v in get(c16, "alphas"))
    audit.check("Threshold range", len(alphas) >= 2 and len(set(alphas)) == len(alphas)
                and all(math.isfinite(v) and 0 < v < 1 for v in alphas))
    audit.check("Reference threshold included", any(same(alpha, a) for a in alphas))

    bins14 = integers(get(c14, "core_bins"), "Step14 core")
    subj14 = integers(get(c14, "complete_case_subjects"), "Step14 subjects")
    dts = integers(get(c14, "dt_values"), "Step14 lags")
    audit.check("Five requested lags", dts == [1, 2, 3, 4, 5])
    audit.check("Step15 frozen core", bins14 == integers(get(c15, "temporal_sensitivity.frozen_step13_core_bins"), "Step15 core"))
    audit.check("Step16 inherited core", bins14 == integers(get(c16, "temporal.frozen_step13_14_core_bins"), "Step16 inherited core"))
    audit.check("Step15 subjects", subj14 == integers(get(c15, "temporal_sensitivity.subject_universe"), "Step15 subjects"))
    audit.check("Step16 subjects", subj14 == integers(get(c16, "temporal.subject_universe"), "Step16 subjects"))
    audit.check("Step15 lag set", dts == integers(get(c15, "temporal_sensitivity.dt_values"), "Step15 lags"))
    audit.check("Step16 lag set", dts == integers(get(c16, "temporal.dt_values"), "Step16 lags"))
    common = integers(get(c16, "temporal.cross_alpha_common_core_bins"), "Cross-alpha core")
    audit.check("Common core is contiguous subset", set(common).issubset(bins14)
                and common == list(range(common[0], common[-1] + 1)))
    audit.check("Step16 sensitivity estimator",
                get(c16, "temporal.estimator") == "common_core_pooled_covariance_equal_bin")

    core14 = load_csv(step14, "02_core_bin_definitions.csv", audit)
    require_columns(core14, ["bin", "x_left", "x_right", "complete_case_subjects"], "Step14 core")
    numeric(core14, ["bin"], "Step14 core", integer=True)
    numeric(core14, ["x_left", "x_right"], "Step14 core")
    unique(core14, ["bin"], "Step14 core", audit)
    audit.check("Step14 core table matches config", sorted(core14["bin"].tolist()) == bins14)
    for s in core14["complete_case_subjects"].astype(str).unique():
        audit.check("Step14 core-table subjects", integers(s, "core subjects") == subj14)
    core16 = load_csv(step16, "11_cross_alpha_temporal_common_core.csv", audit)
    require_columns(core16, ["bin", "selected_cross_alpha_common_core"], "Step16 common core")
    numeric(core16, ["bin"], "Step16 common core", integer=True)
    flag = core16["selected_cross_alpha_common_core"].astype(str).str.lower().str.strip()
    audit.check("Recognized common-core flags", flag.isin(["true", "false", "1", "0", "1.0", "0.0"]).all())
    selected = core16.loc[flag.isin(["true", "1", "1.0"])].copy()
    audit.check("Step16 selected core matches config", sorted(selected["bin"].tolist()) == common)
    for coord in ["x_left", "x_right", "x_center"]:
        if coord in core16 and coord in core14:
            joined = selected[["bin", coord]].merge(core14[["bin", coord]], on="bin", validate="one_to_one")
            audit.check(f"Common-grid coordinates {coord}", same(joined[coord + "_x"], joined[coord + "_y"]))

    # A: stored reference plus the two prespecified matched analyses.
    anch = load_csv(step14, "09_anchored_core_slope_bootstrap.csv", audit)
    cols = ["metric", "anchor_type", "subject_matching", "n_subjects_analysis",
            "observed_slope_per_week", "bootstrap_q025", "bootstrap_q975",
            "step13_observed_slope_per_week", "step13_bootstrap_q025", "step13_bootstrap_q975"]
    require_columns(anch, cols, "Step14 anchors")
    a = anch.loc[anch["metric"].eq("cross_cov")].copy()
    audit.check("Cross-covariance anchor rows present", not a.empty)
    refs = ["step13_observed_slope_per_week", "step13_bootstrap_q025", "step13_bootstrap_q975"]
    numeric(a, refs, "Step14 reference")
    reference = {k: float(a.iloc[0][k]) for k in refs}
    for k in refs:
        constant(a, k, reference[k], f"Repeated stored reference {k}", audit)
        if "step13_reference" in c14:
            audit.check(f"Step14 config reference {k}", same(reference[k], get(c14, "step13_reference." + k)))
    audit.check("Reference interval order", reference[refs[1]] <= reference[refs[2]])
    rows = [{"source_table": "Step14/09_anchored_core_slope_bootstrap.csv",
             "source_csv_row": int(a.iloc[0]["source_csv_row"]), "role": "PRIMARY",
             "label": "Primary analysis", "subject_matching": "reference",
             "n_subjects": len(subj14), "analysis_subject_ids": ";".join(map(str, subj14)),
             "estimate": reference[refs[0]], "ci_lower": reference[refs[1]],
             "ci_upper": reference[refs[2]], "status": "stored_reference", "is_estimable": True}]
    for anchor, label in [("common_start", "Common start"), ("common_end", "Common end")]:
        z = a.loc[a["anchor_type"].eq(anchor) & a["subject_matching"].eq("matched_all_dt")].copy()
        audit.check(f"One matched {anchor} row", len(z) == 1)
        numeric(z, ["n_subjects_analysis"], anchor, integer=True)
        numeric(z, ["observed_slope_per_week"], anchor, allow_nan=True)
        validate_ci(z, "bootstrap_q025", "bootstrap_q975", anchor, allow_nan=True)
        r = z.iloc[0]
        n = int(r["n_subjects_analysis"])
        audit.check(f"{anchor} subject count", 0 <= n <= len(subj14))
        ids = str(r.get("analysis_subject_ids", ""))
        if ids and ids != "nan" and n:
            ids_int = integers(ids, f"{anchor} subjects")
            audit.check(f"{anchor} subject IDs", len(ids_int) == n and set(ids_int).issubset(subj14))
        if "n_observed_dt" in z and np.isfinite(r["observed_slope_per_week"]):
            audit.check(f"{anchor} uses all lags", int(r["n_observed_dt"]) == len(dts))
        status = str(r.get("test_status", "not_recorded"))
        estimable = bool(np.isfinite([r["observed_slope_per_week"], r["bootstrap_q025"], r["bootstrap_q975"]]).all())
        if not estimable:
            audit.note(f"{anchor} unavailable", f"Retained in panel A as not estimable; status={status}, n={n}")
        rows.append({"source_table": "Step14/09_anchored_core_slope_bootstrap.csv",
                     "source_csv_row": int(r["source_csv_row"]), "role": anchor,
                     "label": label, "subject_matching": "matched_all_dt", "n_subjects": n,
                     "analysis_subject_ids": ids, "estimate": r["observed_slope_per_week"],
                     "ci_lower": r["bootstrap_q025"], "ci_upper": r["bootstrap_q975"],
                     "status": status, "is_estimable": estimable})
    panel_a = pd.DataFrame(rows)

    # Optional independent Step-13 reference check; never import a new estimate.
    c13 = None
    if step13 is not None:
        c13 = load_json(step13, "00_run_config.json", audit)
        audit.check("Step13 primary estimator", get(c13, "primary_estimator") == "equal_bin_equal_subject")
        audit.check("Step13 primary metric", get(c13, "primary_metric") == "cross_cov")
        audit.check("Step13 subject universe", integers(get(c13, "complete_case.subjects"), "Step13 subjects") == subj14)
        audit.check("Step13 selected core", integers(get(c13, "core_selection.selected_bins"), "Step13 core") == bins14)
        s13 = load_csv(step13, "05_temporal_slope_summary.csv", audit)
        require_columns(s13, ["metric", "estimator", "observed_slope_per_week", "bootstrap_q025", "bootstrap_q975"], "Step13 slopes")
        z13 = s13.loc[s13["metric"].eq("cross_cov") & s13["estimator"].eq("equal_bin_equal_subject")]
        audit.check("One Step13 primary slope", len(z13) == 1)
        for old, new in zip(refs, ["observed_slope_per_week", "bootstrap_q025", "bootstrap_q975"]):
            audit.check(f"Step13 independent reference {new}", same(reference[old], z13.iloc[0][new]),
                        "Mismatch suggests Step14 used a different Step13 run. Inputs are not reconciled automatically.")
    else:
        audit.note("Independent Step13 reference check not performed", "Supply --step13-dir to compare against Figure 6 source outputs.")

    # B: three already-combined reciprocal-fold curves; retain NaN rows as gaps.
    forward = load_csv(step15, "03_forward_domain_by_bin.csv", audit)
    require_columns(forward, ["mode", "bin", "x_center", "mean_dx", "bootstrap_q025", "bootstrap_q975"], "Step15 forward")
    b = forward.loc[forward["mode"].isin(["PRIMARY", "T0PLUS", "TT"])].copy()
    audit.check("All three forward domains", set(b["mode"]) == {"PRIMARY", "T0PLUS", "TT"})
    numeric(b, ["bin"], "Step15 forward", integer=True)
    numeric(b, ["x_center"], "Step15 forward")
    numeric(b, ["mean_dx"], "Step15 forward", allow_nan=True)
    validate_ci(b, "bootstrap_q025", "bootstrap_q975", "Step15 forward", allow_nan=True)
    unique(b, ["mode", "bin"], "Step15 forward", audit)
    for mode, g in b.groupby("mode"):
        audit.check(f"Finite forward points {mode}", np.isfinite(g["mean_dx"]).sum() >= 2)
        xs = g.sort_values("bin")["x_center"].to_numpy(float)
        audit.check(f"Ordered abundance coordinates {mode}", np.all(np.diff(xs) > 0))
    audit.check("One forward coordinate per bin", b.groupby("bin")["x_center"].nunique().max() == 1)
    if "dt" in b:
        constant(b, "dt", 1, "Forward table lag matches configuration", audit)
    if "alpha" in b:
        constant(b, "alpha", alpha, "Forward table threshold matches configuration", audit)
    b["reference_alpha"] = alpha
    b["dt"] = 1
    b["point_displayed"] = np.isfinite(b["mean_dx"])
    b["ci_displayed"] = b["mode"].eq("PRIMARY") & b["point_displayed"] & np.isfinite(b["bootstrap_q025"]) & np.isfinite(b["bootstrap_q975"])
    panel_b = source(b.sort_values(["mode", "bin"]), "Step15/03_forward_domain_by_bin.csv")

    # C/D: precomputed cross-alpha common-core values; never refit slopes here.
    profiles = load_csv(step16, "12_temporal_profiles_cross_alpha.csv", audit)
    slopes = load_csv(step16, "13_temporal_slopes_cross_alpha.csv", audit)
    pcols = ["mode", "alpha", "dt", "cross_cov", "n_common_core_bins", "n_subjects_universe"]
    scols = ["mode", "alpha", "observed_slope_per_week", "bootstrap_q025", "bootstrap_q975",
             "n_common_core_bins", "n_subjects_universe"]
    require_columns(profiles, pcols, "Step16 temporal profiles")
    require_columns(slopes, scols, "Step16 temporal slopes")
    modes = ["PRIMARY", "NON_FF", "TT"]
    profiles = profiles.loc[profiles["mode"].isin(modes)].copy()
    slopes = slopes.loc[slopes["mode"].isin(modes)].copy()
    numeric(profiles, ["alpha", "cross_cov"], "Step16 profiles")
    numeric(profiles, ["dt", "n_common_core_bins", "n_subjects_universe"], "Step16 profiles", integer=True)
    numeric(slopes, ["alpha", "observed_slope_per_week"], "Step16 slopes")
    validate_ci(slopes, "bootstrap_q025", "bootstrap_q975", "Step16 slopes")
    unique(profiles, ["mode", "alpha", "dt"], "Step16 profiles", audit)
    unique(slopes, ["mode", "alpha"], "Step16 slopes", audit)
    for label, d in [("profiles", profiles), ("slopes", slopes)]:
        audit.check(f"All temporal domains in {label}", set(d["mode"]) == set(modes))
        constant(d, "n_common_core_bins", len(common), f"Common core count in {label}", audit)
        constant(d, "n_subjects_universe", len(subj14), f"Subject count in {label}", audit)
        for mode in modes:
            audit.check(f"Threshold set in {label}/{mode}", same(sorted(d.loc[d["mode"].eq(mode), "alpha"].unique()), alphas))
    for (mode, av), g in profiles.groupby(["mode", "alpha"]):
        audit.check(f"Profile lags {mode}/{av:g}", sorted(g["dt"].tolist()) == dts)
    for dt in dts:
        g = profiles.loc[profiles["mode"].eq("PRIMARY") & profiles["dt"].eq(dt)]
        constant(g.copy(), "cross_cov", float(g.iloc[0]["cross_cov"]), f"PRIMARY amplitude invariant dt={dt}", audit)
    primary_d = slopes.loc[slopes["mode"].eq("PRIMARY")]
    for col in ["observed_slope_per_week", "bootstrap_q025", "bootstrap_q975"]:
        constant(primary_d.copy(), col, float(primary_d.iloc[0][col]), f"PRIMARY slope/reference invariant {col}", audit)
    c = profiles.loc[profiles["dt"].eq(1)].copy().sort_values(["mode", "alpha"])
    c["ci_displayed"] = False
    c["display_as_horizontal_reference"] = c["mode"].eq("PRIMARY")
    c["uncertainty_display"] = "Point estimate only; no amplitude CI inferred by this plotter"
    panel_c = source(c, "Step16/12_temporal_profiles_cross_alpha.csv")
    slopes = slopes.sort_values(["mode", "alpha"])
    slopes["display_as_horizontal_reference"] = slopes["mode"].eq("PRIMARY")
    panel_d = source(slopes, "Step16/13_temporal_slopes_cross_alpha.csv")
    audit.note("Panel C uncertainty", "The inspected Step16 profile schema exports no amplitude CI; panel C is descriptive and has no error bars.")
    audit.note("Different temporal estimators", "A: Step13 equal-bin/equal-subject reference. C/D: Step16 common-core pooled covariance, equal-bin average. They are not substituted for one another.")

    # Compare forward values from two analysis runs when Step16 supplies them.
    extra_path = step16 / "04_forward_domain_cross_alpha_by_bin.csv"
    if extra_path.is_file():
        f16 = load_csv(step16, extra_path.name, audit)
        require_columns(f16, ["alpha", "mode", "bin", "x_center", "mean_dx"], "Step16 forward cross-check")
        f16 = f16.loc[np.isclose(pd.to_numeric(f16["alpha"], errors="raise"), alpha, rtol=0, atol=1e-12)
                      & f16["mode"].isin(["PRIMARY", "T0PLUS", "TT"])].copy()
        unique(f16, ["mode", "bin"], "Step16 reference forward", audit)
        comp = b.merge(f16, on=["mode", "bin"], how="outer", suffixes=("_15", "_16"), indicator=True, validate="one_to_one")
        audit.check("Step15/16 reference forward keys", comp["_merge"].eq("both").all())
        for col in ["x_center", "mean_dx"]:
            audit.check(f"Step15/16 reference forward {col}", same(comp[col + "_15"], comp[col + "_16"]),
                        "Different values may indicate mixed runs or different settings.")
    else:
        audit.note("Cross-run forward-value check unavailable", f"Optional table not supplied: {extra_path}")
    # Paths are documentary, not hashes: differences are reported, not equated
    # with a scientific mismatch (datasets may have been relocated).
    for key in ["step11_dir", "step13_dir"]:
        vals = {str(cfg[key]) for cfg in [c14, c15, c16] if cfg.get(key)}
        if len(vals) > 1:
            audit.note(f"Recorded {key} differs", " | ".join(sorted(vals)))
    if c15.get("transitions") != c16.get("transitions"):
        audit.note("Recorded transition paths differ", "Check the input manifests; path equality alone does not establish dataset identity.")
    audit.note("Provenance scope", "Hashes identify the supplied files. Schema/core/value checks are not a new audit of the full upstream analysis.")
    meta = {"reference_alpha": alpha, "alphas": alphas, "dt_values": dts,
            "step13_core_bins": bins14, "step16_common_core_bins": common,
            "complete_case_subjects": subj14, "n_subjects": len(subj14),
            "source_configs": {"step14": c14, "step15": c15, "step16": c16, "step13": c13},
            "step13_independently_checked": step13 is not None,
            "panel_c_has_confidence_intervals": False}
    return {"A": panel_a, "B": panel_b, "C": panel_c, "D": panel_d}, meta



def load_supplementary(step14: Path, step15: Path, step16: Path, meta: Dict[str, Any],
                       audit: Audit) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """Select finalized Step-14/15/16 source data for Supplementary Figures 9–11.

    This function performs schema/value checks and source selection only. It does
    not recompute permutation tests, bootstrap intervals, differences, variances,
    covariance, temporal slopes or operational classes.
    """
    # ------------------------------ S9 / Step 14 ------------------------------
    core = load_csv(step14, "06_core_interval_position_tests.csv", audit)
    require_columns(core, ["metric", "dt", "test_status",
                           "q_bh_rms_position_structure"], "Step14 core position tests")
    numeric(core, ["dt"], "Step14 core position tests", integer=True)
    numeric(core, ["q_bh_rms_position_structure"], "Step14 core position tests", allow_nan=True)
    s9a = core.loc[
        core["metric"].isin(SUPP9_METRICS)
        & core["test_status"].astype(str).eq("ok")
        & np.isfinite(core["q_bh_rms_position_structure"])
    ].copy()
    audit.check("S9A three metrics", set(s9a["metric"]) == set(SUPP9_METRICS))
    audit.check("S9A evaluable lags", sorted(s9a["dt"].unique().tolist()) == [1, 2, 3])
    unique(s9a, ["metric", "dt"], "S9A core q values", audit)
    s9a = source(s9a.sort_values(["metric", "dt"]), "Step14/06_core_interval_position_tests.csv")

    binpos = load_csv(step14, "07_binwise_interval_position_tests.csv", audit)
    require_columns(binpos, ["metric", "dt", "bin", "x_center", "test_status",
                             "q_bh_rms_position_structure"], "Step14 binwise position tests")
    numeric(binpos, ["dt", "bin"], "Step14 binwise position tests", integer=True)
    numeric(binpos, ["x_center", "q_bh_rms_position_structure"],
            "Step14 binwise position tests", allow_nan=True)
    s9b = binpos.loc[
        binpos["metric"].eq("cross_cov")
        & binpos["test_status"].astype(str).eq("ok")
        & np.isfinite(binpos["q_bh_rms_position_structure"])
    ].copy()
    audit.check("S9B evaluable lags", sorted(s9b["dt"].unique().tolist()) == [1, 2, 3])
    unique(s9b, ["dt", "bin"], "S9B binwise q values", audit)
    s9b = source(s9b.sort_values(["dt", "bin"]), "Step14/07_binwise_interval_position_tests.csv")

    binsl = load_csv(step14, "12_binwise_anchor_slope_bootstrap.csv", audit)
    require_columns(binsl, ["anchor_type", "subject_matching", "metric", "bin", "x_center",
                            "test_status", "n_subjects_analysis", "observed_slope_per_week",
                            "bootstrap_q025", "bootstrap_q975"], "Step14 binwise anchored slopes")
    numeric(binsl, ["bin", "n_subjects_analysis"], "Step14 binwise anchored slopes", integer=True)
    numeric(binsl, ["x_center", "observed_slope_per_week"],
            "Step14 binwise anchored slopes", allow_nan=True)
    validate_ci(binsl, "bootstrap_q025", "bootstrap_q975", "Step14 binwise anchored slopes", allow_nan=True)
    s9_cd: Dict[str, pd.DataFrame] = {}
    for anchor, panel in [("common_start", "C"), ("common_end", "D")]:
        d = binsl.loc[
            binsl["anchor_type"].eq(anchor)
            & binsl["subject_matching"].eq("matched_all_dt")
            & binsl["metric"].eq("cross_cov")
            & binsl["test_status"].astype(str).eq("ok")
        ].copy()
        audit.check(f"S9{panel} has core bins", len(d) == len(meta["step13_core_bins"]),
                    f"found={len(d)}, expected={len(meta['step13_core_bins'])}")
        unique(d, ["bin"], f"S9{panel} anchored bin slopes", audit)
        s9_cd[panel] = source(d.sort_values("bin"), "Step14/12_binwise_anchor_slope_bootstrap.csv")

    # ------------------------------ S10 / Step 15 -----------------------------
    classcomp = load_csv(step15, "01_operational_class_composition_by_dt.csv", audit)
    require_columns(classcomp, ["dt", "obs_class_step15", "fraction_within_dt", "alpha"],
                    "Step15 class composition")
    numeric(classcomp, ["dt"], "Step15 class composition", integer=True)
    numeric(classcomp, ["fraction_within_dt", "alpha"], "Step15 class composition")
    constant(classcomp.copy(), "alpha", meta["reference_alpha"], "S10A reference alpha", audit)
    audit.check("S10A classes", set(classcomp["obs_class_step15"]) == {"TT", "TF", "FT", "FF"})
    unique(classcomp, ["dt", "obs_class_step15"], "S10A class composition", audit)
    s10a = source(classcomp.sort_values(["obs_class_step15", "dt"]),
                  "Step15/01_operational_class_composition_by_dt.csv")

    crossing = load_csv(step15, "02_boundary_crossing_by_dt_bin.csv", audit)
    require_columns(crossing, ["dt", "bin", "x_center", "crossing_imbalance", "alpha"],
                    "Step15 boundary crossing")
    numeric(crossing, ["dt", "bin"], "Step15 boundary crossing", integer=True)
    numeric(crossing, ["x_center", "crossing_imbalance", "alpha"],
            "Step15 boundary crossing", allow_nan=True)
    constant(crossing.copy(), "alpha", meta["reference_alpha"], "S10B reference alpha", audit)
    unique(crossing, ["dt", "bin"], "S10B crossing imbalance", audit)
    s10b = source(crossing.sort_values(["dt", "bin"]), "Step15/02_boundary_crossing_by_dt_bin.csv")

    fdiff = load_csv(step15, "05_forward_difference_vs_primary.csv", audit)
    require_columns(fdiff, ["contrast", "bin", "x_center", "difference",
                            "bootstrap_difference_q025", "bootstrap_difference_q975"],
                    "Step15 forward differences")
    numeric(fdiff, ["bin"], "Step15 forward differences", integer=True)
    numeric(fdiff, ["x_center", "difference"], "Step15 forward differences", allow_nan=True)
    validate_ci(fdiff, "bootstrap_difference_q025", "bootstrap_difference_q975",
                "Step15 forward differences", allow_nan=True)
    wanted_contrasts = set(SUPP10_CONTRASTS)
    s10c = fdiff.loc[fdiff["contrast"].isin(wanted_contrasts)].copy()
    audit.check("S10C contrasts", set(s10c["contrast"]) == wanted_contrasts)
    unique(s10c, ["contrast", "bin"], "S10C forward differences", audit)
    s10c = source(s10c.sort_values(["contrast", "bin"]), "Step15/05_forward_difference_vs_primary.csv")

    fluctdiff = load_csv(step15, "08_fluctuation_difference_vs_primary.csv", audit)
    require_columns(fluctdiff, ["mode", "metric", "dt", "bin", "x_center", "difference",
                              "bootstrap_difference_q025", "bootstrap_difference_q975"],
                    "Step15 fluctuation differences")
    numeric(fluctdiff, ["dt", "bin"], "Step15 fluctuation differences", integer=True)
    numeric(fluctdiff, ["x_center", "difference"], "Step15 fluctuation differences", allow_nan=True)
    validate_ci(fluctdiff, "bootstrap_difference_q025", "bootstrap_difference_q975",
                "Step15 fluctuation differences", allow_nan=True)
    s10d = fluctdiff.loc[
        fluctdiff["metric"].eq("cross_cov")
        & fluctdiff["dt"].eq(1)
        & fluctdiff["mode"].isin(["NON_FF", "TT"])
    ].copy()
    audit.check("S10D modes", set(s10d["mode"]) == {"NON_FF", "TT"})
    unique(s10d, ["mode", "bin"], "S10D covariance differences", audit)
    s10d = source(s10d.sort_values(["mode", "bin"]), "Step15/08_fluctuation_difference_vs_primary.csv")

    fluct = load_csv(step15, "07_fluctuation_domain_by_bin_dt.csv", audit)
    req = ["mode", "dt", "bin", "x_center", "same_var_mean", "replicate_specific_excess",
           "same_var_mean_ci025", "same_var_mean_ci975",
           "replicate_specific_excess_ci025", "replicate_specific_excess_ci975"]
    require_columns(fluct, req, "Step15 fluctuation domains")
    numeric(fluct, ["dt", "bin"], "Step15 fluctuation domains", integer=True)
    numeric(fluct, ["x_center", "same_var_mean", "replicate_specific_excess"],
            "Step15 fluctuation domains", allow_nan=True)
    for lo, hi, label in [("same_var_mean_ci025", "same_var_mean_ci975", "same_var"),
                          ("replicate_specific_excess_ci025", "replicate_specific_excess_ci975", "excess")]:
        validate_ci(fluct, lo, hi, f"Step15 {label} intervals", allow_nan=True)
    base_fluct = fluct.loc[fluct["dt"].eq(1) & fluct["mode"].isin(["PRIMARY", "NON_FF", "TT"])].copy()
    audit.check("S10E/F modes", set(base_fluct["mode"]) == {"PRIMARY", "NON_FF", "TT"})
    unique(base_fluct, ["mode", "bin"], "S10E/F fluctuation profiles", audit)
    s10e = source(base_fluct.sort_values(["mode", "bin"]), "Step15/07_fluctuation_domain_by_bin_dt.csv")
    s10f = s10e.copy()

    temporal = load_csv(step15, "12_temporal_domain_by_dt.csv", audit)
    require_columns(temporal, ["mode", "dt", "value", "n_common_core_bins", "n_subjects_universe"],
                    "Step15 temporal profiles")
    numeric(temporal, ["dt", "n_common_core_bins", "n_subjects_universe"],
            "Step15 temporal profiles", integer=True)
    numeric(temporal, ["value"], "Step15 temporal profiles")
    temporal = temporal.loc[temporal["mode"].isin(["PRIMARY", "NON_FF", "TT"])].copy()
    audit.check("S10G modes", set(temporal["mode"]) == {"PRIMARY", "NON_FF", "TT"})
    for mode, g in temporal.groupby("mode"):
        audit.check(f"S10G lags {mode}", sorted(g["dt"].tolist()) == meta["dt_values"])
    unique(temporal, ["mode", "dt"], "S10G temporal profiles", audit)
    s10g = source(temporal.sort_values(["mode", "dt"]), "Step15/12_temporal_domain_by_dt.csv")

    # ------------------------------ S11 / Step 16 -----------------------------
    ttret = load_csv(step16, "03_tt_retention_by_alpha.csv", audit)
    require_columns(ttret, ["alpha", "alpha_tag", "reference_alpha", "n_tt",
                            "reference_n_tt", "n_tt_intersection_with_reference",
                            "fraction_reference_tt_retained",
                            "fraction_current_tt_shared_with_reference"],
                    "Step16 TT retention")
    numeric(ttret, ["alpha", "reference_alpha", "n_tt", "reference_n_tt",
                    "n_tt_intersection_with_reference",
                    "fraction_reference_tt_retained",
                    "fraction_current_tt_shared_with_reference"], "Step16 TT retention")
    constant(ttret.copy(), "reference_alpha", meta["reference_alpha"],
             "S11A reference alpha", audit)
    audit.check("S11A threshold set", same(sorted(ttret["alpha"].unique()), meta["alphas"]))
    unique(ttret, ["alpha"], "S11A TT retention", audit)
    s11a = source(ttret.sort_values("alpha"), "Step16/03_tt_retention_by_alpha.csv")

    fcommon = load_csv(step16, "05b_forward_cross_alpha_common_support.csv", audit)
    require_columns(fcommon, ["bin", "x_center", "selected_cross_alpha_forward_common_support"],
                    "Step16 forward common support")
    numeric(fcommon, ["bin", "x_center"], "Step16 forward common support")
    flag = fcommon["selected_cross_alpha_forward_common_support"].astype(str).str.strip().str.lower()
    audit.check("S11 forward common-support flags",
                flag.isin(["true", "false", "1", "0", "1.0", "0.0"]).all())
    selected_f = fcommon.loc[flag.isin(["true", "1", "1.0"])].copy()
    selected_forward_bins = sorted(pd.to_numeric(selected_f["bin"], errors="raise").astype(int).tolist())
    audit.check("S11 forward common support nonempty", len(selected_forward_bins) >= 2)
    audit.check("S11 forward common support contiguous",
                selected_forward_bins == list(range(selected_forward_bins[0], selected_forward_bins[-1] + 1)))
    meta["step16_forward_common_support_bins"] = selected_forward_bins

    fwd = load_csv(step16, "04_forward_domain_cross_alpha_by_bin.csv", audit)
    require_columns(fwd, ["mode", "bin", "x_center", "mean_dx", "bootstrap_q025",
                          "bootstrap_q975", "alpha", "alpha_tag"], "Step16 forward cross-alpha")
    numeric(fwd, ["bin"], "Step16 forward cross-alpha", integer=True)
    numeric(fwd, ["x_center", "mean_dx", "alpha"], "Step16 forward cross-alpha", allow_nan=True)
    validate_ci(fwd, "bootstrap_q025", "bootstrap_q975", "Step16 forward cross-alpha", allow_nan=True)
    s11_forward: Dict[str, pd.DataFrame] = {}
    for mode, panel in [("T0PLUS", "B"), ("TT", "C"), ("TF", "D")]:
        d = fwd.loc[fwd["mode"].eq(mode) & fwd["bin"].isin(selected_forward_bins)].copy()
        audit.check(f"S11{panel} threshold set", same(sorted(d["alpha"].unique()), meta["alphas"]))
        for av, g in d.groupby("alpha"):
            audit.check(f"S11{panel} common bins alpha={av:g}", sorted(g["bin"].tolist()) == selected_forward_bins)
        unique(d, ["alpha", "bin"], f"S11{panel} {mode} forward", audit)
        s11_forward[panel] = source(d.sort_values(["alpha", "bin"]),
                                    "Step16/04_forward_domain_cross_alpha_by_bin.csv")

    fsdiff = load_csv(step16, "07_forward_slope_difference_vs_reference_alpha.csv", audit)
    require_columns(fsdiff, ["alpha", "reference_alpha", "mode", "observed_slope_difference",
                             "bootstrap_slope_difference_q025", "bootstrap_slope_difference_q975"],
                    "Step16 forward slope differences")
    numeric(fsdiff, ["alpha", "reference_alpha", "observed_slope_difference"],
            "Step16 forward slope differences")
    validate_ci(fsdiff, "bootstrap_slope_difference_q025", "bootstrap_slope_difference_q975",
                "Step16 forward slope differences")
    constant(fsdiff.copy(), "reference_alpha", meta["reference_alpha"],
             "S11E reference alpha", audit)
    s11e = fsdiff.loc[fsdiff["mode"].isin(["T0PLUS", "TT", "TF"])].copy()
    audit.check("S11E modes", set(s11e["mode"]) == {"T0PLUS", "TT", "TF"})
    nonref_alphas = [a for a in meta["alphas"] if not same(a, meta["reference_alpha"])]
    for mode, g in s11e.groupby("mode"):
        audit.check(f"S11E non-reference thresholds {mode}", same(sorted(g["alpha"].unique()), sorted(nonref_alphas)))
    unique(s11e, ["mode", "alpha"], "S11E forward slope differences", audit)
    s11e = source(s11e.sort_values(["mode", "alpha"]),
                  "Step16/07_forward_slope_difference_vs_reference_alpha.csv")

    fluct16 = load_csv(step16, "08_fluctuation_domain_cross_alpha_by_bin_dt.csv", audit)
    require_columns(fluct16, ["mode", "dt", "bin", "x_center", "meets_min_support",
                              "cross_cov", "alpha", "alpha_tag"], "Step16 fluctuation domains")
    numeric(fluct16, ["dt", "bin"], "Step16 fluctuation domains", integer=True)
    numeric(fluct16, ["x_center", "cross_cov", "alpha"], "Step16 fluctuation domains", allow_nan=True)
    support_flag = fluct16["meets_min_support"].astype(str).str.strip().str.lower()
    audit.check("S11F support flags", support_flag.isin(["true", "false", "1", "0", "1.0", "0.0"]).all())
    s11f = fluct16.loc[
        fluct16["mode"].eq("NON_FF") & fluct16["dt"].eq(1)
        & support_flag.isin(["true", "1", "1.0"])
    ].copy()
    audit.check("S11F threshold set", same(sorted(s11f["alpha"].unique()), meta["alphas"]))
    unique(s11f, ["alpha", "bin"], "S11F NON_FF one-week covariance", audit)
    s11f = source(s11f.sort_values(["alpha", "bin"]),
                  "Step16/08_fluctuation_domain_cross_alpha_by_bin_dt.csv")

    fdiff16 = load_csv(step16, "09_fluctuation_difference_vs_reference_alpha.csv", audit)
    require_columns(fdiff16, ["alpha", "reference_alpha", "mode", "metric", "dt", "bin",
                              "x_center", "difference", "bootstrap_difference_q025",
                              "bootstrap_difference_q975"], "Step16 fluctuation differences")
    numeric(fdiff16, ["alpha", "reference_alpha", "dt", "bin"],
            "Step16 fluctuation differences", allow_nan=False)
    numeric(fdiff16, ["x_center", "difference"], "Step16 fluctuation differences", allow_nan=True)
    validate_ci(fdiff16, "bootstrap_difference_q025", "bootstrap_difference_q975",
                "Step16 fluctuation differences", allow_nan=True)
    s11g = fdiff16.loc[
        fdiff16["mode"].eq("TT") & fdiff16["metric"].eq("cross_cov") & fdiff16["dt"].eq(1)
    ].copy()
    audit.check("S11G non-reference thresholds", same(sorted(s11g["alpha"].unique()), sorted(nonref_alphas)))
    unique(s11g, ["alpha", "bin"], "S11G TT covariance differences", audit)
    s11g = source(s11g.sort_values(["alpha", "bin"]),
                  "Step16/09_fluctuation_difference_vs_reference_alpha.csv")

    temporal16 = load_csv(step16, "12_temporal_profiles_cross_alpha.csv", audit)
    require_columns(temporal16, ["mode", "dt", "cross_cov", "n_common_core_bins",
                                 "n_subjects_universe", "alpha", "alpha_tag"],
                    "Step16 temporal profiles for S11")
    numeric(temporal16, ["dt", "n_common_core_bins", "n_subjects_universe"],
            "Step16 temporal profiles for S11", integer=True)
    numeric(temporal16, ["cross_cov", "alpha"], "Step16 temporal profiles for S11")
    s11h = temporal16.loc[temporal16["mode"].eq("TT")].copy()
    audit.check("S11H threshold set", same(sorted(s11h["alpha"].unique()), meta["alphas"]))
    constant(s11h.copy(), "n_common_core_bins", len(meta["step16_common_core_bins"]),
             "S11H common core count", audit)
    constant(s11h.copy(), "n_subjects_universe", meta["n_subjects"],
             "S11H subject count", audit)
    for av, g in s11h.groupby("alpha"):
        audit.check(f"S11H lags alpha={av:g}", sorted(g["dt"].tolist()) == meta["dt_values"])
    unique(s11h, ["alpha", "dt"], "S11H TT temporal profiles", audit)
    s11h = source(s11h.sort_values(["alpha", "dt"]),
                  "Step16/12_temporal_profiles_cross_alpha.csv")

    audit.note("Supplementary plotting scope",
               "Supplementary Figures 9–11 select already-computed Step14/15/16 outputs only; no scientific statistic is recomputed.")
    s9 = {"A": s9a, "B": s9b, "C": s9_cd["C"], "D": s9_cd["D"]}
    s10 = {"A": s10a, "B": s10b, "C": s10c, "D": s10d,
           "E": s10e, "F": s10f, "G": s10g}
    s11 = {"A": s11a, "B": s11_forward["B"], "C": s11_forward["C"],
           "D": s11_forward["D"], "E": s11e, "F": s11f, "G": s11g, "H": s11h}
    return s9, s10, s11


# ============================== PLOTTING ====================================
def configure_style(font_override: Optional[str], dpi: int) -> str:
    requested = font_override or FONT["family"]
    try:
        font_manager.findfont(font_manager.FontProperties(family=requested), fallback_to_default=False)
        family = requested
    except ValueError:
        family = "DejaVu Sans"
        print(f"[NOTE] Font {requested!r} unavailable; using {family}.", file=sys.stderr)
    mpl.rcParams.update({"font.family": family, "font.size": FONT["size"],
        "axes.labelsize": FONT["axis"], "xtick.labelsize": FONT["tick"],
        "ytick.labelsize": FONT["tick"], "legend.fontsize": FONT["legend"],
        "axes.linewidth": LINE_WIDTHS["axis"], "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.dpi": dpi, "text.usetex": False})
    return family


def style_for(key: str) -> Dict[str, Any]:
    cfg = SERIES[key]
    colours = mpl.rcParams["axes.prop_cycle"].by_key().get("color", ["C0"])
    return {"color": colours[cfg["cycle_index"] % len(colours)],
            "marker": cfg["marker"], "linestyle": cfg["linestyle"]}


def new_panel(panel: str):
    fig, ax = plt.subplots(figsize=FIGURE_SIZES[panel])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.grid(False)
    return fig, ax


def neutral(ax) -> Any:
    return ax.spines["bottom"].get_edgecolor()


def draw_box(ax, panel: str, default_text: str) -> None:
    cfg = BOX_CONFIG[panel]
    if not cfg["enabled"]:
        return
    text = default_text if cfg.get("text") is None else str(cfg["text"])
    if text:
        ax.text(float(cfg["x"]), float(cfg["y"]), text, transform=ax.transAxes,
                ha=cfg["ha"], va=cfg["va"], fontsize=FONT["annotation"],
                bbox={"boxstyle": "round,pad=0.25", "facecolor": mpl.rcParams["figure.facecolor"],
                      "edgecolor": "none", "alpha": 0.85}, zorder=10)


def legend(ax, panel: str) -> None:
    cfg = LEGEND_CONFIG[panel]
    if cfg["enabled"]:
        kw = {"loc": cfg["loc"], "ncol": cfg["ncol"], "frameon": False, "handlelength": 2}
        if cfg["bbox_to_anchor"] is not None:
            kw["bbox_to_anchor"] = cfg["bbox_to_anchor"]
        ax.legend(**kw)


def apply_supp_legend(ax: plt.Axes, key: str) -> None:
    """Apply user-editable legend placement to one supplementary panel.

    Only legend visibility/position/columns are changed; plotted data and
    scientific quantities are untouched.
    """
    if key not in SUPP_LEGEND_CONFIG:
        raise InputError(f"Missing supplementary legend configuration: {key}")
    cfg = SUPP_LEGEND_CONFIG[key]
    current = ax.get_legend()
    if not cfg.get("enabled", True):
        if current is not None:
            current.remove()
        return

    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        if current is not None:
            current.remove()
        return

    ncol = int(cfg.get("ncol", 1))
    if ncol < 1:
        raise InputError(f"{key}.ncol must be >= 1")
    kw = {"loc": cfg.get("loc", "best"), "ncol": ncol,
          "frameon": False, "handlelength": 2}
    if cfg.get("bbox_to_anchor") is not None:
        kw["bbox_to_anchor"] = cfg["bbox_to_anchor"]
    ax.legend(handles, labels, **kw)


def intervals(ax, x, y, lo, hi, color, horizontal=False) -> None:
    # Draw endpoints directly. This also represents intervals that do not
    # contain the point estimate without creating invalid negative yerr.
    for xx, yy, lower, upper in zip(x, y, lo, hi):
        if np.isfinite([xx, yy, lower, upper]).all():
            if horizontal:
                ax.plot([lower, upper], [yy, yy], color=color, linewidth=LINE_WIDTHS["ci"], zorder=2)
                ax.plot([lower, upper], [yy, yy], "|", color=color, markersize=5, markeredgewidth=0.9)
            else:
                ax.plot([xx, xx], [lower, upper], color=color, linewidth=LINE_WIDTHS["ci"], zorder=2)
                ax.plot([xx, xx], [lower, upper], "_", color=color, markersize=5, markeredgewidth=0.9)


def draw_a(data: pd.DataFrame, meta: Dict[str, Any]):
    fig, ax = new_panel("A")
    ax.axvline(0, color=neutral(ax), linestyle=":", linewidth=LINE_WIDTHS["reference"], alpha=0.6)
    ys = np.arange(len(data))[::-1]
    labels = []
    for yy, r in zip(ys, data.itertuples(index=False)):
        st = style_for(r.role)
        labels.append((f"Primary analysis\n(n={r.n_subjects})" if r.role == "PRIMARY" else
                       f"{r.label}, matched\n(n={r.n_subjects})"))
        if r.is_estimable:
            intervals(ax, [r.estimate], [yy], [r.ci_lower], [r.ci_upper], st["color"], horizontal=True)
            ax.plot(r.estimate, yy, marker=st["marker"], linestyle="None", color=st["color"], markersize=MARKER_SIZE + 0.5, zorder=4)
        else:
            ax.text(0.5, yy, "Not estimable", transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=7)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.65, len(data) - 0.35)
    ax.set_xlabel("Temporal covariance slope\n(squared log-frequency units per week)")
    ax.margins(x=0.12)
    draw_box(ax, "A", f"Empirical core: {len(meta['step13_core_bins'])} bins")
    return fig, ax


def draw_b(data: pd.DataFrame, meta: Dict[str, Any]):
    fig, ax = new_panel("B")
    ax.axhline(0, color=neutral(ax), linestyle=":", linewidth=LINE_WIDTHS["reference"], alpha=0.6)
    coords = data[["bin", "x_center"]].drop_duplicates().set_index("bin")["x_center"]
    for mode in ["PRIMARY", "T0PLUS", "TT"]:
        g = data.loc[data["mode"].eq(mode)].set_index("bin").sort_index()
        # Missing abundance bins are gaps, not straight-line interpolation.
        full_index = range(int(g.index.min()), int(g.index.max()) + 1)
        g = g.reindex(full_index)
        g["x_center"] = coords.reindex(full_index).to_numpy()
        st = style_for(mode)
        ax.plot(g["x_center"], g["mean_dx"], **st, markersize=MARKER_SIZE,
                linewidth=LINE_WIDTHS["primary" if mode == "PRIMARY" else "curve"], label=mode)
        if mode == "PRIMARY":
            ax.fill_between(g["x_center"].to_numpy(float), g["bootstrap_q025"].to_numpy(float),
                            g["bootstrap_q975"].to_numpy(float),
                            where=np.isfinite(g["mean_dx"].to_numpy(float)),
                            color=st["color"], alpha=CI_ALPHA, linewidth=0)
    ax.set_xlabel(r"Initial observed log-frequency, $x_0$")
    ax.set_ylabel("Mean replicate-decoupled\ndisplacement")
    legend(ax, "B")
    draw_box(ax, "B", rf"$\Delta t=1$ week; $\alpha={meta['reference_alpha']:g}$")
    return fig, ax


def alpha_axes(ax, meta: Dict[str, Any]):
    aa = np.asarray(meta["alphas"], dtype=float)
    pad = max(float(np.ptp(aa)) * 0.06, 0.001)
    ax.set_xlim(float(aa.min() - pad), float(aa.max() + pad))
    ax.set_xticks(aa)
    ax.set_xticklabels([f"{x:g}" for x in aa])
    ax.set_xlabel(r"Operational observation threshold, $\alpha$")
    ax.axvline(meta["reference_alpha"], color=neutral(ax), linestyle=":", linewidth=0.6, alpha=0.35, zorder=0)


def draw_c(data: pd.DataFrame, meta: Dict[str, Any]):
    fig, ax = new_panel("C")
    for mode in ["PRIMARY", "NON_FF", "TT"]:
        g = data.loc[data["mode"].eq(mode)].sort_values("alpha")
        st = style_for(mode)
        if mode == "PRIMARY":
            ax.axhline(float(g.iloc[0]["cross_cov"]), color=st["color"], linestyle="--",
                       linewidth=LINE_WIDTHS["primary"], label=mode)
        else:
            ax.plot(g["alpha"], g["cross_cov"], **st, markersize=MARKER_SIZE, linewidth=LINE_WIDTHS["curve"], label=mode)
    ax.set_ylabel("Cross-replicate covariance\nat one week")
    # Include zero without truncating signed estimates.
    low, high = ax.get_ylim()
    ax.set_ylim(min(0.0, low), max(0.0, high))
    alpha_axes(ax, meta)
    legend(ax, "C")
    draw_box(ax, "C", f"Shared core: {len(meta['step16_common_core_bins'])} bins")
    return fig, ax


def draw_d(data: pd.DataFrame, meta: Dict[str, Any]):
    fig, ax = new_panel("D")
    ax.axhline(0, color=neutral(ax), linestyle=":", linewidth=LINE_WIDTHS["reference"], alpha=0.6)
    for mode in ["PRIMARY", "NON_FF", "TT"]:
        g = data.loc[data["mode"].eq(mode)].sort_values("alpha")
        st = style_for(mode)
        if mode == "PRIMARY":
            r = g.iloc[0]
            ax.axhspan(float(r["bootstrap_q025"]), float(r["bootstrap_q975"]), color=st["color"], alpha=CI_ALPHA, linewidth=0)
            ax.axhline(float(r["observed_slope_per_week"]), color=st["color"], linestyle="--", linewidth=LINE_WIDTHS["primary"], label=mode)
        else:
            intervals(ax, g["alpha"], g["observed_slope_per_week"], g["bootstrap_q025"], g["bootstrap_q975"], st["color"])
            ax.plot(g["alpha"], g["observed_slope_per_week"], **st, linewidth=LINE_WIDTHS["curve"], markersize=MARKER_SIZE, label=mode)
    alpha_axes(ax, meta)
    ax.set_ylabel("Temporal covariance slope\n(squared log-frequency units per week)")
    legend(ax, "D")
    draw_box(ax, "D", f"{len(meta['step16_common_core_bins'])} shared bins; {meta['n_subjects']} subjects")
    return fig, ax



def cycle_colour(index: int) -> str:
    colours = mpl.rcParams["axes.prop_cycle"].by_key().get("color", ["C0"])
    return colours[index % len(colours)]


def prep_supp_ax(ax: plt.Axes, letter: Optional[str] = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)
    ax.grid(False)
    if letter:
        ax.text(-0.13, 1.05, letter, transform=ax.transAxes, ha="left", va="bottom",
                fontweight="bold", fontsize=FONT["panel_label"], clip_on=False)


def draw_s9_axis(ax: plt.Axes, panel: str, data: pd.DataFrame, show_letter: bool = True) -> None:
    prep_supp_ax(ax, panel if show_letter else None)
    if panel == "A":
        for metric, cfg in SUPP9_METRICS.items():
            g = data.loc[data["metric"].eq(metric)].sort_values("dt")
            ax.plot(g["dt"], g["q_bh_rms_position_structure"], marker=cfg["marker"],
                    color=cycle_colour(cfg["cycle_index"]), linewidth=LINE_WIDTHS["curve"],
                    markersize=MARKER_SIZE, label=cfg["label"])
        ax.axhline(0.05, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.65)
        ax.set_xlabel("Temporal lag (weeks)")
        ax.set_ylabel(r"BH-adjusted $q_{\mathrm{RMS}}$")
        ax.set_xticks([1, 2, 3])
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False, loc="upper left", handlelength=2)
    elif panel == "B":
        for i, dt in enumerate([1, 2, 3]):
            g = data.loc[data["dt"].eq(dt)].sort_values("bin")
            ax.plot(g["x_center"], g["q_bh_rms_position_structure"], marker=["o", "s", "^"][i],
                    color=cycle_colour(i), linewidth=LINE_WIDTHS["curve"], markersize=MARKER_SIZE,
                    label=f"{dt} week" if dt == 1 else f"{dt} weeks")
        ax.axhline(0.05, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.65)
        ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
        ax.set_ylabel(r"BH-adjusted $q_{\mathrm{RMS}}$")
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False, loc="upper left", handlelength=2)
    elif panel in ("C", "D"):
        g = data.sort_values("bin")
        color = cycle_colour(0 if panel == "C" else 1)
        intervals(ax, g["x_center"], g["observed_slope_per_week"],
                  g["bootstrap_q025"], g["bootstrap_q975"], color)
        ax.plot(g["x_center"], g["observed_slope_per_week"], marker="o", linestyle="None",
                color=color, markersize=MARKER_SIZE)
        ax.axhline(0, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.65)
        ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
        ax.set_ylabel("Anchored cross-covariance slope\n(per week)")
    else:
        raise ValueError(panel)


def draw_s10_axis(ax: plt.Axes, panel: str, data: pd.DataFrame, show_letter: bool = True) -> None:
    prep_supp_ax(ax, panel if show_letter else None)
    if panel == "A":
        order = ["TT", "TF", "FT", "FF"]
        for cls in order:
            g = data.loc[data["obs_class_step15"].eq(cls)].sort_values("dt")
            cfg = SUPP10_CLASSES[cls]
            ax.plot(g["dt"], g["fraction_within_dt"], marker=cfg["marker"],
                    color=cycle_colour(cfg["cycle_index"]), linewidth=LINE_WIDTHS["curve"],
                    markersize=MARKER_SIZE, label=cls)
        ax.set_xlabel("Temporal lag (weeks)")
        ax.set_ylabel("Fraction of transitions")
        ax.set_xticks([1, 2, 3, 4, 5])
        ax.set_ylim(0, 1.05)
        ax.legend(frameon=False, loc="upper right", ncol=1)
    elif panel == "B":
        for i, dt in enumerate([1, 2, 3, 4, 5]):
            g = data.loc[data["dt"].eq(dt)].sort_values("bin")
            ax.plot(g["x_center"], g["crossing_imbalance"], marker=["o", "s", "^", "v", "D"][i],
                    color=cycle_colour(i), linewidth=LINE_WIDTHS["curve"], markersize=MARKER_SIZE,
                    label=f"{dt} week" if dt == 1 else f"{dt} weeks")
        ax.axhline(0, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.65)
        ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
        ax.set_ylabel(r"Crossing imbalance, $(FT-TF)/(FT+TF)$")
        ax.set_ylim(-1.05, 1.05)
        ax.legend(frameon=False, loc="lower right", handlelength=2)
    elif panel == "C":
        for contrast, cfg in SUPP10_CONTRASTS.items():
            g = data.loc[data["contrast"].eq(contrast)].sort_values("bin")
            color = cycle_colour(cfg["cycle_index"])
            x = g["x_center"].to_numpy(float)
            y = g["difference"].to_numpy(float)
            lo = g["bootstrap_difference_q025"].to_numpy(float)
            hi = g["bootstrap_difference_q975"].to_numpy(float)
            ax.fill_between(x, lo, hi, color=color, alpha=0.12, linewidth=0)
            ax.plot(x, y, marker=cfg["marker"], color=color, linewidth=LINE_WIDTHS["curve"],
                    markersize=MARKER_SIZE, label=contrast)
        ax.axhline(0, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.65)
        ax.set_xlabel(r"Initial observed log-frequency, $x_0$")
        ax.set_ylabel("Difference in mean displacement")
        ax.legend(frameon=False, loc="lower right", handlelength=2)
    elif panel == "D":
        for mode in ["NON_FF", "TT"]:
            g = data.loc[data["mode"].eq(mode)].sort_values("bin")
            cfg = SERIES[mode]
            color = cycle_colour(cfg["cycle_index"])
            x = g["x_center"].to_numpy(float)
            y = g["difference"].to_numpy(float)
            lo = g["bootstrap_difference_q025"].to_numpy(float)
            hi = g["bootstrap_difference_q975"].to_numpy(float)
            ax.fill_between(x, lo, hi, color=color, alpha=0.12, linewidth=0)
            ax.plot(x, y, marker=cfg["marker"], color=color, linewidth=LINE_WIDTHS["curve"],
                    markersize=MARKER_SIZE, label=f"{mode} - PRIMARY")
        ax.axhline(0, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.65)
        ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
        ax.set_ylabel("Difference in cross-replicate covariance")
        ax.legend(frameon=False, loc="upper right", handlelength=2)
    elif panel in ("E", "F"):
        metric = "same_var_mean" if panel == "E" else "replicate_specific_excess"
        lo_col = "same_var_mean_ci025" if panel == "E" else "replicate_specific_excess_ci025"
        hi_col = "same_var_mean_ci975" if panel == "E" else "replicate_specific_excess_ci975"
        for mode in ["PRIMARY", "NON_FF", "TT"]:
            g = data.loc[data["mode"].eq(mode)].sort_values("bin")
            cfg = SERIES[mode]
            color = cycle_colour(cfg["cycle_index"])
            x = g["x_center"].to_numpy(float)
            y = g[metric].to_numpy(float)
            ax.plot(x, y, marker=cfg["marker"], color=color,
                    linewidth=LINE_WIDTHS["primary" if mode == "PRIMARY" else "curve"],
                    markersize=MARKER_SIZE, label=mode)
            if mode == "PRIMARY":
                lo = g[lo_col].to_numpy(float)
                hi = g[hi_col].to_numpy(float)
                ax.fill_between(x, lo, hi, color=color, alpha=CI_ALPHA, linewidth=0)
        ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
        ax.set_ylabel("Mean within-replicate variance" if panel == "E" else "Replicate-specific excess")
        ax.legend(frameon=False, loc="upper right", handlelength=2)
    elif panel == "G":
        for mode in ["PRIMARY", "NON_FF", "TT"]:
            g = data.loc[data["mode"].eq(mode)].sort_values("dt")
            cfg = SERIES[mode]
            ax.plot(g["dt"], g["value"], marker=cfg["marker"], color=cycle_colour(cfg["cycle_index"]),
                    linewidth=LINE_WIDTHS["primary" if mode == "PRIMARY" else "curve"],
                    markersize=MARKER_SIZE, label=mode)
        ax.set_xlabel("Temporal lag (weeks)")
        ax.set_ylabel("Common-core cross-replicate covariance")
        ax.set_xticks([1, 2, 3, 4, 5])
        ax.legend(frameon=False, loc="best", handlelength=2)
    else:
        raise ValueError(panel)


def alpha_plot_order(meta: Dict[str, Any]) -> List[float]:
    return sorted([float(a) for a in meta["alphas"]], reverse=True)


def alpha_style(alpha: float, meta: Dict[str, Any]) -> Dict[str, Any]:
    alpha = float(alpha)
    ref = float(meta["reference_alpha"])
    if same(alpha, ref):
        return {"color": "#222222", "marker": "s", "linestyle": "-"}
    nonref = [a for a in alpha_plot_order(meta) if not same(a, ref)]
    palette = [cycle_colour(0), cycle_colour(2), cycle_colour(1), cycle_colour(4), cycle_colour(5)]
    markers = ["o", "^", "D", "v", "P"]
    idx = min(nonref.index(alpha), len(palette) - 1)
    return {"color": palette[idx], "marker": markers[idx], "linestyle": "-"}


def alpha_label(alpha: float, meta: Dict[str, Any]) -> str:
    text = rf"$\alpha={float(alpha):g}$"
    if same(alpha, meta["reference_alpha"]):
        text += " (reference)"
    return text


def draw_s11_axis(ax: plt.Axes, panel: str, data: pd.DataFrame, meta: Dict[str, Any], show_letter: bool = True) -> None:
    prep_supp_ax(ax, panel if show_letter else None)
    ref = float(meta["reference_alpha"])
    if panel == "A":
        g = data.sort_values("alpha")
        ax.plot(g["alpha"], g["fraction_reference_tt_retained"],
                marker="^", color=cycle_colour(0), linewidth=LINE_WIDTHS["curve"],
                markersize=MARKER_SIZE)
        ax.axvline(ref, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.65)
        ax.set_xlabel(r"Operational observation threshold, $\alpha$")
        ax.set_ylabel("Fraction of reference TT retained")
        ax.set_xticks(sorted(g["alpha"].unique()))
        ax.set_ylim(0, 1.05)
    elif panel in ("B", "C", "D"):
        mode = {"B": "T0PLUS", "C": "TT", "D": "TF"}[panel]
        ax.axhline(0, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.6)
        for av in alpha_plot_order(meta):
            g = data.loc[np.isclose(data["alpha"].astype(float), av, rtol=0, atol=1e-12)].sort_values("bin")
            st = alpha_style(av, meta)
            ax.plot(g["x_center"], g["mean_dx"], color=st["color"], marker=st["marker"],
                    linestyle="-", linewidth=LINE_WIDTHS["curve"], markersize=MARKER_SIZE,
                    label=alpha_label(av, meta))
        ax.set_xlabel(r"Initial observed log-frequency, $x_0$")
        ax.set_ylabel(f"{mode} mean replicate-decoupled\ndisplacement")
        ax.legend(frameon=False, loc="best", handlelength=2)
    elif panel == "E":
        ax.axhline(0, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.6)
        ax.axvline(ref, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.65)
        mode_labels = {"T0PLUS": "T0PLUS", "TT": "TT", "TF": "TF"}
        for mode in ["T0PLUS", "TT", "TF"]:
            g = data.loc[data["mode"].eq(mode)].sort_values("alpha")
            st = style_for(mode)
            intervals(ax, g["alpha"], g["observed_slope_difference"],
                      g["bootstrap_slope_difference_q025"], g["bootstrap_slope_difference_q975"],
                      st["color"])
            ax.plot(g["alpha"], g["observed_slope_difference"], marker=st["marker"],
                    color=st["color"], linewidth=LINE_WIDTHS["curve"], markersize=MARKER_SIZE,
                    label=mode_labels[mode])
        ax.set_xlabel(r"Operational observation threshold, $\alpha$")
        ax.set_ylabel(r"Forward-slope difference vs $\alpha=0.05$")
        ax.set_xticks(meta["alphas"])
        ax.legend(frameon=False, loc="best", handlelength=2)
    elif panel == "F":
        ax.axhline(0, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.6)
        for av in alpha_plot_order(meta):
            g = data.loc[np.isclose(data["alpha"].astype(float), av, rtol=0, atol=1e-12)].sort_values("bin")
            st = alpha_style(av, meta)
            ax.plot(g["x_center"], g["cross_cov"], color=st["color"], marker=st["marker"],
                    linewidth=LINE_WIDTHS["curve"], markersize=MARKER_SIZE, label=alpha_label(av, meta))
        ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
        ax.set_ylabel("NON_FF cross-replicate covariance")
        ax.legend(frameon=False, loc="best", handlelength=2)
    elif panel == "G":
        ax.axhline(0, color=neutral(ax), linestyle="--", linewidth=LINE_WIDTHS["reference"], alpha=0.6)
        for av in [a for a in alpha_plot_order(meta) if not same(a, ref)]:
            g = data.loc[np.isclose(data["alpha"].astype(float), av, rtol=0, atol=1e-12)].sort_values("bin")
            st = alpha_style(av, meta)
            x = g["x_center"].to_numpy(float)
            y = g["difference"].to_numpy(float)
            lo = g["bootstrap_difference_q025"].to_numpy(float)
            hi = g["bootstrap_difference_q975"].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y)
            ci = ok & np.isfinite(lo) & np.isfinite(hi)
            if ci.any():
                ax.fill_between(x, lo, hi, where=ci, color=st["color"], alpha=0.12, linewidth=0)
            ax.plot(x, y, color=st["color"], marker=st["marker"], linewidth=LINE_WIDTHS["curve"],
                    markersize=MARKER_SIZE, label=alpha_label(av, meta))
        ax.set_xlabel(r"Transition-centred latent log-frequency, $x_{\mathrm{mid}}$")
        ax.set_ylabel(r"TT cross-covariance difference vs $\alpha=0.05$")
        ax.legend(frameon=False, loc="best", handlelength=2)
    elif panel == "H":
        for av in alpha_plot_order(meta):
            g = data.loc[np.isclose(data["alpha"].astype(float), av, rtol=0, atol=1e-12)].sort_values("dt")
            st = alpha_style(av, meta)
            ax.plot(g["dt"], g["cross_cov"], color=st["color"], marker=st["marker"],
                    linewidth=LINE_WIDTHS["curve"], markersize=MARKER_SIZE, label=alpha_label(av, meta))
        ax.set_xlabel("Temporal lag (weeks)")
        ax.set_ylabel("TT common-core cross-replicate covariance")
        ax.set_xticks(meta["dt_values"])
        ax.legend(frameon=False, loc="best", handlelength=2)
    else:
        raise ValueError(panel)


def save_supp_individual(series: str, panels: Dict[str, pd.DataFrame], out: Path,
                         formats: Sequence[str], dpi: int, meta: Optional[Dict[str, Any]] = None) -> List[str]:
    if series == "S9":
        stems, size = SUPP9_STEMS, SUPP_FIGURE_SIZES["S9_panel"]
    elif series == "S10":
        stems, size = SUPP10_STEMS, SUPP_FIGURE_SIZES["S10_panel"]
    elif series == "S11":
        stems, size = SUPP11_STEMS, SUPP_FIGURE_SIZES["S11_panel"]
        if meta is None:
            raise ValueError("S11 plotting requires metadata")
    else:
        raise ValueError(series)
    names: List[str] = []
    for panel, data in panels.items():
        fig, ax = plt.subplots(figsize=size)
        if series == "S9":
            draw_s9_axis(ax, panel, panels[panel], show_letter=SHOW_SUPP_INDIVIDUAL_PANEL_LABELS)
        elif series == "S10":
            draw_s10_axis(ax, panel, panels[panel], show_letter=SHOW_SUPP_INDIVIDUAL_PANEL_LABELS)
        else:
            draw_s11_axis(ax, panel, panels[panel], meta, show_letter=SHOW_SUPP_INDIVIDUAL_PANEL_LABELS)
        key = f"{series}{panel}"
        apply_axis_config(ax, SUPP_AXIS_CONFIG[key], key)
        apply_supp_legend(ax, key)
        fig.tight_layout(pad=0.9)
        for fmt in formats:
            name = stems[panel] + "." + fmt
            fig.savefig(out / name, dpi=dpi)
            names.append(name)
        plt.close(fig)
    return names


def save_supp9_composite(panels: Dict[str, pd.DataFrame], out: Path,
                         formats: Sequence[str], dpi: int) -> List[str]:
    fig, axes = plt.subplots(2, 2, figsize=SUPP_FIGURE_SIZES["S9_composite"])
    for ax, panel in zip(axes.flat, "ABCD"):
        draw_s9_axis(ax, panel, panels[panel], show_letter=SHOW_SUPP_COMPOSITE_PANEL_LABELS)
        key = f"S9{panel}"
        apply_axis_config(ax, SUPP_AXIS_CONFIG[key], key)
        apply_supp_legend(ax, key)
    fig.tight_layout(pad=1.15, w_pad=1.4, h_pad=1.6)
    names = []
    for fmt in formats:
        name = "SupplementaryFigure9." + fmt
        fig.savefig(out / name, dpi=dpi)
        names.append(name)
    plt.close(fig)
    return names


def save_supp10_composite(panels: Dict[str, pd.DataFrame], out: Path,
                          formats: Sequence[str], dpi: int) -> List[str]:
    fig = plt.figure(figsize=SUPP_FIGURE_SIZES["S10_composite"])
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1.02])
    axes = {
        "A": fig.add_subplot(gs[0, 0]), "B": fig.add_subplot(gs[0, 1]),
        "C": fig.add_subplot(gs[1, 0]), "D": fig.add_subplot(gs[1, 1]),
        "E": fig.add_subplot(gs[2, 0]), "F": fig.add_subplot(gs[2, 1]),
        "G": fig.add_subplot(gs[3, :]),
    }
    for panel in "ABCDEFG":
        draw_s10_axis(axes[panel], panel, panels[panel], show_letter=SHOW_SUPP_COMPOSITE_PANEL_LABELS)
        key = f"S10{panel}"
        apply_axis_config(axes[panel], SUPP_AXIS_CONFIG[key], key)
        apply_supp_legend(axes[panel], key)
    fig.tight_layout(pad=1.05, w_pad=1.4, h_pad=1.55)
    names = []
    for fmt in formats:
        name = "SupplementaryFigure10." + fmt
        fig.savefig(out / name, dpi=dpi)
        names.append(name)
    plt.close(fig)
    return names


def save_supp11_composite(panels: Dict[str, pd.DataFrame], meta: Dict[str, Any], out: Path,
                          formats: Sequence[str], dpi: int) -> List[str]:
    fig, axes = plt.subplots(4, 2, figsize=SUPP_FIGURE_SIZES["S11_composite"])
    for ax, panel in zip(axes.flat, "ABCDEFGH"):
        draw_s11_axis(ax, panel, panels[panel], meta, show_letter=SHOW_SUPP_COMPOSITE_PANEL_LABELS)
        key = f"S11{panel}"
        apply_axis_config(ax, SUPP_AXIS_CONFIG[key], key)
        apply_supp_legend(ax, key)
    fig.tight_layout(pad=1.05, w_pad=1.4, h_pad=1.55)
    names: List[str] = []
    for fmt in formats:
        name = "SupplementaryFigure11." + fmt
        fig.savefig(out / name, dpi=dpi)
        names.append(name)
    plt.close(fig)
    return names


def apply_axis_config(ax: plt.Axes, cfg: Dict[str, Any], label: str,
                      audit: Optional[Audit] = None) -> None:
    """Apply user-editable axis limits and tick positions/spacing.

    Explicit ``xticks``/``yticks`` take precedence over ``xtick_step``/
    ``ytick_step``.  This function only changes display geometry; it never
    alters, filters, interpolates or recomputes plotted values.
    """
    for coord in ["x", "y"]:
        bounds = cfg.get(coord + "lim")
        if bounds is not None:
            if len(bounds) != 2:
                raise InputError(f"{label}.{coord}lim must contain exactly two values")
            lo, hi = float(bounds[0]), float(bounds[1])
            if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
                raise InputError(f"{label}.{coord}lim must be a finite increasing pair")
            getattr(ax, "set_" + coord + "lim")(lo, hi)
            if audit is not None:
                audit.note(f"Custom {label} {coord}-limits",
                           f"({lo}, {hi}); manually verify that intended data/CI remain visible.")

        explicit = cfg.get(coord + "ticks")
        step = cfg.get(coord + "tick_step")
        if explicit is not None:
            vals = np.asarray(explicit, dtype=float)
            if vals.ndim != 1 or len(vals) == 0 or not np.isfinite(vals).all():
                raise InputError(f"{label}.{coord}ticks must be a nonempty sequence of finite values")
            if np.any(np.diff(vals) <= 0):
                raise InputError(f"{label}.{coord}ticks must be strictly increasing")
            getattr(ax, "set_" + coord + "ticks")(vals.tolist())
            if audit is not None:
                audit.note(f"Custom {label} {coord}-ticks", str(vals.tolist()))
        elif step is not None:
            step = float(step)
            if not math.isfinite(step) or step <= 0:
                raise InputError(f"{label}.{coord}tick_step must be a finite positive number")
            getattr(ax, coord + "axis").set_major_locator(MultipleLocator(step))
            if audit is not None:
                audit.note(f"Custom {label} {coord}-tick step", str(step))

    for axis in [ax.xaxis, ax.yaxis]:
        if isinstance(axis.get_major_formatter(), ScalarFormatter):
            axis.get_major_formatter().set_useOffset(False)


def finish_panel(fig, ax, panel: str, out: Path, formats: Sequence[str], dpi: int, audit: Audit) -> List[str]:
    apply_axis_config(ax, AXIS_CONFIG[panel], f"Figure8{panel}", audit)
    if SHOW_TITLES:
        ax.set_title(TITLES[panel], fontsize=FONT["axis"])
    fig.tight_layout(pad=0.9, rect=(0.015, 0.015, 0.995, 0.92 if SHOW_PANEL_LABELS else 0.99))
    if SHOW_PANEL_LABELS:
        fig.text(0.025, 0.975, panel, ha="left", va="top", fontweight="bold", fontsize=FONT["panel_label"])
    names = []
    for fmt in formats:
        name = STEMS[panel] + "." + fmt
        fig.savefig(out / name, dpi=dpi)  # fixed physical canvas; no bbox_inches='tight'
        names.append(name)
    plt.close(fig)
    return names


# =============================== EXPORT / CLI ===============================
def figure_note(meta: Dict[str, Any]) -> str:
    return f'''# Figure 8 | Interval and observation-domain sensitivity

(A) Stored Step-13 primary temporal slope and matched common-start/common-end
slopes from Step 14. Horizontal intervals are the supplied 95% subject-bootstrap
intervals. Labels report the corresponding subject counts. These are different
subject sets; the display is not a paired test of their slope differences.

(B) One-week observed replicate-decoupled forward profiles for PRIMARY, T0PLUS
and TT at alpha={meta['reference_alpha']:g}, using the Step-15 output coordinates.
The shaded region is the supplied PRIMARY 95% subject-bootstrap interval only.
TF is omitted from the main display but remains part of the source analysis.
Undefined points and missing bins are not interpolated.

(C) One-week cross-replicate covariance for PRIMARY, NON_FF and TT across
thresholds, read directly from Step-16 common-core temporal profiles. This is
a descriptive point-estimate display: amplitude confidence limits are not
exported by the inspected upstream profile schema and are not constructed here.

(D) Step-16 temporal slopes on the cross-threshold common core. Intervals are
the supplied 95% subject-bootstrap limits. The PRIMARY horizontal line/band
shows the threshold-invariant estimate/interval once; its four original rows
remain in the source CSV. Dotted horizontal reference: zero slope.

The empirical temporal core has {len(meta['step13_core_bins'])} bins. The Step-16
common core has {len(meta['step16_common_core_bins'])} bins. The fixed subject
universe is {', '.join(map(str, meta['complete_case_subjects']))}. C/D use the
Step-16 pooled-covariance/equal-bin sensitivity estimator, not the Step-13
within-subject equal-bin/equal-subject primary estimator used as A's reference.
T/F are operational observation labels, not biological presence/absence.

## Provenance

The program selects already-computed results; it does not estimate covariance,
refit slopes, reclassify states, recompute bootstrap intervals or rerun simulations.
`00_input_manifest.csv` records the exact input hashes; `00_plot_checks.csv`
records the checks and any unverified provenance limitations. `source_csv_row`
is one-based in the original CSV (header is row 1). Coordinates/stats are not
rounded in the exported CSVs. The optional Step-13 independent check was
{'performed' if meta['step13_independently_checked'] else 'not requested'}.

Figures are exported independently for final assembly. No Supplementary Excel
or manuscript text is changed by this plotting program.
'''



def supplementary_note(meta: Dict[str, Any]) -> str:
    return f'''# Supplementary Figures 9, 10 and 11

## Supplementary Figure 9 | Interval-position and anchored temporal controls

(A) BH-adjusted q values for the core-level RMS calendar-position statistic for
cross-replicate covariance, mean within-replicate variance and replicate-specific
excess at the evaluable lags (1–3 weeks). The dashed line marks q=0.05.

(B) Abundance-resolved BH-adjusted q values for the same RMS calendar-position
statistic applied to cross-replicate covariance at lags 1–3 weeks.

(C,D) Abundance-resolved temporal slopes for matched common-start and common-end
controls, respectively. Points are the stored slopes and vertical intervals are
the supplied 95% subject-bootstrap intervals. The dashed line marks zero slope.

## Supplementary Figure 10 | Fixed-threshold operational-domain sensitivity

All panels use the Step-15 reference operational threshold alpha={meta['reference_alpha']:g}.

(A) TT, TF, FT and FF fractions across temporal lags.

(B) Boundary-crossing imbalance (FT-TF)/(FT+TF) across abundance for lags 1–5.

(C) Differences between one-week forward profiles for T0PLUS, TT and PRIMARY.
Lines are stored point differences; shaded regions are supplied paired
subject-bootstrap intervals.

(D) One-week cross-replicate covariance differences for NON_FF and TT relative
to PRIMARY, with supplied paired subject-bootstrap intervals.

(E,F) One-week within-replicate variance and replicate-specific excess for
PRIMARY, NON_FF and TT. The PRIMARY shaded region is its supplied subject-bootstrap
interval; other curves are displayed as point estimates to limit overplotting.

(G) Common-core cross-replicate covariance over lags 1–5 for PRIMARY, NON_FF and
TT, read directly from the Step-15 temporal-domain profile table.

## Supplementary Figure 11 | Cross-threshold operational robustness

(A) Fraction of the reference-alpha TT population retained at each operational
threshold. The vertical dashed line marks the reference alpha.

(B–D) T0PLUS, TT and TF one-week replicate-decoupled forward profiles,
respectively, shown on the abundance support shared by all three modes and all
thresholds. Curves are stored point estimates; no rebinning or interpolation is
performed.

(E) Forward-slope differences relative to the reference alpha for T0PLUS, TT and
TF. Points are stored slope differences and vertical intervals are the supplied
paired subject-bootstrap limits.

(F) NON_FF one-week cross-replicate covariance profiles across operational
thresholds, restricted to bins satisfying the stored minimum-support flag.

(G) TT one-week cross-replicate covariance differences relative to the reference
alpha. Shaded regions are supplied paired bootstrap intervals where finite.

(H) TT common-core cross-replicate covariance over lags 1–5 across all operational
thresholds. The common temporal core contains {len(meta['step16_common_core_bins'])}
bins and the fixed subject universe contains {meta['n_subjects']} participants.

## Provenance

These figures select existing Step-14, Step-15 and Step-16 outputs. No permutation test,
bootstrap, covariance, variance, forward difference or temporal slope is recomputed.
The exported source CSVs preserve the original values and record the source table
and one-based original CSV row where applicable. Input hashes and all validation
checks are recorded with the main Figure 8 plotting outputs.
'''


def dump_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--step14-dir", type=Path, required=True)
    p.add_argument("--step15-dir", type=Path, required=True)
    p.add_argument("--step16-dir", type=Path, required=True)
    p.add_argument("--step13-dir", type=Path, default=None,
                   help="Recommended independent reference/core/cohort cross-check.")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--formats", nargs="+", choices=["pdf", "png", "svg"], default=["pdf", "png"])
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--font-family", default=None)
    p.add_argument("--validate-only", action="store_true", help="Validate and export source data/checks without figures.")
    p.add_argument("--overwrite", action="store_true", help="Replace only files owned by an earlier run of this plotter.")
    p.add_argument("--debug", action="store_true", help="Show traceback on error.")
    return p


def check_output(out: Path, roots: Sequence[Path], overwrite: bool) -> List[str]:
    for root in roots:
        if out == root or out in root.parents:
            raise InputError("Output must not be an analysis input folder or its ancestor.")
    if out.exists() and not out.is_dir():
        raise InputError(f"Output path is not a directory: {out}")
    if not out.exists() or not any(out.iterdir()):
        return []
    if not overwrite:
        raise InputError(f"Output folder is nonempty: {out}. Use a new folder, or --overwrite for this plotter's own files.")
    config = out / "00_plot_config.json"
    if not config.is_file():
        raise InputError("--overwrite requires an existing output owned by this plotter.")
    old = json.loads(config.read_text(encoding="utf-8"))
    if old.get("owner") != OWNER:
        raise InputError("Refusing to overwrite a folder created by another program.")
    names = old.get("managed_files", [])
    if not isinstance(names, list) or not all(isinstance(x, str) for x in names):
        raise InputError("Invalid managed-file list in previous plot config.")
    for name in names:
        q = Path(name)
        if q.is_absolute() or ".." in q.parts:
            raise InputError("Unsafe previous managed-file path")
    return names


def run(args: argparse.Namespace, audit: Audit) -> Path:
    if args.dpi <= 0:
        raise InputError("--dpi must be positive")
    roots = [getattr(args, f"step{s}_dir").expanduser().resolve(strict=True) for s in (14, 15, 16)]
    if any(not p.is_dir() for p in roots):
        raise InputError("All step inputs must be directories")
    step13 = args.step13_dir.expanduser().resolve(strict=True) if args.step13_dir else None
    all_roots = roots + ([step13] if step13 else [])
    out = args.outdir.expanduser().resolve()
    previous = check_output(out, all_roots, args.overwrite)
    panels, meta = load_panels(*roots, step13, audit)
    supp9, supp10, supp11 = load_supplementary(roots[0], roots[1], roots[2], meta, audit)
    # All data checks occur before creating the output/staging directory.
    family = configure_style(args.font_family, args.dpi)
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".figure8-staging-", dir=out.parent))
    try:
        (stage / "figure_source_data").mkdir()
        for panel, data in panels.items():
            data.to_csv(stage / "figure_source_data" / (STEMS[panel] + "_source_data.csv"), index=False)
        for panel, data in supp9.items():
            data.to_csv(stage / "figure_source_data" / (SUPP9_STEMS[panel] + "_source_data.csv"), index=False)
        for panel, data in supp10.items():
            data.to_csv(stage / "figure_source_data" / (SUPP10_STEMS[panel] + "_source_data.csv"), index=False)
        for panel, data in supp11.items():
            data.to_csv(stage / "figure_source_data" / (SUPP11_STEMS[panel] + "_source_data.csv"), index=False)
        if not args.validate_only:
            for panel, draw in zip("ABCD", [draw_a, draw_b, draw_c, draw_d]):
                fig, ax = draw(panels[panel], meta)
                finish_panel(fig, ax, panel, stage, args.formats, args.dpi, audit)
            save_supp_individual("S9", supp9, stage, args.formats, args.dpi, meta)
            save_supp9_composite(supp9, stage, args.formats, args.dpi)
            save_supp_individual("S10", supp10, stage, args.formats, args.dpi, meta)
            save_supp10_composite(supp10, stage, args.formats, args.dpi)
            save_supp_individual("S11", supp11, stage, args.formats, args.dpi, meta)
            save_supp11_composite(supp11, meta, stage, args.formats, args.dpi)
        (stage / "README_figure.md").write_text(figure_note(meta), encoding="utf-8")
        (stage / "README_supplementary.md").write_text(supplementary_note(meta), encoding="utf-8")
        # Detect a concurrent input change, without claiming full upstream validation.
        for item in audit.inputs:
            audit.check(f"Input unchanged: {item['filename']}", digest(Path(item["path"])) == item["sha256"])
        pd.DataFrame(audit.inputs).to_csv(stage / "00_input_manifest.csv", index=False)
        pd.DataFrame(audit.rows).to_csv(stage / "00_plot_checks.csv", index=False)
        records = [{"relative_path": str(p.relative_to(stage)), "size_bytes": p.stat().st_size,
                    "sha256": digest(p)} for p in sorted(stage.rglob("*")) if p.is_file()]
        pd.DataFrame(records).to_csv(stage / "figure_manifest.csv", index=False)
        managed = sorted(str(p.relative_to(stage)) for p in stage.rglob("*") if p.is_file()) + ["00_plot_config.json"]
        config = {"owner": OWNER, "script_version": VERSION, "script_sha256": digest(Path(__file__)),
                  "generated_utc": datetime.now(timezone.utc).isoformat(),
                  "command": shlex.join(sys.argv), "python": platform.python_version(),
                  "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": mpl.__version__,
                  "font_family_used": family, "formats": list(dict.fromkeys(args.formats)),
                  "validate_only": args.validate_only, "analysis_recomputed": False,
                  "metadata": meta, "settings": {"FONT": FONT, "FIGURE_SIZES": FIGURE_SIZES,
                  "SUPP_FIGURE_SIZES": SUPP_FIGURE_SIZES, "LINE_WIDTHS": LINE_WIDTHS,
                  "SERIES": SERIES, "SUPP9_METRICS": SUPP9_METRICS,
                  "SUPP10_CLASSES": SUPP10_CLASSES, "SUPP10_CONTRASTS": SUPP10_CONTRASTS,
                  "SUPP11_STEMS": SUPP11_STEMS, "AXIS_CONFIG": AXIS_CONFIG,
                  "SUPP_AXIS_CONFIG": SUPP_AXIS_CONFIG, "LEGEND_CONFIG": LEGEND_CONFIG,
                  "SUPP_LEGEND_CONFIG": SUPP_LEGEND_CONFIG,
                  "BOX_CONFIG": BOX_CONFIG, "SHOW_PANEL_LABELS": SHOW_PANEL_LABELS,
                  "SHOW_TITLES": SHOW_TITLES, "CI_ALPHA": CI_ALPHA},
                  "managed_files": managed}
        dump_json(stage / "00_plot_config.json", config)
        if out.exists():
            # Never replace unrelated contents. On a re-run, remove only stale
            # files explicitly registered by an earlier run of this program.
            for name in managed:
                target = out / name
                if target.exists() and name not in previous:
                    raise InputError(f"Output collision with an unowned file: {target}")
            for name in managed:
                dst = out / name
                dst.parent.mkdir(parents=True, exist_ok=True)
                os.replace(stage / name, dst)
            for name in set(previous) - set(managed):
                old_file = out / name
                if old_file.is_file():
                    old_file.unlink()
        else:
            os.replace(stage, out)
        return out
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parser().parse_args(argv)
    audit = Audit()
    try:
        out = run(args, audit)
    except Exception as exc:
        if args.debug:
            raise
        print(f"[ERROR] {type(exc).__name__}: {exc}", file=sys.stderr)
        print("No input analysis was modified. Correct the reported source/configuration issue before retrying.", file=sys.stderr)
        return 2
    passes = sum(r["status"] == "PASS" for r in audit.rows)
    print(f"[DONE] {'Validation/source export' if args.validate_only else 'Figure 8 + Supplementary Figures 9–11 plotting'}: {out}")
    print(f"[CHECKS] {passes} passed; no scientific statistics were recomputed.")
    print("[NOTE] Figure 8C contains point estimates only; its upstream table has no amplitude confidence limits.")
    if not args.validate_only:
        print("[OUTPUT] SupplementaryFigure9, SupplementaryFigure10 and SupplementaryFigure11 were written as assembled figures and separate panels.")
    for row in audit.rows:
        if row["status"] == "NOTE":
            print(f"[NOTE] {row['check']}: {row['detail']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
