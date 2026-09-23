#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
15-detectability_boundary_sensitivity.py
====================================

Canonical ClonoDynamics Step 15 after the replicate-resolved dynamics pivot.

SCIENTIFIC QUESTION
-------------------
Steps 9–14 deliberately define the manuscript-primary longitudinal dynamics
WITHOUT an operational TT filter:

    forward drift
        observed replicate-decoupled AB/BA estimator;

    shared fluctuations
        cross_cov = Cov(dx_observed_rep1, dx_observed_rep2
                        | xmid_latent, dt, common4);

    temporal scaling
        fixed abundance core, genuine longitudinal data.

Step 15 asks whether those conclusions depend materially on the fixed
operational observation-domain classification at the reference threshold
alpha = 0.05.

This is therefore an OBSERVATION-DOMAIN SENSITIVITY analysis.

It does NOT redefine the primary estimands and it does NOT interpret T/F as
biological presence/absence.

REFERENCE OPERATIONAL CLASSIFICATION
------------------------------------
The transition-level endpoint p-values carried by Step 5 are classified at the
fixed reference threshold:

    T0 = p_value_t0 < alpha
    T1 = p_value_t1 < alpha

giving:

    TT  : T0 and T1
    TF  : T0 and not T1
    FT  : not T0 and T1
    FF  : not T0 and not T1

The default alpha is 0.05.

Step 15 keeps alpha fixed.

Step 16 is the separate threshold-robustness analysis that varies alpha.

PRIMARY SENSITIVITY BRANCHES
----------------------------

1. Operational-class composition and boundary localization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 15 reports TT/TF/FT/FF composition by temporal lag and by the frozen
Step-11 xmid_latent abundance grid.

Boundary-crossing descriptors are descriptive:

    crossing fraction among NON_FF
        (TF + FT) / (TT + TF + FT)

    crossing imbalance
        (FT - TF) / (TF + FT)

These are not transition rates and do not imply biological appearance,
disappearance, birth or extinction.

2. Forward-drift observation-domain sensitivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The exact Step-9 observed replicate-decoupled geometry is retained:

    AB
        condition on x_observed_rep1_t0
        displacement = dx_observed_rep2
        support = forward_ab_eligible

    BA
        condition on x_observed_rep2_t0
        displacement = dx_observed_rep1
        support = forward_ba_eligible

The Step-9 fixed observed-abundance bin grid is reused.

Four operational domains are evaluated:

    PRIMARY
        the manuscript-primary Step-9 support, with no operational-class filter;

    T0PLUS
        operationally T at t0, t1 unrestricted = TT + TF;

    TT
        operationally T at both endpoints;

    TF
        operationally T at t0 and F at t1
        [descriptive component of T0PLUS].

AB and BA are combined with exact equal fold weight AFTER binning, exactly as
in final Step 9. A combined bin is valid only when BOTH folds satisfy the
support threshold; Step 15 never falls back to a single-fold estimate.

Biological uncertainty is estimated by one joint subject-cluster bootstrap,
using the same subject resample across folds and operational domains.

The primary sensitivity contrasts are:

    T0PLUS - PRIMARY
    TT - PRIMARY
    T0PLUS - TT

3. Shared-fluctuation observation-domain sensitivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The exact Step-11 displacement and conditioning variables are retained:

    dx_A = dx_observed_rep1
    dx_B = dx_observed_rep2
    conditioning = xmid_latent
    base support = common4

The Step-11 fixed xmid_latent bin grid is reused.

Three operational domains are evaluated within common4:

    PRIMARY
        all manuscript-primary common4 transitions;

    NON_FF
        TT + TF + FT;

    TT
        TT only.

For every subject x dt x abundance-bin cell, sufficient statistics are built
for:

    cross_cov
    same_var_mean
    replicate_specific_excess

and pooled Step-11-style bin estimates are obtained by aggregating those
sufficient statistics across biological subjects.

The same joint subject-cluster bootstrap is propagated across PRIMARY, NON_FF
and TT so that domain differences have paired bootstrap uncertainty. Within
each bootstrap draw, the reconstructed cell must still satisfy the declared
`min_n` and `min_subjects` support thresholds before contributing.

4. Temporal-scaling observation-domain sensitivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Temporal sensitivity is restricted to:

    the exact Step-13/14 complete-case subject universe;
    the exact Step-13/14 abundance core;
    dt = 1–5.

Because TT/NON_FF operational filtering can remove support from individual
subject-bin cells, Step 15 does NOT silently redefine the primary Step-13
equal-subject estimator.

Instead, it constructs a deliberately matched sensitivity estimand:

    - start from the frozen Step-13 abundance core;
    - retain only the largest contiguous bin subset having valid pooled
      covariance support for PRIMARY, NON_FF and TT at every requested lag;
    - use the same fixed subject universe from Step 14;
    - compute pooled covariance in each retained bin;
    - average retained bins equally within lag;
    - fit a signed temporal slope to each operational domain;
    - use the same subject-cluster bootstrap draws across all domains;
    - in every bootstrap draw require the SAME fixed common core to satisfy the
      support thresholds at every lag; otherwise that draw does not contribute
      a temporal slope.

This matched-core pooled-covariance sensitivity is reported alongside, but is
not substituted for, the Step-13 primary equal-bin/equal-subject result.

The key robustness question is not whether every estimator has the same exact
negative slope. It is whether positive temporal accumulation appears after
operational-domain restriction.

PIPELINE BOUNDARY
-----------------
Step 15 consumes:

    Step 5  longitudinal_transitions.parquet
    Step 9  forward bin grid and reference outputs
    Step 11 fluctuation bin grid and reference outputs
    Step 13 temporal reference
    Step 14 fixed-core / complete-subject contract

Step 15 does NOT:

    refit latent abundance;
    subtract the pseudo Step-12 technical null;
    change alpha;
    redefine Step-9 or Step-11 primary estimands;
    call operational T/F biological presence/absence.

OUTPUTS
-------
00_run_config.json
00_run_signature.json
00_analysis_summary.csv
00_pipeline_manifest.csv
00_input_manifest.csv
README_outputs.md

01_operational_class_composition_by_dt.csv
02_boundary_crossing_by_dt_bin.csv

03_forward_domain_by_bin.csv
04_forward_domain_summary.csv
05_forward_difference_vs_primary.csv
06_forward_primary_audit.csv

07_fluctuation_domain_by_bin_dt.csv
08_fluctuation_difference_vs_primary.csv
09_fluctuation_domain_summary.csv
10_fluctuation_primary_audit.csv

11_temporal_common_core.csv
12_temporal_domain_by_dt.csv
13_temporal_domain_slope_summary.csv
14_temporal_support_audit.csv

PRIMARY RUN
-----------
python3 ./code/dynamics_4/15-detectability_boundary_sensitivity.py \
    --transitions \
    ./results_4/5-longitudinal_transition_assembly/longitudinal_transitions.parquet \
    --step9-dir \
    ./results_4/9-observed_replicate_decoupled_forward_drift \
    --step11-dir \
    ./results_4/11-cross_replicate_fluctuation_dynamics \
    --step13-dir \
    ./results_4/13-temporal_fluctuation_scaling \
    --step14-dir \
    ./results_4/14-interval_position_structure \
    --outdir \
    ./results_4/15-detectability_boundary_sensitivity \
    --alpha 0.05 \
    --min-n 50 \
    --min-subjects 2 \
    --n-bootstrap 2000 \
    --seed 123
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

SCRIPT_VERSION = "v6-forward-safe-bootstrap-2026-09-18"

FORWARD_MODES = (
    "PRIMARY",
    "T0PLUS",
    "TT",
    "TF",
)

FLUCTUATION_MODES = (
    "PRIMARY",
    "NON_FF",
    "TT",
)

METRICS = (
    "cross_cov",
    "same_var_mean",
    "replicate_specific_excess",
)

OBS_CLASSES = (
    "TT",
    "TF",
    "FT",
    "FF",
)


# =============================================================================
# Generic helpers
# =============================================================================

def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return None if np.isnan(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def read_json(path: Path) -> Dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(path)
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def read_csv_required(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError(f"{label} is empty: {path}")
    return frame


def require_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
    label: str,
) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(
            f"{label}: missing columns {missing}; available={list(frame.columns)}"
        )


def finite_numeric(values) -> np.ndarray:
    x = pd.to_numeric(
        pd.Series(values),
        errors="coerce",
    ).to_numpy(float)
    return x[np.isfinite(x)]


def bool_mask_pandas(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).to_numpy(bool)
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "t", "yes", "y"])
        .to_numpy(bool)
    )


def parse_subject_string(text: object) -> List[int]:
    raw = str(text or "").strip()
    if not raw:
        return []
    return sorted(
        set(
            int(float(x.strip()))
            for x in raw.replace(",", ";").split(";")
            if x.strip()
        )
    )


def largest_contiguous_run(values: Iterable[int]) -> List[int]:
    vals = sorted(set(int(x) for x in values))
    if not vals:
        return []

    runs: List[List[int]] = []
    current = [vals[0]]

    for value in vals[1:]:
        if value == current[-1] + 1:
            current.append(value)
        else:
            runs.append(current)
            current = [value]

    runs.append(current)
    runs.sort(
        key=lambda r: (len(r), -r[0]),
        reverse=True,
    )
    return runs[0]


def safe_linear_fit(
    x_values: np.ndarray,
    y_values: np.ndarray,
) -> Dict[str, float]:
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)

    if int(ok.sum()) < 2:
        return {
            "n_points": int(ok.sum()),
            "intercept": np.nan,
            "slope": np.nan,
        }

    xx = x[ok]
    yy = y[ok]

    if np.ptp(xx) <= 0:
        return {
            "n_points": int(ok.sum()),
            "intercept": np.nan,
            "slope": np.nan,
        }

    slope, intercept = np.polyfit(
        xx,
        yy,
        1,
    )

    return {
        "n_points": int(ok.sum()),
        "intercept": float(intercept),
        "slope": float(slope),
    }


def percentile_summary(
    values: np.ndarray,
) -> Tuple[float, float, float, int]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]

    if x.size == 0:
        return np.nan, np.nan, np.nan, 0

    return (
        float(np.quantile(x, 0.025)),
        float(np.median(x)),
        float(np.quantile(x, 0.975)),
        int(x.size),
    )


def stable_sha256(payload: Mapping[str, object]) -> str:
    import hashlib
    raw = json.dumps(
        json_safe(dict(payload)),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def path_identity(path: Path) -> Dict[str, object]:
    p = path.expanduser().resolve(strict=True)
    st = p.stat()
    return {
        "path": str(p),
        "size_bytes": int(st.st_size),
        "mtime_ns": int(st.st_mtime_ns),
    }


def write_manifest(outdir: Path) -> None:
    rows = []
    for path in sorted(outdir.rglob("*")):
        if path.is_file() and path.name != "00_pipeline_manifest.csv":
            rows.append(
                {
                    "relative_path": str(path.relative_to(outdir)),
                    "suffix": path.suffix.lower(),
                    "size_bytes": int(path.stat().st_size),
                }
            )
    pd.DataFrame(rows).to_csv(
        outdir / "00_pipeline_manifest.csv",
        index=False,
    )


# =============================================================================
# Upstream validation / fixed grids
# =============================================================================

def validate_upstream(
    step9_config: Mapping[str, object],
    step9_signature: Mapping[str, object],
    step11_config: Mapping[str, object],
    step11_signature: Mapping[str, object],
    step13_config: Mapping[str, object],
    step13_signature: Mapping[str, object],
    step14_config: Mapping[str, object],
    step14_signature: Mapping[str, object],
) -> None:
    # Step 9 final exact AB/BA contract.
    step9_support = step9_config.get("primary_support", {})
    if not isinstance(step9_support, dict):
        step9_support = {}

    if step9_support.get("operational_TT_filter") is not False:
        raise ValueError(
            "Step 15 requires final Step 9 without operational TT filtering."
        )

    step9_est = step9_config.get("primary_estimand", {})
    if (
        step9_est.get("AB_BA_combination_policy")
        != "exact_equal_fold_weight_after_binning"
        or step9_est.get("AB_BA_row_count_weighting") is not False
    ):
        raise ValueError(
            "Step 15 requires final Step-9 exact equal-weight AB/BA contract."
        )

    s9_payload = step9_signature.get("signature_payload", {})
    if (
        s9_payload.get("AB_BA_combination_policy")
        != "exact_equal_fold_weight_after_binning"
        or s9_payload.get("AB_BA_row_count_weighting") is not False
    ):
        raise ValueError(
            "Step-9 run signature does not certify final AB/BA weighting."
        )

    # Step 11 final fluctuation contract.
    if str(step11_config.get("primary_support")) != "common4":
        raise ValueError("Step 15 requires Step 11 primary_support='common4'.")
    if str(step11_config.get("primary_conditioning")) != "xmid_latent":
        raise ValueError(
            "Step 15 requires Step 11 primary_conditioning='xmid_latent'."
        )
    if step11_config.get("operational_class_filter") is not None:
        raise ValueError(
            "Step 15 requires Step 11 without operational-class filtering."
        )

    s11_est = step11_config.get("primary_estimands", {})
    if s11_est.get("cross_cov_truncated_at_zero") is not False:
        raise ValueError(
            "Step 15 requires signed, untruncated Step-11 cross_cov."
        )

    s11_payload = step11_signature.get("signature_payload", {})
    if (
        s11_payload.get("primary_support") != "common4"
        or s11_payload.get("primary_conditioning") != "xmid_latent"
        or s11_payload.get("cross_cov_truncated_at_zero") is not False
    ):
        raise ValueError(
            "Step-11 run signature does not certify final fluctuation contract."
        )

    # Step 13 frozen biological cohort/core.
    if str(step13_config.get("primary_metric")) != "cross_cov":
        raise ValueError("Step 15 requires Step 13 primary_metric='cross_cov'.")
    if str(step13_config.get("primary_estimator")) != "equal_bin_equal_subject":
        raise ValueError(
            "Step 15 requires finalized Step-13 primary estimator."
        )
    complete = step13_config.get("complete_case", {})
    if (
        complete.get("defined_by_metric") != "cross_cov"
        or complete.get("technical_comparators_cannot_redefine_cohort") is not True
    ):
        raise ValueError(
            "Step 15 requires final Step-13 cross_cov-defined frozen cohort."
        )
    s13_payload = step13_signature.get("signature_payload", {})
    if s13_payload.get("primary_complete_case_defined_by") != "cross_cov":
        raise ValueError(
            "Step-13 run signature does not certify cross_cov-defined cohort."
        )

    # Step 14 final inherited contract.
    if str(step14_config.get("primary_support")) != "common4":
        raise ValueError("Step 15 requires Step 14 based on common4.")
    if str(step14_config.get("primary_conditioning")) != "xmid_latent":
        raise ValueError("Step 15 requires Step 14 based on xmid_latent.")
    if step14_config.get("step13_core_redefined") is not False:
        raise ValueError("Step 15 requires Step 14 to inherit Step-13 core.")
    if (
        step14_config.get("step13_complete_case_subject_universe_redefined")
        is not False
    ):
        raise ValueError(
            "Step 15 requires Step 14 to inherit Step-13 subject universe."
        )

    if (
        step14_config.get("analysis_signature")
        != step14_signature.get("analysis_signature")
    ):
        raise ValueError(
            "Step-14 run config/signature analysis signatures disagree."
        )


def normalize_edges(
    frame: pd.DataFrame,
    label: str,
) -> pd.DataFrame:
    require_columns(
        frame,
        ["bin", "x_left", "x_right", "x_center"],
        label,
    )

    d = frame[
        ["bin", "x_left", "x_right", "x_center"]
    ].drop_duplicates("bin").copy()

    for c in ["bin", "x_left", "x_right", "x_center"]:
        d[c] = pd.to_numeric(
            d[c],
            errors="raise",
        )

    d["bin"] = d["bin"].astype(int)

    d = d.sort_values("bin").reset_index(drop=True)

    if not np.all(
        np.diff(d["bin"].to_numpy(int)) == 1
    ):
        raise ValueError(
            f"{label}: bin IDs must be contiguous."
        )

    return d


def load_step14_contract(
    step14_dir: Path,
) -> Tuple[pd.DataFrame, List[int], List[int], List[int]]:
    contract = read_csv_required(
        step14_dir / "13_step15_contract.csv",
        "Step-14 Step-15 contract",
    )

    require_columns(
        contract,
        [
            "bin",
            "x_left",
            "x_right",
            "x_center",
            "complete_case_subjects",
            "n_complete_case_subjects",
            "dt_values",
            "primary_metric",
            "primary_support",
            "primary_conditioning",
            "step14_analysis_signature",
        ],
        "Step-14 Step-15 contract",
    )

    core = normalize_edges(
        contract[
            ["bin", "x_left", "x_right", "x_center"]
        ],
        "Step-14 core",
    )

    subject_values = (
        contract["complete_case_subjects"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    if len(subject_values) != 1:
        raise ValueError(
            "Inconsistent complete_case_subjects in Step-14 contract."
        )

    subjects = parse_subject_string(
        subject_values[0]
    )

    dt_values_raw = (
        contract["dt_values"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    if len(dt_values_raw) != 1:
        raise ValueError(
            "Inconsistent dt_values in Step-14 contract."
        )

    dts = parse_subject_string(
        dt_values_raw[0]
    )

    sig_values = (
        contract["step14_analysis_signature"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    if len(sig_values) != 1:
        raise ValueError(
            "Inconsistent step14_analysis_signature in Step-14 contract."
        )

    return (
        core,
        subjects,
        dts,
        sorted(core["bin"].astype(int).tolist()),
    )


# =============================================================================
# Transition-table helpers: pandas / PyArrow only
# =============================================================================

TRANSITION_COLUMNS = [
    "subject",
    "dt",
    "p_value_t0",
    "p_value_t1",
    "xmid_latent",
    "common4",
    "forward_ab_eligible",
    "forward_ba_eligible",
    "x_observed_rep1_t0",
    "x_observed_rep2_t0",
    "dx_observed_rep1",
    "dx_observed_rep2",
]


def read_transition_table(
    path: Path,
) -> pd.DataFrame:
    """
    Read only the Step-15 columns into pandas.

    Polars is intentionally not used in Step 15. On the target Apple-silicon
    workstation the Polars execution engine crashed natively during the
    forward group-by, including when streaming was disabled. The transition
    table contains only ~4.4 million rows and the machine has ample RAM, so a
    selected-column pandas/PyArrow load is the more robust execution path.
    """
    suffix = path.suffix.lower()

    if suffix in {".parquet", ".pq"}:
        try:
            frame = pd.read_parquet(
                path,
                columns=TRANSITION_COLUMNS,
                engine="pyarrow",
                use_threads=False,
            )
        except TypeError:
            frame = pd.read_parquet(
                path,
                columns=TRANSITION_COLUMNS,
                engine="pyarrow",
            )
    elif suffix == ".csv":
        frame = pd.read_csv(
            path,
            usecols=TRANSITION_COLUMNS,
        )
    else:
        raise ValueError(
            f"Unsupported transition table format: {path}"
        )

    validate_transition_schema(
        frame
    )

    # Normalize numeric columns once.
    numeric_cols = [
        "subject",
        "dt",
        "p_value_t0",
        "p_value_t1",
        "xmid_latent",
        "x_observed_rep1_t0",
        "x_observed_rep2_t0",
        "dx_observed_rep1",
        "dx_observed_rep2",
    ]

    for col in numeric_cols:
        frame[col] = pd.to_numeric(
            frame[col],
            errors="coerce",
        )

    frame = frame.dropna(
        subset=[
            "subject",
            "dt",
        ]
    ).copy()

    frame["subject"] = frame[
        "subject"
    ].astype(int)
    frame["dt"] = frame[
        "dt"
    ].astype(int)

    for col in [
        "common4",
        "forward_ab_eligible",
        "forward_ba_eligible",
    ]:
        frame[col] = bool_mask_pandas(
            frame[col]
        )

    return frame


def validate_transition_schema(
    frame: pd.DataFrame,
) -> None:
    missing = [
        c
        for c in TRANSITION_COLUMNS
        if c not in frame.columns
    ]

    if missing:
        raise ValueError(
            "Step-5 transition table is missing columns required by "
            f"Step 15: {missing}"
        )


def add_operational_columns(
    frame: pd.DataFrame,
    alpha: float,
) -> pd.DataFrame:
    """
    Add compact fixed-alpha operational state columns.

    obs_class_code:
        -1 = missing endpoint p-value
         0 = TT
         1 = TF
         2 = FT
         3 = FF
    """
    out = frame.copy()

    p0 = pd.to_numeric(
        out["p_value_t0"],
        errors="coerce",
    ).to_numpy(float)
    p1 = pd.to_numeric(
        out["p_value_t1"],
        errors="coerce",
    ).to_numpy(float)

    finite = (
        np.isfinite(p0)
        & np.isfinite(p1)
    )

    t0 = finite & (
        p0 < float(alpha)
    )
    t1 = finite & (
        p1 < float(alpha)
    )

    code = np.full(
        len(out),
        -1,
        dtype=np.int8,
    )
    code[
        finite & t0 & t1
    ] = 0
    code[
        finite & t0 & ~t1
    ] = 1
    code[
        finite & ~t0 & t1
    ] = 2
    code[
        finite & ~t0 & ~t1
    ] = 3

    out["op_t0"] = t0
    out["op_t1"] = t1
    out["obs_class_code"] = code

    return out


def assign_bins_numpy(
    values,
    edges: pd.DataFrame,
) -> np.ndarray:
    """
    Assign zero-based contiguous Step-bin IDs using the exact supplied edges.
    """
    x = pd.to_numeric(
        pd.Series(values),
        errors="coerce",
    ).to_numpy(float)

    rights = edges[
        "x_right"
    ].to_numpy(float)
    left0 = float(
        edges[
            "x_left"
        ].iloc[0]
    )
    last_right = float(
        rights[-1]
    )

    idx = np.searchsorted(
        rights,
        x,
        side="right",
    )

    # Include the exact upper edge in the last bin.
    at_last = (
        np.isfinite(x)
        & np.isclose(
            x,
            last_right,
            rtol=0.0,
            atol=1e-12,
        )
    )
    idx[
        at_last
    ] = len(
        rights
    ) - 1

    valid = (
        np.isfinite(x)
        & (x >= left0)
        & (x <= last_right)
        & (idx >= 0)
        & (idx < len(rights))
    )

    out = np.full(
        len(x),
        -1,
        dtype=np.int32,
    )
    out[
        valid
    ] = idx[
        valid
    ].astype(
        np.int32
    )

    return out



def calibrate_forward_assignment_edges(
    base_df: pd.DataFrame,
    reported_edges: pd.DataFrame,
    step9_primary: pd.DataFrame,
    forward_dt: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Recover assignment boundaries that reproduce the finalized Step-9
    fold-specific bin counts exactly.

    Why this is needed
    ------------------
    Step 9 was executed with Polars and its in-memory float boundaries.
    Step 15 v4 reads the saved CSV edges back into pandas. For a few highly
    discrete observed-frequency values lying numerically at quantile
    boundaries, CSV round-trip / backend comparison semantics can move an
    entire tied block into the adjacent bin.

    This does not change the scientific bin geometry. We therefore:
      1. retain the reported Step-9 x_left/x_right/x_center coordinates;
      2. derive an internal assignment threshold inside the open interval
         compatible with the finalized Step-9 AB and BA cumulative counts;
      3. require exact reproduction of every Step-9 fold-specific bin count.

    The returned assignment edges are used only to assign transitions. All
    reported coordinates remain the canonical Step-9 CSV coordinates.
    """
    require_columns(
        step9_primary,
        [
            "representation",
            "bin",
            "n_AB",
            "n_BA",
        ],
        "Step-9 primary forward curves",
    )

    ref = (
        step9_primary[
            step9_primary[
                "representation"
            ].astype(str).eq(
                "cross_combined"
            )
        ]
        .copy()
        .sort_values("bin")
        .reset_index(drop=True)
    )

    if len(ref) != len(reported_edges):
        raise ValueError(
            "Step-9 cross_combined row count does not match forward bin grid."
        )

    canonical = np.concatenate(
        [
            [
                float(
                    reported_edges[
                        "x_left"
                    ].iloc[0]
                )
            ],
            reported_edges[
                "x_right"
            ].to_numpy(float),
        ]
    )

    n_bins = len(
        reported_edges
    )

    fold_specs = [
        (
            "AB",
            "forward_ab_eligible",
            "x_observed_rep1_t0",
            "n_AB",
        ),
        (
            "BA",
            "forward_ba_eligible",
            "x_observed_rep2_t0",
            "n_BA",
        ),
    ]

    fold_intervals: Dict[
        str,
        List[
            Tuple[
                float,
                float,
            ]
        ],
    ] = {}

    fold_values: Dict[
        str,
        np.ndarray,
    ] = {}

    for (
        fold,
        eligibility,
        x_col,
        n_col,
    ) in fold_specs:
        x = pd.to_numeric(
            base_df[
                x_col
            ],
            errors="coerce",
        ).to_numpy(float)

        mask = (
            base_df[
                "dt"
            ].eq(
                int(
                    forward_dt
                )
            ).to_numpy(bool)
            & base_df[
                eligibility
            ].to_numpy(bool)
            & np.isfinite(x)
            & (
                x
                >= canonical[0]
            )
            & (
                x
                <= canonical[-1]
            )
        )

        values = np.sort(
            x[
                mask
            ]
        )

        expected = (
            pd.to_numeric(
                ref[
                    n_col
                ],
                errors="raise",
            )
            .astype(int)
            .to_numpy()
        )

        if int(
            expected.sum()
        ) != int(
            values.size
        ):
            raise ValueError(
                f"Cannot calibrate {fold}: in-range transition count "
                f"{values.size} != Step-9 expected {expected.sum()}."
            )

        fold_values[
            fold
        ] = values

        intervals = []

        cumulative = np.cumsum(
            expected
        )

        for k in range(
            1,
            n_bins,
        ):
            c = int(
                cumulative[
                    k - 1
                ]
            )

            if (
                c <= 0
                or c >= values.size
            ):
                raise ValueError(
                    f"Invalid cumulative Step-9 count at {fold} boundary {k}: {c}"
                )

            lower = float(
                values[
                    c - 1
                ]
            )
            upper = float(
                values[
                    c
                ]
            )

            if not (
                lower < upper
            ):
                raise ValueError(
                    f"Step-9 expected count would split an exact tie at "
                    f"{fold} boundary {k}: lower={lower}, upper={upper}."
                )

            # Any boundary e satisfying lower < e <= upper reproduces
            # the finalized Step-9 cumulative count for this fold.
            intervals.append(
                (
                    lower,
                    upper,
                )
            )

        fold_intervals[
            fold
        ] = intervals

    calibrated = canonical.copy()
    audit_rows = []

    for k in range(
        1,
        n_bins,
    ):
        lower = max(
            fold_intervals[
                "AB"
            ][
                k - 1
            ][0],
            fold_intervals[
                "BA"
            ][
                k - 1
            ][0],
        )

        upper = min(
            fold_intervals[
                "AB"
            ][
                k - 1
            ][1],
            fold_intervals[
                "BA"
            ][
                k - 1
            ][1],
        )

        if not (
            lower < upper
        ):
            raise ValueError(
                f"No common AB/BA assignment interval for forward boundary {k}: "
                f"({lower}, {upper}]"
            )

        original = float(
            canonical[
                k
            ]
        )

        if (
            lower
            < original
            <= upper
        ):
            chosen = original
        else:
            chosen = float(
                lower
                + 0.5
                * (
                    upper
                    - lower
                )
            )

            if not (
                lower
                < chosen
                <= upper
            ):
                chosen = float(
                    np.nextafter(
                        lower,
                        np.inf,
                    )
                )

            if not (
                lower
                < chosen
                <= upper
            ):
                chosen = upper

        calibrated[
            k
        ] = chosen

        audit_rows.append(
            {
                "boundary_index": int(
                    k
                ),
                "reported_boundary": original,
                "assignment_boundary": chosen,
                "difference": (
                    chosen
                    - original
                ),
                "allowed_lower_open": lower,
                "allowed_upper_closed": upper,
                "changed": bool(
                    chosen
                    != original
                ),
            }
        )

    assignment = reported_edges.copy()

    assignment[
        "x_left"
    ] = calibrated[:-1]
    assignment[
        "x_right"
    ] = calibrated[1:]
    assignment[
        "x_center"
    ] = 0.5 * (
        calibrated[:-1]
        + calibrated[1:]
    )

    # Hard validation: reproduce the finalized Step-9 AB/BA counts exactly.
    for (
        fold,
        eligibility,
        x_col,
        n_col,
    ) in fold_specs:
        x = pd.to_numeric(
            base_df[
                x_col
            ],
            errors="coerce",
        ).to_numpy(float)

        mask = (
            base_df[
                "dt"
            ].eq(
                int(
                    forward_dt
                )
            ).to_numpy(bool)
            & base_df[
                eligibility
            ].to_numpy(bool)
            & np.isfinite(x)
        )

        bins = assign_bins_numpy(
            x[
                mask
            ],
            assignment,
        )

        observed = np.bincount(
            bins[
                bins >= 0
            ],
            minlength=n_bins,
        ).astype(int)

        expected = (
            pd.to_numeric(
                ref[
                    n_col
                ],
                errors="raise",
            )
            .astype(int)
            .to_numpy()
        )

        if not np.array_equal(
            observed,
            expected,
        ):
            bad = np.flatnonzero(
                observed
                != expected
            )
            raise ValueError(
                f"Forward bin-membership calibration failed for {fold}; "
                f"mismatched bins={bad.tolist()}."
            )

    return (
        assignment,
        pd.DataFrame(
            audit_rows
        ),
    )


def operational_class_label(
    codes: pd.Series,
) -> pd.Series:
    mapping = {
        0: "TT",
        1: "TF",
        2: "FT",
        3: "FF",
        -1: "MISSING",
    }
    return codes.map(
        mapping
    )


# =============================================================================
# Operational composition / boundary descriptors
# =============================================================================

def build_operational_composition(
    base_df: pd.DataFrame,
    alpha: float,
    fluct_edges: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    classified = add_operational_columns(
        base_df,
        alpha,
    )

    overall_source = classified[
        classified[
            "obs_class_code"
        ].ge(0)
    ][
        [
            "dt",
            "subject",
            "obs_class_code",
        ]
    ].copy()

    overall = (
        overall_source.groupby(
            [
                "dt",
                "obs_class_code",
            ],
            as_index=False,
            sort=True,
        )
        .agg(
            n_transitions=(
                "subject",
                "size",
            ),
            n_subjects=(
                "subject",
                "nunique",
            ),
        )
    )

    overall[
        "obs_class_step15"
    ] = operational_class_label(
        overall[
            "obs_class_code"
        ]
    )

    totals = (
        overall.groupby(
            "dt",
            as_index=False,
        )[
            "n_transitions"
        ]
        .sum()
        .rename(
            columns={
                "n_transitions":
                    "n_transitions_dt_total"
            }
        )
    )

    overall = overall.merge(
        totals,
        on="dt",
        how="left",
        validate="many_to_one",
    )

    overall[
        "fraction_within_dt"
    ] = (
        overall[
            "n_transitions"
        ]
        / overall[
            "n_transitions_dt_total"
        ]
    )
    overall[
        "alpha"
    ] = float(alpha)

    # Abundance-resolved boundary localization.
    bin_ids = assign_bins_numpy(
        classified[
            "xmid_latent"
        ],
        fluct_edges,
    )

    binned = classified.loc[
        (
            bin_ids >= 0
        )
        & classified[
            "obs_class_code"
        ].ge(0),
        [
            "dt",
            "subject",
            "obs_class_code",
        ],
    ].copy()

    binned[
        "bin"
    ] = bin_ids[
        (
            bin_ids >= 0
        )
        & classified[
            "obs_class_code"
        ].ge(0).to_numpy(bool)
    ]

    grouped = (
        binned.groupby(
            [
                "dt",
                "bin",
                "obs_class_code",
            ],
            as_index=False,
            sort=True,
        )
        .agg(
            n_transitions=(
                "subject",
                "size",
            ),
            n_subjects=(
                "subject",
                "nunique",
            ),
        )
    )

    rows = []

    for (
        dt,
        bin_id,
    ), g in grouped.groupby(
        [
            "dt",
            "bin",
        ],
        sort=True,
    ):
        counts = {
            0: 0,
            1: 0,
            2: 0,
            3: 0,
        }

        for r in g.itertuples(
            index=False
        ):
            counts[
                int(
                    r.obs_class_code
                )
            ] = int(
                r.n_transitions
            )

        n_total = int(
            sum(
                counts.values()
            )
        )
        n_crossing = int(
            counts[1]
            + counts[2]
        )
        n_nonff = int(
            counts[0]
            + counts[1]
            + counts[2]
        )

        edge = fluct_edges[
            fluct_edges[
                "bin"
            ].eq(
                int(bin_id)
            )
        ].iloc[0]

        rows.append(
            {
                "dt": int(dt),
                "bin": int(bin_id),
                "x_left": float(
                    edge[
                        "x_left"
                    ]
                ),
                "x_right": float(
                    edge[
                        "x_right"
                    ]
                ),
                "x_center": float(
                    edge[
                        "x_center"
                    ]
                ),
                "n_total": n_total,
                "n_TT": counts[0],
                "n_TF": counts[1],
                "n_FT": counts[2],
                "n_FF": counts[3],
                "n_crossing": n_crossing,
                "n_NON_FF": n_nonff,
                "crossing_fraction_all": (
                    n_crossing
                    / n_total
                    if n_total
                    else np.nan
                ),
                "crossing_fraction_non_ff": (
                    n_crossing
                    / n_nonff
                    if n_nonff
                    else np.nan
                ),
                "crossing_imbalance": (
                    (
                        counts[2]
                        - counts[1]
                    )
                    / n_crossing
                    if n_crossing
                    else np.nan
                ),
                "alpha": float(alpha),
            }
        )

    overall = overall.drop(
        columns=[
            "obs_class_code",
        ]
    )

    return (
        overall.sort_values(
            [
                "dt",
                "obs_class_step15",
            ]
        ).reset_index(
            drop=True
        ),
        pd.DataFrame(
            rows
        ).sort_values(
            [
                "dt",
                "bin",
            ]
        ).reset_index(
            drop=True
        ),
    )


# =============================================================================
# Forward sensitivity
# =============================================================================

def forward_mode_mask_pandas(
    frame: pd.DataFrame,
    mode: str,
) -> np.ndarray:
    if mode == "PRIMARY":
        return np.ones(
            len(frame),
            dtype=bool,
        )
    if mode == "T0PLUS":
        return frame[
            "op_t0"
        ].to_numpy(bool)
    if mode == "TT":
        return (
            frame[
                "op_t0"
            ].to_numpy(bool)
            & frame[
                "op_t1"
            ].to_numpy(bool)
        )
    if mode == "TF":
        return (
            frame[
                "op_t0"
            ].to_numpy(bool)
            & ~frame[
                "op_t1"
            ].to_numpy(bool)
        )
    raise ValueError(
        mode
    )


def build_forward_sufficient(
    base_df: pd.DataFrame,
    alpha: float,
    forward_dt: int,
    forward_edges: pd.DataFrame,
) -> pd.DataFrame:
    classified = add_operational_columns(
        base_df,
        alpha,
    )

    fold_specs = [
        (
            "AB",
            "forward_ab_eligible",
            "x_observed_rep1_t0",
            "dx_observed_rep2",
        ),
        (
            "BA",
            "forward_ba_eligible",
            "x_observed_rep2_t0",
            "dx_observed_rep1",
        ),
    ]

    pieces = []

    for (
        fold,
        eligibility,
        x_col,
        dx_col,
    ) in fold_specs:
        mask = (
            classified[
                "dt"
            ].eq(
                int(
                    forward_dt
                )
            ).to_numpy(bool)
            & classified[
                eligibility
            ].to_numpy(bool)
        )

        x = pd.to_numeric(
            classified[
                x_col
            ],
            errors="coerce",
        ).to_numpy(float)

        dx = pd.to_numeric(
            classified[
                dx_col
            ],
            errors="coerce",
        ).to_numpy(float)

        mask &= (
            np.isfinite(x)
            & np.isfinite(dx)
        )

        work = classified.loc[
            mask,
            [
                "subject",
                "op_t0",
                "op_t1",
            ],
        ].copy()

        work[
            "x"
        ] = x[
            mask
        ]
        work[
            "dx"
        ] = dx[
            mask
        ]

        bins = assign_bins_numpy(
            work[
                "x"
            ],
            forward_edges,
        )

        work[
            "bin"
        ] = bins
        work = work[
            work[
                "bin"
            ].ge(0)
        ].copy()

        merged = None

        for mode in FORWARD_MODES:
            mm = forward_mode_mask_pandas(
                work,
                mode,
            )

            sub = work.loc[
                mm,
                [
                    "subject",
                    "bin",
                    "dx",
                ],
            ]

            if sub.empty:
                grouped = pd.DataFrame(
                    columns=[
                        "subject",
                        "bin",
                        f"n_{mode}",
                        f"sum_dx_{mode}",
                    ]
                )
            else:
                grouped = (
                    sub.groupby(
                        [
                            "subject",
                            "bin",
                        ],
                        as_index=False,
                        sort=False,
                    )
                    .agg(
                        **{
                            f"n_{mode}": (
                                "dx",
                                "size",
                            ),
                            f"sum_dx_{mode}": (
                                "dx",
                                "sum",
                            ),
                        }
                    )
                )

            if merged is None:
                merged = grouped
            else:
                merged = merged.merge(
                    grouped,
                    on=[
                        "subject",
                        "bin",
                    ],
                    how="outer",
                )

        assert merged is not None

        for mode in FORWARD_MODES:
            for col in [
                f"n_{mode}",
                f"sum_dx_{mode}",
            ]:
                merged[
                    col
                ] = pd.to_numeric(
                    merged[
                        col
                    ],
                    errors="coerce",
                ).fillna(
                    0.0
                )

        merged[
            "subject"
        ] = pd.to_numeric(
            merged[
                "subject"
            ],
            errors="raise",
        ).astype(int)
        merged[
            "bin"
        ] = pd.to_numeric(
            merged[
                "bin"
            ],
            errors="raise",
        ).astype(int)
        merged[
            "fold"
        ] = fold

        pieces.append(
            merged
        )

    return pd.concat(
        pieces,
        ignore_index=True,
    )


def build_forward_arrays(
    suff: pd.DataFrame,
    subjects: Sequence[int],
    n_bins: int,
) -> Tuple[
    Dict[Tuple[str, str], np.ndarray],
    Dict[Tuple[str, str], np.ndarray],
]:
    s_index = {
        int(s): i
        for i, s in enumerate(
            subjects
        )
    }

    counts = {}
    sums = {}

    for mode in FORWARD_MODES:
        for fold in ["AB", "BA"]:
            counts[
                (mode, fold)
            ] = np.zeros(
                (
                    len(subjects),
                    int(n_bins),
                ),
                dtype=float,
            )

            sums[
                (mode, fold)
            ] = np.zeros(
                (
                    len(subjects),
                    int(n_bins),
                ),
                dtype=float,
            )

    for row in suff.itertuples(
        index=False
    ):
        subject = int(
            row.subject
        )
        if subject not in s_index:
            continue

        i = s_index[
            subject
        ]
        b = int(
            row.bin
        )
        fold = str(
            row.fold
        )

        for mode in FORWARD_MODES:
            counts[
                (mode, fold)
            ][i, b] = float(
                getattr(
                    row,
                    f"n_{mode}",
                )
            )
            sums[
                (mode, fold)
            ][i, b] = float(
                getattr(
                    row,
                    f"sum_dx_{mode}",
                )
            )

    return counts, sums


def combine_fold_means(
    mean_ab: np.ndarray,
    mean_ba: np.ndarray,
    valid_ab: np.ndarray,
    valid_ba: np.ndarray,
) -> np.ndarray:
    """
    Final Step-9-compatible AB/BA combination.

    A combined bin exists only when BOTH folds satisfy support. There is no
    single-fold fallback.
    """
    out = np.full(
        mean_ab.shape,
        np.nan,
        dtype=float,
    )

    both = valid_ab & valid_ba
    out[both] = 0.5 * (
        mean_ab[both]
        + mean_ba[both]
    )
    return out


def forward_point_curves(
    counts: Mapping[
        Tuple[str, str],
        np.ndarray,
    ],
    sums: Mapping[
        Tuple[str, str],
        np.ndarray,
    ],
    forward_edges: pd.DataFrame,
    min_n: int,
) -> pd.DataFrame:
    rows = []

    for mode in FORWARD_MODES:
        n_ab = counts[
            (mode, "AB")
        ].sum(
            axis=0
        )
        n_ba = counts[
            (mode, "BA")
        ].sum(
            axis=0
        )

        sum_ab = sums[
            (mode, "AB")
        ].sum(
            axis=0
        )
        sum_ba = sums[
            (mode, "BA")
        ].sum(
            axis=0
        )

        mean_ab = np.divide(
            sum_ab,
            n_ab,
            out=np.full_like(
                sum_ab,
                np.nan,
            ),
            where=n_ab > 0,
        )
        mean_ba = np.divide(
            sum_ba,
            n_ba,
            out=np.full_like(
                sum_ba,
                np.nan,
            ),
            where=n_ba > 0,
        )

        valid_ab = n_ab >= int(
            min_n
        )
        valid_ba = n_ba >= int(
            min_n
        )

        combined = combine_fold_means(
            mean_ab,
            mean_ba,
            valid_ab,
            valid_ba,
        )

        for b, edge in forward_edges.iterrows():
            bin_id = int(
                edge["bin"]
            )
            rows.append(
                {
                    "mode": mode,
                    "bin": bin_id,
                    "x_left": float(
                        edge["x_left"]
                    ),
                    "x_right": float(
                        edge["x_right"]
                    ),
                    "x_center": float(
                        edge["x_center"]
                    ),
                    "n_AB": int(
                        n_ab[bin_id]
                    ),
                    "n_BA": int(
                        n_ba[bin_id]
                    ),
                    "valid_AB": bool(
                        valid_ab[bin_id]
                    ),
                    "valid_BA": bool(
                        valid_ba[bin_id]
                    ),
                    "mean_dx_AB": float(
                        mean_ab[bin_id]
                    ) if np.isfinite(
                        mean_ab[bin_id]
                    ) else np.nan,
                    "mean_dx_BA": float(
                        mean_ba[bin_id]
                    ) if np.isfinite(
                        mean_ba[bin_id]
                    ) else np.nan,
                    "mean_dx": float(
                        combined[bin_id]
                    ) if np.isfinite(
                        combined[bin_id]
                    ) else np.nan,
                }
            )

    return pd.DataFrame(rows)


def bootstrap_forward(
    counts: Mapping[
        Tuple[str, str],
        np.ndarray,
    ],
    sums: Mapping[
        Tuple[str, str],
        np.ndarray,
    ],
    subjects: Sequence[int],
    forward_edges: pd.DataFrame,
    min_n: int,
    n_bootstrap: int,
    seed: int,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    np.ndarray,
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
]:
    rng = np.random.default_rng(
        int(seed)
    )

    n_subjects = len(
        subjects
    )
    n_bins = len(
        forward_edges
    )

    draws = rng.integers(
        0,
        n_subjects,
        size=(
            int(n_bootstrap),
            n_subjects,
        ),
    )

    multiplicities = np.stack(
        [
            np.bincount(
                draw,
                minlength=n_subjects,
            )
            for draw in draws
        ],
        axis=0,
    ).astype(float)

    boot_curves: Dict[
        str,
        np.ndarray
    ] = {
        mode: np.full(
            (
                int(n_bootstrap),
                n_bins,
            ),
            np.nan,
            dtype=float,
        )
        for mode in FORWARD_MODES
    }

    def weighted_subject_sum(
        weights: np.ndarray,
        values: np.ndarray,
        label: str,
    ) -> np.ndarray:
        """
        Stable subject-bootstrap aggregation.

        This deliberately avoids BLAS-backed matrix multiplication. On some
        Apple-silicon NumPy/Accelerate builds, very small finite float64
        matrices can emit spurious divide/overflow/invalid warnings inside
        `matmul`. The explicit broadcasted weighted sum is tiny here
        (n_bootstrap x n_subjects x n_bins), deterministic, and numerically
        equivalent for this application.
        """
        w = np.asarray(weights, dtype=np.float64)
        v = np.asarray(values, dtype=np.float64)

        if w.ndim != 2 or v.ndim != 2:
            raise ValueError(
                f"{label}: expected 2D weights/values; "
                f"got {w.shape} and {v.shape}."
            )
        if w.shape[1] != v.shape[0]:
            raise ValueError(
                f"{label}: incompatible shapes {w.shape} and {v.shape}."
            )
        if not np.isfinite(w).all():
            raise ValueError(
                f"{label}: non-finite bootstrap multiplicities."
            )
        if not np.isfinite(v).all():
            bad = int(np.size(v) - np.isfinite(v).sum())
            raise ValueError(
                f"{label}: {bad} non-finite sufficient-statistic values "
                "before bootstrap aggregation."
            )

        out = np.sum(
            w[:, :, None] * v[None, :, :],
            axis=1,
            dtype=np.float64,
        )

        if not np.isfinite(out).all():
            bad = int(np.size(out) - np.isfinite(out).sum())
            raise FloatingPointError(
                f"{label}: bootstrap aggregation produced {bad} "
                "non-finite values."
            )

        return out

    for mode in FORWARD_MODES:
        n_ab = weighted_subject_sum(
            multiplicities,
            counts[(mode, "AB")],
            f"{mode}/AB counts",
        )
        n_ba = weighted_subject_sum(
            multiplicities,
            counts[(mode, "BA")],
            f"{mode}/BA counts",
        )

        sum_ab = weighted_subject_sum(
            multiplicities,
            sums[(mode, "AB")],
            f"{mode}/AB displacement sums",
        )
        sum_ba = weighted_subject_sum(
            multiplicities,
            sums[(mode, "BA")],
            f"{mode}/BA displacement sums",
        )

        mean_ab = np.divide(
            sum_ab,
            n_ab,
            out=np.full_like(
                sum_ab,
                np.nan,
            ),
            where=n_ab > 0,
        )
        mean_ba = np.divide(
            sum_ba,
            n_ba,
            out=np.full_like(
                sum_ba,
                np.nan,
            ),
            where=n_ba > 0,
        )

        valid_ab = (
            n_ab >= int(
                min_n
            )
        )
        valid_ba = (
            n_ba >= int(
                min_n
            )
        )

        both = valid_ab & valid_ba

        combined = np.full_like(
            mean_ab,
            np.nan,
        )

        combined[both] = 0.5 * (
            mean_ab[both]
            + mean_ba[both]
        )

        boot_curves[
            mode
        ] = combined

    ci_rows = []

    for mode in FORWARD_MODES:
        for _, edge in forward_edges.iterrows():
            bin_id = int(
                edge["bin"]
            )

            values = boot_curves[
                mode
            ][:, bin_id]

            lo, med, hi, n_valid = (
                percentile_summary(
                    values
                )
            )

            ci_rows.append(
                {
                    "mode": mode,
                    "bin": bin_id,
                    "bootstrap_median": med,
                    "bootstrap_q025": lo,
                    "bootstrap_q975": hi,
                    "n_bootstrap_valid": n_valid,
                }
            )

    x = forward_edges[
        "x_center"
    ].to_numpy(float)

    slope_rows = []
    boot_slope_arrays: Dict[
        str,
        np.ndarray
    ] = {}

    for mode in FORWARD_MODES:
        point_n_ab = counts[
            (mode, "AB")
        ].sum(
            axis=0
        )
        point_n_ba = counts[
            (mode, "BA")
        ].sum(
            axis=0
        )

        point_sum_ab = sums[
            (mode, "AB")
        ].sum(
            axis=0
        )
        point_sum_ba = sums[
            (mode, "BA")
        ].sum(
            axis=0
        )

        point_mean_ab = np.divide(
            point_sum_ab,
            point_n_ab,
            out=np.full_like(
                point_sum_ab,
                np.nan,
            ),
            where=point_n_ab > 0,
        )
        point_mean_ba = np.divide(
            point_sum_ba,
            point_n_ba,
            out=np.full_like(
                point_sum_ba,
                np.nan,
            ),
            where=point_n_ba > 0,
        )

        point_curve = combine_fold_means(
            point_mean_ab,
            point_mean_ba,
            point_n_ab >= int(
                min_n
            ),
            point_n_ba >= int(
                min_n
            ),
        )

        point_fit = safe_linear_fit(
            x,
            point_curve,
        )

        boot_slopes = np.full(
            int(n_bootstrap),
            np.nan,
            dtype=float,
        )

        for b in range(
            int(n_bootstrap)
        ):
            fit = safe_linear_fit(
                x,
                boot_curves[
                    mode
                ][b],
            )
            boot_slopes[b] = fit[
                "slope"
            ]

        boot_slope_arrays[
            mode
        ] = boot_slopes.copy()

        lo, med, hi, n_valid = (
            percentile_summary(
                boot_slopes
            )
        )

        slope_rows.append(
            {
                "mode": mode,
                "observed_slope": point_fit[
                    "slope"
                ],
                "observed_intercept": point_fit[
                    "intercept"
                ],
                "n_valid_bins": point_fit[
                    "n_points"
                ],
                "bootstrap_median_slope": med,
                "bootstrap_q025": lo,
                "bootstrap_q975": hi,
                "n_bootstrap_valid": n_valid,
                "bootstrap_positive_slope_fraction": (
                    float(
                        np.mean(
                            finite_numeric(
                                boot_slopes
                            )
                            > 0
                        )
                    )
                    if n_valid
                    else np.nan
                ),
            }
        )

    slope_summary = pd.DataFrame(
        slope_rows
    )

    primary_boot = boot_slope_arrays[
        "PRIMARY"
    ]
    primary_observed = float(
        slope_summary.loc[
            slope_summary[
                "mode"
            ].eq(
                "PRIMARY"
            ),
            "observed_slope",
        ].iloc[0]
    )

    for mode in [
        "T0PLUS",
        "TT",
        "TF",
    ]:
        diff = (
            boot_slope_arrays[
                mode
            ]
            - primary_boot
        )

        lo, med, hi, n_valid = (
            percentile_summary(
                diff
            )
        )

        mask = slope_summary[
            "mode"
        ].eq(
            mode
        )

        observed_mode = float(
            slope_summary.loc[
                mask,
                "observed_slope",
            ].iloc[0]
        )

        slope_summary.loc[
            mask,
            "observed_slope_difference_vs_PRIMARY",
        ] = (
            observed_mode
            - primary_observed
        )
        slope_summary.loc[
            mask,
            "bootstrap_slope_difference_vs_PRIMARY_median",
        ] = med
        slope_summary.loc[
            mask,
            "bootstrap_slope_difference_vs_PRIMARY_q025",
        ] = lo
        slope_summary.loc[
            mask,
            "bootstrap_slope_difference_vs_PRIMARY_q975",
        ] = hi
        slope_summary.loc[
            mask,
            "n_bootstrap_slope_difference_valid",
        ] = n_valid

    return (
        pd.DataFrame(ci_rows),
        slope_summary,
        draws,
        boot_curves,
        boot_slope_arrays,
    )


def build_forward_differences(
    curves: pd.DataFrame,
    boot_curves: Mapping[
        str,
        np.ndarray,
    ],
) -> pd.DataFrame:
    rows = []

    contrast_specs = [
        (
            "T0PLUS",
            "PRIMARY",
            "T0PLUS-PRIMARY",
        ),
        (
            "TT",
            "PRIMARY",
            "TT-PRIMARY",
        ),
        (
            "T0PLUS",
            "TT",
            "T0PLUS-TT",
        ),
    ]

    lookup = curves.set_index(
        [
            "mode",
            "bin",
        ]
    )

    bins = sorted(
        curves[
            "bin"
        ].astype(int).unique()
    )

    for (
        mode_a,
        mode_b,
        contrast,
    ) in contrast_specs:
        for bin_id in bins:
            key_a = (
                mode_a,
                int(bin_id),
            )
            key_b = (
                mode_b,
                int(bin_id),
            )

            if (
                key_a not in lookup.index
                or key_b not in lookup.index
            ):
                continue

            a = lookup.loc[
                key_a
            ]
            b = lookup.loc[
                key_b
            ]

            diff_boot = (
                boot_curves[
                    mode_a
                ][:, int(bin_id)]
                - boot_curves[
                    mode_b
                ][:, int(bin_id)]
            )

            lo, med, hi, n_valid = (
                percentile_summary(
                    diff_boot
                )
            )

            value_a = (
                float(
                    a[
                        "mean_dx"
                    ]
                )
                if pd.notna(
                    a[
                        "mean_dx"
                    ]
                )
                else np.nan
            )
            value_b = (
                float(
                    b[
                        "mean_dx"
                    ]
                )
                if pd.notna(
                    b[
                        "mean_dx"
                    ]
                )
                else np.nan
            )

            rows.append(
                {
                    "contrast": contrast,
                    "mode_a": mode_a,
                    "mode_b": mode_b,
                    "bin": int(
                        bin_id
                    ),
                    "x_left": float(
                        a[
                            "x_left"
                        ]
                    ),
                    "x_right": float(
                        a[
                            "x_right"
                        ]
                    ),
                    "x_center": float(
                        a[
                            "x_center"
                        ]
                    ),
                    "value_a": value_a,
                    "value_b": value_b,
                    "difference": (
                        value_a
                        - value_b
                        if np.isfinite(
                            value_a
                        )
                        and np.isfinite(
                            value_b
                        )
                        else np.nan
                    ),
                    "bootstrap_difference_median": med,
                    "bootstrap_difference_q025": lo,
                    "bootstrap_difference_q975": hi,
                    "n_bootstrap_difference_valid": n_valid,
                }
            )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# Fluctuation sufficient statistics
# =============================================================================

def fluctuation_mode_mask_pandas(
    frame: pd.DataFrame,
    mode: str,
) -> np.ndarray:
    if mode == "PRIMARY":
        return np.ones(
            len(frame),
            dtype=bool,
        )
    if mode == "NON_FF":
        return (
            frame[
                "obs_class_code"
            ].to_numpy(
                np.int8
            )
            >= 0
        ) & (
            frame[
                "obs_class_code"
            ].to_numpy(
                np.int8
            )
            != 3
        )
    if mode == "TT":
        return (
            frame[
                "obs_class_code"
            ].to_numpy(
                np.int8
            )
            == 0
        )
    raise ValueError(
        mode
    )


def build_fluctuation_sufficient(
    base_df: pd.DataFrame,
    alpha: float,
    fluct_edges: pd.DataFrame,
) -> pd.DataFrame:
    classified = add_operational_columns(
        base_df,
        alpha,
    )

    mask = classified[
        "common4"
    ].to_numpy(bool)

    x = pd.to_numeric(
        classified[
            "xmid_latent"
        ],
        errors="coerce",
    ).to_numpy(float)
    dx_a = pd.to_numeric(
        classified[
            "dx_observed_rep1"
        ],
        errors="coerce",
    ).to_numpy(float)
    dx_b = pd.to_numeric(
        classified[
            "dx_observed_rep2"
        ],
        errors="coerce",
    ).to_numpy(float)

    mask &= (
        np.isfinite(x)
        & np.isfinite(
            dx_a
        )
        & np.isfinite(
            dx_b
        )
    )

    work = classified.loc[
        mask,
        [
            "subject",
            "dt",
            "obs_class_code",
        ],
    ].copy()

    work[
        "xmid_latent"
    ] = x[
        mask
    ]
    work[
        "dxA"
    ] = dx_a[
        mask
    ]
    work[
        "dxB"
    ] = dx_b[
        mask
    ]

    work[
        "bin"
    ] = assign_bins_numpy(
        work[
            "xmid_latent"
        ],
        fluct_edges,
    )

    work = work[
        work[
            "bin"
        ].ge(0)
    ].copy()

    merged = None

    for mode in FLUCTUATION_MODES:
        mm = fluctuation_mode_mask_pandas(
            work,
            mode,
        )

        sub = work.loc[
            mm,
            [
                "subject",
                "dt",
                "bin",
                "dxA",
                "dxB",
            ],
        ].copy()

        if sub.empty:
            grouped = pd.DataFrame(
                columns=[
                    "subject",
                    "dt",
                    "bin",
                    f"n_{mode}",
                    f"sumA_{mode}",
                    f"sumB_{mode}",
                    f"sumAB_{mode}",
                    f"sumA2_{mode}",
                    f"sumB2_{mode}",
                ]
            )
        else:
            a = sub[
                "dxA"
            ].to_numpy(float)
            b = sub[
                "dxB"
            ].to_numpy(float)

            sub[
                "_AB"
            ] = a * b
            sub[
                "_A2"
            ] = a * a
            sub[
                "_B2"
            ] = b * b

            grouped = (
                sub.groupby(
                    [
                        "subject",
                        "dt",
                        "bin",
                    ],
                    as_index=False,
                    sort=False,
                )
                .agg(
                    **{
                        f"n_{mode}": (
                            "dxA",
                            "size",
                        ),
                        f"sumA_{mode}": (
                            "dxA",
                            "sum",
                        ),
                        f"sumB_{mode}": (
                            "dxB",
                            "sum",
                        ),
                        f"sumAB_{mode}": (
                            "_AB",
                            "sum",
                        ),
                        f"sumA2_{mode}": (
                            "_A2",
                            "sum",
                        ),
                        f"sumB2_{mode}": (
                            "_B2",
                            "sum",
                        ),
                    }
                )
            )

        if merged is None:
            merged = grouped
        else:
            merged = merged.merge(
                grouped,
                on=[
                    "subject",
                    "dt",
                    "bin",
                ],
                how="outer",
            )

    assert merged is not None

    for mode in FLUCTUATION_MODES:
        for name in [
            "n",
            "sumA",
            "sumB",
            "sumAB",
            "sumA2",
            "sumB2",
        ]:
            col = f"{name}_{mode}"
            merged[
                col
            ] = pd.to_numeric(
                merged[
                    col
                ],
                errors="coerce",
            ).fillna(
                0.0
            )

    for col in [
        "subject",
        "dt",
        "bin",
    ]:
        merged[
            col
        ] = pd.to_numeric(
            merged[
                col
            ],
            errors="raise",
        ).astype(int)

    return merged


def covariance_metrics_from_sufficient(
    n: np.ndarray,
    sum_a: np.ndarray,
    sum_b: np.ndarray,
    sum_ab: np.ndarray,
    sum_a2: np.ndarray,
    sum_b2: np.ndarray,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    n = np.asarray(
        n,
        dtype=float,
    )

    with np.errstate(
        divide="ignore",
        invalid="ignore",
    ):
        centered_ab = (
            sum_ab
            - (
                sum_a
                * sum_b
                / n
            )
        )
        centered_a2 = (
            sum_a2
            - (
                sum_a ** 2
                / n
            )
        )
        centered_b2 = (
            sum_b2
            - (
                sum_b ** 2
                / n
            )
        )

        denom = (
            n - 1.0
        )

        cross_cov = np.divide(
            centered_ab,
            denom,
            out=np.full_like(
                centered_ab,
                np.nan,
                dtype=float,
            ),
            where=denom > 0,
        )

        var_a = np.divide(
            centered_a2,
            denom,
            out=np.full_like(
                centered_a2,
                np.nan,
                dtype=float,
            ),
            where=denom > 0,
        )

        var_b = np.divide(
            centered_b2,
            denom,
            out=np.full_like(
                centered_b2,
                np.nan,
                dtype=float,
            ),
            where=denom > 0,
        )

    same_var_mean = 0.5 * (
        var_a + var_b
    )

    excess = (
        same_var_mean
        - cross_cov
    )

    return (
        cross_cov,
        same_var_mean,
        excess,
    )


def build_fluctuation_arrays(
    suff: pd.DataFrame,
    subjects: Sequence[int],
    dts: Sequence[int],
    n_bins: int,
) -> Dict[
    str,
    Dict[str, np.ndarray],
]:
    s_index = {
        int(s): i
        for i, s in enumerate(
            subjects
        )
    }
    d_index = {
        int(dt): j
        for j, dt in enumerate(
            dts
        )
    }

    arrays = {}

    for mode in FLUCTUATION_MODES:
        arrays[mode] = {
            name: np.zeros(
                (
                    len(subjects),
                    len(dts),
                    int(n_bins),
                ),
                dtype=float,
            )
            for name in [
                "n",
                "sumA",
                "sumB",
                "sumAB",
                "sumA2",
                "sumB2",
            ]
        }

    for row in suff.itertuples(
        index=False
    ):
        subject = int(
            row.subject
        )
        dt = int(
            row.dt
        )
        b = int(
            row.bin
        )

        if (
            subject not in s_index
            or dt not in d_index
        ):
            continue

        i = s_index[
            subject
        ]
        j = d_index[
            dt
        ]

        for mode in FLUCTUATION_MODES:
            for name in [
                "n",
                "sumA",
                "sumB",
                "sumAB",
                "sumA2",
                "sumB2",
            ]:
                arrays[mode][name][
                    i,
                    j,
                    b,
                ] = float(
                    getattr(
                        row,
                        f"{name}_{mode}",
                    )
                )

    return arrays


def aggregate_fluctuation(
    arrays: Mapping[
        str,
        np.ndarray,
    ],
    multiplicity: np.ndarray,
) -> Dict[str, np.ndarray]:
    aggregated = {
        name: np.tensordot(
            multiplicity,
            arrays[name],
            axes=(0, 0),
        )
        for name in arrays
    }

    cross, same, excess = (
        covariance_metrics_from_sufficient(
            aggregated["n"],
            aggregated["sumA"],
            aggregated["sumB"],
            aggregated["sumAB"],
            aggregated["sumA2"],
            aggregated["sumB2"],
        )
    )

    return {
        **aggregated,
        "cross_cov": cross,
        "same_var_mean": same,
        "replicate_specific_excess": excess,
    }


def fluctuation_point_table(
    arrays: Mapping[
        str,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    subjects: Sequence[int],
    dts: Sequence[int],
    fluct_edges: pd.DataFrame,
    min_n: int,
    min_subjects: int,
) -> pd.DataFrame:
    rows = []

    ones = np.ones(
        len(subjects),
        dtype=float,
    )

    for mode in FLUCTUATION_MODES:
        agg = aggregate_fluctuation(
            arrays[mode],
            ones,
        )

        n_subjects = (
            arrays[mode]["n"]
            > 0
        ).sum(
            axis=0
        )

        for j, dt in enumerate(
            dts
        ):
            for _, edge in fluct_edges.iterrows():
                b = int(
                    edge["bin"]
                )

                n = float(
                    agg["n"][
                        j,
                        b,
                    ]
                )
                ns = int(
                    n_subjects[
                        j,
                        b,
                    ]
                )

                valid = (
                    n >= int(
                        min_n
                    )
                    and ns >= int(
                        min_subjects
                    )
                )

                rows.append(
                    {
                        "mode": mode,
                        "dt": int(dt),
                        "bin": b,
                        "x_left": float(
                            edge["x_left"]
                        ),
                        "x_right": float(
                            edge["x_right"]
                        ),
                        "x_center": float(
                            edge["x_center"]
                        ),
                        "n": int(
                            n
                        ),
                        "n_subjects": ns,
                        "meets_min_support": bool(
                            valid
                        ),
                        "cross_cov": (
                            float(
                                agg[
                                    "cross_cov"
                                ][
                                    j,
                                    b,
                                ]
                            )
                            if valid
                            else np.nan
                        ),
                        "same_var_mean": (
                            float(
                                agg[
                                    "same_var_mean"
                                ][
                                    j,
                                    b,
                                ]
                            )
                            if valid
                            else np.nan
                        ),
                        "replicate_specific_excess": (
                            float(
                                agg[
                                    "replicate_specific_excess"
                                ][
                                    j,
                                    b,
                                ]
                            )
                            if valid
                            else np.nan
                        ),
                    }
                )

    return pd.DataFrame(rows)


def bootstrap_fluctuation(
    arrays: Mapping[
        str,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    subjects: Sequence[int],
    dts: Sequence[int],
    fluct_edges: pd.DataFrame,
    min_n: int,
    min_subjects: int,
    n_bootstrap: int,
    seed: int,
) -> Tuple[
    pd.DataFrame,
    Dict[str, Dict[str, np.ndarray]],
    np.ndarray,
]:
    rng = np.random.default_rng(
        int(seed)
    )

    n_subjects_total = len(
        subjects
    )

    draws = rng.integers(
        0,
        n_subjects_total,
        size=(
            int(n_bootstrap),
            n_subjects_total,
        ),
    )

    boot = {
        mode: {
            metric: np.full(
                (
                    int(n_bootstrap),
                    len(dts),
                    len(fluct_edges),
                ),
                np.nan,
                dtype=float,
            )
            for metric in METRICS
        }
        for mode in FLUCTUATION_MODES
    }

    point_support = {}

    for mode in FLUCTUATION_MODES:
        n_point = arrays[mode]["n"].sum(
            axis=0
        )
        ns_point = (
            arrays[mode]["n"]
            > 0
        ).sum(
            axis=0
        )
        point_support[mode] = (
            (n_point >= int(min_n))
            & (
                ns_point
                >= int(min_subjects)
            )
        )

    for b in range(
        int(n_bootstrap)
    ):
        multiplicity = np.bincount(
            draws[b],
            minlength=n_subjects_total,
        ).astype(float)

        for mode in FLUCTUATION_MODES:
            agg = aggregate_fluctuation(
                arrays[mode],
                multiplicity,
            )

            sampled_subject_support = (
                (arrays[mode]["n"] > 0)
                & (multiplicity[:, None, None] > 0)
            ).sum(axis=0)

            draw_support = (
                point_support[mode]
                & (agg["n"] >= int(min_n))
                & (
                    sampled_subject_support
                    >= int(min_subjects)
                )
            )

            for metric in METRICS:
                values = agg[
                    metric
                ].copy()
                values[
                    ~draw_support
                ] = np.nan
                boot[
                    mode
                ][
                    metric
                ][b] = values

    rows = []

    for mode in FLUCTUATION_MODES:
        for j, dt in enumerate(
            dts
        ):
            for _, edge in fluct_edges.iterrows():
                b = int(
                    edge["bin"]
                )

                for metric in METRICS:
                    lo, med, hi, n_valid = (
                        percentile_summary(
                            boot[
                                mode
                            ][
                                metric
                            ][
                                :,
                                j,
                                b,
                            ]
                        )
                    )

                    rows.append(
                        {
                            "mode": mode,
                            "metric": metric,
                            "dt": int(dt),
                            "bin": b,
                            "bootstrap_median": med,
                            "bootstrap_q025": lo,
                            "bootstrap_q975": hi,
                            "n_bootstrap_valid": n_valid,
                        }
                    )

    return (
        pd.DataFrame(rows),
        boot,
        draws,
    )


def wide_fluctuation_bootstrap_ci(
    long_ci: pd.DataFrame,
) -> pd.DataFrame:
    if long_ci.empty:
        return long_ci

    index_cols = [
        "mode",
        "dt",
        "bin",
    ]

    pieces = []

    for metric in METRICS:
        g = long_ci[
            long_ci[
                "metric"
            ].astype(str).eq(
                metric
            )
        ][
            index_cols
            + [
                "bootstrap_median",
                "bootstrap_q025",
                "bootstrap_q975",
                "n_bootstrap_valid",
            ]
        ].copy()

        g = g.rename(
            columns={
                "bootstrap_median":
                    f"{metric}_bootstrap_median",
                "bootstrap_q025":
                    f"{metric}_ci025",
                "bootstrap_q975":
                    f"{metric}_ci975",
                "n_bootstrap_valid":
                    f"{metric}_n_bootstrap_valid",
            }
        )

        pieces.append(
            g
        )

    out = pieces[0]

    for piece in pieces[1:]:
        out = out.merge(
            piece,
            on=index_cols,
            how="outer",
            validate="one_to_one",
        )

    return out


# =============================================================================
# Temporal matched-domain sensitivity
# =============================================================================

def temporal_support_common_core(
    arrays: Mapping[
        str,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    all_subjects: Sequence[int],
    temporal_subjects: Sequence[int],
    dts: Sequence[int],
    frozen_core_bins: Sequence[int],
    min_n: int,
    min_subjects: int,
) -> Tuple[
    List[int],
    pd.DataFrame,
]:
    subject_index = {
        int(s): i
        for i, s in enumerate(
            all_subjects
        )
    }

    idx = [
        subject_index[
            int(s)
        ]
        for s in temporal_subjects
        if int(s) in subject_index
    ]

    if len(idx) != len(
        temporal_subjects
    ):
        raise ValueError(
            "Step-14 temporal subject universe is not fully present in transitions."
        )

    rows = []
    qualifying = []

    for b in frozen_core_bins:
        bin_ok = True

        row: Dict[str, object] = {
            "bin": int(b),
        }

        for mode in FLUCTUATION_MODES:
            mode_arrays = {
                k: v[
                    idx,
                    :,
                    :,
                ]
                for k, v in arrays[
                    mode
                ].items()
            }

            n_total = mode_arrays[
                "n"
            ].sum(
                axis=0
            )
            n_subjects = (
                mode_arrays[
                    "n"
                ]
                > 0
            ).sum(
                axis=0
            )

            for j, dt in enumerate(
                dts
            ):
                n = int(
                    n_total[
                        j,
                        int(b),
                    ]
                )
                ns = int(
                    n_subjects[
                        j,
                        int(b),
                    ]
                )

                valid = (
                    n >= int(
                        min_n
                    )
                    and ns >= int(
                        min_subjects
                    )
                )

                row[
                    f"{mode}_n_dt{dt}"
                ] = n
                row[
                    f"{mode}_n_subjects_dt{dt}"
                ] = ns
                row[
                    f"{mode}_valid_dt{dt}"
                ] = bool(
                    valid
                )

                bin_ok = (
                    bin_ok
                    and valid
                )

        row[
            "valid_all_modes_all_dt"
        ] = bool(
            bin_ok
        )

        if bin_ok:
            qualifying.append(
                int(b)
            )

        rows.append(row)

    selected = largest_contiguous_run(
        qualifying
    )

    support = pd.DataFrame(rows)
    support[
        "selected_common_temporal_core"
    ] = support[
        "bin"
    ].isin(
        selected
    )

    return (
        selected,
        support,
    )


def temporal_sensitivity(
    arrays: Mapping[
        str,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    all_subjects: Sequence[int],
    temporal_subjects: Sequence[int],
    dts: Sequence[int],
    common_core_bins: Sequence[int],
    n_bootstrap: int,
    seed: int,
    min_n: int,
    min_subjects: int,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
]:
    if not common_core_bins:
        return (
            pd.DataFrame(),
            pd.DataFrame(),
        )

    subject_index = {
        int(s): i
        for i, s in enumerate(
            all_subjects
        )
    }
    idx = [
        subject_index[
            int(s)
        ]
        for s in temporal_subjects
    ]

    rng = np.random.default_rng(
        int(seed)
        + 15001
    )

    n_subjects = len(
        temporal_subjects
    )

    draws = rng.integers(
        0,
        n_subjects,
        size=(
            int(n_bootstrap),
            n_subjects,
        ),
    )

    dt_array = np.asarray(
        dts,
        dtype=float,
    )

    curve_rows = []
    slope_rows = []

    boot_slopes = {
        mode: np.full(
            int(n_bootstrap),
            np.nan,
            dtype=float,
        )
        for mode in FLUCTUATION_MODES
    }

    for mode in FLUCTUATION_MODES:
        mode_arrays = {
            k: v[
                idx,
                :,
                :,
            ]
            for k, v in arrays[
                mode
            ].items()
        }

        point = aggregate_fluctuation(
            mode_arrays,
            np.ones(
                n_subjects,
                dtype=float,
            ),
        )

        point_cross = point[
            "cross_cov"
        ][
            :,
            list(
                common_core_bins
            ),
        ]

        if not np.isfinite(point_cross).all():
            raise RuntimeError(
                f"{mode}: non-finite point covariance inside frozen common "
                "temporal sensitivity core."
            )

        point_curve = np.mean(
            point_cross,
            axis=1,
        )

        for j, dt in enumerate(
            dts
        ):
            curve_rows.append(
                {
                    "mode": mode,
                    "dt": int(dt),
                    "value": float(
                        point_curve[j]
                    ),
                    "n_common_core_bins": int(
                        len(
                            common_core_bins
                        )
                    ),
                    "n_subjects_universe": int(
                        n_subjects
                    ),
                }
            )

        fit = safe_linear_fit(
            dt_array - 1.0,
            point_curve,
        )

        boot_curves = np.full(
            (
                int(n_bootstrap),
                len(dts),
            ),
            np.nan,
            dtype=float,
        )

        for b in range(
            int(n_bootstrap)
        ):
            multiplicity = np.bincount(
                draws[b],
                minlength=n_subjects,
            ).astype(float)

            agg = aggregate_fluctuation(
                mode_arrays,
                multiplicity,
            )

            sampled_subject_support = (
                (mode_arrays["n"] > 0)
                & (multiplicity[:, None, None] > 0)
            ).sum(axis=0)

            selected = list(
                common_core_bins
            )

            valid_cells = (
                (agg["n"][:, selected] >= int(min_n))
                & (
                    sampled_subject_support[:, selected]
                    >= int(min_subjects)
                )
                & np.isfinite(
                    agg["cross_cov"][:, selected]
                )
            )

            # Fixed-core sensitivity: every selected bin must remain valid at
            # every requested lag in this bootstrap draw.
            if not bool(
                np.all(valid_cells)
            ):
                continue

            boot_curves[
                b,
                :,
            ] = np.mean(
                agg["cross_cov"][:, selected],
                axis=1,
            )

            if not np.isfinite(
                boot_curves[b, :]
            ).all():
                continue

            boot_fit = safe_linear_fit(
                dt_array - 1.0,
                boot_curves[
                    b,
                    :,
                ],
            )
            if boot_fit["n_points"] == len(dts):
                boot_slopes[
                    mode
                ][b] = boot_fit[
                    "slope"
                ]

        lo, med, hi, n_valid = (
            percentile_summary(
                boot_slopes[
                    mode
                ]
            )
        )

        slope_rows.append(
            {
                "mode": mode,
                "estimator": (
                    "matched_core_pooled_covariance_equal_bin"
                ),
                "observed_slope_per_week": fit[
                    "slope"
                ],
                "observed_intercept_at_dt1": fit[
                    "intercept"
                ],
                "n_dt": fit[
                    "n_points"
                ],
                "n_common_core_bins": int(
                    len(
                        common_core_bins
                    )
                ),
                "n_subjects_universe": int(
                    n_subjects
                ),
                "bootstrap_median_slope": med,
                "bootstrap_q025": lo,
                "bootstrap_q975": hi,
                "n_bootstrap_valid": n_valid,
                "bootstrap_positive_slope_fraction": (
                    float(
                        np.mean(
                            finite_numeric(
                                boot_slopes[
                                    mode
                                ]
                            )
                            > 0
                        )
                    )
                    if n_valid
                    else np.nan
                ),
            }
        )

    summary = pd.DataFrame(
        slope_rows
    )

    primary_boot = boot_slopes[
        "PRIMARY"
    ]

    for mode in [
        "NON_FF",
        "TT",
    ]:
        diff = (
            boot_slopes[
                mode
            ]
            - primary_boot
        )

        lo, med, hi, n_valid = (
            percentile_summary(
                diff
            )
        )

        mask = summary[
            "mode"
        ].eq(
            mode
        )

        summary.loc[
            mask,
            "bootstrap_slope_difference_vs_PRIMARY_median",
        ] = med
        summary.loc[
            mask,
            "bootstrap_slope_difference_vs_PRIMARY_q025",
        ] = lo
        summary.loc[
            mask,
            "bootstrap_slope_difference_vs_PRIMARY_q975",
        ] = hi
        summary.loc[
            mask,
            "n_bootstrap_difference_valid",
        ] = n_valid

        observed_primary = float(
            summary.loc[
                summary[
                    "mode"
                ].eq(
                    "PRIMARY"
                ),
                "observed_slope_per_week",
            ].iloc[0]
        )

        observed_mode = float(
            summary.loc[
                mask,
                "observed_slope_per_week",
            ].iloc[0]
        )

        summary.loc[
            mask,
            "observed_slope_difference_vs_PRIMARY",
        ] = (
            observed_mode
            - observed_primary
        )

    return (
        pd.DataFrame(
            curve_rows
        ),
        summary,
    )


# =============================================================================
# Reference audits and compact summaries
# =============================================================================

def audit_forward_primary(
    recomputed: pd.DataFrame,
    step9_primary: pd.DataFrame,
) -> pd.DataFrame:
    ref = step9_primary[
        step9_primary[
            "representation"
        ].astype(str).eq(
            "cross_combined"
        )
    ][
        [
            "bin",
            "mean_dx",
        ]
    ].rename(
        columns={
            "mean_dx":
                "step9_mean_dx",
        }
    )

    cur = recomputed[
        recomputed[
            "mode"
        ].astype(str).eq(
            "PRIMARY"
        )
    ][
        [
            "bin",
            "mean_dx",
        ]
    ].rename(
        columns={
            "mean_dx":
                "step15_recomputed_mean_dx",
        }
    )

    out = cur.merge(
        ref,
        on="bin",
        how="inner",
        validate="one_to_one",
    )

    out[
        "difference"
    ] = (
        out[
            "step15_recomputed_mean_dx"
        ]
        - out[
            "step9_mean_dx"
        ]
    )

    return out


def audit_fluctuation_primary(
    recomputed: pd.DataFrame,
    step11_primary: pd.DataFrame,
) -> pd.DataFrame:
    ref = step11_primary[
        [
            "dt",
            "bin",
            "cross_cov",
            "same_var_mean",
            "replicate_specific_excess",
        ]
    ].rename(
        columns={
            "cross_cov":
                "step11_cross_cov",
            "same_var_mean":
                "step11_same_var_mean",
            "replicate_specific_excess":
                "step11_replicate_specific_excess",
        }
    )

    cur = recomputed[
        recomputed[
            "mode"
        ].astype(str).eq(
            "PRIMARY"
        )
    ][
        [
            "dt",
            "bin",
            "cross_cov",
            "same_var_mean",
            "replicate_specific_excess",
        ]
    ].rename(
        columns={
            "cross_cov":
                "step15_cross_cov",
            "same_var_mean":
                "step15_same_var_mean",
            "replicate_specific_excess":
                "step15_replicate_specific_excess",
        }
    )

    out = cur.merge(
        ref,
        on=[
            "dt",
            "bin",
        ],
        how="inner",
        validate="one_to_one",
    )

    for metric in METRICS:
        out[
            f"{metric}_difference"
        ] = (
            out[
                f"step15_{metric}"
            ]
            - out[
                f"step11_{metric}"
            ]
        )

    return out


def build_forward_difference_table(
    curves: pd.DataFrame,
    boot_curves: Mapping[
        str,
        np.ndarray,
    ],
) -> pd.DataFrame:
    return build_forward_differences(
        curves,
        boot_curves,
    )


def build_fluctuation_difference_table(
    point: pd.DataFrame,
    boot: Mapping[
        str,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    dts: Sequence[int],
    fluct_edges: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    primary = point[
        point[
            "mode"
        ].eq(
            "PRIMARY"
        )
    ].set_index(
        [
            "dt",
            "bin",
        ]
    )

    for mode in [
        "NON_FF",
        "TT",
    ]:
        current = point[
            point[
                "mode"
            ].eq(
                mode
            )
        ].set_index(
            [
                "dt",
                "bin",
            ]
        )

        for j, dt in enumerate(
            dts
        ):
            for _, edge in fluct_edges.iterrows():
                b = int(
                    edge["bin"]
                )

                key = (
                    int(dt),
                    b,
                )

                for metric in METRICS:
                    cur = (
                        current.loc[
                            key,
                            metric,
                        ]
                        if key in current.index
                        else np.nan
                    )

                    ref = (
                        primary.loc[
                            key,
                            metric,
                        ]
                        if key in primary.index
                        else np.nan
                    )

                    diff_boot = (
                        boot[
                            mode
                        ][
                            metric
                        ][
                            :,
                            j,
                            b,
                        ]
                        - boot[
                            "PRIMARY"
                        ][
                            metric
                        ][
                            :,
                            j,
                            b,
                        ]
                    )

                    lo, med, hi, n_valid = (
                        percentile_summary(
                            diff_boot
                        )
                    )

                    rows.append(
                        {
                            "mode": mode,
                            "contrast":
                                f"{mode}-PRIMARY",
                            "metric": metric,
                            "dt": int(dt),
                            "bin": b,
                            "x_left": float(
                                edge["x_left"]
                            ),
                            "x_right": float(
                                edge["x_right"]
                            ),
                            "x_center": float(
                                edge["x_center"]
                            ),
                            "mode_value": cur,
                            "primary_value": ref,
                            "difference": (
                                cur - ref
                                if pd.notna(
                                    cur
                                )
                                and pd.notna(
                                    ref
                                )
                                else np.nan
                            ),
                            "bootstrap_difference_median": med,
                            "bootstrap_difference_q025": lo,
                            "bootstrap_difference_q975": hi,
                            "n_bootstrap_difference_valid": n_valid,
                        }
                    )

    return pd.DataFrame(rows)


def summarize_fluctuation_domains(
    point: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for (
        mode,
        metric,
        dt,
    ), g in point.melt(
        id_vars=[
            "mode",
            "dt",
            "bin",
            "x_center",
            "n",
            "n_subjects",
            "meets_min_support",
        ],
        value_vars=list(
            METRICS
        ),
        var_name="metric",
        value_name="value",
    ).groupby(
        [
            "mode",
            "metric",
            "dt",
        ],
        sort=True,
    ):
        values = finite_numeric(
            g["value"]
        )

        rows.append(
            {
                "mode": mode,
                "metric": metric,
                "dt": int(dt),
                "n_valid_bins": int(
                    len(
                        values
                    )
                ),
                "median_across_bins": (
                    float(
                        np.median(
                            values
                        )
                    )
                    if values.size
                    else np.nan
                ),
                "mean_across_bins": (
                    float(
                        np.mean(
                            values
                        )
                    )
                    if values.size
                    else np.nan
                ),
                "min_across_bins": (
                    float(
                        np.min(
                            values
                        )
                    )
                    if values.size
                    else np.nan
                ),
                "max_across_bins": (
                    float(
                        np.max(
                            values
                        )
                    )
                    if values.size
                    else np.nan
                ),
            }
        )

    return pd.DataFrame(rows)


def build_analysis_summary(
    composition: pd.DataFrame,
    forward_summary: pd.DataFrame,
    fluct_summary: pd.DataFrame,
    temporal_summary: pd.DataFrame,
    forward_audit: pd.DataFrame,
    fluct_audit: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []

    for r in forward_summary.itertuples(
        index=False
    ):
        rows.append(
            {
                "section": "forward",
                "mode": r.mode,
                "metric": "mean_dx_slope_vs_abundance",
                "dt": 1,
                "value": r.observed_slope,
                "lower": r.bootstrap_q025,
                "upper": r.bootstrap_q975,
                "n": r.n_valid_bins,
                "note": "observed replicate-decoupled AB/BA; fixed Step-9 bins",
            }
        )

    fs = fluct_summary[
        fluct_summary[
            "metric"
        ].eq(
            "cross_cov"
        )
    ]

    for r in fs.itertuples(
        index=False
    ):
        rows.append(
            {
                "section": "fluctuation",
                "mode": r.mode,
                "metric": "median_cross_cov_across_bins",
                "dt": int(r.dt),
                "value": r.median_across_bins,
                "lower": np.nan,
                "upper": np.nan,
                "n": r.n_valid_bins,
                "note": "common4; fixed Step-11 xmid bins",
            }
        )

    if not temporal_summary.empty:
        for r in temporal_summary.itertuples(
            index=False
        ):
            rows.append(
                {
                    "section": "temporal",
                    "mode": r.mode,
                    "metric": "matched_core_cross_cov_slope_per_week",
                    "dt": np.nan,
                    "value": r.observed_slope_per_week,
                    "lower": r.bootstrap_q025,
                    "upper": r.bootstrap_q975,
                    "n": r.n_common_core_bins,
                    "note": "fixed Step-14 subjects; matched operational-domain core",
                }
            )

    max_forward_audit = (
        float(
            np.nanmax(
                np.abs(
                    pd.to_numeric(
                        forward_audit[
                            "difference"
                        ],
                        errors="coerce",
                    )
                )
            )
        )
        if not forward_audit.empty
        else np.nan
    )

    fluct_diffs = []
    if not fluct_audit.empty:
        for metric in METRICS:
            col = (
                f"{metric}_difference"
            )
            if col in fluct_audit.columns:
                fluct_diffs.extend(
                    finite_numeric(
                        fluct_audit[
                            col
                        ]
                    ).tolist()
                )

    max_fluct_audit = (
        float(
            np.max(
                np.abs(
                    fluct_diffs
                )
            )
        )
        if fluct_diffs
        else np.nan
    )

    rows.extend(
        [
            {
                "section": "audit",
                "mode": "PRIMARY",
                "metric": "max_abs_forward_recompute_difference_vs_step9",
                "dt": 1,
                "value": max_forward_audit,
                "lower": np.nan,
                "upper": np.nan,
                "n": len(
                    forward_audit
                ),
                "note": "should be numerical-level only",
            },
            {
                "section": "audit",
                "mode": "PRIMARY",
                "metric": "max_abs_fluctuation_recompute_difference_vs_step11",
                "dt": np.nan,
                "value": max_fluct_audit,
                "lower": np.nan,
                "upper": np.nan,
                "n": len(
                    fluct_audit
                ),
                "note": "should be numerical-level only",
            },
        ]
    )

    return pd.DataFrame(rows)


def write_readme(outdir: Path) -> None:
    text = """# Step 15 — fixed-threshold operational observation-domain sensitivity

Reference operational threshold
-------------------------------
alpha = 0.05 by default.

T/F is defined from Step-5 endpoint p-values:
T when p_endpoint < alpha, F otherwise.

T/F is an operational observation label and must not be described as biological
presence/absence.

Primary analyses remain unchanged
---------------------------------
Forward:
observed replicate-decoupled AB/BA estimator from Step 9.

Fluctuations:
cross_cov = Cov(dx_observed_rep1, dx_observed_rep2 | xmid_latent, dt, common4).

Step 15 does not replace these primary estimands.

Forward domains
---------------
The Step-9 reported abundance coordinates are retained exactly. Internal
assignment boundaries are calibrated, where necessary, inside numerically
equivalent boundary intervals so that the finalized Step-9 AB and BA bin
counts are reproduced exactly after the pandas/PyArrow backend change.

PRIMARY  : no operational-class filter
T0PLUS   : T at t0, t1 unrestricted = TT + TF
TT       : T at both endpoints
TF       : descriptive T0PLUS component

AB and BA are combined only when both folds meet support; no single-fold
fallback is allowed.

Fluctuation domains
-------------------
PRIMARY  : all common4
NON_FF   : common4 and class != FF
TT       : common4 and class == TT

Temporal sensitivity
--------------------
The Step-14 biological-subject universe and Step-13 abundance core are frozen.
Within that core, the largest contiguous bin subset with valid pooled support
for PRIMARY, NON_FF and TT at every lag is used as one common sensitivity core.

The temporal sensitivity estimator is:
matched_core_pooled_covariance_equal_bin.

Within every subject-bootstrap draw the full fixed common operational core must
remain supported at every lag; draws with incomplete support do not contribute
a temporal slope.

Forward bootstrap sufficient statistics are aggregated with an explicit
float64 subject-weighted sum rather than BLAS-backed matrix multiplication.
This is numerically equivalent for the small subject x bin matrices and avoids
platform-specific NumPy/Accelerate warnings on Apple silicon.

This is a sensitivity estimator and is not substituted for the Step-13
equal-bin/equal-subject primary estimator.

Interpretation boundary
-----------------------
Step 15 asks whether the main conclusions materially depend on operational
observation-domain selection at fixed alpha.

Step 16 varies alpha and is a separate threshold-robustness analysis.
"""
    (
        outdir
        / "README_outputs.md"
    ).write_text(
        text,
        encoding="utf-8",
    )


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Step 15: fixed-alpha operational observation-domain sensitivity "
            "for current replicate-resolved forward and fluctuation estimands."
        ),
    )

    p.add_argument(
        "--transitions",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--step9-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--step11-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--step13-dir",
        required=True,
        type=Path,
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
        "--alpha",
        type=float,
        default=0.05,
    )
    p.add_argument(
        "--forward-dt",
        type=int,
        default=1,
    )
    p.add_argument(
        "--min-n",
        type=int,
        default=50,
    )
    p.add_argument(
        "--min-subjects",
        type=int,
        default=2,
    )
    p.add_argument(
        "--n-bootstrap",
        type=int,
        default=2000,
    )
    p.add_argument(
        "--seed",
        type=int,
        default=123,
    )

    return p


def validate_args(args: argparse.Namespace) -> None:
    if not (
        0.0 < float(
            args.alpha
        ) < 1.0
    ):
        raise ValueError(
            "--alpha must be in (0,1)"
        )
    if int(
        args.forward_dt
    ) < 1:
        raise ValueError(
            "--forward-dt must be >=1"
        )
    if int(
        args.min_n
    ) < 2:
        raise ValueError(
            "--min-n must be >=2"
        )
    if int(
        args.min_subjects
    ) < 2:
        raise ValueError(
            "--min-subjects must be >=2"
        )
    if int(
        args.n_bootstrap
    ) < 1:
        raise ValueError(
            "--n-bootstrap must be >=1"
        )


# =============================================================================
# Main
# =============================================================================

def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(
        argv
    )
    validate_args(args)

    transitions = (
        args.transitions
        .expanduser()
        .resolve(strict=True)
    )
    step9_dir = (
        args.step9_dir
        .expanduser()
        .resolve(strict=True)
    )
    step11_dir = (
        args.step11_dir
        .expanduser()
        .resolve(strict=True)
    )
    step13_dir = (
        args.step13_dir
        .expanduser()
        .resolve(strict=True)
    )
    step14_dir = (
        args.step14_dir
        .expanduser()
        .resolve(strict=True)
    )
    outdir = ensure_dir(
        args.outdir
        .expanduser()
        .resolve()
    )

    step9_config_path = (
        step9_dir / "00_run_config.json"
    )
    step9_signature_path = (
        step9_dir / "00_run_signature.json"
    )
    step11_config_path = (
        step11_dir / "00_run_config.json"
    )
    step11_signature_path = (
        step11_dir / "00_run_signature.json"
    )
    step13_config_path = (
        step13_dir / "00_run_config.json"
    )
    step13_signature_path = (
        step13_dir / "00_run_signature.json"
    )
    step14_config_path = (
        step14_dir / "00_run_config.json"
    )
    step14_signature_path = (
        step14_dir / "00_run_signature.json"
    )

    step9_config = read_json(
        step9_config_path
    )
    step9_signature = read_json(
        step9_signature_path
    )
    step11_config = read_json(
        step11_config_path
    )
    step11_signature = read_json(
        step11_signature_path
    )
    step13_config = read_json(
        step13_config_path
    )
    step13_signature = read_json(
        step13_signature_path
    )
    step14_config = read_json(
        step14_config_path
    )
    step14_signature = read_json(
        step14_signature_path
    )

    validate_upstream(
        step9_config,
        step9_signature,
        step11_config,
        step11_signature,
        step13_config,
        step13_signature,
        step14_config,
        step14_signature,
    )

    step9_min_n = int(
        step9_config.get("binning", {}).get(
            "min_n",
            args.min_n,
        )
    )
    step11_min_n = int(
        step11_config.get("binning", {}).get(
            "min_n_per_pooled_cell",
            args.min_n,
        )
    )

    if (
        int(args.min_n) != step9_min_n
        or int(args.min_n) != step11_min_n
    ):
        raise ValueError(
            "For the primary fixed-alpha Step-15 run, --min-n must match "
            f"both finalized Step 9 and Step 11 thresholds. "
            f"requested={args.min_n}, Step9={step9_min_n}, "
            f"Step11={step11_min_n}."
        )

    forward_edges = normalize_edges(
        read_csv_required(
            step9_dir
            / "01_forward_bin_edges.csv",
            "Step-9 forward bin edges",
        ),
        "Step-9 forward bin edges",
    )

    step9_primary = read_csv_required(
        step9_dir
        / "02_primary_forward_by_bin.csv",
        "Step-9 primary forward curves",
    )

    step11_primary = read_csv_required(
        step11_dir
        / "02_fluctuation_by_bin_dt.csv",
        "Step-11 primary fluctuation curves",
    )

    fluct_edges = normalize_edges(
        step11_primary[
            [
                "bin",
                "x_left",
                "x_right",
                "x_center",
            ]
        ],
        "Step-11 fluctuation bin grid",
    )

    (
        frozen_core,
        temporal_subjects,
        temporal_dts,
        frozen_core_bins,
    ) = load_step14_contract(
        step14_dir
    )

    step13_slope = read_csv_required(
        step13_dir
        / "05_temporal_slope_summary.csv",
        "Step-13 temporal slope summary",
    )

    print(
        "[INFO] Loading selected Step-5 transition columns with pandas/PyArrow..."
    )
    base_df = read_transition_table(
        transitions
    )

    memory_gib = (
        float(
            base_df.memory_usage(
                index=True,
                deep=True,
            ).sum()
        )
        / (1024.0 ** 3)
    )
    print(
        f"[INFO] Transition working table: {len(base_df):,} rows; "
        f"{memory_gib:.2f} GiB in pandas"
    )

    all_subjects = sorted(
        int(x)
        for x in base_df[
            "subject"
        ].dropna().unique()
    )

    all_dts = sorted(
        int(x)
        for x in base_df[
            "dt"
        ].dropna().unique()
    )

    (
        forward_assignment_edges,
        forward_bin_calibration,
    ) = calibrate_forward_assignment_edges(
        base_df,
        forward_edges,
        step9_primary,
        int(
            args.forward_dt
        ),
    )

    n_changed_boundaries = int(
        forward_bin_calibration[
            "changed"
        ].sum()
    )

    max_boundary_shift = float(
        np.max(
            np.abs(
                forward_bin_calibration[
                    "difference"
                ].to_numpy(float)
            )
        )
    )

    print(
        "[INFO] Step-9 forward bin-membership calibration: "
        f"{n_changed_boundaries} internal boundaries adjusted; "
        f"max |shift|={max_boundary_shift:.3e}"
    )

    missing_dts = sorted(
        set(
            temporal_dts
        )
        - set(
            all_dts
        )
    )
    if missing_dts:
        raise ValueError(
            f"Step-14 temporal lags missing from transitions: {missing_dts}"
        )

    print(
        "[INFO] Execution backend: pandas/PyArrow (Polars disabled for Step 15)"
    )

    print(
        "[1/6] Operational class composition and boundary localization..."
    )

    composition, crossing = (
        build_operational_composition(
            base_df,
            float(
                args.alpha
            ),
            fluct_edges,
        )
    )

    print(
        "[2/6] Forward observation-domain sensitivity..."
    )

    forward_suff = (
        build_forward_sufficient(
            base_df,
            float(
                args.alpha
            ),
            int(
                args.forward_dt
            ),
            forward_assignment_edges,
        )
    )

    forward_counts, forward_sums = (
        build_forward_arrays(
            forward_suff,
            all_subjects,
            len(
                forward_edges
            ),
        )
    )

    forward_curves = (
        forward_point_curves(
            forward_counts,
            forward_sums,
            forward_edges,
            int(
                args.min_n
            ),
        )
    )

    (
        forward_boot_ci,
        forward_summary,
        forward_draws,
        forward_boot_curves,
        forward_boot_slope_arrays,
    ) = bootstrap_forward(
        forward_counts,
        forward_sums,
        all_subjects,
        forward_edges,
        int(
            args.min_n
        ),
        int(
            args.n_bootstrap
        ),
        int(
            args.seed
        ),
    )

    forward_curves = forward_curves.merge(
        forward_boot_ci,
        on=[
            "mode",
            "bin",
        ],
        how="left",
        validate="one_to_one",
    )

    forward_diff = (
        build_forward_difference_table(
            forward_curves,
            forward_boot_curves,
        )
    )

    forward_audit = (
        audit_forward_primary(
            forward_curves,
            step9_primary,
        )
    )

    print(
        "[3/6] Shared-fluctuation observation-domain sensitivity..."
    )

    fluct_suff = (
        build_fluctuation_sufficient(
            base_df,
            float(
                args.alpha
            ),
            fluct_edges,
        )
    )

    fluct_arrays = (
        build_fluctuation_arrays(
            fluct_suff,
            all_subjects,
            temporal_dts,
            len(
                fluct_edges
            ),
        )
    )

    fluct_point = (
        fluctuation_point_table(
            fluct_arrays,
            all_subjects,
            temporal_dts,
            fluct_edges,
            int(
                args.min_n
            ),
            int(
                args.min_subjects
            ),
        )
    )

    (
        fluct_boot_ci,
        fluct_boot,
        fluct_draws,
    ) = bootstrap_fluctuation(
        fluct_arrays,
        all_subjects,
        temporal_dts,
        fluct_edges,
        int(
            args.min_n
        ),
        int(
            args.min_subjects
        ),
        int(
            args.n_bootstrap
        ),
        int(
            args.seed
        )
        + 1000,
    )

    fluct_boot_ci_wide = (
        wide_fluctuation_bootstrap_ci(
            fluct_boot_ci
        )
    )

    fluct_point = fluct_point.merge(
        fluct_boot_ci_wide,
        on=[
            "mode",
            "dt",
            "bin",
        ],
        how="left",
        validate="one_to_one",
    )

    fluct_diff = (
        build_fluctuation_difference_table(
            fluct_point,
            fluct_boot,
            temporal_dts,
            fluct_edges,
        )
    )

    fluct_summary = (
        summarize_fluctuation_domains(
            fluct_point
        )
    )

    fluct_audit = (
        audit_fluctuation_primary(
            fluct_point,
            step11_primary[
                step11_primary[
                    "dt"
                ].astype(int).isin(
                    temporal_dts
                )
            ].copy(),
        )
    )

    print(
        "[4/6] Matched operational-domain temporal sensitivity..."
    )

    common_temporal_core, support_audit = (
        temporal_support_common_core(
            fluct_arrays,
            all_subjects,
            temporal_subjects,
            temporal_dts,
            frozen_core_bins,
            int(
                args.min_n
            ),
            int(
                args.min_subjects
            ),
        )
    )

    temporal_curve, temporal_summary = (
        temporal_sensitivity(
            fluct_arrays,
            all_subjects,
            temporal_subjects,
            temporal_dts,
            common_temporal_core,
            int(
                args.n_bootstrap
            ),
            int(
                args.seed
            ),
            int(
                args.min_n
            ),
            int(
                args.min_subjects
            ),
        )
    )

    if not temporal_summary.empty:
        ref = step13_slope[
            step13_slope[
                "estimator"
            ].astype(str).eq(
                "equal_bin_equal_subject"
            )
            & step13_slope[
                "metric"
            ].astype(str).eq(
                "cross_cov"
            )
        ]

        if len(ref) == 1:
            r = ref.iloc[0]
            temporal_summary[
                "step13_primary_slope_per_week"
            ] = float(
                r[
                    "observed_slope_per_week"
                ]
            )
            temporal_summary[
                "step13_primary_bootstrap_q025"
            ] = float(
                r[
                    "bootstrap_q025"
                ]
            )
            temporal_summary[
                "step13_primary_bootstrap_q975"
            ] = float(
                r[
                    "bootstrap_q975"
                ]
            )

    print(
        "[5/6] Compact summaries and audits..."
    )

    analysis_summary = (
        build_analysis_summary(
            composition,
            forward_summary,
            fluct_summary,
            temporal_summary,
            forward_audit,
            fluct_audit,
        )
    )

    print(
        "[6/6] Writing outputs..."
    )

    composition.to_csv(
        outdir
        / "01_operational_class_composition_by_dt.csv",
        index=False,
    )

    crossing.to_csv(
        outdir
        / "02_boundary_crossing_by_dt_bin.csv",
        index=False,
    )

    forward_curves.to_csv(
        outdir
        / "03_forward_domain_by_bin.csv",
        index=False,
    )

    forward_summary.to_csv(
        outdir
        / "04_forward_domain_summary.csv",
        index=False,
    )

    forward_diff.to_csv(
        outdir
        / "05_forward_difference_vs_primary.csv",
        index=False,
    )

    forward_audit.to_csv(
        outdir
        / "06_forward_primary_audit.csv",
        index=False,
    )

    forward_bin_calibration.to_csv(
        outdir
        / "06a_forward_bin_assignment_calibration.csv",
        index=False,
    )

    fluct_point.to_csv(
        outdir
        / "07_fluctuation_domain_by_bin_dt.csv",
        index=False,
    )

    fluct_diff.to_csv(
        outdir
        / "08_fluctuation_difference_vs_primary.csv",
        index=False,
    )

    fluct_summary.to_csv(
        outdir
        / "09_fluctuation_domain_summary.csv",
        index=False,
    )

    fluct_audit.to_csv(
        outdir
        / "10_fluctuation_primary_audit.csv",
        index=False,
    )

    common_core_table = frozen_core.copy()
    common_core_table[
        "selected_common_operational_temporal_core"
    ] = common_core_table[
        "bin"
    ].isin(
        common_temporal_core
    )

    common_core_table.to_csv(
        outdir
        / "11_temporal_common_core.csv",
        index=False,
    )

    temporal_curve.to_csv(
        outdir
        / "12_temporal_domain_by_dt.csv",
        index=False,
    )

    temporal_summary.to_csv(
        outdir
        / "13_temporal_domain_slope_summary.csv",
        index=False,
    )

    support_audit.to_csv(
        outdir
        / "14_temporal_support_audit.csv",
        index=False,
    )

    analysis_summary.to_csv(
        outdir
        / "00_analysis_summary.csv",
        index=False,
    )

    signature_payload = {
        "script_version": SCRIPT_VERSION,
        "transitions": path_identity(transitions),
        "step9_analysis_signature": step9_config.get(
            "analysis_signature"
        ),
        "step11_analysis_signature": step11_config.get(
            "analysis_signature"
        ),
        "step13_analysis_signature": step13_config.get(
            "analysis_signature"
        ),
        "step14_analysis_signature": step14_config.get(
            "analysis_signature"
        ),
        "parameters": {
            "alpha": float(args.alpha),
            "forward_dt": int(args.forward_dt),
            "min_n": int(args.min_n),
            "min_subjects": int(args.min_subjects),
            "n_bootstrap": int(args.n_bootstrap),
            "seed": int(args.seed),
        },
        "temporal_subjects": [
            int(x) for x in temporal_subjects
        ],
        "frozen_step13_core_bins": [
            int(x) for x in frozen_core_bins
        ],
        "common_operational_core_bins": [
            int(x) for x in common_temporal_core
        ],
        "forward_fold_combination":
            "exact equal AB/BA weight; both folds required",
        "bootstrap_cell_support_rechecked_each_draw": True,
        "temporal_bootstrap_requires_complete_fixed_core_all_dt": True,
    }
    analysis_signature = stable_sha256(
        signature_payload
    )

    (
        outdir / "00_run_signature.json"
    ).write_text(
        json.dumps(
            {
                "analysis_signature": analysis_signature,
                "signature_payload": signature_payload,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    input_manifest = pd.DataFrame(
        [
            {"input": "transitions", **path_identity(transitions)},
            {"input": "step9_run_config", **path_identity(step9_config_path)},
            {"input": "step9_run_signature", **path_identity(step9_signature_path)},
            {"input": "step11_run_config", **path_identity(step11_config_path)},
            {"input": "step11_run_signature", **path_identity(step11_signature_path)},
            {"input": "step13_run_config", **path_identity(step13_config_path)},
            {"input": "step13_run_signature", **path_identity(step13_signature_path)},
            {"input": "step14_run_config", **path_identity(step14_config_path)},
            {"input": "step14_run_signature", **path_identity(step14_signature_path)},
            {
                "input": "step14_step15_contract",
                **path_identity(step14_dir / "13_step15_contract.csv"),
            },
        ]
    )
    input_manifest.to_csv(
        outdir / "00_input_manifest.csv",
        index=False,
    )

    run_config = {
        "script": Path(
            __file__
        ).name,
        "script_version": SCRIPT_VERSION,
        "analysis_signature": analysis_signature,
        "transitions": str(
            transitions
        ),
        "step9_dir": str(
            step9_dir
        ),
        "step11_dir": str(
            step11_dir
        ),
        "step13_dir": str(
            step13_dir
        ),
        "step14_dir": str(
            step14_dir
        ),
        "outdir": str(
            outdir
        ),
        "reference_operational_alpha": float(
            args.alpha
        ),
        "operational_definition": {
            "T": "p_endpoint < alpha",
            "F": "p_endpoint >= alpha",
            "interpretation": (
                "operational observation state; not biological presence/absence"
            ),
        },
        "primary_estimands_unchanged": True,
        "step12_used_computationally": False,
        "forward": {
            "dt": int(
                args.forward_dt
            ),
            "conditioning": {
                "AB":
                    "x_observed_rep1_t0",
                "BA":
                    "x_observed_rep2_t0",
            },
            "displacement": {
                "AB":
                    "dx_observed_rep2",
                "BA":
                    "dx_observed_rep1",
            },
            "base_support": {
                "AB":
                    "forward_ab_eligible",
                "BA":
                    "forward_ba_eligible",
            },
            "operational_domains": list(
                FORWARD_MODES
            ),
            "bin_grid": (
                "reported coordinates from exact Step-9 01_forward_bin_edges.csv"
            ),
            "bin_membership_calibration": {
                "purpose": (
                    "reproduce finalized Step-9 AB/BA bin counts exactly after "
                    "CSV/backend boundary round-trip"
                ),
                "n_internal_boundaries_adjusted": int(
                    n_changed_boundaries
                ),
                "max_abs_assignment_boundary_shift": float(
                    max_boundary_shift
                ),
                "reported_coordinates_changed": False,
            },
            "fold_combination": (
                "exact 0.5 AB + 0.5 BA after binning; both folds required"
            ),
            "single_fold_fallback": False,
        },
        "fluctuations": {
            "base_support": "common4",
            "conditioning": "xmid_latent",
            "dx_A":
                "dx_observed_rep1",
            "dx_B":
                "dx_observed_rep2",
            "operational_domains": list(
                FLUCTUATION_MODES
            ),
            "bin_grid": (
                "exact Step-11 fixed xmid grid"
            ),
        },
        "temporal_sensitivity": {
            "subject_universe": [
                int(x)
                for x in temporal_subjects
            ],
            "frozen_step13_core_bins": [
                int(x)
                for x in frozen_core_bins
            ],
            "common_operational_core_bins": [
                int(x)
                for x in common_temporal_core
            ],
            "dt_values": [
                int(x)
                for x in temporal_dts
            ],
            "estimator": (
                "matched_core_pooled_covariance_equal_bin"
            ),
            "interpretive_role": (
                "observation-domain sensitivity only; does not replace Step-13 primary estimator"
            ),
        },
        "support_thresholds": {
            "min_n": int(
                args.min_n
            ),
            "min_subjects": int(
                args.min_subjects
            ),
        },
        "bootstrap": {
            "unit":
                "biological subject",
            "n_bootstrap": int(
                args.n_bootstrap
            ),
            "same_subject_draw_across_operational_domains":
                True,
            "forward_bootstrap_subject_aggregation":
                "explicit float64 broadcasted weighted sum; BLAS matmul disabled",
            "cell_support_rechecked_within_each_draw": True,
            "temporal_slope_requires_complete_fixed_core_all_dt": True,
            "seed": int(
                args.seed
            ),
        },
    }

    (
        outdir
        / "00_run_config.json"
    ).write_text(
        json.dumps(
            json_safe(
                run_config
            ),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    write_readme(
        outdir
    )
    write_manifest(
        outdir
    )

    print(
        "\n[DONE] Step 15 operational observation-domain sensitivity"
    )
    print(
        "Reference alpha          :",
        float(
            args.alpha
        ),
    )
    print(
        "Biological subjects      :",
        all_subjects,
    )
    print(
        "Temporal subject universe:",
        temporal_subjects,
    )
    print(
        "Frozen Step-13 core bins :",
        frozen_core_bins,
    )
    print(
        "Common operational core  :",
        common_temporal_core,
    )

    if not forward_audit.empty:
        print(
            "Forward primary audit max |diff|:",
            float(
                np.nanmax(
                    np.abs(
                        forward_audit[
                            "difference"
                        ].to_numpy(float)
                    )
                )
            ),
        )

    if not fluct_audit.empty:
        audit_values = []
        for metric in METRICS:
            audit_values.extend(
                finite_numeric(
                    fluct_audit[
                        f"{metric}_difference"
                    ]
                ).tolist()
            )
        print(
            "Fluctuation primary audit max |diff|:",
            (
                float(
                    np.max(
                        np.abs(
                            audit_values
                        )
                    )
                )
                if audit_values
                else np.nan
            ),
        )

    print(
        "\nOutputs:"
    )
    for filename in [
        "00_analysis_summary.csv",
        "00_run_signature.json",
        "00_input_manifest.csv",
        "01_operational_class_composition_by_dt.csv",
        "02_boundary_crossing_by_dt_bin.csv",
        "03_forward_domain_by_bin.csv",
        "04_forward_domain_summary.csv",
        "05_forward_difference_vs_primary.csv",
        "06_forward_primary_audit.csv",
        "06a_forward_bin_assignment_calibration.csv",
        "07_fluctuation_domain_by_bin_dt.csv",
        "08_fluctuation_difference_vs_primary.csv",
        "09_fluctuation_domain_summary.csv",
        "10_fluctuation_primary_audit.csv",
        "11_temporal_common_core.csv",
        "12_temporal_domain_by_dt.csv",
        "13_temporal_domain_slope_summary.csv",
        "14_temporal_support_audit.csv",
    ]:
        print(
            " -",
            outdir / filename,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
