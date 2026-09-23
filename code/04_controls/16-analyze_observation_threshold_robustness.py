#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
16-analyze_observation_threshold_robustness.py
===============================================

Canonical ClonoDynamics Step 16 after the replicate-resolved dynamics pivot.

SCIENTIFIC QUESTION
-------------------
Step 15 showed that fixing the operational observation threshold at alpha=0.05
and changing the included operational transition domain can substantially
reshape abundance-resolved forward drift and fluctuation amplitude, while the
absence of positive temporal accumulation remains robust.

Step 16 asks whether those Step-15 sensitivity conclusions themselves depend
on the NUMERICAL operational threshold.

The canonical threshold sweep is:

    alpha = 0.10, 0.05, 0.025, 0.01

with:

    reference alpha = 0.05

Only the operational T/F classification changes.

Latent-frequency inference is NOT refit.
The observed replicate-resolved dynamical estimands are NOT redefined.

OPERATIONAL CLASSIFICATION
--------------------------
For each endpoint:

    T(alpha) = 1[p_endpoint < alpha]
    F(alpha) = 1[p_endpoint >= alpha]

and transitions are assigned to:

    TT, TF, FT, FF.

T/F is an OPERATIONAL OBSERVATION STATE.
It must not be interpreted as biological presence/absence.

CURRENT DYNAMICAL ARCHITECTURE
------------------------------
Step 16 reuses the finalized Step-15 v6 implementation directly, including
exact equal-weight AB/BA support and per-bootstrap support rechecking.

Forward drift
~~~~~~~~~~~~~
Observed replicate-decoupled estimator inherited from Step 9:

    AB:
        condition on x_observed_rep1_t0
        displacement = dx_observed_rep2

    BA:
        condition on x_observed_rep2_t0
        displacement = dx_observed_rep1

Operational domains:

    PRIMARY
        no operational-domain restriction; alpha invariant

    T0PLUS
        TT + TF

    TT

    TF
        descriptive boundary-crossing component of T0PLUS

The exact Step-9 reported abundance coordinates are preserved. Internal
assignment boundaries are calibrated by the Step-15 implementation so that the
finalized Step-9 AB/BA bin memberships are reproduced exactly.

Shared fluctuations
~~~~~~~~~~~~~~~~~~~
Current Step-11 estimand:

    cross_cov =
        Cov(dx_observed_rep1, dx_observed_rep2
            | xmid_latent, dt, common4)

with complementary:

    same_var_mean

    replicate_specific_excess =
        same_var_mean - cross_cov

Operational domains:

    PRIMARY
        all common4 transitions; alpha invariant

    NON_FF
        TT + TF + FT

    TT

Temporal scaling
~~~~~~~~~~~~~~~~
Step 16 freezes:

    - the Step-14 complete-case biological-subject universe;
    - the Step-13/14 abundance core;
    - temporal lags 1–5 weeks.

For each alpha, support is evaluated for PRIMARY, NON_FF and TT at every lag.
The Step-16 temporal comparison then uses ONE cross-alpha common operational
core:

    largest contiguous subset of the frozen Step-13/14 core that has valid
    support for all operational domains, all alpha values and all temporal lags.

Temporal covariance profiles are averaged equally across this common abundance
core. Signed temporal slopes are then fit using:

    M(dt) = intercept_at_dt1 + slope_per_week * (dt - 1)

with one joint biological-subject bootstrap.

CROSS-ALPHA PAIRING
-------------------
Bootstrap seeds and subject universes are fixed across alpha.

Therefore Step 16 computes paired cross-alpha differences for:

    - forward abundance-bin curves;
    - forward abundance slopes;
    - fluctuation curves;
    - temporal slopes.

All non-reference thresholds are compared against alpha=0.05 using identical
biological-subject bootstrap draws.

For common-support forward slopes, a bootstrap slope is retained only when all
cross-alpha common forward bins are finite in that draw. For temporal slopes, a
bootstrap slope is retained only when the entire cross-alpha common operational
core remains valid at every lag in that draw.

PRIMARY is expected to be alpha invariant by construction.

INTERPRETATION BOUNDARY
-----------------------
Step 16 is a NUMERICAL THRESHOLD ROBUSTNESS analysis.

It does NOT:

    refit latent abundance;
    threshold posterior-predictive detectability;
    use T/F as biological presence/absence;
    change the primary Step-9/11 estimands;
    subtract the Step-12 pseudo technical null;
    replace the Step-13 primary temporal estimator.

Step 16 is a sensitivity branch built on the validated current Step-15
implementation.

INPUTS
------
--transitions
    Final Step-5 longitudinal_transitions.parquet.

--step9-dir
    Final Step-9 output directory.

--step11-dir
    Final Step-11 output directory.

--step13-dir
    Final Step-13 output directory.

--step14-dir
    Final Step-14 output directory.

--step15-script
    Current validated Step-15 analyzer. The default is the sibling
    15-detectability_boundary_sensitivity.py.

Optional:
--step15-reference-dir
    Final alpha=0.05 Step-15 output directory, used only as an external
    reference audit. It is not used to calculate the cross-alpha results.

OUTPUTS
-------
00_run_config.json
00_run_signature.json
00_analysis_summary.csv
00_pipeline_manifest.csv
00_input_manifest.csv
README_outputs.md

01_class_composition_by_alpha_dt.csv
02_class_composition_pooled_by_alpha.csv
03_tt_retention_by_alpha.csv

04_forward_domain_cross_alpha_by_bin.csv
05_forward_slope_cross_alpha.csv
05a_forward_slope_native_support_diagnostic.csv
05b_forward_cross_alpha_common_support.csv
06_forward_curve_difference_vs_reference_alpha.csv
07_forward_slope_difference_vs_reference_alpha.csv

08_fluctuation_domain_cross_alpha_by_bin_dt.csv
09_fluctuation_difference_vs_reference_alpha.csv

10_temporal_support_by_alpha_bin.csv
11_cross_alpha_temporal_common_core.csv
12_temporal_profiles_cross_alpha.csv
13_temporal_slopes_cross_alpha.csv
14_temporal_slope_difference_vs_reference_alpha.csv

15_reference_alpha_internal_audit.csv
16_reference_step15_output_audit.csv          [when --step15-reference-dir]
17_support_summary_by_alpha.csv

PRIMARY RUN
-----------
python3 ./code/dynamics_4/16-analyze_observation_threshold_robustness.py \
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
    --step15-script \
    ./code/dynamics_4/15-detectability_boundary_sensitivity.py \
    --step15-reference-dir \
    ./results_4/15-detectability_boundary_sensitivity \
    --outdir \
    ./results_4/16-observation_threshold_robustness \
    --alphas 0.10 0.05 0.025 0.01 \
    --reference-alpha 0.05 \
    --min-n 50 \
    --min-subjects 2 \
    --n-bootstrap 2000 \
    --seed 123
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from types import ModuleType
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SCRIPT_VERSION = "v3-step15v6-fixed-bootstrap-support-2026-09-18"

DEFAULT_ALPHAS = (
    0.10,
    0.05,
    0.025,
    0.01,
)

REFERENCE_ALPHA = 0.05

FORWARD_SENSITIVITY_MODES = (
    "T0PLUS",
    "TT",
    "TF",
)

FLUCTUATION_SENSITIVITY_MODES = (
    "NON_FF",
    "TT",
)

TEMPORAL_MODES = (
    "PRIMARY",
    "NON_FF",
    "TT",
)

METRICS = (
    "cross_cov",
    "same_var_mean",
    "replicate_specific_excess",
)


# =============================================================================
# Generic helpers
# =============================================================================

def ensure_dir(
    path: Path,
) -> Path:
    path.mkdir(
        parents=True,
        exist_ok=True,
    )
    return path


def alpha_tag(
    alpha: float,
) -> str:
    text = (
        f"{float(alpha):.8g}"
        .rstrip("0")
        .rstrip(".")
    )
    return (
        text
        .replace("-", "m")
        .replace(".", "p")
    )


def json_safe(
    value,
):
    if isinstance(
        value,
        Path,
    ):
        return str(
            value
        )

    if isinstance(
        value,
        np.integer,
    ):
        return int(
            value
        )

    if isinstance(
        value,
        np.floating,
    ):
        return (
            None
            if np.isnan(
                value
            )
            else float(
                value
            )
        )

    if isinstance(
        value,
        np.bool_,
    ):
        return bool(
            value
        )

    if isinstance(
        value,
        dict,
    ):
        return {
            str(k):
                json_safe(v)
            for k, v
            in value.items()
        }

    if isinstance(
        value,
        (
            list,
            tuple,
        ),
    ):
        return [
            json_safe(v)
            for v in value
        ]

    return value


def finite_numeric(
    values,
) -> np.ndarray:
    x = pd.to_numeric(
        pd.Series(
            values
        ),
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    return x[
        np.isfinite(
            x
        )
    ]


def percentile_summary(
    values: np.ndarray,
) -> Tuple[
    float,
    float,
    float,
    int,
]:
    x = np.asarray(
        values,
        dtype=float,
    )

    x = x[
        np.isfinite(
            x
        )
    ]

    if x.size == 0:
        return (
            np.nan,
            np.nan,
            np.nan,
            0,
        )

    return (
        float(
            np.quantile(
                x,
                0.025,
            )
        ),
        float(
            np.median(
                x
            )
        ),
        float(
            np.quantile(
                x,
                0.975,
            )
        ),
        int(
            x.size
        ),
    )


def largest_contiguous_run(
    values: Iterable[int],
) -> List[int]:
    vals = sorted(
        set(
            int(x)
            for x in values
        )
    )

    if not vals:
        return []

    runs: List[
        List[int]
    ] = []

    current = [
        vals[0]
    ]

    for value in vals[
        1:
    ]:
        if (
            value
            == current[-1]
            + 1
        ):
            current.append(
                value
            )
        else:
            runs.append(
                current
            )
            current = [
                value
            ]

    runs.append(
        current
    )

    runs.sort(
        key=lambda run: (
            len(
                run
            ),
            -run[0],
        ),
        reverse=True,
    )

    return runs[0]


def write_manifest(
    outdir: Path,
) -> None:
    rows = []

    for path in sorted(
        outdir.rglob("*")
    ):
        if (
            path.is_file()
            and path.name
            != "00_pipeline_manifest.csv"
        ):
            rows.append(
                {
                    "relative_path":
                        str(
                            path.relative_to(
                                outdir
                            )
                        ),
                    "suffix":
                        path.suffix.lower(),
                    "size_bytes":
                        int(
                            path.stat().st_size
                        ),
                }
            )

    pd.DataFrame(
        rows
    ).to_csv(
        outdir
        / "00_pipeline_manifest.csv",
        index=False,
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


def load_python_module(
    path: Path,
) -> ModuleType:
    spec = (
        importlib.util
        .spec_from_file_location(
            "clonodynamics_step15_current",
            path,
        )
    )

    if (
        spec is None
        or spec.loader is None
    ):
        raise ImportError(
            f"Could not import Step-15 module from {path}"
        )

    module = (
        importlib.util
        .module_from_spec(
            spec
        )
    )

    spec.loader.exec_module(
        module
    )

    required = [
        "read_transition_table",
        "read_json",
        "validate_upstream",
        "normalize_edges",
        "load_step14_contract",
        "calibrate_forward_assignment_edges",
        "build_operational_composition",
        "build_forward_sufficient",
        "build_forward_arrays",
        "forward_point_curves",
        "bootstrap_forward",
        "build_fluctuation_sufficient",
        "build_fluctuation_arrays",
        "fluctuation_point_table",
        "bootstrap_fluctuation",
        "temporal_support_common_core",
        "aggregate_fluctuation",
        "safe_linear_fit",
    ]

    missing = [
        name
        for name in required
        if not hasattr(
            module,
            name,
        )
    ]

    if missing:
        raise ValueError(
            "The supplied Step-15 script is not the current validated "
            f"implementation. Missing functions: {missing}"
        )

    version = str(
        getattr(
            module,
            "SCRIPT_VERSION",
            "",
        )
    )
    if not version.startswith("v6-"):
        raise ValueError(
            "Step 16 v3 requires the finalized Step-15 v6 implementation; "
            f"received SCRIPT_VERSION={version!r}."
        )

    return module


# =============================================================================
# Cross-alpha composition and retention
# =============================================================================

def pooled_class_composition(
    class_by_dt: pd.DataFrame,
) -> pd.DataFrame:
    pooled = (
        class_by_dt.groupby(
            [
                "alpha",
                "obs_class_step15",
            ],
            as_index=False,
            sort=True,
        )[
            "n_transitions"
        ]
        .sum()
    )

    totals = (
        pooled.groupby(
            "alpha",
            as_index=False,
        )[
            "n_transitions"
        ]
        .sum()
        .rename(
            columns={
                "n_transitions":
                    "n_total"
            }
        )
    )

    pooled = pooled.merge(
        totals,
        on="alpha",
        how="left",
        validate="many_to_one",
    )

    pooled[
        "fraction"
    ] = (
        pooled[
            "n_transitions"
        ]
        / pooled[
            "n_total"
        ]
    )

    pooled[
        "alpha_tag"
    ] = pooled[
        "alpha"
    ].map(
        alpha_tag
    )

    return pooled.sort_values(
        [
            "alpha",
            "obs_class_step15",
        ],
        ascending=[
            False,
            True,
        ],
    ).reset_index(
        drop=True
    )


def tt_mask(
    frame: pd.DataFrame,
) -> np.ndarray:
    return (
        frame[
            "op_t0"
        ].to_numpy(bool)
        & frame[
            "op_t1"
        ].to_numpy(bool)
    )


def build_tt_retention(
    step15: ModuleType,
    base_df: pd.DataFrame,
    alphas: Sequence[float],
    reference_alpha: float,
) -> pd.DataFrame:
    reference = (
        step15
        .add_operational_columns(
            base_df,
            float(
                reference_alpha
            ),
        )
    )

    ref_tt = tt_mask(
        reference
    )

    n_ref = int(
        ref_tt.sum()
    )

    rows = []

    for alpha in alphas:
        classified = (
            step15
            .add_operational_columns(
                base_df,
                float(
                    alpha
                ),
            )
        )

        current_tt = tt_mask(
            classified
        )

        intersection = (
            current_tt
            & ref_tt
        )

        n_current = int(
            current_tt.sum()
        )
        n_intersection = int(
            intersection.sum()
        )

        rows.append(
            {
                "alpha":
                    float(
                        alpha
                    ),
                "alpha_tag":
                    alpha_tag(
                        alpha
                    ),
                "reference_alpha":
                    float(
                        reference_alpha
                    ),
                "n_tt":
                    n_current,
                "reference_n_tt":
                    n_ref,
                "n_tt_intersection_with_reference":
                    n_intersection,
                "fraction_reference_tt_retained":
                    (
                        n_intersection
                        / n_ref
                        if n_ref
                        else np.nan
                    ),
                "fraction_current_tt_shared_with_reference":
                    (
                        n_intersection
                        / n_current
                        if n_current
                        else np.nan
                    ),
            }
        )

    return pd.DataFrame(
        rows
    ).sort_values(
        "alpha",
        ascending=False,
    ).reset_index(
        drop=True
    )


# =============================================================================
# Cross-alpha paired forward differences
# =============================================================================

def combine_forward_curves(
    curves_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
) -> pd.DataFrame:
    frames = []

    for alpha, frame in curves_by_alpha.items():
        d = frame.copy()
        d[
            "alpha"
        ] = float(
            alpha
        )
        d[
            "alpha_tag"
        ] = alpha_tag(
            alpha
        )

        frames.append(
            d
        )

    return pd.concat(
        frames,
        ignore_index=True,
    ).sort_values(
        [
            "mode",
            "alpha",
            "bin",
        ],
        ascending=[
            True,
            False,
            True,
        ],
    ).reset_index(
        drop=True
    )


def combine_forward_slopes(
    slopes_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
) -> pd.DataFrame:
    frames = []

    for alpha, frame in slopes_by_alpha.items():
        d = frame.copy()
        d[
            "alpha"
        ] = float(
            alpha
        )
        d[
            "alpha_tag"
        ] = alpha_tag(
            alpha
        )

        frames.append(
            d
        )

    return pd.concat(
        frames,
        ignore_index=True,
    ).sort_values(
        [
            "mode",
            "alpha",
        ],
        ascending=[
            True,
            False,
        ],
    ).reset_index(
        drop=True
    )



def resolve_forward_cross_alpha_common_support(
    curves_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    forward_edges: pd.DataFrame,
) -> Tuple[
    List[int],
    pd.DataFrame,
]:
    """
    Define one forward-abundance support shared by all non-primary operational
    domains and all alpha values.

    The original Step-16 v1 slope summaries allowed each alpha/domain to use
    its own finite bin subset. At stringent alpha values this truncated the
    lower-abundance bins, so cross-alpha slope differences mixed threshold
    effects with changing x-support.

    Step 16 v2 fixes this by using the largest contiguous set of Step-9 bins
    for which T0PLUS, TT and TF all have finite point estimates at every alpha.

    PRIMARY is alpha invariant and remains available on its native Step-9 grid;
    threshold-sensitivity slopes for T0PLUS/TT/TF use this common support.
    """
    modes = list(
        FORWARD_SENSITIVITY_MODES
    )

    qualifying = []

    rows = []

    for _, edge in forward_edges.iterrows():
        bin_id = int(
            edge[
                "bin"
            ]
        )

        valid_all = True
        row = {
            "bin":
                bin_id,
            "x_left":
                float(
                    edge[
                        "x_left"
                    ]
                ),
            "x_right":
                float(
                    edge[
                        "x_right"
                    ]
                ),
            "x_center":
                float(
                    edge[
                        "x_center"
                    ]
                ),
        }

        for alpha in sorted(
            curves_by_alpha.keys(),
            reverse=True,
        ):
            frame = curves_by_alpha[
                alpha
            ]

            for mode in modes:
                g = frame[
                    frame[
                        "mode"
                    ].astype(str).eq(
                        mode
                    )
                    & frame[
                        "bin"
                    ].astype(int).eq(
                        bin_id
                    )
                ]

                finite = bool(
                    len(g) == 1
                    and pd.notna(
                        g.iloc[0][
                            "mean_dx"
                        ]
                    )
                    and np.isfinite(
                        float(
                            g.iloc[0][
                                "mean_dx"
                            ]
                        )
                    )
                )

                row[
                    f"{mode}_alpha_{alpha_tag(alpha)}_valid"
                ] = finite

                valid_all = (
                    valid_all
                    and finite
                )

        row[
            "valid_all_modes_all_alpha"
        ] = bool(
            valid_all
        )

        if valid_all:
            qualifying.append(
                bin_id
            )

        rows.append(
            row
        )

    selected = largest_contiguous_run(
        qualifying
    )

    support = pd.DataFrame(
        rows
    )

    support[
        "selected_cross_alpha_forward_common_support"
    ] = support[
        "bin"
    ].isin(
        selected
    )

    return (
        selected,
        support,
    )


def forward_slopes_on_common_support(
    step15: ModuleType,
    curves_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    boot_curves_by_alpha: Mapping[
        float,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    common_bins: Sequence[int],
) -> Tuple[
    Dict[
        float,
        pd.DataFrame,
    ],
    Dict[
        float,
        Dict[
            str,
            np.ndarray,
        ],
    ],
]:
    """
    Refit forward abundance slopes on one common bin set across alpha/domain.

    The same bootstrap draw index is already aligned across alpha, so paired
    slope differences remain valid after restricting every curve to the same
    x-support.
    """
    if not common_bins:
        raise ValueError(
            "Cross-alpha forward common support is empty."
        )

    support_bins = [
        int(x)
        for x in common_bins
    ]

    slopes_by_alpha: Dict[
        float,
        pd.DataFrame,
    ] = {}

    boot_slopes_by_alpha: Dict[
        float,
        Dict[
            str,
            np.ndarray,
        ],
    ] = {}

    for alpha in sorted(
        curves_by_alpha.keys(),
        reverse=True,
    ):
        frame = curves_by_alpha[
            alpha
        ]

        boot_curves = boot_curves_by_alpha[
            alpha
        ]

        rows = []
        boot_dict: Dict[
            str,
            np.ndarray,
        ] = {}

        for mode in FORWARD_SENSITIVITY_MODES:
            g = (
                frame[
                    frame[
                        "mode"
                    ].astype(str).eq(
                        mode
                    )
                    & frame[
                        "bin"
                    ].astype(int).isin(
                        support_bins
                    )
                ]
                .sort_values(
                    "bin"
                )
            )

            if len(g) != len(
                support_bins
            ):
                raise ValueError(
                    f"Forward common-support row mismatch for alpha={alpha}, "
                    f"mode={mode}: {len(g)} vs {len(support_bins)}."
                )

            x = g[
                "x_center"
            ].to_numpy(
                dtype=float
            )
            y = g[
                "mean_dx"
            ].to_numpy(
                dtype=float
            )

            if not np.isfinite(y).all():
                raise ValueError(
                    f"Forward common-support point curve incomplete for "
                    f"alpha={alpha}, mode={mode}."
                )

            point_fit = (
                step15
                .safe_linear_fit(
                    x,
                    y,
                )
            )

            n_boot = (
                boot_curves[
                    mode
                ].shape[0]
            )

            boot_slopes = np.full(
                n_boot,
                np.nan,
                dtype=float,
            )

            for b in range(
                n_boot
            ):
                y_boot = (
                    boot_curves[
                        mode
                    ][
                        b,
                        support_bins,
                    ]
                )

                # Cross-alpha slope inference uses one fixed x-support.
                # A draw with any missing common-support bin does not
                # contribute a slope.
                if not np.isfinite(
                    y_boot
                ).all():
                    continue

                fit = (
                    step15
                    .safe_linear_fit(
                        x,
                        y_boot,
                    )
                )

                if fit["n_points"] == len(
                    support_bins
                ):
                    boot_slopes[
                        b
                    ] = fit[
                        "slope"
                    ]

            boot_dict[
                mode
            ] = boot_slopes

            lo, med, hi, n_valid = (
                percentile_summary(
                    boot_slopes
                )
            )

            rows.append(
                {
                    "mode":
                        mode,
                    "support":
                        "cross_alpha_forward_common_support",
                    "support_bins":
                        ",".join(
                            str(x)
                            for x
                            in support_bins
                        ),
                    "n_valid_bins":
                        int(
                            len(
                                support_bins
                            )
                        ),
                    "x_min":
                        float(
                            g[
                                "x_left"
                            ].min()
                        ),
                    "x_max":
                        float(
                            g[
                                "x_right"
                            ].max()
                        ),
                    "observed_slope":
                        point_fit[
                            "slope"
                        ],
                    "observed_intercept":
                        point_fit[
                            "intercept"
                        ],
                    "bootstrap_median_slope":
                        med,
                    "bootstrap_q025":
                        lo,
                    "bootstrap_q975":
                        hi,
                    "n_bootstrap_valid":
                        n_valid,
                    "bootstrap_positive_slope_fraction":
                        (
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

        slopes = pd.DataFrame(
            rows
        )

        # Paired difference of each common-support mode vs the invariant
        # manuscript PRIMARY slope evaluated on the same common bin set.
        primary = (
            frame[
                frame[
                    "mode"
                ].astype(str).eq(
                    "PRIMARY"
                )
                & frame[
                    "bin"
                ].astype(int).isin(
                    support_bins
                )
            ]
            .sort_values(
                "bin"
            )
        )

        x_primary = primary[
            "x_center"
        ].to_numpy(
            dtype=float
        )
        y_primary = primary[
            "mean_dx"
        ].to_numpy(
            dtype=float
        )

        primary_fit = (
            step15
            .safe_linear_fit(
                x_primary,
                y_primary,
            )
        )

        primary_boot_slopes = np.full(
            boot_curves[
                "PRIMARY"
            ].shape[0],
            np.nan,
            dtype=float,
        )

        for b in range(
            primary_boot_slopes.size
        ):
            y_boot = (
                boot_curves[
                    "PRIMARY"
                ][
                    b,
                    support_bins,
                ]
            )

            if not np.isfinite(
                y_boot
            ).all():
                continue

            fit = (
                step15
                .safe_linear_fit(
                    x_primary,
                    y_boot,
                )
            )
            if fit["n_points"] == len(
                support_bins
            ):
                primary_boot_slopes[
                    b
                ] = fit[
                    "slope"
                ]

        for mode in FORWARD_SENSITIVITY_MODES:
            mask = slopes[
                "mode"
            ].eq(
                mode
            )

            observed_mode = float(
                slopes.loc[
                    mask,
                    "observed_slope",
                ].iloc[0]
            )

            diff = (
                boot_dict[
                    mode
                ]
                - primary_boot_slopes
            )

            lo, med, hi, n_valid = (
                percentile_summary(
                    diff
                )
            )

            slopes.loc[
                mask,
                "observed_slope_difference_vs_PRIMARY_same_support",
            ] = (
                observed_mode
                - float(
                    primary_fit[
                        "slope"
                    ]
                )
            )
            slopes.loc[
                mask,
                "bootstrap_slope_difference_vs_PRIMARY_same_support_median",
            ] = med
            slopes.loc[
                mask,
                "bootstrap_slope_difference_vs_PRIMARY_same_support_q025",
            ] = lo
            slopes.loc[
                mask,
                "bootstrap_slope_difference_vs_PRIMARY_same_support_q975",
            ] = hi
            slopes.loc[
                mask,
                "n_bootstrap_slope_difference_vs_PRIMARY_same_support_valid",
            ] = n_valid

        slopes[
            "primary_same_support_observed_slope"
        ] = float(
            primary_fit[
                "slope"
            ]
        )

        slopes[
            "alpha"
        ] = float(
            alpha
        )
        slopes[
            "alpha_tag"
        ] = alpha_tag(
            alpha
        )

        slopes_by_alpha[
            float(
                alpha
            )
        ] = slopes

        boot_slopes_by_alpha[
            float(
                alpha
            )
        ] = boot_dict

    return (
        slopes_by_alpha,
        boot_slopes_by_alpha,
    )


def forward_curve_difference_vs_reference(
    curves_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    boot_curves_by_alpha: Mapping[
        float,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    reference_alpha: float,
) -> pd.DataFrame:
    ref_curves = curves_by_alpha[
        float(
            reference_alpha
        )
    ].set_index(
        [
            "mode",
            "bin",
        ]
    )

    ref_boot = boot_curves_by_alpha[
        float(
            reference_alpha
        )
    ]

    rows = []

    for alpha in sorted(
        curves_by_alpha.keys(),
        reverse=True,
    ):
        if math.isclose(
            float(alpha),
            float(
                reference_alpha
            ),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            continue

        current = curves_by_alpha[
            alpha
        ].set_index(
            [
                "mode",
                "bin",
            ]
        )

        for mode in FORWARD_SENSITIVITY_MODES:
            for (
                mode_key,
                bin_id,
            ), row in current.loc[
                mode
            ].reset_index().assign(
                mode=mode
            ).set_index(
                [
                    "mode",
                    "bin",
                ]
            ).iterrows():
                key = (
                    mode_key,
                    int(
                        bin_id
                    ),
                )

                if key not in ref_curves.index:
                    continue

                ref_row = ref_curves.loc[
                    key
                ]

                value = (
                    float(
                        row[
                            "mean_dx"
                        ]
                    )
                    if pd.notna(
                        row[
                            "mean_dx"
                        ]
                    )
                    else np.nan
                )

                ref_value = (
                    float(
                        ref_row[
                            "mean_dx"
                        ]
                    )
                    if pd.notna(
                        ref_row[
                            "mean_dx"
                        ]
                    )
                    else np.nan
                )

                current_boot = (
                    boot_curves_by_alpha[
                        alpha
                    ][
                        mode
                    ][:, int(
                        bin_id
                    )]
                )

                reference_boot = (
                    ref_boot[
                        mode
                    ][:, int(
                        bin_id
                    )]
                )

                difference_boot = (
                    current_boot
                    - reference_boot
                )

                lo, med, hi, n_valid = (
                    percentile_summary(
                        difference_boot
                    )
                )

                rows.append(
                    {
                        "alpha":
                            float(
                                alpha
                            ),
                        "alpha_tag":
                            alpha_tag(
                                alpha
                            ),
                        "reference_alpha":
                            float(
                                reference_alpha
                            ),
                        "mode":
                            mode,
                        "bin":
                            int(
                                bin_id
                            ),
                        "x_left":
                            float(
                                row[
                                    "x_left"
                                ]
                            ),
                        "x_right":
                            float(
                                row[
                                    "x_right"
                                ]
                            ),
                        "x_center":
                            float(
                                row[
                                    "x_center"
                                ]
                            ),
                        "value":
                            value,
                        "reference_value":
                            ref_value,
                        "difference":
                            (
                                value
                                - ref_value
                                if np.isfinite(
                                    value
                                )
                                and np.isfinite(
                                    ref_value
                                )
                                else np.nan
                            ),
                        "bootstrap_difference_median":
                            med,
                        "bootstrap_difference_q025":
                            lo,
                        "bootstrap_difference_q975":
                            hi,
                        "n_bootstrap_difference_valid":
                            n_valid,
                    }
                )

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "mode",
            "alpha",
            "bin",
        ],
        ascending=[
            True,
            False,
            True,
        ],
    ).reset_index(
        drop=True
    )


def forward_slope_difference_vs_reference(
    slope_summary_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    boot_slope_arrays_by_alpha: Mapping[
        float,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    reference_alpha: float,
) -> pd.DataFrame:
    ref_summary = (
        slope_summary_by_alpha[
            float(
                reference_alpha
            )
        ]
        .set_index(
            "mode"
        )
    )

    ref_boot = (
        boot_slope_arrays_by_alpha[
            float(
                reference_alpha
            )
        ]
    )

    rows = []

    for alpha in sorted(
        slope_summary_by_alpha.keys(),
        reverse=True,
    ):
        if math.isclose(
            float(
                alpha
            ),
            float(
                reference_alpha
            ),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            continue

        current = (
            slope_summary_by_alpha[
                alpha
            ]
            .set_index(
                "mode"
            )
        )

        for mode in FORWARD_SENSITIVITY_MODES:
            if (
                mode
                not in current.index
                or mode
                not in ref_summary.index
            ):
                continue

            observed = float(
                current.loc[
                    mode,
                    "observed_slope",
                ]
            )

            reference = float(
                ref_summary.loc[
                    mode,
                    "observed_slope",
                ]
            )

            diff_boot = (
                boot_slope_arrays_by_alpha[
                    alpha
                ][
                    mode
                ]
                - ref_boot[
                    mode
                ]
            )

            lo, med, hi, n_valid = (
                percentile_summary(
                    diff_boot
                )
            )

            rows.append(
                {
                    "alpha":
                        float(
                            alpha
                        ),
                    "alpha_tag":
                        alpha_tag(
                            alpha
                        ),
                    "reference_alpha":
                        float(
                            reference_alpha
                        ),
                    "mode":
                        mode,
                    "observed_slope":
                        observed,
                    "reference_observed_slope":
                        reference,
                    "observed_slope_difference":
                        observed
                        - reference,
                    "bootstrap_slope_difference_median":
                        med,
                    "bootstrap_slope_difference_q025":
                        lo,
                    "bootstrap_slope_difference_q975":
                        hi,
                    "n_bootstrap_difference_valid":
                        n_valid,
                }
            )

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "mode",
            "alpha",
        ],
        ascending=[
            True,
            False,
        ],
    ).reset_index(
        drop=True
    )


# =============================================================================
# Cross-alpha fluctuation differences
# =============================================================================

def combine_fluctuation_points(
    points_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
) -> pd.DataFrame:
    frames = []

    for alpha, frame in points_by_alpha.items():
        d = frame.copy()
        d[
            "alpha"
        ] = float(
            alpha
        )
        d[
            "alpha_tag"
        ] = alpha_tag(
            alpha
        )

        frames.append(
            d
        )

    return pd.concat(
        frames,
        ignore_index=True,
    ).sort_values(
        [
            "mode",
            "alpha",
            "dt",
            "bin",
        ],
        ascending=[
            True,
            False,
            True,
            True,
        ],
    ).reset_index(
        drop=True
    )


def fluctuation_difference_vs_reference(
    points_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    boot_by_alpha: Mapping[
        float,
        Mapping[
            str,
            Mapping[
                str,
                np.ndarray,
            ],
        ],
    ],
    dts: Sequence[int],
    fluct_edges: pd.DataFrame,
    reference_alpha: float,
) -> pd.DataFrame:
    ref_point = (
        points_by_alpha[
            float(
                reference_alpha
            )
        ]
        .set_index(
            [
                "mode",
                "dt",
                "bin",
            ]
        )
    )

    ref_boot = (
        boot_by_alpha[
            float(
                reference_alpha
            )
        ]
    )

    rows = []

    for alpha in sorted(
        points_by_alpha.keys(),
        reverse=True,
    ):
        if math.isclose(
            float(
                alpha
            ),
            float(
                reference_alpha
            ),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            continue

        current = (
            points_by_alpha[
                alpha
            ]
            .set_index(
                [
                    "mode",
                    "dt",
                    "bin",
                ]
            )
        )

        for mode in FLUCTUATION_SENSITIVITY_MODES:
            for j, dt in enumerate(
                dts
            ):
                for _, edge in fluct_edges.iterrows():
                    bin_id = int(
                        edge[
                            "bin"
                        ]
                    )

                    key = (
                        mode,
                        int(
                            dt
                        ),
                        bin_id,
                    )

                    if (
                        key
                        not in current.index
                        or key
                        not in ref_point.index
                    ):
                        continue

                    cur_row = current.loc[
                        key
                    ]
                    ref_row = ref_point.loc[
                        key
                    ]

                    for metric in METRICS:
                        value = (
                            float(
                                cur_row[
                                    metric
                                ]
                            )
                            if pd.notna(
                                cur_row[
                                    metric
                                ]
                            )
                            else np.nan
                        )

                        ref_value = (
                            float(
                                ref_row[
                                    metric
                                ]
                            )
                            if pd.notna(
                                ref_row[
                                    metric
                                ]
                            )
                            else np.nan
                        )

                        diff_boot = (
                            boot_by_alpha[
                                alpha
                            ][
                                mode
                            ][
                                metric
                            ][
                                :,
                                j,
                                bin_id,
                            ]
                            - ref_boot[
                                mode
                            ][
                                metric
                            ][
                                :,
                                j,
                                bin_id,
                            ]
                        )

                        lo, med, hi, n_valid = (
                            percentile_summary(
                                diff_boot
                            )
                        )

                        rows.append(
                            {
                                "alpha":
                                    float(
                                        alpha
                                    ),
                                "alpha_tag":
                                    alpha_tag(
                                        alpha
                                    ),
                                "reference_alpha":
                                    float(
                                        reference_alpha
                                    ),
                                "mode":
                                    mode,
                                "metric":
                                    metric,
                                "dt":
                                    int(
                                        dt
                                    ),
                                "bin":
                                    bin_id,
                                "x_left":
                                    float(
                                        edge[
                                            "x_left"
                                        ]
                                    ),
                                "x_right":
                                    float(
                                        edge[
                                            "x_right"
                                        ]
                                    ),
                                "x_center":
                                    float(
                                        edge[
                                            "x_center"
                                        ]
                                    ),
                                "value":
                                    value,
                                "reference_value":
                                    ref_value,
                                "difference":
                                    (
                                        value
                                        - ref_value
                                        if np.isfinite(
                                            value
                                        )
                                        and np.isfinite(
                                            ref_value
                                        )
                                        else np.nan
                                    ),
                                "bootstrap_difference_median":
                                    med,
                                "bootstrap_difference_q025":
                                    lo,
                                "bootstrap_difference_q975":
                                    hi,
                                "n_bootstrap_difference_valid":
                                    n_valid,
                            }
                        )

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "mode",
            "metric",
            "alpha",
            "dt",
            "bin",
        ],
        ascending=[
            True,
            True,
            False,
            True,
            True,
        ],
    ).reset_index(
        drop=True
    )


# =============================================================================
# Cross-alpha common temporal core and temporal bootstrap
# =============================================================================

def build_cross_alpha_temporal_core(
    support_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    frozen_core: pd.DataFrame,
) -> Tuple[
    List[int],
    pd.DataFrame,
]:
    frames = []

    valid_sets = []

    for alpha, support in support_by_alpha.items():
        d = support.copy()
        d[
            "alpha"
        ] = float(
            alpha
        )
        d[
            "alpha_tag"
        ] = alpha_tag(
            alpha
        )

        frames.append(
            d
        )

        valid = set(
            int(x)
            for x in d.loc[
                d[
                    "valid_all_modes_all_dt"
                ].astype(bool),
                "bin",
            ].tolist()
        )

        valid_sets.append(
            valid
        )

    common = (
        set.intersection(
            *valid_sets
        )
        if valid_sets
        else set()
    )

    selected = (
        largest_contiguous_run(
            common
        )
    )

    combined = pd.concat(
        frames,
        ignore_index=True,
    )

    combined[
        "selected_cross_alpha_common_core"
    ] = combined[
        "bin"
    ].isin(
        selected
    )

    geometry = frozen_core[
        frozen_core[
            "bin"
        ].isin(
            selected
        )
    ][
        [
            "bin",
            "x_left",
            "x_right",
            "x_center",
        ]
    ].copy()

    return (
        selected,
        combined,
    )


def temporal_bootstrap_for_alpha(
    step15: ModuleType,
    fluct_arrays: Mapping[
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
    Dict[
        str,
        np.ndarray,
    ],
]:
    if not common_core_bins:
        raise ValueError(
            "Cross-alpha common temporal core is empty."
        )

    subject_index = {
        int(subject):
            i
        for i, subject
        in enumerate(
            all_subjects
        )
    }

    idx = [
        subject_index[
            int(subject)
        ]
        for subject
        in temporal_subjects
    ]

    n_subjects = len(
        temporal_subjects
    )

    rng = np.random.default_rng(
        int(
            seed
        )
        + 15001
    )

    draws = rng.integers(
        0,
        n_subjects,
        size=(
            int(
                n_bootstrap
            ),
            n_subjects,
        ),
    )

    dt_array = np.asarray(
        dts,
        dtype=float,
    )

    curve_rows = []
    slope_rows = []
    boot_slopes: Dict[
        str,
        np.ndarray,
    ] = {}

    for mode in TEMPORAL_MODES:
        mode_arrays = {
            key:
                value[
                    idx,
                    :,
                    :,
                ]
            for key, value
            in fluct_arrays[
                mode
            ].items()
        }

        point = (
            step15
            .aggregate_fluctuation(
                mode_arrays,
                np.ones(
                    n_subjects,
                    dtype=float,
                ),
            )
        )

        selected = list(
            common_core_bins
        )

        point_cross = point[
            "cross_cov"
        ][
            :,
            selected,
        ]

        point_n = point[
            "n"
        ][
            :,
            selected,
        ]

        point_subject_support = (
            mode_arrays[
                "n"
            ][
                :,
                :,
                selected,
            ]
            > 0
        ).sum(
            axis=0
        )

        point_valid = (
            np.isfinite(
                point_cross
            )
            & (
                point_n
                >= int(
                    min_n
                )
            )
            & (
                point_subject_support
                >= int(
                    min_subjects
                )
            )
        )

        if not bool(
            np.all(
                point_valid
            )
        ):
            raise RuntimeError(
                f"{mode}: cross-alpha temporal common core is not fully "
                "supported in the point estimate."
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
                    "mode":
                        mode,
                    "dt":
                        int(
                            dt
                        ),
                    "cross_cov":
                        float(
                            point_curve[
                                j
                            ]
                        ),
                    "n_common_core_bins":
                        int(
                            len(
                                common_core_bins
                            )
                        ),
                    "n_subjects_universe":
                        int(
                            n_subjects
                        ),
                }
            )

        fit = (
            step15
            .safe_linear_fit(
                dt_array
                - 1.0,
                point_curve,
            )
        )

        slopes = np.full(
            int(
                n_bootstrap
            ),
            np.nan,
            dtype=float,
        )

        for b in range(
            int(
                n_bootstrap
            )
        ):
            multiplicity = np.bincount(
                draws[
                    b
                ],
                minlength=n_subjects,
            ).astype(
                float
            )

            agg = (
                step15
                .aggregate_fluctuation(
                    mode_arrays,
                    multiplicity,
                )
            )

            sampled_subject_support = (
                (
                    mode_arrays[
                        "n"
                    ]
                    > 0
                )
                & (
                    multiplicity[
                        :,
                        None,
                        None,
                    ]
                    > 0
                )
            ).sum(
                axis=0
            )

            selected_cross = agg[
                "cross_cov"
            ][
                :,
                selected,
            ]

            selected_n = agg[
                "n"
            ][
                :,
                selected,
            ]

            selected_subjects = (
                sampled_subject_support[
                    :,
                    selected,
                ]
            )

            valid = (
                np.isfinite(
                    selected_cross
                )
                & (
                    selected_n
                    >= int(
                        min_n
                    )
                )
                & (
                    selected_subjects
                    >= int(
                        min_subjects
                    )
                )
            )

            # Fixed cross-alpha common core: a bootstrap draw contributes only
            # if every selected bin remains valid at every requested lag.
            if not bool(
                np.all(
                    valid
                )
            ):
                continue

            profile = np.mean(
                selected_cross,
                axis=1,
            )

            if not np.isfinite(
                profile
            ).all():
                continue

            boot_fit = (
                step15
                .safe_linear_fit(
                    dt_array
                    - 1.0,
                    profile,
                )
            )

            if boot_fit[
                "n_points"
            ] == len(
                dts
            ):
                slopes[
                    b
                ] = boot_fit[
                    "slope"
                ]

        boot_slopes[
            mode
        ] = slopes

        lo, med, hi, n_valid = (
            percentile_summary(
                slopes
            )
        )

        slope_rows.append(
            {
                "mode":
                    mode,
                "observed_slope_per_week":
                    fit[
                        "slope"
                    ],
                "observed_intercept_at_dt1":
                    fit[
                        "intercept"
                    ],
                "n_dt":
                    fit[
                        "n_points"
                    ],
                "n_common_core_bins":
                    int(
                        len(
                            common_core_bins
                        )
                    ),
                "n_subjects_universe":
                    int(
                        n_subjects
                    ),
                "bootstrap_median_slope":
                    med,
                "bootstrap_q025":
                    lo,
                "bootstrap_q975":
                    hi,
                "n_bootstrap_valid":
                    n_valid,
                "bootstrap_positive_slope_fraction":
                    (
                        float(
                            np.mean(
                                finite_numeric(
                                    slopes
                                )
                                > 0
                            )
                        )
                        if n_valid
                        else np.nan
                    ),
            }
        )

    return (
        pd.DataFrame(
            curve_rows
        ),
        pd.DataFrame(
            slope_rows
        ),
        boot_slopes,
    )


def combine_temporal_profiles(
    profiles_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
) -> pd.DataFrame:
    frames = []

    for alpha, frame in profiles_by_alpha.items():
        d = frame.copy()
        d[
            "alpha"
        ] = float(
            alpha
        )
        d[
            "alpha_tag"
        ] = alpha_tag(
            alpha
        )

        frames.append(
            d
        )

    return pd.concat(
        frames,
        ignore_index=True,
    ).sort_values(
        [
            "mode",
            "alpha",
            "dt",
        ],
        ascending=[
            True,
            False,
            True,
        ],
    ).reset_index(
        drop=True
    )


def combine_temporal_slopes(
    slopes_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
) -> pd.DataFrame:
    frames = []

    for alpha, frame in slopes_by_alpha.items():
        d = frame.copy()
        d[
            "alpha"
        ] = float(
            alpha
        )
        d[
            "alpha_tag"
        ] = alpha_tag(
            alpha
        )

        frames.append(
            d
        )

    return pd.concat(
        frames,
        ignore_index=True,
    ).sort_values(
        [
            "mode",
            "alpha",
        ],
        ascending=[
            True,
            False,
        ],
    ).reset_index(
        drop=True
    )


def temporal_slope_difference_vs_reference(
    slopes_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    boot_slopes_by_alpha: Mapping[
        float,
        Mapping[
            str,
            np.ndarray,
        ],
    ],
    reference_alpha: float,
) -> pd.DataFrame:
    reference = (
        slopes_by_alpha[
            float(
                reference_alpha
            )
        ]
        .set_index(
            "mode"
        )
    )

    reference_boot = (
        boot_slopes_by_alpha[
            float(
                reference_alpha
            )
        ]
    )

    rows = []

    for alpha in sorted(
        slopes_by_alpha.keys(),
        reverse=True,
    ):
        if math.isclose(
            float(
                alpha
            ),
            float(
                reference_alpha
            ),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            continue

        current = (
            slopes_by_alpha[
                alpha
            ]
            .set_index(
                "mode"
            )
        )

        for mode in TEMPORAL_MODES:
            observed = float(
                current.loc[
                    mode,
                    "observed_slope_per_week",
                ]
            )

            ref_observed = float(
                reference.loc[
                    mode,
                    "observed_slope_per_week",
                ]
            )

            diff_boot = (
                boot_slopes_by_alpha[
                    alpha
                ][
                    mode
                ]
                - reference_boot[
                    mode
                ]
            )

            lo, med, hi, n_valid = (
                percentile_summary(
                    diff_boot
                )
            )

            rows.append(
                {
                    "alpha":
                        float(
                            alpha
                        ),
                    "alpha_tag":
                        alpha_tag(
                            alpha
                        ),
                    "reference_alpha":
                        float(
                            reference_alpha
                        ),
                    "mode":
                        mode,
                    "observed_slope":
                        observed,
                    "reference_observed_slope":
                        ref_observed,
                    "observed_slope_difference":
                        observed
                        - ref_observed,
                    "bootstrap_slope_difference_median":
                        med,
                    "bootstrap_slope_difference_q025":
                        lo,
                    "bootstrap_slope_difference_q975":
                        hi,
                    "n_bootstrap_difference_valid":
                        n_valid,
                }
            )

    return pd.DataFrame(
        rows
    ).sort_values(
        [
            "mode",
            "alpha",
        ],
        ascending=[
            True,
            False,
        ],
    ).reset_index(
        drop=True
    )


# =============================================================================
# Audits and summaries
# =============================================================================

def max_abs_numeric(
    frame: pd.DataFrame,
    columns: Sequence[str],
) -> float:
    values = []

    for column in columns:
        if column not in frame.columns:
            continue

        values.extend(
            finite_numeric(
                frame[
                    column
                ]
            ).tolist()
        )

    return (
        float(
            np.max(
                np.abs(
                    values
                )
            )
        )
        if values
        else np.nan
    )


def build_reference_internal_audit(
    forward_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    fluct_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    reference_alpha: float,
) -> pd.DataFrame:
    rows = []

    ref_forward = forward_by_alpha[
        float(
            reference_alpha
        )
    ]

    primary_forward = ref_forward[
        ref_forward[
            "mode"
        ].astype(str).eq(
            "PRIMARY"
        )
    ]

    rows.append(
        {
            "audit":
                "reference_alpha_forward_primary_present",
            "value":
                int(
                    len(
                        primary_forward
                    )
                ),
            "expected":
                20,
            "absolute_difference":
                abs(
                    len(
                        primary_forward
                    )
                    - 20
                ),
        }
    )

    ref_fluct = fluct_by_alpha[
        float(
            reference_alpha
        )
    ]

    primary_fluct = ref_fluct[
        ref_fluct[
            "mode"
        ].astype(str).eq(
            "PRIMARY"
        )
    ]

    rows.append(
        {
            "audit":
                "reference_alpha_fluctuation_primary_rows",
            "value":
                int(
                    len(
                        primary_fluct
                    )
                ),
            "expected":
                int(
                    len(
                        primary_fluct
                    )
                ),
            "absolute_difference":
                0,
        }
    )

    return pd.DataFrame(
        rows
    )


def audit_against_step15_reference_outputs(
    step15_reference_dir: Optional[Path],
    reference_alpha: float,
    forward: pd.DataFrame,
    fluctuation: pd.DataFrame,
    temporal_slopes: pd.DataFrame,
) -> pd.DataFrame:
    if step15_reference_dir is None:
        return pd.DataFrame()

    rows = []

    # Forward.
    path = (
        step15_reference_dir
        / "03_forward_domain_by_bin.csv"
    )
    if path.is_file():
        ref = pd.read_csv(
            path
        )

        joined = (
            forward.merge(
                ref,
                on=[
                    "mode",
                    "bin",
                ],
                how="inner",
                suffixes=(
                    "_step16",
                    "_step15",
                ),
            )
        )

        diff = (
            pd.to_numeric(
                joined[
                    "mean_dx_step16"
                ],
                errors="coerce",
            )
            - pd.to_numeric(
                joined[
                    "mean_dx_step15"
                ],
                errors="coerce",
            )
        )

        rows.append(
            {
                "component":
                    "forward_domain_by_bin",
                "n_compared":
                    int(
                        np.isfinite(
                            diff
                        ).sum()
                    ),
                "max_abs_difference":
                    float(
                        np.nanmax(
                            np.abs(
                                diff
                            )
                        )
                    ),
            }
        )

    # Fluctuation.
    path = (
        step15_reference_dir
        / "07_fluctuation_domain_by_bin_dt.csv"
    )
    if path.is_file():
        ref = pd.read_csv(
            path
        )

        joined = (
            fluctuation.merge(
                ref,
                on=[
                    "mode",
                    "dt",
                    "bin",
                ],
                how="inner",
                suffixes=(
                    "_step16",
                    "_step15",
                ),
            )
        )

        diffs = []

        for metric in METRICS:
            diff = (
                pd.to_numeric(
                    joined[
                        f"{metric}_step16"
                    ],
                    errors="coerce",
                )
                - pd.to_numeric(
                    joined[
                        f"{metric}_step15"
                    ],
                    errors="coerce",
                )
            )
            diffs.extend(
                finite_numeric(
                    diff
                ).tolist()
            )

        rows.append(
            {
                "component":
                    "fluctuation_domain_by_bin_dt",
                "n_compared":
                    int(
                        len(
                            diffs
                        )
                    ),
                "max_abs_difference":
                    (
                        float(
                            np.max(
                                np.abs(
                                    diffs
                                )
                            )
                        )
                        if diffs
                        else np.nan
                    ),
            }
        )

    # Temporal slopes are expected to differ if Step 16 uses a narrower
    # cross-alpha common core, so this is reported but not treated as identity.
    path = (
        step15_reference_dir
        / "13_temporal_domain_slope_summary.csv"
    )
    if path.is_file():
        ref = pd.read_csv(
            path
        )

        joined = (
            temporal_slopes.merge(
                ref[
                    [
                        "mode",
                        "observed_slope_per_week",
                    ]
                ],
                on="mode",
                how="inner",
                suffixes=(
                    "_step16_common_core",
                    "_step15_threshold_specific_core",
                ),
            )
        )

        for row in joined.itertuples(
            index=False
        ):
            rows.append(
                {
                    "component":
                        (
                            "temporal_slope_"
                            + str(
                                row.mode
                            )
                        ),
                    "n_compared":
                        1,
                    "max_abs_difference":
                        abs(
                            float(
                                row.observed_slope_per_week_step16_common_core
                            )
                            - float(
                                row.observed_slope_per_week_step15_threshold_specific_core
                            )
                        ),
                    "note":
                        "not an identity audit: Step 16 uses the cross-alpha common core",
                }
            )

    return pd.DataFrame(
        rows
    )


def support_summary(
    class_table: pd.DataFrame,
    support_by_alpha: Mapping[
        float,
        pd.DataFrame,
    ],
    cross_alpha_core: Sequence[int],
) -> pd.DataFrame:
    rows = []

    for alpha in sorted(
        support_by_alpha.keys(),
        reverse=True,
    ):
        class_alpha = (
            class_table[
                np.isclose(
                    class_table[
                        "alpha"
                    ].to_numpy(float),
                    float(
                        alpha
                    ),
                    rtol=0.0,
                    atol=1e-12,
                )
            ]
        )

        counts = (
            class_alpha.groupby(
                "obs_class_step15"
            )[
                "n_transitions"
            ]
            .sum()
            .to_dict()
        )

        support = support_by_alpha[
            alpha
        ]

        rows.append(
            {
                "alpha":
                    float(
                        alpha
                    ),
                "alpha_tag":
                    alpha_tag(
                        alpha
                    ),
                "n_TT":
                    int(
                        counts.get(
                            "TT",
                            0,
                        )
                    ),
                "n_TF":
                    int(
                        counts.get(
                            "TF",
                            0,
                        )
                    ),
                "n_FT":
                    int(
                        counts.get(
                            "FT",
                            0,
                        )
                    ),
                "n_FF":
                    int(
                        counts.get(
                            "FF",
                            0,
                        )
                    ),
                "n_frozen_core_bins_valid_all_modes_all_dt":
                    int(
                        support[
                            "valid_all_modes_all_dt"
                        ].astype(bool).sum()
                    ),
                "n_cross_alpha_common_core_bins":
                    int(
                        len(
                            cross_alpha_core
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    ).sort_values(
        "alpha",
        ascending=False,
    ).reset_index(
        drop=True
    )


def build_analysis_summary(
    pooled_composition: pd.DataFrame,
    forward_slopes: pd.DataFrame,
    temporal_slopes: pd.DataFrame,
    forward_slope_diff: pd.DataFrame,
    temporal_slope_diff: pd.DataFrame,
    cross_alpha_core: Sequence[int],
) -> pd.DataFrame:
    rows = []

    for row in pooled_composition.itertuples(
        index=False
    ):
        rows.append(
            {
                "section":
                    "class_composition",
                "alpha":
                    float(
                        row.alpha
                    ),
                "domain":
                    str(
                        row.obs_class_step15
                    ),
                "metric":
                    "fraction",
                "value":
                    float(
                        row.fraction
                    ),
                "lower":
                    np.nan,
                "upper":
                    np.nan,
                "note":
                    "pooled across temporal lags",
            }
        )

    for row in forward_slopes.itertuples(
        index=False
    ):
        if str(
            row.mode
        ) not in FORWARD_SENSITIVITY_MODES:
            continue

        rows.append(
            {
                "section":
                    "forward_slope",
                "alpha":
                    float(
                        row.alpha
                    ),
                "domain":
                    str(
                        row.mode
                    ),
                "metric":
                    "slope_vs_abundance",
                "value":
                    float(
                        row.observed_slope
                    ),
                "lower":
                    float(
                        row.bootstrap_q025
                    ),
                "upper":
                    float(
                        row.bootstrap_q975
                    ),
                "note":
                    "observed replicate-decoupled forward drift",
            }
        )

    for row in forward_slope_diff.itertuples(
        index=False
    ):
        rows.append(
            {
                "section":
                    "forward_slope_difference_vs_reference_alpha",
                "alpha":
                    float(
                        row.alpha
                    ),
                "domain":
                    str(
                        row.mode
                    ),
                "metric":
                    "paired_slope_difference",
                "value":
                    float(
                        row.bootstrap_slope_difference_median
                    ),
                "lower":
                    float(
                        row.bootstrap_slope_difference_q025
                    ),
                "upper":
                    float(
                        row.bootstrap_slope_difference_q975
                    ),
                "note":
                    "paired biological-subject bootstrap vs alpha=0.05",
            }
        )

    for row in temporal_slopes.itertuples(
        index=False
    ):
        rows.append(
            {
                "section":
                    "temporal_slope",
                "alpha":
                    float(
                        row.alpha
                    ),
                "domain":
                    str(
                        row.mode
                    ),
                "metric":
                    "cross_cov_slope_per_week",
                "value":
                    float(
                        row.observed_slope_per_week
                    ),
                "lower":
                    float(
                        row.bootstrap_q025
                    ),
                "upper":
                    float(
                        row.bootstrap_q975
                    ),
                "note":
                    (
                        f"cross-alpha common operational core; "
                        f"{len(cross_alpha_core)} bins"
                    ),
            }
        )

    for row in temporal_slope_diff.itertuples(
        index=False
    ):
        rows.append(
            {
                "section":
                    "temporal_slope_difference_vs_reference_alpha",
                "alpha":
                    float(
                        row.alpha
                    ),
                "domain":
                    str(
                        row.mode
                    ),
                "metric":
                    "paired_temporal_slope_difference",
                "value":
                    float(
                        row.bootstrap_slope_difference_median
                    ),
                "lower":
                    float(
                        row.bootstrap_slope_difference_q025
                    ),
                "upper":
                    float(
                        row.bootstrap_slope_difference_q975
                    ),
                "note":
                    "paired biological-subject bootstrap vs alpha=0.05",
            }
        )

    return pd.DataFrame(
        rows
    )


def write_readme(
    outdir: Path,
) -> None:
    text = """# Step 16 — operational observation-threshold robustness

Scientific scope
----------------
Step 16 varies only the numerical operational threshold:

    alpha = 0.10, 0.05, 0.025, 0.01

with alpha=0.05 as the reference.

T(alpha) = 1[p_endpoint < alpha]

T/F is an operational observation state and must not be interpreted as
biological presence/absence.

Current estimands
-----------------
Forward:
observed replicate-decoupled AB/BA estimator inherited from Step 9.

Fluctuations:
cross_cov = Cov(dx_observed_rep1, dx_observed_rep2 | xmid_latent, dt, common4).

Step 16 reuses the validated Step-15 implementation for all sufficient
statistics, support rules and biological-subject bootstrap calculations.

Cross-alpha pairing
-------------------
All alpha values use identical subject universes, seeds and bootstrap draws.

Forward abundance slopes are additionally evaluated on one common contiguous
Step-9 abundance support valid for T0PLUS, TT and TF at every alpha. This
prevents threshold-dependent loss of low-abundance bins from being mistaken
for a threshold effect on slope.

Native-support forward slopes are retained only as descriptive diagnostics.
Within the paired bootstrap, a common-support forward slope is retained only
when every selected common-support bin is finite in that draw.

Paired differences are therefore reported for:
- forward curves;
- common-support forward abundance slopes;
- fluctuation curves;
- temporal slopes.

Temporal support
----------------
The Step-14 complete-case subject universe and frozen Step-13 abundance core
are retained.

Step 16 then identifies one largest contiguous abundance subset having valid
PRIMARY, NON_FF and TT support for every alpha and every temporal lag.

All threshold-specific temporal profiles and slopes are evaluated on this one
cross-alpha common operational core. A bootstrap draw contributes a temporal
slope only when the entire common core satisfies `min_n` and `min_subjects` at
every lag.

Interpretation
--------------
Step 16 asks whether the Step-15 observation-domain conclusions are robust to
the numerical value used to define T/F.

It does not refit latent abundance, threshold posterior-predictive detectability,
subtract the pseudo technical null, or redefine the primary dynamical estimands.
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
            "Step 16: robustness of current operational-domain sensitivity "
            "analyses to the numerical endpoint observation threshold."
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
        "--step15-script",
        type=Path,
        default=(
            Path(
                __file__
            ).resolve().parent
            / "15-detectability_boundary_sensitivity.py"
        ),
    )
    p.add_argument(
        "--step15-reference-dir",
        type=Path,
        default=None,
        help=(
            "Optional finalized alpha=0.05 Step-15 output directory used only "
            "for an external replication audit."
        ),
    )
    p.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    p.add_argument(
        "--alphas",
        type=float,
        nargs="+",
        default=list(
            DEFAULT_ALPHAS
        ),
    )
    p.add_argument(
        "--reference-alpha",
        type=float,
        default=REFERENCE_ALPHA,
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


def validate_args(
    args: argparse.Namespace,
) -> List[float]:
    alphas = sorted(
        set(
            float(x)
            for x in args.alphas
        ),
        reverse=True,
    )

    if not alphas:
        raise ValueError(
            "At least one alpha is required."
        )

    if any(
        (
            not math.isfinite(
                x
            )
            or x <= 0.0
            or x >= 1.0
        )
        for x in alphas
    ):
        raise ValueError(
            "All alpha values must be finite and in (0,1)."
        )

    if not any(
        math.isclose(
            x,
            float(
                args.reference_alpha
            ),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for x in alphas
    ):
        raise ValueError(
            "--reference-alpha must be included in --alphas."
        )

    if int(
        args.forward_dt
    ) < 1:
        raise ValueError(
            "--forward-dt must be >=1."
        )

    if int(
        args.min_n
    ) < 2:
        raise ValueError(
            "--min-n must be >=2."
        )

    if int(
        args.min_subjects
    ) < 2:
        raise ValueError(
            "--min-subjects must be >=2."
        )

    if int(
        args.n_bootstrap
    ) < 1:
        raise ValueError(
            "--n-bootstrap must be >=1."
        )

    return alphas


# =============================================================================
# Main
# =============================================================================

def main(
    argv: Optional[
        Sequence[str]
    ] = None,
) -> int:
    args = (
        build_parser()
        .parse_args(
            argv
        )
    )

    alphas = validate_args(
        args
    )

    transitions = (
        args.transitions
        .expanduser()
        .resolve(
            strict=True
        )
    )

    step9_dir = (
        args.step9_dir
        .expanduser()
        .resolve(
            strict=True
        )
    )

    step11_dir = (
        args.step11_dir
        .expanduser()
        .resolve(
            strict=True
        )
    )

    step13_dir = (
        args.step13_dir
        .expanduser()
        .resolve(
            strict=True
        )
    )

    step14_dir = (
        args.step14_dir
        .expanduser()
        .resolve(
            strict=True
        )
    )

    step15_script = (
        args.step15_script
        .expanduser()
        .resolve(
            strict=True
        )
    )

    step15_reference_dir = (
        args.step15_reference_dir
        .expanduser()
        .resolve(
            strict=True
        )
        if args.step15_reference_dir
        is not None
        else None
    )

    outdir = ensure_dir(
        args.outdir
        .expanduser()
        .resolve()
    )

    step15 = load_python_module(
        step15_script
    )

    step15_reference_config = None
    step15_reference_signature = None

    if step15_reference_dir is not None:
        step15_reference_config = step15.read_json(
            step15_reference_dir
            / "00_run_config.json"
        )
        step15_reference_signature = step15.read_json(
            step15_reference_dir
            / "00_run_signature.json"
        )

        if not math.isclose(
            float(
                step15_reference_config.get(
                    "reference_operational_alpha",
                    np.nan,
                )
            ),
            float(
                args.reference_alpha
            ),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise ValueError(
                "Step-15 reference directory alpha does not match "
                "--reference-alpha."
            )

        if (
            step15_reference_config.get(
                "analysis_signature"
            )
            != step15_reference_signature.get(
                "analysis_signature"
            )
        ):
            raise ValueError(
                "Step-15 reference run config/signature disagree."
            )

    step9_config_path = (
        step9_dir
        / "00_run_config.json"
    )
    step9_signature_path = (
        step9_dir
        / "00_run_signature.json"
    )
    step11_config_path = (
        step11_dir
        / "00_run_config.json"
    )
    step11_signature_path = (
        step11_dir
        / "00_run_signature.json"
    )
    step13_config_path = (
        step13_dir
        / "00_run_config.json"
    )
    step13_signature_path = (
        step13_dir
        / "00_run_signature.json"
    )
    step14_config_path = (
        step14_dir
        / "00_run_config.json"
    )
    step14_signature_path = (
        step14_dir
        / "00_run_signature.json"
    )

    step9_config = (
        step15
        .read_json(
            step9_config_path
        )
    )
    step9_signature = (
        step15
        .read_json(
            step9_signature_path
        )
    )

    step11_config = (
        step15
        .read_json(
            step11_config_path
        )
    )
    step11_signature = (
        step15
        .read_json(
            step11_signature_path
        )
    )

    step13_config = (
        step15
        .read_json(
            step13_config_path
        )
    )
    step13_signature = (
        step15
        .read_json(
            step13_signature_path
        )
    )

    step14_config = (
        step15
        .read_json(
            step14_config_path
        )
    )
    step14_signature = (
        step15
        .read_json(
            step14_signature_path
        )
    )

    step15.validate_upstream(
        step9_config,
        step9_signature,
        step11_config,
        step11_signature,
        step13_config,
        step13_signature,
        step14_config,
        step14_signature,
    )

    # Step 16 invokes Step-15 functions directly and therefore must reproduce
    # the primary Step-15 threshold-contract validation performed in its main().
    step9_min_n = int(
        step9_config.get(
            "binning",
            {},
        ).get(
            "min_n",
            args.min_n,
        )
    )
    step11_min_n = int(
        step11_config.get(
            "binning",
            {},
        ).get(
            "min_n_per_pooled_cell",
            args.min_n,
        )
    )

    if (
        int(args.min_n) != step9_min_n
        or int(args.min_n) != step11_min_n
    ):
        raise ValueError(
            "Step 16 --min-n must match finalized Step 9 and Step 11. "
            f"requested={args.min_n}, Step9={step9_min_n}, "
            f"Step11={step11_min_n}."
        )

    forward_edges = (
        step15
        .normalize_edges(
            pd.read_csv(
                step9_dir
                / "01_forward_bin_edges.csv"
            ),
            "Step-9 forward bin edges",
        )
    )

    step9_primary = pd.read_csv(
        step9_dir
        / "02_primary_forward_by_bin.csv"
    )

    step11_primary = pd.read_csv(
        step11_dir
        / "02_fluctuation_by_bin_dt.csv"
    )

    fluct_edges = (
        step15
        .normalize_edges(
            step11_primary[
                [
                    "bin",
                    "x_left",
                    "x_right",
                    "x_center",
                ]
            ],
            "Step-11 fluctuation grid",
        )
    )

    (
        frozen_core,
        temporal_subjects,
        temporal_dts,
        frozen_core_bins,
    ) = (
        step15
        .load_step14_contract(
            step14_dir
        )
    )

    print(
        "[1/7] Loading selected Step-5 transition columns..."
    )

    base_df = (
        step15
        .read_transition_table(
            transitions
        )
    )

    memory_gib = (
        float(
            base_df
            .memory_usage(
                index=True,
                deep=True,
            )
            .sum()
        )
        / (
            1024.0
            ** 3
        )
    )

    print(
        f"[INFO] Transition working table: "
        f"{len(base_df):,} rows; "
        f"{memory_gib:.2f} GiB"
    )

    all_subjects = sorted(
        int(x)
        for x in base_df[
            "subject"
        ].dropna().unique()
    )

    print(
        "[2/7] Calibrating Step-9 forward bin membership..."
    )

    (
        forward_assignment_edges,
        forward_calibration,
    ) = (
        step15
        .calibrate_forward_assignment_edges(
            base_df,
            forward_edges,
            step9_primary,
            int(
                args.forward_dt
            ),
        )
    )

    print(
        "[3/7] Cross-alpha operational composition, forward and fluctuation analyses..."
    )

    class_tables = []

    forward_curves_by_alpha: Dict[
        float,
        pd.DataFrame,
    ] = {}
    forward_slopes_by_alpha: Dict[
        float,
        pd.DataFrame,
    ] = {}
    forward_boot_curves_by_alpha: Dict[
        float,
        Mapping[
            str,
            np.ndarray,
        ],
    ] = {}
    forward_boot_slopes_by_alpha: Dict[
        float,
        Mapping[
            str,
            np.ndarray,
        ],
    ] = {}

    fluct_points_by_alpha: Dict[
        float,
        pd.DataFrame,
    ] = {}
    fluct_boot_by_alpha: Dict[
        float,
        Mapping[
            str,
            Mapping[
                str,
                np.ndarray,
            ],
        ],
    ] = {}
    fluct_arrays_by_alpha: Dict[
        float,
        Mapping,
    ] = {}
    support_by_alpha: Dict[
        float,
        pd.DataFrame,
    ] = {}

    for alpha in alphas:
        print(
            f"  [alpha={alpha:g}]"
        )

        composition, _ = (
            step15
            .build_operational_composition(
                base_df,
                float(
                    alpha
                ),
                fluct_edges,
            )
        )

        composition[
            "alpha"
        ] = float(
            alpha
        )
        composition[
            "alpha_tag"
        ] = alpha_tag(
            alpha
        )

        class_tables.append(
            composition
        )

        forward_suff = (
            step15
            .build_forward_sufficient(
                base_df,
                float(
                    alpha
                ),
                int(
                    args.forward_dt
                ),
                forward_assignment_edges,
            )
        )

        forward_counts, forward_sums = (
            step15
            .build_forward_arrays(
                forward_suff,
                all_subjects,
                len(
                    forward_edges
                ),
            )
        )

        forward_curves = (
            step15
            .forward_point_curves(
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
            forward_slope_summary,
            _,
            forward_boot_curves,
            forward_boot_slope_arrays,
        ) = (
            step15
            .bootstrap_forward(
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
        )

        forward_curves = (
            forward_curves.merge(
                forward_boot_ci,
                on=[
                    "mode",
                    "bin",
                ],
                how="left",
                validate="one_to_one",
            )
        )

        forward_curves_by_alpha[
            float(
                alpha
            )
        ] = forward_curves

        forward_slopes_by_alpha[
            float(
                alpha
            )
        ] = forward_slope_summary

        forward_boot_curves_by_alpha[
            float(
                alpha
            )
        ] = forward_boot_curves

        forward_boot_slopes_by_alpha[
            float(
                alpha
            )
        ] = forward_boot_slope_arrays

        fluct_suff = (
            step15
            .build_fluctuation_sufficient(
                base_df,
                float(
                    alpha
                ),
                fluct_edges,
            )
        )

        fluct_arrays = (
            step15
            .build_fluctuation_arrays(
                fluct_suff,
                all_subjects,
                temporal_dts,
                len(
                    fluct_edges
                ),
            )
        )

        fluct_point = (
            step15
            .fluctuation_point_table(
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
            _,
        ) = (
            step15
            .bootstrap_fluctuation(
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
        )

        fluct_points_by_alpha[
            float(
                alpha
            )
        ] = fluct_point

        fluct_boot_by_alpha[
            float(
                alpha
            )
        ] = fluct_boot

        fluct_arrays_by_alpha[
            float(
                alpha
            )
        ] = fluct_arrays

        _, support = (
            step15
            .temporal_support_common_core(
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

        support_by_alpha[
            float(
                alpha
            )
        ] = support

    class_by_alpha_dt = pd.concat(
        class_tables,
        ignore_index=True,
    )

    pooled_composition = pooled_class_composition(
        class_by_alpha_dt
    )

    tt_retention = build_tt_retention(
        step15,
        base_df,
        alphas,
        float(
            args.reference_alpha
        ),
    )

    print(
        "[4/7] Building paired cross-alpha forward and fluctuation contrasts..."
    )

    forward_cross_alpha = combine_forward_curves(
        forward_curves_by_alpha
    )

    # Native-support slope tables are retained only as descriptive diagnostics.
    # They are not valid for cross-alpha slope robustness because stringent
    # alpha values lose low-abundance bins.
    forward_slope_native_support = combine_forward_slopes(
        forward_slopes_by_alpha
    )

    (
        forward_common_bins,
        forward_common_support,
    ) = resolve_forward_cross_alpha_common_support(
        forward_curves_by_alpha,
        forward_edges,
    )

    if not forward_common_bins:
        raise ValueError(
            "No cross-alpha forward common support exists."
        )

    print(
        "[FORWARD COMMON SUPPORT]",
        forward_common_bins,
    )

    (
        forward_common_slopes_by_alpha,
        forward_common_boot_slopes_by_alpha,
    ) = forward_slopes_on_common_support(
        step15,
        forward_curves_by_alpha,
        forward_boot_curves_by_alpha,
        forward_common_bins,
    )

    forward_slope_cross_alpha = pd.concat(
        [
            frame
            for frame
            in forward_common_slopes_by_alpha.values()
        ],
        ignore_index=True,
    ).sort_values(
        [
            "mode",
            "alpha",
        ],
        ascending=[
            True,
            False,
        ],
    ).reset_index(
        drop=True
    )

    forward_curve_diff = (
        forward_curve_difference_vs_reference(
            forward_curves_by_alpha,
            forward_boot_curves_by_alpha,
            float(
                args.reference_alpha
            ),
        )
    )

    forward_slope_diff = (
        forward_slope_difference_vs_reference(
            forward_common_slopes_by_alpha,
            forward_common_boot_slopes_by_alpha,
            float(
                args.reference_alpha
            ),
        )
    )

    fluct_cross_alpha = (
        combine_fluctuation_points(
            fluct_points_by_alpha
        )
    )

    fluct_diff = (
        fluctuation_difference_vs_reference(
            fluct_points_by_alpha,
            fluct_boot_by_alpha,
            temporal_dts,
            fluct_edges,
            float(
                args.reference_alpha
            ),
        )
    )

    print(
        "[5/7] Resolving one cross-alpha temporal common core..."
    )

    (
        cross_alpha_core,
        temporal_support,
    ) = (
        build_cross_alpha_temporal_core(
            support_by_alpha,
            frozen_core,
        )
    )

    if not cross_alpha_core:
        raise ValueError(
            "No cross-alpha common operational temporal core exists."
        )

    print(
        "[CORE]",
        cross_alpha_core,
    )

    temporal_profiles_by_alpha = {}
    temporal_slopes_by_alpha = {}
    temporal_boot_slopes_by_alpha = {}

    for alpha in alphas:
        (
            profiles,
            slopes,
            boot_slopes,
        ) = (
            temporal_bootstrap_for_alpha(
                step15,
                fluct_arrays_by_alpha[
                    float(
                        alpha
                    )
                ],
                all_subjects,
                temporal_subjects,
                temporal_dts,
                cross_alpha_core,
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

        temporal_profiles_by_alpha[
            float(
                alpha
            )
        ] = profiles

        temporal_slopes_by_alpha[
            float(
                alpha
            )
        ] = slopes

        temporal_boot_slopes_by_alpha[
            float(
                alpha
            )
        ] = boot_slopes

    temporal_profiles = (
        combine_temporal_profiles(
            temporal_profiles_by_alpha
        )
    )

    temporal_slopes = (
        combine_temporal_slopes(
            temporal_slopes_by_alpha
        )
    )

    temporal_slope_diff = (
        temporal_slope_difference_vs_reference(
            temporal_slopes_by_alpha,
            temporal_boot_slopes_by_alpha,
            float(
                args.reference_alpha
            ),
        )
    )

    print(
        "[6/7] Audits and compact summaries..."
    )

    internal_audit = (
        build_reference_internal_audit(
            forward_curves_by_alpha,
            fluct_points_by_alpha,
            float(
                args.reference_alpha
            ),
        )
    )

    external_audit = (
        audit_against_step15_reference_outputs(
            step15_reference_dir,
            float(
                args.reference_alpha
            ),
            forward_curves_by_alpha[
                float(
                    args.reference_alpha
                )
            ],
            fluct_points_by_alpha[
                float(
                    args.reference_alpha
                )
            ],
            temporal_slopes_by_alpha[
                float(
                    args.reference_alpha
                )
            ],
        )
    )

    support = support_summary(
        class_by_alpha_dt,
        support_by_alpha,
        cross_alpha_core,
    )

    analysis_summary = (
        build_analysis_summary(
            pooled_composition,
            forward_slope_cross_alpha,
            temporal_slopes,
            forward_slope_diff,
            temporal_slope_diff,
            cross_alpha_core,
        )
    )

    print(
        "[7/7] Writing outputs..."
    )

    class_by_alpha_dt.to_csv(
        outdir
        / "01_class_composition_by_alpha_dt.csv",
        index=False,
    )

    pooled_composition.to_csv(
        outdir
        / "02_class_composition_pooled_by_alpha.csv",
        index=False,
    )

    tt_retention.to_csv(
        outdir
        / "03_tt_retention_by_alpha.csv",
        index=False,
    )

    forward_cross_alpha.to_csv(
        outdir
        / "04_forward_domain_cross_alpha_by_bin.csv",
        index=False,
    )

    forward_slope_cross_alpha.to_csv(
        outdir
        / "05_forward_slope_cross_alpha.csv",
        index=False,
    )

    forward_slope_native_support.to_csv(
        outdir
        / "05a_forward_slope_native_support_diagnostic.csv",
        index=False,
    )

    forward_common_support.to_csv(
        outdir
        / "05b_forward_cross_alpha_common_support.csv",
        index=False,
    )

    forward_curve_diff.to_csv(
        outdir
        / "06_forward_curve_difference_vs_reference_alpha.csv",
        index=False,
    )

    forward_slope_diff.to_csv(
        outdir
        / "07_forward_slope_difference_vs_reference_alpha.csv",
        index=False,
    )

    fluct_cross_alpha.to_csv(
        outdir
        / "08_fluctuation_domain_cross_alpha_by_bin_dt.csv",
        index=False,
    )

    fluct_diff.to_csv(
        outdir
        / "09_fluctuation_difference_vs_reference_alpha.csv",
        index=False,
    )

    temporal_support.to_csv(
        outdir
        / "10_temporal_support_by_alpha_bin.csv",
        index=False,
    )

    frozen_core_out = frozen_core.copy()
    frozen_core_out[
        "selected_cross_alpha_common_core"
    ] = frozen_core_out[
        "bin"
    ].isin(
        cross_alpha_core
    )

    frozen_core_out.to_csv(
        outdir
        / "11_cross_alpha_temporal_common_core.csv",
        index=False,
    )

    temporal_profiles.to_csv(
        outdir
        / "12_temporal_profiles_cross_alpha.csv",
        index=False,
    )

    temporal_slopes.to_csv(
        outdir
        / "13_temporal_slopes_cross_alpha.csv",
        index=False,
    )

    temporal_slope_diff.to_csv(
        outdir
        / "14_temporal_slope_difference_vs_reference_alpha.csv",
        index=False,
    )

    internal_audit.to_csv(
        outdir
        / "15_reference_alpha_internal_audit.csv",
        index=False,
    )

    if not external_audit.empty:
        external_audit.to_csv(
            outdir
            / "16_reference_step15_output_audit.csv",
            index=False,
        )

    support.to_csv(
        outdir
        / "17_support_summary_by_alpha.csv",
        index=False,
    )

    forward_calibration.to_csv(
        outdir
        / "18_forward_bin_assignment_calibration.csv",
        index=False,
    )

    analysis_summary.to_csv(
        outdir
        / "00_analysis_summary.csv",
        index=False,
    )

    signature_payload = {
        "script_version": SCRIPT_VERSION,
        "transitions": path_identity(
            transitions
        ),
        "step15_script": path_identity(
            step15_script
        ),
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
        "step15_reference_analysis_signature": (
            step15_reference_config.get(
                "analysis_signature"
            )
            if step15_reference_config
            is not None
            else None
        ),
        "parameters": {
            "alphas": [
                float(x)
                for x in alphas
            ],
            "reference_alpha": float(
                args.reference_alpha
            ),
            "forward_dt": int(
                args.forward_dt
            ),
            "min_n": int(
                args.min_n
            ),
            "min_subjects": int(
                args.min_subjects
            ),
            "n_bootstrap": int(
                args.n_bootstrap
            ),
            "seed": int(
                args.seed
            ),
        },
        "forward_common_support_bins": [
            int(x)
            for x in forward_common_bins
        ],
        "cross_alpha_temporal_common_core_bins": [
            int(x)
            for x in cross_alpha_core
        ],
        "forward_bootstrap_requires_complete_common_support": True,
        "temporal_bootstrap_requires_complete_common_core_all_dt": True,
        "same_subject_draws_across_alpha": True,
    }

    analysis_signature = stable_sha256(
        signature_payload
    )

    (
        outdir
        / "00_run_signature.json"
    ).write_text(
        json.dumps(
            {
                "analysis_signature":
                    analysis_signature,
                "signature_payload":
                    signature_payload,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    input_manifest_rows = [
        {
            "input":
                "transitions",
            **path_identity(
                transitions
            ),
        },
        {
            "input":
                "step15_script",
            **path_identity(
                step15_script
            ),
        },
        {
            "input":
                "step9_run_config",
            **path_identity(
                step9_config_path
            ),
        },
        {
            "input":
                "step9_run_signature",
            **path_identity(
                step9_signature_path
            ),
        },
        {
            "input":
                "step11_run_config",
            **path_identity(
                step11_config_path
            ),
        },
        {
            "input":
                "step11_run_signature",
            **path_identity(
                step11_signature_path
            ),
        },
        {
            "input":
                "step13_run_config",
            **path_identity(
                step13_config_path
            ),
        },
        {
            "input":
                "step13_run_signature",
            **path_identity(
                step13_signature_path
            ),
        },
        {
            "input":
                "step14_run_config",
            **path_identity(
                step14_config_path
            ),
        },
        {
            "input":
                "step14_run_signature",
            **path_identity(
                step14_signature_path
            ),
        },
    ]

    if step15_reference_dir is not None:
        input_manifest_rows.extend(
            [
                {
                    "input":
                        "step15_reference_run_config",
                    **path_identity(
                        step15_reference_dir
                        / "00_run_config.json"
                    ),
                },
                {
                    "input":
                        "step15_reference_run_signature",
                    **path_identity(
                        step15_reference_dir
                        / "00_run_signature.json"
                    ),
                },
            ]
        )

    pd.DataFrame(
        input_manifest_rows
    ).to_csv(
        outdir
        / "00_input_manifest.csv",
        index=False,
    )

    run_config = {
        "script":
            Path(
                __file__
            ).name,
        "script_version":
            SCRIPT_VERSION,
        "analysis_signature":
            analysis_signature,
        "scientific_scope":
            "operational_observation_threshold_robustness",
        "classification_rule":
            "T(alpha)=1[p_endpoint<alpha]",
        "interpretation":
            (
                "operational observation state; "
                "not biological presence/absence"
            ),
        "alphas":
            [
                float(
                    x
                )
                for x in alphas
            ],
        "reference_alpha":
            float(
                args.reference_alpha
            ),
        "latent_states_refitted":
            False,
        "posterior_predictive_detectability_thresholded":
            False,
        "primary_estimands_changed":
            False,
        "step12_used_computationally":
            False,
        "transitions":
            str(
                transitions
            ),
        "step9_dir":
            str(
                step9_dir
            ),
        "step11_dir":
            str(
                step11_dir
            ),
        "step13_dir":
            str(
                step13_dir
            ),
        "step14_dir":
            str(
                step14_dir
            ),
        "step15_script":
            str(
                step15_script
            ),
        "step15_script_version":
            getattr(
                step15,
                "SCRIPT_VERSION",
                None,
            ),
        "step15_reference_dir":
            (
                str(
                    step15_reference_dir
                )
                if step15_reference_dir
                is not None
                else None
            ),
        "forward": {
            "dt":
                int(
                    args.forward_dt
                ),
            "modes":
                [
                    "PRIMARY",
                    "T0PLUS",
                    "TT",
                    "TF",
                ],
            "primary_alpha_invariant":
                True,
            "bin_coordinates":
                "final Step-9 coordinates",
            "bin_membership":
                "Step-15 calibrated assignment boundaries",
            "cross_alpha_slope_support":
                "largest contiguous Step-9 bin set valid for T0PLUS, TT and TF at every alpha",
            "cross_alpha_common_support_bins":
                [
                    int(
                        x
                    )
                    for x
                    in forward_common_bins
                ],
            "native_support_slopes_role":
                "descriptive diagnostic only; not used for threshold-comparison inference",
            "bootstrap_common_support_policy":
                "slope retained only when every common-support bin is finite in that draw",
        },
        "fluctuations": {
            "base_support":
                "common4",
            "conditioning":
                "xmid_latent",
            "metrics":
                list(
                    METRICS
                ),
            "modes":
                [
                    "PRIMARY",
                    "NON_FF",
                    "TT",
                ],
            "primary_alpha_invariant":
                True,
        },
        "temporal": {
            "subject_universe":
                [
                    int(
                        x
                    )
                    for x
                    in temporal_subjects
                ],
            "dt_values":
                [
                    int(
                        x
                    )
                    for x
                    in temporal_dts
                ],
            "frozen_step13_14_core_bins":
                [
                    int(
                        x
                    )
                    for x
                    in frozen_core_bins
                ],
            "cross_alpha_common_core_bins":
                [
                    int(
                        x
                    )
                    for x
                    in cross_alpha_core
                ],
            "estimator":
                "common_core_pooled_covariance_equal_bin",
            "role":
                (
                    "cross-alpha threshold sensitivity; "
                    "does not replace Step-13 primary estimator"
                ),
            "bootstrap_common_core_policy":
                (
                    "slope retained only when the full cross-alpha common "
                    "core satisfies min_n/min_subjects at every lag in that draw"
                ),
        },
        "support_thresholds": {
            "min_n":
                int(
                    args.min_n
                ),
            "min_subjects":
                int(
                    args.min_subjects
                ),
        },
        "bootstrap": {
            "unit":
                "biological subject",
            "n_bootstrap":
                int(
                    args.n_bootstrap
                ),
            "seed":
                int(
                    args.seed
                ),
            "same_draws_across_alpha":
                True,
            "same_draws_across_operational_domains":
                True,
            "forward_complete_common_support_required_per_draw":
                True,
            "temporal_complete_common_core_all_dt_required_per_draw":
                True,
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
        "\n[DONE] Step 16 operational observation-threshold robustness"
    )
    print(
        "Alphas                  :",
        alphas,
    )
    print(
        "Reference alpha         :",
        float(
            args.reference_alpha
        ),
    )
    print(
        "Forward common support  :",
        forward_common_bins,
    )
    print(
        "Frozen Step-13/14 core :",
        frozen_core_bins,
    )
    print(
        "Cross-alpha common core :",
        cross_alpha_core,
    )
    print(
        "Temporal subjects       :",
        temporal_subjects,
    )

    if not external_audit.empty:
        identity_rows = external_audit[
            external_audit[
                "component"
            ].isin(
                [
                    "forward_domain_by_bin",
                    "fluctuation_domain_by_bin_dt",
                ]
            )
        ]

        if not identity_rows.empty:
            print(
                "Reference Step-15 audit max |diff|:",
                float(
                    identity_rows[
                        "max_abs_difference"
                    ].max()
                ),
            )

    print(
        "\nOutputs:"
    )

    for filename in [
        "00_analysis_summary.csv",
        "00_run_signature.json",
        "00_input_manifest.csv",
        "01_class_composition_by_alpha_dt.csv",
        "02_class_composition_pooled_by_alpha.csv",
        "03_tt_retention_by_alpha.csv",
        "04_forward_domain_cross_alpha_by_bin.csv",
        "05_forward_slope_cross_alpha.csv",
        "05a_forward_slope_native_support_diagnostic.csv",
        "05b_forward_cross_alpha_common_support.csv",
        "06_forward_curve_difference_vs_reference_alpha.csv",
        "07_forward_slope_difference_vs_reference_alpha.csv",
        "08_fluctuation_domain_cross_alpha_by_bin_dt.csv",
        "09_fluctuation_difference_vs_reference_alpha.csv",
        "10_temporal_support_by_alpha_bin.csv",
        "11_cross_alpha_temporal_common_core.csv",
        "12_temporal_profiles_cross_alpha.csv",
        "13_temporal_slopes_cross_alpha.csv",
        "14_temporal_slope_difference_vs_reference_alpha.csv",
        "15_reference_alpha_internal_audit.csv",
        "16_reference_step15_output_audit.csv",
        "17_support_summary_by_alpha.csv",
        "18_forward_bin_assignment_calibration.csv",
    ]:
        path = (
            outdir
            / filename
        )

        if path.exists():
            print(
                " -",
                path,
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
