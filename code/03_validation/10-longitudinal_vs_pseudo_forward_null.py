#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
10-longitudinal_vs_pseudo_forward_null.py
=========================================

Canonical ClonoDynamics Step 10 under the final replicate-resolved framework.

PURPOSE
-------
Compare the genuine longitudinal observed replicate-decoupled forward profile
from Step 9 with the empirical pseudo-longitudinal technical-null ensemble from
Step 7.

The SAME forward estimand is used in both datasets:

    AB:
        x_A(t0) -> Delta x_B

    BA:
        x_B(t0) -> Delta x_A

    cross_combined:
        equal AB/BA fold weight after binning.

No TT/TF/FT/FF restriction is introduced by Step 10.

TWO COMPLEMENTARY GEOMETRIES
----------------------------

A. ABSOLUTE observed log-frequency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Each dataset retains its native Step-7 / Step-9 abundance-bin geometry.

The already-estimated genuine and pseudo curves are compared only on their
shared ABSOLUTE observed-log-frequency support. No extrapolation is allowed.

This branch answers:

    At the same observable abundance scale represented in both datasets,
    does the genuine longitudinal forward profile differ from the technical
    null?

Because the single-donor pseudo reference can have a narrower dynamic range
than the 10-subject longitudinal cohort, this comparison may cover only part of
the genuine restoring-like profile.


B. WITHIN-UNIT abundance percentile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A second branch compares RELATIVE abundance geometry across the complete
supported repertoire range.

Percentiles are assigned at the CLONOTYPE-TRANSITION level from the initial
observed log-frequency:

    longitudinal:
        within biological subject x fold

    pseudo:
        within randomization configuration x fold

AB and BA are percentile-ranked separately because they use different
conditioning measurements.

Percentiles are mid-rank empirical-CDF coordinates:

    p = (average_rank - 0.5) / n

and are binned on a common fixed [0,1] grid.

The percentile branch does NOT imply that the same percentile corresponds to
the same absolute frequency in the two datasets. It compares the SHAPE of the
abundance-conditioned forward profile after removing between-dataset dynamic-
range differences.

PRIMARY INPUTS
--------------
--longitudinal-dir
    Current Step-9 output directory, expected to contain:

        00_run_config.json
        02_primary_forward_by_bin.csv
        05_subject_forward_by_bin.csv
        07_AB_BA_concordance.csv

    The Step-9 run config provides the canonical longitudinal Step-5 transition
    table used for the percentile branch.

--pseudo-step7-dir
    Current compact-native Step-7 output directory, expected to contain:

        00_run_config.json
        00_run_signature.json
        00_configuration_manifest.csv
        03_configuration_curves_long.parquet

    The Step-7 manifest directly identifies each configuration's
    `step7_cache/step7_dt1.parquet`. No configuration-specific full Step-5
    transition table is required.

RANDOMIZATION UNIT
------------------
Each pseudo `configuration_id` is a randomization realization conditional on
ONE biological source sample. It is never treated as a biological subject.

Pseudo 2.5th-97.5th percentile intervals are randomization intervals, not
biological confidence intervals.

LONGITUDINAL UNCERTAINTY
------------------------
The absolute branch uses the subject-cluster bootstrap confidence intervals
already produced by Step 9.

The percentile branch reconstructs subject x fold x percentile-bin sufficient
statistics and performs the same biological subject-cluster bootstrap:

    default n = inherited from Step 9 (normally 2000)
    default seed = inherited from Step 9 (normally 123)

The same subject draw is used for AB and BA before equal-fold combination.

ABSOLUTE-BRANCH EMPIRICAL NULL
------------------------------
The full Step-7 configuration-level `cross_combined` table first defines
stable pseudo native bins using the fraction of all 2,000 configurations for
which the combined bin satisfies `meets_min_n=True`.

Only bins meeting:

    --absolute-min-pseudo-bin-coverage

are eligible for the absolute comparison. The production default is 0.95.

The shared comparison coordinates are those stable pseudo native-bin centers
that lie inside the supported longitudinal Step-9 range.

For each pseudo configuration, ONLY DIRECTLY VALID `meets_min_n=True`
native-bin values are used. Missing pseudo bins are never interpolated,
imputed, or connected across gaps.

In production, a pseudo configuration enters the curve-level empirical null
only when it has a directly valid value at EVERY stable shared bin. Therefore
the genuine longitudinal metric and every eligible pseudo metric are evaluated
on the same abundance coordinates.

`--absolute-min-grid-coverage < 1.0` is permitted only together with
`--max-pseudo-configs` for structural smoke testing; production inference
requires complete shared-grid coverage.

Curve-level metrics:
    mean_dx
    slope
    centered_mean_abs_drift
    centered_rms_drift
    mean_abs_drift

Prespecified empirical tests use +1 finite-sample correction:

    slope:
        lower tail

    centered_mean_abs_drift:
        upper tail

    centered_rms_drift:
        upper tail

    mean_abs_drift:
        upper tail

    mean_dx:
        two-sided about the pseudo median


PERCENTILE-BRANCH EMPIRICAL NULL
--------------------------------
For each pseudo configuration, AB and BA transition-level percentile curves are
estimated separately and then combined with equal fold weight.

A percentile bin is retained for a fold when:

    n >= --percentile-min-n

A configuration-level combined curve requires at least:

    --percentile-min-valid-bins

valid bins.

Configuration curves are interpolated only within their valid percentile
support onto the fixed percentile-bin centers. Configurations satisfying the
requested grid-coverage rule define the empirical percentile-space null.

The same curve-level metrics and empirical tests are reported as in the
absolute branch.

DIRECTIONAL CROSSOVER
---------------------
For the genuine percentile-space curve the script reports the compact linear
descriptor:

    mean_dx = intercept + slope * percentile

and:

    crossover_percentile = -intercept / slope

when finite.

A subject-bootstrap interval for the crossover percentile is also reported.

This is a compact descriptor of the binned profile, not a claim that the true
forward law is exactly linear.

OUTPUT
------
<outdir>/
    00_run_config.json
    00_input_manifest.csv
    00_analysis_summary.csv
    00_pipeline_manifest.csv
    README_outputs.md

    01_absolute_abundance/
        01_native_curves.csv
        02_shared_grid.csv
        03_pseudo_null_metrics.csv
        04_empirical_tests.csv
        05_support_summary.csv
        06_comparison_summary.csv
        07_pseudo_native_bin_support.csv

    02_percentile_geometry/
        01_longitudinal_percentile_curve.csv
        02_pseudo_configuration_percentile_curves.parquet
        03_pseudo_percentile_envelope.csv
        04_pseudo_null_metrics.csv
        05_empirical_tests.csv
        06_pointwise_comparison.csv
        07_longitudinal_bootstrap_linear_descriptors.csv
        08_support_summary.csv
        09_comparison_summary.csv
        00_percentile_cache_signature.json

TESTING
-------
Use:

    --max-pseudo-configs 10

to smoke-test the analysis on the first 10 configurations. The full Step-7
ensemble is still checked before subsetting, and the frozen full-ensemble
Step-7 envelope continues to define absolute-abundance support. In test mode,
one eligible pseudo configuration is sufficient for structural execution, but
the inferential minimum remains unchanged. Therefore empirical P-values are
reported as NA whenever the requested inferential null size is not reached.
Test-mode P-values must not be reported.

PRODUCTION RUN
--------------
Production absolute inference uses:

    --absolute-min-pseudo-bin-coverage 0.95
    --absolute-min-grid-coverage 1.0

so missing pseudo bins are never imputed and every eligible configuration is
summarized on the same directly observed shared native-bin coordinates.

python3 10-longitudinal_vs_pseudo_forward_null.py \
    --longitudinal-dir \
    ./results_4/9-observed_replicate_decoupled_forward_drift \
    --pseudo-step7-dir \
    ./results_4/7-pseudo_forward_technical_null_validation \
    --outdir \
    ./results_4/10-longitudinal_vs_pseudo_forward_null \
    --expected-pseudo-configs 2000 \
    --percentile-bins 20 \
    --percentile-min-n 50 \
    --percentile-min-valid-bins 15 \
    --seed 123
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import polars as pl
except Exception as exc:
    raise SystemExit(
        "ERROR: polars is required. Activate the ClonoDynamics environment."
    ) from exc


SCRIPT_VERSION = "v5-no-pseudo-gap-interpolation-2026-09-17"
PERCENTILE_CACHE_VERSION = "v4-compact-step7-step9v3-2026-09-17"

PRIMARY_REPRESENTATION = "cross_combined"
PERCENTILE_METRICS = (
    "mean_dx",
    "slope",
    "centered_mean_abs_drift",
    "centered_rms_drift",
    "mean_abs_drift",
)

PERCENTILE_FOLDS = {
    "AB": {
        "x": "x_observed_rep1_t0",
        "dx": "dx_observed_rep2",
        "support": "forward_ab_eligible",
    },
    "BA": {
        "x": "x_observed_rep2_t0",
        "dx": "dx_observed_rep1",
        "support": "forward_ba_eligible",
    },
}

REQUIRED_TRANSITION_COLUMNS = [
    "subject",
    "dt",
    "x_observed_rep1_t0",
    "x_observed_rep2_t0",
    "dx_observed_rep1",
    "dx_observed_rep2",
    "forward_ab_eligible",
    "forward_ba_eligible",
]


# =============================================================================
# Generic helpers
# =============================================================================

def ensure_dir(path: Path) -> Path:
    path = Path(path)
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


def stable_sha256(payload: Mapping[str, object]) -> str:
    text = json.dumps(
        json_safe(dict(payload)),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_json(path: Path, required: bool = True) -> Dict[str, object]:
    if not path.is_file():
        if required:
            raise FileNotFoundError(path)
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_csv(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    d = pd.read_csv(path)
    if d.empty:
        raise ValueError(f"{label} is empty: {path}")
    return d



def read_table(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    if path.suffix.lower() in {".parquet", ".pq"}:
        d = pd.read_parquet(path)
    elif path.suffix.lower() == ".csv":
        d = pd.read_csv(path)
    else:
        raise ValueError(
            f"{label}: unsupported table format {path.suffix}"
        )
    if d.empty:
        raise ValueError(f"{label} is empty: {path}")
    return d


def bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    if pd.api.types.is_numeric_dtype(values):
        return pd.to_numeric(
            values, errors="coerce"
        ).fillna(0).ne(0)
    return (
        values.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "t", "1", "yes", "y"])
    )


def path_identity(path: Path) -> Dict[str, object]:
    p = path.expanduser().resolve(strict=True)
    st = p.stat()
    return {
        "path": str(p),
        "size_bytes": int(st.st_size),
        "mtime_ns": int(st.st_mtime_ns),
    }


def require_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
    label: str,
) -> None:
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise ValueError(
            f"{label}: missing required columns {missing}; "
            f"available={list(frame.columns)}"
        )


def schema_names(path: Path) -> List[str]:
    return list(pl.scan_parquet(path).collect_schema().names())


def collect_frame(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except Exception:
        return lf.collect()


def parse_bool_expr(column: str) -> pl.Expr:
    s = pl.col(column)
    return (
        pl.when(s.cast(pl.Boolean, strict=False).is_not_null())
        .then(s.cast(pl.Boolean, strict=False).fill_null(False))
        .otherwise(
            s.cast(pl.String, strict=False)
            .str.strip_chars()
            .str.to_lowercase()
            .is_in(["true", "t", "1", "yes", "y"])
            .fill_null(False)
        )
    )


def finite_numeric(values) -> np.ndarray:
    x = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    return x[np.isfinite(x)]


def safe_linear_descriptor(
    x: np.ndarray,
    y: np.ndarray,
) -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)

    if int(ok.sum()) < 2:
        return {
            "slope": np.nan,
            "intercept": np.nan,
            "x_zero": np.nan,
        }

    xv = x[ok]
    yv = y[ok]
    if np.ptp(xv) <= 0:
        return {
            "slope": np.nan,
            "intercept": np.nan,
            "x_zero": np.nan,
        }

    slope, intercept = np.polyfit(xv, yv, 1)
    slope = float(slope)
    intercept = float(intercept)
    x_zero = (
        float(-intercept / slope)
        if np.isfinite(slope) and abs(slope) > 1e-12
        else np.nan
    )
    return {
        "slope": slope,
        "intercept": intercept,
        "x_zero": x_zero,
    }


def curve_metrics(
    x: np.ndarray,
    y: np.ndarray,
) -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]

    if len(y) < 3:
        return {
            "n_points": int(len(y)),
            "mean_dx": np.nan,
            "slope": np.nan,
            "centered_mean_abs_drift": np.nan,
            "centered_rms_drift": np.nan,
            "mean_abs_drift": np.nan,
        }

    mean_y = float(np.mean(y))
    centered = y - mean_y
    linear = safe_linear_descriptor(x, y)

    return {
        "n_points": int(len(y)),
        "mean_dx": mean_y,
        "slope": linear["slope"],
        "centered_mean_abs_drift": float(
            np.mean(np.abs(centered))
        ),
        "centered_rms_drift": float(
            np.sqrt(np.mean(centered * centered))
        ),
        "mean_abs_drift": float(
            np.mean(np.abs(y))
        ),
    }


def summarize_values(values: np.ndarray) -> Dict[str, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "sd": np.nan,
            "q005": np.nan,
            "q025": np.nan,
            "median": np.nan,
            "q975": np.nan,
            "q995": np.nan,
        }

    q = np.quantile(
        x,
        [0.005, 0.025, 0.5, 0.975, 0.995],
    )
    return {
        "n": int(x.size),
        "mean": float(np.mean(x)),
        "sd": (
            float(np.std(x, ddof=1))
            if x.size > 1
            else 0.0
        ),
        "q005": float(q[0]),
        "q025": float(q[1]),
        "median": float(q[2]),
        "q975": float(q[3]),
        "q995": float(q[4]),
    }


def empirical_p_lower(
    null_values: np.ndarray,
    observed: float,
) -> float:
    x = np.asarray(null_values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0 or not np.isfinite(observed):
        return np.nan
    return float(
        (1 + np.sum(x <= float(observed)))
        / (x.size + 1)
    )


def empirical_p_upper(
    null_values: np.ndarray,
    observed: float,
) -> float:
    x = np.asarray(null_values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0 or not np.isfinite(observed):
        return np.nan
    return float(
        (1 + np.sum(x >= float(observed)))
        / (x.size + 1)
    )


def empirical_p_two_sided_about_median(
    null_values: np.ndarray,
    observed: float,
) -> float:
    x = np.asarray(null_values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0 or not np.isfinite(observed):
        return np.nan

    med = float(np.median(x))
    d_obs = abs(float(observed) - med)
    d_null = np.abs(x - med)

    return float(
        (1 + np.sum(d_null >= d_obs))
        / (x.size + 1)
    )


def empirical_tests(
    observed_metrics: Mapping[str, float],
    pseudo_metrics: pd.DataFrame,
    branch: str,
    min_null_configurations: int,
) -> pd.DataFrame:
    specs = [
        (
            "slope",
            "lower",
            "genuine slope more negative than pseudo null",
        ),
        (
            "centered_mean_abs_drift",
            "upper",
            "genuine centered mean-absolute structure exceeds pseudo null",
        ),
        (
            "centered_rms_drift",
            "upper",
            "genuine centered RMS structure exceeds pseudo null",
        ),
        (
            "mean_abs_drift",
            "upper",
            "genuine mean absolute displacement exceeds pseudo null",
        ),
        (
            "mean_dx",
            "two_sided_about_median",
            "genuine mean displacement differs from pseudo-null median",
        ),
    ]

    rows = []
    for metric, tail, interpretation in specs:
        null = pd.to_numeric(
            pseudo_metrics[metric],
            errors="coerce",
        ).to_numpy(float)
        null = null[np.isfinite(null)]
        observed = float(observed_metrics.get(metric, np.nan))
        summary = summarize_values(null)

        if int(summary["n"]) < int(min_null_configurations):
            p = np.nan
            status = "insufficient_null_configurations"
        else:
            if tail == "lower":
                p = empirical_p_lower(null, observed)
            elif tail == "upper":
                p = empirical_p_upper(null, observed)
            else:
                p = empirical_p_two_sided_about_median(
                    null,
                    observed,
                )
            status = "complete"

        rows.append(
            {
                "branch": branch,
                "metric": metric,
                "genuine_value": observed,
                "pseudo_n": int(summary["n"]),
                "pseudo_mean": summary["mean"],
                "pseudo_sd": summary["sd"],
                "pseudo_q005": summary["q005"],
                "pseudo_q025": summary["q025"],
                "pseudo_median": summary["median"],
                "pseudo_q975": summary["q975"],
                "pseudo_q995": summary["q995"],
                "tail": tail,
                "empirical_p": p,
                "status": status,
                "interpretation": interpretation,
            }
        )

    return pd.DataFrame(rows)


def interpolate_inside_support(
    x: np.ndarray,
    y: np.ndarray,
    grid: np.ndarray,
) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    grid = np.asarray(grid, dtype=float)

    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]

    if len(x) < 2:
        return np.full(grid.shape, np.nan)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    ux, inverse = np.unique(x, return_inverse=True)
    if len(ux) != len(x):
        sums = np.zeros(len(ux), dtype=float)
        counts = np.zeros(len(ux), dtype=float)
        for i, j in enumerate(inverse):
            sums[j] += y[i]
            counts[j] += 1
        x = ux
        y = sums / counts

    out = np.full(grid.shape, np.nan)
    inside = (grid >= np.min(x)) & (grid <= np.max(x))
    out[inside] = np.interp(
        grid[inside],
        x,
        y,
    )
    return out


def write_manifest(outdir: Path) -> None:
    rows = []
    for path in sorted(outdir.rglob("*")):
        if path.is_file():
            rows.append(
                {
                    "relative_path": str(
                        path.relative_to(outdir)
                    ),
                    "suffix": path.suffix.lower(),
                    "size_bytes": int(path.stat().st_size),
                }
            )
    pd.DataFrame(rows).to_csv(
        outdir / "00_pipeline_manifest.csv",
        index=False,
    )


# =============================================================================
# Input resolution / validation
# =============================================================================

def resolve_step7_curves_path(pseudo_dir: Path) -> Path:
    parquet = pseudo_dir / "03_configuration_curves_long.parquet"
    csv = pseudo_dir / "03_configuration_curves_long.csv"
    if parquet.is_file():
        return parquet
    if csv.is_file():
        return csv
    raise FileNotFoundError(
        "Step-7 configuration curves not found as Parquet or CSV under "
        f"{pseudo_dir}"
    )


def resolve_inputs(
    longitudinal_dir: Path,
    pseudo_step7_dir: Path,
) -> Dict[str, object]:
    long_dir = longitudinal_dir.expanduser().resolve(strict=True)
    pseudo_dir = pseudo_step7_dir.expanduser().resolve(strict=True)

    step7_curves = resolve_step7_curves_path(pseudo_dir)

    paths = {
        "step9_run_config": long_dir / "00_run_config.json",
        "step9_run_signature": long_dir / "00_run_signature.json",
        "step9_primary": long_dir / "02_primary_forward_by_bin.csv",
        "step7_run_config": pseudo_dir / "00_run_config.json",
        "step7_run_signature": pseudo_dir / "00_run_signature.json",
        "step7_manifest": pseudo_dir / "00_configuration_manifest.csv",
        "step7_curves": step7_curves,
    }

    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Required final Step-9/Step-7 files missing:\n  "
            + "\n  ".join(missing)
        )

    step9_cfg = read_json(paths["step9_run_config"])
    step7_cfg = read_json(paths["step7_run_config"])
    step9_sig = read_json(paths["step9_run_signature"])
    step7_sig = read_json(paths["step7_run_signature"])

    long_transitions = Path(str(step9_cfg.get("transitions", ""))).expanduser()
    if not long_transitions.is_file():
        raise FileNotFoundError(
            "Step-9 run config does not resolve to the longitudinal "
            f"transition table: {long_transitions}"
        )

    step9_primary = step9_cfg.get("primary_estimand", {})
    if (
        step9_primary.get("AB") != "x_observed_rep1_t0 -> dx_observed_rep2"
        or step9_primary.get("BA") != "x_observed_rep2_t0 -> dx_observed_rep1"
    ):
        raise ValueError(
            "Step 9 does not expose the expected observed AB/BA estimand."
        )
    if (
        step9_primary.get("AB_BA_combination_policy")
        != "exact_equal_fold_weight_after_binning"
    ):
        raise ValueError(
            "Step 9 is not the final exact equal-weight AB/BA revision."
        )
    if step9_primary.get("AB_BA_row_count_weighting") is not False:
        raise ValueError(
            "Step 9 must explicitly disable AB/BA row-count weighting."
        )

    step9_support = step9_cfg.get("primary_support", {})
    if step9_support.get("operational_TT_filter") is True:
        raise ValueError(
            "Final Step 10 requires Step 9 without a primary TT restriction."
        )

    step7_primary = step7_cfg.get("primary_forward", {})
    if not isinstance(step7_primary, dict):
        raise ValueError("Invalid Step-7 primary_forward metadata.")

    ab = step7_primary.get("AB", {})
    ba = step7_primary.get("BA", {})
    if (
        not isinstance(ab, dict)
        or ab.get("conditioning") != "x_observed_rep1_t0"
        or ab.get("displacement") != "dx_observed_rep2"
        or ab.get("eligibility") != "forward_ab_eligible"
    ):
        raise ValueError(
            "Step 7 AB contract is not the expected compact-native estimator."
        )
    if (
        not isinstance(ba, dict)
        or ba.get("conditioning") != "x_observed_rep2_t0"
        or ba.get("displacement") != "dx_observed_rep1"
        or ba.get("eligibility") != "forward_ba_eligible"
    ):
        raise ValueError(
            "Step 7 BA contract is not the expected compact-native estimator."
        )
    if step7_primary.get("row_count_weighting") is not False:
        raise ValueError(
            "Step 7 must explicitly disable AB/BA row-count weighting."
        )

    sig_payload = step7_sig.get("signature_payload", {})
    if (
        sig_payload.get("AB_BA_combination_policy")
        != "exact_equal_fold_weight_after_binning"
        or sig_payload.get("AB_BA_row_count_weighting") is not False
    ):
        raise ValueError(
            "Step-7 run signature does not certify the final equal-weight "
            "AB/BA contract."
        )

    params = sig_payload.get("parameters", {})
    step7_dt = params.get("dt") if isinstance(params, dict) else None

    manifest = read_csv(
        paths["step7_manifest"],
        "Step-7 configuration manifest",
    )
    require_columns(
        manifest,
        ["configuration_id", "compact"],
        "Step-7 configuration manifest",
    )
    manifest["configuration_id"] = manifest["configuration_id"].astype(str)
    manifest["compact"] = manifest["compact"].astype(str)

    compact_map: Dict[str, Path] = {}
    for row in manifest.itertuples(index=False):
        cid = str(row.configuration_id)
        path = Path(str(row.compact)).expanduser()
        if not path.is_file():
            raise FileNotFoundError(
                f"{cid}: Step-7 compact cache not found: {path}"
            )
        compact_map[cid] = path.resolve(strict=True)

    return {
        "longitudinal_dir": long_dir,
        "pseudo_step7_dir": pseudo_dir,
        "paths": paths,
        "step9_config": step9_cfg,
        "step7_config": step7_cfg,
        "step9_signature": step9_sig,
        "step7_signature": step7_sig,
        "step7_dt": step7_dt,
        "longitudinal_transitions": long_transitions.resolve(strict=True),
        "pseudo_compact_map": compact_map,
        "step7_manifest": manifest,
    }

def validate_transition_schema(
    path: Path,
    label: str,
) -> None:
    available = set(schema_names(path))
    missing = [
        c for c in REQUIRED_TRANSITION_COLUMNS
        if c not in available
    ]
    if missing:
        raise ValueError(
            f"{label}: transition table is missing {missing}"
        )


# =============================================================================
# Absolute-abundance branch
# =============================================================================

def load_real_absolute_curve(
    path: Path,
) -> pd.DataFrame:
    d = read_csv(path, "Step-9 primary forward curve")
    require_columns(
        d,
        [
            "representation",
            "x_center",
            "mean_dx",
            "mean_dx_ci025",
            "mean_dx_ci975",
            "meets_min_n",
        ],
        "Step-9 primary forward curve",
    )

    d = d[
        d["representation"].astype(str).eq(
            PRIMARY_REPRESENTATION
        )
    ].copy()
    d = d[
        d["meets_min_n"].astype(str)
        .str.lower()
        .isin(["true", "1"])
    ].copy()

    for c in [
        "x_center",
        "mean_dx",
        "mean_dx_ci025",
        "mean_dx_ci975",
    ]:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    d = d[
        np.isfinite(d["x_center"])
        & np.isfinite(d["mean_dx"])
    ].sort_values("x_center")

    if len(d) < 3:
        raise ValueError(
            "Too few valid Step-9 cross_combined bins."
        )
    return d.reset_index(drop=True)


def load_pseudo_absolute_curves(
    path: Path,
) -> pd.DataFrame:
    """
    Load final Step-7 configuration-level cross_combined curves.

    Only `meets_min_n=True` points are considered valid for the empirical
    technical null. The full table is retained so stable-bin coverage can be
    estimated across all configurations before any smoke-test subsetting.
    """
    d = read_table(path, "Step-7 configuration curves")
    require_columns(
        d,
        [
            "configuration_id",
            "section",
            "representation",
            "bin",
            "x_center",
            "mean_dx",
            "meets_min_n",
        ],
        "Step-7 configuration curves",
    )

    d = d[
        d["section"].astype(str).eq("primary_forward")
        & d["representation"].astype(str).eq(PRIMARY_REPRESENTATION)
    ].copy()

    d["configuration_id"] = d["configuration_id"].astype(str)
    d["bin"] = pd.to_numeric(d["bin"], errors="coerce")
    d["x_center"] = pd.to_numeric(d["x_center"], errors="coerce")
    d["mean_dx"] = pd.to_numeric(d["mean_dx"], errors="coerce")
    d["meets_min_n_bool"] = bool_series(d["meets_min_n"])

    d = d[
        np.isfinite(d["bin"])
        & np.isfinite(d["x_center"])
        & np.isfinite(d["mean_dx"])
    ].copy()
    if d.empty:
        raise ValueError(
            "No finite Step-7 cross_combined configuration curves."
        )

    d["bin"] = d["bin"].astype(int)
    return d.sort_values(
        ["configuration_id", "bin"]
    ).reset_index(drop=True)


def derive_stable_pseudo_absolute_support(
    all_curves: pd.DataFrame,
    full_ensemble_n: int,
    min_bin_coverage_fraction: float,
) -> pd.DataFrame:
    """
    Define stable absolute support using the fraction of the FULL Step-7
    ensemble with a valid cross_combined estimate in each native bin.
    """
    rows = []

    for bin_id, g in all_curves.groupby("bin", sort=True):
        g_valid = g[g["meets_min_n_bool"].astype(bool)].copy()
        x_center = finite_numeric(g["x_center"])
        y = finite_numeric(g_valid["mean_dx"])

        n_present = int(g["configuration_id"].astype(str).nunique())
        n_valid = int(g_valid["configuration_id"].astype(str).nunique())
        coverage = float(n_valid / float(full_ensemble_n))

        rows.append(
            {
                "bin": int(bin_id),
                "x_center": (
                    float(np.median(x_center))
                    if x_center.size else np.nan
                ),
                "n_configurations_present": n_present,
                "n_configurations_valid": n_valid,
                "valid_coverage_fraction": coverage,
                "pseudo_median_dx": (
                    float(np.median(y)) if y.size else np.nan
                ),
                "pseudo_q025": (
                    float(np.quantile(y, 0.025)) if y.size else np.nan
                ),
                "pseudo_q975": (
                    float(np.quantile(y, 0.975)) if y.size else np.nan
                ),
            }
        )

    support = pd.DataFrame(rows)
    support = support[
        np.isfinite(support["x_center"])
        & (
            support["valid_coverage_fraction"]
            >= float(min_bin_coverage_fraction)
        )
    ].sort_values("bin")

    if len(support) < 3:
        raise ValueError(
            "Fewer than three Step-7 native bins satisfy the full-ensemble "
            f"valid-coverage threshold {min_bin_coverage_fraction:.3f}."
        )

    return support.reset_index(drop=True)

def absolute_branch(
    real_curve: pd.DataFrame,
    pseudo_curves_all: pd.DataFrame,
    stable_support: pd.DataFrame,
    selected_ids: Sequence[str],
    outdir: Path,
    args: argparse.Namespace,
    execution_min_null_configurations: int,
    inferential_min_null_configurations: int,
) -> Dict[str, object]:
    outdir = ensure_dir(outdir)
    n_selected = len(selected_ids)

    real_support_min = float(real_curve["x_center"].min())
    real_support_max = float(real_curve["x_center"].max())

    stable = stable_support[
        (stable_support["x_center"] >= real_support_min)
        & (stable_support["x_center"] <= real_support_max)
    ].copy()
    if len(stable) < int(args.absolute_min_grid_points):
        raise ValueError(
            "Too few stable pseudo native bins lie inside the supported "
            "longitudinal Step-9 abundance range: "
            f"stable_shared={len(stable)}, "
            f"required={int(args.absolute_min_grid_points)}."
        )

    grid = stable["x_center"].to_numpy(float)
    shared_min = float(np.min(grid))
    shared_max = float(np.max(grid))

    real_y = interpolate_inside_support(
        real_curve["x_center"].to_numpy(float),
        real_curve["mean_dx"].to_numpy(float),
        grid,
    )
    real_lo = interpolate_inside_support(
        real_curve["x_center"].to_numpy(float),
        real_curve["mean_dx_ci025"].to_numpy(float),
        grid,
    )
    real_hi = interpolate_inside_support(
        real_curve["x_center"].to_numpy(float),
        real_curve["mean_dx_ci975"].to_numpy(float),
        grid,
    )
    if np.sum(np.isfinite(real_y)) < 3:
        raise ValueError(
            "Longitudinal curve has insufficient support on stable pseudo bins."
        )

    selected = pseudo_curves_all[
        pseudo_curves_all["configuration_id"].astype(str).isin(selected_ids)
    ].copy()
    grouped = {
        cid: g
        for cid, g in selected.groupby("configuration_id", sort=False)
    }

    matrix = np.full((n_selected, len(grid)), np.nan, dtype=float)
    metric_rows = []

    stable_bins = stable["bin"].astype(int).tolist()
    bin_to_pos = {
        int(bin_id): int(i)
        for i, bin_id in enumerate(stable_bins)
    }

    for i, cid in enumerate(selected_ids):
        g = grouped.get(cid)
        if g is None:
            continue

        # Direct native-bin observations only. There is deliberately no
        # interpolation, imputation, or connection across missing pseudo bins.
        g = g[
            g["meets_min_n_bool"].astype(bool)
            & g["bin"].isin(stable_bins)
        ].copy()

        if g["bin"].duplicated().any():
            duplicates = (
                g.loc[g["bin"].duplicated(keep=False), "bin"]
                .astype(int)
                .unique()
                .tolist()
            )
            raise ValueError(
                f"{cid}: duplicate Step-7 cross_combined native bins: "
                f"{duplicates[:10]}"
            )

        y = np.full(len(grid), np.nan, dtype=float)
        for row in g.itertuples(index=False):
            bin_id = int(row.bin)
            if bin_id not in bin_to_pos:
                continue
            value = float(row.mean_dx)
            if np.isfinite(value):
                y[bin_to_pos[bin_id]] = value

        finite = np.isfinite(y)
        coverage_fraction = float(np.mean(finite))
        n_finite = int(finite.sum())

        eligible = (
            coverage_fraction >= float(args.absolute_min_grid_coverage)
            and n_finite >= int(args.absolute_min_grid_points)
        )
        if not eligible:
            continue

        matrix[i, :] = y
        metrics = curve_metrics(grid, y)
        metric_rows.append(
            {
                "configuration_id": cid,
                "grid_coverage_fraction": coverage_fraction,
                "n_grid_points_finite": n_finite,
                "n_grid_points_required": int(len(grid)),
                "complete_shared_grid": bool(n_finite == len(grid)),
                **metrics,
            }
        )

    pseudo_metrics = pd.DataFrame(metric_rows)
    if len(pseudo_metrics) < int(execution_min_null_configurations):
        raise ValueError(
            "Too few pseudo configurations satisfy absolute stable-bin support: "
            f"eligible={len(pseudo_metrics)}, required_for_execution="
            f"{execution_min_null_configurations}."
        )

    pseudo_median = np.nanmedian(matrix, axis=0)
    pseudo_q025 = np.nanquantile(matrix, 0.025, axis=0)
    pseudo_q975 = np.nanquantile(matrix, 0.975, axis=0)
    pseudo_n_point = np.sum(np.isfinite(matrix), axis=0)

    shared = pd.DataFrame(
        {
            "bin": stable["bin"].to_numpy(int),
            "x": grid,
            "longitudinal_mean_dx": real_y,
            "longitudinal_ci025": real_lo,
            "longitudinal_ci975": real_hi,
            "pseudo_median_dx": pseudo_median,
            "pseudo_q025": pseudo_q025,
            "pseudo_q975": pseudo_q975,
            "pseudo_n_configurations": pseudo_n_point,
            "full_ensemble_valid_coverage_fraction":
                stable["valid_coverage_fraction"].to_numpy(float),
        }
    )
    shared["longitudinal_minus_pseudo"] = (
        shared["longitudinal_mean_dx"] - shared["pseudo_median_dx"]
    )
    shared["longitudinal_outside_pseudo_95"] = (
        (shared["longitudinal_mean_dx"] < shared["pseudo_q025"])
        | (shared["longitudinal_mean_dx"] > shared["pseudo_q975"])
    )

    observed_metrics = curve_metrics(
        shared["x"].to_numpy(float),
        shared["longitudinal_mean_dx"].to_numpy(float),
    )
    tests = empirical_tests(
        observed_metrics,
        pseudo_metrics,
        branch="absolute_abundance",
        min_null_configurations=int(inferential_min_null_configurations),
    )

    native_real = real_curve[
        ["x_center", "mean_dx", "mean_dx_ci025", "mean_dx_ci975"]
    ].copy()
    native_real.insert(0, "dataset", "longitudinal")

    pseudo_native = stable_support.copy()
    pseudo_native.insert(0, "dataset", "pseudo")
    pseudo_native = pseudo_native.rename(
        columns={
            "pseudo_median_dx": "mean_dx",
            "pseudo_q025": "q025",
            "pseudo_q975": "q975",
        }
    )
    native = pd.concat(
        [native_real, pseudo_native],
        ignore_index=True,
        sort=False,
    )

    stable_support_export = stable_support.copy()
    stable_support_export["inside_longitudinal_native_support"] = (
        (stable_support_export["x_center"] >= real_support_min)
        & (stable_support_export["x_center"] <= real_support_max)
    )
    stable_support_export["used_in_shared_absolute_analysis"] = (
        stable_support_export["bin"].astype(int).isin(
            stable["bin"].astype(int)
        )
    )
    stable_support_export["minimum_valid_coverage_required"] = float(
        args.absolute_min_pseudo_bin_coverage
    )

    support = pd.DataFrame(
        [
            {
                "dataset": "longitudinal",
                "native_support_min": real_support_min,
                "native_support_max": real_support_max,
                "shared_support_min": shared_min,
                "shared_support_max": shared_max,
                "n_native_bins": int(len(real_curve)),
                "n_stable_shared_bins": int(len(stable)),
                "n_selected_pseudo_configurations": n_selected,
                "n_eligible_pseudo_configurations": int(len(pseudo_metrics)),
            },
            {
                "dataset": "pseudo",
                "native_support_min": float(stable_support["x_center"].min()),
                "native_support_max": float(stable_support["x_center"].max()),
                "shared_support_min": shared_min,
                "shared_support_max": shared_max,
                "n_native_bins": int(len(stable_support)),
                "n_stable_shared_bins": int(len(stable)),
                "n_selected_pseudo_configurations": n_selected,
                "n_eligible_pseudo_configurations": int(len(pseudo_metrics)),
            },
        ]
    )

    n_outside = int(
        shared["longitudinal_outside_pseudo_95"].astype(bool).sum()
    )
    summary = pd.DataFrame(
        [
            {
                "branch": "absolute_abundance",
                **{f"longitudinal_{k}": v for k, v in observed_metrics.items()},
                "shared_support_min": shared_min,
                "shared_support_max": shared_max,
                "n_shared_grid_points": int(len(shared)),
                "n_outside_pseudo_95": n_outside,
                "fraction_outside_pseudo_95": float(n_outside / len(shared)),
                "n_pseudo_selected": n_selected,
                "n_pseudo_eligible": int(len(pseudo_metrics)),
            }
        ]
    )

    native.to_csv(outdir / "01_native_curves.csv", index=False)
    shared.to_csv(outdir / "02_shared_grid.csv", index=False)
    pseudo_metrics.to_csv(outdir / "03_pseudo_null_metrics.csv", index=False)
    tests.to_csv(outdir / "04_empirical_tests.csv", index=False)
    support.to_csv(outdir / "05_support_summary.csv", index=False)
    summary.to_csv(outdir / "06_comparison_summary.csv", index=False)
    stable_support_export.to_csv(
        outdir / "07_pseudo_native_bin_support.csv",
        index=False,
    )

    return {
        "shared": shared,
        "tests": tests,
        "summary": summary,
        "observed_metrics": observed_metrics,
        "pseudo_metrics": pseudo_metrics,
    }


# =============================================================================
# Percentile geometry helpers
# =============================================================================

def midrank_percentile(x: np.ndarray) -> np.ndarray:
    """
    Mid-rank ECDF coordinate:
        p = (average_rank - 0.5) / n
    with average ranks for ties.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)

    if n == 0:
        return np.asarray([], dtype=float)

    order = np.argsort(
        x,
        kind="mergesort",
    )
    xs = x[order]

    starts = np.r_[
        0,
        np.flatnonzero(xs[1:] != xs[:-1]) + 1,
    ]
    ends = np.r_[
        starts[1:],
        n,
    ]

    mid_p = (starts + ends) / (2.0 * n)
    repeated = np.repeat(
        mid_p,
        ends - starts,
    )

    out = np.empty(n, dtype=float)
    out[order] = repeated
    return out


def percentile_bin_index(
    percentile: np.ndarray,
    n_bins: int,
) -> np.ndarray:
    p = np.asarray(percentile, dtype=float)
    idx = np.floor(
        p * int(n_bins)
    ).astype(np.int32)
    return np.clip(
        idx,
        0,
        int(n_bins) - 1,
    )


def fold_percentile_stats_from_arrays(
    x: np.ndarray,
    dx: np.ndarray,
    n_bins: int,
    min_n: int,
) -> pd.DataFrame:
    x = np.asarray(x, dtype=float)
    dx = np.asarray(dx, dtype=float)
    ok = np.isfinite(x) & np.isfinite(dx)
    x = x[ok]
    dx = dx[ok]

    if len(x) == 0:
        return pd.DataFrame()

    p = midrank_percentile(x)
    b = percentile_bin_index(
        p,
        n_bins=int(n_bins),
    )

    rows = []
    for bin_id in range(int(n_bins)):
        mask = b == bin_id
        n = int(np.sum(mask))
        if n == 0:
            continue

        xx = x[mask]
        yy = dx[mask]
        rows.append(
            {
                "bin": int(bin_id),
                "percentile_left": float(
                    bin_id / n_bins
                ),
                "percentile_right": float(
                    (bin_id + 1) / n_bins
                ),
                "percentile_center": float(
                    (bin_id + 0.5) / n_bins
                ),
                "n": n,
                "mean_initial_x": float(
                    np.mean(xx)
                ),
                "median_initial_x": float(
                    np.median(xx)
                ),
                "mean_dx": float(
                    np.mean(yy)
                ),
                "median_dx": float(
                    np.median(yy)
                ),
                "p_dx_gt0": float(
                    np.mean(yy > 0)
                ),
                "sd_dx": (
                    float(np.std(yy, ddof=1))
                    if n > 1
                    else np.nan
                ),
                "meets_min_n": bool(
                    n >= int(min_n)
                ),
            }
        )

    return pd.DataFrame(rows)


def combine_percentile_folds(
    ab: pd.DataFrame,
    ba: pd.DataFrame,
    min_n: int,
) -> pd.DataFrame:
    if ab.empty or ba.empty:
        return pd.DataFrame()

    a = ab.rename(
        columns={
            "n": "n_AB",
            "mean_initial_x": "mean_initial_x_AB",
            "median_initial_x": "median_initial_x_AB",
            "mean_dx": "mean_dx_AB",
            "median_dx": "median_dx_AB",
            "p_dx_gt0": "p_dx_gt0_AB",
            "meets_min_n": "valid_AB",
        }
    )

    b = ba.rename(
        columns={
            "n": "n_BA",
            "mean_initial_x": "mean_initial_x_BA",
            "median_initial_x": "median_initial_x_BA",
            "mean_dx": "mean_dx_BA",
            "median_dx": "median_dx_BA",
            "p_dx_gt0": "p_dx_gt0_BA",
            "meets_min_n": "valid_BA",
        }
    )

    keep_a = [
        "bin",
        "percentile_left",
        "percentile_right",
        "percentile_center",
        "n_AB",
        "mean_initial_x_AB",
        "median_initial_x_AB",
        "mean_dx_AB",
        "median_dx_AB",
        "p_dx_gt0_AB",
        "valid_AB",
    ]
    keep_b = [
        "bin",
        "n_BA",
        "mean_initial_x_BA",
        "median_initial_x_BA",
        "mean_dx_BA",
        "median_dx_BA",
        "p_dx_gt0_BA",
        "valid_BA",
    ]

    m = a[keep_a].merge(
        b[keep_b],
        on="bin",
        how="inner",
    )

    m["meets_min_n"] = (
        m["valid_AB"].astype(bool)
        & m["valid_BA"].astype(bool)
    )
    m["n_min_fold"] = np.minimum(
        m["n_AB"],
        m["n_BA"],
    )
    m["mean_dx"] = 0.5 * (
        m["mean_dx_AB"]
        + m["mean_dx_BA"]
    )
    m["median_dx"] = 0.5 * (
        m["median_dx_AB"]
        + m["median_dx_BA"]
    )
    m["p_dx_gt0"] = 0.5 * (
        m["p_dx_gt0_AB"]
        + m["p_dx_gt0_BA"]
    )
    m["mean_initial_x"] = 0.5 * (
        m["mean_initial_x_AB"]
        + m["mean_initial_x_BA"]
    )

    return m.sort_values("bin").reset_index(drop=True)


def read_fold_transition_rows(
    path: Path,
    dt: int,
    fold: str,
) -> pd.DataFrame:
    spec = PERCENTILE_FOLDS[fold]
    x_col = str(spec["x"])
    dx_col = str(spec["dx"])
    support_col = str(spec["support"])

    lf = (
        pl.scan_parquet(path)
        .filter(
            pl.col("dt").cast(
                pl.Int64,
                strict=False,
            ) == int(dt)
        )
        .select(
            [
                pl.col("subject")
                .cast(pl.String, strict=False)
                .alias("subject"),
                pl.col(x_col)
                .cast(pl.Float64, strict=False)
                .alias("x"),
                pl.col(dx_col)
                .cast(pl.Float64, strict=False)
                .alias("dx"),
                parse_bool_expr(support_col)
                .alias("eligible"),
            ]
        )
        .filter(
            pl.col("eligible")
            & pl.col("x").is_finite()
            & pl.col("dx").is_finite()
        )
        .select(
            ["subject", "x", "dx"]
        )
    )

    d = collect_frame(lf).to_pandas()
    if d.empty:
        raise ValueError(
            f"No eligible rows for fold {fold} in {path}"
        )
    return d


def longitudinal_percentile_fold(
    transitions: Path,
    dt: int,
    fold: str,
    n_bins: int,
    min_n: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d = read_fold_transition_rows(
        transitions,
        dt=dt,
        fold=fold,
    )

    percentile = np.full(
        len(d),
        np.nan,
        dtype=float,
    )

    # Percentile is defined within biological subject x fold.
    for _, idx in d.groupby(
        "subject",
        sort=False,
    ).groups.items():
        pos = np.asarray(list(idx), dtype=int)
        percentile[pos] = midrank_percentile(
            d.loc[pos, "x"].to_numpy(float)
        )

    d["percentile"] = percentile
    d["bin"] = percentile_bin_index(
        percentile,
        n_bins=n_bins,
    )

    point_rows = []
    for bin_id, g in d.groupby(
        "bin",
        sort=True,
    ):
        yy = g["dx"].to_numpy(float)
        xx = g["x"].to_numpy(float)

        point_rows.append(
            {
                "fold": fold,
                "bin": int(bin_id),
                "percentile_left": float(
                    int(bin_id) / n_bins
                ),
                "percentile_right": float(
                    (int(bin_id) + 1) / n_bins
                ),
                "percentile_center": float(
                    (int(bin_id) + 0.5) / n_bins
                ),
                "n": int(len(g)),
                "n_subjects": int(
                    g["subject"].nunique()
                ),
                "mean_initial_x": float(
                    np.mean(xx)
                ),
                "median_initial_x": float(
                    np.median(xx)
                ),
                "mean_dx": float(
                    np.mean(yy)
                ),
                "median_dx": float(
                    np.median(yy)
                ),
                "p_dx_gt0": float(
                    np.mean(yy > 0)
                ),
                "sd_dx": float(
                    np.std(yy, ddof=1)
                ) if len(yy) > 1 else np.nan,
                "meets_min_n": bool(
                    len(g) >= int(min_n)
                ),
            }
        )

    # Subject x bin sufficient statistics for biological bootstrap.
    suff = (
        d.groupby(
            ["subject", "bin"],
            as_index=False,
            sort=True,
        )
        .agg(
            n=("dx", "size"),
            sum_dx=("dx", "sum"),
            n_pos=("dx", lambda x: int(np.sum(np.asarray(x) > 0))),
        )
    )
    suff["fold"] = fold

    return pd.DataFrame(point_rows), suff


def bootstrap_longitudinal_percentile(
    suff_ab: pd.DataFrame,
    suff_ba: pd.DataFrame,
    n_bins: int,
    n_bootstrap: int,
    seed: int,
    min_valid_bins: int,
    fixed_valid_bins: Sequence[int],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    subjects = sorted(
        set(suff_ab["subject"].astype(str))
        | set(suff_ba["subject"].astype(str))
    )
    if len(subjects) < 2:
        raise ValueError(
            "Need at least two biological subjects for percentile bootstrap."
        )

    fixed_bins = sorted(set(int(b) for b in fixed_valid_bins))
    if len(fixed_bins) < int(min_valid_bins):
        raise ValueError(
            "Point-estimate percentile curve has too few fixed valid bins "
            "for bootstrap linear descriptors."
        )

    si = {subject: i for i, subject in enumerate(subjects)}

    def arrays(frame: pd.DataFrame):
        n = np.zeros((len(subjects), n_bins), dtype=float)
        s = np.zeros_like(n)
        p = np.zeros_like(n)
        for row in frame.itertuples(index=False):
            i = si[str(row.subject)]
            b = int(row.bin)
            n[i, b] = float(row.n)
            s[i, b] = float(row.sum_dx)
            p[i, b] = float(row.n_pos)
        return n, s, p

    nA, sA, pA = arrays(suff_ab)
    nB, sB, pB = arrays(suff_ba)
    rng = np.random.default_rng(int(seed))
    curve_rows = []
    linear_rows = []
    centers = (np.arange(n_bins, dtype=float) + 0.5) / float(n_bins)

    for boot in range(int(n_bootstrap)):
        draw = rng.integers(0, len(subjects), size=len(subjects))
        mult = np.bincount(
            draw, minlength=len(subjects)
        ).astype(float)[:, None]

        n_ab = np.sum(mult * nA, axis=0)
        s_ab = np.sum(mult * sA, axis=0)
        p_ab = np.sum(mult * pA, axis=0)
        n_ba = np.sum(mult * nB, axis=0)
        s_ba = np.sum(mult * sB, axis=0)
        p_ba = np.sum(mult * pB, axis=0)

        mean_ab = np.divide(
            s_ab, n_ab, out=np.full(n_bins, np.nan), where=n_ab > 0
        )
        mean_ba = np.divide(
            s_ba, n_ba, out=np.full(n_bins, np.nan), where=n_ba > 0
        )
        ppos_ab = np.divide(
            p_ab, n_ab, out=np.full(n_bins, np.nan), where=n_ab > 0
        )
        ppos_ba = np.divide(
            p_ba, n_ba, out=np.full(n_bins, np.nan), where=n_ba > 0
        )

        combined = 0.5 * (mean_ab + mean_ba)
        combined_p = 0.5 * (ppos_ab + ppos_ba)

        for b in fixed_bins:
            curve_rows.append(
                {
                    "bootstrap": int(boot),
                    "bin": int(b),
                    "percentile_center": float(centers[b]),
                    "mean_dx": float(combined[b]),
                    "p_dx_gt0": float(combined_p[b]),
                }
            )

        valid_bins = [b for b in fixed_bins if np.isfinite(combined[b])]
        if len(valid_bins) >= int(min_valid_bins):
            linear = safe_linear_descriptor(
                centers[valid_bins], combined[valid_bins]
            )
            linear_rows.append(
                {
                    "bootstrap": int(boot),
                    "n_bins_valid": int(len(valid_bins)),
                    "slope": linear["slope"],
                    "intercept": linear["intercept"],
                    "crossover_percentile": linear["x_zero"],
                }
            )

    return pd.DataFrame(curve_rows), pd.DataFrame(linear_rows)

def longitudinal_percentile_curve(
    transitions: Path,
    dt: int,
    n_bins: int,
    min_n: int,
    n_bootstrap: int,
    bootstrap_seed: int,
    min_valid_bins: int,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    ab, suff_ab = longitudinal_percentile_fold(
        transitions,
        dt=dt,
        fold="AB",
        n_bins=n_bins,
        min_n=min_n,
    )
    ba, suff_ba = longitudinal_percentile_fold(
        transitions,
        dt=dt,
        fold="BA",
        n_bins=n_bins,
        min_n=min_n,
    )

    combined = combine_percentile_folds(
        ab,
        ba,
        min_n=min_n,
    )
    combined = combined[
        combined["meets_min_n"].astype(bool)
    ].copy()

    if len(combined) < int(min_valid_bins):
        raise ValueError(
            "Longitudinal percentile curve has too few valid bins."
        )

    fixed_valid_bins = combined["bin"].astype(int).tolist()

    boot_curves, boot_linear = bootstrap_longitudinal_percentile(
        suff_ab,
        suff_ba,
        n_bins=n_bins,
        n_bootstrap=n_bootstrap,
        seed=bootstrap_seed,
        min_valid_bins=min_valid_bins,
        fixed_valid_bins=fixed_valid_bins,
    )

    ci_rows = []
    for bin_id, g in boot_curves.groupby(
        "bin",
        sort=True,
    ):
        x = finite_numeric(g["mean_dx"])
        p = finite_numeric(g["p_dx_gt0"])
        ci_rows.append(
            {
                "bin": int(bin_id),
                "mean_dx_ci025": (
                    float(np.quantile(x, 0.025))
                    if x.size else np.nan
                ),
                "mean_dx_ci975": (
                    float(np.quantile(x, 0.975))
                    if x.size else np.nan
                ),
                "p_dx_gt0_ci025": (
                    float(np.quantile(p, 0.025))
                    if p.size else np.nan
                ),
                "p_dx_gt0_ci975": (
                    float(np.quantile(p, 0.975))
                    if p.size else np.nan
                ),
                "n_bootstrap_valid": int(
                    x.size
                ),
            }
        )

    combined = combined.merge(
        pd.DataFrame(ci_rows),
        on="bin",
        how="left",
    )

    return combined, boot_linear, boot_curves


def pseudo_configuration_percentile_curve(
    path: Path,
    configuration_id: str,
    dt: int,
    n_bins: int,
    min_n: int,
) -> pd.DataFrame:
    validate_transition_schema(
        path,
        f"pseudo {configuration_id}",
    )

    fold_tables = {}
    for fold, spec in PERCENTILE_FOLDS.items():
        x_col = str(spec["x"])
        dx_col = str(spec["dx"])
        support_col = str(spec["support"])

        lf = (
            pl.scan_parquet(path)
            .filter(
                pl.col("dt").cast(
                    pl.Int64,
                    strict=False,
                ) == int(dt)
            )
            .select(
                [
                    pl.col(x_col)
                    .cast(pl.Float64, strict=False)
                    .alias("x"),
                    pl.col(dx_col)
                    .cast(pl.Float64, strict=False)
                    .alias("dx"),
                    parse_bool_expr(support_col)
                    .alias("eligible"),
                ]
            )
            .filter(
                pl.col("eligible")
                & pl.col("x").is_finite()
                & pl.col("dx").is_finite()
            )
            .select(["x", "dx"])
        )

        d = collect_frame(lf)
        x = d["x"].to_numpy().astype(float)
        dx = d["dx"].to_numpy().astype(float)

        fold_tables[fold] = fold_percentile_stats_from_arrays(
            x,
            dx,
            n_bins=n_bins,
            min_n=min_n,
        )

    combined = combine_percentile_folds(
        fold_tables["AB"],
        fold_tables["BA"],
        min_n=min_n,
    )

    if combined.empty:
        return combined

    combined.insert(
        0,
        "configuration_id",
        configuration_id,
    )
    return combined


def percentile_cache_signature(
    step7_cfg: Mapping[str, object],
    selected_ids: Sequence[str],
    args: argparse.Namespace,
) -> Dict[str, object]:
    payload = {
        "script_version": PERCENTILE_CACHE_VERSION,
        "step7_analysis_signature": step7_cfg.get(
            "analysis_signature"
        ),
        "selected_configuration_ids": list(
            selected_ids
        ),
        "dt": int(args.dt),
        "percentile_bins": int(
            args.percentile_bins
        ),
        "percentile_min_n": int(
            args.percentile_min_n
        ),
    }
    return {
        "signature": stable_sha256(payload),
        "payload": payload,
    }


def build_or_reuse_pseudo_percentile_curves(
    compact_map: Mapping[str, Path],
    selected_ids: Sequence[str],
    step7_cfg: Mapping[str, object],
    outdir: Path,
    args: argparse.Namespace,
) -> pd.DataFrame:
    curves_path = (
        outdir / "02_pseudo_configuration_percentile_curves.parquet"
    )
    signature_path = (
        outdir / "00_percentile_cache_signature.json"
    )

    current = percentile_cache_signature(
        step7_cfg,
        selected_ids,
        args,
    )

    if (
        not args.restart_percentile_cache
        and curves_path.is_file()
        and signature_path.is_file()
    ):
        previous = read_json(signature_path, required=False)
        if previous.get("signature") == current["signature"]:
            print(
                "[INFO] reusing current percentile pseudo cache:",
                curves_path,
            )
            return pd.read_parquet(curves_path)

    if curves_path.exists():
        curves_path.unlink()

    rows = []
    for i, cid in enumerate(selected_ids, start=1):
        if cid not in compact_map:
            raise KeyError(
                f"{cid}: absent from Step-7 compact manifest"
            )
        path = Path(compact_map[cid])
        curve = pseudo_configuration_percentile_curve(
            path,
            configuration_id=cid,
            dt=int(args.dt),
            n_bins=int(args.percentile_bins),
            min_n=int(args.percentile_min_n),
        )
        if not curve.empty:
            rows.append(curve)

        if i == 1 or i % 50 == 0 or i == len(selected_ids):
            print(
                f"[percentile pseudo compact] {i}/{len(selected_ids)}"
            )

    if not rows:
        raise RuntimeError(
            "No pseudo percentile configuration curves were generated."
        )

    out = pd.concat(rows, ignore_index=True, sort=False)
    out.to_parquet(curves_path, index=False)
    signature_path.write_text(
        json.dumps(current, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return out

def percentile_branch(
    longitudinal_curve: pd.DataFrame,
    longitudinal_boot_linear: pd.DataFrame,
    pseudo_curves: pd.DataFrame,
    outdir: Path,
    args: argparse.Namespace,
    execution_min_null_configurations: int,
    inferential_min_null_configurations: int,
) -> Dict[str, object]:
    outdir = ensure_dir(outdir)

    centers = (
        np.arange(
            int(args.percentile_bins),
            dtype=float,
        )
        + 0.5
    ) / float(args.percentile_bins)

    # Configuration eligibility and interpolation to the fixed percentile grid.
    selected_ids = sorted(
        pseudo_curves["configuration_id"]
        .astype(str)
        .unique()
        .tolist()
    )

    matrix = np.full(
        (
            len(selected_ids),
            len(centers),
        ),
        np.nan,
        dtype=float,
    )
    metric_rows = []

    grouped = {
        cid: g
        for cid, g in pseudo_curves.groupby(
            "configuration_id",
            sort=False,
        )
    }

    for i, cid in enumerate(selected_ids):
        g = grouped[cid].copy()
        g = g[
            g["meets_min_n"].astype(bool)
        ].sort_values("percentile_center")

        if len(g) < int(
            args.percentile_min_valid_bins
        ):
            continue

        y = interpolate_inside_support(
            g["percentile_center"].to_numpy(float),
            g["mean_dx"].to_numpy(float),
            centers,
        )
        finite = np.isfinite(y)
        coverage = float(np.mean(finite))
        n_finite = int(finite.sum())

        eligible = (
            coverage
            >= float(args.percentile_min_grid_coverage)
            and n_finite
            >= int(args.percentile_min_grid_points)
        )
        if not eligible:
            continue

        matrix[i, :] = y
        metrics = curve_metrics(
            centers,
            y,
        )
        linear = safe_linear_descriptor(
            centers,
            y,
        )

        metric_rows.append(
            {
                "configuration_id": cid,
                "grid_coverage_fraction": coverage,
                "n_grid_points_finite": n_finite,
                **metrics,
                "intercept": linear["intercept"],
                "crossover_percentile": linear[
                    "x_zero"
                ],
            }
        )

    pseudo_metrics = pd.DataFrame(metric_rows)
    if len(pseudo_metrics) < int(
        execution_min_null_configurations
    ):
        raise ValueError(
            "Too few pseudo configurations satisfy percentile-grid support: "
            f"eligible={len(pseudo_metrics)}, required_for_execution="
            f"{execution_min_null_configurations}."
        )

    envelope = pd.DataFrame(
        {
            "bin": np.arange(
                len(centers),
                dtype=int,
            ),
            "percentile_center": centers,
            "pseudo_median_dx": np.nanmedian(
                matrix,
                axis=0,
            ),
            "pseudo_q025": np.nanquantile(
                matrix,
                0.025,
                axis=0,
            ),
            "pseudo_q975": np.nanquantile(
                matrix,
                0.975,
                axis=0,
            ),
            "pseudo_n_configurations": np.sum(
                np.isfinite(matrix),
                axis=0,
            ),
        }
    )

    real = longitudinal_curve.copy()
    real = real.sort_values("bin")
    pointwise = real.merge(
        envelope,
        on=[
            "bin",
            "percentile_center",
        ],
        how="inner",
    )
    pointwise["longitudinal_minus_pseudo"] = (
        pointwise["mean_dx"]
        - pointwise["pseudo_median_dx"]
    )
    pointwise["longitudinal_outside_pseudo_95"] = (
        (
            pointwise["mean_dx"]
            < pointwise["pseudo_q025"]
        )
        |
        (
            pointwise["mean_dx"]
            > pointwise["pseudo_q975"]
        )
    )

    real_metrics = curve_metrics(
        pointwise["percentile_center"].to_numpy(float),
        pointwise["mean_dx"].to_numpy(float),
    )
    real_linear = safe_linear_descriptor(
        pointwise["percentile_center"].to_numpy(float),
        pointwise["mean_dx"].to_numpy(float),
    )

    tests = empirical_tests(
        real_metrics,
        pseudo_metrics,
        branch="percentile_geometry",
        min_null_configurations=int(
            inferential_min_null_configurations
        ),
    )

    crossover_boot = finite_numeric(
        longitudinal_boot_linear.get(
            "crossover_percentile",
            pd.Series(dtype=float),
        )
    )
    slope_boot = finite_numeric(
        longitudinal_boot_linear.get(
            "slope",
            pd.Series(dtype=float),
        )
    )

    crossover_q025 = (
        float(np.quantile(crossover_boot, 0.025))
        if crossover_boot.size else np.nan
    )
    crossover_q975 = (
        float(np.quantile(crossover_boot, 0.975))
        if crossover_boot.size else np.nan
    )
    slope_q025 = (
        float(np.quantile(slope_boot, 0.025))
        if slope_boot.size else np.nan
    )
    slope_q975 = (
        float(np.quantile(slope_boot, 0.975))
        if slope_boot.size else np.nan
    )

    n_outside = int(
        pointwise[
            "longitudinal_outside_pseudo_95"
        ].astype(bool).sum()
    )

    summary = pd.DataFrame(
        [
            {
                "branch": "percentile_geometry",
                **{
                    f"longitudinal_{k}": v
                    for k, v in real_metrics.items()
                },
                "longitudinal_intercept": real_linear[
                    "intercept"
                ],
                "longitudinal_crossover_percentile": real_linear[
                    "x_zero"
                ],
                "longitudinal_slope_bootstrap_q025": slope_q025,
                "longitudinal_slope_bootstrap_q975": slope_q975,
                "longitudinal_crossover_bootstrap_q025": crossover_q025,
                "longitudinal_crossover_bootstrap_q975": crossover_q975,
                "n_percentile_bins": int(
                    len(pointwise)
                ),
                "n_outside_pseudo_95": n_outside,
                "fraction_outside_pseudo_95": float(
                    n_outside / len(pointwise)
                ),
                "n_pseudo_eligible": int(
                    len(pseudo_metrics)
                ),
            }
        ]
    )

    support = pd.DataFrame(
        [
            {
                "dataset": "longitudinal",
                "percentile_definition": (
                    "transition-level midrank within biological subject x fold"
                ),
                "n_bins": int(
                    len(longitudinal_curve)
                ),
                "min_fold_n_per_bin": int(
                    min(
                        longitudinal_curve["n_AB"].min(),
                        longitudinal_curve["n_BA"].min(),
                    )
                ),
                "n_pseudo_configurations_eligible": np.nan,
            },
            {
                "dataset": "pseudo",
                "percentile_definition": (
                    "transition-level midrank within configuration x fold"
                ),
                "n_bins": int(
                    args.percentile_bins
                ),
                "min_fold_n_per_bin": int(
                    args.percentile_min_n
                ),
                "n_pseudo_configurations_eligible": int(
                    len(pseudo_metrics)
                ),
            },
        ]
    )

    longitudinal_curve.to_csv(
        outdir / "01_longitudinal_percentile_curve.csv",
        index=False,
    )
    # 02 is the cache written upstream.
    envelope.to_csv(
        outdir / "03_pseudo_percentile_envelope.csv",
        index=False,
    )
    pseudo_metrics.to_csv(
        outdir / "04_pseudo_null_metrics.csv",
        index=False,
    )
    tests.to_csv(
        outdir / "05_empirical_tests.csv",
        index=False,
    )
    pointwise.to_csv(
        outdir / "06_pointwise_comparison.csv",
        index=False,
    )
    longitudinal_boot_linear.to_csv(
        outdir / "07_longitudinal_bootstrap_linear_descriptors.csv",
        index=False,
    )
    support.to_csv(
        outdir / "08_support_summary.csv",
        index=False,
    )
    summary.to_csv(
        outdir / "09_comparison_summary.csv",
        index=False,
    )

    return {
        "pointwise": pointwise,
        "tests": tests,
        "summary": summary,
        "pseudo_metrics": pseudo_metrics,
        "real_metrics": real_metrics,
        "real_linear": real_linear,
    }


# =============================================================================
# README
# =============================================================================

def write_readme(outdir: Path) -> None:
    text = """# Step 10 — genuine longitudinal versus pseudo forward technical null

Step 10 compares the same observed replicate-decoupled AB/BA forward estimator
in the genuine longitudinal cohort and the 1x12 pseudo-longitudinal technical
reference.

No TT/TF/FT/FF filter is introduced.

## 01_absolute_abundance

Native Step-9 and Step-7 abundance geometries are retained. Comparison is
restricted to their shared absolute observed-log-frequency support and no
extrapolation is allowed.

This branch provides the direct quantitative comparison at the same absolute
abundance scale represented by both datasets.

## 02_percentile_geometry

Initial observed abundance is converted to a transition-level midrank percentile:

- longitudinal: within biological subject x fold;
- pseudo: within randomization configuration x fold.

AB and BA are transformed separately, then estimated separately, and combined
with equal fold weight only after binning.

This coordinate compares relative abundance geometry across the full supported
range. Equal percentile does NOT imply equal absolute frequency across datasets.

## Inference

Longitudinal confidence intervals are biological subject-cluster bootstrap
intervals.

Pseudo intervals are randomization intervals across arbitrary 1x12
configurations from one biological source sample and are not biological
confidence intervals.

For the absolute empirical null, Step-7 configuration curves are used on
their fixed ensemble-wide abundance grid. Only configuration/bin cells with
`meets_min_n=True` contribute to the technical null. Stable absolute support is
defined from the fraction of the full 2,000-configuration ensemble with a valid
combined estimate in each native Step-7 bin (production default: >=95%).

Missing pseudo native bins are NEVER interpolated or imputed. In production, a
pseudo configuration contributes a curve-level null metric only if every
stable shared native bin is directly valid. Thus the genuine curve and every
eligible pseudo curve are summarized on exactly the same abundance coordinates.

The percentile pseudo branch reads the compact Step-7 `step7_dt1.parquet`
files directly. It does not require or reconstruct deleted full Step-5 pseudo
transition tables.

In `--max-pseudo-configs` smoke-test mode, the pipeline is allowed to
complete when at least one pseudo configuration is eligible. The inferential
minimum is not relaxed: empirical P-values remain NA until the requested
minimum number of eligible configurations is reached. Smoke-test results are
therefore structural diagnostics only.

Empirical P-values use a finite-sample +1 correction.

The main directional test is the lower-tail slope test. Centered magnitude
metrics and mean displacement provide complementary curve-level comparisons.
"""
    (outdir / "README_outputs.md").write_text(
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
            "Step 10: compare genuine observed replicate-decoupled forward "
            "dynamics with the pseudo 1x12 technical null on absolute and "
            "within-unit percentile abundance geometries."
        ),
    )

    p.add_argument(
        "--longitudinal-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--pseudo-step7-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    p.add_argument(
        "--expected-pseudo-configs",
        type=int,
        default=2000,
    )
    p.add_argument(
        "--max-pseudo-configs",
        type=int,
        default=None,
        help="Testing only: use the first N pseudo configurations.",
    )

    p.add_argument(
        "--dt",
        type=int,
        default=1,
    )

    # Absolute branch.
    p.add_argument(
        "--absolute-min-pseudo-bin-coverage",
        type=float,
        default=0.95,
        help=(
            "Minimum fraction of the full Step-7 ensemble with "
            "meets_min_n=True in a native cross_combined bin."
        ),
    )
    p.add_argument(
        "--absolute-min-grid-coverage",
        type=float,
        default=1.0,
        help=(
            "Fraction of stable shared native bins that must be directly valid "
            "within a pseudo configuration. Production inference requires 1.0; "
            "lower values are allowed only with --max-pseudo-configs."
        ),
    )
    p.add_argument(
        "--absolute-min-grid-points",
        type=int,
        default=11,
    )

    # Percentile branch.
    p.add_argument(
        "--percentile-bins",
        type=int,
        default=20,
    )
    p.add_argument(
        "--percentile-min-n",
        type=int,
        default=50,
    )
    p.add_argument(
        "--percentile-min-valid-bins",
        type=int,
        default=15,
    )
    p.add_argument(
        "--percentile-min-grid-coverage",
        type=float,
        default=0.80,
    )
    p.add_argument(
        "--percentile-min-grid-points",
        type=int,
        default=15,
    )

    p.add_argument(
        "--min-null-configurations",
        type=int,
        default=30,
    )

    p.add_argument(
        "--n-bootstrap",
        type=int,
        default=None,
        help="Percentile longitudinal bootstrap count; default inherits Step 9.",
    )
    p.add_argument(
        "--bootstrap-seed",
        type=int,
        default=None,
        help="Percentile longitudinal bootstrap seed; default inherits Step 9.",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Reserved reproducibility seed for deterministic future extensions.",
    )

    p.add_argument(
        "--restart-percentile-cache",
        action="store_true",
    )

    return p


def validate_args(args: argparse.Namespace) -> None:
    if int(args.dt) < 1:
        raise ValueError("--dt must be >=1")
    if int(args.expected_pseudo_configs) < 0:
        raise ValueError("--expected-pseudo-configs must be >=0")
    if (
        args.max_pseudo_configs is not None
        and int(args.max_pseudo_configs) < 1
    ):
        raise ValueError("--max-pseudo-configs must be >=1")
    if int(args.percentile_bins) < 5:
        raise ValueError("--percentile-bins must be >=5")
    if int(args.percentile_min_n) < 1:
        raise ValueError("--percentile-min-n must be >=1")
    if int(args.percentile_min_valid_bins) < 3:
        raise ValueError("--percentile-min-valid-bins must be >=3")
    if int(args.percentile_min_grid_points) < 3:
        raise ValueError("--percentile-min-grid-points must be >=3")
    if int(args.min_null_configurations) < 1:
        raise ValueError("--min-null-configurations must be >=1")

    for name in [
        "absolute_min_pseudo_bin_coverage",
        "absolute_min_grid_coverage",
        "percentile_min_grid_coverage",
    ]:
        value = float(getattr(args, name))
        if not (0 < value <= 1):
            raise ValueError(
                f"--{name.replace('_','-')} must be in (0,1]"
            )

    if (
        args.max_pseudo_configs is None
        and not np.isclose(
            float(args.absolute_min_grid_coverage),
            1.0,
            rtol=0.0,
            atol=1e-12,
        )
    ):
        raise ValueError(
            "Production absolute inference requires "
            "--absolute-min-grid-coverage 1.0 so every pseudo null curve "
            "uses exactly the same directly valid stable native bins. "
            "Values <1.0 are allowed only with --max-pseudo-configs."
        )


# =============================================================================
# Main
# =============================================================================

def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    args = build_parser().parse_args(argv)
    validate_args(args)

    inputs = resolve_inputs(
        args.longitudinal_dir,
        args.pseudo_step7_dir,
    )

    outdir = ensure_dir(
        args.outdir.expanduser().resolve()
    )
    absolute_dir = ensure_dir(
        outdir / "01_absolute_abundance"
    )
    percentile_dir = ensure_dir(
        outdir / "02_percentile_geometry"
    )

    validate_transition_schema(
        inputs["longitudinal_transitions"],
        "longitudinal",
    )

    step9_cfg = inputs["step9_config"]
    step7_cfg = inputs["step7_config"]

    step9_boot = step9_cfg.get(
        "bootstrap",
        {},
    )
    inherited_n_boot = int(
        step9_boot.get("n_bootstrap", 2000)
    )
    inherited_seed = int(
        step9_boot.get("seed", 123)
    )

    n_bootstrap = (
        int(args.n_bootstrap)
        if args.n_bootstrap is not None
        else inherited_n_boot
    )
    bootstrap_seed = (
        int(args.bootstrap_seed)
        if args.bootstrap_seed is not None
        else inherited_seed
    )

    # Discover / validate the pseudo configuration ensemble from final Step 7.
    step7_curves_all = load_pseudo_absolute_curves(
        inputs["paths"]["step7_curves"]
    )

    all_ids = sorted(
        step7_curves_all["configuration_id"]
        .astype(str)
        .unique()
        .tolist()
    )
    manifest_ids = sorted(inputs["pseudo_compact_map"].keys())
    if all_ids != manifest_ids:
        raise ValueError(
            "Step-7 configuration curves and compact manifest contain "
            "different configuration ID sets."
        )

    if int(args.expected_pseudo_configs) > 0:
        if len(all_ids) != int(args.expected_pseudo_configs):
            raise ValueError(
                f"Expected {args.expected_pseudo_configs} Step-7 configurations; "
                f"found {len(all_ids)}."
            )

    if inputs.get("step7_dt") is not None:
        if int(inputs["step7_dt"]) != int(args.dt):
            raise ValueError(
                f"Step 10 requested dt={args.dt}, but compact Step 7 was "
                f"generated for dt={inputs['step7_dt']}."
            )

    if args.max_pseudo_configs is not None:
        selected_ids = all_ids[
            : int(args.max_pseudo_configs)
        ]
    else:
        selected_ids = all_ids

    testing_subset = args.max_pseudo_configs is not None

    # Separate "can the smoke test execute?" from "is there enough null support
    # for inferential reporting?". In test mode, one eligible configuration is
    # sufficient to exercise the complete pipeline, but empirical P-values still
    # require the user-requested inferential minimum and therefore remain NA when
    # that threshold is not met.
    inferential_min_null_configurations = int(args.min_null_configurations)
    execution_min_null_configurations = (
        1 if testing_subset
        else inferential_min_null_configurations
    )

    if testing_subset:
        print(
            "[TEST MODE] structural smoke test: at least 1 eligible pseudo "
            "configuration is required for execution."
        )
        print(
            "[TEST MODE] empirical P-values still require",
            inferential_min_null_configurations,
            "eligible configurations and will be NA below that threshold."
        )

    print(
        "[INFO] Step 10",
        SCRIPT_VERSION,
    )
    print(
        "[INFO] genuine Step 9:",
        inputs["longitudinal_dir"],
    )
    print(
        "[INFO] pseudo Step 7:",
        inputs["pseudo_step7_dir"],
    )
    print(
        "[INFO] pseudo configurations discovered:",
        len(all_ids),
    )
    print(
        "[INFO] pseudo configurations selected:",
        len(selected_ids),
    )
    print(
        "[INFO] TT restriction: NO",
    )
    print(
        "[INFO] absolute + within-unit percentile geometries",
    )
    print(
        "[INFO] pseudo percentile source: Step-7 compact caches",
    )
    print(
        "[INFO] pseudo full Step-5 reconstruction: NO",
    )
    print(
        "[INFO] absolute pseudo missing-bin interpolation: NO",
    )
    print(
        "[INFO] absolute stable-bin coverage threshold:",
        float(args.absolute_min_pseudo_bin_coverage),
    )
    print(
        "[INFO] absolute per-configuration direct-grid coverage:",
        float(args.absolute_min_grid_coverage),
    )

    # -----------------------------------------------------------------
    # Absolute abundance comparison from frozen Step-9/Step-7 curves.
    # -----------------------------------------------------------------
    real_absolute = load_real_absolute_curve(
        inputs["paths"]["step9_primary"]
    )
    stable_pseudo_support = derive_stable_pseudo_absolute_support(
        step7_curves_all,
        full_ensemble_n=len(all_ids),
        min_bin_coverage_fraction=float(
            args.absolute_min_pseudo_bin_coverage
        ),
    )

    absolute = absolute_branch(
        real_absolute,
        step7_curves_all,
        stable_pseudo_support,
        selected_ids,
        absolute_dir,
        args,
        execution_min_null_configurations,
        inferential_min_null_configurations,
    )

    # -----------------------------------------------------------------
    # Percentile longitudinal curve.
    # -----------------------------------------------------------------
    real_percentile, real_boot_linear, _ = (
        longitudinal_percentile_curve(
            inputs["longitudinal_transitions"],
            dt=int(args.dt),
            n_bins=int(args.percentile_bins),
            min_n=int(args.percentile_min_n),
            n_bootstrap=int(n_bootstrap),
            bootstrap_seed=int(bootstrap_seed),
            min_valid_bins=int(
                args.percentile_min_valid_bins
            ),
        )
    )

    # -----------------------------------------------------------------
    # Percentile pseudo configuration curves.
    # -----------------------------------------------------------------
    pseudo_percentile = (
        build_or_reuse_pseudo_percentile_curves(
            inputs["pseudo_compact_map"],
            selected_ids,
            step7_cfg,
            percentile_dir,
            args,
        )
    )

    percentile = percentile_branch(
        real_percentile,
        real_boot_linear,
        pseudo_percentile,
        percentile_dir,
        args,
        execution_min_null_configurations,
        inferential_min_null_configurations,
    )

    # -----------------------------------------------------------------
    # Top-level comparison summary.
    # -----------------------------------------------------------------
    abs_summary = absolute["summary"].iloc[0].to_dict()
    pct_summary = percentile["summary"].iloc[0].to_dict()

    def p_for(
        tests: pd.DataFrame,
        metric: str,
    ) -> float:
        d = tests[
            tests["metric"].astype(str).eq(metric)
        ]
        if d.empty:
            return np.nan
        return float(d.iloc[0]["empirical_p"])

    analysis_summary = pd.DataFrame(
        [
            {
                "branch": "absolute_abundance",
                "longitudinal_slope": abs_summary.get(
                    "longitudinal_slope"
                ),
                "longitudinal_mean_dx": abs_summary.get(
                    "longitudinal_mean_dx"
                ),
                "empirical_p_slope": p_for(
                    absolute["tests"],
                    "slope",
                ),
                "empirical_p_mean_dx": p_for(
                    absolute["tests"],
                    "mean_dx",
                ),
                "n_outside_pseudo_95": abs_summary.get(
                    "n_outside_pseudo_95"
                ),
                "n_comparison_points": abs_summary.get(
                    "n_shared_grid_points"
                ),
                "support_min": abs_summary.get(
                    "shared_support_min"
                ),
                "support_max": abs_summary.get(
                    "shared_support_max"
                ),
                "crossover_coordinate": np.nan,
                "crossover_bootstrap_q025": np.nan,
                "crossover_bootstrap_q975": np.nan,
            },
            {
                "branch": "percentile_geometry",
                "longitudinal_slope": pct_summary.get(
                    "longitudinal_slope"
                ),
                "longitudinal_mean_dx": pct_summary.get(
                    "longitudinal_mean_dx"
                ),
                "empirical_p_slope": p_for(
                    percentile["tests"],
                    "slope",
                ),
                "empirical_p_mean_dx": p_for(
                    percentile["tests"],
                    "mean_dx",
                ),
                "n_outside_pseudo_95": pct_summary.get(
                    "n_outside_pseudo_95"
                ),
                "n_comparison_points": pct_summary.get(
                    "n_percentile_bins"
                ),
                "support_min": 0.0,
                "support_max": 1.0,
                "crossover_coordinate": pct_summary.get(
                    "longitudinal_crossover_percentile"
                ),
                "crossover_bootstrap_q025": pct_summary.get(
                    "longitudinal_crossover_bootstrap_q025"
                ),
                "crossover_bootstrap_q975": pct_summary.get(
                    "longitudinal_crossover_bootstrap_q975"
                ),
            },
        ]
    )
    analysis_summary.to_csv(
        outdir / "00_analysis_summary.csv",
        index=False,
    )

    # -----------------------------------------------------------------
    # Input manifest and run config.
    # -----------------------------------------------------------------
    input_manifest = pd.DataFrame(
        [
            {
                "source": "Step9",
                "input": key,
                "path": str(path),
                "exists": bool(path.exists()),
            }
            for key, path in inputs["paths"].items()
            if key.startswith("step9_")
        ]
        +
        [
            {
                "source": "Step7",
                "input": key,
                "path": str(path),
                "exists": bool(path.exists()),
            }
            for key, path in inputs["paths"].items()
            if key.startswith("step7_")
        ]
        +
        [
            {
                "source": "Step9",
                "input": "longitudinal_transitions",
                "path": str(
                    inputs["longitudinal_transitions"]
                ),
                "exists": True,
            },
            {
                "source": "Step7",
                "input": "pseudo_compact_manifest",
                "path": str(
                    inputs["paths"]["step7_manifest"]
                ),
                "exists": True,
            },
        ]
    )
    input_manifest.to_csv(
        outdir / "00_input_manifest.csv",
        index=False,
    )

    run_config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "longitudinal_dir": str(
            inputs["longitudinal_dir"]
        ),
        "pseudo_step7_dir": str(
            inputs["pseudo_step7_dir"]
        ),
        "outdir": str(outdir),
        "dt": int(args.dt),
        "forward_estimand": {
            "AB": (
                "x_observed_rep1_t0 -> dx_observed_rep2 "
                "on forward_ab_eligible"
            ),
            "BA": (
                "x_observed_rep2_t0 -> dx_observed_rep1 "
                "on forward_ba_eligible"
            ),
            "combined": (
                "exact 0.5*AB + 0.5*BA after binning"
            ),
            "AB_BA_combination_policy":
                "exact_equal_fold_weight_after_binning",
            "AB_BA_row_count_weighting": False,
            "operational_TT_filter": False,
        },
        "absolute_abundance": {
            "native_bins_retained": True,
            "stable_pseudo_native_centers_used": True,
            "pseudo_points_require_meets_min_n": True,
            "pseudo_missing_bin_interpolation": False,
            "pseudo_missing_bin_imputation": False,
            "production_requires_complete_shared_grid": True,
            "shared_support_only": True,
            "extrapolation": False,
            "min_pseudo_native_bin_coverage": float(
                args.absolute_min_pseudo_bin_coverage
            ),
            "min_configuration_grid_coverage": float(
                args.absolute_min_grid_coverage
            ),
            "min_configuration_grid_points": int(
                args.absolute_min_grid_points
            ),
        },
        "percentile_geometry": {
            "definition": (
                "transition-level midrank ECDF of initial observed "
                "log-frequency"
            ),
            "longitudinal_unit": (
                "biological subject x fold"
            ),
            "pseudo_unit": (
                "randomization configuration x fold"
            ),
            "equal_percentile_not_equal_absolute_frequency": True,
            "n_bins": int(
                args.percentile_bins
            ),
            "min_n_per_fold_bin": int(
                args.percentile_min_n
            ),
            "min_valid_bins_per_configuration": int(
                args.percentile_min_valid_bins
            ),
            "min_configuration_grid_coverage": float(
                args.percentile_min_grid_coverage
            ),
            "min_configuration_grid_points": int(
                args.percentile_min_grid_points
            ),
        },
        "pseudo_null": {
            "source_biological_subjects": 1,
            "configuration_is_biological_subject": False,
            "n_configurations_discovered": int(
                len(all_ids)
            ),
            "n_configurations_selected": int(
                len(selected_ids)
            ),
            "randomization_interval_not_biological_confidence_interval": True,
            "min_null_configurations_for_inference": int(
                inferential_min_null_configurations
            ),
            "min_null_configurations_for_execution": int(
                execution_min_null_configurations
            ),
            "testing_subset_mode": bool(testing_subset),
            "test_mode_empirical_p_values_reportable": False if testing_subset else True,
        },
        "longitudinal_bootstrap": {
            "unit": "biological subject",
            "n_bootstrap": int(
                n_bootstrap
            ),
            "seed": int(
                bootstrap_seed
            ),
            "same_subject_draw_for_AB_BA": True,
        },
        "empirical_tests": {
            "finite_sample_plus_one_correction": True,
            "slope": "lower",
            "centered_mean_abs_drift": "upper",
            "centered_rms_drift": "upper",
            "mean_abs_drift": "upper",
            "mean_dx": "two_sided_about_pseudo_median",
        },
        "step9_script_version": step9_cfg.get(
            "script_version"
        ),
        "step7_script_version": step7_cfg.get(
            "script_version"
        ),
        "step9_analysis_signature": step9_cfg.get(
            "analysis_signature"
        ),
        "step7_analysis_signature": step7_cfg.get(
            "analysis_signature"
        ),
        "pseudo_percentile_source":
            "Step-7 compact step7_dt1.parquet files",
        "percentile_cache_contract_version":
            PERCENTILE_CACHE_VERSION,
    }

    signature_payload = {
        "script_version": SCRIPT_VERSION,
        "step9_analysis_signature": step9_cfg.get("analysis_signature"),
        "step7_analysis_signature": step7_cfg.get("analysis_signature"),
        "step9_primary": path_identity(inputs["paths"]["step9_primary"]),
        "step7_curves": path_identity(inputs["paths"]["step7_curves"]),
        "step7_manifest": path_identity(inputs["paths"]["step7_manifest"]),
        "selected_configuration_ids": list(selected_ids),
        "parameters": {
            "dt": int(args.dt),
            "absolute_min_pseudo_bin_coverage": float(
                args.absolute_min_pseudo_bin_coverage
            ),
            "absolute_min_grid_coverage": float(
                args.absolute_min_grid_coverage
            ),
            "absolute_min_grid_points": int(
                args.absolute_min_grid_points
            ),
            "percentile_bins": int(args.percentile_bins),
            "percentile_min_n": int(args.percentile_min_n),
            "percentile_min_valid_bins": int(
                args.percentile_min_valid_bins
            ),
            "percentile_min_grid_coverage": float(
                args.percentile_min_grid_coverage
            ),
            "percentile_min_grid_points": int(
                args.percentile_min_grid_points
            ),
            "n_bootstrap": int(n_bootstrap),
            "bootstrap_seed": int(bootstrap_seed),
            "min_null_configurations": int(
                inferential_min_null_configurations
            ),
        },
    }
    analysis_signature = stable_sha256(signature_payload)
    run_config["analysis_signature"] = analysis_signature

    (outdir / "00_run_signature.json").write_text(
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

    (outdir / "00_run_config.json").write_text(
        json.dumps(
            json_safe(run_config),
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    write_readme(outdir)
    write_manifest(outdir)

    # -----------------------------------------------------------------
    # Console summary.
    # -----------------------------------------------------------------
    print("\n[DONE] Step 10 longitudinal versus pseudo forward null")
    print("\n[ABSOLUTE ABUNDANCE]")
    print(
        "shared support:",
        f"{abs_summary.get('shared_support_min'):.6f}",
        "to",
        f"{abs_summary.get('shared_support_max'):.6f}",
    )
    print(
        "longitudinal slope:",
        f"{abs_summary.get('longitudinal_slope'):.6f}",
    )
    print(
        "outside pseudo 95%:",
        f"{int(abs_summary.get('n_outside_pseudo_95'))}/"
        f"{int(abs_summary.get('n_shared_grid_points'))}",
    )
    print(
        "empirical p slope:",
        f"{p_for(absolute['tests'], 'slope'):.6g}",
    )
    print(
        "empirical p mean_dx:",
        f"{p_for(absolute['tests'], 'mean_dx'):.6g}",
    )

    print("\n[PERCENTILE GEOMETRY]")
    print(
        "longitudinal slope:",
        f"{pct_summary.get('longitudinal_slope'):.6f}",
    )
    print(
        "crossover percentile:",
        f"{pct_summary.get('longitudinal_crossover_percentile'):.6f}",
    )
    print(
        "outside pseudo 95%:",
        f"{int(pct_summary.get('n_outside_pseudo_95'))}/"
        f"{int(pct_summary.get('n_percentile_bins'))}",
    )
    print(
        "eligible pseudo configurations:",
        f"{int(pct_summary.get('n_pseudo_eligible'))}",
    )
    print(
        "empirical p slope:",
        f"{p_for(percentile['tests'], 'slope'):.6g}",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
