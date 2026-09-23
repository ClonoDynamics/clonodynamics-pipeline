#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibrate accumulation doses using only the full simulated-frequency oracle.

Purpose and workflow position
-----------------------------
After 01_freeze_empirical_calibration.py, choose q_bio for a dose series before
synthetic counts are generated. This is not measurement-noise calibration and
does not evaluate recovery by the downstream analysis.

Inputs
------
--calibration-dir contains 00_manifest.json, 02_analysis_sampling_design.csv,
11_subject_union_references.csv and the reference files indexed by that table.
References supply f_reference_mean_repertoire. The nominal time grid is retained;
only analysis-eligible visit pairs enter the oracle temporal summaries.
The default sigma_fast is 0.20, seed is 20260918 and --target-Rs is
0,0.25,0.50,1.00. A zero target fixes q_bio=0; it does not suppress fast variation.

Computation
-----------
The log-frequency generator uses a subject reference, independent visit-specific
Gaussian shocks and a Gaussian random walk whose variance increments are
q_bio*delta_t, followed by compositional normalization. The fast shocks and unit
random-walk paths use the same per-subject/per-time seeds as
03_generate_support_conditioned_repertoires.py. Changing q_bio rescales the same
unit path by sqrt(q_bio).
For each eligible interval, oracle displacement variance is calculated across
ALL reference clonotypes (ddof=0). Interval values are averaged within subject
and lag, then equally across the subjects available at each lag. The effect is
R=(V5-V1)/V1. Positive doses use bracketing/bisection (default absolute R tolerance
5e-4, maximum 22 bisection iterations). The best candidate is reported if the
iteration limit is reached; the code does not raise solely for unmet tolerance.

Outputs under --out-dir
-----------------------
00_qbio_oracle_calibration.json: parameters, effect definition and scenarios.
01_scenario_table.csv: target/realized effects, q_bio, oracle slopes and seed.
02_oracle_temporal_profiles.csv: oracle variance at each lag and scenario.
03_search_trace.csv: bracketing and bisection evaluations.
04_manifest.json: script/environment information and output hashes.

Interpretation and safety
-------------------------
Clonotype-specific positivity masks, synthetic counts, common4 selection,
latent midpoint binning, real covariance, real slope and recovered results do
not select q_bio. Eligibility of VISITS is nevertheless inherited from the
empirical design. This full oracle is not the same support-specific estimand as
the recovered covariance on common4 and the complete-case abundance core.
Inputs are read-only and an existing output directory is refused. Memory scales
with reference-clone count and the nominal time grid.
Dependencies: Python >=3.9 with postponed annotations, NumPy and pandas.

Example
-------
python 02_calibrate_oracle_doses.py --calibration-dir CALIBRATION_DIR \
    --out-dir NEW_ORACLE_DIR --sigma-fast 0.20 --seed 20260918

Provenance
----------
Documentation-only edition of 8-calibrate_qbio_oracle_support_conditioned_v5.py.
The computational body and internal versions are unchanged. The payload's
historical generator_match string is retained as a provenance label, not a
runtime import or executable path. See SCRIPT_MAP.csv.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

VERSION = "5.0.0-oracle-qbio-support-conditioned-2026-09-19"
SCHEMA = "clonodynamics_oracle_qbio_support_conditioned_v5"


def require(c: bool, m: str) -> None:
    if not c:
        raise ValueError(m)


def sha(path: Path, bs: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(bs), b""):
            h.update(b)
    return h.hexdigest()


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if math.isfinite(float(x)) else None
    return x


def write_json(path: Path, obj) -> None:
    Path(path).write_text(
        json.dumps(clean(obj), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def stable_log_softmax(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, np.float64)
    m = float(np.max(x))
    return x - (m + math.log(float(np.exp(x - m).sum(dtype=np.float64))))


def parse_bool(s: pd.Series, label: str) -> pd.Series:
    if s.dtype == bool:
        return s.copy()
    mp = {
        "true": True, "false": False,
        "1": True, "0": False,
        "1.0": True, "0.0": False,
    }
    x = s.astype(str).str.strip().str.lower().map(mp)
    require(x.notna().all(), f"Unrecognized boolean in {label}")
    return x.astype(bool)


def parse_targets(text: str) -> List[float]:
    vals = [float(x.strip()) for x in text.split(",") if x.strip()]
    require(vals, "target-Rs is empty")
    require(all(np.isfinite(v) and v >= 0 for v in vals),
            "target-Rs must be finite and >=0")
    out = []
    for v in vals:
        if v not in out:
            out.append(v)
    return out


# These seed functions MUST match generator v5.1 exactly.
def fast_visit_seed(master: int, subject: int, time: int):
    return np.random.SeedSequence(
        [int(master), int(subject), 117031, 501, int(time)]
    )


def rw_increment_seed(master: int, subject: int, time: int):
    return np.random.SeedSequence(
        [int(master), int(subject), 117031, 502, int(time)]
    )


def build_subject_cache(
    cal: Path,
    design_all: pd.DataFrame,
    design_eligible: pd.DataFrame,
    refs: pd.DataFrame,
    seed: int,
) -> Dict[int, dict]:
    times = sorted(int(x) for x in design_all["time"].unique())
    caches: Dict[int, dict] = {}

    for subject in sorted(int(x) for x in refs["subject"].unique()):
        rr = refs.loc[refs["subject"] == subject]
        require(len(rr) == 1, f"Expected one reference row for subject {subject}")
        d = pd.read_csv(
            cal / str(rr.iloc[0]["reference_file"]),
            usecols=["f_reference_mean_repertoire"],
        )
        f = d["f_reference_mean_repertoire"].to_numpy(np.float64)
        require(np.all(np.isfinite(f) & (f > 0)),
                f"Invalid reference frequencies subject {subject}")
        f /= f.sum(dtype=np.float64)
        logf = np.log(f)

        # Frozen common-random-number fields.
        zfast = np.empty((len(times), len(f)), np.float64)
        wunit = np.zeros((len(times), len(f)), np.float64)
        rw = np.zeros(len(f), np.float64)
        prev = times[0]

        for j, t in enumerate(times):
            zfast[j] = np.random.default_rng(
                fast_visit_seed(seed, subject, t)
            ).normal(0.0, 1.0, size=len(f))

            if j > 0:
                dt = float(t - prev)
                require(dt > 0, "times not strictly increasing")
                zr = np.random.default_rng(
                    rw_increment_seed(seed, subject, t)
                ).normal(0.0, 1.0, size=len(f))
                rw = rw + math.sqrt(dt) * zr
                prev = t
            wunit[j] = rw

        observed = set(
            int(x) for x in design_eligible.loc[
                design_eligible["subject"] == subject, "time"
            ]
        )
        require(len(observed) >= 2,
                f"Subject {subject} has <2 analysis-eligible visits")

        caches[subject] = {
            "times": times,
            "logf": logf,
            "zfast": zfast,
            "wunit": wunit,
            "observed": observed,
        }

    return caches


def oracle_profile(caches: Dict[int, dict], sigma_fast: float, q_bio: float):
    require(np.isfinite(sigma_fast) and sigma_fast >= 0, "invalid sigma_fast")
    require(np.isfinite(q_bio) and q_bio >= 0, "invalid q_bio")

    sqrt_q = math.sqrt(float(q_bio))
    rows = []

    for subject, c in caches.items():
        times = c["times"]
        bio = np.empty_like(c["zfast"])
        for j, t in enumerate(times):
            bio[j] = stable_log_softmax(
                c["logf"]
                + float(sigma_fast) * c["zfast"][j]
                + sqrt_q * c["wunit"][j]
            )

        observed = c["observed"]
        for a in range(len(times) - 1):
            for b in range(a + 1, len(times)):
                if times[a] not in observed or times[b] not in observed:
                    continue
                dx = bio[b] - bio[a]
                rows.append({
                    "subject": subject,
                    "time0": int(times[a]),
                    "time1": int(times[b]),
                    "dt": int(times[b] - times[a]),
                    "var_dx_bio": float(np.var(dx, ddof=0)),
                })

    intervals = pd.DataFrame(rows)
    require(not intervals.empty, "No analysis-eligible oracle intervals")

    by_subject_dt = (
        intervals.groupby(["subject", "dt"], as_index=False)["var_dx_bio"]
        .mean()
        .rename(columns={"var_dx_bio": "var_subject_mean"})
    )
    profile = (
        by_subject_dt.groupby("dt", as_index=False)["var_subject_mean"]
        .mean()
        .rename(columns={"var_subject_mean": "var_equal_subject"})
        .sort_values("dt")
        .reset_index(drop=True)
    )
    nsub = (
        by_subject_dt.groupby("dt")["subject"]
        .nunique()
        .rename("n_subjects")
        .reset_index()
    )
    profile = profile.merge(nsub, on="dt", how="left", validate="one_to_one")

    v = dict(zip(profile["dt"].astype(int), profile["var_equal_subject"].astype(float)))
    require(1 in v and 5 in v, "Oracle profile lacks dt=1 or dt=5")
    require(v[1] > 0, "V1 is zero; R is undefined. Use sigma_fast>0.")

    R = float((v[5] - v[1]) / v[1])
    slope = float(np.polyfit(
        profile["dt"].to_numpy(float),
        profile["var_equal_subject"].to_numpy(float),
        1,
    )[0])
    return R, slope, v, profile, by_subject_dt


def calibrate_one(
    caches: Dict[int, dict],
    sigma_fast: float,
    target_R: float,
    q_min: float,
    q_max: float | None,
    iterations: int,
    tolerance: float,
):
    # R=0 is defined as the plateau control q_bio=0; do not tune finite-sample
    # noise to force the realized R exactly to zero.
    if math.isclose(target_R, 0.0, abs_tol=1e-15):
        R, slope, v, prof, bysd = oracle_profile(caches, sigma_fast, 0.0)
        trace = pd.DataFrame([{
            "iteration": 0,
            "target_R": target_R,
            "q_bio": 0.0,
            "R_oracle": R,
            "oracle_slope_per_week": slope,
            "status": "fixed_plateau_q0",
        }])
        return 0.0, R, slope, v, prof, bysd, trace

    q_lo = float(q_min)
    require(q_lo >= 0, "q_min must be >=0")

    if q_max is None:
        # Infinite-clone approximation:
        # Var(dx|dt) ~ 2*sigma_fast^2 + q*dt,
        # target_R ~ 4q / (2*sigma_fast^2 + q).
        denom = max(4.0 - target_R, 1e-9)
        approx = (2.0 * target_R * sigma_fast**2) / denom
        q_hi = max(4.0 * approx, 0.01)
    else:
        q_hi = float(q_max)
    require(q_hi > q_lo, "q_max must exceed q_min")

    trace_rows = []
    Rlo, slo, _, _, _ = oracle_profile(caches, sigma_fast, q_lo)
    Rhi, shi, _, _, _ = oracle_profile(caches, sigma_fast, q_hi)
    trace_rows.extend([
        {"iteration": -2, "target_R": target_R, "q_bio": q_lo,
         "R_oracle": Rlo, "oracle_slope_per_week": slo, "status": "lower_bracket"},
        {"iteration": -1, "target_R": target_R, "q_bio": q_hi,
         "R_oracle": Rhi, "oracle_slope_per_week": shi, "status": "upper_bracket"},
    ])

    # Expand upper bracket if needed.
    expand = 0
    while Rhi < target_R and expand < 12:
        q_hi *= 2.0
        Rhi, shi, _, _, _ = oracle_profile(caches, sigma_fast, q_hi)
        trace_rows.append({
            "iteration": -10 - expand,
            "target_R": target_R,
            "q_bio": q_hi,
            "R_oracle": Rhi,
            "oracle_slope_per_week": shi,
            "status": "expanded_upper_bracket",
        })
        expand += 1
    require(Rhi >= target_R,
            f"Could not bracket target R={target_R}; R(q_hi)={Rhi}")

    best = None
    for it in range(iterations):
        q = 0.5 * (q_lo + q_hi)
        R, slope, _, _, _ = oracle_profile(caches, sigma_fast, q)
        err = abs(R - target_R)
        trace_rows.append({
            "iteration": it,
            "target_R": target_R,
            "q_bio": q,
            "R_oracle": R,
            "oracle_slope_per_week": slope,
            "abs_R_error": err,
            "status": "bisection",
        })
        if best is None or err < best[0]:
            best = (err, q, R, slope)
        if err <= tolerance:
            break
        if R < target_R:
            q_lo = q
        else:
            q_hi = q

    require(best is not None, "Bisection produced no candidate")
    _, q_final, _, _ = best
    R, slope, v, prof, bysd = oracle_profile(caches, sigma_fast, q_final)
    trace = pd.DataFrame(trace_rows)
    return q_final, R, slope, v, prof, bysd, trace


def scenario_label(target_R: float) -> str:
    return "R" + f"{target_R:.2f}".replace(".", "p")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calibration-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--sigma-fast", type=float, default=0.20)
    ap.add_argument("--target-Rs", default="0,0.25,0.50,1.00")
    ap.add_argument("--seed", type=int, default=20260918)
    ap.add_argument("--q-min", type=float, default=0.0)
    ap.add_argument("--q-max", type=float, default=None)
    ap.add_argument("--iterations", type=int, default=22)
    ap.add_argument("--tolerance", type=float, default=5e-4)
    args = ap.parse_args()

    require(np.isfinite(args.sigma_fast) and args.sigma_fast > 0,
            "sigma_fast must be >0 because R uses V1 in the denominator")
    require(args.iterations >= 1, "iterations must be >=1")
    require(args.tolerance > 0, "tolerance must be >0")

    cal = args.calibration_dir.resolve()
    out = args.out_dir.resolve()
    if out.exists():
        raise FileExistsError(f"Refusing overwrite: {out}")
    out.mkdir(parents=True)

    design = pd.read_csv(cal / "02_analysis_sampling_design.csv")
    refs = pd.read_csv(cal / "11_subject_union_references.csv")
    require("analysis_eligible_pair" in design.columns,
            "02_analysis_sampling_design.csv lacks analysis_eligible_pair")
    design["analysis_eligible_pair"] = parse_bool(
        design["analysis_eligible_pair"], "analysis_eligible_pair"
    )
    eligible = design.loc[design["analysis_eligible_pair"]].copy()

    targets = parse_targets(args.target_Rs)
    caches = build_subject_cache(cal, design, eligible, refs, args.seed)

    scenario_rows = []
    profile_rows = []
    trace_parts = []

    for target_R in targets:
        q, R, slope, v, profile, bysd, trace = calibrate_one(
            caches=caches,
            sigma_fast=float(args.sigma_fast),
            target_R=float(target_R),
            q_min=float(args.q_min),
            q_max=args.q_max,
            iterations=int(args.iterations),
            tolerance=float(args.tolerance),
        )
        label = scenario_label(target_R)
        scenario_rows.append({
            "scenario_label": label,
            "target_R": float(target_R),
            "sigma_fast": float(args.sigma_fast),
            "q_bio": float(q),
            "realized_R": float(R),
            "oracle_slope_per_week": float(slope),
            "oracle_V1": float(v[1]),
            "oracle_V5": float(v[5]),
            "abs_R_error": float(abs(R - target_R)),
            "seed": int(args.seed),
        })
        pp = profile.copy()
        pp.insert(0, "scenario_label", label)
        pp.insert(1, "target_R", float(target_R))
        pp.insert(2, "q_bio", float(q))
        profile_rows.append(pp)

        tt = trace.copy()
        tt.insert(0, "scenario_label", label)
        trace_parts.append(tt)

        print(
            f"[ORACLE] {label}: target_R={target_R:.3f} "
            f"q_bio={q:.10f} realized_R={R:.6f} "
            f"slope={slope:.8f} V1={v[1]:.6f} V5={v[5]:.6f}",
            flush=True,
        )

    scen = pd.DataFrame(scenario_rows)
    profiles = pd.concat(profile_rows, ignore_index=True)
    traces = pd.concat(trace_parts, ignore_index=True)

    scen.to_csv(out / "01_scenario_table.csv", index=False)
    profiles.to_csv(out / "02_oracle_temporal_profiles.csv", index=False)
    traces.to_csv(out / "03_search_trace.csv", index=False)

    payload = {
        "schema_version": SCHEMA,
        "script_version": VERSION,
        "calibration_dir": str(cal),
        "calibration_manifest_sha256": sha(cal / "00_manifest.json"),
        "seed": int(args.seed),
        "sigma_fast": float(args.sigma_fast),
        "target_Rs": targets,
        "effect_definition": "R=(V5-V1)/V1 using equal-subject full biological oracle",
        "oracle_weighting": "interval mean within subject and lag, then equal-subject mean by lag",
        "common_random_numbers": {
            "fast_visit_shocks": True,
            "unit_random_walk_path": True,
            "generator_match": "7b-generate_empirical_support_conditioned_v5_1.py",
        },
        "anti_leakage": {
            "empirical_support_used_to_select_q_bio": False,
            "synthetic_counts_used_to_select_q_bio": False,
            "real_cross_cov_used": False,
            "real_temporal_slope_used": False,
            "recovered_pipeline_output_used": False,
        },
        "scenarios": scenario_rows,
        "note": (
            "target_R=0 is the plateau control q_bio=0 and is not numerically tuned; "
            "its finite-sample realized R may differ slightly from zero."
        ),
    }
    write_json(out / "00_qbio_oracle_calibration.json", payload)

    write_json(out / "04_manifest.json", {
        "schema_version": SCHEMA,
        "script_version": VERSION,
        "script_sha256": sha(Path(__file__).resolve()),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "outputs": [
            {"path": p.name, "sha256": sha(p)}
            for p in sorted(out.iterdir())
            if p.is_file() and p.name != "04_manifest.json"
        ],
    })

    print("\nOracle-only q_bio calibration completed.")
    print(f"  sigma_fast={args.sigma_fast:.6f}")
    for r in scenario_rows:
        print(
            f"  {r['scenario_label']}: target_R={r['target_R']:.2f}; "
            f"q_bio={r['q_bio']:.10f}; realized_R={r['realized_R']:.6f}; "
            f"slope={r['oracle_slope_per_week']:.8f}"
        )
    print(f"  wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
