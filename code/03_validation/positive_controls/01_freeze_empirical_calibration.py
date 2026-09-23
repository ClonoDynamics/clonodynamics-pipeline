#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Freeze the empirical design and reference repertoires for synthetic controls.

Purpose and workflow position
-----------------------------
Prepare the input bundle used by 02_calibrate_oracle_doses.py and
03_generate_support_conditioned_repertoires.py. This is calibration of the
empirical reference/input structure, NOT oracle dose calibration or simulation.
It separates physically observed replicate pairs from pairs eligible after
state inference, and records the real abundance grid and temporal-analysis core.

Inputs
------
--step1-characteristics: repertoire_characteristics.csv.
--step2-qc: latent_pair_qc.csv, including success/included_in_downstream flags.
--step2-params: latent_params_by_pair.csv, including k, gamma, fmin and grid_size.
--source-dir: corrected repertoire files containing aaSeqCDR3 and readCount
(or --count-column). Supported formats include text, Parquet and Feather.
--transition-report: the real Step-5 transition_assembly_report.md.
--step11-dir: 00_run_config.json, 01_fluctuation_bin_edges.csv and
06_support_by_dt.csv from the real longitudinal analysis.
--step13-dir: 00_run_config.json, 01_core_bin_selection.csv,
02_complete_case_subjects.csv and 05_temporal_slope_summary.csv.

Main operations
---------------
Repeated aaSeqCDR3 entries are summed WITHIN each repertoire. Positive counts
are checked against Step-1 depth/richness. Each subject's reference includes
the union of all raw observed clonotypes, including raw files from a later
excluded pair. f_reference_mean_repertoire is the equal-repertoire mean of
within-repertoire frequencies, with zero contributions when a clone is absent;
f_reference_pooled_all_counts is also exported as a distinct reference.
Step-2 k is copied to kappa; it is not re-estimated here. The observed support,
Step-11 grid and Step-13 core/cohort are recorded without tuning a simulated
signal to the real covariance or temporal slope.

Outputs under --out-dir
-----------------------
00_manifest.json
01_raw_sampling_design.csv
02_analysis_sampling_design.csv
03_step2_pair_parameters.csv
04_observation_structure_targets.json
05_cohort_and_analysis_summary.json
06_step11_bin_edges.csv; 06_step11_run_config.json
07_step11_support_by_dt.csv
08_step13_core_selection.csv; 08_step13_run_config.json
09_step13_complete_case_subjects.csv
10_step13_temporal_slope_summary.csv
11_subject_union_references.csv
12_source_repertoire_audit.csv
13_validation_contract.json
subject_union_repertoires/subject_<subject>.csv.gz

Safety, scope and execution
---------------------------
Inputs are read-only. An existing output directory is refused; a new bundle is
built in a staging directory before being moved into place. Default expected
counts describe the current empirical dataset and can be changed explicitly.
Dependencies: Python >=3.9, NumPy and pandas; a compatible columnar-file engine
is additionally needed for Parquet/Feather input. No simulation or model fit is
performed. Real temporal estimates copied into the bundle are documentary only.
Run `python 01_freeze_empirical_calibration.py --help` for the full CLI.

Provenance
----------
Documentation-only edition of 0-build_validation_calibration_v5.py. The supplied
original was retrieved from the file Library and its SHA-256 matched the
calibration manifest. Computational code, CLI, output names, schema and internal
version constants are unchanged; see SCRIPT_MAP.csv for original/current hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

SCRIPT_VERSION = "5.0.0-corrected-data-freeze-2026-09-18"
SCHEMA_VERSION = "clonodynamics_validation_calibration_v5"

DEFAULT_TIMES = [1, 2, 3, 4, 5, 6]
SUPPORTED_SUFFIXES = (
    ".tsv", ".tsv.gz", ".csv", ".csv.gz", ".txt", ".txt.gz",
    ".parquet", ".pq", ".feather", ".arrow",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(block_size), b""):
            h.update(block)
    return h.hexdigest()


def json_clean(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): json_clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [json_clean(v) for v in x]
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if math.isfinite(float(x)) else None
    return x


def write_json(path: Path, obj: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(json_clean(obj), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def parse_bool_series(series: pd.Series, label: str) -> pd.Series:
    if series.dtype == bool:
        return series.copy()
    mapping = {
        "true": True, "1": True, "1.0": True, "yes": True,
        "false": False, "0": False, "0.0": False, "no": False,
    }
    vals = series.astype(str).str.strip().str.lower()
    require(vals.isin(mapping).all(), f"{label}: unrecognized boolean values.")
    return vals.map(mapping).astype(bool)


def parse_time_values(text: str) -> List[int]:
    vals = [int(x.strip()) for x in text.split(",") if x.strip()]
    require(vals == sorted(set(vals)) and len(vals) >= 2,
            "--time-values must be an increasing duplicate-free list.")
    return vals


def read_repertoire(path: Path, count_column: str) -> pd.DataFrame:
    name = path.name.lower()
    cols = ["aaSeqCDR3", count_column]
    if name.endswith((".parquet", ".pq")):
        df = pd.read_parquet(path, columns=cols)
    elif name.endswith((".feather", ".arrow")):
        df = pd.read_feather(path, columns=cols)
    else:
        sep = "\t" if name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
        df = pd.read_csv(path, sep=sep, usecols=cols, keep_default_na=False,
                         dtype={"aaSeqCDR3": str})

    require(not df.empty, f"Empty repertoire: {path}")
    require(df["aaSeqCDR3"].notna().all(), f"Missing aaSeqCDR3 in {path}")
    ids = df["aaSeqCDR3"].astype(str)
    require(((ids.str.len() > 0) & (ids == ids.str.strip())).all(),
            f"Blank or whitespace-padded aaSeqCDR3 in {path}")

    counts = pd.to_numeric(df[count_column], errors="raise").to_numpy(dtype=float)
    require(np.all(np.isfinite(counts)), f"Non-finite {count_column} in {path}")
    require(np.all((counts >= 0) & (counts < 2**53) & (counts == np.floor(counts))),
            f"{count_column} must contain exact non-negative integers in {path}")

    x = pd.DataFrame({
        "aaSeqCDR3": ids.to_numpy(),
        "readCount": counts.astype(np.int64),
    })
    # This matches the analysis convention: aggregate repeated AA CDR3 rows within
    # each repertoire only, never across repertoire files.
    x = x.groupby("aaSeqCDR3", as_index=False, sort=False)["readCount"].sum()
    x = x.loc[x["readCount"] > 0].copy()
    require(not x.empty, f"No positive clonotypes after aggregation: {path}")
    return x


def strip_data_suffix(name: str) -> str:
    out = Path(name).name
    # Repeated stripping handles .tsv.gz etc.
    changed = True
    while changed:
        changed = False
        low = out.lower()
        for suf in (".gz", ".tsv", ".csv", ".txt", ".parquet", ".pq", ".feather", ".arrow"):
            if low.endswith(suf):
                out = out[:-len(suf)]
                changed = True
                break
    return out


class SourceResolver:
    def __init__(self, root: Path):
        self.root = root.resolve()
        require(self.root.is_dir(), f"Source repertoire directory not found: {self.root}")
        self.by_name: Dict[str, List[Path]] = {}
        self.by_stem: Dict[str, List[Path]] = {}
        for p in self.root.rglob("*"):
            if not p.is_file():
                continue
            low = p.name.lower()
            if not any(low.endswith(s) for s in SUPPORTED_SUFFIXES):
                continue
            self.by_name.setdefault(p.name, []).append(p.resolve())
            self.by_stem.setdefault(strip_data_suffix(p.name), []).append(p.resolve())

    def resolve(self, requested: str) -> Path:
        base = Path(str(requested)).name
        direct = self.by_name.get(base, [])
        if len(direct) == 1:
            return direct[0]
        if len(direct) > 1:
            raise RuntimeError(f"Ambiguous exact filename {base}: {direct}")
        stem = strip_data_suffix(base)
        hits = self.by_stem.get(stem, [])
        if len(hits) == 1:
            return hits[0]
        if not hits:
            raise FileNotFoundError(
                f"Could not locate source repertoire {requested!r} below {self.root}"
            )
        raise RuntimeError(
            f"Ambiguous equivalent-format repertoire for {requested!r}: {hits}"
        )


def parse_transition_report(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")

    def scalar(name: str, cast):
        m = re.search(rf"(?m)^-\s*{re.escape(name)}:\s*(.+?)\s*$", text)
        if not m:
            return None
        raw = m.group(1).strip()
        try:
            return cast(raw)
        except Exception:
            return raw

    def json_value(name: str):
        m = re.search(rf"(?m)^-\s*{re.escape(name)}:\s*(\{{.*\}}|\[.*\])\s*$", text)
        if not m:
            return None
        try:
            return json.loads(m.group(1))
        except json.JSONDecodeError:
            return None

    return {
        "script": scalar("script", str),
        "script_version": scalar("script_version", str),
        "prepared_state_rows": scalar("prepared_state_rows", int),
        "min_dt": scalar("min_dt", int),
        "max_dt": scalar("max_dt", int),
        "transition_rows": scalar("transition_rows", int),
        "n_subjects": scalar("n_subjects", int),
        "n_subject_clone_pairs": scalar("n_subject_clone_pairs", int),
        "n_subject_interval_pairs": scalar("n_subject_interval_pairs", int),
        "fraction_common4": scalar("fraction_common4", float),
        "fraction_forward_ab_eligible": scalar("fraction_forward_ab_eligible", float),
        "fraction_forward_ba_eligible": scalar("fraction_forward_ba_eligible", float),
        "counts_by_dt": json_value("counts_by_dt"),
        "common4_fraction_by_dt": json_value("common4_fraction_by_dt"),
    }


def validate_step11_contract(config: Dict[str, Any]) -> None:
    require(config.get("primary_support") == "common4",
            "Step 11 contract mismatch: primary_support != common4")
    require(config.get("primary_conditioning") == "xmid_latent",
            "Step 11 contract mismatch: primary_conditioning != xmid_latent")
    est = config.get("primary_estimands", {})
    text = json.dumps(est).lower()
    require("cross_cov" in text and "dx_observed_rep1" in text and "dx_observed_rep2" in text,
            "Step 11 contract mismatch: cross-replicate observed displacement covariance not found.")


def validate_step13_contract(config: Dict[str, Any]) -> None:
    require(config.get("primary_metric") == "cross_cov",
            "Step 13 contract mismatch: primary_metric != cross_cov")
    require(config.get("primary_estimator") == "equal_bin_equal_subject",
            "Step 13 contract mismatch: primary_estimator != equal_bin_equal_subject")
    require(config.get("complete_case", {}).get("requires_all_selected_bins_at_all_selected_dt") is True,
            "Step 13 contract mismatch: expected fixed complete-case geometry.")


def build_design(step1: pd.DataFrame, qc: pd.DataFrame, times: List[int]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    req1 = {"subject", "time", "replica", "file", "depth", "n_clonotypes"}
    req_qc = {
        "subject", "time", "pair_id", "success", "included_in_downstream",
        "message", "depth_rep1", "depth_rep2", "n_clones_rep1", "n_clones_rep2",
        "n_clones_union", "n_clones_intersection",
    }
    require(req1.issubset(step1.columns),
            f"Step 1 missing columns: {sorted(req1-set(step1.columns))}")
    require(req_qc.issubset(qc.columns),
            f"Step 2 QC missing columns: {sorted(req_qc-set(qc.columns))}")

    for c in ["subject", "time", "replica", "depth", "n_clonotypes"]:
        step1[c] = pd.to_numeric(step1[c], errors="raise").astype(int)
    qc["subject"] = pd.to_numeric(qc["subject"], errors="raise").astype(int)
    qc["time"] = pd.to_numeric(qc["time"], errors="raise").astype(int)
    qc["success"] = parse_bool_series(qc["success"], "Step2.success")
    qc["included_in_downstream"] = parse_bool_series(
        qc["included_in_downstream"], "Step2.included_in_downstream"
    )

    require(not step1.duplicated(["subject", "time", "replica"]).any(),
            "Step 1 has duplicate subject/time/replica keys.")
    require(set(step1["replica"].unique()) <= {1, 2},
            "Step 1 contains replica labels other than 1/2.")
    require(not qc.duplicated(["subject", "time"]).any(),
            "Step 2 QC has duplicate subject/time keys.")

    subjects = sorted(int(x) for x in step1["subject"].unique())
    idx = step1.set_index(["subject", "time", "replica"])
    qidx = qc.set_index(["subject", "time"])

    raw_rows: List[dict] = []
    analysis_rows: List[dict] = []

    for s in subjects:
        for t in times:
            row = {"subject": s, "time": t}
            observed = []
            for r in (1, 2):
                key = (s, t, r)
                exists = key in idx.index
                observed.append(exists)
                row[f"observed_rep{r}"] = bool(exists)
                row[f"source_file_rep{r}"] = str(idx.loc[key, "file"]) if exists else ""
                row[f"depth_rep{r}"] = int(idx.loc[key, "depth"]) if exists else np.nan
                row[f"n_clonotypes_rep{r}"] = int(idx.loc[key, "n_clonotypes"]) if exists else np.nan

            n_obs = sum(observed)
            row["raw_pair_complete"] = n_obs == 2
            row["raw_sampling_status"] = {
                0: "not_observed",
                1: "incomplete_pair",
                2: "complete_pair",
            }[n_obs]
            row["pair_id"] = f"{s}_{t}"
            raw_rows.append(dict(row))

            a = dict(row)
            if (s, t) in qidx.index:
                q = qidx.loc[(s, t)]
                if isinstance(q, pd.DataFrame):
                    raise ValueError(f"Duplicate Step 2 QC key {(s,t)}")
                a["step2_success"] = bool(q["success"])
                a["included_in_downstream"] = bool(q["included_in_downstream"])
                a["step2_message"] = str(q["message"])
                a["n_clones_union_step2"] = int(q["n_clones_union"])
                a["n_clones_intersection_step2"] = int(q["n_clones_intersection"])
                if n_obs == 2:
                    require(int(q["depth_rep1"]) == int(a["depth_rep1"]),
                            f"Depth mismatch Step1/Step2 for {(s,t,1)}")
                    require(int(q["depth_rep2"]) == int(a["depth_rep2"]),
                            f"Depth mismatch Step1/Step2 for {(s,t,2)}")
                    require(int(q["n_clones_rep1"]) == int(a["n_clonotypes_rep1"]),
                            f"Richness mismatch Step1/Step2 for {(s,t,1)}")
                    require(int(q["n_clones_rep2"]) == int(a["n_clonotypes_rep2"]),
                            f"Richness mismatch Step1/Step2 for {(s,t,2)}")
            else:
                a["step2_success"] = False
                a["included_in_downstream"] = False
                a["step2_message"] = "NO_STEP2_PAIR"
                a["n_clones_union_step2"] = np.nan
                a["n_clones_intersection_step2"] = np.nan

            a["analysis_eligible_pair"] = bool(
                a["raw_pair_complete"]
                and a["step2_success"]
                and a["included_in_downstream"]
            )
            if not a["raw_pair_complete"]:
                a["analysis_exclusion_reason"] = "raw_pair_not_complete"
            elif not a["step2_success"]:
                a["analysis_exclusion_reason"] = "step2_fit_failure"
            elif not a["included_in_downstream"]:
                a["analysis_exclusion_reason"] = "step2_not_included_downstream"
            else:
                a["analysis_exclusion_reason"] = ""
            analysis_rows.append(a)

    raw = pd.DataFrame(raw_rows)
    analysis = pd.DataFrame(analysis_rows)
    return raw, analysis


def build_subject_union_templates(
    step1: pd.DataFrame,
    raw_design: pd.DataFrame,
    analysis_design: pd.DataFrame,
    source_dir: Path,
    count_column: str,
    out: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    resolver = SourceResolver(source_dir)
    refdir = out / "subject_union_repertoires"
    refdir.mkdir()

    subjects = sorted(int(x) for x in step1["subject"].unique())
    summary_rows: List[dict] = []
    audit_rows: List[dict] = []

    for s in subjects:
        s1 = step1.loc[step1["subject"] == s].sort_values(["time", "replica"]).copy()
        require(not s1.empty, f"No Step 1 repertoires for subject {s}")

        sum_frac: Optional[pd.Series] = None
        sum_count: Optional[pd.Series] = None
        n_rep_detected: Optional[pd.Series] = None
        n_tp_raw: Optional[pd.Series] = None
        n_tp_analysis: Optional[pd.Series] = None

        total_depth = 0
        max_richness = 0
        raw_files = 0
        positive_entries = 0
        odd_entries = 0
        singleton_entries = 0
        gcd_all = 0

        # Cache per timepoint clone sets to count union-positive timepoint states
        # only once per clone/time, irrespective of A/B.
        raw_tp_sets: Dict[int, List[pd.Index]] = {}
        analysis_tp_sets: Dict[int, List[pd.Index]] = {}

        for row in s1.itertuples(index=False):
            path = resolver.resolve(str(row.file))
            d = read_repertoire(path, count_column=count_column)
            ser = d.set_index("aaSeqCDR3")["readCount"].astype(np.int64)

            depth = int(ser.to_numpy(dtype=np.int64).sum(dtype=np.int64))
            richness = int(len(ser))
            require(depth == int(row.depth),
                    f"Source depth mismatch for {row.file}: {depth} != Step1 {int(row.depth)}")
            require(richness == int(row.n_clonotypes),
                    f"Source richness mismatch for {row.file}: {richness} != Step1 {int(row.n_clonotypes)}")

            frac = ser.astype(np.float64) / float(depth)
            sum_frac = frac.copy() if sum_frac is None else sum_frac.add(frac, fill_value=0.0)
            sum_count = ser.astype(np.float64).copy() if sum_count is None else sum_count.add(
                ser.astype(np.float64), fill_value=0.0
            )
            det = pd.Series(np.ones(len(ser), dtype=np.float64), index=ser.index)
            n_rep_detected = det.copy() if n_rep_detected is None else n_rep_detected.add(
                det, fill_value=0.0
            )

            t = int(row.time)
            raw_tp_sets.setdefault(t, []).append(ser.index)
            elig = analysis_design.loc[
                (analysis_design["subject"] == s)
                & (analysis_design["time"] == t),
                "analysis_eligible_pair",
            ]
            is_analysis_time = bool(elig.iloc[0]) if len(elig) == 1 else False
            if is_analysis_time:
                analysis_tp_sets.setdefault(t, []).append(ser.index)

            vals = ser.to_numpy(np.int64)
            positive_entries += len(vals)
            odd_entries += int(np.count_nonzero(vals % 2))
            singleton_entries += int(np.count_nonzero(vals == 1))
            this_gcd = int(np.gcd.reduce(vals))
            gcd_all = this_gcd if gcd_all == 0 else math.gcd(gcd_all, this_gcd)
            total_depth += depth
            max_richness = max(max_richness, richness)
            raw_files += 1

            audit_rows.append({
                "subject": s,
                "time": t,
                "replica": int(row.replica),
                "requested_file": str(row.file),
                "resolved_file": str(path),
                "source_sha256": sha256_file(path),
                "depth": depth,
                "n_observed_clonotypes": richness,
                "positive_count_entries": int(len(vals)),
                "odd_positive_count_entries": int(np.count_nonzero(vals % 2)),
                "singleton_entries": int(np.count_nonzero(vals == 1)),
                "gcd_positive_counts": this_gcd,
                "analysis_eligible_timepoint": is_analysis_time,
            })

        assert sum_frac is not None and sum_count is not None and n_rep_detected is not None

        def accumulate_timepoint_sets(tp_sets: Dict[int, List[pd.Index]]) -> pd.Series:
            acc: Optional[pd.Series] = None
            for _t, idxs in sorted(tp_sets.items()):
                u = idxs[0]
                for ix in idxs[1:]:
                    u = u.union(ix, sort=False)
                one = pd.Series(np.ones(len(u), dtype=np.float64), index=u)
                acc = one.copy() if acc is None else acc.add(one, fill_value=0.0)
            if acc is None:
                return pd.Series(dtype=np.float64)
            return acc

        n_tp_raw = accumulate_timepoint_sets(raw_tp_sets)
        n_tp_analysis = accumulate_timepoint_sets(analysis_tp_sets)

        union_index = sum_frac.index
        n_tp_raw = n_tp_raw.reindex(union_index, fill_value=0.0)
        n_tp_analysis = n_tp_analysis.reindex(union_index, fill_value=0.0)

        ref = pd.DataFrame({
            "aaSeqCDR3": union_index.astype(str),
            "total_readCount_all_raw_repertoires": np.rint(
                sum_count.reindex(union_index).to_numpy(np.float64)
            ).astype(np.int64),
            "n_repertoires_detected_raw": np.rint(
                n_rep_detected.reindex(union_index).to_numpy(np.float64)
            ).astype(np.int16),
            "n_timepoints_observed_raw": np.rint(
                n_tp_raw.to_numpy(np.float64)
            ).astype(np.int8),
            "n_timepoints_observed_analysis": np.rint(
                n_tp_analysis.to_numpy(np.float64)
            ).astype(np.int8),
            "f_reference_mean_repertoire": (
                sum_frac / float(raw_files)
            ).reindex(union_index).to_numpy(np.float64),
            "f_reference_pooled_all_counts": (
                sum_count / float(total_depth)
            ).reindex(union_index).to_numpy(np.float64),
        })

        for c in ("f_reference_mean_repertoire", "f_reference_pooled_all_counts"):
            vals = ref[c].to_numpy(np.float64)
            require(np.all(np.isfinite(vals) & (vals > 0)),
                    f"Invalid reference probabilities for subject {s}, {c}")
            ref[c] = vals / vals.sum(dtype=np.float64)
            require(np.isclose(ref[c].sum(), 1.0, rtol=0, atol=1e-10),
                    f"Reference normalization failed for subject {s}, {c}")

        outfile = refdir / f"subject_{s}.csv.gz"
        ref.sort_values("aaSeqCDR3").to_csv(
            outfile, index=False, compression={"method": "gzip", "mtime": 0}
        )

        n_union = len(ref)
        raw_tp_positive = int(ref["n_timepoints_observed_raw"].to_numpy(dtype=np.int64).sum(dtype=np.int64))
        analysis_tp_positive = int(
            ref["n_timepoints_observed_analysis"].to_numpy(dtype=np.int64).sum(dtype=np.int64)
        )
        summary_rows.append({
            "subject": s,
            "n_raw_source_repertoires": raw_files,
            "n_raw_observed_timepoints": int(
                raw_design.loc[
                    (raw_design["subject"] == s) & raw_design["raw_pair_complete"],
                    "time",
                ].nunique()
            ),
            "n_analysis_eligible_timepoints": int(
                analysis_design.loc[
                    (analysis_design["subject"] == s)
                    & analysis_design["analysis_eligible_pair"],
                    "time",
                ].nunique()
            ),
            "n_subject_union_clonotypes": n_union,
            "max_single_repertoire_richness": max_richness,
            "union_to_max_richness_ratio": n_union / max_richness,
            "sum_union_positive_states_raw": raw_tp_positive,
            "sum_union_positive_states_analysis_reference_universe": analysis_tp_positive,
            "total_raw_depth": total_depth,
            "positive_count_entries": positive_entries,
            "odd_positive_count_entries": odd_entries,
            "singleton_entries": singleton_entries,
            "gcd_all_positive_counts": gcd_all,
            "reference_file": str(outfile.relative_to(out)),
            "reference_sha256": sha256_file(outfile),
        })
        print(
            f"[REFERENCE] subject={s} union={n_union:,} "
            f"raw_repertoires={raw_files} analysis_timepoints="
            f"{summary_rows[-1]['n_analysis_eligible_timepoints']}",
            flush=True,
        )

    return pd.DataFrame(summary_rows), pd.DataFrame(audit_rows)


def build(args: argparse.Namespace, out: Path) -> Dict[str, Any]:
    times = parse_time_values(args.time_values)

    step1_path = args.step1_characteristics.resolve(strict=True)
    qc_path = args.step2_qc.resolve(strict=True)
    params_path = args.step2_params.resolve(strict=True)
    report_path = args.transition_report.resolve(strict=True)
    step11 = args.step11_dir.resolve(strict=True)
    step13 = args.step13_dir.resolve(strict=True)
    source_dir = args.source_dir.resolve(strict=True)

    step1 = pd.read_csv(step1_path)
    qc = pd.read_csv(qc_path)
    params = pd.read_csv(params_path)

    require(not params.duplicated(["subject", "time"]).any(),
            "Step 2 params contain duplicate subject/time keys.")
    require({"subject", "time", "success", "k", "gamma", "fmin", "grid_size"}.issubset(params.columns),
            "Step 2 params missing required calibration columns.")

    raw_design, analysis_design = build_design(step1, qc, times)

    # Current corrected-data expectations are defaults only; they can be changed
    # explicitly when the dataset changes again.
    n_subjects = int(step1["subject"].nunique())
    n_repertoires = int(len(step1))
    n_raw_pairs = int(raw_design["raw_pair_complete"].sum())
    n_analysis_pairs = int(analysis_design["analysis_eligible_pair"].sum())

    require(n_subjects == args.expected_subjects,
            f"Found {n_subjects} subjects; expected {args.expected_subjects}")
    require(n_repertoires == args.expected_repertoires,
            f"Found {n_repertoires} Step-1 repertoires; expected {args.expected_repertoires}")
    require(n_raw_pairs == args.expected_raw_pairs,
            f"Found {n_raw_pairs} raw complete pairs; expected {args.expected_raw_pairs}")
    require(n_analysis_pairs == args.expected_analysis_pairs,
            f"Found {n_analysis_pairs} analysis-eligible pairs; expected {args.expected_analysis_pairs}")

    # Scheduled interval geometry from analysis-eligible visits.
    interval_pairs = 0
    for s, g in analysis_design.groupby("subject"):
        n = int(g["analysis_eligible_pair"].sum())
        interval_pairs += math.comb(n, 2)
    require(interval_pairs == args.expected_interval_pairs,
            f"Analysis design implies {interval_pairs} within-subject interval pairs; "
            f"expected {args.expected_interval_pairs}")

    raw_design.to_csv(out / "01_raw_sampling_design.csv", index=False)
    analysis_design.to_csv(out / "02_analysis_sampling_design.csv", index=False)

    # Freeze Step-2 parameter table as-is, but add explicit downstream eligibility.
    p = params.copy()
    p["subject"] = pd.to_numeric(p["subject"], errors="raise").astype(int)
    p["time"] = pd.to_numeric(p["time"], errors="raise").astype(int)
    p["success"] = parse_bool_series(p["success"], "Step2 params.success")
    p = p.merge(
        analysis_design[["subject", "time", "analysis_eligible_pair", "analysis_exclusion_reason"]],
        on=["subject", "time"], how="left", validate="one_to_one",
    )
    p["kappa"] = pd.to_numeric(p["k"], errors="coerce")
    p["kappa_source"] = "Step2 latent_params_by_pair:k"
    p["calibration_status"] = "observation_model_calibration_parameter_not_ground_truth"
    p.to_csv(out / "03_step2_pair_parameters.csv", index=False)

    # Observation targets from eligible Step-2 cells.
    qc_eligible = qc.loc[parse_bool_series(qc["included_in_downstream"], "included_in_downstream")].copy()
    positive_rep_states = int(
        (pd.to_numeric(qc_eligible["n_clones_rep1"]) +
         pd.to_numeric(qc_eligible["n_clones_rep2"])).sum()
    )
    union_states = int(pd.to_numeric(qc_eligible["n_clones_union"]).sum())
    both_states = int(pd.to_numeric(qc_eligible["n_clones_intersection"]).sum())

    transition_report = parse_transition_report(report_path)
    require(transition_report.get("transition_rows") is not None,
            "Could not parse transition_rows from transition report.")
    require(transition_report.get("n_subject_interval_pairs") is not None,
            "Could not parse n_subject_interval_pairs from transition report.")
    require(int(transition_report["n_subject_interval_pairs"]) == interval_pairs,
            "Transition report interval-pair count disagrees with frozen analysis design.")
    if transition_report.get("prepared_state_rows") is not None:
        require(int(transition_report["prepared_state_rows"]) == union_states,
                "Step-5 prepared_state_rows disagrees with Step-2 union-positive target.")

    # Final Step 11.
    s11_cfg_path = step11 / "00_run_config.json"
    s11_edges_path = step11 / "01_fluctuation_bin_edges.csv"
    s11_support_path = step11 / "06_support_by_dt.csv"
    for f in (s11_cfg_path, s11_edges_path, s11_support_path):
        require(f.is_file(), f"Required Step-11 file not found: {f}")
    s11_cfg = json.loads(s11_cfg_path.read_text(encoding="utf-8"))
    validate_step11_contract(s11_cfg)
    edges = pd.read_csv(s11_edges_path)
    support = pd.read_csv(s11_support_path)
    require({"dt", "n_all_transitions", "n_common4"}.issubset(support.columns),
            "Step 11 support table missing required columns.")
    step11_transition_rows = int(pd.to_numeric(support["n_all_transitions"]).sum())
    common4_rows = int(pd.to_numeric(support["n_common4"]).sum())
    require(step11_transition_rows == int(transition_report["transition_rows"]),
            "Step-11 total transitions disagree with Step-5 report.")

    if transition_report.get("fraction_common4") is not None:
        expected_fraction = common4_rows / step11_transition_rows
        require(
            math.isclose(expected_fraction, float(transition_report["fraction_common4"]),
                         rel_tol=0, abs_tol=1e-12),
            "Step-11 exact common4 fraction disagrees with Step-5 report.",
        )

    shutil.copy2(s11_edges_path, out / "06_step11_bin_edges.csv")
    shutil.copy2(s11_support_path, out / "07_step11_support_by_dt.csv")
    write_json(out / "06_step11_run_config.json", s11_cfg)

    # Final Step 13.
    s13_cfg_path = step13 / "00_run_config.json"
    s13_core_path = step13 / "01_core_bin_selection.csv"
    s13_cc_path = step13 / "02_complete_case_subjects.csv"
    s13_slope_path = step13 / "05_temporal_slope_summary.csv"
    for f in (s13_cfg_path, s13_core_path, s13_cc_path, s13_slope_path):
        require(f.is_file(), f"Required Step-13 file not found: {f}")
    s13_cfg = json.loads(s13_cfg_path.read_text(encoding="utf-8"))
    validate_step13_contract(s13_cfg)
    core = pd.read_csv(s13_core_path)
    cc = pd.read_csv(s13_cc_path)
    slopes = pd.read_csv(s13_slope_path)

    require("selected_for_primary_core" in core.columns,
            "Step 13 core table lacks selected_for_primary_core.")
    selected_mask = parse_bool_series(core["selected_for_primary_core"],
                                      "selected_for_primary_core")
    selected_bins = sorted(pd.to_numeric(core.loc[selected_mask, "bin"]).astype(int).tolist())
    cfg_bins = [int(x) for x in s13_cfg["core_selection"]["selected_bins"]]
    require(selected_bins == cfg_bins,
            "Step-13 selected core bins disagree between table and run config.")

    require("selected_for_primary_complete_case" in cc.columns,
            "Step 13 complete-case table lacks primary selection flag.")
    cc_mask = parse_bool_series(
        cc["selected_for_primary_complete_case"],
        "selected_for_primary_complete_case",
    )
    primary_subjects = sorted(
        pd.to_numeric(cc.loc[cc_mask, "subject"]).astype(int).tolist()
    )
    cfg_subjects = sorted(int(x) for x in s13_cfg["complete_case"]["subjects"])
    require(primary_subjects == cfg_subjects,
            "Step-13 complete-case subjects disagree between table and run config.")

    primary_slope_rows = slopes.loc[
        (slopes["estimator"].astype(str) == "equal_bin_equal_subject")
        & (slopes["metric"].astype(str) == "cross_cov")
    ]
    require(len(primary_slope_rows) == 1,
            "Expected exactly one Step-13 primary equal-bin/equal-subject cross_cov row.")
    primary_slope = primary_slope_rows.iloc[0]

    shutil.copy2(s13_core_path, out / "08_step13_core_selection.csv")
    shutil.copy2(s13_cc_path, out / "09_step13_complete_case_subjects.csv")
    shutil.copy2(s13_slope_path, out / "10_step13_temporal_slope_summary.csv")
    write_json(out / "08_step13_run_config.json", s13_cfg)

    # Build subject-union templates from ALL raw observed repertoires.
    subject_summary, source_audit = build_subject_union_templates(
        step1=step1,
        raw_design=raw_design,
        analysis_design=analysis_design,
        source_dir=source_dir,
        count_column=args.count_column,
        out=out,
    )
    subject_summary.to_csv(out / "11_subject_union_references.csv", index=False)
    source_audit.to_csv(out / "12_source_repertoire_audit.csv", index=False)

    # Cross-check raw source audit counts against Step 1 totals.
    require(len(source_audit) == n_repertoires,
            f"Audited {len(source_audit)} source repertoires; expected {n_repertoires}")
    require(int(source_audit["depth"].sum()) == int(pd.to_numeric(step1["depth"]).sum()),
            "Source-audit total depth disagrees with Step 1.")

    # Freeze targets.
    observation_targets = {
        "source": "corrected_post_bedtools_real_pipeline",
        "raw_observed_repertoires": n_repertoires,
        "raw_complete_pairs": n_raw_pairs,
        "analysis_eligible_pairs": n_analysis_pairs,
        "analysis_eligible_repertoires": 2 * n_analysis_pairs,
        "analysis_within_subject_interval_pairs": interval_pairs,
        "positive_repertoire_states": positive_rep_states,
        "union_timepoint_states": union_states,
        "both_timepoint_states": both_states,
        "step5_transition_rows": step11_transition_rows,
        "common4_rows": common4_rows,
        "common4_fraction": common4_rows / step11_transition_rows,
        "counts_by_dt": {
            str(int(r.dt)): int(r.n_all_transitions)
            for r in support.itertuples(index=False)
        },
        "common4_by_dt": {
            str(int(r.dt)): int(r.n_common4)
            for r in support.itertuples(index=False)
        },
        "transition_report_n_subject_clone_pairs": transition_report.get("n_subject_clone_pairs"),
        "transition_report_prepared_state_rows": transition_report.get("prepared_state_rows"),
    }
    write_json(out / "04_observation_structure_targets.json", observation_targets)

    # Explicit synthetic-validation contract.
    contract = {
        "schema_version": SCHEMA_VERSION,
        "script_version": SCRIPT_VERSION,
        "reference_universe": {
            "definition": "subject-level union of all raw observed corrected repertoire files",
            "uses_all_raw_observed_repertoires": True,
            "includes_step2_failed_pair_in_reference_if_raw_files_exist": True,
            "temporal_order_used": False,
            "never_observed_tail_inferred": False,
            "primary_frequency": "f_reference_mean_repertoire",
            "sensitivity_frequency": "f_reference_pooled_all_counts",
        },
        "analysis_geometry": {
            "analysis_eligibility_source": "Step2 included_in_downstream AND successful raw complete pair",
            "step2_fit_failure_is_biological_absence": False,
            "analysis_eligible_pairs": n_analysis_pairs,
            "within_subject_interval_pairs": interval_pairs,
        },
        "observation_model_guidance": {
            "legacy_even_count_lattice_retained": False,
            "recommended_count_scale": 1,
            "corrected_counts_may_be_odd": True,
            "calibrate_observation_structure_without_real_cross_cov": True,
            "targets_file": "04_observation_structure_targets.json",
        },
        "step11_frozen_geometry": {
            "support": "common4",
            "conditioning": "xmid_latent",
            "estimand": "Cov(dx_observed_rep1, dx_observed_rep2)",
            "dt_values": [int(x) for x in s11_cfg.get("dt_values", [])],
            "n_bins": int(len(edges)),
            "x_min": float(edges["x_left"].min()),
            "x_max": float(edges["x_right"].max()),
            "bin_edges_file": "06_step11_bin_edges.csv",
        },
        "step13_frozen_primary": {
            "metric": "cross_cov",
            "estimator": "equal_bin_equal_subject",
            "selected_bins": selected_bins,
            "core_x_min": float(s13_cfg["core_selection"]["x_min"]),
            "core_x_max": float(s13_cfg["core_selection"]["x_max"]),
            "complete_case_subjects_real_reference": primary_subjects,
            "n_complete_case_subjects_real_reference": len(primary_subjects),
            "real_primary_slope_per_week_documentary_only": float(
                primary_slope["observed_slope_per_week"]
            ),
            "real_primary_bootstrap_q025_documentary_only": float(
                primary_slope["bootstrap_q025"]
            ),
            "real_primary_bootstrap_q975_documentary_only": float(
                primary_slope["bootstrap_q975"]
            ),
            "real_temporal_result_must_not_be_used_to_calibrate_positive_control": True,
        },
        "next_validation_stage": {
            "null_observation_calibration_parameters": [
                "reference_shape_or_alpha",
                "sigma_fast",
                "sigma_visibility",
            ],
            "forbidden_calibration_targets": [
                "real_cross_cov",
                "real_temporal_slope",
                "positive_control_recovery",
            ],
            "q_bio_selection": "oracle-only after observation-model freeze",
        },
    }
    write_json(out / "13_validation_contract.json", contract)

    # Compact cohort summary.
    raw_missing = raw_design.loc[~raw_design["raw_pair_complete"], ["subject", "time"]]
    analysis_excluded = analysis_design.loc[
        ~analysis_design["analysis_eligible_pair"],
        ["subject", "time", "raw_sampling_status", "step2_success",
         "included_in_downstream", "analysis_exclusion_reason", "step2_message"],
    ]
    cohort_summary = {
        "subjects": n_subjects,
        "scheduled_subject_time_cells": n_subjects * len(times),
        "raw_observed_repertoires": n_repertoires,
        "raw_complete_pairs": n_raw_pairs,
        "raw_missing_or_incomplete_cells": raw_missing.to_dict("records"),
        "analysis_eligible_pairs": n_analysis_pairs,
        "analysis_excluded_cells": analysis_excluded.to_dict("records"),
        "analysis_interval_pairs": interval_pairs,
        "step2_successful_pairs": int(parse_bool_series(qc["success"], "success").sum()),
        "step13_primary_core_bins": selected_bins,
        "step13_primary_complete_case_subjects": primary_subjects,
    }
    write_json(out / "05_cohort_and_analysis_summary.json", cohort_summary)

    # Input provenance and output hashes.
    input_files = [
        {"role": "step1_characteristics", "path": str(step1_path), "sha256": sha256_file(step1_path)},
        {"role": "step2_qc", "path": str(qc_path), "sha256": sha256_file(qc_path)},
        {"role": "step2_params", "path": str(params_path), "sha256": sha256_file(params_path)},
        {"role": "transition_report", "path": str(report_path), "sha256": sha256_file(report_path)},
        {"role": "step11_run_config", "path": str(s11_cfg_path), "sha256": sha256_file(s11_cfg_path)},
        {"role": "step11_edges", "path": str(s11_edges_path), "sha256": sha256_file(s11_edges_path)},
        {"role": "step11_support", "path": str(s11_support_path), "sha256": sha256_file(s11_support_path)},
        {"role": "step13_run_config", "path": str(s13_cfg_path), "sha256": sha256_file(s13_cfg_path)},
        {"role": "step13_core", "path": str(s13_core_path), "sha256": sha256_file(s13_core_path)},
        {"role": "step13_complete_case", "path": str(s13_cc_path), "sha256": sha256_file(s13_cc_path)},
        {"role": "step13_slope_summary", "path": str(s13_slope_path), "sha256": sha256_file(s13_slope_path)},
    ]
    source_hashes = source_audit[
        ["subject", "time", "replica", "resolved_file", "source_sha256"]
    ].to_dict("records")

    outputs = []
    for pth in sorted(out.rglob("*")):
        if pth.is_file() and pth.name != "00_manifest.json":
            outputs.append({
                "path": str(pth.relative_to(out)),
                "sha256": sha256_file(pth),
            })

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "script_version": SCRIPT_VERSION,
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "input_files": input_files,
        "source_repertoire_files": source_hashes,
        "outputs": outputs,
        "calibration_scope": (
            "corrected-data observation-structure freeze and subject-union reference "
            "construction; no synthetic dynamics or recovery fitting"
        ),
        "real_temporal_crosscov_used_for_observation_calibration": False,
        "legacy_count_scale_2_assumption": False,
        "recommended_synthetic_count_scale": 1,
    }
    write_json(out / "00_manifest.json", manifest)

    return {
        "subjects": n_subjects,
        "raw_repertoires": n_repertoires,
        "raw_pairs": n_raw_pairs,
        "analysis_pairs": n_analysis_pairs,
        "interval_pairs": interval_pairs,
        **observation_targets,
        "step13_core_bins": selected_bins,
        "step13_complete_case_subjects": primary_subjects,
    }


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--step1-characteristics", required=True, type=Path)
    p.add_argument("--step2-qc", required=True, type=Path)
    p.add_argument("--step2-params", required=True, type=Path)
    p.add_argument("--source-dir", required=True, type=Path,
                   help="Directory containing the corrected repertoire files used by Step 1.")
    p.add_argument("--transition-report", required=True, type=Path)
    p.add_argument("--step11-dir", required=True, type=Path)
    p.add_argument("--step13-dir", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--count-column", default="readCount")
    p.add_argument("--time-values", default="1,2,3,4,5,6")
    p.add_argument("--expected-subjects", type=int, default=10)
    p.add_argument("--expected-repertoires", type=int, default=114)
    p.add_argument("--expected-raw-pairs", type=int, default=57)
    p.add_argument("--expected-analysis-pairs", type=int, default=56)
    p.add_argument("--expected-interval-pairs", type=int, default=130)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parser().parse_args(argv)
    target = args.out_dir.expanduser().resolve()
    if target.exists():
        print(f"ERROR: Output path already exists; nothing changed: {target}", file=sys.stderr)
        return 2

    stage: Optional[Path] = None
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix=".validation-v5-", dir=str(target.parent)))
        result = build(args, stage)
        require(not target.exists(), f"Output appeared during execution: {target}")
        os.rename(stage, target)
        stage = None
        print("\n[DONE] ClonoDynamics validation calibration v5")
        print(f"  output: {target}")
        print(f"  raw repertoires: {result['raw_repertoires']}")
        print(f"  raw complete pairs: {result['raw_pairs']}")
        print(f"  analysis-eligible pairs: {result['analysis_pairs']}")
        print(f"  interval pairs: {result['interval_pairs']}")
        print(f"  positive repertoire states: {result['positive_repertoire_states']:,}")
        print(f"  union timepoint states: {result['union_timepoint_states']:,}")
        print(f"  BOTH timepoint states: {result['both_timepoint_states']:,}")
        print(f"  Step-5 transitions: {result['step5_transition_rows']:,}")
        print(f"  common4: {result['common4_rows']:,} ({result['common4_fraction']:.6%})")
        print(f"  Step-13 core bins: {result['step13_core_bins']}")
        print(f"  Step-13 complete-case subjects: {result['step13_complete_case_subjects']}")
        print("  recommended synthetic count_scale: 1")
        return 0
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    finally:
        if stage is not None and stage.exists():
            shutil.rmtree(stage)


if __name__ == "__main__":
    raise SystemExit(main())
