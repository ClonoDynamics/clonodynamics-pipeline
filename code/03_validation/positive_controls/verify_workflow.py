#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_workflow.py
==================

Read-only verification helper for the ClonoDynamics support-conditioned
synthetic positive-control workflow under the reorganized repository layout.

Expected repository layout
--------------------------

code/
    01_core/
        1-repertoire_characterization.py
        2-multirepresentation_clonotype_state_inference.py
        noiseK_latent.py
        4-multirepresentation_trajectory_assembly.py
        5-longitudinal_transition_assembly.py

    03_validation/
        11-cross_replicate_fluctuation_dynamics.py
        13-temporal_fluctuation_scaling.py
        positive_controls/
            01_freeze_empirical_calibration.py
            02_calibrate_oracle_doses.py
            03_generate_support_conditioned_repertoires.py
            04_summarize_temporal_recovery.py
            05_plot_synthetic_validation.py
            verify_workflow.py

The historical ``pipeline_reference/13-temporal_fluctuation_scaling.py`` copy
is deliberately NOT used.  The positive-control workflow must reuse the same
canonical Step 13 that is used by the real longitudinal analysis.

Scientific scope
----------------
This helper does not fit latent states, generate synthetic counts, alter support,
or replace any production analysis.  It checks:

* source-code versions, hashes and CLI contracts;
* the corrected empirical filename design (114 repertoires, 57 raw pairs,
  subjects 1..10, nominal times 1..6);
* retention of the intended 56 analysis-eligible subject/time pairs in each
  synthetic Step-2 run;
* equality of the frozen Step-11 abundance grid and support across scenarios;
* equality of the frozen Step-13 abundance core and complete-case cohort;
* consistency of oracle-selected sigma_fast / q_bio / seed with the generator;
* exact preservation of empirical depth and richness targets;
* exact reconstruction of the primary Step-13 joint subject bootstrap;
* reuse of the same subject-bootstrap draws across the four dose scenarios.

The completed-run audit is a contract/provenance verification.  It is not a
byte-for-byte historical software reproduction and it is not a multi-realization
power analysis.

Typical commands
----------------

Preflight before starting a validation workspace::

    python verify_workflow.py preflight \
        --data /path/to/dataset_longitudinal \
        --code-root /path/to/code

Freeze source-code provenance::

    python verify_workflow.py source-inventory \
        --code-root /path/to/code \
        --out /path/to/work/source_code_inventory.json

Verify that the frozen code has not changed::

    python verify_workflow.py check-inventory \
        /path/to/work/source_code_inventory.json

Other orchestration helpers::

    python verify_workflow.py fitted-pairs CALIBRATION_DIR SYNTHETIC_STEP2_DIR
    python verify_workflow.py core-bins CALIBRATION_DIR
    python verify_workflow.py scenario-rows ORACLE_DIR
    python verify_workflow.py completed WORK_DIR

Python >= 3.9.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


VERIFY_VERSION = "v2-split-core-validation-layout-2026-09-23"


# -----------------------------------------------------------------------------
# Canonical code contracts
# -----------------------------------------------------------------------------

CORE_SPECS: Mapping[str, Tuple[Optional[str], Sequence[str]]] = {
    "1-repertoire_characterization.py": (
        "v2-production-descriptive-heavy-tail-2026-09-09",
        ("--data_dir", "--out_dir", "--file_sep", "--pattern"),
    ),
    "2-multirepresentation_clonotype_state_inference.py": (
        "v9-replicate-resolved-state-contract-2026-09-09",
        (
            "--data-dir",
            "--results-dir",
            "--pattern",
            "--file-sep",
            "--gamma-init",
            "--k-init",
            "--grid-size",
            "--chunk-size",
            "--alpha",
            "--null",
            "--tail",
            "--min-pairs-per-subject",
            "--intermediate-format",
            "--parquet-compression",
            "--n-jobs",
            "--n-jobs-profile",
        ),
    ),
    "4-multirepresentation_trajectory_assembly.py": (
        "v12-multirepresentation-trajectory-contract-2026-09-09",
        (
            "--per-clone",
            "--outdir",
            "--reference-alpha",
            "--output-format",
            "--parquet-compression",
            "--streaming",
        ),
    ),
    "5-longitudinal_transition_assembly.py": (
        "v4-multirepresentation-transition-contract-2026-09-09",
        (
            "--trajectories",
            "--out",
            "--min-dt",
            "--max-dt",
            "--parquet-compression",
        ),
    ),
}

VALIDATION_SPECS: Mapping[str, Tuple[Optional[str], Sequence[str]]] = {
    "11-cross_replicate_fluctuation_dynamics.py": (
        "v2-fixed-support-bootstrap-signature-2026-09-17",
        (
            "--transitions",
            "--outdir",
            "--dataset-label",
            "--dt-values",
            "--n-bins",
            "--binning",
            "--q-low",
            "--q-high",
            "--bin-edges-file",
            "--min-n",
            "--n-bootstrap",
            "--bootstrap-seed",
        ),
    ),
    "13-temporal_fluctuation_scaling.py": (
        "v3-final-step11v2-fixed-primary-cohort-2026-09-17",
        (
            "--step11-dir",
            "--outdir",
            "--dt-values",
            "--core-min-subject-fraction",
            "--min-subject-cell-n",
            "--n-bootstrap",
            "--bootstrap-seed",
            "--core-bins",
        ),
    ),
}

POSITIVE_CONTROL_SPECS: Mapping[str, Tuple[Optional[str], Sequence[str]]] = {
    "01_freeze_empirical_calibration.py": (
        "5.0.0-corrected-data-freeze-2026-09-18",
        (
            "--step1-characteristics",
            "--step2-qc",
            "--step2-params",
            "--source-dir",
            "--transition-report",
            "--step11-dir",
            "--step13-dir",
            "--out-dir",
            "--expected-subjects",
            "--expected-repertoires",
            "--expected-raw-pairs",
            "--expected-analysis-pairs",
            "--expected-interval-pairs",
        ),
    ),
    "02_calibrate_oracle_doses.py": (
        "5.0.0-oracle-qbio-support-conditioned-2026-09-19",
        (
            "--calibration-dir",
            "--out-dir",
            "--sigma-fast",
            "--target-Rs",
            "--seed",
            "--iterations",
            "--tolerance",
        ),
    ),
    "03_generate_support_conditioned_repertoires.py": (
        "5.1.0-empirical-support-conditioned-common-random-numbers-2026-09-19",
        (
            "--calibration-dir",
            "--out-dir",
            "--scenario",
            "--sigma-fast",
            "--q-bio",
            "--seed",
            "--source-dir",
        ),
    ),
    "04_summarize_temporal_recovery.py": (
        None,
        (
            "--oracle-dir",
            "--dose-root",
            "--outdir",
            "--no-png",
            "--no-pdf",
            "--no-svg",
        ),
    ),
    "05_plot_synthetic_validation.py": (
        None,
        (
            "--summary-dir",
            "--outdir",
            "--include-supplementary",
            "--pdf",
            "--png",
        ),
    ),
}

SCENARIOS = ("R0p00", "R0p25", "R0p50", "R1p00")
PRIMARY_ESTIMATOR = "equal_bin_equal_subject"

REPERTOIRE_PATTERN = re.compile(
    r"^(?P<subject>\d+)_(?P<time>\d+)-(?P<replica>[12])"
    r"\.(?:tsv|csv|parquet|pq|feather|arrow)$"
)


@dataclass(frozen=True)
class CodeLayout:
    code_root: Path
    core_dir: Path
    validation_dir: Path
    positive_dir: Path


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)



def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for part in iter(lambda: handle.read(chunk_size), b""):
            h.update(part)
    return h.hexdigest()



def read_json(path: Path):
    return json.loads(Path(path).read_text(encoding="utf-8"))



def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")



def declared_version(path: Path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8-sig"))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        if any(name in {"SCRIPT_VERSION", "VERSION"} for name in names):
            try:
                return ast.literal_eval(node.value)
            except (ValueError, TypeError):
                pass
    return None



def infer_code_root_from_script() -> Optional[Path]:
    """Infer code/ when this file lives in 03_validation/positive_controls/."""
    here = Path(__file__).resolve()
    parent = here.parent
    if parent.name in {"positive_controls", "synthetic_validation_methods"}:
        validation = parent.parent
        if validation.name == "03_validation":
            return validation.parent
    return None



def resolve_layout(
    code_root: Optional[Path] = None,
    core_dir: Optional[Path] = None,
    validation_dir: Optional[Path] = None,
    positive_dir: Optional[Path] = None,
) -> CodeLayout:
    inferred_root = infer_code_root_from_script()

    if code_root is not None:
        root = Path(code_root).expanduser().resolve(strict=True)
    elif inferred_root is not None:
        root = inferred_root.expanduser().resolve(strict=True)
    elif core_dir is not None and validation_dir is not None:
        root = Path(core_dir).expanduser().resolve(strict=True).parent
    else:
        raise ValueError(
            "Cannot infer repository code root. Supply --code-root, or place "
            "verify_workflow.py under code/03_validation/positive_controls/."
        )

    core = (
        Path(core_dir).expanduser().resolve(strict=True)
        if core_dir is not None
        else (root / "01_core").resolve(strict=True)
    )
    validation = (
        Path(validation_dir).expanduser().resolve(strict=True)
        if validation_dir is not None
        else (root / "03_validation").resolve(strict=True)
    )

    if positive_dir is not None:
        positive = Path(positive_dir).expanduser().resolve(strict=True)
    else:
        preferred = validation / "positive_controls"
        legacy = validation / "synthetic_validation_methods"
        if preferred.is_dir():
            positive = preferred.resolve(strict=True)
        elif legacy.is_dir():
            # Transitional compatibility while the repository is being cleaned.
            positive = legacy.resolve(strict=True)
        else:
            raise FileNotFoundError(
                f"Neither {preferred} nor {legacy} exists."
            )

    for path, label in [
        (root, "code root"),
        (core, "01_core"),
        (validation, "03_validation"),
        (positive, "positive_controls"),
    ]:
        require(path.is_dir(), f"Missing {label} directory: {path}")

    return CodeLayout(root, core, validation, positive)



def source_entries(layout: CodeLayout) -> List[Dict[str, object]]:
    entries: List[Dict[str, object]] = []

    for name, (expected, _flags) in CORE_SPECS.items():
        entries.append(
            {
                "role": "core",
                "name": name,
                "path": layout.core_dir / name,
                "expected_version": expected,
            }
        )

    entries.append(
        {
            "role": "core_dependency",
            "name": "noiseK_latent.py",
            "path": layout.core_dir / "noiseK_latent.py",
            "expected_version": None,
        }
    )

    for name, (expected, _flags) in VALIDATION_SPECS.items():
        entries.append(
            {
                "role": "validation",
                "name": name,
                "path": layout.validation_dir / name,
                "expected_version": expected,
            }
        )

    for name, (expected, _flags) in POSITIVE_CONTROL_SPECS.items():
        entries.append(
            {
                "role": "positive_control",
                "name": name,
                "path": layout.positive_dir / name,
                "expected_version": expected,
            }
        )

    entries.append(
        {
            "role": "audit_helper",
            "name": Path(__file__).name,
            "path": Path(__file__).resolve(),
            "expected_version": VERIFY_VERSION,
        }
    )
    return entries



def validate_source_files(layout: CodeLayout) -> None:
    for entry in source_entries(layout):
        path = Path(entry["path"])
        require(path.is_file(), f"Missing source file: {path}")
        ast.parse(path.read_text(encoding="utf-8-sig"))

        expected = entry["expected_version"]
        if entry["role"] == "audit_helper":
            found = VERIFY_VERSION
        else:
            found = declared_version(path)
        if expected is not None:
            require(
                found == expected,
                f"{path.name}: expected version {expected!r}; found {found!r}",
            )



def help_text(path: Path, cwd: Path) -> str:
    result = subprocess.run(
        [sys.executable, str(path), "--help"],
        capture_output=True,
        text=True,
        timeout=90,
        cwd=str(cwd),
    )
    require(
        result.returncode == 0,
        f"--help failed for {path}:\n{result.stderr}",
    )
    return (result.stdout or "") + "\n" + (result.stderr or "")



def validate_cli_contracts(layout: CodeLayout) -> None:
    checks: List[Tuple[Path, Sequence[str]]] = []
    checks.extend(
        (layout.core_dir / name, flags)
        for name, (_version, flags) in CORE_SPECS.items()
    )
    checks.extend(
        (layout.validation_dir / name, flags)
        for name, (_version, flags) in VALIDATION_SPECS.items()
    )
    checks.extend(
        (layout.positive_dir / name, flags)
        for name, (_version, flags) in POSITIVE_CONTROL_SPECS.items()
    )

    for path, flags in checks:
        text = help_text(path, layout.code_root)
        missing = [flag for flag in flags if flag not in text]
        require(
            not missing,
            f"{path.name}: missing expected CLI flags {missing}",
        )



def validate_noise_module(layout: CodeLayout) -> None:
    # Avoid a previously imported module shadowing the selected 01_core file.
    old_module = sys.modules.pop("noiseK_latent", None)
    sys.path.insert(0, str(layout.core_dir.resolve()))
    try:
        noise = importlib.import_module("noiseK_latent")
        require(
            Path(noise.__file__).resolve()
            == (layout.core_dir / "noiseK_latent.py").resolve(),
            "Another noiseK_latent module shadows code/01_core/noiseK_latent.py.",
        )
        require(
            callable(getattr(noise, "fit_noiseK_latent_powerlaw", None)),
            "Missing noiseK_latent.fit_noiseK_latent_powerlaw.",
        )
    finally:
        if sys.path and sys.path[0] == str(layout.core_dir.resolve()):
            sys.path.pop(0)
        sys.modules.pop("noiseK_latent", None)
        if old_module is not None:
            sys.modules["noiseK_latent"] = old_module


# -----------------------------------------------------------------------------
# Preflight / provenance
# -----------------------------------------------------------------------------


def validate_empirical_filename_design(data_dir: Path) -> Dict[str, object]:
    pairs: Dict[Tuple[int, int], set] = {}
    n_files = 0

    for path in sorted(data_dir.iterdir()):
        match = REPERTOIRE_PATTERN.fullmatch(path.name) if path.is_file() else None
        if match is None:
            continue

        subject = int(match["subject"])
        time = int(match["time"])
        replicate = int(match["replica"])

        key = (subject, time)
        require(
            replicate not in pairs.setdefault(key, set()),
            f"Duplicate repertoire representation: {(subject, time, replicate)}",
        )
        pairs[key].add(replicate)
        n_files += 1

    require(
        n_files == 114,
        f"Expected 114 corrected empirical files; found {n_files}. "
        "Check paths/names; do not invent missing files.",
    )
    require(
        len(pairs) == 57 and all(reps == {1, 2} for reps in pairs.values()),
        "Expected 57 complete raw replicate pairs.",
    )
    require(
        {subject for subject, _time in pairs} == set(range(1, 11)),
        "Expected empirical subject IDs 1..10.",
    )
    require(
        {time for _subject, time in pairs} == set(range(1, 7)),
        "Expected nominal time IDs 1..6.",
    )

    return {
        "n_repertoires": n_files,
        "n_raw_pairs": len(pairs),
        "subjects": sorted({subject for subject, _ in pairs}),
        "nominal_times": sorted({time for _, time in pairs}),
    }



def preflight(data_dir: Path, layout: CodeLayout) -> None:
    require(sys.version_info >= (3, 9), "Python >=3.9 is required.")
    require(data_dir.is_dir(), f"Dataset directory not found: {data_dir}")

    validate_source_files(layout)
    validate_cli_contracts(layout)

    for package in ("numpy", "pandas", "scipy", "polars", "pyarrow", "matplotlib"):
        module = importlib.import_module(package)
        print(f"{package}: {getattr(module, '__version__', 'unknown')}")

    validate_noise_module(layout)
    empirical = validate_empirical_filename_design(data_dir)

    print("\nRepository layout")
    print("  code root        :", layout.code_root)
    print("  core             :", layout.core_dir)
    print("  validation       :", layout.validation_dir)
    print("  positive controls:", layout.positive_dir)
    print("\nEmpirical design")
    print("  repertoires      :", empirical["n_repertoires"])
    print("  raw pairs        :", empirical["n_raw_pairs"])
    print("  subjects         :", empirical["subjects"])
    print("  nominal times    :", empirical["nominal_times"])
    print(
        "\nPreflight passed: canonical source versions, CLI interfaces, "
        "noise-model import and empirical filename design are consistent."
    )



def source_inventory(layout: CodeLayout, out: Path) -> None:
    validate_source_files(layout)

    records = []
    for entry in source_entries(layout):
        path = Path(entry["path"]).resolve(strict=True)
        version = VERIFY_VERSION if entry["role"] == "audit_helper" else declared_version(path)
        records.append(
            {
                "role": entry["role"],
                "name": entry["name"],
                "path": str(path),
                "sha256": sha256_file(path),
                "declared_version": version,
                "expected_version": entry["expected_version"],
            }
        )

    payload = {
        "verify_workflow_version": VERIFY_VERSION,
        "code_root": str(layout.code_root),
        "core_dir": str(layout.core_dir),
        "validation_dir": str(layout.validation_dir),
        "positive_controls_dir": str(layout.positive_dir),
        "files": records,
    }
    write_json(out, payload)
    print(f"Wrote source inventory: {out}")



def check_inventory(path: Path) -> None:
    payload = read_json(path)
    records = payload.get("files", payload if isinstance(payload, list) else [])
    require(isinstance(records, list) and records, f"No source records in {path}")

    for record in records:
        source = Path(record["path"])
        require(source.is_file(), f"Source missing after initialization: {source}")
        actual = sha256_file(source)
        require(
            actual == record["sha256"],
            f"Source changed after initialization: {source}",
        )

    print(f"Source inventory unchanged: {len(records)} files verified.")


# -----------------------------------------------------------------------------
# Scientific audit helpers retained from the validated workflow
# -----------------------------------------------------------------------------


def booleans(series):
    converted = (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map(
            {
                "true": True,
                "false": False,
                "1": True,
                "0": False,
                "1.0": True,
                "0.0": False,
            }
        )
    )
    require(converted.notna().all(), "Invalid Boolean field.")
    return converted.astype(bool)



def fitted_pairs(calibration_dir: Path, step2_dir: Path) -> None:
    import pandas as pd

    design_path = calibration_dir / "02_analysis_sampling_design.csv"
    qc_path = step2_dir / "latent_pair_qc.csv"
    require(design_path.is_file(), f"Missing calibration design: {design_path}")
    require(qc_path.is_file(), f"Missing synthetic Step-2 QC table: {qc_path}")

    design = pd.read_csv(design_path)
    expected = set(
        map(
            tuple,
            design.loc[
                booleans(design["analysis_eligible_pair"]),
                ["subject", "time"],
            ]
            .astype(int)
            .to_numpy(),
        )
    )

    qc = pd.read_csv(qc_path)
    active = booleans(qc["included_in_downstream"]) & booleans(qc["success"])
    actual = set(
        map(
            tuple,
            qc.loc[active, ["subject", "time"]].astype(int).to_numpy(),
        )
    )

    require(
        actual == expected,
        "Synthetic Step 2 changed the eligible visit set: "
        f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}. "
        "Stop; do not drop visits silently.",
    )
    print(f"All {len(expected)} intended synthetic pairs retained.")



def reconstruct_primary(step13_dir: Path):
    import numpy as np
    import pandas as pd

    cfg = read_json(step13_dir / "00_run_config.json")
    subjects = list(map(int, cfg["complete_case"]["subjects"]))
    dt_values = list(map(int, cfg["dt_values"]))

    core = pd.read_csv(step13_dir / "03_subject_core_metric_by_dt.csv")
    core = core.loc[
        (core["estimator"] == PRIMARY_ESTIMATOR)
        & (core["metric"] == "cross_cov")
    ]

    matrix = (
        core.pivot(index="subject", columns="dt", values="value")
        .reindex(index=subjects, columns=dt_values)
        .to_numpy(float)
    )
    require(
        np.isfinite(matrix).all(),
        f"Incomplete primary subject-by-lag table: {step13_dir}",
    )

    n_bootstrap = int(cfg["bootstrap"]["n_bootstrap"])
    seed = int(cfg["bootstrap"]["seed"])
    draws = np.random.default_rng(seed).integers(
        0,
        len(subjects),
        size=(n_bootstrap, len(subjects)),
    )

    curves = matrix[draws].mean(axis=1)
    x = np.asarray(dt_values, dtype=float)
    x -= x.mean()
    slopes = curves @ x / (x @ x)

    boot = pd.read_csv(step13_dir / "10_joint_subject_bootstrap.csv")
    boot = boot.loc[
        (boot["estimator"] == PRIMARY_ESTIMATOR)
        & (boot["metric"] == "cross_cov")
    ].sort_values("bootstrap")

    require(
        boot["bootstrap"].tolist() == list(range(n_bootstrap)),
        f"Missing/repeated bootstrap IDs: {step13_dir}",
    )
    require(
        np.allclose(
            slopes,
            boot["slope_per_week"].to_numpy(float),
            rtol=1e-10,
            atol=1e-12,
        ),
        f"Primary bootstrap slope reconstruction failed: {step13_dir}",
    )

    for j, dt in enumerate(dt_values):
        require(
            np.allclose(
                curves[:, j],
                boot[f"value_dt{dt}"].to_numpy(float),
                rtol=1e-10,
                atol=1e-12,
            ),
            f"Bootstrap profile mismatch at dt={dt}: {step13_dir}",
        )

    max_error = float(
        np.max(
            np.abs(
                slopes - boot["slope_per_week"].to_numpy(float)
            )
        )
    )
    return cfg, subjects, draws, slopes, max_error



def completed(work_dir: Path) -> None:
    import numpy as np
    import pandas as pd

    calibration = work_dir / "calibration"
    oracle_dir = work_dir / "oracle"
    scenarios_root = work_dir / "scenarios"

    oracle_path = oracle_dir / "01_scenario_table.csv"
    require(oracle_path.is_file(), f"Missing oracle table: {oracle_path}")
    oracle = pd.read_csv(oracle_path).set_index("scenario_label")

    ref_cfg = read_json(calibration / "08_step13_run_config.json")
    ref_edges = pd.read_csv(calibration / "06_step11_bin_edges.csv")
    ref_support = pd.read_csv(calibration / "07_step11_support_by_dt.csv").sort_values("dt")

    checks: List[Dict[str, object]] = []
    first_draws = None
    first_subjects = None
    bootstrap_slopes: Dict[str, object] = {}

    for scenario in SCENARIOS:
        root = scenarios_root / scenario
        require(root.is_dir(), f"Missing scenario directory: {root}")

        fitted_pairs(calibration, root / "2-clonotype_state_inference")

        step11_dir = root / "11-cross_replicate_fluctuation_dynamics"
        step13_dir = root / "13-temporal_fluctuation_scaling"
        cfg11 = read_json(step11_dir / "00_run_config.json")
        cfg13 = read_json(step13_dir / "00_run_config.json")

        require(
            cfg13["upstream_step11_analysis_signature"] == cfg11["analysis_signature"],
            f"{scenario}: Step-11/Step-13 signature mismatch",
        )
        require(
            cfg11["primary_support"] == "common4"
            and cfg11["primary_conditioning"] == "xmid_latent",
            f"{scenario}: incorrect Step-11 estimand support/conditioning",
        )
        require(
            cfg11.get("operational_class_filter") is None,
            f"{scenario}: unexpected operational-class filtering",
        )
        require(
            cfg11["primary_estimands"]["cross_cov_truncated_at_zero"] is False,
            f"{scenario}: cross covariance was truncated",
        )

        edges = pd.read_csv(step11_dir / "01_fluctuation_bin_edges.csv")
        require(
            np.allclose(
                edges[["x_left", "x_right"]],
                ref_edges[["x_left", "x_right"]],
                rtol=0,
                atol=1e-12,
            ),
            f"{scenario}: Step-11 abundance-bin edges changed",
        )

        support = pd.read_csv(step11_dir / "06_support_by_dt.csv").sort_values("dt")
        for column in (
            "dt",
            "n_all_transitions",
            "n_common4",
            "n_subjects",
            "n_subject_intervals",
        ):
            require(
                np.array_equal(
                    support[column].to_numpy(),
                    ref_support[column].to_numpy(),
                ),
                f"{scenario}: support mismatch in {column}",
            )

        require(
            cfg13["core_selection"]["selected_bins"]
            == ref_cfg["core_selection"]["selected_bins"],
            f"{scenario}: temporal core changed",
        )
        require(
            cfg13["complete_case"]["subjects"]
            == ref_cfg["complete_case"]["subjects"],
            f"{scenario}: complete-case subject set changed",
        )
        require(
            int(cfg13["bootstrap"]["seed"]) == 123
            and int(cfg13["bootstrap"]["n_bootstrap"]) == 2000,
            f"{scenario}: unexpected Step-13 bootstrap setup",
        )

        generator_manifest = read_json(root / "00_manifest.json")
        for key in ("sigma_fast", "q_bio"):
            require(
                np.isclose(
                    float(generator_manifest["biology"][key]),
                    float(oracle.loc[scenario, key]),
                    rtol=1e-12,
                    atol=1e-15,
                ),
                f"{scenario}: oracle/generator {key} mismatch",
            )
        require(
            int(generator_manifest["seed"])
            == int(oracle.loc[scenario, "seed"]),
            f"{scenario}: oracle/generator seed mismatch",
        )
        require(
            all(generator_manifest["support"]["frozen_support_audit"].values()),
            f"{scenario}: generator frozen-support audit failure",
        )

        sampling = pd.read_csv(root / "01_generated_sampling_summary.csv")
        require(
            len(sampling) == 112,
            f"{scenario}: expected 112 generated synthetic repertoires; found {len(sampling)}",
        )
        require(
            (sampling["nominal_depth"] == sampling["realized_depth"]).all(),
            f"{scenario}: empirical sequencing depth was not preserved",
        )
        require(
            (
                sampling["target_empirical_richness"]
                == sampling["synthetic_richness"]
            ).all(),
            f"{scenario}: empirical richness target was not preserved",
        )

        cfg, subjects, draws, slopes, error = reconstruct_primary(step13_dir)

        if first_draws is None:
            first_draws = draws
            first_subjects = subjects
        require(
            np.array_equal(first_draws, draws),
            f"{scenario}: subject-bootstrap draws differ from the other scenarios",
        )
        require(
            list(first_subjects) == list(subjects),
            f"{scenario}: bootstrap subject ordering differs from the reference scenario",
        )

        bootstrap_slopes[scenario] = slopes
        checks.append(
            {
                "scenario": scenario,
                "status": "PASS",
                "n_primary_bootstrap": len(slopes),
                "max_abs_reconstruction_error": error,
                "n_subjects": len(subjects),
                "n_core_bins": len(cfg["core_selection"]["selected_bins"]),
            }
        )

    require(first_draws is not None, "No scenario bootstrap draws reconstructed.")
    require(first_subjects is not None, "No scenario subject set reconstructed.")

    verification = work_dir / "verification"
    require(
        not verification.exists(),
        f"Verification output already exists: {verification}. "
        "Delete/rewrite it explicitly before rerunning this audit.",
    )
    verification.mkdir(parents=True)

    pd.DataFrame(checks).to_csv(
        verification / "verification_checks.csv",
        index=False,
    )

    subject_ids = np.asarray(first_subjects)
    draw_frame = pd.DataFrame(
        subject_ids[first_draws],
        columns=[f"draw_{i + 1}" for i in range(len(first_subjects))],
    )
    draw_frame.insert(0, "bootstrap", np.arange(len(first_draws)))
    draw_frame.to_csv(
        verification / "bootstrap_subject_draws.csv",
        index=False,
    )

    paired = pd.DataFrame({"bootstrap": np.arange(len(first_draws))})
    for scenario in SCENARIOS[1:]:
        paired[f"{scenario}_minus_R0p00"] = (
            bootstrap_slopes[scenario] - bootstrap_slopes["R0p00"]
        )
    paired.to_csv(
        verification / "verified_paired_increments.csv",
        index=False,
    )

    write_json(
        verification / "verification_summary.json",
        {
            "status": "PASS",
            "verify_workflow_version": VERIFY_VERSION,
            "checks": checks,
            "scope": (
                "support, frozen abundance geometry/core/cohort, generator-oracle "
                "contract and primary Step-13 bootstrap reconstruction; not a "
                "historical byte-for-byte reproduction"
            ),
        },
    )

    print(
        "All four scenarios passed support/core/cohort checks, generator-oracle "
        "consistency and primary bootstrap reconstruction. Summarization may proceed."
    )


# -----------------------------------------------------------------------------
# Small machine-readable helpers used by the orchestrator
# -----------------------------------------------------------------------------


def core_bins(calibration_dir: Path) -> List[int]:
    cfg = read_json(calibration_dir / "08_step13_run_config.json")
    bins = [int(x) for x in cfg["core_selection"]["selected_bins"]]
    require(bool(bins), "Frozen empirical Step-13 core is empty.")
    return bins



def scenario_rows(oracle_dir: Path) -> List[Tuple[str, float, float, int]]:
    import pandas as pd

    table = pd.read_csv(oracle_dir / "01_scenario_table.csv").sort_values("target_R")
    require(
        tuple(table["scenario_label"].astype(str)) == SCENARIOS,
        "Unexpected oracle scenario labels/order.",
    )
    rows = []
    for row in table.itertuples(index=False):
        rows.append(
            (
                str(row.scenario_label),
                float(row.sigma_fast),
                float(row.q_bio),
                int(row.seed),
            )
        )
    return rows


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def add_layout_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--code-root",
        type=Path,
        default=None,
        help=(
            "Repository code/ directory containing 01_core and 03_validation. "
            "Normally inferred automatically when this script lives under "
            "03_validation/positive_controls/."
        ),
    )
    parser.add_argument("--core-dir", type=Path, default=None)
    parser.add_argument("--validation-dir", type=Path, default=None)
    parser.add_argument("--positive-dir", type=Path, default=None)



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Verify the ClonoDynamics support-conditioned synthetic positive-control "
            "workflow under the split 01_core / 03_validation repository layout."
        ),
    )
    parser.add_argument("--version", action="version", version=VERIFY_VERSION)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("preflight", help="Check code contracts and empirical dataset layout.")
    p.add_argument("--data", required=True, type=Path)
    add_layout_arguments(p)

    p = sub.add_parser("source-inventory", help="Freeze code paths, versions and SHA-256 hashes.")
    p.add_argument("--out", required=True, type=Path)
    add_layout_arguments(p)

    p = sub.add_parser("check-inventory", help="Verify that a frozen source inventory has not changed.")
    p.add_argument("inventory", type=Path)

    p = sub.add_parser("fitted-pairs", help="Verify retention of the frozen 56 analysis-eligible visits.")
    p.add_argument("calibration_dir", type=Path)
    p.add_argument("step2_dir", type=Path)

    p = sub.add_parser("core-bins", help="Print the frozen empirical Step-13 core-bin IDs.")
    p.add_argument("calibration_dir", type=Path)

    p = sub.add_parser("scenario-rows", help="Print oracle scenario rows for orchestration.")
    p.add_argument("oracle_dir", type=Path)

    p = sub.add_parser("completed", help="Audit all completed synthetic scenarios and write verification outputs.")
    p.add_argument("work_dir", type=Path)

    return parser



def layout_from_args(args: argparse.Namespace) -> CodeLayout:
    return resolve_layout(
        code_root=getattr(args, "code_root", None),
        core_dir=getattr(args, "core_dir", None),
        validation_dir=getattr(args, "validation_dir", None),
        positive_dir=getattr(args, "positive_dir", None),
    )



def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "preflight":
        data = args.data.expanduser().resolve(strict=True)
        preflight(data, layout_from_args(args))

    elif args.command == "source-inventory":
        out = args.out.expanduser().resolve()
        source_inventory(layout_from_args(args), out)

    elif args.command == "check-inventory":
        check_inventory(args.inventory.expanduser().resolve(strict=True))

    elif args.command == "fitted-pairs":
        fitted_pairs(
            args.calibration_dir.expanduser().resolve(strict=True),
            args.step2_dir.expanduser().resolve(strict=True),
        )

    elif args.command == "core-bins":
        bins = core_bins(args.calibration_dir.expanduser().resolve(strict=True))
        print(" ".join(map(str, bins)))

    elif args.command == "scenario-rows":
        for label, sigma, q_bio, seed in scenario_rows(
            args.oracle_dir.expanduser().resolve(strict=True)
        ):
            print(label, repr(sigma), repr(q_bio), seed, sep="\t")

    elif args.command == "completed":
        completed(args.work_dir.expanduser().resolve(strict=True))

    else:  # pragma: no cover
        raise ValueError(f"Unknown command: {args.command}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
