#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orchestrator_longitudinal_vs_pseudo.py
======================================

Interactive, dependency-aware orchestrator for the cross-dataset
ClonoDynamics comparisons between genuine longitudinal and pseudo-longitudinal
reference data.

Managed public Steps
--------------------
    Step 10
    10-longitudinal_vs_pseudo_forward_drift.py

    Step 12
    12-longitudinal_vs_pseudo_fluctuation_baseline.py

These are deliberately kept outside the genuine-longitudinal orchestrator
because both require inputs from TWO independently processed result trees.

Scientific dependency graph
----------------------------

Forward-drift comparison:

    LONGITUDINAL                              PSEUDO-LONGITUDINAL

    Step 9                                    Step 7
    replicate-decoupled                      pseudo conditioning /
    forward drift                            ordering benchmark
        |                                         |
        +------------------+----------------------+
                           |
                           v
                        STEP 10
              longitudinal vs pseudo
                   forward drift


Fluctuation-baseline comparison:

    LONGITUDINAL                              PSEUDO-LONGITUDINAL

    Step 11                                   Step 8
    fluctuation dynamics                     representation assessment
        |                                         |
        |                                         |
        |                                    Step 4 trajectories
        |                                         |
        +------------------+----------------------+
                           |
                           v
                        STEP 12
              longitudinal vs pseudo
                fluctuation baseline


Exact computational dependencies
--------------------------------
Step 10 requires:

    longitudinal Step-9 output directory containing:
        forward_drift_by_bin.csv
        forward_drift_by_fold.csv
        forward_drift_fold_concordance.csv
        forward_drift_subject_by_bin.csv
        forward_drift_loso_by_bin.csv

    pseudo Step-7 output directory containing:
        03_ordering_robustness/
            pseudo_reverse_pair_envelope.csv
            pseudo_reverse_pair_curves.csv

Step 10 compares the already-estimated longitudinal and pseudo
replicate-decoupled forward-drift curves on their common absolute x0 support.
It does not percentile-normalize or force common bin edges.

Step 12 requires:

    longitudinal Step-11 output directory containing:
        10_binned_dynamics_long.csv/parquet

    pseudo Step-8 output directory containing:
        02_conditioning_bin_edges.csv
        03_pooled_representation_curves.csv   [descriptive MAD reference]

    pseudo Step-4:
        latent_trajectories_long.parquet

Step 12 constructs arbitrary pseudo-time orderings from the pseudo trajectories
using the fixed Step-8 xstar bin geometry, then compares genuine longitudinal
dt=1 sample Var(dx|xstar) and MSD with the pseudo technical baseline on common
absolute xstar support.

Step 12 can reuse its own previously generated:
    02_pseudo_directed_pair_sufficient.parquet
during an incomplete resume instead of rebuilding those sufficient statistics
from the pseudo trajectories.

Execution choices
-----------------
Interactive mode offers:

    [1] Steps 10 and 12
    [2] Step 10 only
    [3] Step 12 only

Result roots
------------
Three roots are conceptually distinct:

    1. longitudinal results root
    2. pseudo-longitudinal results root
    3. comparison results root

The first two are READ-ONLY upstream trees for this orchestrator.

The comparison root receives:

    10-longitudinal_vs-pseudo-forward-drift/
    12-longitudinal-vs-pseudo-fluctuation-baseline/

    logs/
        step10.log
        step12.log

    00_orchestrator_longitudinal_vs_pseudo_config.json
    00_orchestrator_longitudinal_vs_pseudo_manifest.csv

    orchestrator_runs_longitudinal_vs_pseudo/

If --comparison-root is omitted, a default is proposed/derived under the common
ancestor of the longitudinal and pseudo result trees:

    <common-ancestor>/clonodynamics_longitudinal_vs_pseudo_results

Root normalization
------------------
For convenience:

    .../<dataset>_clonodynamics_results

is accepted when it contains the standard child:

    longitudinal/

or:

    pseudo-longitudinal/

The orchestrator automatically descends into the appropriate child.

Provenance safeguards
---------------------
When the Steps-1-6 orchestrator configs are available, this orchestrator checks:

    longitudinal root:
        dataset_type == longitudinal

    pseudo root:
        dataset_type == pseudo-longitudinal

It also checks cross-dataset compatibility when metadata are available:

    observability_alpha
        must match across the two upstream pipelines;

    common upstream script hashes
        step2
        step4
        step5
        noise_model.

Both Steps-1-6 provenance schemas are supported:

    v1 complete-run:
        script_sha256

    v2 selective-run:
        selected_script_sha256.

For v2 selective runs, the active config may contain only the most recently
selected Step. The comparison orchestrator therefore reconstructs the latest
known Step-2/4/5/noise-model hashes from immutable:

    orchestrator_runs/*_config.json

and uses the active config as the final/highest-priority record.

A mismatch stops the comparison by default because it would mix data processed
under different latent-inference/transition definitions.

Use:
    --allow-provenance-mismatch
only for an intentional sensitivity analysis.

Primary production settings
----------------------------
Step 10:
    representation = cross_combined
    shared-grid points = 101

Step 12:
    longitudinal dt = 1
    pseudo configurations = 2000
    seed = 123
    pseudo min n per bin = 200
    pseudo min subjects = 3
    stable pseudo-bin fraction = 0.95
    integrated grid-coverage fraction = 0.95
    shared-grid points = 101

Step 12 does NOT subtract the pseudo baseline before Step 13 temporal fitting.
It is a technical-baseline comparison, not a preprocessing correction.

Fresh / resume
--------------
Fresh:
    deletes ONLY the selected Step-10/12 comparison output directories and
    selected logs. Longitudinal and pseudo upstream result trees are never
    modified.

Resume:
    skips a selected comparison when its completion markers are already
    present.

    If Step 12 is incomplete but its
    02_pseudo_directed_pair_sufficient.parquet
    exists, the orchestrator reuses it automatically.

Fail-fast
---------
Each Step runs as a subprocess using the same Python executable used for this
orchestrator.

stdout/stderr are streamed both to terminal and a Step-specific log.

A non-zero return code or missing completion marker stops the selected workflow.

Interactive use
---------------
    python3 orchestrator_longitudinal_vs_pseudo.py

Dry run
-------
    python3 orchestrator_longitudinal_vs_pseudo.py --dry-run

Non-interactive example
-----------------------
    python3 orchestrator_longitudinal_vs_pseudo.py \
        --workflow all \
        --longitudinal-root /path/to/longitudinal \
        --pseudo-root /path/to/pseudo-longitudinal \
        --comparison-root /path/to/comparison-results \
        --yes \
        --non-interactive
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


VERSION = "v2-longitudinal-vs-pseudo-provenance-2026-08-25"

STEP10_SCRIPT = "10-longitudinal_vs_pseudo_forward_drift.py"
STEP12_SCRIPT = "12-longitudinal_vs_pseudo_fluctuation_baseline.py"

STEP10_OUTDIR = "10-longitudinal-vs-pseudo-forward-drift"
STEP12_OUTDIR = "12-longitudinal-vs-pseudo-fluctuation-baseline"

LONG_STEP9_DIR = "9-replicate_decoupled_forward_drift"
LONG_STEP11_DIR = "11-longitudinal_fluctuation_dynamics"

PSEUDO_STEP4_DIR = "4-latent_trajectory_construction"
PSEUDO_STEP7_DIR = "7-pseudo_longitudinal_conditioning_benchmark"
PSEUDO_STEP8_DIR = "8-pseudo_longitudinal_representation_assessment"

TRAJECTORY_BASENAME = "latent_trajectories_long.parquet"

STEP9_MARKERS = (
    "forward_drift_by_bin.csv",
    "forward_drift_by_fold.csv",
    "forward_drift_fold_concordance.csv",
    "forward_drift_subject_by_bin.csv",
    "forward_drift_loso_by_bin.csv",
)

STEP7_MARKERS = (
    "03_ordering_robustness/pseudo_reverse_pair_envelope.csv",
    "03_ordering_robustness/pseudo_reverse_pair_curves.csv",
)

STEP11_MARKERS = (
    "00_analysis_metadata.json",
)

# Either CSV or Parquet is accepted for 10_binned_dynamics_long.
STEP11_TABLE_STEMS = (
    "10_binned_dynamics_long.csv",
    "10_binned_dynamics_long.parquet",
)

STEP8_MARKERS = (
    "02_conditioning_bin_edges.csv",
)

# This file is optional for Step 12 in the source script but is expected from
# the manuscript-primary Step-8 production run and supplies descriptive MAD.
STEP8_REFERENCE = "03_pooled_representation_curves.csv"

STEP10_COMPLETION = (
    "00_run_config.json",
    "04_empirical_tests.csv",
    "09_comparison_summary.csv",
    "10_longitudinal_vs_pseudo_report.json",
)

STEP12_COMPLETION = (
    "00_run_config.json",
    "07_shared_grid_comparison.csv",
    "09_empirical_tests.csv",
    "10_comparison_summary.csv",
    "12_fluctuation_baseline_report.json",
)

STEP12_SUFFICIENT_PARQUET = "02_pseudo_directed_pair_sufficient.parquet"


# =============================================================================
# Generic utilities
# =============================================================================


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def timestamp_slug() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def shell_join(command: Sequence[str]) -> str:
    try:
        return shlex.join([str(x) for x in command])
    except AttributeError:
        return " ".join(shlex.quote(str(x)) for x in command)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def normalize_workflow(value: str) -> Tuple[int, ...]:
    v = str(value).strip().lower().replace(" ", "")
    mapping = {
        "1": (10, 12),
        "all": (10, 12),
        "10,12": (10, 12),
        "10-12": (10, 12),
        "2": (10,),
        "10": (10,),
        "step10": (10,),
        "3": (12,),
        "12": (12,),
        "step12": (12,),
    }
    if v not in mapping:
        raise ValueError("workflow must be all, 10, or 12")
    return mapping[v]


def workflow_label(steps: Sequence[int]) -> str:
    values = tuple(int(x) for x in steps)
    if values == (10, 12):
        return "Steps 10 and 12"
    return f"Step {values[0]}"


def prompt_workflow() -> Tuple[int, ...]:
    print("\nLongitudinal-vs-pseudo workflow:")
    print("  [1] Steps 10 and 12")
    print("  [2] Step 10 only — forward drift")
    print("  [3] Step 12 only — fluctuation baseline")

    while True:
        raw = input("Selection [1]: ").strip()
        if not raw:
            return (10, 12)
        try:
            return normalize_workflow(raw)
        except ValueError:
            print("Please enter 1, 2, or 3.")


def step_outdir(comparison_root: Path, step: int) -> Path:
    if int(step) == 10:
        return comparison_root / STEP10_OUTDIR
    if int(step) == 12:
        return comparison_root / STEP12_OUTDIR
    raise ValueError(step)


# =============================================================================
# Upstream root normalization / provenance
# =============================================================================


def looks_like_long_root(path: Path) -> bool:
    return any(
        (path / name).exists()
        for name in (
            "00_orchestrator_config.json",
            LONG_STEP9_DIR,
            LONG_STEP11_DIR,
            "4-latent_trajectory_construction",
            "5-latent_transition_construction",
        )
    )


def looks_like_pseudo_root(path: Path) -> bool:
    return any(
        (path / name).exists()
        for name in (
            "00_orchestrator_config.json",
            PSEUDO_STEP4_DIR,
            PSEUDO_STEP7_DIR,
            PSEUDO_STEP8_DIR,
            "5-latent_transition_construction",
        )
    )


def normalize_longitudinal_root(path: Path) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(path)

    if looks_like_long_root(path):
        return path

    child = path / "longitudinal"
    if child.is_dir() and looks_like_long_root(child):
        print(f"[root] Using longitudinal child:\n  {child}")
        return child.resolve()

    return path


def normalize_pseudo_root(path: Path) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(path)

    if looks_like_pseudo_root(path):
        return path

    child = path / "pseudo-longitudinal"
    if child.is_dir() and looks_like_pseudo_root(child):
        print(f"[root] Using pseudo-longitudinal child:\n  {child}")
        return child.resolve()

    return path


def prompt_existing_root(label: str, normalizer) -> Path:
    while True:
        raw = input(f"\n{label}: ").strip()
        if not raw:
            print("A path is required.")
            continue
        try:
            return normalizer(Path(raw))
        except Exception as exc:
            print(f"Invalid path: {exc}")


def read_pipeline_config(root: Path) -> Optional[Dict[str, object]]:
    path = root / "00_orchestrator_config.json"
    if not path.exists():
        return None
    try:
        return load_json(path)
    except Exception as exc:
        print(f"[warning] Could not read {path}: {exc}")
        return None


def dataset_type_from_config(
    config: Optional[Dict[str, object]],
) -> Optional[str]:
    if not config:
        return None
    value = config.get("dataset_type")
    if value is None:
        return None
    return str(value).strip().lower().replace("_", "-")


def config_alpha(
    config: Optional[Dict[str, object]],
) -> Optional[float]:
    if not config:
        return None
    value = config.get("observability_alpha")
    try:
        alpha = float(value)
    except Exception:
        return None
    return alpha if 0 < alpha < 1 else None


def hashes_from_one_config(
    config: Optional[Dict[str, object]],
) -> Dict[str, str]:
    """
    Extract script hashes from either Steps-1-6 orchestrator schema.

    v1 complete-run schema:
        script_sha256

    v2 selective-run schema:
        selected_script_sha256
    """
    if not config:
        return {}

    merged: Dict[str, str] = {}

    for field in ("script_sha256", "selected_script_sha256"):
        value = config.get(field)
        if not isinstance(value, dict):
            continue
        for key, hash_value in value.items():
            if hash_value:
                merged[str(key)] = str(hash_value)

    return merged


def config_history_paths(root: Path) -> List[Path]:
    """
    Return Steps-1-6 provenance configs in chronological order.

    The v2 selective orchestrator stores immutable per-run configs under:
        orchestrator_runs/*_config.json

    The active 00_orchestrator_config.json is appended last so it has highest
    priority when it contains a hash for the same Step.
    """
    paths: List[Path] = []

    history_dir = root / "orchestrator_runs"
    if history_dir.is_dir():
        paths.extend(
            sorted(
                (
                    p for p in history_dir.glob("*_config.json")
                    if p.is_file()
                ),
                key=lambda p: (p.stat().st_mtime, p.name),
            )
        )

    active = root / "00_orchestrator_config.json"
    if active.is_file():
        paths.append(active)

    return paths


def reconstruct_pipeline_hashes(
    root: Path,
    active_config: Optional[Dict[str, object]],
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, List[str]]]:
    """
    Reconstruct the latest known hash for each Steps-1-6 component.

    This is needed because the selective Steps-1-6 orchestrator intentionally
    records only the scripts selected in each run. The active config can
    therefore contain, for example, only Step 6 after a Step-6-only rerun.

    We walk immutable run configs chronologically and keep the most recent hash
    observed for each component. This corresponds to the latest known execution
    of that Step in the result tree.

    Returns
    -------
    hashes
        latest known hash by component key
    sources
        config file supplying the latest known hash
    observed_values
        all distinct historical hashes observed for each component
    """
    hashes: Dict[str, str] = {}
    sources: Dict[str, str] = {}
    observed_values: Dict[str, List[str]] = {}

    paths = config_history_paths(root)

    # If there is no on-disk active config but one was already supplied by the
    # caller, still inspect it.
    configs: List[Tuple[str, Dict[str, object]]] = []

    for path in paths:
        try:
            configs.append((str(path), load_json(path)))
        except Exception as exc:
            print(f"[warning] Could not read provenance config {path}: {exc}")

    active_path = root / "00_orchestrator_config.json"
    if active_config is not None and not active_path.is_file():
        configs.append(("<active-config-object>", active_config))

    for source, config in configs:
        for key, value in hashes_from_one_config(config).items():
            hashes[key] = value
            sources[key] = source
            seen = observed_values.setdefault(key, [])
            if value not in seen:
                seen.append(value)

    return hashes, sources, observed_values


def validate_cross_dataset_provenance(
    *,
    longitudinal_root: Path,
    pseudo_root: Path,
    long_cfg: Optional[Dict[str, object]],
    pseudo_cfg: Optional[Dict[str, object]],
    allow_mismatch: bool,
) -> Dict[str, object]:
    issues: List[str] = []
    warnings: List[str] = []

    long_type = dataset_type_from_config(long_cfg)
    pseudo_type = dataset_type_from_config(pseudo_cfg)

    if long_type is not None and long_type != "longitudinal":
        issues.append(
            f"longitudinal root config reports dataset_type={long_type!r}"
        )

    if pseudo_type is not None and pseudo_type not in {
        "pseudo-longitudinal",
        "pseudolongitudinal",
    }:
        issues.append(
            f"pseudo root config reports dataset_type={pseudo_type!r}"
        )

    if long_cfg is None:
        warnings.append(
            "longitudinal 00_orchestrator_config.json not available"
        )
    if pseudo_cfg is None:
        warnings.append(
            "pseudo 00_orchestrator_config.json not available"
        )

    long_alpha = config_alpha(long_cfg)
    pseudo_alpha = config_alpha(pseudo_cfg)

    if long_alpha is not None and pseudo_alpha is not None:
        if abs(long_alpha - pseudo_alpha) > 1e-12:
            issues.append(
                "observability_alpha mismatch: "
                f"longitudinal={long_alpha}, pseudo={pseudo_alpha}"
            )
    else:
        warnings.append(
            "observability_alpha could not be compared across both roots"
        )

    (
        long_hashes,
        long_hash_sources,
        long_hash_history,
    ) = reconstruct_pipeline_hashes(longitudinal_root, long_cfg)

    (
        pseudo_hashes,
        pseudo_hash_sources,
        pseudo_hash_history,
    ) = reconstruct_pipeline_hashes(pseudo_root, pseudo_cfg)

    compared_hash_keys: List[str] = []
    mismatched_hash_keys: List[str] = []
    unavailable_hash_keys: List[str] = []

    for key in ("step2", "step4", "step5", "noise_model"):
        if key in long_hashes and key in pseudo_hashes:
            compared_hash_keys.append(key)
            if long_hashes[key] != pseudo_hashes[key]:
                mismatched_hash_keys.append(key)
        else:
            unavailable_hash_keys.append(key)

    if mismatched_hash_keys:
        details = []
        for key in mismatched_hash_keys:
            details.append(
                "{}: longitudinal={} [{}], pseudo={} [{}]".format(
                    key,
                    long_hashes.get(key),
                    long_hash_sources.get(key, "unknown source"),
                    pseudo_hashes.get(key),
                    pseudo_hash_sources.get(key, "unknown source"),
                )
            )
        issues.append(
            "upstream script-hash mismatch for: "
            + ", ".join(mismatched_hash_keys)
            + "\n    "
            + "\n    ".join(details)
        )

    if unavailable_hash_keys:
        warnings.append(
            "provenance hashes unavailable on one or both sides for: "
            + ", ".join(unavailable_hash_keys)
        )

    if compared_hash_keys:
        print(
            "[provenance] matched comparable upstream hash keys available: "
            + ", ".join(compared_hash_keys)
        )
    else:
        warnings.append(
            "no common Step-2/4/5/noise-model hashes were available for comparison"
        )

    for warning in warnings:
        print(f"[warning] {warning}")

    if issues:
        message = (
            "Cross-dataset provenance incompatibility detected:\n"
            + "\n".join(f"  - {issue}" for issue in issues)
        )
        if not allow_mismatch:
            raise ValueError(
                message
                + "\nUse --allow-provenance-mismatch only if this is intentional."
            )
        print("[warning] " + message.replace("\n", "\n[warning] "))

    return {
        "longitudinal_root": str(longitudinal_root),
        "pseudo_root": str(pseudo_root),
        "longitudinal_dataset_type": long_type,
        "pseudo_dataset_type": pseudo_type,
        "longitudinal_observability_alpha": long_alpha,
        "pseudo_observability_alpha": pseudo_alpha,
        "compared_upstream_hash_keys": compared_hash_keys,
        "mismatched_upstream_hash_keys": mismatched_hash_keys,
        "unavailable_upstream_hash_keys": unavailable_hash_keys,
        "longitudinal_latest_upstream_hashes": long_hashes,
        "pseudo_latest_upstream_hashes": pseudo_hashes,
        "longitudinal_hash_sources": long_hash_sources,
        "pseudo_hash_sources": pseudo_hash_sources,
        "longitudinal_historical_hash_values": long_hash_history,
        "pseudo_historical_hash_values": pseudo_hash_history,
        "issues": issues,
        "warnings": warnings,
        "allow_provenance_mismatch": bool(allow_mismatch),
    }


# =============================================================================
# Comparison-root resolution
# =============================================================================


def default_comparison_root(
    longitudinal_root: Path,
    pseudo_root: Path,
) -> Path:
    try:
        common = Path(
            os.path.commonpath(
                [str(longitudinal_root), str(pseudo_root)]
            )
        )
    except Exception:
        common = longitudinal_root.parent

    if str(common) in {"/", ""}:
        common = longitudinal_root.parent

    return (
        common / "clonodynamics_longitudinal_vs_pseudo_results"
    ).resolve()


def prompt_comparison_root(default: Path) -> Path:
    raw = input(
        "\nComparison results root "
        f"[{default}]: "
    ).strip()

    path = Path(raw).expanduser() if raw else default
    return path.resolve()


# =============================================================================
# Code/script validation
# =============================================================================


def resolve_code_dir(raw: Optional[Path]) -> Path:
    if raw is None:
        path = Path(__file__).resolve().parent
    else:
        path = raw.expanduser().resolve(strict=True)

    if not path.is_dir():
        raise NotADirectoryError(path)

    return path


def validate_scripts(
    code_dir: Path,
    selected_steps: Sequence[int],
) -> Dict[int, Path]:
    scripts: Dict[int, Path] = {}
    missing: List[Path] = []

    if 10 in selected_steps:
        scripts[10] = code_dir / STEP10_SCRIPT
    if 12 in selected_steps:
        scripts[12] = code_dir / STEP12_SCRIPT

    for path in scripts.values():
        if not path.is_file():
            missing.append(path)

    if missing:
        raise FileNotFoundError(
            "Required comparison scripts are missing from --code-dir:\n"
            + "\n".join(f"  - {path}" for path in missing)
        )

    return scripts


# =============================================================================
# Input discovery helpers
# =============================================================================


def validate_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"{label} missing/empty:\n  {path}")
    return path


def validate_directory_markers(
    path: Path,
    *,
    label: str,
    markers: Sequence[str],
) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(f"{label}: {path}")

    missing = [
        marker
        for marker in markers
        if not (path / marker).is_file()
        or (path / marker).stat().st_size <= 0
    ]

    if missing:
        raise FileNotFoundError(
            f"{label} incomplete; missing/empty markers:\n"
            + "\n".join(f"  - {marker}" for marker in missing)
            + f"\nDirectory:\n  {path}"
        )

    return path


def validate_step11_dir(path: Path) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(path)

    for marker in STEP11_MARKERS:
        if not (path / marker).is_file():
            raise FileNotFoundError(
                f"Step-11 directory missing {marker}:\n  {path}"
            )

    if not any((path / name).is_file() for name in STEP11_TABLE_STEMS):
        raise FileNotFoundError(
            "Step-11 directory lacks 10_binned_dynamics_long.csv/parquet:\n"
            f"  {path}"
        )

    return path


def validate_step8_dir(path: Path) -> Path:
    path = path.expanduser().resolve(strict=True)
    if not path.is_dir():
        raise NotADirectoryError(path)

    for marker in STEP8_MARKERS:
        if not (path / marker).is_file():
            raise FileNotFoundError(
                f"Step-8 directory missing {marker}:\n  {path}"
            )

    # Optional in the source code, but expected for the production reference.
    reference = path / STEP8_REFERENCE
    if not reference.is_file():
        print(
            "[warning] Step-8 descriptive reference is absent: "
            f"{reference}. Step 12 can still run; MAD reference will be unavailable."
        )

    return path


def candidate_directories(
    root: Path,
    *,
    markers: Sequence[str],
) -> List[Path]:
    parent_sets: Optional[set] = None

    for marker in markers:
        basename = Path(marker).name
        relative_parent = Path(marker).parent

        current: set = set()
        for found in root.rglob(basename):
            if not found.is_file() or found.stat().st_size <= 0:
                continue

            # Reconstruct candidate root so nested marker paths are supported.
            candidate = found.parent
            for _ in relative_parent.parts:
                candidate = candidate.parent

            marker_path = candidate / marker
            if marker_path.is_file() and marker_path.stat().st_size > 0:
                current.add(candidate.resolve())

        if parent_sets is None:
            parent_sets = current
        else:
            parent_sets &= current

    return sorted(
        parent_sets or set(),
        key=lambda p: (len(p.parts), str(p)),
    )


def file_candidates(root: Path, basename: str) -> List[Path]:
    return sorted(
        {
            path.resolve()
            for path in root.rglob(basename)
            if path.is_file() and path.stat().st_size > 0
        },
        key=lambda p: (len(p.parts), str(p)),
    )


def prompt_dir_choice(
    label: str,
    candidates: Sequence[Path],
    validator,
) -> Path:
    print(f"\nRequired upstream directory: {label}")

    if candidates:
        print("Candidates found:")
        for i, path in enumerate(candidates, 1):
            print(f"  [{i}] {path}")
        print("  [0] enter another path")

        while True:
            raw = input("Selection: ").strip()
            try:
                choice = int(raw)
            except ValueError:
                print("Enter a candidate number.")
                continue

            if choice == 0:
                break

            if 1 <= choice <= len(candidates):
                return validator(candidates[choice - 1])

            print("Invalid candidate number.")

    while True:
        raw = input(f"Path for {label}: ").strip()
        if not raw:
            print("A path is required.")
            continue

        try:
            return validator(Path(raw))
        except Exception as exc:
            print(exc)


def prompt_file_choice(
    label: str,
    candidates: Sequence[Path],
) -> Path:
    print(f"\nRequired upstream input: {label}")

    if candidates:
        print("Candidates found:")
        for i, path in enumerate(candidates, 1):
            print(f"  [{i}] {path}")
        print("  [0] enter another path")

        while True:
            raw = input("Selection: ").strip()
            try:
                choice = int(raw)
            except ValueError:
                print("Enter a candidate number.")
                continue

            if choice == 0:
                break

            if 1 <= choice <= len(candidates):
                return candidates[choice - 1]

            print("Invalid candidate number.")

    while True:
        raw = input(f"Path for {label}: ").strip()
        if not raw:
            print("A path is required.")
            continue
        try:
            return validate_file(Path(raw), label)
        except Exception as exc:
            print(exc)


def resolve_marker_dir(
    *,
    root: Path,
    canonical: Path,
    label: str,
    markers: Sequence[str],
    override: Optional[Path],
    interactive: bool,
    validator=None,
) -> Path:
    if validator is None:
        validator = lambda p: validate_directory_markers(
            p, label=label, markers=markers
        )

    if override is not None:
        return validator(override)

    if canonical.is_dir():
        try:
            return validator(canonical)
        except Exception:
            pass

    candidates = candidate_directories(root, markers=markers)

    valid: List[Path] = []
    for candidate in candidates:
        try:
            valid.append(validator(candidate))
        except Exception:
            continue

    if len(valid) == 1:
        print(f"[input] {label}:\n  {valid[0]}")
        return valid[0]

    if interactive:
        return prompt_dir_choice(label, valid, validator)

    if not valid:
        raise FileNotFoundError(
            f"Could not resolve {label}. Canonical path:\n  {canonical}"
        )

    raise RuntimeError(
        f"Multiple candidates found for {label}; provide an explicit override:\n"
        + "\n".join(f"  - {path}" for path in valid)
    )


def resolve_file(
    *,
    root: Path,
    canonical: Path,
    basename: str,
    label: str,
    override: Optional[Path],
    interactive: bool,
) -> Path:
    if override is not None:
        return validate_file(override, label)

    if canonical.is_file():
        return validate_file(canonical, label)

    candidates = file_candidates(root, basename)

    if len(candidates) == 1:
        print(f"[input] {label}:\n  {candidates[0]}")
        return candidates[0]

    if interactive:
        return prompt_file_choice(label, candidates)

    if not candidates:
        raise FileNotFoundError(
            f"Could not resolve {label}. Canonical path:\n  {canonical}"
        )

    raise RuntimeError(
        f"Multiple candidates found for {label}; provide an explicit override:\n"
        + "\n".join(f"  - {path}" for path in candidates)
    )


def resolve_inputs(
    *,
    selected_steps: Sequence[int],
    longitudinal_root: Path,
    pseudo_root: Path,
    long_step9_override: Optional[Path],
    pseudo_step7_override: Optional[Path],
    long_step11_override: Optional[Path],
    pseudo_step8_override: Optional[Path],
    pseudo_trajectories_override: Optional[Path],
    interactive: bool,
) -> Dict[str, Path]:
    inputs: Dict[str, Path] = {}

    if 10 in selected_steps:
        inputs["longitudinal_step9_dir"] = resolve_marker_dir(
            root=longitudinal_root,
            canonical=longitudinal_root / LONG_STEP9_DIR,
            label="Longitudinal Step-9 forward-drift directory",
            markers=STEP9_MARKERS,
            override=long_step9_override,
            interactive=interactive,
        )

        inputs["pseudo_step7_dir"] = resolve_marker_dir(
            root=pseudo_root,
            canonical=pseudo_root / PSEUDO_STEP7_DIR,
            label="Pseudo Step-7 conditioning/ordering benchmark directory",
            markers=STEP7_MARKERS,
            override=pseudo_step7_override,
            interactive=interactive,
        )

    if 12 in selected_steps:
        inputs["longitudinal_step11_dir"] = resolve_marker_dir(
            root=longitudinal_root,
            canonical=longitudinal_root / LONG_STEP11_DIR,
            label="Longitudinal Step-11 fluctuation-dynamics directory",
            markers=STEP11_MARKERS,
            override=long_step11_override,
            interactive=interactive,
            validator=validate_step11_dir,
        )

        inputs["pseudo_step8_dir"] = resolve_marker_dir(
            root=pseudo_root,
            canonical=pseudo_root / PSEUDO_STEP8_DIR,
            label="Pseudo Step-8 representation-assessment directory",
            markers=STEP8_MARKERS,
            override=pseudo_step8_override,
            interactive=interactive,
            validator=validate_step8_dir,
        )

        inputs["pseudo_trajectories"] = resolve_file(
            root=pseudo_root,
            canonical=(
                pseudo_root
                / PSEUDO_STEP4_DIR
                / TRAJECTORY_BASENAME
            ),
            basename=TRAJECTORY_BASENAME,
            label="Pseudo Step-4 latent_trajectories_long.parquet",
            override=pseudo_trajectories_override,
            interactive=interactive,
        )

    return inputs


# =============================================================================
# Completion / fresh / resume
# =============================================================================


def required_completion(
    comparison_root: Path,
    step: int,
) -> List[Path]:
    base = step_outdir(comparison_root, step)
    rels = STEP10_COMPLETION if step == 10 else STEP12_COMPLETION
    return [base / rel for rel in rels]


def step_complete(
    comparison_root: Path,
    step: int,
) -> bool:
    for path in required_completion(comparison_root, step):
        if not path.is_file() or path.stat().st_size <= 0:
            return False
    return True


def selected_outputs_exist(
    comparison_root: Path,
    selected_steps: Sequence[int],
) -> bool:
    for step in selected_steps:
        path = step_outdir(comparison_root, step)
        if path.is_dir() and any(path.iterdir()):
            return True
    return False


def clear_selected_outputs(
    comparison_root: Path,
    selected_steps: Sequence[int],
) -> None:
    for step in selected_steps:
        path = step_outdir(comparison_root, step)
        if path.exists():
            shutil.rmtree(path)

        log = comparison_root / "logs" / f"step{step}.log"
        if log.exists():
            log.unlink()


def choose_run_mode(
    *,
    comparison_root: Path,
    selected_steps: Sequence[int],
    restart: bool,
    resume: bool,
    interactive: bool,
) -> str:
    if restart and resume:
        raise ValueError("Use only one of --restart or --resume.")

    if not selected_outputs_exist(comparison_root, selected_steps):
        return "fresh"

    if restart:
        return "fresh"

    if resume:
        return "resume"

    if not interactive:
        raise RuntimeError(
            "Selected comparison outputs already exist. "
            "Use --restart or --resume."
        )

    print(f"\nExisting outputs detected for {workflow_label(selected_steps)}:")
    print("  [1] fresh   — delete only selected comparison outputs")
    print("  [2] resume  — skip complete comparisons / reuse Step-12 sufficient stats")
    print("  [3] abort")

    while True:
        raw = input("Selection [2]: ").strip().lower()

        if raw in {"", "2", "resume", "r"}:
            return "resume"

        if raw in {"1", "fresh", "f"}:
            return "fresh"

        if raw in {"3", "abort", "a", "q", "quit"}:
            return "abort"

        print("Please enter 1, 2, or 3.")


# =============================================================================
# Command construction
# =============================================================================


def build_step10_command(
    *,
    script: Path,
    comparison_root: Path,
    inputs: Dict[str, Path],
    grid_points: int,
) -> List[str]:
    return [
        str(Path(sys.executable).resolve()),
        "-u",
        str(script),
        "--longitudinal-dir",
        str(inputs["longitudinal_step9_dir"]),
        "--pseudo-benchmark-dir",
        str(inputs["pseudo_step7_dir"]),
        "--outdir",
        str(step_outdir(comparison_root, 10)),
        "--representation",
        "cross_combined",
        "--n-grid-points",
        str(int(grid_points)),
    ]


def step12_resume_sufficient(
    comparison_root: Path,
    mode: str,
) -> Optional[Path]:
    if mode != "resume":
        return None

    path = (
        step_outdir(comparison_root, 12)
        / STEP12_SUFFICIENT_PARQUET
    )

    if path.is_file() and path.stat().st_size > 0:
        return path.resolve()

    return None


def build_step12_command(
    *,
    script: Path,
    comparison_root: Path,
    inputs: Dict[str, Path],
    mode: str,
    n_configurations: int,
    seed: int,
    pseudo_min_n: int,
    pseudo_min_subjects: int,
    min_config_fraction: float,
    min_grid_fraction: float,
    grid_points: int,
) -> Tuple[List[str], str]:
    sufficient = step12_resume_sufficient(comparison_root, mode)

    command = [
        str(Path(sys.executable).resolve()),
        "-u",
        str(script),
        "--longitudinal-dir",
        str(inputs["longitudinal_step11_dir"]),
        "--pseudo-representation-dir",
        str(inputs["pseudo_step8_dir"]),
        "--outdir",
        str(step_outdir(comparison_root, 12)),
    ]

    if sufficient is not None:
        command.extend(
            [
                "--pseudo-sufficient",
                str(sufficient),
            ]
        )
        pseudo_source = "existing_step12_sufficient"
    else:
        command.extend(
            [
                "--pseudo-trajectories",
                str(inputs["pseudo_trajectories"]),
            ]
        )
        pseudo_source = "pseudo_step4_trajectories"

    command.extend(
        [
            "--longitudinal-dt",
            "1",
            "--n-configurations",
            str(int(n_configurations)),
            "--seed",
            str(int(seed)),
            "--pseudo-min-n",
            str(int(pseudo_min_n)),
            "--pseudo-min-subjects",
            str(int(pseudo_min_subjects)),
            "--min-config-fraction",
            repr(float(min_config_fraction)),
            "--min-grid-fraction",
            repr(float(min_grid_fraction)),
            "--grid-points",
            str(int(grid_points)),
        ]
    )

    return command, pseudo_source


# =============================================================================
# Logging / manifests
# =============================================================================


def run_and_tee(
    command: Sequence[str],
    *,
    log_path: Path,
    cwd: Path,
    env: Dict[str, str],
) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)

    print("\n$ " + shell_join(command))
    print(f"[log] {log_path}")

    with log_path.open(
        "w",
        encoding="utf-8",
        buffering=1,
    ) as log:
        process = subprocess.Popen(
            [str(x) for x in command],
            cwd=str(cwd),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert process.stdout is not None

        for line in process.stdout:
            print(line, end="")
            log.write(line)

        return int(process.wait())


def write_manifest(
    path: Path,
    rows: List[Dict[str, object]],
) -> None:
    fields = [
        "step",
        "script",
        "status",
        "started_at",
        "finished_at",
        "returncode",
        "log_path",
        "command",
        "message",
    ]

    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            writer.writerow(
                {field: row.get(field, "") for field in fields}
            )


# =============================================================================
# Display
# =============================================================================


def print_preflight(
    *,
    selected_steps: Sequence[int],
    longitudinal_root: Path,
    pseudo_root: Path,
    comparison_root: Path,
    inputs: Dict[str, Path],
    provenance: Dict[str, object],
    n_configurations: int,
    seed: int,
) -> None:
    print("\n" + "=" * 86)
    print(
        "ClonoDynamics longitudinal vs pseudo — "
        + workflow_label(selected_steps)
    )
    print("=" * 86)
    print(f"Longitudinal root       : {longitudinal_root}")
    print(f"Pseudo root             : {pseudo_root}")
    print(f"Comparison results root : {comparison_root}")
    print(f"Selected comparisons    : {','.join(map(str, selected_steps))}")

    long_alpha = provenance.get("longitudinal_observability_alpha")
    pseudo_alpha = provenance.get("pseudo_observability_alpha")
    print(f"Longitudinal alpha      : {long_alpha}")
    print(f"Pseudo alpha            : {pseudo_alpha}")

    compared = provenance.get("compared_upstream_hash_keys") or []
    print(
        "Upstream hash checks    : "
        + (", ".join(compared) if compared else "not available")
    )

    print("\nResolved upstream inputs:")
    display_order = (
        "longitudinal_step9_dir",
        "pseudo_step7_dir",
        "longitudinal_step11_dir",
        "pseudo_step8_dir",
        "pseudo_trajectories",
    )
    for key in display_order:
        if key in inputs:
            print(f"  {key:26s}: {inputs[key]}")

    if 10 in selected_steps:
        print("\nStep 10:")
        print("  representation         : cross_combined")
        print("  grid points            : 101")
        print("  rank normalization     : disabled")
        print("  forced common bins     : disabled")

    if 12 in selected_steps:
        print("\nStep 12:")
        print("  longitudinal dt        : 1")
        print(f"  pseudo configurations  : {n_configurations}")
        print("  pseudo min n           : 200")
        print("  pseudo min subjects    : 3")
        print("  min config fraction    : 0.95")
        print("  min grid fraction      : 0.95")
        print("  grid points            : 101")
        print(f"  seed                   : {seed}")

    print("=" * 86)


# =============================================================================
# CLI
# =============================================================================


def parse_args(
    argv: Optional[Sequence[str]] = None,
) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Dependency-aware ClonoDynamics orchestrator for "
            "Steps 10 and 12 longitudinal-vs-pseudo comparisons."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--workflow",
        choices=["all", "10", "12"],
        default=None,
        help="Run both comparisons, Step 10 only, or Step 12 only.",
    )

    parser.add_argument(
        "--longitudinal-root",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--pseudo-root",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--comparison-root",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--code-dir",
        type=Path,
        default=None,
    )

    # Explicit upstream overrides.
    parser.add_argument("--step9-dir", type=Path, default=None)
    parser.add_argument("--step7-dir", type=Path, default=None)
    parser.add_argument("--step11-dir", type=Path, default=None)
    parser.add_argument("--step8-dir", type=Path, default=None)
    parser.add_argument("--pseudo-trajectories", type=Path, default=None)

    parser.add_argument(
        "--n-configurations",
        type=int,
        default=2000,
        help="Step-12 pseudo-time ordering configurations.",
    )
    parser.add_argument("--seed", type=int, default=123)

    parser.add_argument(
        "--allow-provenance-mismatch",
        action="store_true",
        help=(
            "Allow comparison despite detectable alpha/upstream-script mismatch. "
            "Use only intentionally."
        ),
    )

    parser.add_argument(
        "--restart",
        action="store_true",
        help="Delete selected comparison outputs and rerun.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Skip completed comparisons; Step 12 reuses its existing sufficient "
            "table when available."
        ),
    )

    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--yes", action="store_true")
    parser.add_argument("--non-interactive", action="store_true")

    return parser.parse_args(argv)


# =============================================================================
# Main
# =============================================================================


def main(
    argv: Optional[Sequence[str]] = None,
) -> int:
    print(
        "ClonoDynamics longitudinal-vs-pseudo orchestrator | "
        f"{VERSION}"
    )

    args = parse_args(argv)
    interactive = not bool(args.non_interactive)

    # Workflow.
    if args.workflow is not None:
        selected_steps = normalize_workflow(args.workflow)
    elif interactive:
        selected_steps = prompt_workflow()
    else:
        raise ValueError("--workflow is required in non-interactive mode.")

    # Longitudinal root.
    if args.longitudinal_root is not None:
        longitudinal_root = normalize_longitudinal_root(
            args.longitudinal_root
        )
    elif interactive:
        longitudinal_root = prompt_existing_root(
            "Existing longitudinal results root",
            normalize_longitudinal_root,
        )
    else:
        raise ValueError(
            "--longitudinal-root is required in non-interactive mode."
        )

    # Pseudo root.
    if args.pseudo_root is not None:
        pseudo_root = normalize_pseudo_root(args.pseudo_root)
    elif interactive:
        pseudo_root = prompt_existing_root(
            "Existing pseudo-longitudinal results root",
            normalize_pseudo_root,
        )
    else:
        raise ValueError(
            "--pseudo-root is required in non-interactive mode."
        )

    if longitudinal_root == pseudo_root:
        raise ValueError(
            "Longitudinal and pseudo result roots must be different."
        )

    # Upstream provenance.
    long_cfg = read_pipeline_config(longitudinal_root)
    pseudo_cfg = read_pipeline_config(pseudo_root)

    provenance = validate_cross_dataset_provenance(
        longitudinal_root=longitudinal_root,
        pseudo_root=pseudo_root,
        long_cfg=long_cfg,
        pseudo_cfg=pseudo_cfg,
        allow_mismatch=bool(args.allow_provenance_mismatch),
    )

    # Comparison result root.
    proposed_comparison_root = default_comparison_root(
        longitudinal_root,
        pseudo_root,
    )

    if args.comparison_root is not None:
        comparison_root = args.comparison_root.expanduser().resolve()
    elif interactive:
        comparison_root = prompt_comparison_root(
            proposed_comparison_root
        )
    else:
        comparison_root = proposed_comparison_root

    # Never allow comparison writes inside either upstream root.
    for upstream_root, label in (
        (longitudinal_root, "longitudinal"),
        (pseudo_root, "pseudo"),
    ):
        try:
            comparison_root.relative_to(upstream_root)
        except ValueError:
            pass
        else:
            raise ValueError(
                "Comparison results root must not be inside the "
                f"{label} upstream result tree:\n  {comparison_root}"
            )

    # Code/scripts.
    code_dir = resolve_code_dir(args.code_dir)
    scripts = validate_scripts(code_dir, selected_steps)

    # Upstream analysis artifacts.
    inputs = resolve_inputs(
        selected_steps=selected_steps,
        longitudinal_root=longitudinal_root,
        pseudo_root=pseudo_root,
        long_step9_override=args.step9_dir,
        pseudo_step7_override=args.step7_dir,
        long_step11_override=args.step11_dir,
        pseudo_step8_override=args.step8_dir,
        pseudo_trajectories_override=args.pseudo_trajectories,
        interactive=interactive,
    )

    if int(args.n_configurations) < 1:
        raise ValueError("--n-configurations must be >= 1")

    print_preflight(
        selected_steps=selected_steps,
        longitudinal_root=longitudinal_root,
        pseudo_root=pseudo_root,
        comparison_root=comparison_root,
        inputs=inputs,
        provenance=provenance,
        n_configurations=int(args.n_configurations),
        seed=int(args.seed),
    )

    # Fresh/resume.
    mode = choose_run_mode(
        comparison_root=comparison_root,
        selected_steps=selected_steps,
        restart=bool(args.restart),
        resume=bool(args.resume),
        interactive=interactive,
    )

    if mode == "abort":
        print("Aborted.")
        return 0

    if mode == "fresh":
        clear_selected_outputs(
            comparison_root,
            selected_steps,
        )

    comparison_root.mkdir(parents=True, exist_ok=True)
    (comparison_root / "logs").mkdir(
        parents=True,
        exist_ok=True,
    )

    run_dir = (
        comparison_root
        / "orchestrator_runs_longitudinal_vs_pseudo"
    )
    run_dir.mkdir(parents=True, exist_ok=True)

    # Commands.
    commands: Dict[int, List[str]] = {}
    step12_pseudo_source: Optional[str] = None

    if 10 in selected_steps:
        commands[10] = build_step10_command(
            script=scripts[10],
            comparison_root=comparison_root,
            inputs=inputs,
            grid_points=101,
        )

    if 12 in selected_steps:
        command12, step12_pseudo_source = build_step12_command(
            script=scripts[12],
            comparison_root=comparison_root,
            inputs=inputs,
            mode=mode,
            n_configurations=int(args.n_configurations),
            seed=int(args.seed),
            pseudo_min_n=200,
            pseudo_min_subjects=3,
            min_config_fraction=0.95,
            min_grid_fraction=0.95,
            grid_points=101,
        )
        commands[12] = command12

    # Reproducibility config.
    config = {
        "orchestrator": Path(__file__).name,
        "orchestrator_version": VERSION,
        "generated_at": now_iso(),
        "selected_steps": list(selected_steps),
        "workflow": workflow_label(selected_steps),
        "longitudinal_root": str(longitudinal_root),
        "pseudo_root": str(pseudo_root),
        "comparison_root": str(comparison_root),
        "code_dir": str(code_dir),
        "cross_dataset_provenance": provenance,
        "resolved_inputs": {
            key: str(value)
            for key, value in inputs.items()
        },
        "script_sha256": {
            str(step): sha256_file(path)
            for step, path in scripts.items()
        },
        "step10": {
            "selected": 10 in selected_steps,
            "representation": "cross_combined",
            "n_grid_points": 101,
            "absolute_x0_support": True,
            "rank_or_percentile_normalization": False,
            "common_bin_edges_forced": False,
        },
        "step12": {
            "selected": 12 in selected_steps,
            "longitudinal_dt": 1,
            "n_configurations": int(args.n_configurations),
            "seed": int(args.seed),
            "pseudo_min_n": 200,
            "pseudo_min_subjects": 3,
            "min_config_fraction": 0.95,
            "min_grid_fraction": 0.95,
            "grid_points": 101,
            "pseudo_pair_source": step12_pseudo_source,
            "primary_metric": "sample Var(dx_latent | xstar_latent), TT",
            "complementary_metric": "MSD",
            "mad_role": "descriptive fixed pseudo reference",
            "temporal_scaling_rule": (
                "pseudo baseline is not subtracted before Step 13 fitting"
            ),
        },
        "python_executable": str(Path(sys.executable).resolve()),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
    }

    config_text = json.dumps(
        config,
        indent=2,
        ensure_ascii=False,
    )

    config_path = (
        comparison_root
        / "00_orchestrator_longitudinal_vs_pseudo_config.json"
    )
    config_path.write_text(
        config_text,
        encoding="utf-8",
    )

    run_config_path = (
        run_dir
        / (
            f"{timestamp_slug()}_"
            f"{'-'.join(map(str, selected_steps))}_config.json"
        )
    )
    run_config_path.write_text(
        config_text,
        encoding="utf-8",
    )

    # Plan.
    print("\nExecution plan:")

    for step in selected_steps:
        state = (
            "SKIP (complete)"
            if mode == "resume"
            and step_complete(comparison_root, step)
            else "RUN"
        )

        print(f"\nStep {step}: {state}")
        print("  " + shell_join(commands[step]))

        if step == 12 and state == "RUN":
            if step12_pseudo_source == "existing_step12_sufficient":
                print(
                    "  [resume] Reusing existing Step-12 "
                    "02_pseudo_directed_pair_sufficient.parquet"
                )
            else:
                print(
                    "  [source] Building pseudo directed-pair sufficient "
                    "statistics from Step-4 trajectories"
                )

    if args.dry_run:
        print(
            "\n[DRY RUN] Cross-dataset dependencies resolved and "
            "commands built. Nothing executed."
        )
        return 0

    if not args.yes:
        answer = input(
            f"\nRun {workflow_label(selected_steps)} now? [Y/n]: "
        ).strip().lower()

        if answer not in {
            "",
            "y",
            "yes",
            "s",
            "si",
            "sì",
        }:
            print("Aborted.")
            return 0

    # Environment.
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    old_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        str(code_dir)
        if not old_pythonpath
        else str(code_dir)
        + os.pathsep
        + old_pythonpath
    )

    manifest_path = (
        comparison_root
        / "00_orchestrator_longitudinal_vs_pseudo_manifest.csv"
    )

    run_manifest_path = (
        run_dir
        / (
            f"{timestamp_slug()}_"
            f"{'-'.join(map(str, selected_steps))}_manifest.csv"
        )
    )

    rows: List[Dict[str, object]] = []

    # Execution.
    for step in selected_steps:
        script_name = (
            STEP10_SCRIPT
            if step == 10
            else STEP12_SCRIPT
        )

        log_path = (
            comparison_root
            / "logs"
            / f"step{step}.log"
        )

        if (
            mode == "resume"
            and step_complete(comparison_root, step)
        ):
            row = {
                "step": step,
                "script": script_name,
                "status": "skipped_complete",
                "started_at": now_iso(),
                "finished_at": now_iso(),
                "returncode": 0,
                "log_path": str(log_path),
                "command": shell_join(commands[step]),
                "message": "Required completion markers already present.",
            }

            rows.append(row)
            write_manifest(manifest_path, rows)
            write_manifest(run_manifest_path, rows)

            print(f"\n[SKIP] Step {step}: complete.")
            continue

        step_outdir(comparison_root, step).mkdir(
            parents=True,
            exist_ok=True,
        )

        started = now_iso()

        print("\n" + "=" * 86)
        print(f"STEP {step} — {script_name}")
        print("=" * 86)

        returncode = run_and_tee(
            commands[step],
            log_path=log_path,
            cwd=code_dir,
            env=env,
        )

        finished = now_iso()
        status = (
            "completed"
            if returncode == 0
            else "failed"
        )
        message = ""

        if (
            returncode == 0
            and not step_complete(comparison_root, step)
        ):
            status = "failed_output_check"

            missing = [
                str(path)
                for path in required_completion(
                    comparison_root,
                    step,
                )
                if not path.is_file()
                or path.stat().st_size <= 0
            ]

            message = (
                "Process exited with code 0 but required completion "
                "markers are missing/empty: "
                + "; ".join(missing)
            )

        row = {
            "step": step,
            "script": script_name,
            "status": status,
            "started_at": started,
            "finished_at": finished,
            "returncode": returncode,
            "log_path": str(log_path),
            "command": shell_join(commands[step]),
            "message": message,
        }

        rows.append(row)
        write_manifest(manifest_path, rows)
        write_manifest(run_manifest_path, rows)

        if status != "completed":
            print(f"\n[FAIL] Step {step} stopped the workflow.")
            print(f"See log:\n  {log_path}")

            if message:
                print(message)

            return (
                returncode
                if returncode != 0
                else 2
            )

        print(
            f"[OK] Step {step} completed and required outputs "
            "were verified."
        )

    print("\n" + "=" * 86)
    print(
        "ClonoDynamics longitudinal-vs-pseudo comparisons "
        "completed successfully."
    )
    print(f"Workflow : {workflow_label(selected_steps)}")
    print(f"Results  : {comparison_root}")
    print(f"Config   : {config_path}")
    print(f"Manifest : {manifest_path}")
    print("=" * 86)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
