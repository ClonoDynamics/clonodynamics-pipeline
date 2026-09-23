#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00-generate_pseudo_design.py
============================

Generate the canonical ClonoDynamics pseudo-longitudinal 1x12 randomization
manifests consumed by the current pseudo-reference builders.

BIOLOGICAL DESIGN
-----------------
Input consists of exactly 12 independently processed TCRbeta RepSeq
measurements derived from ONE biological source sample.

The production design contains 2,000 randomized configurations. Each
configuration contains:

    6 pseudo-timepoints x 2 technical replicates

and uses every one of the 12 source measurements exactly once.

RANDOMIZATION GEOMETRY
----------------------
For every configuration:

    1. randomly permute the 12 source measurements;
    2. split the permutation into six consecutive pairs;
    3. canonicalize the within-pair A/B orientation by source ID;
    4. preserve the randomized ORDER of the six pairs as pseudo-times 1..6.

Within-pair A/B swaps do not define scientifically distinct configurations in
the current pseudo-reference contract because the downstream forward benchmark
uses reciprocal AB/BA orientations symmetrically and the fluctuation benchmark
is symmetric in replicate labels.

The number of distinct configurations under this equivalence is:

    12! / 2^6 = 7,484,400.

The historical production randomization used by ClonoDynamics is preserved:

    n_configurations = 2000
    seed             = 12345

INPUT
-----
The input directory must contain exactly 12 supported repertoire tables. Each
must contain at least:

    aaSeqCDR3
    readCount

Supported formats:

    .tsv .txt .csv .parquet .pq .feather .arrow

OUTPUT
------
<outdir>/
    source_measurements.tsv
    pseudo_configurations.tsv
    pseudo_configurations_wide.tsv
    00_pseudo_1x12_config.json

By default only the manifests are written. Optional canonical raw symlinks can
be materialized with --materialize-symlinks, but they are NOT required by the
current multirepresentation pseudo builder, which reads the manifests directly.

Typical production run:

    python3 00-generate_pseudo_design.py /path/to/12_repseq_tables \
        --outdir /path/to/pseudo_reference_root \
        --n-configs 2000 \
        --seed 12345
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import shutil
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


SCRIPT_VERSION = "v1-canonical-design-manifests-2026-09-22"
MAX_CONFIGURATIONS = 7_484_400

TEXT_EXTS = {".tsv", ".txt", ".csv"}
TABLE_EXTS = TEXT_EXTS | {".parquet", ".pq", ".feather", ".arrow"}
REQUIRED_COLUMNS = {"aaSeqCDR3", "readCount"}


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            h.update(block)
    return h.hexdigest()


def discover_files(input_dir: Path) -> List[Path]:
    """Discover exactly the supported repertoire tables, in historical name order."""
    files = [
        path.resolve()
        for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in TABLE_EXTS
    ]
    # Preserve the historical generator convention exactly: case-insensitive
    # lexicographic filename order, not natural-numeric sorting.
    return sorted(files, key=lambda path: path.name.lower())


def columns_of(path: Path) -> List[str]:
    suffix = path.suffix.lower()

    if suffix in {".tsv", ".txt"}:
        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            return handle.readline().rstrip("\r\n").split("\t")

    if suffix == ".csv":
        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            first_line = handle.readline().rstrip("\r\n")
        return next(csv.reader([first_line]))

    try:
        import pandas as pd
    except ImportError as exc:
        raise RuntimeError(
            f"{path.name}: pandas/pyarrow is required to inspect {suffix}."
        ) from exc

    if suffix in {".parquet", ".pq"}:
        try:
            import pyarrow.parquet as pq
            return list(pq.ParquetFile(path).schema.names)
        except Exception:
            return list(pd.read_parquet(path).columns)

    return list(pd.read_feather(path).columns)


def safe_symlink(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    destination.symlink_to(source.resolve())


def pair_key(a: int, b: int) -> Tuple[int, int]:
    return (a, b) if a < b else (b, a)


def configuration_key(
    permutation: Sequence[int],
) -> Tuple[Tuple[int, int], ...]:
    """
    Canonical configuration signature.

    Pair order is retained because it defines pseudo-time. Within-pair A/B
    orientation is removed by sorting each adjacent pair.
    """
    if len(permutation) != 12:
        raise ValueError("A pseudo-1x12 permutation must contain exactly 12 entries")
    return tuple(
        pair_key(int(permutation[j]), int(permutation[j + 1]))
        for j in range(0, 12, 2)
    )


def generate_unique_configurations(
    n_configs: int,
    seed: int,
) -> Tuple[List[Tuple[Tuple[int, int], ...]], int]:
    if n_configs < 1:
        raise ValueError("--n-configs must be >= 1")
    if n_configs > MAX_CONFIGURATIONS:
        raise ValueError(
            f"Requested {n_configs:,} configurations but only "
            f"{MAX_CONFIGURATIONS:,} are distinct after within-pair A/B equivalence."
        )

    rng = random.Random(int(seed))
    universe = list(range(12))
    seen = set()
    configurations: List[Tuple[Tuple[int, int], ...]] = []

    attempts = 0
    max_attempts = max(10_000, 100 * int(n_configs))

    while len(configurations) < int(n_configs):
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError(
                f"Could not generate {n_configs} unique configurations after "
                f"{attempts} attempts."
            )

        permutation = universe[:]
        rng.shuffle(permutation)
        key = configuration_key(permutation)
        if key in seen:
            continue
        seen.add(key)
        configurations.append(key)

    return configurations, attempts


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the canonical one-subject x 12-measurement pseudo-longitudinal "
            "randomization design manifests."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing exactly 12 repertoire files.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=(
            "Destination pseudo root. Default: sibling directory named "
            "<input_dir>_pseudo_1x12."
        ),
    )
    parser.add_argument("--n-configs", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--materialize-symlinks",
        action="store_true",
        help=(
            "Also create configurations/Cxxxxxx/raw/ canonical symlinks. "
            "Not required by the current pseudo builder."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing non-empty output directory.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    input_dir = args.input_dir.expanduser().resolve(strict=True)
    if not input_dir.is_dir():
        raise NotADirectoryError(input_dir)

    outdir = (
        args.outdir.expanduser().resolve()
        if args.outdir is not None
        else input_dir.parent / f"{input_dir.name}_pseudo_1x12"
    )

    if outdir.exists() and any(outdir.iterdir()):
        if not args.force:
            raise FileExistsError(
                f"Output directory is not empty: {outdir}\n"
                "Use --force only when you intentionally want to rebuild the design."
            )
        shutil.rmtree(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    files = discover_files(input_dir)
    if len(files) != 12:
        raise ValueError(
            f"Exactly 12 supported repertoire files are required; "
            f"found {len(files)} in {input_dir}."
        )

    source_rows: List[Dict[str, object]] = []
    for index, path in enumerate(files, start=1):
        columns = set(columns_of(path))
        missing = REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(
                f"{path.name}: missing required columns {sorted(missing)}"
            )

        source_rows.append(
            {
                "source_id": f"R{index:02d}",
                "source_index": index,
                "source_file": path.name,
                "source_path": str(path),
                "suffix": path.suffix.lower(),
                "sha256": sha256_file(path),
                "has_readFraction": int("readFraction" in columns),
            }
        )

    source_manifest = outdir / "source_measurements.tsv"
    with source_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(source_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(source_rows)

    configurations, attempts = generate_unique_configurations(
        int(args.n_configs),
        int(args.seed),
    )

    long_rows: List[Dict[str, object]] = []
    wide_rows: List[Dict[str, object]] = []

    for config_index, configuration in enumerate(configurations, start=1):
        cid = f"C{config_index:06d}"
        wide: Dict[str, object] = {
            "configuration_id": cid,
            "configuration_index": config_index,
            "source_subject": 1,
            "randomization_seed": int(args.seed),
        }

        signature_parts: List[str] = []
        used: List[str] = []
        raw_dir = outdir / "configurations" / cid / "raw"

        for pseudo_time, (idx_a, idx_b) in enumerate(configuration, start=1):
            pair_source_rows = [
                source_rows[int(idx_a)],
                source_rows[int(idx_b)],
            ]

            pair_ids = [str(row["source_id"]) for row in pair_source_rows]
            signature_parts.append(f"{pair_ids[0]}-{pair_ids[1]}")

            for replicate, source in enumerate(pair_source_rows, start=1):
                source_id = str(source["source_id"])
                source_file = Path(str(source["source_path"]))
                used.append(source_id)

                canonical_name = (
                    f"1_{pseudo_time}-{replicate}{source_file.suffix.lower()}"
                )

                long_rows.append(
                    {
                        "configuration_id": cid,
                        "configuration_index": config_index,
                        "source_subject": 1,
                        "pseudo_time": pseudo_time,
                        "replicate": replicate,
                        "source_id": source_id,
                        "source_index": int(source["source_index"]),
                        "source_file": str(source["source_file"]),
                        "source_path": str(source["source_path"]),
                        "canonical_filename": canonical_name,
                        "randomization_seed": int(args.seed),
                        "within_pair_orientation_policy": "canonical_source_id",
                    }
                )

                wide[f"T{pseudo_time}_R{replicate}"] = source_id

                if args.materialize_symlinks:
                    safe_symlink(
                        source_file,
                        raw_dir / canonical_name,
                    )

        if len(used) != 12 or len(set(used)) != 12:
            raise RuntimeError(
                f"{cid}: source-use invariant failed; values={used}"
            )

        signature = "|".join(signature_parts)
        wide["configuration_signature"] = signature
        wide["configuration_sha256"] = hashlib.sha256(
            signature.encode("utf-8")
        ).hexdigest()
        wide_rows.append(wide)

    signatures = [str(row["configuration_signature"]) for row in wide_rows]
    if len(set(signatures)) != len(signatures):
        raise RuntimeError("Duplicate canonical configuration signatures detected")

    long_fields = [
        "configuration_id",
        "configuration_index",
        "source_subject",
        "pseudo_time",
        "replicate",
        "source_id",
        "source_index",
        "source_file",
        "source_path",
        "canonical_filename",
        "randomization_seed",
        "within_pair_orientation_policy",
    ]

    configuration_manifest = outdir / "pseudo_configurations.tsv"
    with configuration_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=long_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(long_rows)

    wide_fields = [
        "configuration_id",
        "configuration_index",
        "source_subject",
        "randomization_seed",
    ] + [
        f"T{time}_R{replicate}"
        for time in range(1, 7)
        for replicate in (1, 2)
    ] + [
        "configuration_signature",
        "configuration_sha256",
    ]

    wide_manifest = outdir / "pseudo_configurations_wide.tsv"
    with wide_manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=wide_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(wide_rows)

    config = {
        "script": Path(__file__).name,
        "script_version": SCRIPT_VERSION,
        "design": (
            "1 biological source sample x 12 independent technical RepSeq measurements"
        ),
        "source_subject_count": 1,
        "n_source_measurements": 12,
        "n_configurations": int(args.n_configs),
        "pseudo_timepoints_per_configuration": 6,
        "technical_replicates_per_pseudo_timepoint": 2,
        "pairing_randomized": True,
        "pseudo_time_order_randomized": True,
        "within_pair_AB_orientation_randomized": False,
        "within_pair_AB_orientation_policy": (
            "canonicalized by source ID because downstream AB/BA and covariance "
            "estimands are symmetric under within-pair label exchange"
        ),
        "configuration_equivalence": (
            "pairing plus ordered pseudo-time sequence; within-pair A/B swaps "
            "do not define new configurations"
        ),
        "max_unique_configurations": MAX_CONFIGURATIONS,
        "generation_attempts": int(attempts),
        "seed": int(args.seed),
        "configuration_is_biological_subject": False,
        "input_dir": str(input_dir),
        "output_dir": str(outdir),
        "raw_symlinks_materialized": bool(args.materialize_symlinks),
        "source_file_ordering": "case-insensitive lexicographic filename order",
        "source_manifest": source_manifest.name,
        "configuration_manifest": configuration_manifest.name,
        "configuration_wide_manifest": wide_manifest.name,
    }

    config_path = outdir / "00_pseudo_1x12_config.json"
    config_path.write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("\n[DONE] canonical pseudo 1x12 design")
    print("Biological source samples        : 1")
    print("Technical source measurements    : 12")
    print("Configurations                   :", len(configurations))
    print("Pseudo-times per configuration   : 6")
    print("Replicates per pseudo-time       : 2")
    print("Within-pair A/B orientation      : canonicalized")
    print("configuration_id as biological N : NO")
    print("Random seed                      :", int(args.seed))
    print("Generation attempts              :", int(attempts))
    print("Raw symlinks materialized        :", "YES" if args.materialize_symlinks else "NO")
    print("Output                           :", outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
