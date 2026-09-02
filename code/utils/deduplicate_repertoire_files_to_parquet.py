#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse
import pandas as pd


def read_table(path: Path, sep: str):
    if path.suffix.lower() in [".parquet", ".pq"]:
        return pd.read_parquet(path)
    if path.suffix.lower() in [".tsv", ".txt"]:
        return pd.read_csv(path, sep=sep)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported file type: {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--pattern", default="*")
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted([
        p for p in in_dir.glob(args.pattern)
        if p.suffix.lower() in [".tsv", ".txt", ".csv", ".parquet", ".pq"]
    ])

    if not files:
        raise SystemExit(f"No input files found in {in_dir}")

    report = []

    for p in files:
        df = read_table(p, args.sep)

        if "aaSeqCDR3" not in df.columns:
            raise ValueError(f"{p} missing aaSeqCDR3")
        if "readCount" not in df.columns:
            raise ValueError(f"{p} missing readCount")

        before = len(df)
        unique_before = df["aaSeqCDR3"].nunique()

        df["readCount"] = pd.to_numeric(df["readCount"], errors="coerce").fillna(0)

        out = (
            df.groupby("aaSeqCDR3", as_index=False)
              .agg(readCount=("readCount", "sum"))
        )

        out = out[out["readCount"] > 0].copy()
        out["readFraction"] = out["readCount"] / out["readCount"].sum()

        out_path = out_dir / f"{p.stem}.parquet"
        out.to_parquet(out_path, index=False)

        report.append({
            "file": p.name,
            "out_file": out_path.name,
            "rows_before": before,
            "unique_before": unique_before,
            "duplicates_removed": before - unique_before,
            "rows_after": len(out),
            "total_reads_after": int(out["readCount"].sum()),
        })

        print(f"[OK] {p.name} -> {out_path.name} | {before} -> {len(out)}")

    rep = pd.DataFrame(report)
    rep.to_csv(out_dir / "deduplication_report.csv", index=False)
    print(f"\n[DONE] Report: {out_dir / 'deduplication_report.csv'}")


if __name__ == "__main__":
    main()
