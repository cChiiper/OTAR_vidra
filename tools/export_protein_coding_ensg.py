#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export distinct protein-coding ENSG IDs from Open Targets target parquet files."
    )
    parser.add_argument(
        "--target-parquet",
        required=True,
        help="Parquet glob or dataset path for Open Targets target data.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV path. The file is written without a header.",
    )
    return parser.parse_args(argv)


def protein_coding_ids(target_parquet: str) -> list[str]:
    try:
        import duckdb
    except ImportError as exc:
        raise RuntimeError("duckdb is required; run `make venv` first.") from exc

    con = duckdb.connect(":memory:")
    try:
        rows = con.execute(
            """
            SELECT DISTINCT id
            FROM read_parquet(?)
            WHERE biotype = 'protein_coding'
              AND id IS NOT NULL
            ORDER BY id
            """,
            [target_parquet],
        ).fetchall()
    finally:
        con.close()

    return [row[0] for row in rows if row[0]]


def write_gene_csv(gene_ids: list[str], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerows((gene_id,) for gene_id in gene_ids)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    output = Path(args.output)
    gene_ids = protein_coding_ids(args.target_parquet)
    if not gene_ids:
        raise RuntimeError(f"No protein-coding ENSG IDs found in {args.target_parquet!r}.")

    write_gene_csv(gene_ids, output)
    print(f"Wrote {len(gene_ids)} protein-coding ENSG IDs to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
