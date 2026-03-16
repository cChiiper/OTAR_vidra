#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path


def create_spark_session():
    from pyspark.sql import SparkSession

    return (
        SparkSession.builder
        .appName("headtext")
        .master("local[1]")
        .getOrCreate()
    )


def read_head(input_path: Path, n_lines: int, spark=None) -> str:
    """Return the first ``n_lines`` rows from a Parquet file as JSON lines."""

    if n_lines < 0:
        raise ValueError("n_lines must be >= 0")

    should_stop_spark = spark is None
    spark = spark or create_spark_session()

    try:
        rows = (
            spark.read
            .parquet(str(input_path))
            .limit(n_lines)
            .toJSON()
            .collect()
        )
    finally:
        if should_stop_spark:
            spark.stop()

    return "\n".join(rows) + ("\n" if rows else "")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Write the first N rows of a Parquet file to a text file as JSON lines.",
    )
    parser.add_argument("--input", required=True, type=Path, help="Input Parquet file")
    parser.add_argument("--output", required=True, type=Path, help="Output text file")
    parser.add_argument(
        "--n_lines",
        type=int,
        default=5,
        help="Number of rows to keep from the start of the dataset",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    text = read_head(args.input, args.n_lines)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
