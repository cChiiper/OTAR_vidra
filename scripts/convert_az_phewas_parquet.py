#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import lzma
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

DEFAULT_CHUNK_SIZE = 250_000
DEFAULT_TARGET_FILE_SIZE_MB = 512.0

BINARY_DTYPES = {
    "Gene": "string",
    "Phenotype": "string",
    "CollapsingModel": "string",
    "Type": "string",
    "pValue": "float64",
    "Category": "string",
    "nSamples": "Int64",
    "BinNcases": "Int64",
    "BinQVcases": "Int64",
    "BinNcontrols": "Int64",
    "BinQVcontrols": "Int64",
    "BinCaseFreq": "float64",
    "BinCtrlFreq": "float64",
    "BinOddsRatio": "float64",
    "BinOddsRatioLCI": "float64",
    "BinOddsRatioUCI": "float64",
}

QUANTITATIVE_DTYPES = {
    "Gene": "string",
    "Phenotype": "string",
    "CollapsingModel": "string",
    "Type": "string",
    "pValue": "float64",
    "Category": "string",
    "nSamples": "Int64",
    "YesQV": "Int64",
    "NoQV": "Int64",
    "ConCaseMedian": "float64",
    "ConCtrMedian": "float64",
    "beta": "float64",
    "ConBetaSe": "float64",
    "LCI": "float64",
    "UCI": "float64",
}


@dataclass(frozen=True)
class DatasetConfig:
    name: str
    dtypes: dict[str, str]


@dataclass(frozen=True)
class ConversionJob:
    dataset_name: str
    input_path: Path
    output_dir: Path
    dtypes: dict[str, str]


DATASET_CONFIGS = {
    "binary": DatasetConfig(
        name="binary",
        dtypes=BINARY_DTYPES,
    ),
    "quantitative": DatasetConfig(
        name="quantitative",
        dtypes=QUANTITATIVE_DTYPES,
    ),
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the AstraZeneca 470k PheWAS binary and quantitative CSV.xz "
            "files into snappy parquet datasets. Input and output paths are "
            "supplied explicitly, typically via the Makefile."
        )
    )
    parser.add_argument(
        "--dataset",
        choices=("all", "binary", "quantitative"),
        default="all",
        help="Which dataset to convert. Default converts both datasets.",
    )
    parser.add_argument(
        "--input",
        default="",
        help="Input .csv.xz for single-dataset runs.",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="Output directory for single-dataset runs.",
    )
    parser.add_argument(
        "--binary-input",
        default="",
        help="Binary input path. Used by --dataset all or as a fallback for --dataset binary.",
    )
    parser.add_argument(
        "--binary-output-dir",
        default="",
        help="Binary output directory. Used by --dataset all or as a fallback for --dataset binary.",
    )
    parser.add_argument(
        "--quantitative-input",
        default="",
        help="Quantitative input path. Used by --dataset all or as a fallback for --dataset quantitative.",
    )
    parser.add_argument(
        "--quantitative-output-dir",
        default="",
        help="Quantitative output directory. Used by --dataset all or as a fallback for --dataset quantitative.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Rows per pandas chunk when streaming the compressed CSV.",
    )
    parser.add_argument(
        "--target-file-size-mb",
        type=float,
        default=DEFAULT_TARGET_FILE_SIZE_MB,
        help=(
            "Approximate uncompressed Arrow size per parquet part. "
            "Larger values produce fewer, larger files."
        ),
    )
    return parser


def resolve_input_path(path_like: str) -> Path:
    return Path(path_like).expanduser().resolve()


def resolve_output_dir(path_like: str) -> Path:
    output_dir = Path(path_like).expanduser()
    if not output_dir.is_absolute():
        output_dir = Path.cwd() / output_dir
    return output_dir


def require_path_argument(path_like: str, *, dataset_name: str, option_name: str, kind: str) -> str:
    text = str(path_like).strip()
    if not text:
        raise ValueError(f"{dataset_name} dataset requires {kind} via {option_name}")
    return text


def build_job(dataset_name: str, *, input_override: str = "", output_override: str = "") -> ConversionJob:
    config = DATASET_CONFIGS[dataset_name]
    input_option = "--input/--binary-input" if dataset_name == "binary" else "--input/--quantitative-input"
    output_option = "--output-dir/--binary-output-dir" if dataset_name == "binary" else "--output-dir/--quantitative-output-dir"
    input_path = resolve_input_path(
        require_path_argument(
            input_override,
            dataset_name=dataset_name,
            option_name=input_option,
            kind="an input path",
        )
    )
    output_dir = resolve_output_dir(
        require_path_argument(
            output_override,
            dataset_name=dataset_name,
            option_name=output_option,
            kind="an output directory",
        )
    )
    return ConversionJob(
        dataset_name=dataset_name,
        input_path=input_path,
        output_dir=output_dir,
        dtypes=config.dtypes,
    )


def resolve_jobs(args: argparse.Namespace) -> list[ConversionJob]:
    if args.dataset == "all":
        if args.input or args.output_dir:
            raise ValueError("--input and --output-dir can only be used with --dataset binary or quantitative")
        return [
            build_job("binary", input_override=args.binary_input, output_override=args.binary_output_dir),
            build_job(
                "quantitative",
                input_override=args.quantitative_input,
                output_override=args.quantitative_output_dir,
            ),
        ]

    input_override = args.input or getattr(args, f"{args.dataset}_input")
    output_override = args.output_dir or getattr(args, f"{args.dataset}_output_dir")
    return [build_job(args.dataset, input_override=input_override, output_override=output_override)]


def required_columns_for(dtypes: dict[str, str]) -> list[str]:
    return list(dtypes)


def read_header_columns(input_path: Path) -> list[str]:
    with lzma.open(input_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        return next(reader, [])


def ensure_required_columns(
    header_columns: list[str],
    required_columns: list[str],
    dataset_name: str,
    input_path: Path,
) -> None:
    if not header_columns:
        raise ValueError(f"{dataset_name} input has no header row: {input_path}")
    missing = [column for column in required_columns if column not in header_columns]
    if missing:
        raise ValueError(f"{dataset_name} input is missing required columns {missing}: {input_path}")


def iter_input_chunks(input_path: Path, chunk_size: int, dtypes: dict[str, str]):
    columns = required_columns_for(dtypes)
    return pd.read_csv(
        input_path,
        compression="xz",
        chunksize=chunk_size,
        dtype=dtypes,
        usecols=columns,
        keep_default_na=True,
    )


def remove_existing_path(path: Path) -> None:
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def part_filename(part_index: int) -> str:
    return f"part-{part_index:05d}.snappy.parquet"


def target_file_size_bytes(size_mb: float) -> int:
    if size_mb <= 0:
        raise ValueError(f"target_file_size_mb must be positive, got {size_mb}")
    return max(1, int(size_mb * 1024 * 1024))


def row_count_for_bytes(table: pa.Table, byte_budget: int) -> int:
    if table.num_rows <= 1:
        return table.num_rows
    bytes_per_row = max(table.nbytes / table.num_rows, 1.0)
    return max(1, min(table.num_rows, int(byte_budget / bytes_per_row)))


def write_dataset(
    input_path: Path,
    output_dir: Path,
    chunk_size: int,
    approx_target_file_size_mb: float,
    dataset_name: str,
    dtypes: dict[str, str],
) -> tuple[int, int]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")

    header_columns = read_header_columns(input_path)
    ensure_required_columns(header_columns, required_columns_for(dtypes), dataset_name, input_path)
    target_bytes = target_file_size_bytes(approx_target_file_size_mb)

    remove_existing_path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    writer: pq.ParquetWriter | None = None
    schema: pa.Schema | None = None
    current_file_bytes = 0
    part_count = 0
    total_rows = 0

    try:
        for chunk in iter_input_chunks(input_path, chunk_size, dtypes):
            if chunk.empty:
                continue

            total_rows += len(chunk)
            table = pa.Table.from_pandas(chunk, preserve_index=False)
            if schema is None:
                schema = table.schema
            elif table.schema != schema:
                table = table.cast(schema, safe=False)

            rows_written = 0
            while rows_written < table.num_rows:
                if writer is None:
                    writer = pq.ParquetWriter(
                        output_dir / part_filename(part_count),
                        schema,
                        compression="snappy",
                    )
                    current_file_bytes = 0
                    part_count += 1

                remaining = table.slice(rows_written)
                budget = max(1, target_bytes - current_file_bytes)
                rows_this_write = row_count_for_bytes(remaining, budget)
                piece = remaining.slice(0, rows_this_write)
                writer.write_table(piece)
                current_file_bytes += max(piece.nbytes, piece.num_rows)
                rows_written += piece.num_rows

                if current_file_bytes >= target_bytes:
                    writer.close()
                    writer = None
                    current_file_bytes = 0

        if total_rows == 0 or schema is None:
            raise ValueError(f"Input CSV has no data rows: {input_path}")

        if writer is not None:
            writer.close()
            writer = None
        return total_rows, part_count
    except Exception:
        shutil.rmtree(output_dir, ignore_errors=True)
        raise
    finally:
        if writer is not None:
            writer.close()


def finalize_dataset(output_dir: Path) -> None:
    (output_dir / "_SUCCESS").touch()


def move_completed_dataset(temp_output_dir: Path, output_dir: Path) -> None:
    remove_existing_path(output_dir)
    temp_output_dir.rename(output_dir)


def convert_job(job: ConversionJob, chunk_size: int, target_file_size_mb: float) -> tuple[int, int]:
    job.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_output_dir = Path(tempfile.mkdtemp(prefix=f"{job.output_dir.name}.build.", dir=job.output_dir.parent))

    try:
        total_rows, part_count = write_dataset(
            input_path=job.input_path,
            output_dir=temp_output_dir,
            chunk_size=chunk_size,
            approx_target_file_size_mb=target_file_size_mb,
            dataset_name=job.dataset_name,
            dtypes=job.dtypes,
        )
        finalize_dataset(temp_output_dir)
        move_completed_dataset(temp_output_dir, job.output_dir)
        return total_rows, part_count
    except Exception:
        shutil.rmtree(temp_output_dir, ignore_errors=True)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        jobs = resolve_jobs(args)
    except ValueError as exc:
        parser.error(str(exc))

    for job in jobs:
        total_rows, part_count = convert_job(job, args.chunk_size, args.target_file_size_mb)
        print(
            f"{job.dataset_name}: wrote {total_rows} rows from {job.input_path} "
            f"into {part_count} parquet part files under {job.output_dir}",
            flush=True,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
