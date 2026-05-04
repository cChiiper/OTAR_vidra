import csv
import lzma
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "convert_az_phewas_parquet.py"
)


def write_fixture_csv_xz(path: Path, rows: list[dict[str, object]]) -> pd.DataFrame:
    with lzma.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return pd.DataFrame(rows)


def assert_parquet_dataset_matches(output_dir: Path, expected_df: pd.DataFrame) -> None:
    assert output_dir.is_dir()
    assert (output_dir / "_SUCCESS").is_file()

    part_paths = sorted(output_dir.glob("part-*.snappy.parquet"))
    assert part_paths
    assert all(re.fullmatch(r"part-\d{5}\.snappy\.parquet", path.name) for path in part_paths)

    actual_df = pd.concat(
        [pq.read_table(path).to_pandas() for path in part_paths],
        ignore_index=True,
    ).sort_values(["Gene", "Phenotype", "CollapsingModel"]).reset_index(drop=True)

    expected_df = expected_df.sort_values(["Gene", "Phenotype", "CollapsingModel"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(actual_df, expected_df, check_dtype=False)


def test_cli_converts_binary_dataset_in_single_dataset_mode(tmp_path):
    rows = [
        {
            "Gene": "A1BG",
            "Phenotype": "120000#Example phenotype one",
            "CollapsingModel": "flexnonsynmtr",
            "Type": "Binary",
            "pValue": 0.0667,
            "Category": "Example category",
            "nSamples": 120745,
            "BinNcases": 43153,
            "BinQVcases": 117,
            "BinNcontrols": 77592,
            "BinQVcontrols": 258,
            "BinCaseFreq": 0.00271128310893797,
            "BinCtrlFreq": 0.0033250850603155,
            "BinOddsRatio": 0.8149,
            "BinOddsRatioLCI": 0.6548,
            "BinOddsRatioUCI": 1.0142,
        },
        {
            "Gene": "A1BG",
            "Phenotype": "120001#Example phenotype two",
            "CollapsingModel": "raredmg",
            "Type": "Binary",
            "pValue": 0.0479,
            "Category": "Example category",
            "nSamples": 122543,
            "BinNcases": 9148,
            "BinQVcases": 2,
            "BinNcontrols": 113395,
            "BinQVcontrols": 3,
            "BinCaseFreq": 0.000218627022299956,
            "BinCtrlFreq": 2.64561929538339e-05,
            "BinOddsRatio": 8.2653,
            "BinOddsRatioLCI": 1.3809,
            "BinOddsRatioUCI": 49.4728,
        },
        {
            "Gene": "A1CF",
            "Phenotype": "120005#Example phenotype three",
            "CollapsingModel": "rec",
            "Type": "Binary",
            "pValue": 0.0146,
            "Category": "Example category",
            "nSamples": 143003,
            "BinNcases": 5695,
            "BinQVcases": 2,
            "BinNcontrols": 137308,
            "BinQVcontrols": 3,
            "BinCaseFreq": 0.000351185250219491,
            "BinCtrlFreq": 2.18486905351473e-05,
            "BinOddsRatio": 16.0788,
            "BinOddsRatioLCI": 2.6861,
            "BinOddsRatioUCI": 96.2478,
        },
    ]

    input_path = tmp_path / "azphewas-com-470k-phewas-binary.csv.xz"
    output_dir = tmp_path / "azphewas-com-470k-phewas-binary"
    expected_df = write_fixture_csv_xz(input_path, rows)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--dataset",
            "binary",
            "--input",
            str(input_path),
            "--output-dir",
            str(output_dir),
            "--chunk-size",
            "2",
            "--target-file-size-mb",
            "0.000001",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "binary: wrote 3 rows" in result.stdout
    assert_parquet_dataset_matches(output_dir, expected_df)


def test_cli_converts_binary_and_quantitative_datasets_in_all_mode(tmp_path):
    binary_rows = [
        {
            "Gene": "A1BG",
            "Phenotype": "120000#Example phenotype one",
            "CollapsingModel": "flexnonsynmtr",
            "Type": "Binary",
            "pValue": 0.0667,
            "Category": "Example category",
            "nSamples": 120745,
            "BinNcases": 43153,
            "BinQVcases": 117,
            "BinNcontrols": 77592,
            "BinQVcontrols": 258,
            "BinCaseFreq": 0.00271128310893797,
            "BinCtrlFreq": 0.0033250850603155,
            "BinOddsRatio": 0.8149,
            "BinOddsRatioLCI": 0.6548,
            "BinOddsRatioUCI": 1.0142,
        },
        {
            "Gene": "A1CF",
            "Phenotype": "120005#Example phenotype three",
            "CollapsingModel": "rec",
            "Type": "Binary",
            "pValue": 0.0146,
            "Category": "Example category",
            "nSamples": 143003,
            "BinNcases": 5695,
            "BinQVcases": 2,
            "BinNcontrols": 137308,
            "BinQVcontrols": 3,
            "BinCaseFreq": 0.000351185250219491,
            "BinCtrlFreq": 2.18486905351473e-05,
            "BinOddsRatio": 16.0788,
            "BinOddsRatioLCI": 2.6861,
            "BinOddsRatioUCI": 96.2478,
        },
    ]
    quantitative_rows = [
        {
            "Gene": "A1BG",
            "Phenotype": "102#Pulse rate automated reading",
            "CollapsingModel": "UR",
            "Type": "Quantitative",
            "pValue": 0.00729629075022466,
            "Category": "Example category",
            "nSamples": 395230,
            "YesQV": 63,
            "NoQV": 395167,
            "ConCaseMedian": 0.34377956603166,
            "ConCtrMedian": -0.00777789704557489,
            "beta": 0.336764612364056,
            "ConBetaSe": 0.125516834752413,
            "LCI": 0.0907553834019236,
            "UCI": 0.582773841326189,
        },
        {
            "Gene": "A1CF",
            "Phenotype": "200#Standing height",
            "CollapsingModel": "ptv",
            "Type": "Quantitative",
            "pValue": 0.0001,
            "Category": "Example category",
            "nSamples": 395230,
            "YesQV": 20,
            "NoQV": 395210,
            "ConCaseMedian": 1.1,
            "ConCtrMedian": 0.1,
            "beta": 1.0,
            "ConBetaSe": 0.2,
            "LCI": 0.6,
            "UCI": 1.4,
        },
    ]

    binary_input = tmp_path / "azphewas-com-470k-phewas-binary.csv.xz"
    quantitative_input = tmp_path / "azphewas-com-470k-phewas-quantitative.csv.xz"
    binary_output = tmp_path / "azphewas-com-470k-phewas-binary"
    quantitative_output = tmp_path / "azphewas-com-470k-phewas-quantitative"

    expected_binary_df = write_fixture_csv_xz(binary_input, binary_rows)
    expected_quantitative_df = write_fixture_csv_xz(quantitative_input, quantitative_rows)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--dataset",
            "all",
            "--binary-input",
            str(binary_input),
            "--binary-output-dir",
            str(binary_output),
            "--quantitative-input",
            str(quantitative_input),
            "--quantitative-output-dir",
            str(quantitative_output),
            "--chunk-size",
            "2",
            "--target-file-size-mb",
            "0.000001",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "binary: wrote 2 rows" in result.stdout
    assert "quantitative: wrote 2 rows" in result.stdout
    assert_parquet_dataset_matches(binary_output, expected_binary_df)
    assert_parquet_dataset_matches(quantitative_output, expected_quantitative_df)


def test_cli_requires_explicit_paths_when_not_supplied():
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--dataset",
            "binary",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 2
    assert "binary dataset requires an input path via --input/--binary-input" in result.stderr
