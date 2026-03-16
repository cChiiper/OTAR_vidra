from pathlib import Path

import pytest

import tools.headtext as headtext
from tools import read_head


class FakeJsonRows:
    def __init__(self, rows: list[str]):
        self.rows = rows

    def collect(self) -> list[str]:
        return self.rows


class FakeDataFrame:
    def __init__(self, rows: list[str]):
        self.rows = rows
        self.limit_value = None

    def limit(self, n_lines: int):
        self.limit_value = n_lines
        return self

    def toJSON(self) -> FakeJsonRows:
        if self.limit_value is None:
            raise AssertionError("limit() must be called before toJSON()")
        return FakeJsonRows(self.rows[:self.limit_value])


class FakeReader:
    def __init__(self, rows: list[str]):
        self.rows = rows
        self.paths = []
        self.dataframes = []

    def parquet(self, path: str) -> FakeDataFrame:
        self.paths.append(path)
        dataframe = FakeDataFrame(self.rows)
        self.dataframes.append(dataframe)
        return dataframe


class FakeSparkSession:
    def __init__(self, rows: list[str]):
        self.read = FakeReader(rows)
        self.stop_called = False

    def stop(self):
        self.stop_called = True


def test_read_head_returns_first_rows_as_json_lines(tmp_path: Path):
    input_file = tmp_path / "example.snappy.parquet"
    fake_spark = FakeSparkSession([
        '{"fruit":"apple","count":1}',
        '{"fruit":"banana","count":2}',
        '{"fruit":"cherry","count":3}',
    ])

    result = read_head(input_file, 2, spark=fake_spark)

    assert result == '{"fruit":"apple","count":1}\n{"fruit":"banana","count":2}\n'
    assert fake_spark.read.paths == [str(input_file)]
    assert fake_spark.read.dataframes[0].limit_value == 2
    assert fake_spark.stop_called is False


def test_read_head_handles_short_dataset(tmp_path: Path):
    input_file = tmp_path / "short.snappy.parquet"
    fake_spark = FakeSparkSession([
        '{"fruit":"apple","count":1}',
        '{"fruit":"banana","count":2}',
    ])

    result = read_head(input_file, 5, spark=fake_spark)

    assert result == '{"fruit":"apple","count":1}\n{"fruit":"banana","count":2}\n'


def test_read_head_allows_zero_lines(tmp_path: Path):
    input_file = tmp_path / "empty.snappy.parquet"
    fake_spark = FakeSparkSession([
        '{"fruit":"apple","count":1}',
        '{"fruit":"banana","count":2}',
    ])

    result = read_head(input_file, 0, spark=fake_spark)

    assert result == ""


def test_read_head_rejects_negative_lines():
    input_file = Path("invalid.snappy.parquet")

    with pytest.raises(ValueError, match="n_lines must be >= 0"):
        read_head(input_file, -1)


def test_read_head_stops_created_spark_session(monkeypatch, tmp_path: Path):
    input_file = tmp_path / "managed.snappy.parquet"
    fake_spark = FakeSparkSession([
        '{"fruit":"apple","count":1}',
    ])

    monkeypatch.setattr(headtext, "create_spark_session", lambda: fake_spark)

    result = read_head(input_file, 1)

    assert result == '{"fruit":"apple","count":1}\n'
    assert fake_spark.stop_called is True
