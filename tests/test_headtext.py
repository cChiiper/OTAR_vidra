from pathlib import Path

from tools.headtext import read_head


def test_read_head_returns_first_lines(tmp_path: Path):

    input_file = tmp_path / "example.txt"

    input_file.write_text(
        "line1\nline2\nline3\nline4\n",
        encoding="utf-8",
    )

    result = read_head(input_file, 2)

    assert result == "line1\nline2\n"


def test_read_head_handles_short_file(tmp_path: Path):

    input_file = tmp_path / "short.txt"

    input_file.write_text(
        "only1\nonly2\n",
        encoding="utf-8",
    )

    result = read_head(input_file, 5)

    assert result == "only1\nonly2\n"
