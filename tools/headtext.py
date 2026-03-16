#!/usr/bin/env python3

from pathlib import Path
import argparse


def read_head(input_path: Path, n_lines: int) -> str:
    lines = []
    with input_path.open("r", encoding="utf-8") as f:
        for _ in range(n_lines):
            line = f.readline()
            if not line:
                break
            lines.append(line)
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--n_lines", type=int, default=5)
    args = parser.parse_args()

    text = read_head(args.input, args.n_lines)
    args.output.write_text(text, encoding="utf-8")

    print(text, end="")


if __name__ == "__main__":
    main()
