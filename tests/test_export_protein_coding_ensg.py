import importlib.util
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_DIR / "tools" / "export_protein_coding_ensg.py"

spec = importlib.util.spec_from_file_location("export_protein_coding_ensg", SCRIPT_PATH)
exporter = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(exporter)


def test_exports_unique_protein_coding_ids_without_header(tmp_path):
    target_dir = tmp_path / "target"
    target_dir.mkdir()
    pq.write_table(
        pa.table(
            {
                "id": [
                    "ENSG00000000002",
                    "ENSG00000000001",
                    "ENSG00000000002",
                    "ENSG00000000003",
                    None,
                ],
                "biotype": [
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                    "lncRNA",
                    "protein_coding",
                ],
            }
        ),
        target_dir / "part-00000.parquet",
    )

    output = tmp_path / "all_protein_coding_ENSG.csv"

    assert exporter.main(
        [
            "--target-parquet",
            str(target_dir / "*.parquet"),
            "--output",
            str(output),
        ]
    ) == 0

    assert output.read_text(encoding="utf-8") == "ENSG00000000001\nENSG00000000002\n"
