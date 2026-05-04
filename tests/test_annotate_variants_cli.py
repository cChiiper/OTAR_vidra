import gzip
import json
import tarfile
from argparse import Namespace
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from tools import annotate_variants_cli as annotate


def test_plugin_name_and_amino_acid_parsing_helpers():
    assert annotate.plugin_name("AlphaMissense,file=/tmp/foo.tsv.gz") == "AlphaMissense"
    assert annotate.parse_amino_acids("A/V") == ("A", "V")
    assert annotate.parse_amino_acids("") == (None, None)


def test_build_plugin_args_only_enables_plugins_with_module_and_data(tmp_path):
    data_root = tmp_path / "data"
    modules_root = tmp_path / "plugins"
    (data_root / "alphamissense").mkdir(parents=True)
    (data_root / "CADD").mkdir(parents=True)
    modules_root.mkdir(parents=True)

    (data_root / "alphamissense" / "AlphaMissense_hg38.tsv.gz").write_text("x", encoding="utf-8")
    (data_root / "alphamissense" / "AlphaMissense_hg38.tsv.gz.tbi").write_text("idx", encoding="utf-8")
    (modules_root / "AlphaMissense.pm").write_text("1;", encoding="utf-8")
    (data_root / "CADD" / "whole_genome_SNVs.tsv.gz").write_text("x", encoding="utf-8")
    (data_root / "CADD" / "whole_genome_SNVs.tsv.gz.tbi").write_text("idx", encoding="utf-8")
    (modules_root / "CADD.pm").write_text("1;", encoding="utf-8")

    plugin_args, statuses = annotate.build_plugin_args(
        str(data_root),
        str(modules_root),
        {"AlphaMissense", "CADD"},
    )

    assert "--plugin" in plugin_args
    assert any("AlphaMissense" in value for value in plugin_args)
    assert any("CADD,snv=" in value for value in plugin_args)
    assert not any("REVEL" in value for value in plugin_args)
    assert any(value == "Blosum62" for value in plugin_args)
    assert any(status["name"] == "AlphaMissense" and status["required"] for status in statuses)


def test_build_plugin_args_raises_for_missing_required_plugin(tmp_path):
    data_root = tmp_path / "data"
    modules_root = tmp_path / "plugins"
    data_root.mkdir()
    modules_root.mkdir()

    try:
        annotate.build_plugin_args(str(data_root), str(modules_root), {"AlphaMissense"})
    except RuntimeError as exc:
        assert "AlphaMissense" in str(exc)
    else:
        raise AssertionError("Expected required plugin check to fail")


def test_resolve_vep_parallelism_auto_and_limits():
    assert annotate.resolve_vep_parallelism(0, 0, 0, cpu_count=8) == (1, 8, 8)
    assert annotate.resolve_vep_parallelism(10, 0, 0, cpu_count=8) == (4, 2, 8)
    assert annotate.resolve_vep_parallelism(2, 4, 0, cpu_count=8) == (2, 4, 8)
    assert annotate.resolve_vep_parallelism(10, 3, 5, cpu_count=8) == (3, 5, 8)


def test_plan_vep_chunks_handles_zero_single_and_large_chunk_requests():
    single = "1_100_A_G"

    assert annotate.plan_vep_chunks([], 4) == []
    assert annotate.plan_vep_chunks([single], 4) == [[single]]

    variant_ids = [
        "2_300_G_A",
        "1_200_C_T",
        "1_100_A_G",
        "1_100_A_G",
    ]
    chunks = annotate.plan_vep_chunks(variant_ids, 10)

    assert len(chunks) == 3
    assert [variant for chunk in chunks for variant in chunk] == annotate.dedupe_and_sort_variant_ids(variant_ids)


def test_load_manifest_variant_ids_dedupes_and_sorts(tmp_path):
    manifest_dir = tmp_path / "vidra_analysis_ready_manifest"
    manifest_dir.mkdir()

    pq.write_table(
        pa.Table.from_pylist(
            [
                {"variant": "2_300_G_A"},
                {"variant": "1_200_C_T"},
                {"variant": "1_100_A_G"},
                {"variant": "1_100_A_G"},
            ]
        ),
        manifest_dir / "part-00000.parquet",
    )

    assert annotate.load_manifest_variant_ids(str(manifest_dir)) == [
        "1_100_A_G",
        "1_200_C_T",
        "2_300_G_A",
    ]


def test_merge_vep_chunk_outputs_preserves_chunk_order(tmp_path):
    chunk_a = tmp_path / "chunk_a.json"
    chunk_b = tmp_path / "chunk_b.json"
    merged = tmp_path / "merged.json"

    chunk_a.write_text('{"id":"1_100_A_G"}', encoding="utf-8")
    chunk_b.write_text('{"id":"2_200_C_T"}\n', encoding="utf-8")

    annotate.merge_vep_chunk_outputs([chunk_a, chunk_b], merged)

    assert merged.read_text(encoding="utf-8") == '{"id":"1_100_A_G"}\n{"id":"2_200_C_T"}\n'


def test_run_vep_executes_directly_without_nested_docker(tmp_path, monkeypatch):
    cache_root = tmp_path / "cache"
    (cache_root / "homo_sapiens").mkdir(parents=True)
    plugins_root = tmp_path / "plugins"
    plugins_root.mkdir()
    plugin_data_root = tmp_path / "plugin_data"
    plugin_data_root.mkdir()

    seen = {}

    def fake_run(cmd, check):
        seen["cmd"] = cmd
        assert check is True

    monkeypatch.setattr(annotate.subprocess, "run", fake_run)

    args = Namespace(
        vep_dir_cache=str(cache_root),
        vep_plugins_dir=str(plugins_root),
        vep_plugin_data_dir=str(plugin_data_root),
        required_vep_plugins="",
        vep_plugin_args=["--plugin", "Blosum62"],
    )

    annotate.run_vep(
        args,
        tmp_path / "input.tmp",
        tmp_path / "output.json",
        vep_forks=2,
        vep_buffer_size=10000,
    )

    assert seen["cmd"][0] == "vep"
    assert "docker" not in seen["cmd"]
    assert "--everything" not in seen["cmd"]
    assert seen["cmd"][seen["cmd"].index("--sift") + 1] == "b"
    assert seen["cmd"][seen["cmd"].index("--polyphen") + 1] == "b"
    assert "--no_stats" in seen["cmd"]
    assert seen["cmd"][seen["cmd"].index("--fork") + 1] == "2"
    assert seen["cmd"][seen["cmd"].index("--buffer_size") + 1] == "10000"


def test_prepare_runtime_resources_extracts_local_vep_cache_archive(tmp_path):
    archive_root = tmp_path / "archive_root"
    cache_source = archive_root / "homo_sapiens" / "115_GRCh38"
    cache_source.mkdir(parents=True)
    (cache_source / "info.txt").write_text("ok", encoding="utf-8")

    cache_archive = tmp_path / "vep_cache.tar.gz"
    with tarfile.open(cache_archive, "w:gz") as handle:
        handle.add(archive_root / "homo_sapiens", arcname="homo_sapiens")

    args = Namespace(
        vep_dir_cache="",
        vep_cache_archive=str(cache_archive),
        vep_plugins_dir=str(tmp_path / "plugins"),
        vep_plugin_data_dir=str(tmp_path / "plugin_data"),
        foldx_file="",
    )

    runtime_args = annotate.prepare_runtime_resources(args, tmp_path / "runtime")

    extracted = annotate.Path(runtime_args.vep_dir_cache) / "homo_sapiens" / "115_GRCh38" / "info.txt"
    assert extracted.read_text(encoding="utf-8") == "ok"


def test_prepare_runtime_resources_localizes_plugin_sidecars_from_gs(monkeypatch, tmp_path):
    plugin_source = "gs://bucket/alphamissense/AlphaMissense_hg38.tsv.gz"
    copied = []

    def fake_path_exists(path):
        return path in {
            plugin_source,
            f"{plugin_source}.tbi",
        }

    def fake_copy(source, destination):
        copied.append((source, str(destination)))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("x", encoding="utf-8")

    monkeypatch.setattr(annotate, "path_exists", fake_path_exists)
    monkeypatch.setattr(annotate, "copy_file_to_local", fake_copy)

    args = Namespace(
        vep_dir_cache=str(tmp_path / "cache"),
        vep_cache_archive="",
        vep_plugins_dir=str(tmp_path / "plugins"),
        vep_plugin_data_dir="gs://bucket",
        foldx_file="",
    )

    runtime_args = annotate.prepare_runtime_resources(args, tmp_path / "runtime")

    assert runtime_args.vep_plugin_data_dir.endswith("plugin_resources")
    assert (
        plugin_source,
        str(tmp_path / "runtime" / "plugin_resources" / "alphamissense" / "AlphaMissense_hg38.tsv.gz"),
    ) in copied
    assert (
        f"{plugin_source}.tbi",
        str(tmp_path / "runtime" / "plugin_resources" / "alphamissense" / "AlphaMissense_hg38.tsv.gz.tbi"),
    ) in copied


def test_prepare_runtime_resources_generates_missing_alphamissense_index(monkeypatch, tmp_path):
    plugin_source = "gs://bucket/AlphaMissense_hg38.tsv.gz"
    copied = []
    seen = {}

    def fake_path_exists(path):
        return path == plugin_source

    def fake_copy(source, destination):
        copied.append((source, str(destination)))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("x", encoding="utf-8")

    def fake_run(cmd, check):
        seen["cmd"] = cmd
        assert check is True
        index_path = tmp_path / "runtime" / "plugin_resources" / "alphamissense" / "AlphaMissense_hg38.tsv.gz.tbi"
        index_path.write_text("idx", encoding="utf-8")

    monkeypatch.setattr(annotate, "path_exists", fake_path_exists)
    monkeypatch.setattr(annotate, "copy_file_to_local", fake_copy)
    monkeypatch.setattr(annotate.subprocess, "run", fake_run)

    args = Namespace(
        vep_dir_cache=str(tmp_path / "cache"),
        vep_cache_archive="",
        vep_plugins_dir=str(tmp_path / "plugins"),
        vep_plugin_data_dir="gs://bucket",
        foldx_file="",
    )

    annotate.prepare_runtime_resources(args, tmp_path / "runtime")

    assert copied == [
        (
            plugin_source,
            str(tmp_path / "runtime" / "plugin_resources" / "alphamissense" / "AlphaMissense_hg38.tsv.gz"),
        )
    ]
    assert seen["cmd"] == [
        "tabix",
        "-s",
        "1",
        "-b",
        "2",
        "-e",
        "2",
        "-f",
        "-S",
        "1",
        str(tmp_path / "runtime" / "plugin_resources" / "alphamissense" / "AlphaMissense_hg38.tsv.gz"),
    ]


def test_merge_annotation_tables():
    existing = pa.Table.from_pylist(
        [{"variant": "1_100_A_G", "as_blosum62": 0.1}],
        schema=pa.schema([("variant", pa.string()), ("as_blosum62", pa.float64())]),
    )
    new = pa.Table.from_pylist(
        [{"variant": "1_200_C_T", "as_blosum62": 0.2}],
        schema=pa.schema([("variant", pa.string()), ("as_blosum62", pa.float64())]),
    )

    merged = annotate.merge_annotation_tables(existing, new)
    merged_rows = {row["variant"]: row["as_blosum62"] for row in merged.to_pylist()}
    assert merged_rows == {"1_100_A_G": 0.1, "1_200_C_T": 0.2}


def test_foldx_lookup_and_annotation_table_use_matching_keys(tmp_path):
    raw_records = {
        "10_100042447_G_A": {
            "variant": "10_100042447_G_A",
            "most_severe_consequence": "missense_variant",
            "clinicalSignificance": "pathogenic",
            "swissprot": "P12345.2",
            "protein_start": "42",
            "ref_aa": "A",
            "alt_aa": "V",
        }
    }

    keys = annotate.collect_needed_foldx_keys(raw_records)
    assert keys == {("P12345", 42, "A", "V")}

    foldx_file = tmp_path / "foldx.csv.gz"
    with gzip.open(foldx_file, "wt", encoding="utf-8") as handle:
        handle.write("uniprot_accession,uniprot_position,alphafold_fragment_id,alphafold_fragment_position,wild_type,mutated_type,foldx_ddg,plddt\n")
        handle.write("P12345,42,frag,42,A,V,1.75,82.0\n")

    foldx_lookup = annotate.load_foldx_lookup(str(foldx_file), keys)
    table = annotate.build_annotation_table(["10_100042447_G_A"], raw_records, foldx_lookup=foldx_lookup)
    row = table.to_pylist()[0]

    assert row["foldxDdq_raw"] == 1.75
    assert row["as_plddt"] == 0.82
    assert row["as_clinicalSignificance"] > 0.0


def test_validate_runtime_args_rejects_invalid_parallel_settings():
    with pytest.raises(ValueError):
        annotate.validate_runtime_args(Namespace(vep_parallel=-1, vep_forks=0, vep_buffer_size=10000))
    with pytest.raises(ValueError):
        annotate.validate_runtime_args(Namespace(vep_parallel=0, vep_forks=-1, vep_buffer_size=10000))
    with pytest.raises(ValueError):
        annotate.validate_runtime_args(Namespace(vep_parallel=0, vep_forks=0, vep_buffer_size=0))


def test_parallel_and_serial_bulk_annotation_produce_identical_final_table(tmp_path, monkeypatch):
    cache_root = tmp_path / "cache"
    (cache_root / "homo_sapiens").mkdir(parents=True)

    def fake_run(cmd, check):
        assert check is True
        input_path = Path(cmd[cmd.index("-i") + 1])
        output_path = Path(cmd[cmd.index("-o") + 1])

        with input_path.open(encoding="utf-8") as handle, output_path.open("w", encoding="utf-8") as out:
            for line in handle:
                chrom, pos, _, ref, alt = line.strip().split()
                variant_id = f"{chrom}_{pos}_{ref}_{alt}"
                record = {
                    "id": variant_id,
                    "input": f"{chrom} {pos} . {ref} {alt}",
                    "most_severe_consequence": "missense_variant" if alt in {"G", "T"} else "synonymous_variant",
                    "transcript_consequences": [
                        {
                            "canonical": 1,
                            "blosum62": "1.5",
                            "conservation": "2.0",
                            "sift_score": "0.2",
                            "polyphen_score": "0.8",
                            "cadd_phred": "25.0",
                            "am_pathogenicity": "0.4",
                            "revel": "0.7",
                            "primateai_score": "0.1",
                            "swissprot": "P12345.1",
                            "protein_start": str((int(pos) % 50) + 1),
                            "amino_acids": "A/V",
                        }
                    ],
                    "colocated_variants": [
                        {"clin_sig": ["pathogenic" if alt == "T" else "benign"]}
                    ],
                }
                out.write(json.dumps(record) + "\n")

    monkeypatch.setattr(annotate.subprocess, "run", fake_run)

    base_args = {
        "vep_dir_cache": str(cache_root),
        "vep_plugins_dir": str(tmp_path / "plugins"),
        "vep_plugin_data_dir": str(tmp_path / "plugin_data"),
        "vep_plugin_args": ["--plugin", "Blosum62"],
        "required_vep_plugins": "",
        "vep_forks": 1,
        "vep_buffer_size": 10000,
    }
    variant_ids = [
        "1_100_A_G",
        "1_200_C_T",
        "1_100_A_G",
        "2_300_G_A",
        "3_400_T_C",
    ]

    serial_args = Namespace(**base_args, vep_parallel=1)
    serial_json = tmp_path / "serial.json"
    annotate.annotate_variant_ids(serial_args, "bulk_variant_annotation", variant_ids, str(serial_json))
    serial_records = annotate.parse_vep_json_outputs([str(serial_json)])
    serial_table = annotate.build_annotation_table(annotate.dedupe_and_sort_variant_ids(variant_ids), serial_records)

    parallel_args = Namespace(**base_args, vep_parallel=2)
    parallel_json = tmp_path / "parallel.json"
    annotate.annotate_variant_ids(parallel_args, "bulk_variant_annotation", variant_ids, str(parallel_json))
    parallel_records = annotate.parse_vep_json_outputs([str(parallel_json)])
    parallel_table = annotate.build_annotation_table(annotate.dedupe_and_sort_variant_ids(variant_ids), parallel_records)

    assert serial_table.to_pylist() == parallel_table.to_pylist()
