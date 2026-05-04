import csv
import inspect
import io
import sys
import types
from pathlib import Path

import pytest

from tools import prepare_analysis_input as prep

TESTDATA_DIR = Path(__file__).resolve().parents[1] / "testdata"


def _read_dict_rows(path: Path, delimiter: str = ","):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def test_normalize_az_variant_id_converts_dashes_to_underscores():
    assert prep.normalize_az_variant_id("10-100042447-G-A") == "10_100042447_G_A"


def test_normalize_az_variant_id_rejects_malformed_values():
    assert prep.normalize_az_variant_id("bad-value") is None


def test_clinvar_confidence_excludes_current_and_legacy_no_assertion_labels():
    assert prep.is_excluded_clinvar_confidence("no classification provided") is True
    assert prep.is_excluded_clinvar_confidence("no assertion provided") is True
    assert prep.is_excluded_clinvar_confidence("reviewed by expert panel") is False


def test_coloc_filters_use_vidra4_right_study_type_and_blood_biosamples():
    assert prep.is_allowed_right_study_type("eqtl") is True
    assert prep.is_allowed_right_study_type("pqtl") is True
    assert prep.is_allowed_right_study_type("sqtl") is False
    assert prep.is_allowed_blood_biosample_id("UBERON_0000178") is True
    assert prep.is_allowed_blood_biosample_id("EFO_0005292") is False


def test_coding_consequence_filter_uses_so_notation():
    assert prep.is_allowed_coding_consequence("SO:0001583") is True
    assert prep.is_allowed_coding_consequence("SO_0001589") is True
    assert prep.is_allowed_coding_consequence("SO_0001627") is False


def test_burden_policy_allows_450k_and_470k_with_vidra2_method_subset():
    assert prep.burden_row_is_allowed("AstraZeneca PheWAS Portal", "UK Biobank 470k", "ptv") is True
    assert prep.burden_row_is_allowed("AstraZeneca PheWAS Portal", "UK Biobank 450k", "ptv") is True
    assert prep.burden_row_is_allowed("AstraZeneca PheWAS Portal", "UK Biobank 470k", "flexdmg") is False
    assert prep.burden_row_is_allowed("AstraZeneca PheWAS Portal", "UK Biobank 500k", "ptv") is False


def test_analysis_source_cols_exclude_post_burden_fields():
    assert "bO" not in prep.ANALYSIS_SOURCE_COLS
    assert "bOse" not in prep.ANALYSIS_SOURCE_COLS
    assert "bO" in prep.ANALYSIS_OUTPUT_COLS
    assert "bOse" in prep.ANALYSIS_OUTPUT_COLS


def test_resolve_parquet_input_accepts_gs_paths_without_local_checks():
    assert prep.resolve_parquet_input("gs://vidra-2-0/example.parquet") == "gs://vidra-2-0/example.parquet"


def test_path_exists_uses_requester_pays_gcsfs_when_project_is_provided(monkeypatch):
    captured = {}

    class FakeFS:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def exists(self, path):
            captured["exists_path"] = path
            return True

    monkeypatch.setitem(sys.modules, "gcsfs", types.SimpleNamespace(GCSFileSystem=FakeFS))

    assert prep.path_exists("gs://bucket/path/file", gcp_project="open-targets-eu-dev") is True
    assert captured["project"] == "open-targets-eu-dev"
    assert captured["requester_pays"] is True
    assert captured["exists_path"] == "bucket/path/file"


def test_load_gene_ids_from_csv_uses_requester_pays_gcsfs_when_project_is_provided(monkeypatch):
    captured = {}

    class FakeFS:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def open(self, path, mode):
            captured["open_path"] = path
            captured["open_mode"] = mode
            return io.StringIO("gene\nENSG00000141510\n")

    monkeypatch.setitem(sys.modules, "gcsfs", types.SimpleNamespace(GCSFileSystem=FakeFS))

    genes = prep.load_gene_ids_from_csv("gs://bucket/genes.csv", gcp_project="open-targets-eu-dev")
    assert genes == ["ENSG00000141510"]
    assert captured["project"] == "open-targets-eu-dev"
    assert captured["requester_pays"] is True
    assert captured["open_path"] == "bucket/genes.csv"
    assert captured["open_mode"] == "r"


def test_configure_spark_defaults_updates_shared_preparation_defaults():
    original = dict(prep.SPARK_SESSION_DEFAULTS)
    try:
        prep.configure_spark_defaults(
            master="local[4]",
            driver_memory="8g",
            shuffle_partitions=64,
            default_parallelism=16,
        )
        assert prep.SPARK_SESSION_DEFAULTS["master"] == "local[4]"
        assert prep.SPARK_SESSION_DEFAULTS["driver_memory"] == "8g"
        assert prep.SPARK_SESSION_DEFAULTS["shuffle_partitions"] == 64
        assert prep.SPARK_SESSION_DEFAULTS["default_parallelism"] == 16
    finally:
        prep.SPARK_SESSION_DEFAULTS.update(original)


def test_configure_local_spark_runtime_uses_explicit_writable_dirs(tmp_path, monkeypatch):
    home_dir = tmp_path / "home"
    spark_local_dir = tmp_path / "spark_local"

    monkeypatch.setenv("VIDRA_SPARK_USER_HOME", str(home_dir))
    monkeypatch.setenv("VIDRA_SPARK_LOCAL_DIRS", str(spark_local_dir))
    monkeypatch.setenv("HOME", "?")

    runtime = prep.configure_local_spark_runtime()

    assert runtime["home"] == str(home_dir)
    assert runtime["ivy"] == str(home_dir / ".ivy2")
    assert runtime["ivy_cache"] == str(home_dir / ".ivy2" / "cache")
    assert runtime["spark_local"] == str(spark_local_dir)
    assert home_dir.is_dir()
    assert (home_dir / ".ivy2").is_dir()
    assert (home_dir / ".ivy2" / "cache").is_dir()
    assert spark_local_dir.is_dir()


def test_mock_az_fixture_rows_are_valid_for_normalization_and_or_ci_handling():
    rows = _read_dict_rows(TESTDATA_DIR / "mock_az_variants.csv")

    assert rows

    for row in rows:
        normalized = prep.normalize_az_variant_id(row["Variant"])
        assert normalized is not None
        assert normalized.count("_") == 3
        assert row["Model"] == "allelic"
        assert float(row["Odds ratio"]) > 0.0
        assert float(row["Odds ratio LCI"]) > 0.0
        assert float(row["Odds ratio UCI"]) > 0.0


def test_mock_az_supporting_fixtures_cover_fixture_genes_and_phenotypes():
    az_rows = _read_dict_rows(TESTDATA_DIR / "mock_az_variants.csv")
    efo_rows = _read_dict_rows(TESTDATA_DIR / "mock_az_phewas_to_efo.tsv", delimiter="\t")

    gene_rows = []
    with (TESTDATA_DIR / "mock_all_human_protein_coding_genes.csv").open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                gene_rows.append(line.split(","))

    az_genes = {row["Gene"].replace("'", "").strip() for row in az_rows}
    az_phenotypes = {row["Phenotype"] for row in az_rows}
    mapped_genes = {parts[1].strip() for parts in gene_rows if len(parts) >= 2}
    mapped_phenotypes = {
        row["PROPERTY_VALUE"]
        for row in efo_rows
        if "AstraZeneca" in row["STUDY"]
    }

    assert az_genes <= mapped_genes
    assert az_phenotypes <= mapped_phenotypes



def test_parse_args_defaults_step1_debug_outputs_to_empty(monkeypatch):
    monkeypatch.setattr(
        prep.sys,
        "argv",
        [
            "prepare_analysis_input.py",
            "--genes",
            "genes.csv",
            "--colocalisation_threshold",
            "0.7",
            "--coloc_data_dir",
            "coloc",
            "--credible_set_dir",
            "credible",
            "--study_data_dir",
            "study",
            "--variant_data_dir",
            "variant",
            "--burden_evidence_dir",
            "burden",
            "--clinvar_evidence_dir",
            "clinvar",
            "--az_variants_file",
            "az.csv",
            "--az_mapping_file",
            "az_map.tsv",
            "--az_gene_map_file",
            "genes_map.csv",
        ],
    )

    args = prep.parse_args()

    assert args.coloc_output == ""
    assert args.gwas_output == ""
    assert args.qtl_output == ""
    assert args.burden_output == ""
    assert args.az_output == ""
    assert args.az_burden_binary_dir == ""
    assert args.az_burden_quantitative_dir == ""
    assert args.target_data_dir == ""
    assert args.lead_variant_effect_dir == ""


def test_validate_raw_az_burden_args_warns_and_keeps_release_only(monkeypatch):
    messages = []
    args = type(
        "Args",
        (),
        {
            "az_burden_binary_dir": "",
            "az_burden_quantitative_dir": "",
            "target_data_dir": "",
        },
    )()

    monkeypatch.setattr(prep, "log", messages.append)

    prep.validate_raw_az_burden_args(args)

    assert prep.raw_az_burden_inputs_enabled(args) is False
    assert messages == []


def test_validate_step1_input_paths_accepts_present_local_inputs(tmp_path):
    genes_file = tmp_path / "genes.csv"
    genes_file.write_text("gene\nENSG00000141510\n", encoding="utf-8")

    for dirname in [
        "coloc",
        "credible",
        "study",
        "variant",
        "burden",
        "clinvar",
    ]:
        (tmp_path / dirname).mkdir()

    az_file = tmp_path / "az.csv.bz2"
    az_file.write_text("mock", encoding="utf-8")
    az_map = tmp_path / "az_map.tsv"
    az_map.write_text("mock", encoding="utf-8")
    az_gene_map = tmp_path / "az_gene_map.csv"
    az_gene_map.write_text("mock", encoding="utf-8")

    args = type(
        "Args",
        (),
        {
            "gcp_project": "",
            "genes": str(genes_file),
            "coloc_data_dir": str(tmp_path / "coloc"),
            "credible_set_dir": str(tmp_path / "credible"),
            "study_data_dir": str(tmp_path / "study"),
            "variant_data_dir": str(tmp_path / "variant"),
            "burden_evidence_dir": str(tmp_path / "burden"),
            "clinvar_evidence_dir": str(tmp_path / "clinvar"),
            "az_variants_file": str(az_file),
            "az_mapping_file": str(az_map),
            "az_gene_map_file": str(az_gene_map),
            "lead_variant_effect_dir": "",
            "az_burden_binary_dir": "",
            "az_burden_quantitative_dir": "",
            "target_data_dir": "",
        },
    )()

    prep.validate_step1_input_paths(args)


def test_validate_step1_input_paths_fails_before_spark_for_missing_optional_and_required_inputs(tmp_path):
    genes_file = tmp_path / "genes.csv"
    genes_file.write_text("gene\nENSG00000141510\n", encoding="utf-8")
    (tmp_path / "coloc").mkdir()
    (tmp_path / "credible").mkdir()
    (tmp_path / "study").mkdir()
    (tmp_path / "variant").mkdir()
    (tmp_path / "burden").mkdir()
    (tmp_path / "clinvar").mkdir()
    (tmp_path / "target").mkdir()
    az_file = tmp_path / "az.csv.bz2"
    az_file.write_text("mock", encoding="utf-8")
    az_map = tmp_path / "az_map.tsv"
    az_map.write_text("mock", encoding="utf-8")
    az_gene_map = tmp_path / "az_gene_map.csv"
    az_gene_map.write_text("mock", encoding="utf-8")

    args = type(
        "Args",
        (),
        {
            "gcp_project": "",
            "genes": str(genes_file),
            "coloc_data_dir": str(tmp_path / "coloc"),
            "credible_set_dir": str(tmp_path / "credible"),
            "study_data_dir": str(tmp_path / "study"),
            "variant_data_dir": str(tmp_path / "variant"),
            "burden_evidence_dir": str(tmp_path / "burden"),
            "clinvar_evidence_dir": str(tmp_path / "clinvar"),
            "az_variants_file": str(az_file),
            "az_mapping_file": str(az_map),
            "az_gene_map_file": str(az_gene_map),
            "lead_variant_effect_dir": str(tmp_path / "missing_lead"),
            "az_burden_binary_dir": str(tmp_path / "missing_raw_binary"),
            "az_burden_quantitative_dir": str(tmp_path / "missing_raw_quant"),
            "target_data_dir": str(tmp_path / "target"),
        },
    )()

    with pytest.raises(FileNotFoundError, match="Lead variant effect input not found or unreadable"):
        prep.validate_step1_input_paths(args)


def test_validate_step1_input_paths_fails_before_spark_for_missing_raw_az_burden_inputs(tmp_path):
    genes_file = tmp_path / "genes.csv"
    genes_file.write_text("gene\nENSG00000141510\n", encoding="utf-8")

    for dirname in [
        "coloc",
        "credible",
        "study",
        "variant",
        "burden",
        "clinvar",
        "target",
    ]:
        (tmp_path / dirname).mkdir()

    az_file = tmp_path / "az.csv.bz2"
    az_file.write_text("mock", encoding="utf-8")
    az_map = tmp_path / "az_map.tsv"
    az_map.write_text("mock", encoding="utf-8")
    az_gene_map = tmp_path / "az_gene_map.csv"
    az_gene_map.write_text("mock", encoding="utf-8")

    args = type(
        "Args",
        (),
        {
            "gcp_project": "",
            "genes": str(genes_file),
            "coloc_data_dir": str(tmp_path / "coloc"),
            "credible_set_dir": str(tmp_path / "credible"),
            "study_data_dir": str(tmp_path / "study"),
            "variant_data_dir": str(tmp_path / "variant"),
            "burden_evidence_dir": str(tmp_path / "burden"),
            "clinvar_evidence_dir": str(tmp_path / "clinvar"),
            "az_variants_file": str(az_file),
            "az_mapping_file": str(az_map),
            "az_gene_map_file": str(az_gene_map),
            "lead_variant_effect_dir": "",
            "az_burden_binary_dir": str(tmp_path / "missing_raw_binary"),
            "az_burden_quantitative_dir": str(tmp_path / "missing_raw_quant"),
            "target_data_dir": str(tmp_path / "target"),
        },
    )()

    with pytest.raises(FileNotFoundError, match="Raw AZ burden binary input not found or unreadable"):
        prep.validate_step1_input_paths(args)


def test_log_step1_optional_input_modes_reports_disabled_modes(monkeypatch):
    messages = []
    args = type(
        "Args",
        (),
        {
            "lead_variant_effect_dir": "",
            "az_burden_binary_dir": "",
            "az_burden_quantitative_dir": "",
            "target_data_dir": "",
        },
    )()

    monkeypatch.setattr(prep, "log", messages.append)

    prep.log_step1_optional_input_modes(args)

    assert messages == [
        "SOURCE lead_variant_effect: disabled -> using release credible_set + study + variant inputs",
        "SOURCE raw AZ burden: disabled -> using only release gene_burden evidence",
    ]


def test_log_step1_optional_input_modes_reports_enabled_modes(monkeypatch):
    messages = []
    args = type(
        "Args",
        (),
        {
            "lead_variant_effect_dir": "gs://bucket/lead_variant_effect",
            "az_burden_binary_dir": "gs://bucket/az_binary",
            "az_burden_quantitative_dir": "gs://bucket/az_quantitative",
            "target_data_dir": "gs://bucket/target",
        },
    )()

    monkeypatch.setattr(prep, "log", messages.append)

    prep.log_step1_optional_input_modes(args)

    assert messages == [
        "SOURCE lead_variant_effect: enabled -> rescaled statistics from gs://bucket/lead_variant_effect",
        "SOURCE raw AZ burden: enabled -> binary=gs://bucket/az_binary quantitative=gs://bucket/az_quantitative target=gs://bucket/target",
    ]


def test_validate_raw_az_burden_args_requires_complete_pair_and_target():
    missing_quant = type(
        "Args",
        (),
        {
            "az_burden_binary_dir": "binary",
            "az_burden_quantitative_dir": "",
            "target_data_dir": "target",
        },
    )()
    missing_target = type(
        "Args",
        (),
        {
            "az_burden_binary_dir": "binary",
            "az_burden_quantitative_dir": "quant",
            "target_data_dir": "",
        },
    )()

    with pytest.raises(ValueError, match="requires both"):
        prep.validate_raw_az_burden_args(missing_quant)
    with pytest.raises(ValueError, match="requires --target_data_dir"):
        prep.validate_raw_az_burden_args(missing_target)


def test_optional_step1_debug_outputs_only_write_when_explicit(monkeypatch):
    calls = []

    monkeypatch.setattr(prep, "write_parquet_dataset", lambda df, output: calls.append(output))

    prep._maybe_write_section_output("gwas_variants", object(), "")
    prep._maybe_write_section_output("gwas_variants", object(), "debug.parquet")

    assert calls == ["debug.parquet"]


def test_build_analysis_ready_dataset_consumes_dataframes_directly():
    params = inspect.signature(prep.build_analysis_ready_dataset).parameters
    source = inspect.getsource(prep.build_analysis_ready_dataset)

    for name in ["spark", "args", "study_raw", "coloc", "gwas", "qtl", "coding", "clinvar", "burden", "az_ready"]:
        assert name in params

    assert "spark.read.parquet(args.coloc_output)" not in source
    assert "spark.read.parquet(args.gwas_output)" not in source
    assert "spark.read.parquet(args.qtl_output)" not in source
    assert "spark.read.parquet(args.burden_output)" not in source
    assert "spark.read.parquet(args.az_output)" not in source


def test_raw_az_burden_filters_to_requested_target_symbols_before_efo_join():
    source = inspect.getsource(prep._build_raw_az_burden_output_df)

    assert "requested_target_map" in source
    assert '.join(F.broadcast(requested_target_map), "gene_symbol_norm", "inner")' in source
    assert source.index('.join(F.broadcast(requested_target_map), "gene_symbol_norm", "inner")') < source.index(".join(efo_map,")


def test_raw_az_burden_efo_map_keeps_all_phenotype_efo_pairs():
    source = inspect.getsource(prep._load_az_burden_efo_map)

    assert '.dropDuplicates(["diseaseFromSource", "diseaseId"])' in source
    assert '.groupBy("diseaseFromSource")' not in source
    assert 'F.min(F.col("diseaseId"))' not in source


def test_final_burden_dedup_resolves_gene_disease_rows_by_pvalue():
    source = inspect.getsource(prep.build_analysis_ready_dataset)

    assert 'Window.partitionBy("targetId", "disease_key").orderBy(*burden_order)' in source
    assert 'F.col("pValueExponent").asc_nulls_last()' in source
    assert 'F.col("pValueMantissa").asc_nulls_last()' in source


def test_main_does_not_persist_large_raw_step1_inputs():
    source = inspect.getsource(prep.main)

    assert "_cache_frame(spark.read.parquet(resolve_parquet_input(args.credible_set_dir)))" not in source
    assert "_cache_frame(spark.read.parquet(resolve_parquet_input(args.study_data_dir)))" not in source
    assert "_cache_frame(\n            spark.read.parquet(resolve_parquet_input(args.variant_data_dir))" not in source
    assert 'spark.read.parquet(resolve_parquet_input(args.variant_data_dir))' in source
    assert 'F.col("variantId").cast("string").alias("ot_variant")' in source


def test_analysis_ready_uses_vidra2_style_common_and_coding_dedup():
    source = inspect.getsource(prep.build_analysis_ready_dataset)

    assert 'Window.partitionBy("variant", "as_gene", "as_disease", "GsourceLab", "GqtlLab").orderBy(' in source
    assert 'F.col("_gwasPValue").asc_nulls_last()' in source
    assert 'F.col("_qtlPValue").asc_nulls_last()' in source
    assert "common_candidates" in source
    assert "coding_candidates" in source


def test_study_disease_mapping_supports_lead_variant_fallback_disease_ids():
    source = inspect.getsource(prep._study_disease_mapping)

    assert '"traitFromSourceMappedIds"' in source
    assert '"diseaseIds"' in source
    assert '"left_anti"' in source


def test_lead_variant_effect_loader_uses_rescaled_beta_and_standard_error_fallbacks():
    source = inspect.getsource(prep._normalise_lead_variant_effect_df)

    assert "directionOfEffect" in source
    assert "absEstimatedBeta" in source
    assert "estimatedSE" in source
    assert "absZScore" in source
    assert "originalBeta" in source
    assert "originalStandardError" in source
    assert 'F.col("pValue") <= F.lit(LEAD_VARIANT_EFFECT_PVALUE_THRESHOLD)' in source
    assert 'F.abs(F.col("beta")) <= F.lit(LEAD_VARIANT_EFFECT_BETA_ABS_MAX)' in source
    assert 'F.col("alleleFrequency").isNull() | (F.col("alleleFrequency") != F.lit(0.0))' in source


def test_coding_transcript_consequence_rule_uses_gentropy_cutoff_and_protein_coding():
    source = inspect.getsource(prep._coding_transcript_consequence_is_eligible)

    assert "1 - (23/41)" in source
    assert "prep_coding_TRANSCRIPT_CONSEQUENCE_SCORE_MIN" in source
    assert "prep_coding_TRANSCRIPT_CONSEQUENCE_PROTEIN_BIOTYPE" in source
    assert 'transcript_consequence["consequenceScore"]' in source
    assert 'transcript_consequence["biotype"]' in source


def test_release_coding_path_uses_transcript_consequences_for_gene_assignment():
    source = inspect.getsource(prep._build_coding_output_df)

    assert '"transcriptConsequences"' in source
    assert "F.filter(" in source
    assert "_candidateTranscriptConsequences" in source
    assert "_coding_transcript_consequence_is_eligible" in source
    assert 'F.explode("_candidateTranscriptConsequences")' in source
    assert '_transcriptConsequence.targetId' in source
    assert 'F.max("consequenceScore").over(Window.partitionBy("variantId"))' in source
    assert "_variantEffect.targetId" not in source


def test_release_coding_path_collapses_same_target_transcripts_but_keeps_max_score_ties():
    source = inspect.getsource(prep._build_coding_output_df)

    assert '.filter(F.col("consequenceScore") == F.col("maxConsequenceScore"))' in source
    assert '.select("variantId", "geneId", "mostSevereConsequenceId", "variantDescription")' in source
    assert ".dropDuplicates()" in source


def test_lead_variant_effect_loader_exposes_transcript_consequence_target_fields():
    source = inspect.getsource(prep._normalise_lead_variant_effect_df)

    assert 'leadVariantConsequence.mostSevereConsequence.transcriptConsequence.targetId' in source
    assert 'leadVariantConsequence.mostSevereConsequence.transcriptConsequence.biotype' in source
    assert 'leadVariantConsequence.mostSevereConsequence.transcriptConsequence.consequenceScore' in source
    assert '"leadTranscriptTargetId"' in source
    assert '"leadTranscriptBiotype"' in source
    assert '"leadTranscriptConsequenceScore"' in source


def test_lead_coding_path_uses_lead_transcript_consequence_target_assignment():
    source = inspect.getsource(prep._build_coding_output_df_from_lead)

    assert '"leadTranscriptTargetId"' in source
    assert '"leadTranscriptBiotype"' in source
    assert '"leadTranscriptConsequenceScore"' in source
    assert "_leadTranscriptConsequence" in source
    assert "_coding_transcript_consequence_is_eligible" in source
    assert 'F.col("leadTranscriptTargetId").cast("string")' in source
    assert "_variantEffect.targetId" not in source


def test_standalone_prep_coding_entrypoint_delegates_to_shared_helper():
    source = inspect.getsource(prep.prep_coding_main)

    assert "_build_coding_output_df(" in source
    assert "_variantEffect.targetId" not in source


def test_build_coloc_output_df_from_lead_uses_locus_keys_from_raw_coloc():
    source = inspect.getsource(prep._build_coloc_output_df_from_lead)

    assert 'coloc_filtered_df.select("rightStudyLocusId").dropDuplicates(["rightStudyLocusId"])' in source
    assert 'join(right_loci_df, ["rightStudyLocusId"], "inner")' in source
    assert 'coloc_with_right_df = coloc_filtered_df.join(' in source
    assert '["rightStudyLocusId"],' in source
    assert 'coloc_with_right_df.select("leftStudyLocusId").dropDuplicates(["leftStudyLocusId"])' in source
    assert 'dropDuplicates(["leftStudyLocusId"])' in source
    assert 'join(left_loci_df, ["leftStudyLocusId"], "inner")' in source


def test_main_supports_optional_lead_variant_effect_branch():
    source = inspect.getsource(prep.main)

    assert "log_step1_optional_input_modes(args)" in source
    assert "if lead_variant_effect_inputs_enabled(args):" in source
    assert 'spark.read.parquet(resolve_parquet_input(args.lead_variant_effect_dir))' in source
    assert "_normalise_lead_variant_effect_df(lead_variant_effect_raw)" in source
    assert "_build_coloc_output_df_from_lead(" in source
    assert "_build_gwas_output_df_from_lead(" in source
    assert "_build_coding_output_df_from_lead(" in source
    assert 'coding_df = _cache_frame(' in source
    assert "StorageLevel.DISK_ONLY,\n            )\n        else:" in source
    assert "_build_qtl_output_df_from_lead(" in source


def test_clinvar_uses_vidra2_first_significance_without_exploding():
    source = inspect.getsource(prep.prep_clinvar_add_clinical_significance_column)

    assert ".getItem(0)" in source
    assert "F.explode" not in source


def test_analysis_ready_deduplicates_clinvar_by_source_key():
    source = inspect.getsource(prep.build_analysis_ready_dataset)

    assert "clinvar_candidates" in source
    assert "clinvar_window" in source
    assert 'F.col("_clinvarScore").desc_nulls_last()' in source
    assert 'F.col("as_clinicalSignificance").desc_nulls_last()' in source


def test_az_ready_deduplicates_by_lowest_pvalue():
    source = inspect.getsource(prep._deduplicate_az_ready)

    assert 'Window.partitionBy("variant", "as_gene", "as_disease", "GsourceLab", "GqtlLab").orderBy(' in source
    assert 'F.col("pValue").asc_nulls_last()' in source
    assert 'F.col("oddsRatio").desc_nulls_last()' in source


def test_both_az_builders_use_shared_az_dedup_helper():
    legacy_source = inspect.getsource(prep._build_az_ready_dataset)
    current_source = inspect.getsource(prep._build_az_ready_dataset_from_inputs)

    assert "return _deduplicate_az_ready(az_ready)" in legacy_source
    assert "return _deduplicate_az_ready(az_ready)" in current_source


def test_standalone_az_builder_always_uses_release_variant_ids_for_matching():
    source = inspect.getsource(prep._build_az_ready_dataset)

    assert 'spark.read.parquet(resolve_parquet_input(args.variant_data_dir))' in source
    assert 'F.col("variantId").cast("string").alias("ot_variant")' in source
    assert "lead_variant_effect_dir" not in source


def test_main_passes_release_variant_index_to_az_builder_in_all_modes():
    source = inspect.getsource(prep.main)

    assert "az_variant_index = None" in source
    assert 'az_variant_index = (\n                spark.read.parquet(resolve_parquet_input(args.variant_data_dir))' in source
    assert 'az_variant_index = (\n                variant_projection.select(F.col("variantId").alias("ot_variant"))' in source
    assert "_build_az_ready_dataset_from_inputs(" in source
    assert "az_variant_index," in source
