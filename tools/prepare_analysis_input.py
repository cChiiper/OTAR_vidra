#!/usr/bin/env python3
"""VIDRA Step 1: prepare analysis inputs (single independent script).

Step 1 reads these input families:
- colocalisation parquet
  Provides coloc links between GWAS loci and QTL loci. It is filtered to
  `h4 > colocalisation_threshold`, restricted to allowed right-side study types
  (eQTL/pQTL), and optionally to blood biosamples.
- GWAS/QTL/coding summary statistics
  Default mode reads release `credible_set + study + variant`.
  Optional lead mode reads `lead_variant_effect` instead and reconstructs the
  same internal contract from rescaled statistics.
- burden evidence
  Always reads release `evidence_gene_burden`. Optional raw AZ burden can be
  added on top and then collapsed together with release burden by best p-value.
- ClinVar/EVA evidence
  Reads `evidence_eva`, filters to datasource `eva`, excludes low-confidence
  records, and keeps only requested genes.
- AZ rare-variant file
  Reads the CSV-style AZ rare-variant input, filters to allelic model and
  genome-wide significance, maps HGNC gene symbols to Ensembl and AZ
  phenotypes to EFO, then checks whether the variant exists in the active
  release variant universe.

Step 1 always writes the same public handoff outputs:
- `vidra_analysis_ready/`
- `vidra_analysis_ready_manifest/`
- `GWAS_coding_variants_from_cs.parquet`
- `clinvar_variants.parquet`

Important runtime switches:
- `--lead_variant_effect_dir`
  If set, GWAS/QTL/coding summary statistics come from the rescaled
  `lead_variant_effect` dataset instead of the release
  `credible_set + study + variant` path.
- `--az_burden_binary_dir` + `--az_burden_quantitative_dir`
  If both are set, raw AZ burden evidence is added on top of the release
  `evidence_gene_burden` baseline. The release burden source remains required.

Final analysis-ready schema semantics:
- `variant`, `as_gene`, `as_disease`
  The public analysis key used downstream by Step 3.
- `yc`, `ycse`
  Variant-level effect and SE for GWAS-like evidence.
  For common coloc rows this is the GWAS beta/SE.
  For coding GWAS rows this is the coding GWAS beta/SE.
  For AZ rare-variant rows this is `log(odds ratio)` and its SE proxy from the
  OR confidence interval.
  ClinVar rows keep `yc=0` because their signal is not an effect size.
- `xc`, `xcse`
  QTL effect and SE.
  Only common coloc rows carry real QTL values here; coding, ClinVar and AZ
  rare-variant rows use the default `xc=0`, `xcse=0.1`.
- `bO`, `bOse`
  Gene-level burden effect and SE, joined onto every variant-level row by
  `(as_gene, as_disease)`. These come from release burden evidence, optionally
  enriched with raw AZ burden, then collapsed to one best row per gene/disease.
- `as_clinicalSignificance`
  ClinVar pathogenicity-like score on [0, 1].
  Only ClinVar rows carry a real value here; all other sources fill 0.

High-level merge scheme:
- `common_ready`
  Coloc GWAS + QTL rows. Deduplicated per
  `(variant, as_gene, as_disease, GsourceLab, GqtlLab)` using the lowest GWAS
  p-value, then lowest QTL p-value, then stable locus keys.
- `coding_ready`
  Non-coloc coding GWAS rows. Deduplicated on the same source-aware key using
  the lowest GWAS p-value, then stable locus key.
- `clinvar_ready`
  ClinVar/EVA rows. Deduplicated on the same source-aware key using EVA
  `score` first, then VIDRA clinical significance rank.
- `az_ready`
  AZ rare-variant rows. Deduplicated on the same source-aware key using the
  lowest p-value, then largest odds ratio.
- `burden_for_join`
  Gene/disease burden rows. If raw AZ burden is enabled it is unioned with the
  release burden rows, then one row per `(targetId, disease)` is kept by lowest
  p-value, with release evidence preferred on exact ties.

The goal of this file is to keep those decisions explicit. Comments therefore
focus on selection, fallback, and dedup logic rather than restating basic
Spark operations.
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
import tempfile
from pathlib import Path

try:
    from pyspark import StorageLevel
    from pyspark.sql import SparkSession, Window
    from pyspark.sql import functions as F
    from pyspark.sql import types as T
except ModuleNotFoundError:
    StorageLevel = None
    SparkSession = None
    Window = None
    F = None
    T = None

# ============================================================================
# Parameters
# ============================================================================
DEFAULT_COLOC_OUTPUT = "colocalising_variants.parquet"
DEFAULT_GWAS_OUTPUT = "GWAS_variants.parquet"
DEFAULT_QTL_OUTPUT = "QTL_variants.parquet"
DEFAULT_BURDEN_OUTPUT = "burden_tests.parquet"
DEFAULT_CODING_OUTPUT = "GWAS_coding_variants_from_cs.parquet"
DEFAULT_CLINVAR_OUTPUT = "clinvar_variants.parquet"
DEFAULT_AZ_OUTPUT = "AZ_variants.parquet"
DEFAULT_ANALYSIS_READY_DIR = "vidra_analysis_ready"
DEFAULT_ANALYSIS_MANIFEST_DIR = "vidra_analysis_ready_manifest"

# Final public schema consumed by Step 3.
# Source meaning:
# - GsourceLab=0 common coloc, 1 AZ rare variant, 2 ClinVar, 3 coding GWAS
# - GqtlLab=0 eQTL/common, 1 pQTL/common, 2 no-QTL source
#
# Signal meaning:
# - yc/ycse: GWAS-like or rare-variant effect channel
# - xc/xcse: QTL effect channel
# - bO/bOse: burden effect channel joined by gene+disease
# - as_clinicalSignificance: ClinVar-only pathogenicity-like score
ANALYSIS_OUTPUT_COLS = [
    "variant", "as_gene", "as_disease",
    "GsourceLab", "GqtlLab",
    "yc", "ycse", "xc", "xcse", "bO", "bOse",
    "as_clinicalSignificance",
]

ANALYSIS_SOURCE_COLS = [
    "variant", "as_gene", "as_disease",
    "GsourceLab", "GqtlLab",
    "yc", "ycse", "xc", "xcse",
    "as_clinicalSignificance",
]

ANALYSIS_DEFAULT_FILL = {
    "as_clinicalSignificance": 0.0,
    "bO": 0.0,
    "bOse": 2.0,
    "ycse": 0.14,
    "yc": 0.0,
    "xc": 0.0,
    "xcse": 0.1,
}

# ClinVar labels are mapped to [0, 1] by index / 16. Examples:
# not provided=0.0, uncertain significance=9/16, likely pathogenic=15/16,
# pathogenic=1.0. The order mirrors the VIDRA_2 ordinal mapping.
CLINICAL_SIG_ORDER = [
    "not provided",
    "association not found",
    "other",
    "benign",
    "likely benign",
    "low penetrance",
    "confers sensitivity",
    "uncertain risk allele",
    "drug response",
    "uncertain significance",
    "association",
    "affects",
    "likely risk allele",
    "risk factor",
    "established risk allele",
    "likely pathogenic",
    "pathogenic",
]

AZ_PVAL_THRESHOLD = 5e-8
LEAD_VARIANT_EFFECT_PVALUE_THRESHOLD = 5e-8
LEAD_VARIANT_EFFECT_BETA_ABS_MAX = 3.0
AZ_BURDEN_COLLAPSING_MODELS = ["ptv", "ptv5pcnt", "ptvraredmg"]

if T is not None:
    # The 470k AZ rare-variant file keeps the legacy analysis columns used by
    # this pipeline and adds trailing Case/Control MAF columns that we ignore
    # downstream. We still include them here so Spark header validation matches
    # the on-disk CSV exactly and does not silently rely on positional parsing.
    AZ_SCHEMA = T.StructType([
        T.StructField("Variant", T.StringType(), True),
        T.StructField("Variant type", T.StringType(), True),
        T.StructField("Phenotype", T.StringType(), True),
        T.StructField("Category", T.StringType(), True),
        T.StructField("Model", T.StringType(), True),
        T.StructField("Consequence type", T.StringType(), True),
        T.StructField("Gene", T.StringType(), True),
        T.StructField("Transcript", T.StringType(), True),
        T.StructField("cDNA change", T.StringType(), True),
        T.StructField("Amino acid change", T.StringType(), True),
        T.StructField("Exon rank", T.StringType(), True),
        T.StructField("No. cases", T.DoubleType(), True),
        T.StructField("No. AA cases", T.DoubleType(), True),
        T.StructField("No. AB cases", T.DoubleType(), True),
        T.StructField("No. BB cases", T.DoubleType(), True),
        T.StructField("Case AAF", T.DoubleType(), True),
        T.StructField("% AB or BB cases", T.DoubleType(), True),
        T.StructField("% BB cases", T.DoubleType(), True),
        T.StructField("No. controls", T.DoubleType(), True),
        T.StructField("No. AA controls", T.DoubleType(), True),
        T.StructField("No. AB controls", T.DoubleType(), True),
        T.StructField("No. BB controls", T.DoubleType(), True),
        T.StructField("Control AAF", T.DoubleType(), True),
        T.StructField("% AB or BB controls", T.DoubleType(), True),
        T.StructField("% BB controls", T.DoubleType(), True),
        T.StructField("p-value", T.DoubleType(), True),
        T.StructField("Odds ratio", T.DoubleType(), True),
        T.StructField("Odds ratio LCI", T.DoubleType(), True),
        T.StructField("Odds ratio UCI", T.DoubleType(), True),
        T.StructField("Case MAF", T.DoubleType(), True),
        T.StructField("Control MAF", T.DoubleType(), True),
    ])
else:
    AZ_SCHEMA = None


def log(message: str) -> None:
    print(f"[prepare_analysis_input] {message}", flush=True)


def require_pyspark() -> None:
    if StorageLevel is None or SparkSession is None or Window is None or F is None or T is None:
        raise ModuleNotFoundError(
            "pyspark is required to run prepare_analysis_input.py. "
            "Use the VIDRA_5 runtime image or install pyspark locally."
        )


def normalize_az_variant_id(variant: str | None) -> str | None:
    if variant is None:
        return None

    cleaned = str(variant).strip().replace(" ", "_").replace("-", "_")
    parts = cleaned.split("_")
    if len(parts) != 4 or any(not part for part in parts):
        return None
    return "_".join(parts)


def burden_row_is_allowed(project_id: str | None, cohort_id: str | None, statistical_method: str | None) -> bool:
    return (
        project_id in prep_burden_ALLOWED_PROJECT_IDS
        and cohort_id in prep_burden_ALLOWED_COHORT_IDS
        and statistical_method in prep_burden_ALLOWED_STAT_METHODS
    )


def raw_az_burden_inputs_enabled(args) -> bool:
    return bool(str(getattr(args, "az_burden_binary_dir", "") or "").strip()) and bool(
        str(getattr(args, "az_burden_quantitative_dir", "") or "").strip()
    )


def lead_variant_effect_inputs_enabled(args) -> bool:
    return bool(str(getattr(args, "lead_variant_effect_dir", "") or "").strip())


def validate_raw_az_burden_args(args) -> None:
    binary_dir = str(getattr(args, "az_burden_binary_dir", "") or "").strip()
    quantitative_dir = str(getattr(args, "az_burden_quantitative_dir", "") or "").strip()
    target_data_dir = str(getattr(args, "target_data_dir", "") or "").strip()

    if not binary_dir and not quantitative_dir:
        return

    if not binary_dir or not quantitative_dir:
        raise ValueError(
            "Raw AZ burden enrichment requires both --az_burden_binary_dir "
            "and --az_burden_quantitative_dir, or neither."
        )

    if not target_data_dir:
        raise ValueError("Raw AZ burden enrichment requires --target_data_dir.")


def log_step1_optional_input_modes(args) -> None:
    if lead_variant_effect_inputs_enabled(args):
        log(
            "SOURCE lead_variant_effect: enabled -> "
            f"rescaled statistics from {args.lead_variant_effect_dir}"
        )
    else:
        log(
            "SOURCE lead_variant_effect: disabled -> "
            "using release credible_set + study + variant inputs"
        )

    if raw_az_burden_inputs_enabled(args):
        log(
            "SOURCE raw AZ burden: enabled -> "
            f"binary={args.az_burden_binary_dir} "
            f"quantitative={args.az_burden_quantitative_dir} "
            f"target={args.target_data_dir}"
        )
    else:
        log(
            "SOURCE raw AZ burden: disabled -> "
            "using only release gene_burden evidence"
        )


def is_excluded_clinvar_confidence(confidence: str | None) -> bool:
    return str(confidence or "").strip().lower() in {
        value.lower() for value in prep_clinvar_EXCLUDED_CONFIDENCE
    }


def is_allowed_right_study_type(study_type: str | None) -> bool:
    return str(study_type or "").strip().lower() in set(prep_coloc_ALLOWED_RIGHT_STUDY_TYPES)


def is_allowed_blood_biosample_id(biosample_id: str | None) -> bool:
    return str(biosample_id or "").strip() in set(prep_coloc_BLOOD_BIOSAMPLE_IDS)


def normalise_most_severe_consequence_id(value: str | None) -> str:
    return str(value or "").upper().replace(":", "_")


def is_allowed_coding_consequence(value: str | None) -> bool:
    return normalise_most_severe_consequence_id(value) in set(prep_coding_CODING_MOST_SEVERE_CONSEQUENCE_IDS)


def run_inlined_main(section_name: str, main_fn, argv: list[str], expected_output: str) -> None:
    old_argv = sys.argv[:]
    try:
        sys.argv = [section_name, *argv]
        log(f"START {section_name}")
        main_fn()
    finally:
        sys.argv = old_argv

    if str(expected_output).startswith("gs://"):
        log(f"DONE  {section_name} -> {expected_output}")
        return

    output = Path(expected_output)
    if not output.exists():
        raise RuntimeError(f"{section_name} completed but output is missing: {output}")
    log(f"DONE  {section_name} -> {output}")



# ============================================================================
# Inlined Module: spark_pipeline_utils.py
# ============================================================================

"""
Shared Spark/parquet utilities for VIDRA ingest scripts.

These helpers keep IO/session behavior consistent across the per-source Spark jobs.
"""
def resolve_parquet_input(path_like):
    """
    Resolve input that can be:
      - a parquet directory
      - a single parquet file
      - a glob pattern
    """
    if not path_like:
        raise ValueError("Parquet path is empty.")

    if str(path_like).startswith("gs://"):
        return path_like

    if any(token in path_like for token in ["*", "?", "["]):
        if not glob.glob(path_like):
            raise FileNotFoundError("No files matched parquet glob: {0}".format(path_like))
        return path_like

    if os.path.isdir(path_like):
        parquet_glob = os.path.join(path_like.rstrip("/"), "*.parquet")
        if not glob.glob(parquet_glob):
            raise FileNotFoundError("No parquet files found under directory: {0}".format(path_like))
        return parquet_glob

    if os.path.isfile(path_like):
        if not path_like.endswith(".parquet"):
            raise ValueError("Input file is not a parquet file: {0}".format(path_like))
        return path_like

    raise FileNotFoundError("Path not found: {0}".format(path_like))


def _gcsfs_kwargs(gcp_project: str | None = None) -> dict[str, object]:
    kwargs: dict[str, object] = {}
    if gcp_project:
        kwargs["project"] = str(gcp_project)
        kwargs["requester_pays"] = True
    return kwargs


def path_exists(path_like, gcp_project: str | None = None):
    if not path_like:
        return False
    if str(path_like).startswith("gs://"):
        try:
            import gcsfs
        except ImportError as exc:
            raise RuntimeError("Checking gs:// paths requires gcsfs.") from exc
        fs = gcsfs.GCSFileSystem(**_gcsfs_kwargs(gcp_project))
        return fs.exists(str(path_like)[len("gs://"):])
    return Path(path_like).exists()


def _assert_input_path_readable(path_like, label: str, gcp_project: str | None = None) -> None:
    text = str(path_like or "").strip()
    if not text:
        raise ValueError(f"{label} is empty.")

    if not path_exists(text, gcp_project=gcp_project):
        raise FileNotFoundError(f"{label} not found or unreadable: {text}")

    if text.startswith("gs://"):
        return

    path = Path(text)
    if path.is_file() and not os.access(path, os.R_OK):
        raise PermissionError(f"{label} is not readable: {text}")
    if path.is_dir() and not os.access(path, os.R_OK | os.X_OK):
        raise PermissionError(f"{label} directory is not readable: {text}")


def validate_step1_input_paths(args) -> None:
    gcp_project = str(getattr(args, "gcp_project", "") or "").strip()
    required_inputs = [
        ("Gene list", args.genes),
        ("Colocalisation input", args.coloc_data_dir),
        ("Credible set input", args.credible_set_dir),
        ("Study input", args.study_data_dir),
        ("Variant input", args.variant_data_dir),
        ("Burden evidence input", args.burden_evidence_dir),
        ("ClinVar evidence input", args.clinvar_evidence_dir),
        ("AZ rare-variant input", args.az_variants_file),
        ("AZ phenotype mapping input", args.az_mapping_file),
        ("AZ gene mapping input", args.az_gene_map_file),
    ]

    if lead_variant_effect_inputs_enabled(args):
        required_inputs.append(("Lead variant effect input", args.lead_variant_effect_dir))

    if raw_az_burden_inputs_enabled(args):
        required_inputs.extend(
            [
                ("Raw AZ burden binary input", args.az_burden_binary_dir),
                ("Raw AZ burden quantitative input", args.az_burden_quantitative_dir),
                ("Target data input", args.target_data_dir),
            ]
        )

    for label, path_like in required_inputs:
        _assert_input_path_readable(path_like, label, gcp_project=gcp_project)


def ensure_columns_present(df, required_columns, dataset_label):
    """
    Validate presence of required columns in a source DataFrame.
    """
    missing = [c for c in required_columns if c not in df.columns]
    if missing:
        raise ValueError("{0} parquet is missing required columns: {1}".format(dataset_label, missing))


def nested_field_exists(schema, field_path: str) -> bool:
    current = schema
    for part in field_path.split("."):
        if not isinstance(current, T.StructType) or part not in current.fieldNames():
            return False
        current = current[part].dataType
    return True


def ensure_nested_fields_present(df, required_paths, dataset_label):
    missing = [path for path in required_paths if not nested_field_exists(df.schema, path)]
    if missing:
        raise ValueError("{0} parquet is missing required nested fields: {1}".format(dataset_label, missing))


def optional_col(col_name, schema, dtype):
    """
    Return a typed Spark column if present in schema, else NULL.
    """
    if col_name not in schema.fieldNames():
        return F.lit(None).cast(dtype)
    return F.col(col_name).cast(dtype)


def optional_nested_col(field_path, schema, dtype):
    if not nested_field_exists(schema, field_path):
        return F.lit(None).cast(dtype)
    return F.col(field_path).cast(dtype)


def write_parquet_dataset(df, output_path, n_partitions=None):
    """
    Write parquet dataset to a deterministic output path.

    If n_partitions is None, keep Spark's current partitioning. This matches the
    VIDRA_2 style of only forcing partition layout when the downstream consumer
    explicitly benefits from it.
    """
    if not str(output_path).startswith("gs://") and os.path.exists(output_path):
        if os.path.isdir(output_path):
            shutil.rmtree(output_path)
        else:
            os.remove(output_path)

    writer_df = df if n_partitions is None else df.coalesce(n_partitions)
    writer_df.write.mode("overwrite").parquet(output_path)


SPARK_SESSION_DEFAULTS = {
    # Blank values mean "use the ambient Spark environment", which is required
    # for Dataproc Serverless. Local / Batch single-node runs pass explicit
    # values via Nextflow config.
    "master": "",
    "driver_memory": "",
    "shuffle_partitions": None,
    "default_parallelism": None,
    "gcp_project": "",
    "requester_pays_buckets": "",
}


def configure_spark_defaults(
    *,
    master: str | None = None,
    driver_memory: str | None = None,
    shuffle_partitions: int | None = None,
    default_parallelism: int | None = None,
    gcp_project: str | None = None,
    requester_pays_buckets: str | None = None,
):
    if master is not None:
        SPARK_SESSION_DEFAULTS["master"] = master
    if driver_memory is not None:
        SPARK_SESSION_DEFAULTS["driver_memory"] = driver_memory
    if shuffle_partitions is not None:
        SPARK_SESSION_DEFAULTS["shuffle_partitions"] = int(shuffle_partitions)
    if default_parallelism is not None:
        SPARK_SESSION_DEFAULTS["default_parallelism"] = int(default_parallelism)
    if gcp_project is not None:
        SPARK_SESSION_DEFAULTS["gcp_project"] = str(gcp_project or "")
    if requester_pays_buckets is not None:
        SPARK_SESSION_DEFAULTS["requester_pays_buckets"] = str(requester_pays_buckets or "")


def build_local_spark(
    app_name,
    master=None,
    driver_memory=None,
    shuffle_partitions=None,
    default_parallelism=None,
    gcp_project=None,
    requester_pays_buckets=None,
):
    """
    Build a Spark session for either single-node local execution or ambient
    cluster-managed execution.

    When `master` is blank, Spark uses the surrounding environment (for example
    Dataproc Serverless). When `master` is set to `local[N]`, Spark runs
    single-node inside the current container / VM.
    """
    require_pyspark()
    master = SPARK_SESSION_DEFAULTS["master"] if master is None else master
    driver_memory = SPARK_SESSION_DEFAULTS["driver_memory"] if driver_memory is None else driver_memory
    shuffle_partitions = SPARK_SESSION_DEFAULTS["shuffle_partitions"] if shuffle_partitions is None else shuffle_partitions
    default_parallelism = SPARK_SESSION_DEFAULTS["default_parallelism"] if default_parallelism is None else default_parallelism
    gcp_project = SPARK_SESSION_DEFAULTS["gcp_project"] if gcp_project is None else gcp_project
    requester_pays_buckets = SPARK_SESSION_DEFAULTS["requester_pays_buckets"] if requester_pays_buckets is None else requester_pays_buckets

    if master and str(master).startswith("local"):
        runtime_dirs = configure_local_spark_runtime()
        os.environ["HOME"] = runtime_dirs["home"]
        os.environ["SPARK_LOCAL_DIRS"] = runtime_dirs["spark_local"]

        submit_opts = os.environ.get("SPARK_SUBMIT_OPTS", "").strip()
        required_opts = [
            f"-Duser.home={runtime_dirs['home']}",
            f"-Divy.home={runtime_dirs['ivy']}",
            f"-Divy.cache.dir={runtime_dirs['ivy_cache']}",
        ]
        missing_opts = [opt for opt in required_opts if opt not in submit_opts]
        combined_opts = " ".join(missing_opts + ([submit_opts] if submit_opts else []))
        os.environ["SPARK_SUBMIT_OPTS"] = combined_opts.strip()

    builder = SparkSession.builder.appName(app_name)
    if master:
        builder = builder.master(master)

    if driver_memory:
        builder = builder.config("spark.driver.memory", driver_memory)
    if shuffle_partitions:
        builder = builder.config("spark.sql.shuffle.partitions", str(shuffle_partitions))

    builder = (
        builder
        .config("spark.sql.adaptive.enabled", "true")
        .config("spark.sql.adaptive.coalescePartitions.enabled", "true")
        .config("spark.sql.session.timeZone", "UTC")
    )

    if default_parallelism:
        builder = builder.config("spark.default.parallelism", str(default_parallelism))

    if master and str(master).startswith("local"):
        builder = (
            builder
            .config("spark.jars.ivy", runtime_dirs["ivy"])
            .config("spark.local.dir", runtime_dirs["spark_local"])
        )

    if gcp_project and requester_pays_buckets:
        builder = (
            builder
            .config("spark.hadoop.fs.gs.requester.pays.mode", "CUSTOM")
            .config("spark.hadoop.fs.gs.requester.pays.project.id", str(gcp_project))
            .config("spark.hadoop.fs.gs.requester.pays.buckets", str(requester_pays_buckets))
        )

    spark = builder.getOrCreate()
    print(
        "[prepare_analysis_input] Spark session "
        f"app={app_name} master={spark.sparkContext.master} "
        f"defaultParallelism={spark.sparkContext.defaultParallelism} "
        f"shufflePartitions={spark.conf.get('spark.sql.shuffle.partitions')}",
        flush=True,
    )
    return spark


def _coerce_absolute_path(path_like: str | Path) -> Path:
    path = Path(str(path_like)).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    return path


def configure_local_spark_runtime() -> dict[str, str]:
    """Return writable absolute runtime directories for local Spark launches.

    This matters for containerized local runs under a host UID/GID mapping,
    where Java sometimes resolves `user.home` to `?` and Ivy then crashes while
    trying to initialize `?/.ivy2/local`.
    """

    configured_home = str(os.environ.get("VIDRA_SPARK_USER_HOME", "")).strip()
    current_home = str(os.environ.get("HOME", "")).strip()

    if configured_home and configured_home != "?":
        home_dir = _coerce_absolute_path(configured_home)
    elif current_home and current_home != "?":
        home_dir = _coerce_absolute_path(current_home)
    else:
        home_dir = Path(tempfile.gettempdir()) / "vidra5_spark_home"

    ivy_dir = _coerce_absolute_path(os.environ.get("VIDRA_SPARK_IVY_DIR", home_dir / ".ivy2"))
    ivy_cache_dir = ivy_dir / "cache"
    spark_local_dir = _coerce_absolute_path(os.environ.get("VIDRA_SPARK_LOCAL_DIRS", home_dir / "spark_local"))

    for path in (home_dir, ivy_dir, ivy_cache_dir, spark_local_dir):
        path.mkdir(parents=True, exist_ok=True)

    return {
        "home": str(home_dir),
        "ivy": str(ivy_dir),
        "ivy_cache": str(ivy_cache_dir),
        "spark_local": str(spark_local_dir),
    }


def _gs_bucket_name(path_like: str | None) -> str:
    text = str(path_like or "").strip()
    if not text.startswith("gs://"):
        return ""
    return text[len("gs://"):].split("/", 1)[0].strip()


def load_gene_ids_from_csv(path, gcp_project: str | None = None):
    """
    Load IDs from first CSV column, tolerating optional header rows.
    """
    if str(path).startswith("gs://"):
        try:
            import gcsfs
        except ImportError as exc:
            raise RuntimeError(
                "Reading genes from gs:// requires gcsfs."
            ) from exc
        fs = gcsfs.GCSFileSystem(**_gcsfs_kwargs(gcp_project))
        handle = fs.open(path[len("gs://"):], "r")
    else:
        if not os.path.isfile(path):
            raise FileNotFoundError("Genes file not found: {0}".format(path))
        handle = open(path, "r", encoding="utf-8")

    genes = []
    with handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            value = line.split(",")[0].strip().strip('"').strip("'")
            if not value:
                continue

            lower = value.lower()
            if lower in {"targetid", "gene", "ensembl_id", "gene_id"}:
                continue
            genes.append(value)

    genes = sorted(set(genes))
    if not genes:
        raise ValueError("No genes found in input file.")
    return genes



# ============================================================================
# Inlined Module: get_coloc_vars_per_gene_spark.py
# ============================================================================
########################################################################################################
### Script to extract colocalising variants per gene using Spark and release parquet data.           ###
### Script written with GPT-5.3-Codex on the 16.02.2026                                              ###
########################################################################################################
prep_coloc_DEFAULT_OUTPUT_FILE = "colocalising_variants.parquet"
# Filtering switches (reintroduced from working/vidra pre-processing logic)
prep_coloc_FILTER_LEFT_STUDY_TYPE_TO_GWAS = True
prep_coloc_ALLOWED_RIGHT_STUDY_TYPES = ["eqtl", "pqtl"]
prep_coloc_FILTER_QTL_TO_BLOOD_BIOSAMPLE_IDS = False

# Blood biosample allowlist used for QTL branch filtering.
# Matches the legacy allowlist used in working/vidra.
prep_coloc_BLOOD_BIOSAMPLE_IDS = [
    "CL_0000084",      # T cell
    "CL_0000233",      # platelet
    "CL_0000576",      # monocyte
    "CL_0000775",      # neutrophil
    "UBERON_0000178",  # blood
    "UBERON_0001969",  # blood plasma
]

prep_coloc_COLOC_COLUMNS = [
    "leftStudyLocusId",
    "rightStudyLocusId",
    "rightStudyType",
    "chromosome",
    "colocalisationMethod",
    "numberColocalisingVariants",
    "h3",
    "h4",
    "clpp",
    "betaRatioSignAverage",
]

prep_coloc_CREDIBLE_SET_COLUMNS = [
    "studyLocusId",
    "studyType",
    "studyId",
    "variantId",
]

prep_coloc_STUDY_COLUMNS = [
    "studyId",
    "geneId",
    "biosampleId",
]


def prep_coloc_parse_args():
    """
    Parse CLI arguments.

    Args:
        argv[1]: genes CSV (1 gene per row, Ensembl IDs)
        argv[2]: colocalisation H4 threshold
        argv[3]: colocalisation parquet path (dir/file/glob)
        argv[4]: credible_set parquet path (dir/file/glob)
        argv[5]: study parquet path (dir/file/glob)
        argv[6]: output parquet dataset path (optional)
    """
    if len(sys.argv) < 6:
        raise ValueError(
            "Usage: python get_coloc_vars_per_gene_spark.py <genes.csv> <h4_threshold> "
            "<colocalisation_parquet> <credible_set_parquet> <study_parquet>"
        )

    genes_file = sys.argv[1]
    coloc_threshold = float(sys.argv[2])
    coloc_path = resolve_parquet_input(sys.argv[3])
    credible_set_path = resolve_parquet_input(sys.argv[4])
    study_path = resolve_parquet_input(sys.argv[5])

    output_file = sys.argv[6] if len(sys.argv) > 6 else prep_coloc_DEFAULT_OUTPUT_FILE

    return genes_file, coloc_threshold, coloc_path, credible_set_path, study_path, output_file


def prep_coloc_main():
    # 1) Parse CLI inputs and load requested genes.
    genes_file, coloc_threshold, coloc_path, credible_set_path, study_path, output_file = prep_coloc_parse_args()
    genes = load_gene_ids_from_csv(genes_file)

    # 2) Build Spark session and load parquet inputs.
    spark = build_local_spark("coloc_vars_per_gene_spark")
    try:
        coloc_raw = spark.read.parquet(coloc_path)
        credible_set_raw = spark.read.parquet(credible_set_path)
        study_raw = spark.read.parquet(study_path)

        # 3) Validate required source schema columns.
        ensure_columns_present(coloc_raw, prep_coloc_COLOC_COLUMNS, "Colocalisation")
        ensure_columns_present(credible_set_raw, prep_coloc_CREDIBLE_SET_COLUMNS, "Credible set")
        ensure_columns_present(study_raw, prep_coloc_STUDY_COLUMNS, "Study")

        # 4) Keep coloc signals passing threshold and configured right study types.
        coloc_df = coloc_raw.select(*prep_coloc_COLOC_COLUMNS)
        coloc_filtered_df = coloc_df.filter(F.col("h4") > F.lit(coloc_threshold))

        if prep_coloc_ALLOWED_RIGHT_STUDY_TYPES:
            allowed_right_types = [x.lower() for x in prep_coloc_ALLOWED_RIGHT_STUDY_TYPES]
            coloc_filtered_df = coloc_filtered_df.filter(
                F.lower(F.col("rightStudyType")).isin(allowed_right_types)
            )

        # 5) Restrict right-side studies to the input gene list.
        genes_df = spark.createDataFrame([(gene.upper(),) for gene in genes], ["rightGeneId"])

        study_filtered_df = (
            study_raw.select(*prep_coloc_STUDY_COLUMNS)
            .withColumn("rightGeneId", F.upper(F.col("geneId")))
            .join(F.broadcast(genes_df), "rightGeneId", "inner")
            .select(
                F.col("studyId").alias("rightStudyId"),
                F.col("rightGeneId"),
                F.col("biosampleId").alias("rightBiosampleId"),
            )
            .dropna(subset=["rightStudyId", "rightGeneId"])
            .dropDuplicates()
        )

        credible_set_df = credible_set_raw.select(*prep_coloc_CREDIBLE_SET_COLUMNS)

        # 6) Join right-side coloc loci to credible set and study metadata.
        right_loci_df = coloc_filtered_df.select("rightStudyLocusId").dropDuplicates(["rightStudyLocusId"])

        right_credible_set_df = (
            right_loci_df.join(
                credible_set_df.select(
                    F.col("studyLocusId").alias("rightStudyLocusId"),
                    F.col("studyId").alias("rightStudyId"),
                    F.col("variantId").alias("rightVariantId"),
                ),
                "rightStudyLocusId",
                "inner",
            )
            .join(F.broadcast(study_filtered_df), "rightStudyId", "inner")
        )

        coloc_with_right_df = coloc_filtered_df.join(right_credible_set_df, "rightStudyLocusId", "inner")

        # 7) Add left-side locus metadata from credible set.
        left_loci_df = coloc_with_right_df.select("leftStudyLocusId").dropDuplicates(["leftStudyLocusId"])
        left_credible_set_df = (
            credible_set_df.select(
                F.col("studyLocusId").alias("leftStudyLocusId"),
                F.col("studyType").alias("leftStudyType"),
                F.col("studyId").alias("leftStudyId"),
                F.col("variantId").alias("leftVariantId"),
            )
            .dropDuplicates(["leftStudyLocusId"])
            .join(left_loci_df, "leftStudyLocusId", "inner")
        )

        joined_df = coloc_with_right_df.join(left_credible_set_df, "leftStudyLocusId", "inner")

        # Reintroduced filtering from working/vidra.
        if prep_coloc_FILTER_LEFT_STUDY_TYPE_TO_GWAS:
            joined_df = joined_df.filter(F.lower(F.col("leftStudyType")) == F.lit("gwas"))

        if prep_coloc_FILTER_QTL_TO_BLOOD_BIOSAMPLE_IDS:
            joined_df = joined_df.filter(
                F.col("rightBiosampleId").isin(prep_coloc_BLOOD_BIOSAMPLE_IDS)
            )

        # 8) Select final output columns with stable order and deduplicate rows.
        output_df = (
            joined_df.select(
                # Core colocalisation columns.
                "leftStudyLocusId",
                "rightStudyLocusId",
                "rightStudyType",
                "chromosome",
                "colocalisationMethod",
                "numberColocalisingVariants",
                "h3",
                "h4",
                "clpp",
                "betaRatioSignAverage",
                # Additional columns.
                "leftStudyType",
                "leftStudyId",
                "leftVariantId",
                "rightStudyId",
                "rightBiosampleId",
                "rightVariantId",
                "rightGeneId",
            )
            .dropDuplicates()
        )

        # 9) Write deterministic parquet output for downstream modules.
        write_parquet_dataset(output_df, output_file)
    finally:
        try:
            spark.stop()
        except Exception:  # pragma: no cover - defensive shutdown path
            pass

# ============================================================================
# Inlined Module: get_GWAS_vars_spark.py
# ============================================================================
########################################################################################################
### Script to extract GWAS summary statistics for coloc lead variants using Spark parquet sources.   ###
### Script written with GPT-5-Codex on the 19.02.2026                                                ###
########################################################################################################
prep_gwas_DEFAULT_OUTPUT_FILE = "GWAS_variants.parquet"
prep_gwas_COLOC_REQUIRED_COLUMNS = [
    "leftStudyLocusId",
    "leftVariantId",
    "leftStudyType",
]

prep_gwas_CREDIBLE_SET_REQUIRED_COLUMNS = [
    "studyLocusId",
    "variantId",
    "beta",
    "pValueMantissa",
    "pValueExponent",
    "studyType",
    "studyId",
]


def prep_gwas_parse_args():
    """
    Parse CLI arguments.

    Args:
        argv[1]: coloc parquet dataset generated by get_coloc_vars_per_gene_spark.py
        argv[2]: credible_set parquet path (dir/file/glob)
        argv[3]: output parquet dataset path (optional)
    """
    if len(sys.argv) < 3:
        raise ValueError(
            "Usage: python get_GWAS_vars_spark.py <coloc_variants_parquet> <credible_set_parquet> [output_parquet]"
        )

    coloc_path = resolve_parquet_input(sys.argv[1])
    credible_set_path = resolve_parquet_input(sys.argv[2])
    output_path = sys.argv[3] if len(sys.argv) > 3 else prep_gwas_DEFAULT_OUTPUT_FILE

    return coloc_path, credible_set_path, output_path


def prep_gwas_main():
    # 1) Parse CLI inputs.
    coloc_path, credible_set_path, output_path = prep_gwas_parse_args()

    # 2) Build Spark session and load parquet inputs.
    spark = build_local_spark("gwas_vars_spark")
    try:
        coloc_raw = spark.read.parquet(coloc_path)
        credible_set_raw = spark.read.parquet(credible_set_path)

        # 3) Validate required source schema columns.
        ensure_columns_present(coloc_raw, prep_gwas_COLOC_REQUIRED_COLUMNS, "Colocalising variants")
        ensure_columns_present(credible_set_raw, prep_gwas_CREDIBLE_SET_REQUIRED_COLUMNS, "Credible set")

        credible_schema = credible_set_raw.schema

        # 4) Keep only left-side GWAS locus/variant keys from colocalisation output.
        coloc_gwas = (
            coloc_raw.select("leftStudyLocusId", "leftVariantId", "leftStudyType")
            .filter(F.lower(F.col("leftStudyType")) == F.lit("gwas"))
            .dropna(subset=["leftStudyLocusId", "leftVariantId"])
            .dropDuplicates(["leftStudyLocusId", "leftVariantId"])
        )

        # 5) Select required credible-set summary statistic columns.
        credible_selected = (
            credible_set_raw
            .select(
                F.col("studyLocusId").cast("string").alias("leftStudyLocusId"),
                F.col("variantId").cast("string").alias("variantId"),
                optional_col("beta", credible_schema, "double").alias("beta"),
                optional_col("standardError", credible_schema, "double").alias("standardError"),
                optional_col("zScore", credible_schema, "double").alias("zScore"),
                optional_col("pValueMantissa", credible_schema, "double").alias("pValueMantissa"),
                optional_col("pValueExponent", credible_schema, "long").alias("pValueExponent"),
                optional_col("studyType", credible_schema, "string").alias("studyType"),
                optional_col("studyId", credible_schema, "string").alias("studyId"),
            )
        )

        # 6) Match colocalisation locus/variant pairs to credible-set rows.
        joined = (
            coloc_gwas.join(
                credible_selected,
                on=(
                    (coloc_gwas["leftStudyLocusId"] == credible_selected["leftStudyLocusId"])
                    & (coloc_gwas["leftVariantId"] == credible_selected["variantId"])
                ),
                how="inner",
            )
            .select(
                credible_selected["leftStudyLocusId"].alias("studyLocusId"),
                credible_selected["variantId"],
                credible_selected["beta"],
                credible_selected["standardError"],
                credible_selected["zScore"],
                credible_selected["pValueMantissa"],
                credible_selected["pValueExponent"],
                credible_selected["studyType"],
                credible_selected["studyId"],
            )
            .dropDuplicates()
        )

        # 7) Fill missing SE from beta/zScore and keep complete effect-size rows only.
        out_df = (
            joined.withColumn(
                "standardError",
                F.coalesce(
                    F.col("standardError"),
                    F.when(
                        F.col("beta").isNotNull()
                        & F.col("zScore").isNotNull()
                        & (F.col("zScore") != F.lit(0.0)),
                        F.col("beta") / F.col("zScore"),
                    ),
                ),
            )
            # Keep only complete effect-size rows for downstream modelling.
            .filter(F.col("beta").isNotNull() & F.col("standardError").isNotNull())
            .drop("zScore")
        )

        # 8) Write deterministic parquet output for downstream modules.
        write_parquet_dataset(out_df, output_path)
    finally:
        try:
            spark.stop()
        except Exception:  # pragma: no cover - defensive shutdown path
            pass

# ============================================================================
# Inlined Module: get_coding_GWAS_NoColoc_spark.py
# ============================================================================
########################################################################################################
### Script to extract GWAS coding variants per gene from credible set and variant parquet sources.   ###
### Script written with GPT-5-Codex on the 19.02.2026                                                ###
########################################################################################################
prep_coding_DEFAULT_OUTPUT_FILE = "GWAS_coding_variants_from_cs.parquet"
prep_coding_GWAS_PVALUE_THRESHOLD = 5e-8

prep_coding_VARIANT_REQUIRED_COLUMNS = [
    "variantId",
    "transcriptConsequences",
    "mostSevereConsequenceId",
    "variantDescription",
]

prep_coding_CREDIBLE_SET_REQUIRED_COLUMNS = [
    "studyLocusId",
    "variantId",
    "beta",
    "pValueMantissa",
    "pValueExponent",
    "studyType",
    "studyId",
]

prep_coding_EXISTING_GWAS_REQUIRED_COLUMNS = [
    "variantId",
]

prep_coding_CODING_MOST_SEVERE_CONSEQUENCE_IDS = [
    "SO_0001821",  # inframe_insertion
    "SO_0001589",  # frameshift_variant
    "SO_0001587",  # stop_gained
    "SO_0001575",  # splice_donor_variant
    "SO_0001580",  # coding_sequence_variant
    "SO_0001578",  # stop_lost
    "SO_0001567",  # stop_retained_variant
    "SO_0001583",  # missense_variant
    "SO_0001626",  # incomplete_terminal_codon_variant
    "SO_0001818",  # protein_altering_variant
    "SO_0002012",  # start_lost
    "SO_0001574",  # splice_acceptor_variant
    "SO_0001822",  # inframe_deletion
]

prep_coding_TRANSCRIPT_CONSEQUENCE_SCORE_MIN = 1.0 - (23.0 / 41.0)
prep_coding_TRANSCRIPT_CONSEQUENCE_PROTEIN_BIOTYPE = "protein_coding"


def _coding_consequence_match_column(normalised_ids_col: str = "_normalisedConsequenceIds"):
    return F.coalesce(
        *[
            F.when(F.array_contains(F.col(normalised_ids_col), consequence_id), F.lit(consequence_id))
            for consequence_id in prep_coding_CODING_MOST_SEVERE_CONSEQUENCE_IDS
        ]
    )


def _coding_transcript_consequence_is_eligible(transcript_consequence):
    # Keep only transcript consequences above 1 - (23/41) ~= 0.439, which is
    # the gentropy-style VEP ranking cutoff the pipeline now uses for coding
    # target assignment. Gene assignment is then driven by the max remaining
    # consequenceScore, restricted to protein-coding transcripts.
    return (
        transcript_consequence["targetId"].isNotNull()
        & transcript_consequence["consequenceScore"].isNotNull()
        & (transcript_consequence["consequenceScore"] > F.lit(prep_coding_TRANSCRIPT_CONSEQUENCE_SCORE_MIN))
        & (F.lower(transcript_consequence["biotype"]) == F.lit(prep_coding_TRANSCRIPT_CONSEQUENCE_PROTEIN_BIOTYPE))
    )


def _array_struct_has_fields(array_dtype, required_fields) -> bool:
    return (
        isinstance(array_dtype, T.ArrayType)
        and isinstance(array_dtype.elementType, T.StructType)
        and all(field in array_dtype.elementType.fieldNames() for field in required_fields)
    )


def prep_coding_parse_args():
    """
    Parse CLI arguments.

    Args:
        argv[1]: genes CSV (1 gene per row, Ensembl IDs)
        argv[2]: credible_set parquet path (dir/file/glob)
        argv[3]: variant parquet path (dir/file/glob)
        argv[4]: existing GWAS variants parquet path to exclude (dir/file/glob)
        argv[5]: output parquet dataset path (optional)
    """
    if len(sys.argv) < 5:
        raise ValueError(
            "Usage: python get_coding_GWAS_NoColoc_spark.py <genes.csv> <credible_set_parquet> "
            "<variant_parquet> <existing_gwas_variants_parquet> [output_parquet]"
        )

    genes_file = sys.argv[1]
    credible_set_path = resolve_parquet_input(sys.argv[2])
    variant_path = resolve_parquet_input(sys.argv[3])
    existing_gwas_path = resolve_parquet_input(sys.argv[4])
    output_path = sys.argv[5] if len(sys.argv) > 5 else prep_coding_DEFAULT_OUTPUT_FILE

    return genes_file, credible_set_path, variant_path, existing_gwas_path, output_path


def prep_coding_load_genes(genes_file):
    """
    Load Ensembl gene IDs from the first CSV column.
    """
    return load_gene_ids_from_csv(genes_file)


def prep_coding_main():
    # 1) Parse CLI inputs and load requested genes.
    genes_file, credible_set_path, variant_path, existing_gwas_path, output_path = prep_coding_parse_args()
    genes = prep_coding_load_genes(genes_file)

    # 2) Build Spark session and load parquet inputs.
    spark = build_local_spark("coding_gwas_nocoloc_spark")
    try:
        credible_set_raw = spark.read.parquet(credible_set_path)
        variant_raw = spark.read.parquet(variant_path)
        existing_gwas_raw = spark.read.parquet(existing_gwas_path)

        # Keep the standalone compatibility entrypoint aligned with the shared
        # Step 1 helper so transcript-consequence target selection and coding
        # consequence filtering are implemented in one place only.
        out_df = _build_coding_output_df(
            spark,
            genes,
            credible_set_raw,
            variant_raw,
            existing_gwas_raw,
        )

        write_parquet_dataset(out_df, output_path)
    finally:
        try:
            spark.stop()
        except Exception:  # pragma: no cover - defensive shutdown path
            pass

# ============================================================================
# Inlined Module: get_QTL_vars_spark.py
# ============================================================================
########################################################################################################
### Script to extract QTL summary statistics for coloc lead variants using Spark parquet sources.    ###
### Script written with GPT-5-Codex on the 19.02.2026                                                ###
########################################################################################################
prep_qtl_DEFAULT_OUTPUT_FILE = "QTL_variants.parquet"
prep_qtl_COLOC_REQUIRED_COLUMNS = [
    "rightStudyLocusId",
    "rightVariantId",
    "rightGeneId",
    "rightBiosampleId",
]

prep_qtl_CREDIBLE_SET_REQUIRED_COLUMNS = [
    "studyLocusId",
    "variantId",
    "beta",
    "pValueMantissa",
    "pValueExponent",
    "studyType",
    "studyId",
]


def prep_qtl_parse_args():
    """
    Parse CLI arguments.

    Args:
        argv[1]: coloc parquet dataset generated by get_coloc_vars_per_gene_spark.py
        argv[2]: credible_set parquet path (dir/file/glob)
        argv[3]: output parquet dataset path (optional)
    """
    if len(sys.argv) < 3:
        raise ValueError(
            "Usage: python get_QTL_vars_spark.py <coloc_variants_parquet> <credible_set_parquet> [output_parquet]"
        )

    coloc_path = resolve_parquet_input(sys.argv[1])
    credible_set_path = resolve_parquet_input(sys.argv[2])
    output_path = sys.argv[3] if len(sys.argv) > 3 else prep_qtl_DEFAULT_OUTPUT_FILE

    return coloc_path, credible_set_path, output_path


def prep_qtl_main():
    # 1) Parse CLI inputs.
    coloc_path, credible_set_path, output_path = prep_qtl_parse_args()

    # 2) Build Spark session and load parquet inputs.
    spark = build_local_spark("qtl_vars_spark")
    try:
        coloc_raw = spark.read.parquet(coloc_path)
        credible_set_raw = spark.read.parquet(credible_set_path)

        # 3) Validate required source schema columns.
        ensure_columns_present(coloc_raw, prep_qtl_COLOC_REQUIRED_COLUMNS, "Colocalising variants")
        ensure_columns_present(credible_set_raw, prep_qtl_CREDIBLE_SET_REQUIRED_COLUMNS, "Credible set")

        credible_schema = credible_set_raw.schema

        # 4) Keep right-side locus/variant keys from colocalisation output.
        coloc_qtl = (
            coloc_raw.select(
                "rightStudyLocusId",
                "rightVariantId",
                "rightGeneId",
                "rightBiosampleId",
            )
            .dropna(subset=["rightStudyLocusId", "rightVariantId"])
            .dropDuplicates(["rightStudyLocusId", "rightVariantId", "rightGeneId", "rightBiosampleId"])
        )

        # 5) Select required credible-set summary statistic columns.
        credible_selected = (
            credible_set_raw
            .select(
                F.col("studyLocusId").cast("string").alias("rightStudyLocusId"),
                F.col("variantId").cast("string").alias("variantId"),
                optional_col("beta", credible_schema, "double").alias("beta"),
                optional_col("standardError", credible_schema, "double").alias("standardError"),
                optional_col("zScore", credible_schema, "double").alias("zScore"),
                optional_col("pValueMantissa", credible_schema, "double").alias("pValueMantissa"),
                optional_col("pValueExponent", credible_schema, "long").alias("pValueExponent"),
                optional_col("studyType", credible_schema, "string").alias("studyType"),
                optional_col("studyId", credible_schema, "string").alias("studyId"),
            )
        )

        # 6) Match right-side colocalisation locus/variant pairs to credible-set rows.
        joined = (
            coloc_qtl.join(
                credible_selected,
                on=(
                    (coloc_qtl["rightStudyLocusId"] == credible_selected["rightStudyLocusId"])
                    & (coloc_qtl["rightVariantId"] == credible_selected["variantId"])
                ),
                how="inner",
            )
            .select(
                credible_selected["rightStudyLocusId"].alias("studyLocusId"),
                credible_selected["variantId"],
                credible_selected["beta"],
                credible_selected["standardError"],
                credible_selected["zScore"],
                credible_selected["pValueMantissa"],
                credible_selected["pValueExponent"],
                credible_selected["studyType"],
                credible_selected["studyId"],
                coloc_qtl["rightGeneId"].cast("string").alias("geneId"),
                coloc_qtl["rightBiosampleId"].cast("string").alias("biosampleId"),
            )
            .dropDuplicates()
        )

        # 7) Fill missing SE from beta/zScore and keep complete effect-size rows only.
        out_df = (
            joined.withColumn(
                "standardError",
                F.coalesce(
                    F.col("standardError"),
                    F.when(
                        F.col("beta").isNotNull()
                        & F.col("zScore").isNotNull()
                        & (F.col("zScore") != F.lit(0.0)),
                        F.col("beta") / F.col("zScore"),
                    ),
                ),
            )
            # Keep only complete effect-size rows for downstream modelling.
            .filter(F.col("beta").isNotNull() & F.col("standardError").isNotNull())
            .drop("zScore")
        )

        # 8) Write deterministic parquet output for downstream modules.
        write_parquet_dataset(out_df, output_path)
    finally:
        try:
            spark.stop()
        except Exception:  # pragma: no cover - defensive shutdown path
            pass

# ============================================================================
# Inlined Module: get_burden_tests_per_gene_spark.py
# ============================================================================
########################################################################################################
### Script to extract burden tests per gene using Spark and local parquet release data.              ###
### It also keeps other sources of burden tests at the moment.                                       ###
### Script written with GPT-5-Codex on the 18.02.2026                                                ###
########################################################################################################
prep_burden_DEFAULT_OUTPUT_FILE = "burden_tests.parquet"
# Burden filtering switches and allowlists (reintroduced from working/vidra).
prep_burden_FILTER_TO_ALLOWED_PROJECTS = True
prep_burden_FILTER_TO_ALLOWED_COHORTS = True
prep_burden_FILTER_TO_ALLOWED_METHODS = True

# Keep aligned with VIDRA_2 method semantics, but on the current release cohort.
prep_burden_ALLOWED_STAT_METHODS = ["pLoF", "ptv", "ptvraredmg", "ptv5pcnt"]
prep_burden_ALLOWED_COHORT_IDS = ["UK Biobank 470k", "UK Biobank 450k"]
prep_burden_ALLOWED_PROJECT_IDS = ["AstraZeneca PheWAS Portal"]


def prep_burden_load_genes(genes_file):
    """
    Read Ensembl IDs from first CSV column (no header expected).
    """
    return load_gene_ids_from_csv(genes_file)


def prep_burden_parse_args():
    """
    Parse CLI arguments.

    Args:
        argv[1]: genes CSV (Ensembl IDs)
        argv[2]: burden evidence parquet path (dir/file/glob), passed from Nextflow config.
        argv[3]: output parquet dataset path (optional).
    """
    if len(sys.argv) < 3:
        raise ValueError(
            "Usage: python get_burden_tests_per_gene_spark.py <genes.csv> <burden_evidence_parquet>"
        )

    genes_file = sys.argv[1]
    evidence_path = sys.argv[2]

    genes = prep_burden_load_genes(genes_file)
    evidence_path = resolve_parquet_input(evidence_path)

    output_file = sys.argv[3] if len(sys.argv) > 3 else prep_burden_DEFAULT_OUTPUT_FILE

    return genes, evidence_path, output_file


def prep_burden_stringified_or_null(col_name, schema):
    """
    Return string representation of a column, JSON-encoding complex types for CSV compatibility.
    """
    if col_name not in schema.fieldNames():
        return F.lit(None).cast("string")

    source_type = schema[col_name].dataType
    if isinstance(source_type, (T.ArrayType, T.MapType, T.StructType)):
        return F.to_json(F.col(col_name))

    return F.col(col_name).cast("string")


def prep_burden_main():
    # 1) Parse CLI inputs (genes list, burden evidence source, output path).
    genes, evidence_path, output_file = prep_burden_parse_args()

    # 2) Build Spark session and load burden evidence parquet input.
    spark = build_local_spark("burden_tests_per_gene_spark")
    try:
        raw = spark.read.parquet(evidence_path)
        schema = raw.schema

        # 3) Validate minimum required schema.
        if "targetId" not in schema.fieldNames():
            raise ValueError("Input burden evidence parquet is missing required column: targetId")

        # 4) Restrict records to requested genes (matched on targetId/Ensembl ID).
        genes_df = spark.createDataFrame([(gene,) for gene in genes], ["targetId"])
        filtered = raw.join(F.broadcast(genes_df), "targetId", "inner")

        # 5) Keep gene-burden datasource rows when datasourceId is present.
        if "datasourceId" in schema.fieldNames():
            filtered = filtered.filter(F.col("datasourceId") == F.lit("gene_burden"))

        # 5b) Reintroduced burden allowlist filters from working/vidra.
        if prep_burden_FILTER_TO_ALLOWED_PROJECTS and "projectId" in schema.fieldNames():
            filtered = filtered.filter(F.col("projectId").isin(prep_burden_ALLOWED_PROJECT_IDS))
        if prep_burden_FILTER_TO_ALLOWED_COHORTS and "cohortId" in schema.fieldNames():
            filtered = filtered.filter(F.col("cohortId").isin(prep_burden_ALLOWED_COHORT_IDS))
        if prep_burden_FILTER_TO_ALLOWED_METHODS and "statisticalMethod" in schema.fieldNames():
            filtered = filtered.filter(
                F.col("statisticalMethod").isin(prep_burden_ALLOWED_STAT_METHODS)
            )

        # 6) Select output-compatible columns, adding nulls for optional missing fields.
        selected = filtered.select(
            optional_col("datasourceId", schema, "string").alias("datasourceId"),
            optional_col("targetId", schema, "string").alias("targetId"),
            prep_burden_stringified_or_null("allelicRequirements", schema).alias("allelicRequirements"),
            optional_col("ancestryId", schema, "string").alias("ancestryId"),
            optional_col("beta", schema, "double").alias("beta"),
            optional_col("betaConfidenceIntervalLower", schema, "double").alias("betaConfidenceIntervalLower"),
            optional_col("betaConfidenceIntervalUpper", schema, "double").alias("betaConfidenceIntervalUpper"),
            optional_col("cohortId", schema, "string").alias("cohortId"),
            optional_col("diseaseFromSource", schema, "string").alias("diseaseFromSource"),
            optional_col("diseaseFromSourceMappedId", schema, "string").alias("diseaseFromSourceMappedId"),
            optional_col("oddsRatio", schema, "double").alias("oddsRatio"),
            optional_col("oddsRatioConfidenceIntervalLower", schema, "double").alias("oddsRatioConfidenceIntervalLower"),
            optional_col("oddsRatioConfidenceIntervalUpper", schema, "double").alias("oddsRatioConfidenceIntervalUpper"),
            optional_col("pValueExponent", schema, "long").alias("pValueExponent"),
            optional_col("pValueMantissa", schema, "double").alias("pValueMantissa"),
            optional_col("projectId", schema, "string").alias("projectId"),
            optional_col("resourceScore", schema, "double").alias("resourceScore"),
            optional_col("statisticalMethod", schema, "string").alias("statisticalMethod"),
            optional_col("targetFromSourceId", schema, "string").alias("targetFromSourceId"),
            optional_col("score", schema, "double").alias("score"),
            optional_col("studyCases", schema, "long").alias("studyCases"),
            optional_col("studyCasesWithQualifyingVariants", schema, "long").alias("studyCasesWithQualifyingVariants"),
            optional_col("studySampleSize", schema, "long").alias("studySampleSize"),
            F.lit(0).cast("int").alias("burdenSourceRank"),
        )

        # 7) Harmonize OR fields from beta-scale inputs and apply CI lower imputation rule.
        burden_df = (
            selected
            .withColumn(
                "oddsRatio",
                F.coalesce(F.col("oddsRatio"), F.exp(F.col("beta")))
            )
            .withColumn(
                "oddsRatioConfidenceIntervalLower",
                F.coalesce(
                    F.col("oddsRatioConfidenceIntervalLower"),
                    F.exp(F.col("betaConfidenceIntervalLower"))
                ),
            )
            .withColumn(
                "oddsRatioConfidenceIntervalUpper",
                F.coalesce(
                    F.col("oddsRatioConfidenceIntervalUpper"),
                    F.exp(F.col("betaConfidenceIntervalUpper"))
                ),
            )
            # If lower CI is missing but OR and upper CI exist, impute lower as OR^2 / upper.
            # This follows the standard symmetry assumption on the log(OR) scale.
            .withColumn(
                "oddsRatioConfidenceIntervalLower",
                F.when(
                    F.col("oddsRatioConfidenceIntervalLower").isNull()
                    & F.col("oddsRatio").isNotNull()
                    & F.col("oddsRatioConfidenceIntervalUpper").isNotNull()
                    & (F.col("oddsRatioConfidenceIntervalUpper") > F.lit(0.0)),  # Keep upper CI positive.
                    (F.col("oddsRatio") * F.col("oddsRatio")) / F.col("oddsRatioConfidenceIntervalUpper"),
                ).otherwise(F.col("oddsRatioConfidenceIntervalLower")),
            )
            .filter(F.col("oddsRatio").isNotNull())
            .drop("beta", "betaConfidenceIntervalLower", "betaConfidenceIntervalUpper")
        )

        # 8) Write deterministic parquet output for downstream modules.
        write_parquet_dataset(burden_df, output_file)
    finally:
        try:
            spark.stop()
        except Exception:  # pragma: no cover - defensive shutdown path
            pass

# ============================================================================
# Inlined Module: get_clinvar_vars_per_gene_spark.py
# ============================================================================
########################################################################################################
### Script to extract ClinVar variants per gene using Spark for scalable Parquet processing.         ###
### Script written with GPT-5.3-Codex on the 13.02.2026                                              ###
########################################################################################################
prep_clinvar_DEFAULT_OUTPUT_FILE = "clinvar_variants.parquet"

### Confidence filters also containing legacy values to be safe
## These are the possible criterias:
# no assertion criteria provided                        | +0.00   Excluded                     
# no classification provided                            | +0.00   Excluded  USED TO BE "no assertion provided"
# no assertion for the individual variant               | +0.00   Excluded  NOT PRESENT ANYMORE
# criteria provided, conflicting classifications        | +0.02   Excluded  USED TO BE "criteria provided, conflicting interpretation"
# criteria provided, single submitter                   | +0.02   Can be discussed if kept or not | New version kept          
# criteria provided, multiple submitters, no conflicts  | +0.05   Kept
# reviewed by expert panel                              | +0.07   Kept              
# practice guideline                                    | +0.10   Kept                          

prep_clinvar_EXCLUDED_CONFIDENCE = [
    "no assertion criteria provided",
    "no classification provided",
    "no assertion provided", # legacy value, not present in current data but included for safety
    "no assertion for the individual variant", # legacy value, not present in current data but included for safety
    "criteria provided, conflicting classifications",
    "criteria provided, conflicting interpretation", # legacy value, not present in current data but included for safety
]


def prep_clinvar_parse_args():
    """
    Parse CLI arguments.

    Args:
        argv[1]: CSV file with one gene ID per row (no header expected).
        argv[2]: Folder containing multiple Parquet files with ClinVar/eva evidence.
        argv[3]: Output parquet dataset path (optional).
    """
    if len(sys.argv) < 3:
        raise ValueError(
            "Usage: python get_clinvar_vars_per_gene_spark.py <genes.csv> <clinvar_evidence_dir> [output_parquet]"
        )

    genes_file = sys.argv[1]
    evidence_dir = sys.argv[2] if len(sys.argv) > 2 else ""

    if not evidence_dir or evidence_dir.lower() in ("none", "null"):
        raise ValueError("Missing evidence parquet folder path. Pass a valid directory as second argument.")
    evidence_dir = resolve_parquet_input(evidence_dir)

    output_file = sys.argv[3] if len(sys.argv) > 3 else prep_clinvar_DEFAULT_OUTPUT_FILE

    return genes_file, evidence_dir, output_file


def prep_clinvar_load_genes(genes_file):
    """
    Load the genes input channel file into a plain Python list of strings.

    The pipeline passes a 1-column CSV with Ensembl gene IDs. We preserve that
    convention and avoid any extra normalization so behavior stays consistent.
    """
    return load_gene_ids_from_csv(genes_file)


def prep_clinvar_add_clinical_significance_column(df, schema):
    """
    Create a scalar `_clinicalSignificance` column compatible with VIDRA_2 output.

    VIDRA_2 used the first clinical significance value from OT evidence rather
    than expanding one evidence row into multiple rows.
    Supported source layouts:
      1) `clinicalSignificances` as array<string>
      2) `clinicalSignificances.list` as array<string> (nested struct form)
      3) fallback scalar `clinicalSignificances` / `list`
      4) missing fields -> null
    """
    if "clinicalSignificances" in schema.fieldNames():
        clinical_dtype = schema["clinicalSignificances"].dataType
        if isinstance(clinical_dtype, T.ArrayType):
            return df.withColumn("_clinicalSignificance", F.col("clinicalSignificances").getItem(0))
        if isinstance(clinical_dtype, T.StructType) and "list" in clinical_dtype.fieldNames():
            return df.withColumn("_clinicalSignificance", F.col("clinicalSignificances.list").getItem(0))
        return df.withColumn("_clinicalSignificance", F.col("clinicalSignificances").cast("string"))

    if "list" in schema.fieldNames():
        list_dtype = schema["list"].dataType
        if isinstance(list_dtype, T.ArrayType):
            return df.withColumn("_clinicalSignificance", F.col("list").getItem(0))
        return df.withColumn("_clinicalSignificance", F.col("list").cast("string"))

    return df.withColumn("_clinicalSignificance", F.lit(None).cast("string"))


def prep_clinvar_main():
    # 1) Parse CLI inputs (genes list, evidence folder, output path).
    genes_file, evidence_dir, output_file = prep_clinvar_parse_args()
    genes = prep_clinvar_load_genes(genes_file)

    # 2) Build Spark session and load ClinVar evidence parquet files.
    spark = build_local_spark("clinvar_vars_per_gene_spark")
    try:
        raw = spark.read.parquet(evidence_dir)
        schema = raw.schema

        # 3) Restrict input to requested genes via a broadcast semi-join.
        genes_df = spark.createDataFrame([(gene.upper(),) for gene in genes], ["targetId_upper"])
        gene_filtered = (
            raw.withColumn("targetId_upper", F.upper(F.col("targetId")))
            .join(F.broadcast(genes_df), "targetId_upper", "inner")
            .drop("targetId_upper")
        )

        # 4) Apply the same filters as the legacy BigQuery implementation.
        filtered = (
            gene_filtered
            .filter(F.col("datasourceId") == F.lit("eva"))
            .filter(~F.col("confidence").isin(prep_clinvar_EXCLUDED_CONFIDENCE))
            .filter(F.col("variantId").isNotNull())
        )

        # 5) Normalize clinical significance into one scalar field per row.
        with_clin_sig = prep_clinvar_add_clinical_significance_column(filtered, schema)

        # 6) Keep output schema/order identical to the non-Spark implementation.
        out = with_clin_sig.select(
            F.col("_clinicalSignificance").cast("string").alias("clinicalSignificances"),
            optional_col("confidence", schema, "string").alias("confidence"),
            optional_col("diseaseFromSourceId", schema, "string").alias("diseaseFromSourceId"),
            optional_col("diseaseFromSourceMappedId", schema, "string").alias("diseaseFromSourceMappedId"),
            optional_col("studyId", schema, "string").alias("studyId"),
            optional_col("targetFromSourceId", schema, "string").alias("targetFromSourceId"),
            optional_col("targetId", schema, "string").alias("targetId"),
            optional_col("variantFunctionalConsequenceId", schema, "string").alias("variantFunctionalConsequenceId"),
            optional_col("variantHgvsId", schema, "string").alias("variantHgvsId"),
            optional_col("variantId", schema, "string").alias("variantId"),
            optional_col("variantRsId", schema, "string").alias("variantRsId"),
            optional_col("diseaseId", schema, "string").alias("diseaseId"),
            optional_col("score", schema, "double").alias("score"),
            optional_col("datasourceId", schema, "string").alias("datasourceId"),
        )

        # 7) Write deterministic parquet output for downstream modules.
        write_parquet_dataset(out, output_file)
    finally:
        spark.stop()


# ============================================================================
# Analysis-ready assembly (bucket-compatible output)
# ============================================================================
def _join_bucket_path(base_uri: str, name: str) -> str:
    return f"{base_uri.rstrip('/')}/{name.lstrip('/')}"


def _resolve_bucket_base(args) -> str:
    if args.bucket_uri:
        return args.bucket_uri.rstrip("/")
    if args.bucket_name:
        return f"gs://{args.bucket_name}"
    return ""


def _study_disease_mapping(study_raw):
    # Disease mapping is intentionally two-stage:
    # 1) prefer mapped trait/EFO IDs from `traitFromSourceMappedIds`
    # 2) fall back to `diseaseIds` only for studyIds that have no mapped trait
    #
    # This matches the "best available mapped disease" policy we settled on for
    # both the release study dataset and the lead-variant-effect replacement.
    columns = set(study_raw.columns)
    study_id_col = "studyId" if "studyId" in columns else "study_id" if "study_id" in columns else None
    if study_id_col is None:
        raise ValueError("Study dataset must contain 'studyId' or 'study_id'.")
    if "traitFromSourceMappedIds" not in columns and "diseaseIds" not in columns:
        raise ValueError(
            "Study dataset must contain 'traitFromSourceMappedIds' or 'diseaseIds'. "
            f"Available columns sample: {sorted(list(columns))[:20]}"
        )

    def _mapping_from_field(field_name: str):
        if field_name not in columns:
            return None
        dtype = study_raw.schema[field_name].dataType
        if isinstance(dtype, T.ArrayType):
            return (
                study_raw.select(
                    F.col(study_id_col).cast("string").alias("studyId"),
                    F.explode_outer(F.col(field_name)).alias("_mappedDisease"),
                )
                .withColumn("mappedDisease", F.col("_mappedDisease").cast("string"))
                .drop("_mappedDisease")
                .dropna(subset=["studyId", "mappedDisease"])
                .dropDuplicates(["studyId", "mappedDisease"])
            )
        return (
            study_raw.select(
                F.col(study_id_col).cast("string").alias("studyId"),
                F.col(field_name).cast("string").alias("mappedDisease"),
            )
            .dropna(subset=["studyId", "mappedDisease"])
            .dropDuplicates(["studyId", "mappedDisease"])
        )

    trait_mapping = _mapping_from_field("traitFromSourceMappedIds")
    disease_mapping = _mapping_from_field("diseaseIds")

    if trait_mapping is None:
        return disease_mapping
    if disease_mapping is None:
        return trait_mapping

    trait_studies = trait_mapping.select("studyId").distinct()
    disease_fallback = disease_mapping.join(trait_studies, "studyId", "left_anti")
    return trait_mapping.unionByName(disease_fallback, allowMissingColumns=True).dropDuplicates(["studyId", "mappedDisease"])


def _clinvar_significance_expr():
    denom = float(len(CLINICAL_SIG_ORDER) - 1)
    map_items = []
    for idx, category in enumerate(CLINICAL_SIG_ORDER):
        map_items.extend([F.lit(category), F.lit(float(idx) / denom)])
    lookup = F.create_map(*map_items)

    normalised = F.lower(
        F.trim(
            F.regexp_replace(
                F.coalesce(F.col("clinicalSignificances"), F.lit("not provided")),
                "_",
                " ",
            )
        )
    )
    return lookup[normalised]


def _pvalue_from_mantissa_exponent(mantissa_col: str, exponent_col: str):
    return F.when(
        F.col(mantissa_col).isNotNull() & F.col(exponent_col).isNotNull(),
        F.col(mantissa_col).cast("double")
        * F.exp(F.log(F.lit(10.0)) * F.col(exponent_col).cast("double")),
    )


def _load_az_gene_map(spark: SparkSession, path: str):
    raw = (
        spark.read
        .option("header", "false")
        .option("inferSchema", "false")
        .csv(path)
    )
    if "_c0" not in raw.columns or "_c1" not in raw.columns:
        raise ValueError(f"AZ gene map must be a 2-column CSV without header: {path}")

    return (
        raw.select(
            F.col("_c0").cast("string").alias("ensembl_id"),
            F.col("_c1").cast("string").alias("gene_symbol"),
        )
        .dropna(subset=["ensembl_id", "gene_symbol"])
        .dropDuplicates(["ensembl_id", "gene_symbol"])
    )


def _load_az_efo_map(spark: SparkSession, path: str):
    raw = (
        spark.read
        .option("header", "true")
        .option("sep", "\t")
        .csv(path)
    )

    return (
        raw.filter(F.col("STUDY").contains("AstraZeneca"))
        .withColumn("diseaseId", F.element_at(F.split(F.col("SEMANTIC_TAG"), "/"), -1))
        .select(
            F.col("PROPERTY_VALUE").cast("string").alias("diseaseFromSource"),
            F.col("diseaseId").cast("string").alias("diseaseId"),
        )
        .dropna(subset=["diseaseFromSource", "diseaseId"])
        .dropDuplicates(["diseaseFromSource", "diseaseId"])
    )


def _load_az_burden_efo_map(spark: SparkSession, path: str):
    raw = (
        spark.read
        .option("header", "true")
        .option("sep", "\t")
        .csv(path)
    )
    ensure_columns_present(raw, ["PROPERTY_VALUE", "SEMANTIC_TAG"], "AZ phenotype-to-EFO map")

    return (
        raw
        .withColumn("diseaseFromSource", F.trim(F.col("PROPERTY_VALUE").cast("string")))
        .withColumn("diseaseId", F.element_at(F.split(F.col("SEMANTIC_TAG").cast("string"), "/"), -1))
        .dropna(subset=["diseaseFromSource", "diseaseId"])
        # Keep every phenotype->EFO mapping. A phenotype can legitimately map to
        # multiple EFO terms, and the final burden join later resolves duplicate
        # (gene, disease) rows by p-value instead of choosing an arbitrary EFO here.
        .dropDuplicates(["diseaseFromSource", "diseaseId"])
    )


def _load_target_symbol_map(spark: SparkSession, path: str):
    raw = spark.read.parquet(resolve_parquet_input(path))
    ensure_columns_present(raw, ["id", "approvedSymbol"], "target data")

    return (
        raw.select(
            F.col("id").cast("string").alias("targetId"),
            F.upper(F.trim(F.col("approvedSymbol").cast("string"))).alias("gene_symbol_norm"),
        )
        .dropna(subset=["targetId", "gene_symbol_norm"])
        .dropDuplicates(["targetId", "gene_symbol_norm"])
    )


def _normalise_lead_variant_effect_df(raw):
    # The lead-variant-effect path replaces release credible_set/study/variant
    # for GWAS, QTL and coding evidence, but it must still emit the same Step 1
    # internal contract: one row with usable beta/SE per studyLocusId/variant.
    #
    # Decisions encoded here:
    # - use signed rescaled beta when direction + absEstimatedBeta are present
    # - otherwise fall back to originalBeta
    # - prefer estimatedSE, then reconstruct from absZScore, then fall back to
    #   originalStandardError
    # - keep only significant rows (p <= 5e-8)
    # - drop rows with zero allele frequency
    # - drop rows whose final signed beta is outside [-3, 3]
    #
    # The beta bound is applied after all fallback logic so both the rescaled
    # path and originalBeta fallback obey the same sanity filter.
    ensure_columns_present(
        raw,
        ["studyLocusId", "variantId", "studyId", "geneId", "biosampleId"],
        "Lead variant effect",
    )
    ensure_nested_fields_present(
        raw,
        [
            "studyStatistics.studyType",
            "variantStatistics.pValueMantissa",
            "variantStatistics.pValueExponent",
            "leadVariantConsequence.mostSevereConsequence.transcriptConsequence.targetId",
            "leadVariantConsequence.mostSevereConsequence.transcriptConsequence.biotype",
            "leadVariantConsequence.mostSevereConsequence.transcriptConsequence.consequenceScore",
        ],
        "Lead variant effect",
    )

    schema = raw.schema
    lead_df = (
        raw.select(
            F.col("studyLocusId").cast("string").alias("studyLocusId"),
            F.col("variantId").cast("string").alias("variantId"),
            F.col("studyId").cast("string").alias("studyId"),
            optional_col("geneId", schema, "string").alias("geneId"),
            optional_col("biosampleId", schema, "string").alias("biosampleId"),
            optional_nested_col("studyStatistics.studyType", schema, "string").alias("studyType"),
            optional_nested_col("variantStatistics.pValueMantissa", schema, "double").alias("pValueMantissa"),
            optional_nested_col("variantStatistics.pValueExponent", schema, "long").alias("pValueExponent"),
            optional_nested_col("majorLdPopulationAf.alleleFrequency", schema, "double").alias("alleleFrequency"),
            optional_col("originalBeta", schema, "double").alias("originalBeta"),
            optional_col("originalStandardError", schema, "double").alias("originalStandardError"),
            optional_nested_col("rescaledStatistics.directionOfEffect", schema, "double").alias("directionOfEffect"),
            optional_nested_col("rescaledStatistics.absEstimatedBeta", schema, "double").alias("absEstimatedBeta"),
            optional_nested_col("rescaledStatistics.estimatedSE", schema, "double").alias("estimatedSE"),
            optional_nested_col("rescaledStatistics.absZScore", schema, "double").alias("absZScore"),
            optional_nested_col(
                "leadVariantConsequence.mostSevereConsequence.transcriptConsequence.variantFunctionalConsequenceIds",
                schema,
                T.ArrayType(T.StringType()),
            ).alias("mostSevereConsequenceIds"),
            optional_nested_col(
                "leadVariantConsequence.mostSevereConsequence.transcriptConsequence.targetId",
                schema,
                "string",
            ).alias("leadTranscriptTargetId"),
            optional_nested_col(
                "leadVariantConsequence.mostSevereConsequence.transcriptConsequence.biotype",
                schema,
                "string",
            ).alias("leadTranscriptBiotype"),
            optional_nested_col(
                "leadVariantConsequence.mostSevereConsequence.transcriptConsequence.consequenceScore",
                schema,
                "double",
            ).alias("leadTranscriptConsequenceScore"),
            optional_nested_col(
                "leadVariantConsequence.mostSevereConsequence.type",
                schema,
                "string",
            ).alias("variantDescription"),
        )
        .withColumn(
            "_signedRescaledBeta",
            F.when(
                F.col("directionOfEffect").isNotNull()
                & (F.col("directionOfEffect") != F.lit(0.0))
                & F.col("absEstimatedBeta").isNotNull(),
                F.when(F.col("directionOfEffect") > F.lit(0.0), F.lit(1.0)).otherwise(F.lit(-1.0))
                * F.col("absEstimatedBeta"),
            ),
        )
        .withColumn("beta", F.coalesce(F.col("_signedRescaledBeta"), F.col("originalBeta")))
        .withColumn(
            "standardError",
            F.coalesce(
                F.when(F.col("estimatedSE").isNotNull() & (F.col("estimatedSE") > F.lit(0.0)), F.col("estimatedSE")),
                F.when(
                    F.col("beta").isNotNull()
                    & F.col("absZScore").isNotNull()
                    & (F.col("absZScore") > F.lit(0.0)),
                    F.abs(F.col("beta")) / F.col("absZScore"),
                ),
                F.col("originalStandardError"),
            ),
        )
        .withColumn("pValue", _pvalue_from_mantissa_exponent("pValueMantissa", "pValueExponent"))
        .filter(F.col("studyLocusId").isNotNull() & F.col("variantId").isNotNull() & F.col("studyId").isNotNull())
        .filter(F.col("pValue").isNotNull() & (F.col("pValue") <= F.lit(LEAD_VARIANT_EFFECT_PVALUE_THRESHOLD)))
        .filter(F.abs(F.col("beta")) <= F.lit(LEAD_VARIANT_EFFECT_BETA_ABS_MAX))
        .filter(F.col("alleleFrequency").isNull() | (F.col("alleleFrequency") != F.lit(0.0)))
        .filter(F.col("beta").isNotNull() & F.col("standardError").isNotNull())
        .select(
            "studyLocusId",
            "variantId",
            "studyId",
            "geneId",
            "biosampleId",
            "studyType",
            "beta",
            "standardError",
            "pValueMantissa",
            "pValueExponent",
            "mostSevereConsequenceIds",
            "leadTranscriptTargetId",
            "leadTranscriptBiotype",
            "leadTranscriptConsequenceScore",
            "variantDescription",
        )
        .dropDuplicates()
    )
    return lead_df


def _build_coloc_output_df_from_lead(spark: SparkSession, genes, coloc_threshold: float, coloc_raw, lead_variant_effect_df):
    ensure_columns_present(coloc_raw, prep_coloc_COLOC_COLUMNS, "Colocalisation")
    ensure_columns_present(
        lead_variant_effect_df,
        ["studyLocusId", "studyId", "variantId", "studyType", "geneId", "biosampleId"],
        "Lead variant effect",
    )

    coloc_df = coloc_raw.select(*prep_coloc_COLOC_COLUMNS)
    coloc_filtered_df = coloc_df.filter(F.col("h4") > F.lit(coloc_threshold))

    if prep_coloc_ALLOWED_RIGHT_STUDY_TYPES:
        allowed_right_types = [x.lower() for x in prep_coloc_ALLOWED_RIGHT_STUDY_TYPES]
        coloc_filtered_df = coloc_filtered_df.filter(
            F.lower(F.col("rightStudyType")).isin(allowed_right_types)
        )

    genes_df = spark.createDataFrame([(gene.upper(),) for gene in genes], ["rightGeneId"])

    # Raw coloc parquet only exposes locus ids on each side. Variant ids come
    # from the summary-stat source joined behind those loci. The lead-backed
    # path therefore mirrors the release path: first join right loci to lead
    # rows, then join left loci, instead of assuming rightVariantId already
    # exists in raw coloc.
    right_loci_df = coloc_filtered_df.select("rightStudyLocusId").dropDuplicates(["rightStudyLocusId"])
    right_lead_df = (
        lead_variant_effect_df.select(
            F.col("studyLocusId").alias("rightStudyLocusId"),
            F.col("variantId").alias("rightVariantId"),
            F.col("studyId").alias("rightStudyId"),
            F.upper(F.col("geneId")).alias("rightGeneId"),
            F.col("biosampleId").alias("rightBiosampleId"),
        )
        .dropna(subset=["rightStudyLocusId", "rightVariantId", "rightStudyId", "rightGeneId"])
        .join(F.broadcast(genes_df), "rightGeneId", "inner")
        .join(right_loci_df, ["rightStudyLocusId"], "inner")
        .dropDuplicates()
    )

    coloc_with_right_df = coloc_filtered_df.join(
        right_lead_df,
        ["rightStudyLocusId"],
        "inner",
    )

    left_loci_df = coloc_with_right_df.select("leftStudyLocusId").dropDuplicates(["leftStudyLocusId"])
    left_lead_df = lead_variant_effect_df.select(
        F.col("studyLocusId").alias("leftStudyLocusId"),
        F.col("variantId").alias("leftVariantId"),
        F.col("studyType").alias("leftStudyType"),
        F.col("studyId").alias("leftStudyId"),
    )
    if prep_coloc_FILTER_LEFT_STUDY_TYPE_TO_GWAS:
        left_lead_df = left_lead_df.filter(F.lower(F.col("leftStudyType")) == F.lit("gwas"))

    left_lead_df = (
        left_lead_df
        .dropna(subset=["leftStudyLocusId", "leftVariantId", "leftStudyId"])
        .dropDuplicates(["leftStudyLocusId"])
        .join(left_loci_df, ["leftStudyLocusId"], "inner")
    )

    joined_df = coloc_with_right_df.join(
        left_lead_df,
        ["leftStudyLocusId"],
        "inner",
    )

    if prep_coloc_FILTER_QTL_TO_BLOOD_BIOSAMPLE_IDS:
        joined_df = joined_df.filter(F.col("rightBiosampleId").isin(prep_coloc_BLOOD_BIOSAMPLE_IDS))

    return (
        joined_df.select(
            "leftStudyLocusId",
            "rightStudyLocusId",
            "rightStudyType",
            "chromosome",
            "colocalisationMethod",
            "numberColocalisingVariants",
            "h3",
            "h4",
            "clpp",
            "betaRatioSignAverage",
            "leftStudyType",
            "leftStudyId",
            "leftVariantId",
            "rightStudyId",
            "rightBiosampleId",
            "rightVariantId",
            "rightGeneId",
        )
        .dropDuplicates()
    )


def _build_gwas_output_df_from_lead(coloc_raw, lead_variant_effect_df):
    ensure_columns_present(coloc_raw, prep_gwas_COLOC_REQUIRED_COLUMNS, "Colocalising variants")
    ensure_columns_present(
        lead_variant_effect_df,
        ["studyLocusId", "variantId", "beta", "standardError", "pValueMantissa", "pValueExponent", "studyType", "studyId"],
        "Lead variant effect",
    )

    coloc_gwas = (
        coloc_raw.select("leftStudyLocusId", "leftVariantId", "leftStudyType")
        .filter(F.lower(F.col("leftStudyType")) == F.lit("gwas"))
        .dropna(subset=["leftStudyLocusId", "leftVariantId"])
        .dropDuplicates(["leftStudyLocusId", "leftVariantId"])
    )

    lead_selected = lead_variant_effect_df.select(
        F.col("studyLocusId").cast("string").alias("leftStudyLocusId"),
        F.col("variantId").cast("string").alias("variantId"),
        F.col("beta").cast("double").alias("beta"),
        F.col("standardError").cast("double").alias("standardError"),
        F.col("pValueMantissa").cast("double").alias("pValueMantissa"),
        F.col("pValueExponent").cast("long").alias("pValueExponent"),
        F.col("studyType").cast("string").alias("studyType"),
        F.col("studyId").cast("string").alias("studyId"),
    )

    return (
        coloc_gwas.join(
            lead_selected,
            on=(
                (coloc_gwas["leftStudyLocusId"] == lead_selected["leftStudyLocusId"])
                & (coloc_gwas["leftVariantId"] == lead_selected["variantId"])
            ),
            how="inner",
        )
        .select(
            lead_selected["leftStudyLocusId"].alias("studyLocusId"),
            lead_selected["variantId"],
            lead_selected["beta"],
            lead_selected["standardError"],
            lead_selected["pValueMantissa"],
            lead_selected["pValueExponent"],
            lead_selected["studyType"],
            lead_selected["studyId"],
        )
        .dropDuplicates()
    )


def _build_coding_output_df_from_lead(spark: SparkSession, genes, lead_variant_effect_df, existing_gwas_raw):
    ensure_columns_present(existing_gwas_raw, prep_coding_EXISTING_GWAS_REQUIRED_COLUMNS, "Existing GWAS variants")
    ensure_columns_present(
        lead_variant_effect_df,
        [
            "studyLocusId",
            "variantId",
            "beta",
            "standardError",
            "pValueMantissa",
            "pValueExponent",
            "studyType",
            "studyId",
            "leadTranscriptTargetId",
            "leadTranscriptBiotype",
            "leadTranscriptConsequenceScore",
            "variantDescription",
        ],
        "Lead variant effect",
    )

    existing_gwas_variant_ids_df = (
        existing_gwas_raw
        .select(F.col("variantId").cast("string").alias("variantId"))
        .dropna(subset=["variantId"])
        .dropDuplicates()
    )

    # Coding GWAS is "GWAS only" evidence that is not already represented by
    # the coloc/common branch. We therefore anti-join the GWAS variants already
    # selected into `gwas_df` before doing consequence-based coding selection.
    lead_gwas_df = (
        lead_variant_effect_df
        .filter(F.lower(F.col("studyType")) == F.lit("gwas"))
        .dropna(subset=["studyLocusId", "variantId"])
        .join(existing_gwas_variant_ids_df, "variantId", "left_anti")
        .select(
            "studyLocusId",
            "variantId",
            "beta",
            "standardError",
            "pValueMantissa",
            "pValueExponent",
            "studyType",
            "studyId",
            "mostSevereConsequenceIds",
            "leadTranscriptTargetId",
            "leadTranscriptBiotype",
            "leadTranscriptConsequenceScore",
            "variantDescription",
        )
    )

    genes_df = spark.createDataFrame([(gene.upper(),) for gene in genes], ["geneId"])

    normalised_consequence_ids = F.expr(
        "transform(coalesce(mostSevereConsequenceIds, cast(array() as array<string>)), x -> upper(regexp_replace(x, ':', '_')))"
    )
    consequence_match = _coding_consequence_match_column()

    variant_gene_df = (
        lead_gwas_df.select(
            "variantId",
            "mostSevereConsequenceIds",
            "leadTranscriptTargetId",
            "leadTranscriptBiotype",
            "leadTranscriptConsequenceScore",
            "variantDescription",
        )
        .withColumn("_normalisedConsequenceIds", normalised_consequence_ids)
        .withColumn("mostSevereConsequenceId", consequence_match)
        .filter(F.col("mostSevereConsequenceId").isNotNull())
        .withColumn(
            "_leadTranscriptConsequence",
            F.struct(
                F.col("leadTranscriptTargetId").alias("targetId"),
                F.col("leadTranscriptBiotype").alias("biotype"),
                F.col("leadTranscriptConsequenceScore").alias("consequenceScore"),
            ),
        )
        .filter(_coding_transcript_consequence_is_eligible(F.col("_leadTranscriptConsequence")))
        .withColumn("geneId", F.upper(F.col("leadTranscriptTargetId").cast("string")))
        .drop("_leadTranscriptConsequence")
        .dropna(subset=["variantId", "geneId", "leadTranscriptConsequenceScore"])
        .join(F.broadcast(genes_df), "geneId", "inner")
        .select("variantId", "geneId", "mostSevereConsequenceId", "variantDescription")
        .dropDuplicates()
    )

    return (
        lead_gwas_df.alias("lead")
        .join(variant_gene_df.alias("coding"), "variantId", "inner")
        .select(
            F.col("lead.studyLocusId").alias("studyLocusId"),
            F.col("variantId"),
            F.col("lead.beta").alias("beta"),
            F.col("lead.standardError").alias("standardError"),
            F.col("lead.pValueMantissa").alias("pValueMantissa"),
            F.col("lead.pValueExponent").alias("pValueExponent"),
            F.col("lead.studyType").alias("studyType"),
            F.col("lead.studyId").alias("studyId"),
            F.col("coding.geneId").alias("geneId"),
            F.col("coding.mostSevereConsequenceId").alias("mostSevereConsequenceId"),
            F.coalesce(
                F.col("coding.variantDescription"),
                F.col("lead.variantDescription"),
            ).alias("variantDescription"),
        )
        .dropDuplicates()
    )


def _build_qtl_output_df_from_lead(coloc_raw, lead_variant_effect_df):
    ensure_columns_present(coloc_raw, prep_qtl_COLOC_REQUIRED_COLUMNS, "Colocalising variants")
    ensure_columns_present(
        lead_variant_effect_df,
        ["studyLocusId", "variantId", "beta", "standardError", "pValueMantissa", "pValueExponent", "studyType", "studyId", "geneId", "biosampleId"],
        "Lead variant effect",
    )

    coloc_qtl = (
        coloc_raw.select(
            "rightStudyLocusId",
            "rightVariantId",
            "rightGeneId",
            "rightBiosampleId",
        )
        .dropna(subset=["rightStudyLocusId", "rightVariantId"])
        .dropDuplicates(["rightStudyLocusId", "rightVariantId", "rightGeneId", "rightBiosampleId"])
    )

    lead_selected = lead_variant_effect_df.select(
        F.col("studyLocusId").cast("string").alias("rightStudyLocusId"),
        F.col("variantId").cast("string").alias("variantId"),
        F.col("beta").cast("double").alias("beta"),
        F.col("standardError").cast("double").alias("standardError"),
        F.col("pValueMantissa").cast("double").alias("pValueMantissa"),
        F.col("pValueExponent").cast("long").alias("pValueExponent"),
        F.col("studyType").cast("string").alias("studyType"),
        F.col("studyId").cast("string").alias("studyId"),
        F.col("geneId").cast("string").alias("leadGeneId"),
        F.col("biosampleId").cast("string").alias("leadBiosampleId"),
    )

    return (
        coloc_qtl.join(
            lead_selected,
            on=(
                (coloc_qtl["rightStudyLocusId"] == lead_selected["rightStudyLocusId"])
                & (coloc_qtl["rightVariantId"] == lead_selected["variantId"])
            ),
            how="inner",
        )
        .select(
            lead_selected["rightStudyLocusId"].alias("studyLocusId"),
            lead_selected["variantId"],
            lead_selected["beta"],
            lead_selected["standardError"],
            lead_selected["pValueMantissa"],
            lead_selected["pValueExponent"],
            lead_selected["studyType"],
            lead_selected["studyId"],
            F.coalesce(coloc_qtl["rightGeneId"].cast("string"), lead_selected["leadGeneId"]).alias("geneId"),
            F.coalesce(coloc_qtl["rightBiosampleId"].cast("string"), lead_selected["leadBiosampleId"]).alias("biosampleId"),
        )
        .dropDuplicates()
    )


def _deduplicate_az_ready(df):
    # VIDRA_2 effectively kept AZ duplicates late in Step 3. Here we make the
    # rule explicit in Step 1: keep one row per source-aware analysis key using
    # the lowest p-value, then largest odds ratio, then stable string tie-breaks.
    # This preserves source identity while choosing the most informative row.
    az_window = Window.partitionBy("variant", "as_gene", "as_disease", "GsourceLab", "GqtlLab").orderBy(
        F.col("pValue").asc_nulls_last(),
        F.col("oddsRatio").desc_nulls_last(),
        F.col("sourceVariant").asc_nulls_last(),
        F.col("diseaseFromSource").asc_nulls_last(),
    )
    return (
        df
        .withColumn("rn", F.row_number().over(az_window))
        .filter(F.col("rn") == 1)
        .drop("rn")
        .dropDuplicates()
    )


def _build_az_ready_dataset(spark: SparkSession, args) -> T.DataFrame:
    if not path_exists(args.az_variants_file, getattr(args, "gcp_project", "")):
        raise FileNotFoundError(f"AZ variants file not found: {args.az_variants_file}")
    if not path_exists(args.az_mapping_file, getattr(args, "gcp_project", "")):
        raise FileNotFoundError(f"AZ mapping file not found: {args.az_mapping_file}")
    if not path_exists(args.az_gene_map_file, getattr(args, "gcp_project", "")):
        raise FileNotFoundError(f"AZ gene map file not found: {args.az_gene_map_file}")

    genes = load_gene_ids_from_csv(args.genes, getattr(args, "gcp_project", ""))
    genes_df = spark.createDataFrame([(gene,) for gene in genes], ["ensembl_id"])

    gene_map = _load_az_gene_map(spark, args.az_gene_map_file)
    efo_map = _load_az_efo_map(spark, args.az_mapping_file)
    # AZ rare-variant matching always uses the full release OT variant universe.
    # This stays true even when lead_variant_effect is enabled for GWAS/QTL/coding,
    # because AZ evidence should not be narrowed to only lead variants.
    variant_index = (
        spark.read.parquet(resolve_parquet_input(args.variant_data_dir))
        .select(F.col("variantId").cast("string").alias("ot_variant"))
        .dropna(subset=["ot_variant"])
        .dropDuplicates(["ot_variant"])
    )

    az_raw = (
        spark.read
        .option("header", "true")
        .schema(AZ_SCHEMA)
        .csv(args.az_variants_file)
    )

    az_filtered = (
        az_raw
        .filter(F.col("Model") == F.lit("allelic"))
        .filter(F.col("p-value") < F.lit(AZ_PVAL_THRESHOLD))
        .withColumn("Gene", F.regexp_replace(F.col("Gene"), "'", ""))
        .withColumn(
            "variant",
            F.regexp_replace(F.regexp_replace(F.trim(F.col("Variant")), " ", "_"), "-", "_"),
        )
        .filter(F.col("variant").isNotNull())
        .filter(F.col("Odds ratio").isNotNull() & (F.col("Odds ratio") > F.lit(0.0)))
        .filter(F.col("Odds ratio LCI").isNotNull() & (F.col("Odds ratio LCI") > F.lit(0.0)))
        .filter(F.col("Odds ratio UCI").isNotNull() & (F.col("Odds ratio UCI") > F.lit(0.0)))
    )

    az_mapped = (
        az_filtered
        .join(gene_map, az_filtered["Gene"] == gene_map["gene_symbol"], "inner")
        .drop("gene_symbol")
        .join(F.broadcast(genes_df), "ensembl_id", "inner")
        .join(efo_map, F.col("Phenotype") == efo_map["diseaseFromSource"], "left")
        .drop("diseaseFromSource")
        .filter(F.col("diseaseId").isNotNull())
    )

    az_with_match = (
        az_mapped
        .join(variant_index, az_mapped["variant"] == variant_index["ot_variant"], "left")
        .withColumn("is_variant_matched", F.col("ot_variant").isNotNull())
    )

    summary = az_with_match.agg(
        F.count("*").alias("rows"),
        F.sum(F.when(F.col("is_variant_matched"), F.lit(1)).otherwise(F.lit(0))).alias("matched"),
        F.sum(F.when(~F.col("is_variant_matched"), F.lit(1)).otherwise(F.lit(0))).alias("unmatched"),
    ).collect()[0]
    unmatched_sample = [
        row["variant"]
        for row in az_with_match
        .filter(~F.col("is_variant_matched"))
        .select("variant")
        .dropDuplicates()
        .limit(5)
        .collect()
    ]
    log(
        "AZ variant normalization summary against release variant universe: "
        f"rows={summary['rows']}, matched={summary['matched']}, "
        f"unmatched={summary['unmatched']}, sample_unmatched={unmatched_sample}"
    )

    az_ready = (
        az_with_match
        .withColumn("yc", F.log(F.col("Odds ratio")))
        .withColumn(
            "ycse",
            (F.log(F.col("Odds ratio UCI")) - F.log(F.col("Odds ratio LCI"))) / F.lit(3.92),
        )
        .withColumn(
            "ycse",
            F.when(F.col("ycse").isNull() | (F.col("ycse") <= F.lit(0.0)), F.lit(0.14)).otherwise(F.col("ycse")),
        )
        .select(
            F.col("variant").cast("string").alias("variant"),
            F.col("ensembl_id").cast("string").alias("as_gene"),
            F.col("diseaseId").cast("string").alias("as_disease"),
            F.lit(1).cast("int").alias("GsourceLab"),
            F.lit(2).cast("int").alias("GqtlLab"),
            F.col("yc").cast("double").alias("yc"),
            F.col("ycse").cast("double").alias("ycse"),
            F.lit(0.0).cast("double").alias("xc"),
            F.lit(0.1).cast("double").alias("xcse"),
            F.lit(0.0).cast("double").alias("as_clinicalSignificance"),
            F.col("Variant").cast("string").alias("sourceVariant"),
            F.col("Gene").cast("string").alias("geneSymbol"),
            F.col("Phenotype").cast("string").alias("diseaseFromSource"),
            F.col("ot_variant").cast("string").alias("matchedOtVariantId"),
            F.col("is_variant_matched").cast("boolean").alias("is_variant_matched"),
            F.col("p-value").cast("double").alias("pValue"),
            F.col("Odds ratio").cast("double").alias("oddsRatio"),
            F.col("Odds ratio LCI").cast("double").alias("oddsRatioConfidenceIntervalLower"),
            F.col("Odds ratio UCI").cast("double").alias("oddsRatioConfidenceIntervalUpper"),
        )
        .dropna(subset=["variant", "as_gene", "as_disease"])
    )
    return _deduplicate_az_ready(az_ready)




def _maybe_write_section_output(section_name: str, df, output_path: str) -> None:
    if not output_path:
        log(f"DONE  {section_name} -> in_memory")
        return
    write_parquet_dataset(df, output_path)
    log(f"DONE  {section_name} -> {output_path}")


def _build_coloc_output_df(spark: SparkSession, genes, coloc_threshold: float, coloc_raw, credible_set_raw, study_raw):
    ensure_columns_present(coloc_raw, prep_coloc_COLOC_COLUMNS, "Colocalisation")
    ensure_columns_present(credible_set_raw, prep_coloc_CREDIBLE_SET_COLUMNS, "Credible set")
    ensure_columns_present(study_raw, prep_coloc_STUDY_COLUMNS, "Study")

    coloc_df = coloc_raw.select(*prep_coloc_COLOC_COLUMNS)
    coloc_filtered_df = coloc_df.filter(F.col("h4") > F.lit(coloc_threshold))

    if prep_coloc_ALLOWED_RIGHT_STUDY_TYPES:
        allowed_right_types = [x.lower() for x in prep_coloc_ALLOWED_RIGHT_STUDY_TYPES]
        coloc_filtered_df = coloc_filtered_df.filter(
            F.lower(F.col("rightStudyType")).isin(allowed_right_types)
        )

    genes_df = spark.createDataFrame([(gene.upper(),) for gene in genes], ["rightGeneId"])

    study_filtered_df = (
        study_raw.select(*prep_coloc_STUDY_COLUMNS)
        .withColumn("rightGeneId", F.upper(F.col("geneId")))
        .join(F.broadcast(genes_df), "rightGeneId", "inner")
        .select(
            F.col("studyId").alias("rightStudyId"),
            F.col("rightGeneId"),
            F.col("biosampleId").alias("rightBiosampleId"),
        )
        .dropna(subset=["rightStudyId", "rightGeneId"])
        .dropDuplicates()
    )

    credible_set_df = credible_set_raw.select(*prep_coloc_CREDIBLE_SET_COLUMNS)

    right_loci_df = coloc_filtered_df.select("rightStudyLocusId").dropDuplicates(["rightStudyLocusId"])

    right_credible_set_df = (
        right_loci_df.join(
            credible_set_df.select(
                F.col("studyLocusId").alias("rightStudyLocusId"),
                F.col("studyId").alias("rightStudyId"),
                F.col("variantId").alias("rightVariantId"),
            ),
            "rightStudyLocusId",
            "inner",
        )
        .join(F.broadcast(study_filtered_df), "rightStudyId", "inner")
    )

    coloc_with_right_df = coloc_filtered_df.join(right_credible_set_df, "rightStudyLocusId", "inner")

    left_loci_df = coloc_with_right_df.select("leftStudyLocusId").dropDuplicates(["leftStudyLocusId"])
    left_credible_selected_df = credible_set_df.select(
            F.col("studyLocusId").alias("leftStudyLocusId"),
            F.col("studyType").alias("leftStudyType"),
            F.col("studyId").alias("leftStudyId"),
            F.col("variantId").alias("leftVariantId"),
        )
    if prep_coloc_FILTER_LEFT_STUDY_TYPE_TO_GWAS:
        left_credible_selected_df = left_credible_selected_df.filter(
            F.lower(F.col("leftStudyType")) == F.lit("gwas")
        )

    left_credible_set_df = (
        left_credible_selected_df
        .dropDuplicates(["leftStudyLocusId"])
        .join(left_loci_df, "leftStudyLocusId", "inner")
    )

    joined_df = coloc_with_right_df.join(left_credible_set_df, "leftStudyLocusId", "inner")

    if prep_coloc_FILTER_QTL_TO_BLOOD_BIOSAMPLE_IDS:
        joined_df = joined_df.filter(
            F.col("rightBiosampleId").isin(prep_coloc_BLOOD_BIOSAMPLE_IDS)
        )

    return (
        joined_df.select(
            "leftStudyLocusId",
            "rightStudyLocusId",
            "rightStudyType",
            "chromosome",
            "colocalisationMethod",
            "numberColocalisingVariants",
            "h3",
            "h4",
            "clpp",
            "betaRatioSignAverage",
            "leftStudyType",
            "leftStudyId",
            "leftVariantId",
            "rightStudyId",
            "rightBiosampleId",
            "rightVariantId",
            "rightGeneId",
        )
        .dropDuplicates()
    )


def _build_gwas_output_df(coloc_raw, credible_set_raw):
    ensure_columns_present(coloc_raw, prep_gwas_COLOC_REQUIRED_COLUMNS, "Colocalising variants")
    ensure_columns_present(credible_set_raw, prep_gwas_CREDIBLE_SET_REQUIRED_COLUMNS, "Credible set")

    credible_schema = credible_set_raw.schema

    coloc_gwas = (
        coloc_raw.select("leftStudyLocusId", "leftVariantId", "leftStudyType")
        .filter(F.lower(F.col("leftStudyType")) == F.lit("gwas"))
        .dropna(subset=["leftStudyLocusId", "leftVariantId"])
        .dropDuplicates(["leftStudyLocusId", "leftVariantId"])
    )

    credible_selected = (
        credible_set_raw
        .select(
            F.col("studyLocusId").cast("string").alias("leftStudyLocusId"),
            F.col("variantId").cast("string").alias("variantId"),
            optional_col("beta", credible_schema, "double").alias("beta"),
            optional_col("standardError", credible_schema, "double").alias("standardError"),
            optional_col("zScore", credible_schema, "double").alias("zScore"),
            optional_col("pValueMantissa", credible_schema, "double").alias("pValueMantissa"),
            optional_col("pValueExponent", credible_schema, "long").alias("pValueExponent"),
            optional_col("studyType", credible_schema, "string").alias("studyType"),
            optional_col("studyId", credible_schema, "string").alias("studyId"),
        )
    )

    joined = (
        coloc_gwas.join(
            credible_selected,
            on=(
                (coloc_gwas["leftStudyLocusId"] == credible_selected["leftStudyLocusId"])
                & (coloc_gwas["leftVariantId"] == credible_selected["variantId"])
            ),
            how="inner",
        )
        .select(
            credible_selected["leftStudyLocusId"].alias("studyLocusId"),
            credible_selected["variantId"],
            credible_selected["beta"],
            credible_selected["standardError"],
            credible_selected["zScore"],
            credible_selected["pValueMantissa"],
            credible_selected["pValueExponent"],
            credible_selected["studyType"],
            credible_selected["studyId"],
        )
        .dropDuplicates()
    )

    return (
        joined.withColumn(
            "standardError",
            F.coalesce(
                F.col("standardError"),
                F.when(
                    F.col("beta").isNotNull()
                    & F.col("zScore").isNotNull()
                    & (F.col("zScore") != F.lit(0.0)),
                    F.col("beta") / F.col("zScore"),
                ),
            ),
        )
        .filter(F.col("beta").isNotNull() & F.col("standardError").isNotNull())
        .drop("zScore")
    )


def _build_coding_output_df(spark: SparkSession, genes, credible_set_raw, variant_raw, existing_gwas_raw):
    ensure_columns_present(credible_set_raw, prep_coding_CREDIBLE_SET_REQUIRED_COLUMNS, "Credible set")
    ensure_columns_present(variant_raw, prep_coding_VARIANT_REQUIRED_COLUMNS, "Variant")
    ensure_columns_present(existing_gwas_raw, prep_coding_EXISTING_GWAS_REQUIRED_COLUMNS, "Existing GWAS variants")

    variant_schema = variant_raw.schema
    credible_schema = credible_set_raw.schema
    transcript_consequences_dtype = variant_schema["transcriptConsequences"].dataType
    if not _array_struct_has_fields(
        transcript_consequences_dtype,
        ["targetId", "biotype", "consequenceScore"],
    ):
        raise ValueError(
            "Variant parquet column `transcriptConsequences` must be "
            "array<struct<...targetId..., ...biotype..., ...consequenceScore...>>."
        )

    existing_gwas_variant_ids_df = (
        existing_gwas_raw
        .select(F.col("variantId").cast("string").alias("variantId"))
        .dropna(subset=["variantId"])
        .dropDuplicates()
    )

    credible_gwas_df = (
        credible_set_raw.select(
            F.col("studyLocusId").cast("string").alias("studyLocusId"),
            F.col("variantId").cast("string").alias("variantId"),
            optional_col("beta", credible_schema, "double").alias("beta"),
            optional_col("standardError", credible_schema, "double").alias("standardError"),
            optional_col("zScore", credible_schema, "double").alias("zScore"),
            optional_col("pValueMantissa", credible_schema, "double").alias("pValueMantissa"),
            optional_col("pValueExponent", credible_schema, "long").alias("pValueExponent"),
            optional_col("studyType", credible_schema, "string").alias("studyType"),
            optional_col("studyId", credible_schema, "string").alias("studyId"),
        )
        .filter(F.lower(F.col("studyType")) == F.lit("gwas"))
        .dropna(subset=["studyLocusId", "variantId"])
        .withColumn(
            "pValue",
            F.when(
                F.col("pValueMantissa").isNotNull() & F.col("pValueExponent").isNotNull(),
                F.col("pValueMantissa").cast("double")
                * F.pow(F.lit(10.0), F.col("pValueExponent").cast("double")),
            ),
        )
        .filter(F.col("pValue").isNotNull() & (F.col("pValue") <= F.lit(prep_coding_GWAS_PVALUE_THRESHOLD)))
        .drop("pValue")
        .join(existing_gwas_variant_ids_df, "variantId", "left_anti")
    )

    candidate_variant_ids_df = credible_gwas_df.select("variantId").dropDuplicates()
    genes_df = spark.createDataFrame([(gene.upper(),) for gene in genes], ["geneId"])
    empty_transcript_consequences = F.lit([]).cast(transcript_consequences_dtype)
    consequence_match = _coding_consequence_match_column()

    variant_gene_df = (
        variant_raw.select(
            F.col("variantId").cast("string").alias("variantId"),
            F.col("mostSevereConsequenceId").cast("string").alias("mostSevereConsequenceId"),
            F.col("variantDescription").cast("string").alias("variantDescription"),
            F.col("transcriptConsequences"),
        )
        .join(candidate_variant_ids_df, "variantId", "inner")
        .withColumn(
            "mostSevereConsequenceId",
            F.upper(F.regexp_replace(F.col("mostSevereConsequenceId"), ":", "_")),
        )
        .filter(F.col("mostSevereConsequenceId").isin(prep_coding_CODING_MOST_SEVERE_CONSEQUENCE_IDS))
        .withColumn(
            "_candidateTranscriptConsequences",
            F.coalesce(
                F.filter(
                    F.col("transcriptConsequences"),
                    lambda tc: _coding_transcript_consequence_is_eligible(tc),
                ),
                empty_transcript_consequences,
            ),
        )
        .select(
            "variantId",
            "mostSevereConsequenceId",
            "variantDescription",
            F.explode("_candidateTranscriptConsequences").alias("_transcriptConsequence"),
        )
        .withColumn("geneId", F.upper(F.col("_transcriptConsequence.targetId").cast("string")))
        .withColumn("consequenceScore", F.col("_transcriptConsequence.consequenceScore").cast("double"))
        .drop("_transcriptConsequence")
        .dropna(subset=["variantId", "geneId", "consequenceScore"])
        .join(F.broadcast(genes_df), "geneId", "inner")
        .withColumn("maxConsequenceScore", F.max("consequenceScore").over(Window.partitionBy("variantId")))
        .filter(F.col("consequenceScore") == F.col("maxConsequenceScore"))
        .select("variantId", "geneId", "mostSevereConsequenceId", "variantDescription")
        .dropDuplicates()
    )

    joined = (
        credible_gwas_df.join(variant_gene_df, "variantId", "inner")
        .select(
            "studyLocusId",
            "variantId",
            "beta",
            "standardError",
            "zScore",
            "pValueMantissa",
            "pValueExponent",
            "studyType",
            "studyId",
            "geneId",
            "mostSevereConsequenceId",
            "variantDescription",
        )
        .dropDuplicates()
    )

    return (
        joined.withColumn(
            "standardError",
            F.coalesce(
                F.col("standardError"),
                F.when(
                    F.col("beta").isNotNull()
                    & F.col("zScore").isNotNull()
                    & (F.col("zScore") != F.lit(0.0)),
                    F.col("beta") / F.col("zScore"),
                ),
            ),
        )
        .filter(F.col("beta").isNotNull() & F.col("standardError").isNotNull())
        .drop("zScore")
    )


def _build_qtl_output_df(coloc_raw, credible_set_raw):
    ensure_columns_present(coloc_raw, prep_qtl_COLOC_REQUIRED_COLUMNS, "Colocalising variants")
    ensure_columns_present(credible_set_raw, prep_qtl_CREDIBLE_SET_REQUIRED_COLUMNS, "Credible set")

    credible_schema = credible_set_raw.schema

    coloc_qtl = (
        coloc_raw.select(
            "rightStudyLocusId",
            "rightVariantId",
            "rightGeneId",
            "rightBiosampleId",
        )
        .dropna(subset=["rightStudyLocusId", "rightVariantId"])
        .dropDuplicates(["rightStudyLocusId", "rightVariantId", "rightGeneId", "rightBiosampleId"])
    )

    credible_selected = (
        credible_set_raw
        .select(
            F.col("studyLocusId").cast("string").alias("rightStudyLocusId"),
            F.col("variantId").cast("string").alias("variantId"),
            optional_col("beta", credible_schema, "double").alias("beta"),
            optional_col("standardError", credible_schema, "double").alias("standardError"),
            optional_col("zScore", credible_schema, "double").alias("zScore"),
            optional_col("pValueMantissa", credible_schema, "double").alias("pValueMantissa"),
            optional_col("pValueExponent", credible_schema, "long").alias("pValueExponent"),
            optional_col("studyType", credible_schema, "string").alias("studyType"),
            optional_col("studyId", credible_schema, "string").alias("studyId"),
        )
    )

    joined = (
        coloc_qtl.join(
            credible_selected,
            on=(
                (coloc_qtl["rightStudyLocusId"] == credible_selected["rightStudyLocusId"])
                & (coloc_qtl["rightVariantId"] == credible_selected["variantId"])
            ),
            how="inner",
        )
        .select(
            credible_selected["rightStudyLocusId"].alias("studyLocusId"),
            credible_selected["variantId"],
            credible_selected["beta"],
            credible_selected["standardError"],
            credible_selected["zScore"],
            credible_selected["pValueMantissa"],
            credible_selected["pValueExponent"],
            credible_selected["studyType"],
            credible_selected["studyId"],
            coloc_qtl["rightGeneId"].cast("string").alias("geneId"),
            coloc_qtl["rightBiosampleId"].cast("string").alias("biosampleId"),
        )
        .dropDuplicates()
    )

    return (
        joined.withColumn(
            "standardError",
            F.coalesce(
                F.col("standardError"),
                F.when(
                    F.col("beta").isNotNull()
                    & F.col("zScore").isNotNull()
                    & (F.col("zScore") != F.lit(0.0)),
                    F.col("beta") / F.col("zScore"),
                ),
            ),
        )
        .filter(F.col("beta").isNotNull() & F.col("standardError").isNotNull())
        .drop("zScore")
    )


def _build_burden_output_df(spark: SparkSession, genes, raw):
    schema = raw.schema
    if "targetId" not in schema.fieldNames():
        raise ValueError("Input burden evidence parquet is missing required column: targetId")

    filtered = raw

    if "datasourceId" in schema.fieldNames():
        filtered = filtered.filter(F.col("datasourceId") == F.lit("gene_burden"))
    if prep_burden_FILTER_TO_ALLOWED_PROJECTS and "projectId" in schema.fieldNames():
        filtered = filtered.filter(F.col("projectId").isin(prep_burden_ALLOWED_PROJECT_IDS))
    if prep_burden_FILTER_TO_ALLOWED_COHORTS and "cohortId" in schema.fieldNames():
        filtered = filtered.filter(F.col("cohortId").isin(prep_burden_ALLOWED_COHORT_IDS))
    if prep_burden_FILTER_TO_ALLOWED_METHODS and "statisticalMethod" in schema.fieldNames():
        filtered = filtered.filter(F.col("statisticalMethod").isin(prep_burden_ALLOWED_STAT_METHODS))

    genes_df = spark.createDataFrame([(gene,) for gene in genes], ["targetId"])
    filtered = filtered.join(F.broadcast(genes_df), "targetId", "inner")

    selected = filtered.select(
        optional_col("datasourceId", schema, "string").alias("datasourceId"),
        optional_col("targetId", schema, "string").alias("targetId"),
        prep_burden_stringified_or_null("allelicRequirements", schema).alias("allelicRequirements"),
        optional_col("ancestryId", schema, "string").alias("ancestryId"),
        optional_col("beta", schema, "double").alias("beta"),
        optional_col("betaConfidenceIntervalLower", schema, "double").alias("betaConfidenceIntervalLower"),
        optional_col("betaConfidenceIntervalUpper", schema, "double").alias("betaConfidenceIntervalUpper"),
        optional_col("cohortId", schema, "string").alias("cohortId"),
        optional_col("diseaseFromSource", schema, "string").alias("diseaseFromSource"),
        optional_col("diseaseFromSourceMappedId", schema, "string").alias("diseaseFromSourceMappedId"),
        optional_col("oddsRatio", schema, "double").alias("oddsRatio"),
        optional_col("oddsRatioConfidenceIntervalLower", schema, "double").alias("oddsRatioConfidenceIntervalLower"),
        optional_col("oddsRatioConfidenceIntervalUpper", schema, "double").alias("oddsRatioConfidenceIntervalUpper"),
        optional_col("pValueExponent", schema, "long").alias("pValueExponent"),
        optional_col("pValueMantissa", schema, "double").alias("pValueMantissa"),
        optional_col("projectId", schema, "string").alias("projectId"),
        optional_col("resourceScore", schema, "double").alias("resourceScore"),
        optional_col("statisticalMethod", schema, "string").alias("statisticalMethod"),
        optional_col("targetFromSourceId", schema, "string").alias("targetFromSourceId"),
        optional_col("score", schema, "double").alias("score"),
        optional_col("studyCases", schema, "long").alias("studyCases"),
        optional_col("studyCasesWithQualifyingVariants", schema, "long").alias("studyCasesWithQualifyingVariants"),
        optional_col("studySampleSize", schema, "long").alias("studySampleSize"),
        F.lit(0).cast("int").alias("burdenSourceRank"),
    )

    return (
        selected
        .withColumn("oddsRatio", F.coalesce(F.col("oddsRatio"), F.exp(F.col("beta"))))
        .withColumn(
            "oddsRatioConfidenceIntervalLower",
            F.coalesce(F.col("oddsRatioConfidenceIntervalLower"), F.exp(F.col("betaConfidenceIntervalLower"))),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalUpper",
            F.coalesce(F.col("oddsRatioConfidenceIntervalUpper"), F.exp(F.col("betaConfidenceIntervalUpper"))),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalLower",
            F.when(
                F.col("oddsRatioConfidenceIntervalLower").isNull()
                & F.col("oddsRatio").isNotNull()
                & F.col("oddsRatioConfidenceIntervalUpper").isNotNull()
                & (F.col("oddsRatioConfidenceIntervalUpper") > F.lit(0.0)),
                (F.col("oddsRatio") * F.col("oddsRatio")) / F.col("oddsRatioConfidenceIntervalUpper"),
            ).otherwise(F.col("oddsRatioConfidenceIntervalLower")),
        )
        .filter(F.col("oddsRatio").isNotNull())
        .drop("beta", "betaConfidenceIntervalLower", "betaConfidenceIntervalUpper")
    )


def _build_raw_az_burden_output_df(spark: SparkSession, args, genes):
    # Raw AZ burden is opt-in enrichment on top of release burden evidence.
    # We restrict to the requested genes before the phenotype/EFO join so the
    # large raw tables stay as small as possible for both local and gcloud runs.
    target_map = _load_target_symbol_map(spark, args.target_data_dir)
    efo_map = _load_az_burden_efo_map(spark, args.az_mapping_file)
    genes_df = spark.createDataFrame([(gene,) for gene in genes], ["targetId"])
    requested_target_map = (
        target_map
        .join(F.broadcast(genes_df), "targetId", "inner")
        .select("targetId", "gene_symbol_norm")
        .dropDuplicates(["targetId", "gene_symbol_norm"])
    )

    raw = spark.read.option("mergeSchema", "true").parquet(
        resolve_parquet_input(args.az_burden_binary_dir),
        resolve_parquet_input(args.az_burden_quantitative_dir),
    )
    ensure_columns_present(
        raw,
        [
            "Gene",
            "Phenotype",
            "CollapsingModel",
            "Type",
            "pValue",
            "beta",
            "LCI",
            "UCI",
            "BinOddsRatio",
            "BinOddsRatioLCI",
            "BinOddsRatioUCI",
        ],
        "raw AZ burden",
    )

    prepared = (
        raw.select(
            "Gene",
            "Phenotype",
            "CollapsingModel",
            "Type",
            "pValue",
            "beta",
            "LCI",
            "UCI",
            "BinOddsRatio",
            "BinOddsRatioLCI",
            "BinOddsRatioUCI",
        )
        .withColumn("gene_symbol_norm", F.upper(F.trim(F.regexp_replace(F.col("Gene").cast("string"), "'", ""))))
        .withColumn("CollapsingModel", F.lower(F.trim(F.col("CollapsingModel").cast("string"))))
        .withColumn("pValue", F.col("pValue").cast("double"))
        .filter(F.col("CollapsingModel").isin(AZ_BURDEN_COLLAPSING_MODELS))
        .filter(F.col("pValue").isNotNull() & (F.col("pValue") >= F.lit(0.0)))
        # The HGNC -> Ensembl map is applied before phenotype mapping on purpose:
        # it is the cheapest early filter and keeps memory usage reasonable when
        # the optional raw AZ burden source is enabled on large runs.
        .join(F.broadcast(requested_target_map), "gene_symbol_norm", "inner")
        .withColumn("Phenotype", F.trim(F.col("Phenotype").cast("string")))
        .withColumn("Type", F.trim(F.col("Type").cast("string")))
        .join(efo_map, F.col("Phenotype") == efo_map["diseaseFromSource"], "inner")
        .drop("diseaseFromSource")
    )

    quantitative = F.lower(F.col("Type")) == F.lit("quantitative")
    with_effect = (
        prepared
        .withColumn(
            "oddsRatio",
            F.when(quantitative, F.exp(F.col("beta").cast("double")))
            .otherwise(F.col("BinOddsRatio").cast("double")),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalLower",
            F.when(quantitative, F.exp(F.col("LCI").cast("double")))
            .otherwise(F.col("BinOddsRatioLCI").cast("double")),
        )
        .withColumn(
            "oddsRatioConfidenceIntervalUpper",
            F.when(quantitative, F.exp(F.col("UCI").cast("double")))
            .otherwise(F.col("BinOddsRatioUCI").cast("double")),
        )
        .filter(
            F.col("oddsRatio").isNotNull()
            & (F.col("oddsRatio") > F.lit(0.0))
            & F.col("oddsRatioConfidenceIntervalLower").isNotNull()
            & (F.col("oddsRatioConfidenceIntervalLower") > F.lit(0.0))
            & F.col("oddsRatioConfidenceIntervalUpper").isNotNull()
            & (F.col("oddsRatioConfidenceIntervalUpper") > F.lit(0.0))
            & (F.col("oddsRatioConfidenceIntervalUpper") > F.col("oddsRatioConfidenceIntervalLower"))
        )
        .withColumn(
            "_pvalue_for_log",
            F.when(F.col("pValue") > F.lit(0.0), F.col("pValue")).otherwise(F.lit(1e-300)),
        )
        .withColumn("pValueExponent", F.floor(F.log10(F.col("_pvalue_for_log"))).cast("long"))
        .withColumn(
            "pValueMantissa",
            F.when(
                F.col("pValue") > F.lit(0.0),
                F.col("pValue") / F.exp(F.log(F.lit(10.0)) * F.col("pValueExponent").cast("double")),
            )
            .otherwise(F.lit(0.0)),
        )
    )

    return with_effect.select(
        F.lit("raw_az_phewas_burden").cast("string").alias("datasourceId"),
        F.col("targetId").cast("string").alias("targetId"),
        F.lit(None).cast("string").alias("allelicRequirements"),
        F.lit(None).cast("string").alias("ancestryId"),
        F.lit("UK Biobank 470k").cast("string").alias("cohortId"),
        F.col("Phenotype").cast("string").alias("diseaseFromSource"),
        F.col("diseaseId").cast("string").alias("diseaseFromSourceMappedId"),
        F.col("oddsRatio").cast("double").alias("oddsRatio"),
        F.col("oddsRatioConfidenceIntervalLower").cast("double").alias("oddsRatioConfidenceIntervalLower"),
        F.col("oddsRatioConfidenceIntervalUpper").cast("double").alias("oddsRatioConfidenceIntervalUpper"),
        F.col("pValueExponent").cast("long").alias("pValueExponent"),
        F.col("pValueMantissa").cast("double").alias("pValueMantissa"),
        F.lit("AstraZeneca PheWAS Portal").cast("string").alias("projectId"),
        F.lit(None).cast("double").alias("resourceScore"),
        F.col("CollapsingModel").cast("string").alias("statisticalMethod"),
        F.col("Gene").cast("string").alias("targetFromSourceId"),
        F.lit(None).cast("double").alias("score"),
        F.lit(None).cast("long").alias("studyCases"),
        F.lit(None).cast("long").alias("studyCasesWithQualifyingVariants"),
        F.lit(None).cast("long").alias("studySampleSize"),
        F.lit(1).cast("int").alias("burdenSourceRank"),
    )


def _build_clinvar_output_df(spark: SparkSession, genes, raw):
    schema = raw.schema
    filtered = (
        raw
        .filter(F.col("datasourceId") == F.lit("eva"))
        .filter(~F.col("confidence").isin(prep_clinvar_EXCLUDED_CONFIDENCE))
        .filter(F.col("variantId").isNotNull())
    )

    genes_df = spark.createDataFrame([(gene.upper(),) for gene in genes], ["targetId_upper"])
    filtered = (
        filtered.withColumn("targetId_upper", F.upper(F.col("targetId")))
        .join(F.broadcast(genes_df), "targetId_upper", "inner")
        .drop("targetId_upper")
    )

    # We keep the VIDRA_2 decision to take the first clinical significance value
    # from OT EVA evidence rather than exploding one evidence row into many rows.
    with_clin_sig = prep_clinvar_add_clinical_significance_column(filtered, schema)

    return with_clin_sig.select(
        F.col("_clinicalSignificance").cast("string").alias("clinicalSignificances"),
        optional_col("confidence", schema, "string").alias("confidence"),
        optional_col("diseaseFromSourceId", schema, "string").alias("diseaseFromSourceId"),
        optional_col("diseaseFromSourceMappedId", schema, "string").alias("diseaseFromSourceMappedId"),
        optional_col("studyId", schema, "string").alias("studyId"),
        optional_col("targetFromSourceId", schema, "string").alias("targetFromSourceId"),
        optional_col("targetId", schema, "string").alias("targetId"),
        optional_col("variantFunctionalConsequenceId", schema, "string").alias("variantFunctionalConsequenceId"),
        optional_col("variantHgvsId", schema, "string").alias("variantHgvsId"),
        optional_col("variantId", schema, "string").alias("variantId"),
        optional_col("variantRsId", schema, "string").alias("variantRsId"),
        optional_col("diseaseId", schema, "string").alias("diseaseId"),
        optional_col("score", schema, "double").alias("score"),
        optional_col("datasourceId", schema, "string").alias("datasourceId"),
    )


def _build_az_ready_dataset_from_inputs(spark: SparkSession, args, genes, variant_index):
    if not path_exists(args.az_variants_file, getattr(args, "gcp_project", "")):
        raise FileNotFoundError(f"AZ variants file not found: {args.az_variants_file}")
    if not path_exists(args.az_mapping_file, getattr(args, "gcp_project", "")):
        raise FileNotFoundError(f"AZ mapping file not found: {args.az_mapping_file}")
    if not path_exists(args.az_gene_map_file, getattr(args, "gcp_project", "")):
        raise FileNotFoundError(f"AZ gene map file not found: {args.az_gene_map_file}")

    genes_df = spark.createDataFrame([(gene,) for gene in genes], ["ensembl_id"])
    gene_map = _load_az_gene_map(spark, args.az_gene_map_file)
    efo_map = _load_az_efo_map(spark, args.az_mapping_file)

    az_raw = (
        spark.read
        .option("header", "true")
        .schema(AZ_SCHEMA)
        .csv(args.az_variants_file)
    )

    az_filtered = (
        az_raw
        .filter(F.col("Model") == F.lit("allelic"))
        .filter(F.col("p-value") < F.lit(AZ_PVAL_THRESHOLD))
        .withColumn("Gene", F.regexp_replace(F.col("Gene"), "'", ""))
        .withColumn(
            "variant",
            F.regexp_replace(F.regexp_replace(F.trim(F.col("Variant")), " ", "_"), "-", "_"),
        )
        .filter(F.col("variant").isNotNull())
        .filter(F.col("Odds ratio").isNotNull() & (F.col("Odds ratio") > F.lit(0.0)))
        .filter(F.col("Odds ratio LCI").isNotNull() & (F.col("Odds ratio LCI") > F.lit(0.0)))
        .filter(F.col("Odds ratio UCI").isNotNull() & (F.col("Odds ratio UCI") > F.lit(0.0)))
    )

    az_mapped = (
        az_filtered
        .join(gene_map, az_filtered["Gene"] == gene_map["gene_symbol"], "inner")
        .drop("gene_symbol")
        .join(F.broadcast(genes_df), "ensembl_id", "inner")
        .join(efo_map, F.col("Phenotype") == efo_map["diseaseFromSource"], "left")
        .drop("diseaseFromSource")
        .filter(F.col("diseaseId").isNotNull())
    )

    az_with_match = (
        az_mapped
        .join(variant_index, az_mapped["variant"] == variant_index["ot_variant"], "left")
        .withColumn("is_variant_matched", F.col("ot_variant").isNotNull())
    )

    # The summary is intentionally printed because AZ variants are matched
    # against the full release OT `variant.variantId` universe. This makes the
    # match rate stable across release-mode and lead-mode Step 1 runs.
    summary = az_with_match.agg(
        F.count("*").alias("rows"),
        F.sum(F.when(F.col("is_variant_matched"), F.lit(1)).otherwise(F.lit(0))).alias("matched"),
        F.sum(F.when(~F.col("is_variant_matched"), F.lit(1)).otherwise(F.lit(0))).alias("unmatched"),
    ).collect()[0]
    unmatched_sample = [
        row["variant"]
        for row in az_with_match
        .filter(~F.col("is_variant_matched"))
        .select("variant")
        .dropDuplicates()
        .limit(5)
        .collect()
    ]
    log(
        "AZ variant normalization summary against release variant universe: "
        f"rows={summary['rows']}, matched={summary['matched']}, "
        f"unmatched={summary['unmatched']}, sample_unmatched={unmatched_sample}"
    )

    az_ready = (
        az_with_match
        .withColumn("yc", F.log(F.col("Odds ratio")))
        .withColumn(
            "ycse",
            (F.log(F.col("Odds ratio UCI")) - F.log(F.col("Odds ratio LCI"))) / F.lit(3.92),
        )
        .withColumn(
            "ycse",
            F.when(F.col("ycse").isNull() | (F.col("ycse") <= F.lit(0.0)), F.lit(0.14)).otherwise(F.col("ycse")),
        )
        .select(
            F.col("variant").cast("string").alias("variant"),
            F.col("ensembl_id").cast("string").alias("as_gene"),
            F.col("diseaseId").cast("string").alias("as_disease"),
            F.lit(1).cast("int").alias("GsourceLab"),
            F.lit(2).cast("int").alias("GqtlLab"),
            F.col("yc").cast("double").alias("yc"),
            F.col("ycse").cast("double").alias("ycse"),
            F.lit(0.0).cast("double").alias("xc"),
            F.lit(0.1).cast("double").alias("xcse"),
            F.lit(0.0).cast("double").alias("as_clinicalSignificance"),
            F.col("Variant").cast("string").alias("sourceVariant"),
            F.col("Gene").cast("string").alias("geneSymbol"),
            F.col("Phenotype").cast("string").alias("diseaseFromSource"),
            F.col("ot_variant").cast("string").alias("matchedOtVariantId"),
            F.col("is_variant_matched").cast("boolean").alias("is_variant_matched"),
            F.col("p-value").cast("double").alias("pValue"),
            F.col("Odds ratio").cast("double").alias("oddsRatio"),
            F.col("Odds ratio LCI").cast("double").alias("oddsRatioConfidenceIntervalLower"),
            F.col("Odds ratio UCI").cast("double").alias("oddsRatioConfidenceIntervalUpper"),
        )
        .dropna(subset=["variant", "as_gene", "as_disease"])
    )
    return _deduplicate_az_ready(az_ready)

def build_az_dataset(args) -> None:
    log("START az_variants")

    spark = build_local_spark("vidra_az_variants")
    try:
        az_ready = _build_az_ready_dataset(spark, args)
        write_parquet_dataset(az_ready, args.az_output)
        log(f"DONE  az_variants -> {args.az_output}")
    finally:
        try:
            spark.stop()
        except Exception:
            pass




def build_analysis_ready_dataset(
    spark: SparkSession,
    args,
    *,
    study_raw,
    coloc,
    gwas,
    qtl,
    coding,
    clinvar,
    burden,
    az_ready,
) -> None:
    log("START analysis_ready_assembly")

    # Final assembly is source-aware by construction. Each branch writes the
    # same public output schema but keeps its own `(GsourceLab, GqtlLab)` code,
    # and all dedup windows below include those source labels so rows from
    # different public evidence sources are never collapsed together.

    # Enforce the same core schema contract as working/vidra pre-processing.
    ensure_columns_present(
        coloc,
        ["leftStudyLocusId", "leftStudyId", "leftVariantId", "rightStudyLocusId", "rightVariantId", "rightGeneId", "rightStudyType"],
        "colocalising_variants",
    )
    ensure_columns_present(
        gwas,
        ["studyLocusId", "variantId", "beta", "standardError", "studyId"],
        "GWAS_variants",
    )
    ensure_columns_present(
        qtl,
        ["studyLocusId", "variantId", "geneId", "beta", "standardError"],
        "QTL_variants",
    )
    ensure_columns_present(
        coding,
        ["variantId", "geneId", "studyId", "beta", "standardError"],
        "GWAS_coding_variants_from_cs",
    )
    ensure_columns_present(
        clinvar,
        ["variantId", "targetId", "clinicalSignificances", "diseaseFromSourceMappedId", "score"],
        "clinvar_variants",
    )
    ensure_columns_present(
        burden,
        [
            "targetId",
            "diseaseFromSourceMappedId",
            "oddsRatio",
            "oddsRatioConfidenceIntervalLower",
            "oddsRatioConfidenceIntervalUpper",
            "pValueMantissa",
            "pValueExponent",
            "projectId",
            "cohortId",
            "statisticalMethod",
        ],
        "burden_tests",
    )
    ensure_columns_present(
        az_ready,
        ANALYSIS_SOURCE_COLS,
        "AZ_variants",
    )

    # `study_raw` is polymorphic:
    # - release mode: the OT study parquet
    # - lead mode: a minimal projection from lead_variant_effect carrying only
    #   studyId + disease mapping arrays
    # `_study_disease_mapping` hides that difference and enforces the same
    # "mapped trait first, diseaseIds fallback" rule for both modes.
    study_map = _study_disease_mapping(study_raw)

    gwas_schema = gwas.schema
    qtl_schema = qtl.schema
    coding_schema = coding.schema

    gwas_common = gwas.select(
        F.col("studyLocusId").alias("gwasStudyLocusId"),
        F.col("variantId").alias("gwasVariantId"),
        F.col("beta").alias("gwas_beta"),
        F.col("standardError").alias("gwas_se"),
        F.col("studyType").alias("gwas_studyType"),
        F.col("studyId").alias("gwas_studyId"),
        optional_col("pValueMantissa", gwas_schema, "double").alias("gwasPValueMantissa"),
        optional_col("pValueExponent", gwas_schema, "long").alias("gwasPValueExponent"),
    )
    qtl_common = qtl.select(
        F.col("studyLocusId").alias("qtlStudyLocusId"),
        F.col("variantId").alias("qtlVariantId"),
        F.col("geneId").alias("qtlGeneId"),
        F.col("beta").alias("qtl_beta"),
        F.col("standardError").alias("qtl_se"),
        F.col("studyType").alias("qtl_studyType"),
        F.col("studyId").alias("qtl_studyId"),
        optional_col("pValueMantissa", qtl_schema, "double").alias("qtlPValueMantissa"),
        optional_col("pValueExponent", qtl_schema, "long").alias("qtlPValueExponent"),
    )

    common = (
        coloc.select(
            "leftStudyLocusId",
            "leftStudyId",
            "leftVariantId",
            "rightStudyLocusId",
            "rightVariantId",
            "rightGeneId",
            "rightStudyType",
        )
        .join(
            gwas_common,
            (F.col("leftStudyLocusId") == gwas_common["gwasStudyLocusId"])
            & (F.col("leftVariantId") == gwas_common["gwasVariantId"]),
            "inner",
        )
        .join(
            qtl_common,
            (F.col("rightStudyLocusId") == qtl_common["qtlStudyLocusId"])
            & (F.col("rightVariantId") == qtl_common["qtlVariantId"])
            & (F.col("rightGeneId") == qtl_common["qtlGeneId"]),
            "inner",
        )
    )

    common = common.join(
        F.broadcast(study_map.withColumnRenamed("studyId", "mapStudyId")),
        F.col("leftStudyId") == F.col("mapStudyId"),
        "left",
    )

    common_candidates = (
        common.select(
            F.coalesce(F.col("leftVariantId"), F.col("rightVariantId")).cast("string").alias("variant"),
            F.col("rightGeneId").cast("string").alias("as_gene"),
            F.coalesce(F.col("mappedDisease"), F.col("leftStudyId")).cast("string").alias("as_disease"),
            F.lit(0).cast("int").alias("GsourceLab"),
            F.when(
                F.lower(F.coalesce(F.col("rightStudyType"), F.col("qtl_studyType"))).contains("pqtl"),
                F.lit(1),
            ).otherwise(F.lit(0)).cast("int").alias("GqtlLab"),
            F.coalesce(F.col("gwas_beta").cast("double"), F.lit(0.0)).alias("yc"),
            F.coalesce(F.col("gwas_se").cast("double"), F.lit(0.14)).alias("ycse"),
            F.coalesce(F.col("qtl_beta").cast("double"), F.lit(0.0)).alias("xc"),
            F.coalesce(F.col("qtl_se").cast("double"), F.lit(0.1)).alias("xcse"),
            F.lit(0.0).cast("double").alias("as_clinicalSignificance"),
            _pvalue_from_mantissa_exponent("gwasPValueMantissa", "gwasPValueExponent").alias("_gwasPValue"),
            _pvalue_from_mantissa_exponent("qtlPValueMantissa", "qtlPValueExponent").alias("_qtlPValue"),
            F.col("leftStudyLocusId").cast("string").alias("_leftStudyLocusId"),
            F.col("rightStudyLocusId").cast("string").alias("_rightStudyLocusId"),
        )
        .dropna(subset=["variant", "as_gene", "as_disease"])
    )
    # Common/coloc duplicates are ranked by the strongest association evidence:
    # lowest GWAS p-value first, then lowest QTL p-value. Locus ids are only
    # deterministic tie-breakers once statistical evidence is equivalent.
    common_window = Window.partitionBy("variant", "as_gene", "as_disease", "GsourceLab", "GqtlLab").orderBy(
        F.col("_gwasPValue").asc_nulls_last(),
        F.col("_qtlPValue").asc_nulls_last(),
        F.col("_leftStudyLocusId").asc_nulls_last(),
        F.col("_rightStudyLocusId").asc_nulls_last(),
    )
    common_ready = (
        common_candidates
        .withColumn("rn", F.row_number().over(common_window))
        .filter(F.col("rn") == 1)
        .drop("rn", "_gwasPValue", "_qtlPValue", "_leftStudyLocusId", "_rightStudyLocusId")
        .dropDuplicates()
    )

    coding_for_ready = coding.select(
        "variantId",
        "geneId",
        "studyId",
        "beta",
        "standardError",
        optional_col("pValueMantissa", coding_schema, "double").alias("pValueMantissa"),
        optional_col("pValueExponent", coding_schema, "long").alias("pValueExponent"),
        optional_col("studyLocusId", coding_schema, "string").alias("studyLocusId"),
    )

    coding_candidates = (
        coding_for_ready
        .join(
            F.broadcast(study_map.withColumnRenamed("studyId", "mapStudyId")),
            coding_for_ready["studyId"] == F.col("mapStudyId"),
            "left",
        )
        .select(
            F.col("variantId").cast("string").alias("variant"),
            F.col("geneId").cast("string").alias("as_gene"),
            F.coalesce(F.col("mappedDisease"), F.col("studyId")).cast("string").alias("as_disease"),
            F.lit(3).cast("int").alias("GsourceLab"),
            F.lit(2).cast("int").alias("GqtlLab"),
            F.coalesce(F.col("beta").cast("double"), F.lit(0.0)).alias("yc"),
            F.coalesce(F.col("standardError").cast("double"), F.lit(0.14)).alias("ycse"),
            F.lit(0.0).cast("double").alias("xc"),
            F.lit(0.1).cast("double").alias("xcse"),
            F.lit(0.0).cast("double").alias("as_clinicalSignificance"),
            _pvalue_from_mantissa_exponent("pValueMantissa", "pValueExponent").alias("_gwasPValue"),
            F.col("studyLocusId").cast("string").alias("_studyLocusId"),
        )
        .dropna(subset=["variant", "as_gene", "as_disease"])
    )
    # Coding duplicates are simpler: one source, no QTL component, so rank by
    # lowest GWAS p-value and then stable locus id.
    coding_window = Window.partitionBy("variant", "as_gene", "as_disease", "GsourceLab", "GqtlLab").orderBy(
        F.col("_gwasPValue").asc_nulls_last(),
        F.col("_studyLocusId").asc_nulls_last(),
    )
    coding_ready = (
        coding_candidates
        .withColumn("rn", F.row_number().over(coding_window))
        .filter(F.col("rn") == 1)
        .drop("rn", "_gwasPValue", "_studyLocusId")
        .dropDuplicates()
    )

    clinvar_candidates = (
        clinvar
        .select(
            F.col("variantId").cast("string").alias("variant"),
            F.col("targetId").cast("string").alias("as_gene"),
            F.coalesce(
                F.col("diseaseFromSourceMappedId"),
                F.col("diseaseId"),
                F.col("studyId"),
            ).cast("string").alias("as_disease"),
            F.lit(2).cast("int").alias("GsourceLab"),
            F.lit(2).cast("int").alias("GqtlLab"),
            F.lit(0.0).cast("double").alias("yc"),
            F.lit(0.14).cast("double").alias("ycse"),
            F.lit(0.0).cast("double").alias("xc"),
            F.lit(0.1).cast("double").alias("xcse"),
            F.coalesce(
                _clinvar_significance_expr(),
                F.col("score").cast("double"),
                F.lit(0.0),
            ).alias("as_clinicalSignificance"),
            F.col("score").cast("double").alias("_clinvarScore"),
        )
        .dropna(subset=["variant", "as_gene", "as_disease"])
    )
    # ClinVar duplicates are resolved deterministically in Step 1 rather than
    # later by row order. We keep the highest EVA `score` first, then the
    # strongest VIDRA ordinal clinical-significance label.
    clinvar_window = Window.partitionBy("variant", "as_gene", "as_disease", "GsourceLab", "GqtlLab").orderBy(
        F.col("_clinvarScore").desc_nulls_last(),
        F.col("as_clinicalSignificance").desc_nulls_last()
    )
    clinvar_ready = (
        clinvar_candidates
        .withColumn("rn", F.row_number().over(clinvar_window))
        .filter(F.col("rn") == 1)
        .drop("rn", "_clinvarScore")
        .dropDuplicates()
    )

    burden_filled = (
        burden
        .withColumn("disease_key", F.col("diseaseFromSourceMappedId").cast("string"))
        .withColumn("oddsRatio_filled", F.col("oddsRatio").cast("double"))
        .withColumn("oddsRatioCI_Lower_filled", F.col("oddsRatioConfidenceIntervalLower").cast("double"))
        .withColumn("oddsRatioCI_Upper_filled", F.col("oddsRatioConfidenceIntervalUpper").cast("double"))
        .filter(
            F.col("targetId").isNotNull()
            & F.col("disease_key").isNotNull()
            & F.col("oddsRatio_filled").isNotNull()
            & (F.col("oddsRatio_filled") > F.lit(0.0))
            & F.col("oddsRatioCI_Lower_filled").isNotNull()
            & (F.col("oddsRatioCI_Lower_filled") > F.lit(0.0))
            & F.col("oddsRatioCI_Upper_filled").isNotNull()
            & (F.col("oddsRatioCI_Upper_filled") > F.lit(0.0))
        )
        .withColumn("bO", F.log(F.col("oddsRatio_filled")))
        .withColumn(
            "bOse",
            (F.log(F.col("oddsRatioCI_Upper_filled")) - F.log(F.col("oddsRatioCI_Lower_filled"))) / F.lit(3.92),
        )
        .withColumn(
            "bOse",
            F.when(F.col("bOse").isNull() | (F.col("bOse") <= F.lit(0.0)), F.lit(2.0)).otherwise(F.col("bOse")),
        )
    )

    # Burden rows are collapsed per gene/disease after release burden and
    # optional raw AZ burden are unioned. Ordering is:
    # 1) lowest p-value exponent
    # 2) lowest p-value mantissa
    # 3) lowest burdenSourceRank, where release=0 and raw AZ=1
    #
    # This means the best p-value wins overall, while exact p-value ties prefer
    # the release burden row.
    burden_order = [F.col("pValueExponent").asc_nulls_last()] if "pValueExponent" in burden_filled.columns else [F.lit(0)]
    if "pValueMantissa" in burden_filled.columns:
        burden_order.append(F.col("pValueMantissa").asc_nulls_last())
    if "burdenSourceRank" in burden_filled.columns:
        burden_order.append(F.col("burdenSourceRank").asc_nulls_last())

    burden_window = Window.partitionBy("targetId", "disease_key").orderBy(*burden_order)
    burden_for_join = (
        burden_filled
        .withColumn("rn", F.row_number().over(burden_window))
        .filter(F.col("rn") == 1)
        .select(
            F.col("targetId").cast("string").alias("b_gene"),
            F.col("disease_key").cast("string").alias("b_disease"),
            F.col("bO").cast("double").alias("bO"),
            F.col("bOse").cast("double").alias("bOse"),
        )
    )

    # Burden is joined onto the union of variant-centric evidence sources rather
    # than emitted as its own row type. This preserves the original Step 1
    # output contract consumed by Step 3.
    all_variants = (
        az_ready.select(*ANALYSIS_SOURCE_COLS)
        .unionByName(common_ready, allowMissingColumns=True)
        .unionByName(coding_ready, allowMissingColumns=True)
        .unionByName(clinvar_ready, allowMissingColumns=True)
    )

    all_with_burden = (
        all_variants.join(
            burden_for_join,
            (all_variants.as_gene == burden_for_join.b_gene)
            & (all_variants.as_disease == burden_for_join.b_disease),
            "left",
        )
        .withColumn("bO", F.coalesce(F.col("bO"), F.lit(ANALYSIS_DEFAULT_FILL["bO"])))
        .withColumn("bOse", F.coalesce(F.col("bOse"), F.lit(ANALYSIS_DEFAULT_FILL["bOse"])))
        .drop("b_gene", "b_disease")
    )

    for col_name, default_val in ANALYSIS_DEFAULT_FILL.items():
        if col_name in all_with_burden.columns:
            all_with_burden = all_with_burden.fillna({col_name: default_val})

    final_df = (
        all_with_burden
        .filter(F.col("variant").isNotNull())
        .filter(F.col("as_gene").isNotNull())
        .filter(F.col("as_disease").isNotNull())
        .select(*ANALYSIS_OUTPUT_COLS)
        .dropDuplicates()
    )

    # Repartition by gene so downstream per-gene reads hit gene-local files and
    # Nextflow Step 3 can consume partitions predictably.
    n_genes = final_df.select("as_gene").distinct().count()
    if n_genes > 0:
        final_df = final_df.repartition(n_genes, "as_gene")

    bucket_base = _resolve_bucket_base(args)
    if bucket_base:
        analysis_path = _join_bucket_path(bucket_base, args.analysis_ready_dir)
        manifest_path = _join_bucket_path(bucket_base, args.analysis_manifest_dir)
    else:
        analysis_path = args.analysis_ready_dir
        manifest_path = args.analysis_manifest_dir

    spark.conf.set("spark.sql.sources.partitionOverwriteMode", "dynamic")
    final_df.write.mode("overwrite").partitionBy("as_gene").parquet(analysis_path)

    manifest_df = final_df.select("variant").distinct()
    if path_exists(manifest_path, args.gcp_project):
        existing_manifest = spark.read.parquet(manifest_path).select("variant").distinct().cache()
        existing_manifest.count()
        manifest_df = existing_manifest.unionByName(manifest_df, allowMissingColumns=True).dropDuplicates(["variant"])
    manifest_df.coalesce(1).write.mode("overwrite").parquet(manifest_path)

    log(f"DONE  analysis_ready_assembly -> {analysis_path}")
    log(f"DONE  variant_manifest -> {manifest_path}")


# ============================================================================
# CLI
# ============================================================================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare VIDRA analysis-ready source tables using current OT release paths with AZ reintegration."
    )
    parser.add_argument("--genes", required=True, help="CSV file with Ensembl gene IDs")
    parser.add_argument("--colocalisation_threshold", required=True, type=float)
    parser.add_argument("--coloc_data_dir", required=True)
    parser.add_argument("--credible_set_dir", required=True)
    parser.add_argument("--study_data_dir", required=True)
    parser.add_argument("--variant_data_dir", required=True)
    parser.add_argument("--burden_evidence_dir", required=True)
    parser.add_argument("--clinvar_evidence_dir", required=True)
    parser.add_argument("--az_variants_file", required=True)
    parser.add_argument("--az_mapping_file", required=True)
    parser.add_argument("--az_gene_map_file", required=True)
    parser.add_argument("--az_burden_binary_dir", default="")
    parser.add_argument("--az_burden_quantitative_dir", default="")
    parser.add_argument("--target_data_dir", default="")
    parser.add_argument("--lead_variant_effect_dir", default="")
    parser.add_argument("--coloc_output", default="")
    parser.add_argument("--gwas_output", default="")
    parser.add_argument("--qtl_output", default="")
    parser.add_argument("--burden_output", default="")
    parser.add_argument("--coding_output", default=DEFAULT_CODING_OUTPUT)
    parser.add_argument("--clinvar_output", default=DEFAULT_CLINVAR_OUTPUT)
    parser.add_argument("--az_output", default="")
    parser.add_argument("--bucket_name", default="", help="GCS bucket name without gs://")
    parser.add_argument("--bucket_uri", default="", help="Bucket base URI/path override, e.g. gs://bucket or /data/bucket")
    parser.add_argument("--gcp_project", default="", help="Billing project used when Spark reads requester-pays GCS buckets")
    parser.add_argument("--spark_master", default="", help="Spark master for preparation jobs (for example local[4]); leave blank for ambient Spark such as Dataproc")
    parser.add_argument("--spark_driver_memory", default="", help="Spark driver memory for preparation jobs; leave blank to keep cluster defaults")
    parser.add_argument("--spark_shuffle_partitions", type=int, default=0, help="Spark SQL shuffle partitions for preparation jobs; 0 keeps cluster defaults")
    parser.add_argument("--spark_default_parallelism", type=int, default=0, help="Spark default parallelism for preparation jobs; 0 keeps cluster defaults")
    parser.add_argument("--analysis_ready_dir", default=DEFAULT_ANALYSIS_READY_DIR)
    parser.add_argument("--analysis_manifest_dir", default=DEFAULT_ANALYSIS_MANIFEST_DIR)
    parser.add_argument("--skip_analysis_ready", action="store_true", help="Skip final analysis-ready assembly")
    return parser.parse_args()


# ============================================================================
# Main orchestration
# ============================================================================


def main() -> None:
    require_pyspark()
    args = parse_args()
    validate_raw_az_burden_args(args)
    validate_step1_input_paths(args)
    log_step1_optional_input_modes(args)

    configure_spark_defaults(
        master=args.spark_master,
        driver_memory=args.spark_driver_memory,
        shuffle_partitions=args.spark_shuffle_partitions,
        default_parallelism=args.spark_default_parallelism,
        gcp_project=args.gcp_project,
        requester_pays_buckets=",".join(
            sorted(
                bucket
                for bucket in {
                    _gs_bucket_name(args.coloc_data_dir),
                    _gs_bucket_name(args.credible_set_dir),
                    _gs_bucket_name(args.study_data_dir),
                    _gs_bucket_name(args.variant_data_dir),
                    _gs_bucket_name(args.burden_evidence_dir),
                    _gs_bucket_name(args.clinvar_evidence_dir),
                    _gs_bucket_name(args.az_burden_binary_dir),
                    _gs_bucket_name(args.az_burden_quantitative_dir),
                    _gs_bucket_name(args.target_data_dir),
                    _gs_bucket_name(args.lead_variant_effect_dir),
                }
                if bucket
            )
        ),
    )

    spark = build_local_spark("vidra_prepare_analysis_input")
    cached_frames = []

    def _cache_frame(df, storage_level=None):
        cached_frames.append(df.persist(storage_level or StorageLevel.MEMORY_AND_DISK))
        return cached_frames[-1]

    try:
        genes = load_gene_ids_from_csv(args.genes, args.gcp_project)
        coloc_raw = spark.read.parquet(resolve_parquet_input(args.coloc_data_dir))
        study_raw = None
        credible_set_raw = None
        az_variant_index = None
        variant_projection = None
        lead_variant_effect_df = None
        # Optional source switch:
        # - with lead_variant_effect: do not load credible_set/study/full variant projection for
        #   production GWAS/QTL/coding assembly
        # - without it: keep the original release-based path unchanged
        if lead_variant_effect_inputs_enabled(args):
            az_variant_index = (
                spark.read.parquet(resolve_parquet_input(args.variant_data_dir))
                .select(F.col("variantId").cast("string").alias("ot_variant"))
                .dropna(subset=["ot_variant"])
                .dropDuplicates(["ot_variant"])
            )
            if not path_exists(args.lead_variant_effect_dir, args.gcp_project):
                raise FileNotFoundError(f"Lead variant effect data not found: {args.lead_variant_effect_dir}")
            lead_variant_effect_raw = spark.read.parquet(resolve_parquet_input(args.lead_variant_effect_dir))
            study_raw = _cache_frame(
                lead_variant_effect_raw.select(
                    F.col("studyId").cast("string").alias("studyId"),
                    optional_col("traitFromSourceMappedIds", lead_variant_effect_raw.schema, T.ArrayType(T.StringType())).alias("traitFromSourceMappedIds"),
                    optional_col("diseaseIds", lead_variant_effect_raw.schema, T.ArrayType(T.StringType())).alias("diseaseIds"),
                )
            )
            lead_variant_effect_df = _cache_frame(_normalise_lead_variant_effect_df(lead_variant_effect_raw))
        else:
            credible_set_raw = spark.read.parquet(resolve_parquet_input(args.credible_set_dir))
            study_raw = spark.read.parquet(resolve_parquet_input(args.study_data_dir))
            variant_projection = (
                spark.read.parquet(resolve_parquet_input(args.variant_data_dir)).select(
                    F.col("variantId").cast("string").alias("variantId"),
                    F.col("variantEffect"),
                    F.col("mostSevereConsequenceId").cast("string").alias("mostSevereConsequenceId"),
                    F.col("variantDescription").cast("string").alias("variantDescription"),
                )
            )
            az_variant_index = (
                variant_projection.select(F.col("variantId").alias("ot_variant"))
                .dropna(subset=["ot_variant"])
                .dropDuplicates(["ot_variant"])
            )
        burden_raw = spark.read.parquet(resolve_parquet_input(args.burden_evidence_dir))
        clinvar_raw = spark.read.parquet(resolve_parquet_input(args.clinvar_evidence_dir))

        log("START coloc_variants_per_gene")
        if lead_variant_effect_df is not None:
            coloc_df = _cache_frame(
                _build_coloc_output_df_from_lead(
                    spark,
                    genes,
                    args.colocalisation_threshold,
                    coloc_raw,
                    lead_variant_effect_df,
                )
            )
        else:
            coloc_df = _cache_frame(
                _build_coloc_output_df(
                    spark,
                    genes,
                    args.colocalisation_threshold,
                    coloc_raw,
                    credible_set_raw,
                    study_raw,
                )
            )
        _maybe_write_section_output("coloc_variants_per_gene", coloc_df, args.coloc_output)

        log("START gwas_variants")
        if lead_variant_effect_df is not None:
            gwas_df = _build_gwas_output_df_from_lead(coloc_df, lead_variant_effect_df)
        else:
            gwas_df = _build_gwas_output_df(coloc_df, credible_set_raw)
        if args.gwas_output:
            gwas_df = _cache_frame(gwas_df)
        _maybe_write_section_output("gwas_variants", gwas_df, args.gwas_output)

        log("START coding_gwas_non_coloc")
        if lead_variant_effect_df is not None:
            coding_df = _cache_frame(
                _build_coding_output_df_from_lead(
                    spark,
                    genes,
                    lead_variant_effect_df,
                    gwas_df,
                ),
                StorageLevel.DISK_ONLY,
            )
        else:
            coding_df = _cache_frame(
                _build_coding_output_df(
                    spark,
                    genes,
                    credible_set_raw,
                    variant_projection,
                    gwas_df,
                ),
                StorageLevel.DISK_ONLY,
            )
        _maybe_write_section_output("coding_gwas_non_coloc", coding_df, args.coding_output)

        log("START qtl_variants")
        if lead_variant_effect_df is not None:
            qtl_df = _build_qtl_output_df_from_lead(coloc_df, lead_variant_effect_df)
        else:
            qtl_df = _build_qtl_output_df(coloc_df, credible_set_raw)
        if args.qtl_output:
            qtl_df = _cache_frame(qtl_df)
        _maybe_write_section_output("qtl_variants", qtl_df, args.qtl_output)

        log("START burden_tests")
        burden_df = _build_burden_output_df(spark, genes, burden_raw)
        if raw_az_burden_inputs_enabled(args):
            log("START raw_az_burden_tests")
            raw_az_burden_df = _build_raw_az_burden_output_df(spark, args, genes)
            burden_df = burden_df.unionByName(raw_az_burden_df, allowMissingColumns=True)
            log("DONE  raw_az_burden_tests -> in_memory")
        if args.burden_output:
            burden_df = _cache_frame(burden_df)
        _maybe_write_section_output("burden_tests", burden_df, args.burden_output)

        log("START clinvar_variants")
        clinvar_df = _cache_frame(_build_clinvar_output_df(spark, genes, clinvar_raw))
        _maybe_write_section_output("clinvar_variants", clinvar_df, args.clinvar_output)

        log("START az_variants")
        # AZ rare variants always match against the release OT variant universe,
        # independent of whether lead_variant_effect is used for GWAS/QTL/coding.
        az_ready = _build_az_ready_dataset_from_inputs(
            spark,
            args,
            genes,
            az_variant_index,
        )
        if args.az_output:
            az_ready = _cache_frame(az_ready)
        _maybe_write_section_output("az_variants", az_ready, args.az_output)

        if not args.skip_analysis_ready:
            build_analysis_ready_dataset(
                spark,
                args,
                study_raw=study_raw,
                coloc=coloc_df,
                gwas=gwas_df,
                qtl=qtl_df,
                coding=coding_df,
                clinvar=clinvar_df,
                burden=burden_df,
                az_ready=az_ready,
            )

        log("Completed all ingestion sections.")
        persisted_outputs = [args.coding_output, args.clinvar_output]
        optional_outputs = [args.coloc_output, args.gwas_output, args.qtl_output, args.burden_output, args.az_output]
        for output in [*optional_outputs, *persisted_outputs]:
            if output:
                log(f"  - persisted: {output}")
        if not args.skip_analysis_ready:
            bucket_base = _resolve_bucket_base(args)
            if bucket_base:
                log(f"  - persisted: {_join_bucket_path(bucket_base, args.analysis_ready_dir)}")
                log(f"  - persisted: {_join_bucket_path(bucket_base, args.analysis_manifest_dir)}")
            else:
                log(f"  - persisted: {args.analysis_ready_dir}")
                log(f"  - persisted: {args.analysis_manifest_dir}")
    finally:
        for frame in reversed(cached_frames):
            try:
                frame.unpersist()
            except Exception:
                pass
        try:
            spark.stop()
        except Exception:
            pass


if __name__ == "__main__":
    main()
