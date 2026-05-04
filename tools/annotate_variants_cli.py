#!/usr/bin/env python3
"""VIDRA Step 2: annotate manifest variants with offline Ensembl VEP.

Annotation approach:
  1. Read the Step 1 `vidra_analysis_ready_manifest/` parquet dataset and
     collect distinct variant IDs in `CHR_POS_REF_ALT` form.
  2. Reuse existing annotations when available, then run one bulk Ensembl VEP
     pass offline against the remaining missing variants.
  3. Parse the bulk VEP JSON output into the compact annotation parquet
     consumed by `run_bayesian_analysis.py`.

This implementation keeps the permissive VIDRA_4 behavior around missing plugin
data, but preserves the VIDRA_2 encoding semantics for the annotation fields
that feed the Stan models.

Core plugin policy:
  - AlphaMissense, CADD, and REVEL are required because their scores are used
    directly by the current Stan models.
  - Conservation, PrimateAI, FoldX, and the remaining plugins stay optional
    and fall back to neutral defaults when absent.

Pipeline position:
  1. prepare_analysis_input.py  -> vidra_analysis_ready + manifest
  2. THIS SCRIPT               -> variant_annotations/annotations.parquet
  3. run_bayesian_analysis.py  -> joins analysis-ready rows with annotations

Transforms written to the annotation parquet:
  as_blosum62:
    sigmoid `1 / (1 + exp(-x))`; default 1.0 when VEP does not provide a score
  as_sift:
    raw SIFT score; missing values preserved for Step 3 imputation
  as_polyphen:
    raw PolyPhen score; missing values preserved
  as_cadd:
    `phred / 50` clamped to `[0, 1]`; missing values preserved
  as_alphamissense:
    raw AlphaMissense `am_pathogenicity`; missing values preserved
  as_revel:
    raw REVEL score; missing values preserved
  as_primateai:
    raw PrimateAI score; defaults to 0.0 when absent
  as_loftool:
    currently unused by the Stan model; always 0.0
  as_plddt:
    FoldX/ProtVar pLDDT scaled to `[0, 1]`; defaults to 0.0 when absent
  as_conservation:
    GERP RS clamped to the observed BigWig range, scaled to `[0, 1]`, then
    inverted as `1 - scaled`; missing values preserved
  as_consequence:
    ordinal code 0-12 ordered from least to most damaging consequence
  as_clinicalSignificance:
    17 ordered ClinVar-style categories scaled onto `[0, 1]`
  foldxDdq_raw:
    raw FoldX DDG value, left for downstream handling in Step 3
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.fs as pafs
import pyarrow.parquet as pq

# ============================================================================
# Constants and output schema
# ============================================================================
DEFAULT_ANALYSIS_MANIFEST_DIR = "vidra_analysis_ready_manifest"
DEFAULT_VEP_JSON_OUTPUT = "vep_annotations.json"
DEFAULT_ANNOTATIONS_OUTPUT = "variant_annotations/annotations.parquet"
DEFAULT_VEP_PARALLEL = 0
DEFAULT_VEP_FORKS = 0
DEFAULT_VEP_BUFFER_SIZE = 10000
MAX_AUTO_VEP_PARALLEL = 4

PLUGIN_ROOT = os.environ.get("VEP_PLUGIN_RESOURCES_DIR", "/plugin_resources")

VEP_PLUGIN_FILES = [
    {
        "module": "AlphaMissense",
        "plugin_arg": "AlphaMissense,file={path}",
        "local_relative_path": "alphamissense/AlphaMissense_hg38.tsv.gz",
        "candidate_paths": [
            "alphamissense/AlphaMissense_hg38.tsv.gz",
            "AlphaMissense_hg38.tsv.gz",
        ],
        "index_command": ["tabix", "-s", "1", "-b", "2", "-e", "2", "-f", "-S", "1", "{path}"],
    },
    {
        "module": "CADD",
        "plugin_arg": "CADD,snv={path}",
        "local_relative_path": "CADD/whole_genome_SNVs.tsv.gz",
        "candidate_paths": [
            "CADD/whole_genome_SNVs.tsv.gz",
        ],
    },
    {
        "module": "mutfunc",
        "plugin_arg": "mutfunc,db={path}",
        "local_relative_path": "mutfunc/mutfunc_data.db",
        "candidate_paths": [
            "mutfunc/mutfunc_data.db",
        ],
    },
    {
        "module": "EVE",
        "plugin_arg": "EVE,file={path}",
        "local_relative_path": "EVE/eve_plugin/eve_merged.vcf.gz",
        "candidate_paths": [
            "EVE/eve_plugin/eve_merged.vcf.gz",
        ],
    },
    {
        "module": "MaveDB",
        "plugin_arg": "MaveDB,file={path}",
        "local_relative_path": "maveDB/MaveDB_variants.tsv.gz",
        "candidate_paths": [
            "maveDB/MaveDB_variants.tsv.gz",
        ],
    },
    {
        "module": "LoFtool",
        "plugin_arg": "LoFtool,{path}",
        "local_relative_path": "loftool/LoFtool_scores.txt",
        "candidate_paths": [
            "loftool/LoFtool_scores.txt",
        ],
    },
    {
        "module": "PrimateAI",
        "plugin_arg": "PrimateAI,{path}",
        "local_relative_path": "primateAI/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz",
        "candidate_paths": [
            "primateAI/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz",
            "PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz",
        ],
        "index_command": ["tabix", "-s", "1", "-b", "2", "-e", "2", "-f", "{path}"],
    },
    {
        "module": "REVEL",
        "plugin_arg": "REVEL,file={path}",
        "local_relative_path": "revel/new_tabbed_revel_grch38.tsv.gz",
        "candidate_paths": [
            "revel/new_tabbed_revel_grch38.tsv.gz",
            "new_tabbed_revel_grch38.tsv.gz",
        ],
    },
    {
        "module": "pLI",
        "plugin_arg": "pLI,{path}",
        "local_relative_path": "pli/pLI_values.txt",
        "candidate_paths": [
            "pli/pLI_values.txt",
        ],
    },
    {
        "module": "Conservation",
        "plugin_arg": "Conservation,{path}",
        "local_relative_path": "conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
        "candidate_paths": [
            "conservation/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
            "gerp_conservation_scores.homo_sapiens.GRCh38.bw",
        ],
    },
]

OPTIONAL_PLUGIN_SIDECAR_SUFFIXES = (".tbi", ".csi")
PLUGIN_SIDE_CAR_REQUIRED = {"AlphaMissense", "CADD", "REVEL", "PrimateAI"}

# Global min/max from the GRCh38 GERP resource used to scale conservation to
# the `[0, 1]` range before inversion.
GERP_RS_MIN = -12.36
GERP_RS_MAX = 6.18

# Output parquet schema consumed directly by Step 3.
ANNOTATION_FIELDS = [
    ("variant", pa.string()),
    ("as_blosum62", pa.float64()),
    ("as_conservation", pa.float64()),
    ("as_sift", pa.float64()),
    ("as_polyphen", pa.float64()),
    ("as_cadd", pa.float64()),
    ("as_alphamissense", pa.float64()),
    ("as_revel", pa.float64()),
    ("as_primateai", pa.float64()),
    ("as_loftool", pa.float64()),
    ("as_plddt", pa.float64()),
    ("as_consequence", pa.int64()),
    ("as_clinicalSignificance", pa.float64()),
    ("foldxDdq_raw", pa.float64()),
    ("most_severe_consequence", pa.string()),
]
ANNOTATION_SCHEMA = pa.schema(ANNOTATION_FIELDS)

# Consequence ordinal encoding ordered from least to most damaging.
# This keeps the VIDRA_2 / Ensembl severity interpretation while using the
# current OT field names (`mostSevereConsequenceId` -> parsed VEP labels).
CONSEQUENCE_CATEGORIES = [
    "__unknown__",
    "non_coding_transcript_exon_variant",
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
    "synonymous_variant",
    "splice_donor_region_variant",
    "splice_region_variant",
    "missense_variant",
    "start_lost",
    "frameshift_variant",
    "stop_gained",
    "splice_donor_variant",
    "splice_acceptor_variant",
]
_CONSEQUENCE_LOOKUP = {name: idx for idx, name in enumerate(CONSEQUENCE_CATEGORIES)}

# Clinical significance encoding: 17 ordered categories scaled onto `[0, 1]`.
# Convention:
#   0.0  = benign / uninformative
#   1.0  = pathogenic
# Step 3 then uses this as a prior-like feature, matching VIDRA_2 semantics.
CLINICAL_SIG_CATEGORIES = [
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
_CLIN_SIG_LOOKUP = {
    cat: idx / float(len(CLINICAL_SIG_CATEGORIES) - 1)
    for idx, cat in enumerate(CLINICAL_SIG_CATEGORIES)
}
# Collapse underscore / compound labels from VEP or ClinVar-like payloads into
# the canonical 17-category set above.
_CLIN_SIG_NORMALISE = {
    "likely_pathogenic": "likely pathogenic",
    "likely_benign": "likely benign",
    "uncertain_significance": "uncertain significance",
    "not_provided": "not provided",
    "association_not_found": "association not found",
    "low_penetrance": "low penetrance",
    "confers_sensitivity": "confers sensitivity",
    "uncertain_risk_allele": "uncertain risk allele",
    "drug_response": "drug response",
    "likely_risk_allele": "likely risk allele",
    "risk_factor": "risk factor",
    "established_risk_allele": "established risk allele",
    "variant of uncertain significance": "uncertain significance",
    "benign/likely benign": "likely benign",
    "benign/likely_benign": "likely benign",
    "pathogenic/likely pathogenic": "likely pathogenic",
    "pathogenic/likely_pathogenic": "likely pathogenic",
    "conflicting interpretations of pathogenicity": "uncertain significance",
    "conflicting_interpretations_of_pathogenicity": "uncertain significance",
}


def log(message: str) -> None:
    print(f"[annotate_variants_cli] {message}", flush=True)


def log_progress(message: str) -> None:
    log(f"PROGRESS: {message}")


def parse_bool_flag(value) -> bool:
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "f", "no", "n", "off"}:
        return False
    raise ValueError(f"Could not interpret boolean flag value: {value}")


def parse_plugin_list(value) -> set[str]:
    plugins = set()
    for item in str(value or "").split(","):
        name = item.strip()
        if name:
            plugins.add(name)
    return plugins


def plugin_name(plugin_spec: str) -> str:
    return plugin_spec.split(",", 1)[0]


def parse_amino_acids(value) -> tuple[str | None, str | None]:
    if value is None:
        return None, None
    text = str(value).strip()
    if not text or "/" not in text:
        return None, None
    ref_aa, alt_aa = text.split("/", 1)
    ref_aa = ref_aa.strip() or None
    alt_aa = alt_aa.strip() or None
    return ref_aa, alt_aa


def maybe_first(value):
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


# ============================================================================
# Path utilities
# ============================================================================
def is_gs_path(path: str | Path | None) -> bool:
    return str(path or "").startswith("gs://")


def _strip_gs_prefix(path: str) -> str:
    return str(path)[len("gs://"):]


def _get_gcsfs():
    try:
        import gcsfs
    except ImportError as exc:
        raise RuntimeError(
            "Working with gs:// paths requires gcsfs."
        ) from exc
    return gcsfs.GCSFileSystem()


def path_exists(path: str | Path) -> bool:
    if is_gs_path(path):
        fs = _get_gcsfs()
        return fs.exists(_strip_gs_prefix(str(path)))
    return Path(path).exists()


def copy_gs_file_to_local(source: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    fs = _get_gcsfs()
    fs.get(_strip_gs_prefix(source), str(destination))


def copy_file_to_local(source: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if is_gs_path(source):
        copy_gs_file_to_local(source, destination)
        return
    shutil.copy2(source, destination)


def copy_optional_sidecars_to_local(source: str, destination: Path) -> None:
    for suffix in OPTIONAL_PLUGIN_SIDECAR_SUFFIXES:
        sidecar_source = f"{source}{suffix}"
        if path_exists(sidecar_source):
            copy_file_to_local(sidecar_source, Path(f"{destination}{suffix}"))


def has_local_plugin_sidecar(path: Path) -> bool:
    return any(Path(f"{path}{suffix}").exists() for suffix in OPTIONAL_PLUGIN_SIDECAR_SUFFIXES)


def ensure_plugin_index(local_path: Path, plugin_spec: dict) -> None:
    if has_local_plugin_sidecar(local_path):
        return
    index_command = plugin_spec.get("index_command")
    if not index_command:
        return

    cmd = [part.format(path=str(local_path)) for part in index_command]
    log(f"$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def open_parquet_dataset(path_like: str):
    if is_gs_path(path_like):
        fs = _get_gcsfs()
        arrow_fs = pafs.PyFileSystem(pafs.FSSpecHandler(fs))
        return ds.dataset(_strip_gs_prefix(path_like), filesystem=arrow_fs, format="parquet")
    return ds.dataset(path_like, format="parquet")


def join_root_path(root: str, relative_path: str) -> str:
    base = str(root).rstrip("/")
    rel = relative_path.lstrip("/")
    return f"{base}/{rel}" if base else rel


def _join_bucket_path(base_uri: str, name: str) -> str:
    return f"{base_uri.rstrip('/')}/{name.lstrip('/')}"


def _resolve_bucket_base(args) -> str:
    if args.bucket_uri:
        return args.bucket_uri.rstrip("/")
    if args.bucket_name:
        return f"gs://{args.bucket_name}"
    return ""


def _write_table_to_path(table: pa.Table, output_path: str) -> None:
    if is_gs_path(output_path):
        fs = _get_gcsfs()
        path_no_scheme = _strip_gs_prefix(output_path)
        with fs.open(path_no_scheme, "wb") as handle:
            pq.write_table(table, handle)
        return

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        if out.is_dir():
            shutil.rmtree(out)
        else:
            out.unlink()
    pq.write_table(table, out)


# ============================================================================
# Variant preprocessing
# ============================================================================
def parse_variant_id(variant_id: str):
    parts = str(variant_id).strip().split("_")
    if len(parts) != 4:
        return None

    chrom, pos, ref, alt = parts
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if not pos:
        return None
    return chrom, pos, ref, alt


def chrom_sort_key(chrom: str):
    c = str(chrom).upper().replace("CHR", "")
    if c.isdigit():
        return (0, int(c))
    if c == "X":
        return (1, 23)
    if c == "Y":
        return (1, 24)
    if c in {"MT", "M"}:
        return (1, 25)
    return (2, c)


def pos_sort_key(pos: str):
    p = str(pos)
    if p.isdigit():
        return (0, int(p))
    return (1, p)


def tuple_to_variant_id(chrom: str, pos: str, ref: str, alt: str) -> str:
    return f"{chrom}_{pos}_{ref}_{alt}"


def variant_sort_key(variant_id: str):
    parsed = parse_variant_id(variant_id)
    if parsed is None:
        return (3, str(variant_id))
    chrom, pos, ref, alt = parsed
    return (chrom_sort_key(chrom), pos_sort_key(pos), ref, alt)


def dedupe_and_sort_variant_ids(variant_ids) -> list[str]:
    return sorted({str(variant_id).strip() for variant_id in variant_ids if variant_id}, key=variant_sort_key)


def variant_rows_from_ids(variant_ids) -> list[tuple[str, str, str, str]]:
    rows = []
    invalid_variant_count = 0
    invalid_variant_ids = []
    for variant_id in dedupe_and_sort_variant_ids(variant_ids):
        parsed = parse_variant_id(variant_id)
        if parsed is None:
            invalid_variant_count += 1
            if len(invalid_variant_ids) < 10:
                invalid_variant_ids.append(variant_id)
            continue
        rows.append(parsed)

    if invalid_variant_ids:
        log(
            f"WARN: skipped {invalid_variant_count} non-parseable variant value(s) "
            f"while preparing VEP input; sample={invalid_variant_ids[:5]}"
        )

    return rows


def resolve_analysis_manifest_path(args) -> str:
    bucket_base = _resolve_bucket_base(args)
    if bucket_base:
        return _join_bucket_path(bucket_base, args.analysis_manifest_dir)
    return args.analysis_manifest_dir


def load_manifest_variant_ids(manifest_path: str, *, allow_empty: bool = False) -> list[str]:
    dataset = open_parquet_dataset(manifest_path)
    if "variant" not in dataset.schema.names:
        raise ValueError(f"Manifest parquet must contain 'variant': {manifest_path}")

    variant_ids: set[str] = set()
    total_rows = 0
    null_variants = 0
    invalid_variant_ids = 0
    invalid_samples = []

    scanner = dataset.scanner(columns=["variant"])
    for batch in scanner.to_batches():
        for value in batch.column(0).to_pylist():
            total_rows += 1
            if value is None:
                null_variants += 1
                continue
            variant_id = str(value).strip()
            if parse_variant_id(variant_id) is None:
                invalid_variant_ids += 1
                if len(invalid_samples) < 10:
                    invalid_samples.append(variant_id)
                continue
            variant_ids.add(variant_id)

    if not variant_ids:
        if allow_empty:
            return []
        raise ValueError(
            f"No valid variants extracted from {manifest_path}. "
            f"rows={total_rows}, null_variant={null_variants}, "
            f"invalid_variant={invalid_variant_ids}, "
            f"sample_invalid={invalid_samples}"
        )

    if invalid_variant_ids:
        log(
            f"WARN: skipped {invalid_variant_ids} non-parseable manifest variant value(s) "
            f"in {manifest_path}; sample={invalid_samples[:5]}"
        )

    return dedupe_and_sort_variant_ids(variant_ids)


def write_vep_input_file(variant_ids, output_path: Path) -> None:
    deduped = variant_rows_from_ids(variant_ids)
    with output_path.open("w", encoding="utf-8") as handle:
        for chrom, pos, ref, alt in deduped:
            handle.write(f"{chrom} {pos} . {ref} {alt}\n")


# ============================================================================
# VEP execution
# ============================================================================
def resolve_cpu_count(cpu_count: int | None = None) -> int:
    detected = cpu_count if cpu_count is not None else os.cpu_count()
    try:
        value = int(detected)
    except (TypeError, ValueError):
        value = 1
    return max(1, value)


def resolve_vep_parallelism(
    variant_count: int,
    requested_parallel: int,
    requested_forks: int,
    *,
    cpu_count: int | None = None,
) -> tuple[int, int, int]:
    resolved_cpu_count = resolve_cpu_count(cpu_count)
    auto_parallel = min(MAX_AUTO_VEP_PARALLEL, max(1, resolved_cpu_count // 2))
    resolved_parallel = auto_parallel if requested_parallel == 0 else max(1, int(requested_parallel))
    if variant_count > 0:
        resolved_parallel = min(resolved_parallel, variant_count)
    else:
        resolved_parallel = 1

    resolved_forks = (
        max(1, resolved_cpu_count // resolved_parallel)
        if requested_forks == 0
        else max(1, int(requested_forks))
    )
    return resolved_parallel, resolved_forks, resolved_cpu_count


def plan_vep_chunks(variant_ids, resolved_parallel: int):
    ordered_variant_ids = dedupe_and_sort_variant_ids(variant_ids)
    if not ordered_variant_ids:
        return []

    chunk_count = max(1, min(int(resolved_parallel), len(ordered_variant_ids)))
    chunk_size = (len(ordered_variant_ids) + chunk_count - 1) // chunk_count
    return [
        ordered_variant_ids[start:start + chunk_size]
        for start in range(0, len(ordered_variant_ids), chunk_size)
    ]


def resolve_existing_resource(root: str, candidate_paths: list[str]) -> str | None:
    for relative_path in candidate_paths:
        candidate = join_root_path(root, relative_path)
        if path_exists(candidate):
            return candidate
    return None


def prepare_runtime_resources(args, runtime_root: Path) -> argparse.Namespace:
    runtime_args = argparse.Namespace(**vars(args))

    if args.vep_cache_archive:
        source = str(args.vep_cache_archive)
        if not path_exists(source):
            raise FileNotFoundError(f"Configured VEP cache archive was not found: {source}")
        cache_archive = runtime_root / Path(source).name
        cache_root = runtime_root / "vep_cache"
        copy_file_to_local(source, cache_archive)
        cache_root.mkdir(parents=True, exist_ok=True)
        with tarfile.open(cache_archive, "r:*") as archive:
            archive.extractall(path=cache_root)
        runtime_args.vep_dir_cache = str(cache_root)
    elif is_gs_path(args.vep_dir_cache):
        raise RuntimeError(
            "gs:// cache directories are not used directly for VEP. "
            "Provide --vep_cache_archive so the cache can be localized into the container."
        )

    if is_gs_path(args.vep_plugin_data_dir):
        localized_root = runtime_root / "plugin_resources"
        for plugin_spec in VEP_PLUGIN_FILES:
            source = resolve_existing_resource(args.vep_plugin_data_dir, plugin_spec["candidate_paths"])
            if source is None:
                continue
            target = localized_root / plugin_spec["local_relative_path"]
            copy_file_to_local(source, target)
            copy_optional_sidecars_to_local(source, target)
            ensure_plugin_index(target, plugin_spec)
        runtime_args.vep_plugin_data_dir = str(localized_root)

    if args.foldx_file and is_gs_path(args.foldx_file):
        localized_foldx = runtime_root / "foldx" / Path(args.foldx_file).name
        copy_gs_file_to_local(args.foldx_file, localized_foldx)
        runtime_args.foldx_file = str(localized_foldx)

    return runtime_args


def plugin_statuses(plugin_data_root: str, plugins_dir: str, required_plugins: set[str]) -> list[dict]:
    statuses: list[dict] = []
    for plugin_spec in VEP_PLUGIN_FILES:
        module_name = plugin_spec["module"]
        local_path = Path(plugin_data_root) / plugin_spec["local_relative_path"]
        module_path = Path(plugins_dir) / f"{module_name}.pm"
        plugin = plugin_spec["plugin_arg"].format(path=str(local_path))
        requires_sidecar = module_name in PLUGIN_SIDE_CAR_REQUIRED
        sidecar_ok = True
        if requires_sidecar:
            sidecar_ok = has_local_plugin_sidecar(local_path)

        statuses.append(
            {
                "name": module_name,
                "plugin_arg": plugin,
                "required": module_name in required_plugins,
                "module_ok": module_path.is_file(),
                "data_ok": local_path.is_file(),
                "sidecar_ok": sidecar_ok,
                "requires_sidecar": requires_sidecar,
            }
        )
    return statuses


def log_plugin_statuses(statuses: list[dict]) -> None:
    for status in statuses:
        log(
            "PLUGIN "
            f"{status['name']}: "
            f"required={status['required']} "
            f"module={status['module_ok']} "
            f"data={status['data_ok']} "
            f"index={status['sidecar_ok']}"
        )


def build_plugin_args(plugin_data_root: str, plugins_dir: str, required_plugins: set[str]) -> tuple[list[str], list[dict]]:
    plugin_args: list[str] = []
    statuses = plugin_statuses(plugin_data_root, plugins_dir, required_plugins)
    missing_required = []

    for status in statuses:
        plugin = status["plugin_arg"]
        ready = status["module_ok"] and status["data_ok"] and status["sidecar_ok"]

        if ready:
            plugin_args.extend(["--plugin", plugin])
        else:
            level = "ERROR" if status["required"] else "WARN"
            log(
                f"{level}: skipping plugin '{plugin}' "
                f"(module={status['module_ok']}, data={status['data_ok']}, index={status['sidecar_ok']})"
            )
            if status["required"]:
                missing_required.append(status["name"])

    log_plugin_statuses(statuses)
    if missing_required:
        raise RuntimeError(
            "Required VEP plugins are unavailable: "
            + ", ".join(sorted(missing_required))
        )

    plugin_args.extend(["--plugin", "Blosum62"])
    return plugin_args, statuses


def build_vep_command(
    args,
    input_tmp: Path,
    output_json: Path,
    *,
    vep_forks: int,
    vep_buffer_size: int,
) -> list[str]:
    cache_dir = Path(str(args.vep_dir_cache))
    if not (cache_dir / "homo_sapiens").exists():
        raise RuntimeError(
            f"Missing VEP cache directory {(cache_dir / 'homo_sapiens')}. "
            "Mount a local cache root or provide --vep_cache_archive."
        )

    vep_bin = "/opt/vep/src/ensembl-vep/vep"
    if not Path(vep_bin).exists():
        vep_bin = "vep"
    plugin_args = getattr(args, "vep_plugin_args", None)
    if plugin_args is None:
        plugin_args, _ = build_plugin_args(
            args.vep_plugin_data_dir,
            args.vep_plugins_dir,
            parse_plugin_list(getattr(args, "required_vep_plugins", "")),
        )
    cmd = [
        vep_bin,
        "--cache", "--offline", "--species", "homo_sapiens", "--assembly", "GRCh38",
        "--dir_cache", args.vep_dir_cache,
        "--dir_plugins", args.vep_plugins_dir,
        "--json", "--canonical",
        "--sift", "b",
        "--polyphen", "b",
        "--force_overwrite", "--no_stats",
        "--check_existing", "--uniprot",
        "--fork", str(vep_forks),
        "--buffer_size", str(vep_buffer_size),
        "-i", str(input_tmp),
        "-o", str(output_json),
        *plugin_args,
    ]
    return cmd


def run_vep(
    args,
    input_tmp: Path,
    output_json: Path,
    *,
    vep_forks: int,
    vep_buffer_size: int,
) -> None:
    cmd = build_vep_command(
        args,
        input_tmp,
        output_json,
        vep_forks=vep_forks,
        vep_buffer_size=vep_buffer_size,
    )

    log("$ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_vep_chunk(
    args,
    section_name: str,
    chunk_index: int,
    total_chunks: int,
    variant_ids,
    input_tmp: Path,
    output_json: Path,
    *,
    vep_forks: int,
    vep_buffer_size: int,
) -> Path:
    chunk_variant_count = len(variant_ids)
    log_progress(
        f"{section_name}: chunk {chunk_index + 1}/{total_chunks} starting "
        f"({chunk_variant_count:,} variant(s))"
    )
    write_vep_input_file(variant_ids, input_tmp)
    run_vep(
        args,
        input_tmp,
        output_json,
        vep_forks=vep_forks,
        vep_buffer_size=vep_buffer_size,
    )
    if not output_json.exists():
        raise RuntimeError(
            f"{section_name} chunk {chunk_index + 1}/{total_chunks} "
            "completed but did not produce a JSON output file"
        )
    log_progress(
        f"{section_name}: chunk {chunk_index + 1}/{total_chunks} completed "
        f"({chunk_variant_count:,} variant(s))"
    )
    return output_json


def merge_vep_chunk_outputs(chunk_output_paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    wrote_any = False
    previous_chunk_ended_with_newline = True

    with output_path.open("wb") as destination:
        for chunk_output_path in chunk_output_paths:
            if not chunk_output_path.exists():
                raise FileNotFoundError(f"VEP chunk output not found: {chunk_output_path}")

            with chunk_output_path.open("rb") as source:
                first_chunk = source.read(1024 * 1024)
                if not first_chunk:
                    continue

                if wrote_any and not previous_chunk_ended_with_newline:
                    destination.write(b"\n")

                destination.write(first_chunk)
                previous_chunk_ended_with_newline = first_chunk.endswith(b"\n")

                while True:
                    payload = source.read(1024 * 1024)
                    if not payload:
                        break
                    destination.write(payload)
                    previous_chunk_ended_with_newline = payload.endswith(b"\n")

                wrote_any = True


# ============================================================================
# Bulk VEP runner
# ============================================================================
def write_placeholder_json(output_json: str) -> None:
    output_path = Path(output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        if output_path.is_dir():
            raise RuntimeError(f"Expected file output but found directory: {output_path}")
        output_path.unlink()
    output_path.write_text("", encoding="utf-8")


def annotate_variant_ids(args, section_name: str, variant_ids, output_json: str):
    output_path = Path(output_json)
    if output_path.exists():
        if output_path.is_dir():
            raise RuntimeError(f"Expected file output but found directory: {output_path}")
        output_path.unlink()

    ordered_variant_ids = dedupe_and_sort_variant_ids(variant_ids)
    distinct_variant_count = len(ordered_variant_ids)
    if not distinct_variant_count:
        write_placeholder_json(output_json)
        log(f"SKIP  {section_name} -> {output_path} (0 missing variants)")
        log_progress(f"{section_name}: skipped because no missing variants required annotation")
        return []

    resolved_parallel, resolved_forks, resolved_cpu_count = resolve_vep_parallelism(
        distinct_variant_count,
        args.vep_parallel,
        args.vep_forks,
    )
    chunk_variant_ids = plan_vep_chunks(ordered_variant_ids, resolved_parallel)
    chunk_sizes = ", ".join(f"{len(chunk):,}" for chunk in chunk_variant_ids)

    log(f"START {section_name}")
    log_progress(
        f"{section_name}: starting VEP for {distinct_variant_count:,} distinct variant(s); "
        f"chunks={len(chunk_variant_ids)}, forks_per_chunk={resolved_forks}, "
        f"buffer_size={args.vep_buffer_size}, cpu_count={resolved_cpu_count}"
    )
    log_progress(
        f"{section_name}: chunk sizes [{chunk_sizes}]"
    )

    with tempfile.TemporaryDirectory(prefix=f"vidra_{section_name}_") as tmpdir:
        work = Path(tmpdir)
        chunk_specs = [
            (
                chunk_index,
                chunk_ids,
                work / f"input.chunk{chunk_index:04d}.tmp",
                work / f"vars-annotated.chunk{chunk_index:04d}.json",
            )
            for chunk_index, chunk_ids in enumerate(chunk_variant_ids)
        ]

        if len(chunk_specs) == 1:
            _, chunk_ids, input_tmp, json_tmp = chunk_specs[0]
            run_vep_chunk(
                args,
                section_name,
                0,
                1,
                chunk_ids,
                input_tmp,
                json_tmp,
                vep_forks=resolved_forks,
                vep_buffer_size=args.vep_buffer_size,
            )
        else:
            completed_outputs: dict[int, Path] = {}
            failures: list[tuple[int, Exception]] = []
            with ThreadPoolExecutor(max_workers=len(chunk_specs)) as executor:
                futures = {
                    executor.submit(
                        run_vep_chunk,
                        args,
                        section_name,
                        chunk_index,
                        len(chunk_specs),
                        chunk_ids,
                        input_tmp,
                        json_tmp,
                        vep_forks=resolved_forks,
                        vep_buffer_size=args.vep_buffer_size,
                    ): chunk_index
                    for chunk_index, chunk_ids, input_tmp, json_tmp in chunk_specs
                }
                for future in as_completed(futures):
                    chunk_index = futures[future]
                    try:
                        completed_outputs[chunk_index] = future.result()
                    except Exception as exc:
                        failures.append((chunk_index, exc))
                        log(
                            f"ERROR: {section_name} chunk {chunk_index + 1}/{len(chunk_specs)} failed: {exc}"
                        )

            if failures:
                first_chunk_index, first_exception = sorted(failures, key=lambda item: item[0])[0]
                raise RuntimeError(
                    f"{section_name} failed during VEP chunk {first_chunk_index + 1}/{len(chunk_specs)}"
                ) from first_exception

        merge_vep_chunk_outputs(
            [chunk_output_path for _, _, _, chunk_output_path in chunk_specs],
            output_path,
        )

    log(f"DONE  {section_name} -> {output_path}")
    log_progress(
        f"{section_name}: completed VEP for {distinct_variant_count:,} distinct variant(s)"
    )
    return ordered_variant_ids


# ============================================================================
# VEP JSON parsing + transforms
# ============================================================================
def _safe_float(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(x) or math.isinf(x):
        return None
    return x


def _normalise_clin_sig(value: str | None) -> str:
    if not value:
        return "not provided"
    cleaned = value.lower().strip()
    return _CLIN_SIG_NORMALISE.get(cleaned, cleaned)


def _variant_from_vep_item(item: dict):
    candidate = item.get("id")
    parsed = parse_variant_id(str(candidate)) if candidate is not None else None
    if parsed is not None:
        return tuple_to_variant_id(*parsed)

    inp = item.get("input", "")
    parts = str(inp).replace("\t", " ").split()
    if len(parts) >= 5:
        chrom = parts[0]
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        pos = parts[1]
        ref = parts[3]
        alt = parts[4]
        return tuple_to_variant_id(chrom, pos, ref, alt)

    return None


def _best_transcript(item: dict):
    transcripts = item.get("transcript_consequences", [])
    if not isinstance(transcripts, list) or not transcripts:
        return None

    canonical = [tc for tc in transcripts if isinstance(tc, dict) and tc.get("canonical") == 1]
    if not canonical:
        canonical = [tc for tc in transcripts if isinstance(tc, dict)]
    if not canonical:
        return None

    return max(canonical, key=lambda tc: sum(1 for value in tc.values() if value is not None))


def _get_field(mapping: dict, *keys):
    for key in keys:
        if key in mapping:
            return mapping[key]
    lower_map = {k.lower(): v for k, v in mapping.items()}
    for key in keys:
        value = lower_map.get(key.lower())
        if value is not None:
            return value
    return None


def _get_nested_field(mapping: dict, *paths):
    """Get nested field from a mapping using case-insensitive keys."""
    for path in paths:
        cur = mapping
        ok = True
        for key in path:
            if not isinstance(cur, dict):
                ok = False
                break
            cur = _get_field(cur, key)
            if cur is None:
                ok = False
                break
        if ok and cur is not None:
            return cur
    return None


def _first_not_none(*values):
    for value in values:
        if value is not None:
            return value
    return None


def _pick_best_record(current: dict | None, candidate: dict) -> dict:
    if current is None:
        return candidate

    keys = [
        "blosum62", "conservation", "sift_score", "polyphen_score", "cadd_phred",
        "am_pathogenicity", "revel", "primateai_score", "clinicalSignificance",
        "most_severe_consequence", "swissprot", "protein_start", "ref_aa", "alt_aa",
    ]
    cur_score = sum(1 for key in keys if current.get(key) is not None)
    cand_score = sum(1 for key in keys if candidate.get(key) is not None)
    return candidate if cand_score > cur_score else current


def parse_vep_json_outputs(json_paths: list[str]) -> dict[str, dict]:
    by_variant: dict[str, dict] = {}

    for json_path in json_paths:
        path = Path(json_path)
        if not path.exists():
            raise FileNotFoundError(f"VEP JSON output not found: {path}")

        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue

                try:
                    item = json.loads(line)
                except json.JSONDecodeError:
                    continue

                variant = _variant_from_vep_item(item)
                if not variant:
                    continue

                rec = {
                    "variant": variant,
                    "most_severe_consequence": item.get("most_severe_consequence"),
                    "blosum62": None,
                    "conservation": None,
                    "sift_score": None,
                    "polyphen_score": None,
                    "cadd_phred": None,
                    "am_pathogenicity": None,
                    "revel": None,
                    "primateai_score": None,
                    "clinicalSignificance": None,
                    "swissprot": None,
                    "protein_start": None,
                    "ref_aa": None,
                    "alt_aa": None,
                }

                tc = _best_transcript(item)
                if tc:
                    rec["blosum62"] = _first_not_none(
                        _get_field(tc, "blosum62", "Blosum62"),
                        _get_nested_field(tc, ("blosum62", "score"), ("Blosum62", "score")),
                    )
                    rec["conservation"] = _first_not_none(
                        _get_field(tc, "conservation", "Conservation", "gerp", "GERP"),
                        _get_nested_field(
                            tc,
                            ("conservation", "score"),
                            ("Conservation", "score"),
                            ("conservation", "gerp"),
                            ("Conservation", "gerp"),
                            ("conservation", "gerp_rs"),
                            ("Conservation", "gerp_rs"),
                        ),
                    )
                    rec["sift_score"] = _first_not_none(
                        _get_field(tc, "sift_score", "SIFT_score"),
                        _get_nested_field(tc, ("sift", "score"), ("SIFT", "score")),
                    )
                    rec["polyphen_score"] = _first_not_none(
                        _get_field(tc, "polyphen_score", "PolyPhen_score"),
                        _get_nested_field(tc, ("polyphen", "score"), ("PolyPhen", "score")),
                    )
                    rec["cadd_phred"] = _first_not_none(
                        _get_field(tc, "cadd_phred", "CADD_phred", "CADD_PHRED"),
                        _get_nested_field(
                            tc,
                            ("cadd", "phred"),
                            ("CADD", "phred"),
                            ("cadd", "cadd_phred"),
                            ("CADD", "CADD_phred"),
                            ("cadd", "score"),
                            ("CADD", "score"),
                        ),
                    )
                    rec["am_pathogenicity"] = _first_not_none(
                        _get_field(tc, "am_pathogenicity", "AM_pathogenicity"),
                        _get_nested_field(
                            tc,
                            ("alphamissense", "am_pathogenicity"),
                            ("AlphaMissense", "am_pathogenicity"),
                            ("alphamissense", "pathogenicity"),
                            ("AlphaMissense", "pathogenicity"),
                        ),
                    )
                    rec["revel"] = _first_not_none(
                        _get_field(tc, "revel", "revel_score", "REVEL"),
                        _get_nested_field(
                            tc,
                            ("revel", "score"),
                            ("REVEL", "score"),
                            ("revel", "revel_score"),
                            ("REVEL", "revel_score"),
                        ),
                    )
                    rec["primateai_score"] = _first_not_none(
                        _get_field(
                            tc,
                            "primateai",
                            "primateai_score",
                            "PrimateAI_score",
                            "primatai_score",
                            "primatai_pred",
                        ),
                        _get_nested_field(
                            tc,
                            ("primateai", "score"),
                            ("PrimateAI", "score"),
                            ("primateai", "primateai_score"),
                            ("PrimateAI", "PrimateAI_score"),
                            ("primateai", "pred"),
                            ("PrimateAI", "pred"),
                        ),
                    )
                    rec["swissprot"] = maybe_first(
                        _first_not_none(
                            _get_field(tc, "swissprot", "SWISSPROT"),
                            _get_nested_field(tc, ("swissprot",), ("SWISSPROT",)),
                        )
                    )
                    rec["protein_start"] = _first_not_none(
                        _get_field(tc, "protein_start", "proteinStart"),
                        _get_nested_field(tc, ("protein_start",), ("proteinStart",)),
                    )
                    rec["ref_aa"], rec["alt_aa"] = parse_amino_acids(
                        _first_not_none(
                            _get_field(tc, "amino_acids", "aminoAcids"),
                            _get_nested_field(tc, ("amino_acids",), ("aminoAcids",)),
                        )
                    )

                # ClinVar significance from colocated variants.
                best_sig = None
                severity = {
                    "pathogenic": 4,
                    "likely_pathogenic": 3,
                    "pathogenic/likely_pathogenic": 3,
                    "benign": 2,
                    "likely_benign": 2,
                    "benign/likely_benign": 2,
                }
                for coloc in item.get("colocated_variants", []) or []:
                    clin_sig = coloc.get("clin_sig")
                    if not clin_sig:
                        continue
                    if isinstance(clin_sig, (list, tuple, set)):
                        clin_sig_iter = [str(x) for x in clin_sig]
                    else:
                        clin_sig_iter = str(clin_sig).split(",")
                    for sig in clin_sig_iter:
                        sig_norm = _normalise_clin_sig(sig)
                        if best_sig is None or severity.get(sig_norm, 0) > severity.get(best_sig, 0):
                            best_sig = sig_norm
                if best_sig:
                    rec["clinicalSignificance"] = best_sig

                if rec["conservation"] is None:
                    rec["conservation"] = _first_not_none(
                        _get_field(item, "conservation", "Conservation", "gerp", "GERP"),
                        _get_nested_field(
                            item,
                            ("conservation", "score"),
                            ("Conservation", "score"),
                            ("conservation", "gerp"),
                            ("Conservation", "gerp"),
                            ("conservation", "gerp_rs"),
                            ("Conservation", "gerp_rs"),
                        ),
                    )

                by_variant[variant] = _pick_best_record(by_variant.get(variant), rec)

    return by_variant


def load_foldx_lookup(foldx_path: str, needed_keys: set[tuple[str, int, str, str]]) -> dict[tuple[str, int, str, str], tuple[float, float]]:
    lookup: dict[tuple[str, int, str, str], tuple[float, float]] = {}
    if not foldx_path or not Path(foldx_path).is_file():
        return lookup

    opener = gzip.open if str(foldx_path).endswith(".gz") else open
    with opener(foldx_path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            if len(row) < 8:
                continue
            try:
                key = (row[0], int(row[1]), row[4], row[5])
                if key not in needed_keys:
                    continue
                lookup[key] = (float(row[6]), float(row[7]))
            except (TypeError, ValueError):
                continue
    return lookup


def collect_needed_foldx_keys(raw_records: dict[str, dict]) -> set[tuple[str, int, str, str]]:
    keys: set[tuple[str, int, str, str]] = set()
    for record in raw_records.values():
        accession = record.get("swissprot")
        protein_start = record.get("protein_start")
        ref_aa = record.get("ref_aa")
        alt_aa = record.get("alt_aa")
        if not accession or not protein_start or not ref_aa or not alt_aa:
            continue
        try:
            keys.add((str(accession).split(".", 1)[0], int(float(protein_start)), str(ref_aa), str(alt_aa)))
        except (TypeError, ValueError):
                continue
    return keys


def count_foldx_eligible_variants(raw_records: dict[str, dict]) -> int:
    count = 0
    for record in raw_records.values():
        if (
            record.get("swissprot")
            and record.get("protein_start")
            and record.get("ref_aa")
            and record.get("alt_aa")
        ):
            count += 1
    return count


def _transform_record(raw: dict | None, variant: str, foldx_lookup: dict[tuple[str, int, str, str], tuple[float, float]] | None = None) -> dict:
    """Map one VEP/FoldX record onto the compact Step 2 annotation schema.

    Missing-value policy intentionally mirrors the earlier VIDRA pipeline:
    missense-specific scores are usually left as null so Step 3 can do per-gene
    imputation, while a few features keep explicit neutral defaults.
    """
    raw = raw or {}

    blosum = _safe_float(raw.get("blosum62"))
    sift = _safe_float(raw.get("sift_score"))
    polyphen = _safe_float(raw.get("polyphen_score"))
    cadd = _safe_float(raw.get("cadd_phred"))
    alphamissense = _safe_float(raw.get("am_pathogenicity"))
    revel = _safe_float(raw.get("revel"))
    primateai = _safe_float(raw.get("primateai_score"))
    conservation = _safe_float(raw.get("conservation"))

    if blosum is None:
        as_blosum62 = 1.0
    else:
        as_blosum62 = 1.0 / (1.0 + math.exp(-blosum))

    # GERP RS: clamp to the resource-wide range, scale to [0,1], then invert so
    # larger values track the "more damaging" direction used elsewhere.
    if conservation is None:
        as_conservation = None
    else:
        scaled = (conservation - GERP_RS_MIN) / (GERP_RS_MAX - GERP_RS_MIN)
        scaled = max(0.0, min(1.0, scaled))
        as_conservation = 1.0 - scaled

    most_severe = raw.get("most_severe_consequence") or "__unknown__"
    as_consequence = _CONSEQUENCE_LOOKUP.get(most_severe, 0)

    clin_sig = _normalise_clin_sig(raw.get("clinicalSignificance"))
    as_clinical = _CLIN_SIG_LOOKUP.get(clin_sig, 0.0)

    foldx_ddg = None
    plddt = None
    if foldx_lookup:
        accession = raw.get("swissprot")
        protein_start = raw.get("protein_start")
        ref_aa = raw.get("ref_aa")
        alt_aa = raw.get("alt_aa")
        if accession and protein_start and ref_aa and alt_aa:
            try:
                key = (str(accession).split(".", 1)[0], int(float(protein_start)), str(ref_aa), str(alt_aa))
            except (TypeError, ValueError):
                key = None
            if key is not None:
                foldx_match = foldx_lookup.get(key)
                if foldx_match:
                    foldx_ddg, plddt = foldx_match

    # Keep this output close to VIDRA_2:
    # - raw missense scores are preserved where available
    # - neutral defaults are used only for fields that previously behaved that way
    # - Step 3 remains responsible for any downstream inversion / imputation
    return {
        "variant": variant,
        "as_blosum62": float(as_blosum62),
        "as_conservation": as_conservation,
        "as_sift": sift,
        "as_polyphen": polyphen,
        "as_cadd": None if cadd is None else max(0.0, min(1.0, cadd / 50.0)),
        "as_alphamissense": alphamissense,
        "as_revel": revel,
        "as_primateai": 0.0 if primateai is None else primateai,
        "as_loftool": 0.0,
        "as_plddt": 0.0 if plddt is None else max(0.0, min(1.0, float(plddt) / 100.0)),
        "as_consequence": int(as_consequence),
        "as_clinicalSignificance": float(as_clinical),
        "foldxDdq_raw": foldx_ddg,
        "most_severe_consequence": most_severe,
    }


def build_annotation_table(
    variant_ids: list[str],
    raw_records: dict[str, dict],
    foldx_lookup: dict[tuple[str, int, str, str], tuple[float, float]] | None = None,
) -> pa.Table:
    rows = [_transform_record(raw_records.get(variant), variant, foldx_lookup=foldx_lookup) for variant in variant_ids]
    if not rows:
        return pa.Table.from_pylist([], schema=ANNOTATION_SCHEMA)
    return pa.Table.from_pylist(rows, schema=ANNOTATION_SCHEMA)


def load_existing_annotation_table(output_path: str) -> pa.Table | None:
    if not path_exists(output_path):
        return None
    return open_parquet_dataset(output_path).to_table()


def load_existing_annotation_variants(output_path: str) -> set[str]:
    if not path_exists(output_path):
        return set()
    dataset = open_parquet_dataset(output_path)
    if "variant" not in dataset.schema.names:
        raise ValueError(f"Existing annotation dataset is missing 'variant': {output_path}")
    variants = set()
    for batch in dataset.scanner(columns=["variant"]).to_batches():
        variants.update(str(value).strip() for value in batch.column(0).to_pylist() if value)
    return variants


def merge_annotation_tables(existing_table: pa.Table | None, new_table: pa.Table) -> pa.Table:
    if existing_table is None or existing_table.num_rows == 0:
        return new_table
    if new_table.num_rows == 0:
        return existing_table

    frames = [
        existing_table.to_pandas(),
        new_table.to_pandas(),
    ]
    ordered_columns = [name for name, _ in ANNOTATION_FIELDS]
    for frame in frames:
        for column in ordered_columns:
            if column not in frame.columns:
                frame[column] = None
        frame = frame[ordered_columns]

    merged = pd.concat(
        [
            frames[0][ordered_columns],
            frames[1][ordered_columns],
        ],
        ignore_index=True,
    )
    merged = merged.drop_duplicates(subset=["variant"], keep="last")
    return pa.Table.from_pandas(merged, schema=ANNOTATION_SCHEMA, preserve_index=False)


def resolve_annotations_output(args) -> str:
    if args.annotations_output:
        return args.annotations_output

    bucket_base = _resolve_bucket_base(args)
    if bucket_base:
        return _join_bucket_path(
            _join_bucket_path(bucket_base, args.annotations_dir),
            args.annotations_output_name,
        )

    return DEFAULT_ANNOTATIONS_OUTPUT


# ============================================================================
# CLI / Main
# ============================================================================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run VIDRA Step 2 bulk VEP annotation from the analysis-ready manifest."
    )
    parser.add_argument("--analysis_manifest_dir", default=DEFAULT_ANALYSIS_MANIFEST_DIR)
    parser.add_argument("--vep_json_output", default=DEFAULT_VEP_JSON_OUTPUT)
    parser.add_argument("--bucket_name", default="", help="GCS bucket name without gs://")
    parser.add_argument("--bucket_uri", default="", help="Bucket base URI/path override, e.g. gs://bucket or /data/bucket")
    parser.add_argument("--annotations_dir", default="variant_annotations", help="Annotation directory under bucket base")
    parser.add_argument("--annotations_output_name", default="annotations.parquet", help="Output parquet filename inside annotation directory")
    parser.add_argument("--annotations_output", default="", help="Explicit output parquet path; overrides bucket-based output")
    parser.add_argument("--vep_dir_cache", default="", help="VEP cache root that contains homo_sapiens/<version>_GRCh38")
    parser.add_argument("--vep_cache_archive", default="", help="Optional tar/tar.gz archive for the VEP cache; useful for cloud runs")
    parser.add_argument("--vep_plugins_dir", default="/opt/vep/src/ensembl-vep/Plugins", help="Directory containing VEP plugin Perl modules")
    parser.add_argument("--vep_plugin_data_dir", default=PLUGIN_ROOT, help="Directory containing optional VEP plugin data files")
    parser.add_argument("--required_vep_plugins", default="AlphaMissense,CADD,REVEL", help="Comma-separated VEP plugins that must be available for the run to succeed")
    parser.add_argument("--vep_parallel", type=int, default=DEFAULT_VEP_PARALLEL, help="Concurrent VEP processes for the bulk run; use 0 for the auto setting")
    parser.add_argument("--vep_forks", type=int, default=DEFAULT_VEP_FORKS, help="VEP --fork value per process; use 0 for the auto setting")
    parser.add_argument("--vep_buffer_size", type=int, default=DEFAULT_VEP_BUFFER_SIZE, help="VEP --buffer_size value per process")
    parser.add_argument("--reuse_existing_annotations", default="true", help="Reuse and extend an existing annotations parquet when available")
    parser.add_argument("--foldx_file", default="", help="Optional FoldX lookup file; missing files are tolerated")
    return parser.parse_args()


def validate_runtime_args(args) -> None:
    if args.vep_parallel < 0:
        raise ValueError("--vep_parallel must be >= 0")
    if args.vep_forks < 0:
        raise ValueError("--vep_forks must be >= 0")
    if args.vep_buffer_size < 1:
        raise ValueError("--vep_buffer_size must be >= 1")


def main() -> None:
    args = parse_args()
    args.reuse_existing_annotations = parse_bool_flag(args.reuse_existing_annotations)
    args.required_vep_plugins = parse_plugin_list(args.required_vep_plugins)
    validate_runtime_args(args)
    with tempfile.TemporaryDirectory(prefix="vidra_annotation_runtime_") as runtime_tmp:
        log_progress("Preparing local VEP runtime resources")
        runtime_args = prepare_runtime_resources(args, Path(runtime_tmp))
        runtime_args.vep_plugin_args, runtime_args.vep_plugin_statuses = build_plugin_args(
            runtime_args.vep_plugin_data_dir,
            runtime_args.vep_plugins_dir,
            runtime_args.required_vep_plugins,
        )
        log_progress("VEP runtime resources ready")
        manifest_path = resolve_analysis_manifest_path(runtime_args)
        annotations_output = resolve_annotations_output(runtime_args)

        all_variant_ids = load_manifest_variant_ids(manifest_path, allow_empty=True)
        log_progress(
            f"Loaded {len(all_variant_ids):,} distinct variant(s) from manifest {manifest_path}"
        )

        existing_table = None
        existing_variant_ids: set[str] = set()
        if runtime_args.reuse_existing_annotations:
            existing_table = load_existing_annotation_table(annotations_output)
            existing_variant_ids = load_existing_annotation_variants(annotations_output)
            if existing_variant_ids:
                log(f"Existing annotations available: {len(existing_variant_ids)} variant(s)")

        missing_variant_ids = [variant for variant in all_variant_ids if variant not in existing_variant_ids]
        log_progress(
            f"Annotation reuse: existing={len(existing_variant_ids):,}, missing={len(missing_variant_ids):,}"
        )

        if not all_variant_ids:
            write_placeholder_json(runtime_args.vep_json_output)
            table = existing_table if existing_table is not None else build_annotation_table([], {})
            _write_table_to_path(table, annotations_output)
            log("No manifest variants require annotation for this run.")
            log(f"  - {annotations_output}")
            return

        if runtime_args.reuse_existing_annotations and not missing_variant_ids and existing_table is not None:
            write_placeholder_json(runtime_args.vep_json_output)
            log("All required annotation variants already exist; skipping VEP.")
            log(f"  - {annotations_output}")
            return

        annotate_variant_ids(
            runtime_args,
            "bulk_variant_annotation",
            missing_variant_ids,
            runtime_args.vep_json_output,
        )

        log_progress("Parsing bulk VEP JSON output")
        raw_records = parse_vep_json_outputs([runtime_args.vep_json_output])
        needed_foldx_keys = collect_needed_foldx_keys(raw_records)
        foldx_eligible_variants = count_foldx_eligible_variants(raw_records)
        log_progress(
            f"Bulk VEP parse complete: {len(raw_records):,} annotated variant record(s); "
            f"{foldx_eligible_variants:,} coding/protein consequence variant(s) eligible for FoldX lookup"
        )
        if runtime_args.foldx_file and path_exists(runtime_args.foldx_file):
            log_progress(
                f"Loading FoldX lookup for {len(needed_foldx_keys):,} protein consequence key(s)"
            )
            foldx_lookup = load_foldx_lookup(runtime_args.foldx_file, needed_foldx_keys)
            log(f"FoldX lookup loaded: {len(foldx_lookup)} matching entries")
        else:
            foldx_lookup = {}
            if runtime_args.foldx_file:
                log(f"WARN: FoldX file not found at {runtime_args.foldx_file}; using default foldx/plddt values")
            else:
                log("WARN: No FoldX file configured; using default foldx/plddt values")

        new_table = build_annotation_table(missing_variant_ids, raw_records, foldx_lookup=foldx_lookup)
        table = merge_annotation_tables(existing_table, new_table)

        log_progress(f"Writing merged annotation table with {table.num_rows:,} variant(s)")
        _write_table_to_path(table, annotations_output)

        log("Completed bulk annotation run.")
        log_progress(
            f"Annotation complete: wrote {table.num_rows:,} variant(s) to {annotations_output}"
        )
        log(f"  - {runtime_args.vep_json_output}")
        log(f"  - {annotations_output}")


if __name__ == "__main__":
    main()
