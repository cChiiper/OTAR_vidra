#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?}"
: "${HOST_UID:?}"
: "${HOST_GID:?}"
: "${PREPARE_IMAGE_LOCAL:?}"
: "${LOCAL_GENES:?}"
: "${LOCAL_OUTDIR:?}"
: "${LOCAL_RELEASE_ROOT:?}"
: "${LOCAL_VIDRA_DATA_ROOT:?}"
: "${LOCAL_AZ_VARIANTS_FILE:?}"
: "${COLOCALISATION_THRESHOLD:?}"
: "${LOCAL_PREPARE_SPARK_MASTER:?}"
: "${LOCAL_PREPARE_SPARK_DRIVER_MEMORY:?}"
: "${LOCAL_PREPARE_SPARK_SHUFFLE_PARTITIONS:?}"
: "${LOCAL_PREPARE_SPARK_DEFAULT_PARALLELISM:?}"

LOCAL_AZ_BURDEN_BINARY_DIR="${LOCAL_AZ_BURDEN_BINARY_DIR:-}"
LOCAL_AZ_BURDEN_QUANTITATIVE_DIR="${LOCAL_AZ_BURDEN_QUANTITATIVE_DIR:-}"
LOCAL_TARGET_DATA_DIR="${LOCAL_TARGET_DATA_DIR:-${LOCAL_RELEASE_ROOT}/target}"
LOCAL_LEAD_VARIANT_EFFECT_DIR="${LOCAL_LEAD_VARIANT_EFFECT_DIR:-}"
ANALYSIS_READY_DIR="vidra_analysis_ready"
ANALYSIS_MANIFEST_DIR="vidra_analysis_ready_manifest"
LOCAL_PREPARE_RUNTIME_ROOT="${LOCAL_OUTDIR}/.prepare_spark_runtime"
LOCAL_PREPARE_USER_HOME="${LOCAL_PREPARE_RUNTIME_ROOT}/home"
LOCAL_PREPARE_SPARK_LOCAL_DIRS="${LOCAL_PREPARE_RUNTIME_ROOT}/spark_local"
CODING_OUTPUT="${LOCAL_OUTDIR}/GWAS_coding_variants_from_cs.parquet"
CLINVAR_OUTPUT="${LOCAL_OUTDIR}/clinvar_variants.parquet"

mkdir -p "${LOCAL_OUTDIR}"
mkdir -p "${LOCAL_PREPARE_USER_HOME}" "${LOCAL_PREPARE_SPARK_LOCAL_DIRS}"

docker run --rm --platform=linux/amd64 \
    --user root \
    -e "HOST_UID=${HOST_UID}" \
    -e "HOST_GID=${HOST_GID}" \
    -e "PROJECT_DIR=${PROJECT_DIR}" \
    -e "HOME=${LOCAL_PREPARE_USER_HOME}" \
    -e "VIDRA_SPARK_USER_HOME=${LOCAL_PREPARE_USER_HOME}" \
    -e "VIDRA_SPARK_LOCAL_DIRS=${LOCAL_PREPARE_SPARK_LOCAL_DIRS}" \
    -e "SPARK_USER=vidra_local" \
    -e "HADOOP_USER_NAME=vidra_local" \
    -e "USER=vidra_local" \
    -e "LOGNAME=vidra_local" \
    -e "LOCAL_GENES=${LOCAL_GENES}" \
    -e "LOCAL_OUTDIR=${LOCAL_OUTDIR}" \
    -e "LOCAL_RELEASE_ROOT=${LOCAL_RELEASE_ROOT}" \
    -e "LOCAL_VIDRA_DATA_ROOT=${LOCAL_VIDRA_DATA_ROOT}" \
    -e "LOCAL_AZ_VARIANTS_FILE=${LOCAL_AZ_VARIANTS_FILE}" \
    -e "LOCAL_LEAD_VARIANT_EFFECT_DIR=${LOCAL_LEAD_VARIANT_EFFECT_DIR}" \
    -e "LOCAL_AZ_BURDEN_BINARY_DIR=${LOCAL_AZ_BURDEN_BINARY_DIR}" \
    -e "LOCAL_AZ_BURDEN_QUANTITATIVE_DIR=${LOCAL_AZ_BURDEN_QUANTITATIVE_DIR}" \
    -e "LOCAL_TARGET_DATA_DIR=${LOCAL_TARGET_DATA_DIR}" \
    -e "COLOCALISATION_THRESHOLD=${COLOCALISATION_THRESHOLD}" \
    -e "LOCAL_PREPARE_SPARK_MASTER=${LOCAL_PREPARE_SPARK_MASTER}" \
    -e "LOCAL_PREPARE_SPARK_DRIVER_MEMORY=${LOCAL_PREPARE_SPARK_DRIVER_MEMORY}" \
    -e "LOCAL_PREPARE_SPARK_SHUFFLE_PARTITIONS=${LOCAL_PREPARE_SPARK_SHUFFLE_PARTITIONS}" \
    -e "LOCAL_PREPARE_SPARK_DEFAULT_PARALLELISM=${LOCAL_PREPARE_SPARK_DEFAULT_PARALLELISM}" \
    -e "ANALYSIS_READY_DIR=${ANALYSIS_READY_DIR}" \
    -e "ANALYSIS_MANIFEST_DIR=${ANALYSIS_MANIFEST_DIR}" \
    -e "CODING_OUTPUT=${CODING_OUTPUT}" \
    -e "CLINVAR_OUTPUT=${CLINVAR_OUTPUT}" \
    -v "${PROJECT_DIR}:${PROJECT_DIR}" \
    -v "${LOCAL_RELEASE_ROOT}:${LOCAL_RELEASE_ROOT}:ro" \
    -v "${LOCAL_VIDRA_DATA_ROOT}:${LOCAL_VIDRA_DATA_ROOT}:ro" \
    "${PREPARE_IMAGE_LOCAL}" \
    bash -lc '
        cleanup() {
            chown -R "${HOST_UID}:${HOST_GID}" "${LOCAL_OUTDIR}" || true
        }
        trap cleanup EXIT

        raw_az_burden_args=()
        if [[ -n "${LOCAL_AZ_BURDEN_BINARY_DIR}" || -n "${LOCAL_AZ_BURDEN_QUANTITATIVE_DIR}" ]]; then
            raw_az_burden_args+=(
                --az_burden_binary_dir "${LOCAL_AZ_BURDEN_BINARY_DIR}"
                --az_burden_quantitative_dir "${LOCAL_AZ_BURDEN_QUANTITATIVE_DIR}"
                --target_data_dir "${LOCAL_TARGET_DATA_DIR}"
            )
        fi

        lead_variant_effect_args=()
        if [[ -n "${LOCAL_LEAD_VARIANT_EFFECT_DIR}" ]]; then
            lead_variant_effect_args+=(
                --lead_variant_effect_dir "${LOCAL_LEAD_VARIANT_EFFECT_DIR}"
            )
        fi

        /opt/conda/bin/python "${PROJECT_DIR}/tools/prepare_analysis_input.py" \
            --genes "${LOCAL_GENES}" \
            --colocalisation_threshold "${COLOCALISATION_THRESHOLD}" \
            --coloc_data_dir "${LOCAL_RELEASE_ROOT}/colocalisation" \
            --credible_set_dir "${LOCAL_RELEASE_ROOT}/credible_set" \
            --study_data_dir "${LOCAL_RELEASE_ROOT}/study" \
            --variant_data_dir "${LOCAL_RELEASE_ROOT}/variant" \
            --burden_evidence_dir "${LOCAL_RELEASE_ROOT}/evidence_gene_burden" \
            --clinvar_evidence_dir "${LOCAL_RELEASE_ROOT}/evidence_eva" \
            --spark_master "${LOCAL_PREPARE_SPARK_MASTER}" \
            --spark_driver_memory "${LOCAL_PREPARE_SPARK_DRIVER_MEMORY}" \
            --spark_shuffle_partitions "${LOCAL_PREPARE_SPARK_SHUFFLE_PARTITIONS}" \
            --spark_default_parallelism "${LOCAL_PREPARE_SPARK_DEFAULT_PARALLELISM}" \
            --az_variants_file "${LOCAL_AZ_VARIANTS_FILE}" \
            --az_mapping_file "${LOCAL_VIDRA_DATA_ROOT}/az_phewas_to_efo.tsv" \
            --az_gene_map_file "${LOCAL_VIDRA_DATA_ROOT}/all_human_protein_coding_genes.csv" \
            "${lead_variant_effect_args[@]}" \
            "${raw_az_burden_args[@]}" \
            --bucket_uri "${LOCAL_OUTDIR}" \
            --analysis_ready_dir "${ANALYSIS_READY_DIR}" \
            --analysis_manifest_dir "${ANALYSIS_MANIFEST_DIR}" \
            --coding_output "${CODING_OUTPUT}" \
            --clinvar_output "${CLINVAR_OUTPUT}"
    '
