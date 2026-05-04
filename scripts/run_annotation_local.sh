#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?}"
: "${HOST_UID:?}"
: "${HOST_GID:?}"
: "${ANNOTATION_IMAGE_LOCAL:?}"
: "${LOCAL_OUTDIR:?}"
: "${LOCAL_VEP_RESOURCES_DIR:?}"
: "${LOCAL_REQUIRED_VEP_PLUGINS:?}"
: "${ANNOTATION_VEP_PARALLEL:=0}"
: "${ANNOTATION_VEP_FORKS:=0}"
: "${ANNOTATION_VEP_BUFFER_SIZE:=10000}"
: "${SAVE_INTERMEDIATE_FILES:=false}"

ANALYSIS_MANIFEST_DIR="vidra_analysis_ready_manifest"
ANNOTATIONS_DIR="variant_annotations"
ANNOTATIONS_OUTPUT_NAME="annotations.parquet"

case "$(printf '%s' "${SAVE_INTERMEDIATE_FILES}" | tr '[:upper:]' '[:lower:]')" in
    true|1|yes)
        LOCAL_ANNOTATION_INTERMEDIATE_DIR="${LOCAL_OUTDIR}/intermediate/annotate_variants"
        mkdir -p "${LOCAL_ANNOTATION_INTERMEDIATE_DIR}"
        VEP_JSON_OUTPUT="${LOCAL_ANNOTATION_INTERMEDIATE_DIR}/vep_annotations.json"
        ;;
    *)
        VEP_JSON_OUTPUT="/tmp/vep_annotations.json"
        ;;
esac

docker run --rm --platform=linux/amd64 \
    --user "${HOST_UID}:${HOST_GID}" \
    -v "${PROJECT_DIR}:${PROJECT_DIR}" \
    -v "${LOCAL_VEP_RESOURCES_DIR}:${LOCAL_VEP_RESOURCES_DIR}:ro" \
    "${ANNOTATION_IMAGE_LOCAL}" \
    python3 /app/tools/annotate_variants_cli.py \
        --analysis_manifest_dir "${ANALYSIS_MANIFEST_DIR}" \
        --vep_json_output "${VEP_JSON_OUTPUT}" \
        --bucket_uri "${LOCAL_OUTDIR}" \
        --annotations_dir "${ANNOTATIONS_DIR}" \
        --annotations_output_name "${ANNOTATIONS_OUTPUT_NAME}" \
        --vep_dir_cache "${LOCAL_VEP_RESOURCES_DIR}" \
        --vep_plugins_dir "/opt/vep/src/ensembl-vep/Plugins" \
        --vep_plugin_data_dir "${LOCAL_VEP_RESOURCES_DIR}" \
        --required_vep_plugins "${LOCAL_REQUIRED_VEP_PLUGINS}" \
        --vep_parallel "${ANNOTATION_VEP_PARALLEL}" \
        --vep_forks "${ANNOTATION_VEP_FORKS}" \
        --vep_buffer_size "${ANNOTATION_VEP_BUFFER_SIZE}" \
        --foldx_file "${LOCAL_FOLDX_FILE}"
