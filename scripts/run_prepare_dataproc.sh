#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_DIR:?}"
: "${PREPARE_IMAGE_GCLOUD:?}"
: "${GCP_PROJECT:?}"
: "${GCP_REGION:?}"
: "${GCP_SERVICE_ACCOUNT_EMAIL:?}"
: "${DATAPROC_DEPS_BUCKET:?}"
: "${DATAPROC_BATCH_PREFIX:?}"
: "${DATAPROC_STEP1_TTL:?}"
: "${DATAPROC_STEP1_PROPERTIES:?}"
: "${GCLOUD_OUTDIR:?}"
: "${GCLOUD_GENES:?}"
: "${GCLOUD_RELEASE_ROOT:?}"
: "${GCLOUD_VIDRA_DATA_ROOT:?}"
: "${GCLOUD_AZ_VARIANTS_FILE:?}"
: "${COLOCALISATION_THRESHOLD:?}"

GCLOUD_AZ_BURDEN_BINARY_DIR="${GCLOUD_AZ_BURDEN_BINARY_DIR:-}"
GCLOUD_AZ_BURDEN_QUANTITATIVE_DIR="${GCLOUD_AZ_BURDEN_QUANTITATIVE_DIR:-}"
GCLOUD_TARGET_DATA_DIR="${GCLOUD_TARGET_DATA_DIR:-${GCLOUD_RELEASE_ROOT}/target}"
GCLOUD_LEAD_VARIANT_EFFECT_DIR="${GCLOUD_LEAD_VARIANT_EFFECT_DIR:-}"
ANALYSIS_READY_DIR="vidra_analysis_ready"
ANALYSIS_MANIFEST_DIR="vidra_analysis_ready_manifest"
CODING_OUTPUT="${GCLOUD_OUTDIR}/GWAS_coding_variants_from_cs.parquet"
CLINVAR_OUTPUT="${GCLOUD_OUTDIR}/clinvar_variants.parquet"

trace_suffix="$(date +%Y%m%d_%H%M%S)"
metadata_dir="${PROJECT_DIR}/results/launcher_logs"
metadata_output="${metadata_dir}/prepare_analysis_input_${trace_suffix}.json"
mkdir -p "${metadata_dir}"

submit_cmd=(
    python3 "${PROJECT_DIR}/tools/dataproc_batch_submit.py"
    --step_name "prepare_analysis_input"
    --script_path "${PROJECT_DIR}/tools/prepare_analysis_input.py"
    --project "${GCP_PROJECT}"
    --region "${GCP_REGION}"
    --deps_bucket "${DATAPROC_DEPS_BUCKET}"
    --container_image "${PREPARE_IMAGE_GCLOUD}"
    --service_account "${GCP_SERVICE_ACCOUNT_EMAIL}"
    --batch_prefix "${DATAPROC_BATCH_PREFIX}"
    --trace_suffix "${trace_suffix}"
    --ttl "${DATAPROC_STEP1_TTL}"
    --properties "${DATAPROC_STEP1_PROPERTIES}"
    --metadata_output "${metadata_output}"
    --script_arg=--genes
    --script_arg=${GCLOUD_GENES}
    --script_arg=--colocalisation_threshold
    --script_arg=${COLOCALISATION_THRESHOLD}
    --script_arg=--coloc_data_dir
    --script_arg=${GCLOUD_RELEASE_ROOT}/colocalisation
    --script_arg=--credible_set_dir
    --script_arg=${GCLOUD_RELEASE_ROOT}/credible_set
    --script_arg=--study_data_dir
    --script_arg=${GCLOUD_RELEASE_ROOT}/study
    --script_arg=--variant_data_dir
    --script_arg=${GCLOUD_RELEASE_ROOT}/variant
    --script_arg=--burden_evidence_dir
    --script_arg=${GCLOUD_RELEASE_ROOT}/evidence_gene_burden
    --script_arg=--clinvar_evidence_dir
    --script_arg=${GCLOUD_RELEASE_ROOT}/evidence_eva
    --script_arg=--az_variants_file
    --script_arg=${GCLOUD_AZ_VARIANTS_FILE}
    --script_arg=--az_mapping_file
    --script_arg=${GCLOUD_VIDRA_DATA_ROOT}/az_phewas_to_efo.tsv
    --script_arg=--az_gene_map_file
    --script_arg=${GCLOUD_VIDRA_DATA_ROOT}/all_human_protein_coding_genes.csv
    --script_arg=--bucket_uri
    --script_arg=${GCLOUD_OUTDIR}
    --script_arg=--gcp_project
    --script_arg=${GCP_PROJECT}
    --script_arg=--analysis_ready_dir
    --script_arg=${ANALYSIS_READY_DIR}
    --script_arg=--analysis_manifest_dir
    --script_arg=${ANALYSIS_MANIFEST_DIR}
    --script_arg=--coding_output
    --script_arg=${CODING_OUTPUT}
    --script_arg=--clinvar_output
    --script_arg=${CLINVAR_OUTPUT}
)

if [[ -n "${GCLOUD_AZ_BURDEN_BINARY_DIR}" || -n "${GCLOUD_AZ_BURDEN_QUANTITATIVE_DIR}" ]]; then
    submit_cmd+=(
        --script_arg=--az_burden_binary_dir
        --script_arg=${GCLOUD_AZ_BURDEN_BINARY_DIR}
        --script_arg=--az_burden_quantitative_dir
        --script_arg=${GCLOUD_AZ_BURDEN_QUANTITATIVE_DIR}
        --script_arg=--target_data_dir
        --script_arg=${GCLOUD_TARGET_DATA_DIR}
    )
fi

if [[ -n "${GCLOUD_LEAD_VARIANT_EFFECT_DIR}" ]]; then
    submit_cmd+=(
        --script_arg=--lead_variant_effect_dir
        --script_arg=${GCLOUD_LEAD_VARIANT_EFFECT_DIR}
    )
fi

"${submit_cmd[@]}"
