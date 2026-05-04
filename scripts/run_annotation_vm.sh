#!/usr/bin/env bash
set -euo pipefail

: "${GCP_PROJECT:?}"
: "${GCP_ZONE:?}"
: "${GCP_SERVICE_ACCOUNT_EMAIL:?}"
: "${ANNOTATION_IMAGE_GCLOUD:?}"
: "${GCLOUD_OUTDIR:?}"
: "${GCLOUD_VEP_CACHE_ARCHIVE:?}"
: "${GCLOUD_VEP_PLUGIN_DATA_DIR:?}"
: "${GCLOUD_REQUIRED_VEP_PLUGINS:?}"
: "${ANNOTATION_VM_NAME:?}"
: "${ANNOTATION_VM_MACHINE_TYPE:?}"
: "${ANNOTATION_VM_BOOT_DISK_SIZE:?}"
: "${ANNOTATION_VM_BOOT_DISK_TYPE:=pd-ssd}"
: "${ANNOTATION_VM_IMAGE_FAMILY:=debian-12}"
: "${ANNOTATION_VM_IMAGE_PROJECT:=debian-cloud}"
: "${ANNOTATION_VM_AUTO_DELETE:?}"
: "${ANNOTATION_VM_TIMEOUT_SECONDS:=10800}"
: "${ANNOTATION_VEP_PARALLEL:=0}"
: "${ANNOTATION_VEP_FORKS:=0}"
: "${ANNOTATION_VEP_BUFFER_SIZE:=10000}"
: "${GCLOUD_FOLDX_FILE:=}"

ANALYSIS_MANIFEST_DIR="vidra_analysis_ready_manifest"
ANNOTATIONS_DIR="variant_annotations"
ANNOTATIONS_OUTPUT_NAME="annotations.parquet"
RUN_TRACE_SUFFIX="$(date +%Y%m%d_%H%M%S)"
STATUS_ROOT="${GCLOUD_OUTDIR}/pipeline_info/annotation_vm/${RUN_TRACE_SUFFIX}"
VEP_JSON_OUTPUT="/tmp/vep_annotations.json"
VM_SUFFIX="${RUN_TRACE_SUFFIX//_/-}"
VM_NAME="${ANNOTATION_VM_NAME}-${VM_SUFFIX}"
STATUS_URI="${STATUS_ROOT}/status.txt"
PROGRESS_URI="${STATUS_ROOT}/progress.txt"
LOG_URI="${STATUS_ROOT}/annotation.log"

poll_status() {
    local status_uri="$1"
    local timeout_seconds="$2"
    local start="$SECONDS"
    local status

    while (( SECONDS - start < timeout_seconds )); do
        if status="$(gcloud storage cat "${status_uri}" 2>/dev/null | tr -d '\r\n')" && [[ -n "${status}" ]]; then
            printf '%s' "${status}"
            return 0
        fi
        sleep 30
    done

    return 1
}

delete_instance_if_present() {
    local vm_name="$1"
    gcloud compute instances delete "${vm_name}" \
        --project "${GCP_PROJECT}" \
        --zone "${GCP_ZONE}" \
        --quiet >/dev/null 2>&1 || true
}

startup_file="$(mktemp /tmp/vidra5_annotation_vm_XXXXXX)"
cat > "${startup_file}" <<STARTUP
#!/bin/bash
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
STATUS_URI="${STATUS_URI}"
PROGRESS_URI="${PROGRESS_URI}"
LOG_URI="${LOG_URI}"
GCP_PROJECT_VALUE="${GCP_PROJECT}"
GCP_ZONE_VALUE="${GCP_ZONE}"
VM_NAME_VALUE="${VM_NAME}"
ANNOTATION_VM_AUTO_DELETE_VALUE="${ANNOTATION_VM_AUTO_DELETE}"
LOG_FILE="/var/log/vidra5/annotation.log"
PROGRESS_FILE="/var/log/vidra5/progress.txt"
LOG_SYNC_INTERVAL_SECONDS="60"

mkdir -p /var/log/vidra5
exec > >(tee -a "\${LOG_FILE}") 2>&1
: > "\${PROGRESS_FILE}"

log_msg() {
    echo "[\$(date -u '+%Y-%m-%d %H:%M:%S UTC')] \$*"
}

sync_log_snapshot() {
    if [[ -f "\${LOG_FILE}" ]]; then
        gcloud storage cp "\${LOG_FILE}" "\${LOG_URI}" >/dev/null 2>&1 || true
    fi
}

sync_log_loop() {
    while true; do
        sync_log_snapshot
        sleep "\${LOG_SYNC_INTERVAL_SECONDS}"
    done
}

update_progress() {
    local msg="\$1"
    log_msg "PROGRESS: \${msg}"
    printf '%s - %s\n' "\${msg}" "\$(date -u '+%Y-%m-%d %H:%M:%S UTC')" >> "\${PROGRESS_FILE}"
    gcloud storage cp "\${PROGRESS_FILE}" "\${PROGRESS_URI}" >/dev/null 2>&1 || true
}

update_progress_from_annotation_line() {
    local line="\$1"
    case "\${line}" in
        *"[annotate_variants_cli] PROGRESS: "*)
            update_progress "\${line#*PROGRESS: }"
            ;;
    esac
}

cleanup() {
    status=\$?
    if [[ -n "\${LOG_SYNC_PID:-}" ]]; then
        kill "\${LOG_SYNC_PID}" >/dev/null 2>&1 || true
        wait "\${LOG_SYNC_PID}" 2>/dev/null || true
    fi
    sync_log_snapshot
    if [[ -f "\${PROGRESS_FILE}" ]]; then
        gcloud storage cp "\${PROGRESS_FILE}" "\${PROGRESS_URI}" >/dev/null 2>&1 || true
    fi
    if [[ \$status -eq 0 ]]; then
        printf 'SUCCESS\n' | gcloud storage cp - "\${STATUS_URI}" >/dev/null 2>&1 || true
    else
        printf 'FAILURE\n' | gcloud storage cp - "\${STATUS_URI}" >/dev/null 2>&1 || true
    fi
    if [[ "\${ANNOTATION_VM_AUTO_DELETE_VALUE}" == "true" ]]; then
        gcloud compute instances delete "\${VM_NAME_VALUE}" \
            --project "\${GCP_PROJECT_VALUE}" \
            --zone "\${GCP_ZONE_VALUE}" \
            --quiet >/dev/null 2>&1 || true
    fi
}
trap cleanup EXIT

sync_log_loop &
LOG_SYNC_PID=\$!

update_progress "Phase 1/4: Installing system packages"
apt-get update -qq
apt-get install -y -qq docker.io curl ca-certificates gnupg

install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /etc/apt/keyrings/cloud.google.gpg
echo "deb [signed-by=/etc/apt/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" > /etc/apt/sources.list.d/google-cloud-sdk.list
update_progress "Phase 2/4: Installing Google Cloud CLI and starting Docker"
apt-get update -qq
apt-get install -y -qq google-cloud-cli

systemctl start docker
gcloud auth configure-docker europe-west1-docker.pkg.dev --quiet
update_progress "Phase 3/4: Pulling annotation image"
docker pull "${ANNOTATION_IMAGE_GCLOUD}"

update_progress "Phase 4/4: Running annotation container"
docker run --rm "${ANNOTATION_IMAGE_GCLOUD}" \
    python3 /app/tools/annotate_variants_cli.py \
        --analysis_manifest_dir "${ANALYSIS_MANIFEST_DIR}" \
        --vep_json_output "${VEP_JSON_OUTPUT}" \
        --bucket_uri "${GCLOUD_OUTDIR}" \
        --annotations_dir "${ANNOTATIONS_DIR}" \
        --annotations_output_name "${ANNOTATIONS_OUTPUT_NAME}" \
        --vep_cache_archive "${GCLOUD_VEP_CACHE_ARCHIVE}" \
        --vep_plugins_dir "/opt/vep/src/ensembl-vep/Plugins" \
        --vep_plugin_data_dir "${GCLOUD_VEP_PLUGIN_DATA_DIR}" \
        --required_vep_plugins "${GCLOUD_REQUIRED_VEP_PLUGINS}" \
        --vep_parallel "${ANNOTATION_VEP_PARALLEL}" \
        --vep_forks "${ANNOTATION_VEP_FORKS}" \
        --vep_buffer_size "${ANNOTATION_VEP_BUFFER_SIZE}" \
        --foldx_file "${GCLOUD_FOLDX_FILE}" 2>&1 | while IFS= read -r line; do
    printf '%s\n' "\${line}"
    update_progress_from_annotation_line "\${line}"
done
STARTUP

echo "[gcloud_step2] Starting annotation VM ${VM_NAME} on ${ANNOTATION_VM_MACHINE_TYPE} with ${ANNOTATION_VM_BOOT_DISK_SIZE} ${ANNOTATION_VM_BOOT_DISK_TYPE} disk (timeout ${ANNOTATION_VM_TIMEOUT_SECONDS}s)"
echo "[gcloud_step2] Status URI: ${STATUS_URI}"
echo "[gcloud_step2] Progress URI: ${PROGRESS_URI}"
echo "[gcloud_step2] Log URI: ${LOG_URI}"

gcloud compute instances create "${VM_NAME}" \
    --project "${GCP_PROJECT}" \
    --zone "${GCP_ZONE}" \
    --machine-type "${ANNOTATION_VM_MACHINE_TYPE}" \
    --boot-disk-size "${ANNOTATION_VM_BOOT_DISK_SIZE}" \
    --boot-disk-type "${ANNOTATION_VM_BOOT_DISK_TYPE}" \
    --image-family "${ANNOTATION_VM_IMAGE_FAMILY}" \
    --image-project "${ANNOTATION_VM_IMAGE_PROJECT}" \
    --service-account "${GCP_SERVICE_ACCOUNT_EMAIL}" \
    --scopes "https://www.googleapis.com/auth/cloud-platform" \
    --metadata-from-file startup-script="${startup_file}"

rm -f "${startup_file}"

if ! status="$(poll_status "${STATUS_URI}" "${ANNOTATION_VM_TIMEOUT_SECONDS}")"; then
    echo "[gcloud_step2] Annotation VM timed out waiting for ${STATUS_URI}" >&2
    delete_instance_if_present "${VM_NAME}"
    exit 1
fi

if [[ "${status}" != "SUCCESS" ]]; then
    echo "[gcloud_step2] Annotation VM reported ${status}; see ${LOG_URI}" >&2
    delete_instance_if_present "${VM_NAME}"
    exit 1
fi

echo "[gcloud_step2] Annotation VM completed successfully; log: ${LOG_URI}"
