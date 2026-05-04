from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_makefile_and_base_config_use_hybrid_runtime_names():
    base_config = _read(PROJECT_DIR / "conf" / "base.config")
    makefile = _read(PROJECT_DIR / "Makefile")
    make_defaults = _read(PROJECT_DIR / "conf" / "make.defaults.mk")

    assert "prepare_image            = 'vidra5-spark-dataproc:latest'" in base_config
    assert "annotation_image         = 'vidra5-annotation:latest'" in base_config
    assert "bayes_image              = 'vidra5-bayes:latest'" in base_config
    assert "prepare_image_registry   = 'europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra5-spark-dataproc:latest'" in base_config
    assert "annotation_image_registry = 'europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra5-annotation:latest'" in base_config
    assert "bayes_image_registry     = 'europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra5-bayes:latest'" in base_config

    assert "include conf/make.defaults.mk" in makefile
    assert "-include conf/make.local.mk" in makefile
    assert "PREPARE_IMAGE_LOCAL ?= vidra5-spark-dataproc:latest" in make_defaults
    assert "ANNOTATION_IMAGE_LOCAL ?= vidra5-annotation:latest" in make_defaults
    assert "BAYES_IMAGE_LOCAL ?= vidra5-bayes:latest" in make_defaults
    assert "SAVE_INTERMEDIATE_FILES ?= false" in make_defaults
    assert "ANNOTATION_VEP_PARALLEL ?= 4" in make_defaults
    assert "ANNOTATION_VEP_FORKS ?= 2" in make_defaults
    assert "ANNOTATION_VEP_BUFFER_SIZE ?= 5000" in make_defaults
    assert "LOCAL_RELEASE_ROOT ?=" in make_defaults
    assert "OT_RELEASE_RSYNC_BASE ?= rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output" in make_defaults
    assert "OT_RELEASE_RSYNC_FOLDERS ?= colocalisation credible_set study variant evidence_gene_burden evidence_eva target" in make_defaults
    assert "LOCAL_AZ_PHEWAS_INPUT_ROOT ?= $(LOCAL_VIDRA_DATA_ROOT)/TEMP" in make_defaults
    assert "LOCAL_PROTEIN_CODING_ENSG_FILE ?= $(LOCAL_VIDRA_DATA_ROOT)/all_protein_coding_ENSG.csv" in make_defaults
    assert "DATAPROC_STEP1_TTL ?= 14400s" in make_defaults
    assert "ANNOTATION_VM_MACHINE_TYPE ?= e2-standard-8" in make_defaults
    assert "ANNOTATION_VM_BOOT_DISK_SIZE ?= 250GB" in make_defaults
    assert "ANNOTATION_VM_TIMEOUT_SECONDS ?= 18000" in make_defaults
    assert "ANNOTATION_VM_BOOT_DISK_TYPE ?= pd-ssd" in make_defaults
    assert "ANNOTATION_VM_IMAGE_FAMILY ?= debian-12" in make_defaults
    assert "ANNOTATION_VM_IMAGE_PROJECT ?= debian-cloud" in make_defaults
    assert "docker_build_prepare:" in makefile
    assert "docker_build_annotation:" in makefile
    assert "docker_build_bayes:" in makefile
    assert "local_sync_ot_release_inputs:" in makefile
    assert "local_convert_az_phewas_parquet:" in makefile
    assert "local_all_protein_coding_ensg:" in makefile
    assert "tools/export_protein_coding_ensg.py" in makefile
    assert "scripts/convert_az_phewas_parquet.py" in makefile
    assert 'rsync -rpltvz "$(OT_RELEASE_RSYNC_BASE)/$$folder" "$(LOCAL_RELEASE_ROOT)/"; \\' in makefile
    assert "local_step1:" in makefile
    assert "local_step2:" in makefile
    assert "gcloud_step1:" in makefile
    assert "gcloud_step2:" in makefile
    assert "nextflow_local_step3:" in makefile
    assert '--outdir "$(LOCAL_OUTDIR)" --bucket_uri "$(LOCAL_OUTDIR)"' in makefile
    assert "nextflow_gcloud_test_step3:" in makefile
    assert "nextflow_dataproc_test" not in makefile


def test_main_workflow_only_runs_step3_per_gene():
    workflow = _read(PROJECT_DIR / "main.nf")
    per_gene_module = _read(PROJECT_DIR / "modules" / "shared" / "run_bayesian_analysis_per_gene.nf")

    assert "RUN_BAYESIAN_ANALYSIS_PER_GENE" in workflow
    assert "as_gene=*" in workflow
    assert "type: 'dir'" in workflow
    assert "annotations_output_name" in workflow
    assert ".join(requestedGenes)" in workflow
    assert ".combine(requestedGenes)" not in workflow

    assert "SELECT_GENES_TO_RUN" not in workflow
    assert "PREPARE_ANALYSIS_INPUT" not in workflow
    assert "ANNOTATE_VARIANTS" not in workflow
    assert "DATAPROC" not in workflow
    assert "REFRESH_COMPLETED_GENES_MANIFEST" not in workflow

    assert '--analysis_dir "${gene_partition}"' in per_gene_module
    assert '--annotations_file "${annotations_file}"' in per_gene_module
    assert '--gene "${gene}"' in per_gene_module
    assert '--output_dir "${resultsDir}"' in per_gene_module


def test_only_local_test_and_gcloud_profiles_remain():
    nextflow_config = _read(PROJECT_DIR / "nextflow.config")

    assert "local {" in nextflow_config
    assert "test {" in nextflow_config
    assert "gcloud_test {" in nextflow_config
    assert "gcloud_prod {" in nextflow_config
    assert "gcloud_placeholder" not in nextflow_config
    assert "dataproc_test" not in nextflow_config
    assert "dataproc_prod" not in nextflow_config


def test_legacy_nextflow_module_directories_are_gone():
    assert not (PROJECT_DIR / "modules" / "local").exists()
    assert not (PROJECT_DIR / "modules" / "dataproc").exists()


def test_step3_configs_are_small_and_no_longer_reference_step0_or_spark():
    base_config = _read(PROJECT_DIR / "conf" / "base.config")
    local_config = _read(PROJECT_DIR / "conf" / "local.config")
    test_config = _read(PROJECT_DIR / "conf" / "test.config")
    gcloud_base = _read(PROJECT_DIR / "conf" / "gcloud_base.config")
    gcloud_test = _read(PROJECT_DIR / "conf" / "gcloud_test.config")
    gcloud_prod = _read(PROJECT_DIR / "conf" / "gcloud_prod.config")
    reporting_config = _read(PROJECT_DIR / "conf" / "reporting.config")
    gcloud_google = _read(PROJECT_DIR / "conf" / "gcloud_google.config")

    assert "select_genes_script" not in base_config
    assert "skip_completed_genes" not in base_config
    assert "reuse_existing_annotations" not in base_config
    assert "reuse_existing_results" not in base_config
    assert "spark_master" not in base_config
    assert "max_stan_tasks" not in base_config
    assert "manifest_script" not in base_config
    assert "completed_genes_dir" not in base_config
    assert "completed_genes_name" not in base_config

    assert "container     = params.bayes_image" in local_config
    assert "withName: RUN_BAYESIAN_ANALYSIS_PER_GENE" in local_config
    assert "includeConfig 'reporting.config'" in local_config
    assert 'genes                    = "${projectDir}/testdata/genes_ensembl_small.csv"' in test_config
    assert "REFRESH_COMPLETED_GENES_MANIFEST" not in test_config
    assert "includeConfig 'reporting.config'" in test_config
    assert "container     = params.bayes_image_registry" in gcloud_base

    assert "withName: RUN_BAYESIAN_ANALYSIS_PER_GENE" in gcloud_test
    assert "errorStrategy = 'retry'" in gcloud_test
    assert "maxRetries    = 1" in gcloud_test
    assert "REFRESH_COMPLETED_GENES_MANIFEST" not in gcloud_test
    assert "includeConfig 'gcloud_google.config'" in gcloud_test
    assert "includeConfig 'reporting.config'" in gcloud_test

    assert "REFRESH_COMPLETED_GENES_MANIFEST" not in gcloud_prod
    assert "includeConfig 'gcloud_google.config'" in gcloud_prod
    assert "includeConfig 'reporting.config'" in gcloud_prod

    assert 'execution_trace_${params.trace_report_suffix}.txt' in reporting_config
    assert 'execution_report_${params.trace_report_suffix}.html' in reporting_config
    assert 'execution_timeline_${params.trace_report_suffix}.html' in reporting_config
    assert 'pipeline_dag_${params.trace_report_suffix}.html' in reporting_config
    assert 'location = params.gcp_region' in gcloud_google
    assert 'project  = params.gcp_project' in gcloud_google
    assert 'serviceAccountEmail = params.gcp_service_account_email' in gcloud_google


def test_local_prepare_launcher_runs_single_in_memory_step1_path():
    launcher = _read(PROJECT_DIR / "scripts" / "run_prepare_local.sh")
    makefile = _read(PROJECT_DIR / "Makefile")
    make_defaults = _read(PROJECT_DIR / "conf" / "make.defaults.mk")

    assert '--user root' in launcher
    assert '-e "SPARK_USER=vidra_local"' in launcher
    assert '-e "HADOOP_USER_NAME=vidra_local"' in launcher
    assert '-e "PROJECT_DIR=${PROJECT_DIR}"' in launcher
    assert '/opt/conda/bin/python "${PROJECT_DIR}/tools/prepare_analysis_input.py"' in launcher
    assert 'CODING_OUTPUT="${LOCAL_OUTDIR}/GWAS_coding_variants_from_cs.parquet"' in launcher
    assert 'CLINVAR_OUTPUT="${LOCAL_OUTDIR}/clinvar_variants.parquet"' in launcher
    assert 'SAVE_INTERMEDIATE_FILES' not in launcher
    assert 'LOCAL_PREPARE_STAGE_DIR' not in launcher
    assert '--coloc_output' not in launcher
    assert '--gwas_output' not in launcher
    assert '--qtl_output' not in launcher
    assert '--burden_output' not in launcher
    assert '--az_output' not in launcher
    assert 'LOCAL_AZ_VARIANTS_XZ_FILE ?= $(LOCAL_VIDRA_DATA_ROOT)/azphewas-com-470k-exwas-binary.csv.xz' in make_defaults
    assert 'LOCAL_AZ_VARIANTS_FILE ?= $(LOCAL_VIDRA_DATA_ROOT)/azphewas-com-470k-exwas-binary.csv.bz2' in make_defaults
    assert 'LOCAL_LEAD_VARIANT_EFFECT_DIR ?=' in make_defaults
    assert 'LOCAL_AZ_BURDEN_BINARY_XZ_FILE ?= $(LOCAL_AZ_PHEWAS_INPUT_ROOT)/azphewas-com-470k-phewas-binary.csv.xz' in make_defaults
    assert 'LOCAL_AZ_BURDEN_QUANTITATIVE_XZ_FILE ?= $(LOCAL_AZ_PHEWAS_INPUT_ROOT)/azphewas-com-470k-phewas-quantitative.csv.xz' in make_defaults
    assert 'LOCAL_AZ_BURDEN_BINARY_DIR ?= $(LOCAL_VIDRA_DATA_ROOT)/azphewas-com-470k-phewas-binary' in make_defaults
    assert 'LOCAL_AZ_BURDEN_QUANTITATIVE_DIR ?= $(LOCAL_VIDRA_DATA_ROOT)/azphewas-com-470k-phewas-quantitative' in make_defaults
    assert 'local_convert_az_phewas_parquet:' in makefile
    assert '--binary-input "$(LOCAL_AZ_BURDEN_BINARY_XZ_FILE)"' in makefile
    assert '--binary-output-dir "$(LOCAL_AZ_BURDEN_BINARY_DIR)"' in makefile
    assert '--quantitative-input "$(LOCAL_AZ_BURDEN_QUANTITATIVE_XZ_FILE)"' in makefile
    assert '--quantitative-output-dir "$(LOCAL_AZ_BURDEN_QUANTITATIVE_DIR)"' in makefile
    assert 'lead_variant_effect_args=()' in launcher
    assert '--az_variants_file "${LOCAL_AZ_VARIANTS_FILE}"' in launcher
    assert '--lead_variant_effect_dir "${LOCAL_LEAD_VARIANT_EFFECT_DIR}"' in launcher
    assert 'raw_az_burden_args=()' in launcher
    assert '--az_burden_binary_dir "${LOCAL_AZ_BURDEN_BINARY_DIR}"' in launcher
    assert '--target_data_dir "${LOCAL_TARGET_DATA_DIR}"' in launcher
    assert 'chown -R "${HOST_UID}:${HOST_GID}" "${LOCAL_OUTDIR}" || true' in launcher
    assert 'local_convert_az_variants_bz2:' in makefile
    assert 'xz -dc "$(LOCAL_AZ_VARIANTS_XZ_FILE)" | bzip2 -c > "$(LOCAL_AZ_VARIANTS_FILE)"' in makefile
    assert 'local_step1:' in makefile
    assert 'local_step1: $(LOCAL_AZ_VARIANTS_FILE)' not in makefile
    assert 'SAVE_INTERMEDIATE_FILES="$(SAVE_INTERMEDIATE_FILES)" bash scripts/run_prepare_local.sh' not in makefile


def test_gcloud_prepare_launcher_runs_single_in_memory_step1_path():
    launcher = _read(PROJECT_DIR / "scripts" / "run_prepare_dataproc.sh")
    submit_helper = _read(PROJECT_DIR / "tools" / "dataproc_batch_submit.py")
    makefile = _read(PROJECT_DIR / "Makefile")
    make_defaults = _read(PROJECT_DIR / "conf" / "make.defaults.mk")

    assert 'submit_cmd=(' in launcher
    assert '--script_arg=--genes' in launcher
    assert 'EXTRA_SCRIPT_ARGS' not in launcher
    assert 'GCLOUD_PREPARE_STAGE_DIR' not in launcher
    assert '--script_arg=--coloc_output' not in launcher
    assert '--script_arg=--gwas_output' not in launcher
    assert '--script_arg=--qtl_output' not in launcher
    assert '--script_arg=--burden_output' not in launcher
    assert '--script_arg=--az_output' not in launcher
    assert 'GCLOUD_AZ_VARIANTS_FILE ?= $(GCLOUD_VIDRA_DATA_ROOT)/azphewas-com-470k-exwas-binary.csv.bz2' in make_defaults
    assert 'GCLOUD_LEAD_VARIANT_EFFECT_DIR ?= gs://ot_orchestration/gentropy_manuscript/data/26.03/intermediate/lead_variant_effect' in make_defaults
    assert 'GCLOUD_AZ_BURDEN_BINARY_DIR ?= gs://otar000-evidence_input/GeneBurden/data_files/azphewas-com-470k-phewas-binary' in make_defaults
    assert 'GCLOUD_AZ_BURDEN_QUANTITATIVE_DIR ?= gs://otar000-evidence_input/GeneBurden/data_files/azphewas-com-470k-phewas-quantitative' in make_defaults
    assert '--script_arg=${GCLOUD_AZ_VARIANTS_FILE}' in launcher
    assert '--script_arg=--lead_variant_effect_dir' in launcher
    assert '--script_arg=${GCLOUD_LEAD_VARIANT_EFFECT_DIR}' in launcher
    assert '--script_arg=--az_burden_binary_dir' in launcher
    assert '--script_arg=${GCLOUD_AZ_BURDEN_BINARY_DIR}' in launcher
    assert '--script_arg=${GCLOUD_TARGET_DATA_DIR}' in launcher
    assert '"${submit_cmd[@]}"' in launcher
    assert 'SAVE_INTERMEDIATE_FILES="$(SAVE_INTERMEDIATE_FILES)" bash scripts/run_prepare_dataproc.sh' not in makefile
    assert '"--async"' not in submit_helper
    assert '"--wait"' not in submit_helper
    assert 'console_url =' in submit_helper
    assert 'wait_command' in submit_helper


def test_annotation_launcher_uses_optional_intermediate_json_outputs():
    launcher = _read(PROJECT_DIR / "scripts" / "run_annotation_local.sh")
    makefile = _read(PROJECT_DIR / "Makefile")

    assert 'SAVE_INTERMEDIATE_FILES' in launcher
    assert 'ANALYSIS_MANIFEST_DIR="vidra_analysis_ready_manifest"' in launcher
    assert 'ANNOTATION_VEP_PARALLEL' in launcher
    assert 'ANNOTATION_VEP_FORKS' in launcher
    assert 'ANNOTATION_VEP_BUFFER_SIZE' in launcher
    assert 'LOCAL_ANNOTATION_INTERMEDIATE_DIR="${LOCAL_OUTDIR}/intermediate/annotate_variants"' in launcher
    assert 'VEP_JSON_OUTPUT="/tmp/vep_annotations.json"' in launcher
    assert '--analysis_manifest_dir "${ANALYSIS_MANIFEST_DIR}" \\' in launcher
    assert '--vep_json_output "${VEP_JSON_OUTPUT}" \\' in launcher
    assert '--coding_input' not in launcher
    assert '--clinvar_input' not in launcher
    assert '--coding_json_output' not in launcher
    assert '--clinvar_json_output' not in launcher
    assert '--vep_parallel "${ANNOTATION_VEP_PARALLEL}" \\' in launcher
    assert '--vep_forks "${ANNOTATION_VEP_FORKS}" \\' in launcher
    assert '--vep_buffer_size "${ANNOTATION_VEP_BUFFER_SIZE}" \\' in launcher
    assert 'ANNOTATION_VEP_PARALLEL="$(ANNOTATION_VEP_PARALLEL)"' in makefile
    assert 'ANNOTATION_VEP_FORKS="$(ANNOTATION_VEP_FORKS)"' in makefile
    assert 'ANNOTATION_VEP_BUFFER_SIZE="$(ANNOTATION_VEP_BUFFER_SIZE)"' in makefile
    assert 'SAVE_INTERMEDIATE_FILES="$(SAVE_INTERMEDIATE_FILES)" bash scripts/run_annotation_local.sh' in makefile


def test_gcloud_annotation_launcher_waits_for_status_without_vm_retry():
    launcher = _read(PROJECT_DIR / "scripts" / "run_annotation_vm.sh")
    makefile = _read(PROJECT_DIR / "Makefile")

    assert 'poll_status()' in launcher
    assert 'status.txt' in launcher
    assert 'annotation.log' in launcher
    assert 'VM_SUFFIX="${RUN_TRACE_SUFFIX//_/-}"' in launcher
    assert 'VM_NAME="${ANNOTATION_VM_NAME}-${VM_SUFFIX}"' in launcher
    assert 'STATUS_URI="${STATUS_ROOT}/status.txt"' in launcher
    assert 'PROGRESS_URI="${STATUS_ROOT}/progress.txt"' in launcher
    assert 'LOG_URI="${STATUS_ROOT}/annotation.log"' in launcher
    assert 'attempt${attempt}' not in launcher
    assert 'launch_annotation_attempt' not in launcher
    assert 'update_progress()' in launcher
    assert 'Phase 1/4: Installing system packages' in launcher
    assert 'Phase 4/4: Running annotation container' in launcher
    assert 'ANALYSIS_MANIFEST_DIR="vidra_analysis_ready_manifest"' in launcher
    assert 'VEP_JSON_OUTPUT="/tmp/vep_annotations.json"' in launcher
    assert 'ANNOTATION_VEP_PARALLEL' in launcher
    assert 'ANNOTATION_VEP_FORKS' in launcher
    assert 'ANNOTATION_VEP_BUFFER_SIZE' in launcher
    assert '--analysis_manifest_dir "${ANALYSIS_MANIFEST_DIR}" \\' in launcher
    assert '--vep_json_output "${VEP_JSON_OUTPUT}" \\' in launcher
    assert '--coding_input' not in launcher
    assert '--clinvar_input' not in launcher
    assert '--coding_json_output' not in launcher
    assert '--clinvar_json_output' not in launcher
    assert '--vep_parallel "${ANNOTATION_VEP_PARALLEL}" \\' in launcher
    assert '--vep_forks "${ANNOTATION_VEP_FORKS}" \\' in launcher
    assert '--vep_buffer_size "${ANNOTATION_VEP_BUFFER_SIZE}" \\' in launcher
    assert 'SAVE_INTERMEDIATE_FILES' not in launcher
    assert 'Annotation VM completed successfully; log: ${LOG_URI}' in launcher
    assert '--boot-disk-type "${ANNOTATION_VM_BOOT_DISK_TYPE}" \\' in launcher
    assert '--image-family "${ANNOTATION_VM_IMAGE_FAMILY}" \\' in launcher
    assert '--image-project "${ANNOTATION_VM_IMAGE_PROJECT}" \\' in launcher
    assert 'ANNOTATION_VM_BOOT_DISK_TYPE="$(ANNOTATION_VM_BOOT_DISK_TYPE)"' in makefile
    assert 'ANNOTATION_VM_IMAGE_FAMILY="$(ANNOTATION_VM_IMAGE_FAMILY)"' in makefile
    assert 'ANNOTATION_VM_IMAGE_PROJECT="$(ANNOTATION_VM_IMAGE_PROJECT)"' in makefile
    assert 'ANNOTATION_VM_TIMEOUT_SECONDS="$(ANNOTATION_VM_TIMEOUT_SECONDS)"' in makefile
    assert 'ANNOTATION_VEP_PARALLEL="$(ANNOTATION_VEP_PARALLEL)"' in makefile
    assert 'ANNOTATION_VEP_FORKS="$(ANNOTATION_VEP_FORKS)"' in makefile
    assert 'ANNOTATION_VEP_BUFFER_SIZE="$(ANNOTATION_VEP_BUFFER_SIZE)"' in makefile
    assert 'SAVE_INTERMEDIATE_FILES="$(SAVE_INTERMEDIATE_FILES)" bash scripts/run_annotation_vm.sh' not in makefile
    assert 'ANNOTATION_VM_RETRY_MACHINE_TYPE' not in makefile


def test_github_workflow_builds_prepare_annotation_and_bayes_images():
    workflow = _read(PROJECT_DIR / ".github" / "workflows" / "build-image.yaml")

    assert "Build Step 1 spark image" in workflow
    assert "Build annotation image" in workflow
    assert "Build bayes image" in workflow
    assert "./docker/spark-dataproc.Dockerfile" in workflow
    assert "./docker/annotation.Dockerfile" in workflow
    assert "./docker/bayes.Dockerfile" in workflow
    assert "spark-batch.Dockerfile" not in workflow


def test_prepare_image_is_trimmed_to_step1_runtime():
    dockerfile = _read(PROJECT_DIR / "docker" / "spark-dataproc.Dockerfile")
    requirements = _read(PROJECT_DIR / "docker" / "prepare_requirements.txt")

    assert "external Step 1 only" in dockerfile
    assert "prepare_requirements.txt" in dockerfile
    assert "python_requirements.txt" not in dockerfile
    assert "CMDSTAN" not in dockerfile
    assert "cmdstanpy" not in dockerfile
    assert "stan_models" not in dockerfile
    assert requirements.strip() == "gcsfs"
