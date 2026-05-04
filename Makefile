PROJECT_DIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

include conf/make.defaults.mk
-include conf/make.local.mk

.PHONY: \
	clean \
	clean_gcloud_artifacts \
	list_gcloud_artifacts \
	local_convert_az_phewas_parquet \
	local_convert_az_variants_bz2 \
	local_all_protein_coding_ensg \
	local_sync_ot_release_inputs \
	docker_build \
	docker_build_prepare \
	docker_build_annotation \
	docker_build_bayes \
	docker_tag_gcloud \
	docker_push_gcloud \
	gcloud_upload_vep_cache \
	gcloud_upload_testdata \
	venv \
	test \
	local_step1 \
	local_step2 \
	gcloud_step1 \
	gcloud_step2 \
	nextflow_local_step3 \
	nextflow_test_step3 \
	nextflow_gcloud_test_step3 \
	nextflow_gcloud_prod_step3 \
	sync_gcloud_prod_results \
	config_test

clean:
	rm -rf .nextflow* work .pytest_cache
	find . -type d -name "__pycache__" -exec rm -rf {} +

clean_gcloud_artifacts:
	@echo "Mode: $(MODE)"
	@echo "Target image: $(IMAGE)"
	@read -p "Type 'DELETE' to confirm: " confirm && [ "$$confirm" = "DELETE" ] || exit 1
	@if [ "$(MODE)" = "untagged" ]; then \
		gcloud artifacts docker images list $(IMAGE) --include-tags | \
		awk '/sha256:/ && $$0 !~ / latest( |$$)/ {print $$2}' | \
		while read digest; do \
			[ -n "$$digest" ] || continue; \
			echo "Deleting untagged $$digest"; \
			gcloud artifacts docker images delete "$(IMAGE)@$$digest" --quiet; \
		done; \
	elif [ "$(MODE)" = "all" ]; then \
		echo "Deleting image path $(IMAGE) and all of its tags/digests"; \
		gcloud artifacts docker images delete "$(IMAGE)" --delete-tags --quiet; \
	else \
		echo "Invalid MODE='$(MODE)'. Use MODE=untagged or MODE=all"; \
		exit 1; \
	fi

list_gcloud_artifacts:
	gcloud artifacts docker images list $(IMAGE) --include-tags --format="table(package,version,tags)"

docker_build: docker_build_prepare docker_build_annotation docker_build_bayes

docker_build_prepare:
	docker build --no-cache --platform=linux/amd64 -t $(PREPARE_IMAGE_LOCAL) -f docker/spark-dataproc.Dockerfile .

docker_build_annotation:
	docker build --no-cache --platform=linux/amd64 -t $(ANNOTATION_IMAGE_LOCAL) -f docker/annotation.Dockerfile .

docker_build_bayes:
	docker build --no-cache --platform=linux/amd64 -t $(BAYES_IMAGE_LOCAL) -f docker/bayes.Dockerfile .

docker_tag_gcloud:
	docker tag $(PREPARE_IMAGE_LOCAL) $(PREPARE_IMAGE_GCLOUD)
	docker tag $(ANNOTATION_IMAGE_LOCAL) $(ANNOTATION_IMAGE_GCLOUD)
	docker tag $(BAYES_IMAGE_LOCAL) $(BAYES_IMAGE_GCLOUD)

docker_push_gcloud: docker_tag_gcloud
	docker push $(PREPARE_IMAGE_GCLOUD)
	docker push $(ANNOTATION_IMAGE_GCLOUD)
	docker push $(BAYES_IMAGE_GCLOUD)

gcloud_upload_vep_cache:
	@set -eu; \
	tmp_dir="$$(mktemp -d /tmp/vidra5_vep_cache_XXXXXX)"; \
	trap 'rm -rf "$$tmp_dir"' EXIT; \
	cache_archive="$$tmp_dir/homo_sapiens_vep_115_GRCh38.tar.gz"; \
	curl -L "$(VEP_CACHE_SOURCE_URL)" -o "$$cache_archive"; \
	gcloud storage cp "$$cache_archive" "$(VEP_CACHE_GCLOUD_URI)"

gcloud_upload_testdata:
	gcloud storage cp "$(TEST_GENES_FILE)" "$(GCLOUD_TESTDATA_DIR)/genes_ensembl_small.csv"

venv:
	python3 -m venv .venv
	./.venv/bin/pip install --upgrade pip
	./.venv/bin/pip install -r python_requirements.txt pytest duckdb

test:
	./.venv/bin/pytest

local_sync_ot_release_inputs:
	@set -eu; \
	mkdir -p "$(LOCAL_RELEASE_ROOT)"; \
	for folder in $(OT_RELEASE_RSYNC_FOLDERS); do \
		echo "Syncing $$folder into $(LOCAL_RELEASE_ROOT)"; \
		rsync -rpltvz "$(OT_RELEASE_RSYNC_BASE)/$$folder" "$(LOCAL_RELEASE_ROOT)/"; \
	done

local_convert_az_phewas_parquet:
	$(PROJECT_DIR)/.venv/bin/python $(PROJECT_DIR)/scripts/convert_az_phewas_parquet.py \
		--dataset all \
		--binary-input "$(LOCAL_AZ_BURDEN_BINARY_XZ_FILE)" \
		--binary-output-dir "$(LOCAL_AZ_BURDEN_BINARY_DIR)" \
		--quantitative-input "$(LOCAL_AZ_BURDEN_QUANTITATIVE_XZ_FILE)" \
		--quantitative-output-dir "$(LOCAL_AZ_BURDEN_QUANTITATIVE_DIR)"

local_convert_az_variants_bz2:
	@if [ -f "$(LOCAL_AZ_VARIANTS_FILE)" ]; then \
		echo "$(LOCAL_AZ_VARIANTS_FILE) already exists"; \
	else \
		xz -dc "$(LOCAL_AZ_VARIANTS_XZ_FILE)" | bzip2 -c > "$(LOCAL_AZ_VARIANTS_FILE)"; \
	fi

local_all_protein_coding_ensg:
	$(PROJECT_DIR)/.venv/bin/python $(PROJECT_DIR)/tools/export_protein_coding_ensg.py \
		--target-parquet "$(LOCAL_TARGET_DATA_DIR)/*.parquet" \
		--output "$(LOCAL_PROTEIN_CODING_ENSG_FILE)"

local_step1:
	PROJECT_DIR="$(PROJECT_DIR)" HOST_UID="$(HOST_UID)" HOST_GID="$(HOST_GID)" PREPARE_IMAGE_LOCAL="$(PREPARE_IMAGE_LOCAL)" LOCAL_GENES="$(LOCAL_GENES)" LOCAL_OUTDIR="$(LOCAL_OUTDIR)" LOCAL_RELEASE_ROOT="$(LOCAL_RELEASE_ROOT)" LOCAL_VIDRA_DATA_ROOT="$(LOCAL_VIDRA_DATA_ROOT)" LOCAL_AZ_VARIANTS_FILE="$(LOCAL_AZ_VARIANTS_FILE)" LOCAL_LEAD_VARIANT_EFFECT_DIR="$(LOCAL_LEAD_VARIANT_EFFECT_DIR)" LOCAL_AZ_BURDEN_BINARY_DIR="$(LOCAL_AZ_BURDEN_BINARY_DIR)" LOCAL_AZ_BURDEN_QUANTITATIVE_DIR="$(LOCAL_AZ_BURDEN_QUANTITATIVE_DIR)" LOCAL_TARGET_DATA_DIR="$(LOCAL_TARGET_DATA_DIR)" COLOCALISATION_THRESHOLD="$(COLOCALISATION_THRESHOLD)" LOCAL_PREPARE_SPARK_MASTER="$(LOCAL_PREPARE_SPARK_MASTER)" LOCAL_PREPARE_SPARK_DRIVER_MEMORY="$(LOCAL_PREPARE_SPARK_DRIVER_MEMORY)" LOCAL_PREPARE_SPARK_SHUFFLE_PARTITIONS="$(LOCAL_PREPARE_SPARK_SHUFFLE_PARTITIONS)" LOCAL_PREPARE_SPARK_DEFAULT_PARALLELISM="$(LOCAL_PREPARE_SPARK_DEFAULT_PARALLELISM)" bash scripts/run_prepare_local.sh

local_step2:
	PROJECT_DIR="$(PROJECT_DIR)" HOST_UID="$(HOST_UID)" HOST_GID="$(HOST_GID)" ANNOTATION_IMAGE_LOCAL="$(ANNOTATION_IMAGE_LOCAL)" LOCAL_OUTDIR="$(LOCAL_OUTDIR)" LOCAL_VEP_RESOURCES_DIR="$(LOCAL_VEP_RESOURCES_DIR)" LOCAL_FOLDX_FILE="$(LOCAL_FOLDX_FILE)" LOCAL_REQUIRED_VEP_PLUGINS="$(LOCAL_REQUIRED_VEP_PLUGINS)" ANNOTATION_VEP_PARALLEL="$(ANNOTATION_VEP_PARALLEL)" ANNOTATION_VEP_FORKS="$(ANNOTATION_VEP_FORKS)" ANNOTATION_VEP_BUFFER_SIZE="$(ANNOTATION_VEP_BUFFER_SIZE)" SAVE_INTERMEDIATE_FILES="$(SAVE_INTERMEDIATE_FILES)" bash scripts/run_annotation_local.sh

gcloud_step1:
	PROJECT_DIR="$(PROJECT_DIR)" PREPARE_IMAGE_GCLOUD="$(PREPARE_IMAGE_GCLOUD)" GCP_PROJECT="$(GCP_PROJECT)" GCP_REGION="$(GCP_REGION)" GCP_SERVICE_ACCOUNT_EMAIL="$(GCP_SERVICE_ACCOUNT_EMAIL)" DATAPROC_DEPS_BUCKET="$(DATAPROC_DEPS_BUCKET)" DATAPROC_BATCH_PREFIX="$(DATAPROC_BATCH_PREFIX)" DATAPROC_STEP1_TTL="$(DATAPROC_STEP1_TTL)" DATAPROC_STEP1_PROPERTIES="$(DATAPROC_STEP1_PROPERTIES)" GCLOUD_OUTDIR="$(GCLOUD_OUTDIR)" GCLOUD_GENES="$(GCLOUD_GENES)" GCLOUD_RELEASE_ROOT="$(GCLOUD_RELEASE_ROOT)" GCLOUD_VIDRA_DATA_ROOT="$(GCLOUD_VIDRA_DATA_ROOT)" GCLOUD_AZ_VARIANTS_FILE="$(GCLOUD_AZ_VARIANTS_FILE)" GCLOUD_LEAD_VARIANT_EFFECT_DIR="$(GCLOUD_LEAD_VARIANT_EFFECT_DIR)" GCLOUD_AZ_BURDEN_BINARY_DIR="$(GCLOUD_AZ_BURDEN_BINARY_DIR)" GCLOUD_AZ_BURDEN_QUANTITATIVE_DIR="$(GCLOUD_AZ_BURDEN_QUANTITATIVE_DIR)" GCLOUD_TARGET_DATA_DIR="$(GCLOUD_TARGET_DATA_DIR)" COLOCALISATION_THRESHOLD="$(COLOCALISATION_THRESHOLD)" bash scripts/run_prepare_dataproc.sh

gcloud_step2:
	PROJECT_DIR="$(PROJECT_DIR)" GCP_PROJECT="$(GCP_PROJECT)" GCP_ZONE="$(GCP_ZONE)" GCP_SERVICE_ACCOUNT_EMAIL="$(GCP_SERVICE_ACCOUNT_EMAIL)" ANNOTATION_IMAGE_GCLOUD="$(ANNOTATION_IMAGE_GCLOUD)" GCLOUD_OUTDIR="$(GCLOUD_OUTDIR)" GCLOUD_VEP_CACHE_ARCHIVE="$(GCLOUD_VEP_CACHE_ARCHIVE)" GCLOUD_VEP_PLUGIN_DATA_DIR="$(GCLOUD_VEP_PLUGIN_DATA_DIR)" GCLOUD_FOLDX_FILE="$(GCLOUD_FOLDX_FILE)" GCLOUD_REQUIRED_VEP_PLUGINS="$(GCLOUD_REQUIRED_VEP_PLUGINS)" ANNOTATION_VM_NAME="$(ANNOTATION_VM_NAME)" ANNOTATION_VM_MACHINE_TYPE="$(ANNOTATION_VM_MACHINE_TYPE)" ANNOTATION_VM_BOOT_DISK_SIZE="$(ANNOTATION_VM_BOOT_DISK_SIZE)" ANNOTATION_VM_BOOT_DISK_TYPE="$(ANNOTATION_VM_BOOT_DISK_TYPE)" ANNOTATION_VM_IMAGE_FAMILY="$(ANNOTATION_VM_IMAGE_FAMILY)" ANNOTATION_VM_IMAGE_PROJECT="$(ANNOTATION_VM_IMAGE_PROJECT)" ANNOTATION_VM_TIMEOUT_SECONDS="$(ANNOTATION_VM_TIMEOUT_SECONDS)" ANNOTATION_VM_AUTO_DELETE="$(ANNOTATION_VM_AUTO_DELETE)" ANNOTATION_VEP_PARALLEL="$(ANNOTATION_VEP_PARALLEL)" ANNOTATION_VEP_FORKS="$(ANNOTATION_VEP_FORKS)" ANNOTATION_VEP_BUFFER_SIZE="$(ANNOTATION_VEP_BUFFER_SIZE)" bash scripts/run_annotation_vm.sh

nextflow_local_step3:
	HOST_UID="$(HOST_UID)" HOST_GID="$(HOST_GID)" nextflow run main.nf -profile local -resume --outdir "$(LOCAL_OUTDIR)" --bucket_uri "$(LOCAL_OUTDIR)"

nextflow_test_step3:
	HOST_UID="$(HOST_UID)" HOST_GID="$(HOST_GID)" nextflow run main.nf -profile test -resume

nextflow_gcloud_test_step3:
	nextflow run main.nf -profile gcloud_test -resume

nextflow_gcloud_prod_step3:
	nextflow run main.nf -profile gcloud_prod -resume

config_test:
	nextflow config -flat -profile local
	nextflow config -flat -profile test
	nextflow config -flat -profile gcloud_test
	nextflow config -flat -profile gcloud_prod

sync_gcloud_prod_results:
	@set -eu; \
	mkdir -p "$(LOCAL_GCLOUD_BAYES_RESULTS_DIR)"; \
	echo "Checking remote results: $(GCLOUD_BAYES_RESULTS_DIR)/"; \
	gcloud storage ls "$(GCLOUD_BAYES_RESULTS_DIR)/" >/dev/null; \
	echo "Syncing new or changed result files to $(LOCAL_GCLOUD_BAYES_RESULTS_DIR)"; \
	gcloud storage rsync -r "$(GCLOUD_BAYES_RESULTS_DIR)" "$(LOCAL_GCLOUD_BAYES_RESULTS_DIR)"
