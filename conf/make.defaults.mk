# Shared Make defaults for launch targets only.
# This file is Make-specific on purpose and is not consumed by Nextflow.

# Host user mapping for local Docker runs and local Nextflow Step 3.
HOST_UID ?= $(shell id -u)
HOST_GID ?= $(shell id -g)

# Images built locally.
PREPARE_IMAGE_LOCAL ?= vidra5-spark-dataproc:latest
ANNOTATION_IMAGE_LOCAL ?= vidra5-annotation:latest
BAYES_IMAGE_LOCAL ?= vidra5-bayes:latest

# Images pushed to Artifact Registry.
PREPARE_IMAGE_GCLOUD ?= europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra5-spark-dataproc:latest
ANNOTATION_IMAGE_GCLOUD ?= europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra5-annotation:latest
BAYES_IMAGE_GCLOUD ?= europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra5-bayes:latest

# Small checked-in test gene list mirrored into GCS for smoke tests.
GCLOUD_TESTDATA_DIR ?= gs://vidra-2-0/nextflow/testdata
TEST_GENES_FILE ?= testdata/genes_ensembl_small.csv

# Artifact Registry cleanup target.
MODE ?= untagged
IMAGE ?= $(ANNOTATION_IMAGE_GCLOUD)

# Shared local Step 1 / Step 2 defaults.
LOCAL_GENES ?= $(PROJECT_DIR)/testdata/genes_ensembl_small.csv
LOCAL_OUTDIR ?= $(PROJECT_DIR)/results/test
LOCAL_RELEASE_ROOT ?= 
OT_RELEASE_RSYNC_BASE ?= rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output
OT_RELEASE_RSYNC_FOLDERS ?= colocalisation credible_set study variant evidence_gene_burden evidence_eva target
LOCAL_VIDRA_DATA_ROOT ?= 
LOCAL_AZ_PHEWAS_INPUT_ROOT ?= $(LOCAL_VIDRA_DATA_ROOT)/TEMP
LOCAL_AZ_VARIANTS_XZ_FILE ?= $(LOCAL_VIDRA_DATA_ROOT)/azphewas-com-470k-exwas-binary.csv.xz
LOCAL_AZ_VARIANTS_FILE ?= $(LOCAL_VIDRA_DATA_ROOT)/azphewas-com-470k-exwas-binary.csv.bz2
LOCAL_LEAD_VARIANT_EFFECT_DIR ?= 
LOCAL_AZ_BURDEN_BINARY_XZ_FILE ?= $(LOCAL_AZ_PHEWAS_INPUT_ROOT)/azphewas-com-470k-phewas-binary.csv.xz
LOCAL_AZ_BURDEN_QUANTITATIVE_XZ_FILE ?= $(LOCAL_AZ_PHEWAS_INPUT_ROOT)/azphewas-com-470k-phewas-quantitative.csv.xz
LOCAL_AZ_BURDEN_BINARY_DIR ?= $(LOCAL_VIDRA_DATA_ROOT)/azphewas-com-470k-phewas-binary
LOCAL_AZ_BURDEN_QUANTITATIVE_DIR ?= $(LOCAL_VIDRA_DATA_ROOT)/azphewas-com-470k-phewas-quantitative
LOCAL_TARGET_DATA_DIR ?= $(LOCAL_RELEASE_ROOT)/target
LOCAL_PROTEIN_CODING_ENSG_FILE ?= $(LOCAL_VIDRA_DATA_ROOT)/all_protein_coding_ENSG.csv
LOCAL_VEP_RESOURCES_DIR ?= 
LOCAL_REQUIRED_VEP_PLUGINS ?= AlphaMissense
ANNOTATION_VEP_PARALLEL ?= 4
ANNOTATION_VEP_FORKS ?= 2
ANNOTATION_VEP_BUFFER_SIZE ?= 5000
LOCAL_PREPARE_SPARK_MASTER ?= local[2]
LOCAL_PREPARE_SPARK_DRIVER_MEMORY ?= 4g
LOCAL_PREPARE_SPARK_SHUFFLE_PARTITIONS ?= 16
LOCAL_PREPARE_SPARK_DEFAULT_PARALLELISM ?= 4
COLOCALISATION_THRESHOLD ?= 0.7
SAVE_INTERMEDIATE_FILES ?= false

# Shared cloud Step 1 / Step 2 defaults.
GCP_PROJECT ?= open-targets-eu-dev
GCP_REGION ?= europe-west1
GCP_ZONE ?= europe-west1-b
GCP_SERVICE_ACCOUNT_EMAIL ?= 
DATAPROC_DEPS_BUCKET ?= gs://vidra-2-0
DATAPROC_BATCH_PREFIX ?= vidra5
DATAPROC_STEP1_TTL ?= 14400s
DATAPROC_STEP1_PROPERTIES ?= spark.sql.execution.arrow.pyspark.enabled=true
GCLOUD_OUTDIR ?= gs://vidra-2-0/nextflow/results/gcloud_prod
GCLOUD_BAYES_RESULTS_DIR ?= $(GCLOUD_OUTDIR)/vidra_results
LOCAL_GCLOUD_PROD_OUTDIR ?= $(PROJECT_DIR)/results/gcloud_prod
LOCAL_GCLOUD_BAYES_RESULTS_DIR ?= $(LOCAL_GCLOUD_PROD_OUTDIR)/vidra_results
GCLOUD_GENES ?= gs://vidra-2-0/nextflow/reference/all_protein_coding_ENSG.csv
GCLOUD_RELEASE_ROOT ?= gs://open-targets-data-releases/26.03/output
GCLOUD_VIDRA_DATA_ROOT ?= gs://vidra-2-0/raw_data
GCLOUD_AZ_VARIANTS_FILE ?= $(GCLOUD_VIDRA_DATA_ROOT)/azphewas-com-470k-exwas-binary.csv.bz2
GCLOUD_LEAD_VARIANT_EFFECT_DIR ?= gs://ot_orchestration/gentropy_manuscript/data/26.03/intermediate/lead_variant_effect
GCLOUD_AZ_BURDEN_BINARY_DIR ?= gs://otar000-evidence_input/GeneBurden/data_files/azphewas-com-470k-phewas-binary
GCLOUD_AZ_BURDEN_QUANTITATIVE_DIR ?= gs://otar000-evidence_input/GeneBurden/data_files/azphewas-com-470k-phewas-quantitative
GCLOUD_TARGET_DATA_DIR ?= $(GCLOUD_RELEASE_ROOT)/target
GCLOUD_VEP_CACHE_ARCHIVE ?= gs://vidra-2-0/nextflow/reference/vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
GCLOUD_VEP_PLUGIN_DATA_DIR ?= gs://vidra-2-0/raw_data
GCLOUD_FOLDX_FILE ?= gs://vidra-2-0/raw_data/foldx_energy.csv.gz
GCLOUD_REQUIRED_VEP_PLUGINS ?= AlphaMissense,CADD,REVEL
ANNOTATION_VM_NAME ?= vidra5-annotation-vm
ANNOTATION_VM_MACHINE_TYPE ?= e2-standard-8
ANNOTATION_VM_BOOT_DISK_SIZE ?= 250GB
ANNOTATION_VM_TIMEOUT_SECONDS ?= 18000
ANNOTATION_VM_BOOT_DISK_TYPE ?= pd-ssd
ANNOTATION_VM_IMAGE_FAMILY ?= debian-12
ANNOTATION_VM_IMAGE_PROJECT ?= debian-cloud
ANNOTATION_VM_AUTO_DELETE ?= true

VEP_CACHE_SOURCE_URL ?= https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
VEP_CACHE_GCLOUD_URI ?= gs://vidra-2-0/nextflow/reference/vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
