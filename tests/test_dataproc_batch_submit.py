import json

from tools import dataproc_batch_submit as submit


def test_build_batch_id_sanitizes_components():
    batch_id = submit.build_batch_id("VIDRA5", "Run Bayesian Analysis", "dataproc_prod_20260331_101010_001")

    assert batch_id.startswith("vidra5-run-bayesian-analysis-dataproc-prod-20260331-101010-001")
    assert len(batch_id) <= 63


def test_build_submit_command_matches_expected_vidra_shape():
    command = submit.build_submit_command(
        script_path="/repo/tools/run_bayesian_analysis.py",
        batch_id="vidra5-step3",
        project="open-targets-eu-dev",
        region="europe-west1",
        deps_bucket="gs://vidra-2-0",
        container_image="europe-west1-docker.pkg.dev/open-targets-genetics-dev/opentargets/vidra5-spark-dataproc:latest",
        ttl="86400s",
        properties="spark.sql.execution.arrow.pyspark.enabled=true,spark.executor.instances=20,spark.dynamicAllocation.maxExecutors=50",
        script_args=["--bucket_uri", "gs://vidra-2-0/nextflow/results/dataproc_prod"],
    )

    assert command[:6] == [
        "gcloud",
        "dataproc",
        "batches",
        "submit",
        "pyspark",
        "/repo/tools/run_bayesian_analysis.py",
    ]
    assert "--wait" not in command
    assert "--ttl=86400s" in command
    assert "--properties=spark.sql.execution.arrow.pyspark.enabled=true,spark.executor.instances=20,spark.dynamicAllocation.maxExecutors=50" in command
    assert command[-3:] == ["--", "--bucket_uri", "gs://vidra-2-0/nextflow/results/dataproc_prod"]


def test_write_json_creates_pretty_json_file(tmp_path):
    output = tmp_path / "metadata.json"
    submit.write_json(str(output), {"batch_id": "vidra5-step1"})

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload == {"batch_id": "vidra5-step1"}
