#!/usr/bin/env python3
"""Submit a VIDRA PySpark step to Dataproc Serverless and record metadata."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
from pathlib import Path


def sanitize_batch_component(value: str) -> str:
    text = re.sub(r"[^a-z0-9-]+", "-", str(value).strip().lower())
    text = re.sub(r"-{2,}", "-", text).strip("-")
    return text or "vidra"


def build_batch_id(batch_prefix: str, step_name: str, trace_suffix: str) -> str:
    parts = [
        sanitize_batch_component(batch_prefix),
        sanitize_batch_component(step_name),
        sanitize_batch_component(trace_suffix),
    ]
    batch_id = "-".join(part for part in parts if part)
    return batch_id[:63].rstrip("-")


def build_submit_command(
    *,
    script_path: str,
    batch_id: str,
    project: str,
    region: str,
    deps_bucket: str,
    container_image: str,
    service_account: str = "",
    ttl: str = "",
    properties: str = "",
    script_args: list[str] | None = None,
) -> list[str]:
    command = [
        "gcloud",
        "dataproc",
        "batches",
        "submit",
        "pyspark",
        script_path,
        f"--project={project}",
        f"--region={region}",
        f"--deps-bucket={deps_bucket}",
        f"--container-image={container_image}",
        f"--batch={batch_id}",
    ]
    if service_account:
        command.append(f"--service-account={service_account}")
    if ttl:
        command.append(f"--ttl={ttl}")
    if properties:
        command.append(f"--properties={properties}")
    if script_args:
        command.append("--")
        command.extend(script_args)
    return command


def run_json_command(command: list[str]) -> dict:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    stdout = (completed.stdout or "").strip()
    return json.loads(stdout) if stdout else {}


def write_json(path: str, payload: dict) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Submit a VIDRA step to Dataproc Serverless.")
    parser.add_argument("--step_name", required=True)
    parser.add_argument("--script_path", required=True)
    parser.add_argument("--project", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--deps_bucket", required=True)
    parser.add_argument("--container_image", required=True)
    parser.add_argument("--service_account", default="")
    parser.add_argument("--batch_prefix", required=True)
    parser.add_argument("--trace_suffix", required=True)
    parser.add_argument("--metadata_output", required=True)
    parser.add_argument("--ttl", default="")
    parser.add_argument("--properties", default="")
    parser.add_argument("--script_arg", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    batch_id = build_batch_id(args.batch_prefix, args.step_name, args.trace_suffix)
    submit_cmd = build_submit_command(
        script_path=args.script_path,
        batch_id=batch_id,
        project=args.project,
        region=args.region,
        deps_bucket=args.deps_bucket,
        container_image=args.container_image,
        service_account=args.service_account,
        ttl=args.ttl,
        properties=args.properties,
        script_args=args.script_arg,
    )

    console_url = (
        f"https://console.cloud.google.com/dataproc/batches/{args.region}/{batch_id}?project={args.project}"
    )
    describe_cmd = [
        "gcloud",
        "dataproc",
        "batches",
        "describe",
        batch_id,
        f"--project={args.project}",
        f"--region={args.region}",
        "--format=json",
    ]
    wait_cmd = [
        "gcloud",
        "dataproc",
        "batches",
        "wait",
        batch_id,
        f"--project={args.project}",
        f"--region={args.region}",
    ]

    print(f"[dataproc_batch_submit] batch_id={batch_id}", flush=True)
    print(f"[dataproc_batch_submit] console={console_url}", flush=True)
    print(
        "[dataproc_batch_submit] describe=" + " ".join(shlex.quote(part) for part in describe_cmd),
        flush=True,
    )
    print(
        "[dataproc_batch_submit] wait=" + " ".join(shlex.quote(part) for part in wait_cmd),
        flush=True,
    )

    submit_error: dict | None = None
    try:
        subprocess.run(submit_cmd, check=True)
    except subprocess.CalledProcessError as exc:
        submit_error = {
            "returncode": exc.returncode,
            "command": exc.cmd,
        }

    describe_payload = {}
    describe_error = None
    try:
        describe_payload = run_json_command(describe_cmd)
    except subprocess.CalledProcessError as exc:
        describe_error = {
            "returncode": exc.returncode,
            "command": exc.cmd,
        }

    write_json(
        args.metadata_output,
        {
            "step_name": args.step_name,
            "batch_id": batch_id,
            "console_url": console_url,
            "submit_command": submit_cmd,
            "describe_command": describe_cmd,
            "wait_command": wait_cmd,
            "submit_error": submit_error,
            "describe_error": describe_error,
            "describe_payload": describe_payload,
            "runtime_output_uri": describe_payload.get("runtimeInfo", {}).get("outputUri", ""),
            "approximate_usage": describe_payload.get("runtimeInfo", {}).get("approximateUsage", {}),
        },
    )

    if submit_error is not None:
        raise subprocess.CalledProcessError(
            submit_error["returncode"],
            submit_error["command"],
        )


if __name__ == "__main__":
    main()
