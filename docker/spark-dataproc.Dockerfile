# Dataproc runtime image for external Step 1 only.
# Keep the final user as `spark`, matching Dataproc Serverless expectations.
FROM --platform=linux/amd64 europe-west1-docker.pkg.dev/cloud-dataproc/spark/dataproc_2.2:latest

USER root

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH=/app

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY docker/prepare_requirements.txt /tmp/prepare_requirements.txt
RUN pip install --no-cache-dir -r /tmp/prepare_requirements.txt

COPY tools /app/tools
WORKDIR /app

USER spark

CMD ["bash"]
