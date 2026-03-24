FROM python:3.12-slim

WORKDIR /app

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH=/app
ENV PYSPARK_PYTHON=python3
ENV SPARK_LOCAL_IP=127.0.0.1
ENV SPARK_LOCAL_HOSTNAME=localhost

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    default-jre-headless \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir pyspark

COPY tools /app/tools

CMD ["bash"]
