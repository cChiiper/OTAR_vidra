FROM ensemblorg/ensembl-vep:release_115.0

USER root

ARG VEP_PLUGINS_RELEASE=release/115
ARG VEP_PLUGINS_BASE_URL=https://raw.githubusercontent.com/Ensembl/VEP_plugins/${VEP_PLUGINS_RELEASE}

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH=/app
ENV VEP_PLUGIN_RESOURCES_DIR=/plugin_resources

RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    ca-certificates \
    curl \
    procps \
    python3 \
    python3-pip \
    python3-venv \
    samtools \
    tabix \
    && rm -rf /var/lib/apt/lists/*

COPY docker/annotation_requirements.txt /tmp/annotation_requirements.txt
RUN python3 -m pip install --no-cache-dir -r /tmp/annotation_requirements.txt

# Install the core VEP Perl plugin modules required by the VIDRA models. The missing ones will be ignored.
RUN set -eux; \
    plugins_dir="/opt/vep/src/ensembl-vep/Plugins"; \
    mkdir -p "${plugins_dir}"; \
    curl -fsSL "${VEP_PLUGINS_BASE_URL}/AlphaMissense.pm" -o "${plugins_dir}/AlphaMissense.pm"; \
    curl -fsSL "${VEP_PLUGINS_BASE_URL}/CADD.pm" -o "${plugins_dir}/CADD.pm"; \
    curl -fsSL "${VEP_PLUGINS_BASE_URL}/REVEL.pm" -o "${plugins_dir}/REVEL.pm"

COPY tools /app/tools

WORKDIR /app

CMD ["bash"]
