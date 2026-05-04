FROM python:3.11-slim

USER root

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONPATH=/app
ENV CMDSTAN_INSTALL_PATH=/opt/cmdstan

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    build-essential \
    libcurl4-openssl-dev \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY python_requirements.txt /tmp/python_requirements.txt
RUN pip install --no-cache-dir -r /tmp/python_requirements.txt

RUN python3 -c "import cmdstanpy; cmdstanpy.install_cmdstan(dir='/opt/cmdstan', version='2.34.1')"
ENV CMDSTAN=/opt/cmdstan/cmdstan-2.34.1

COPY tools /app/tools
COPY stan_models /opt/vidra/stan_models

WORKDIR /opt/vidra/stan_models
RUN python3 -c "from cmdstanpy import CmdStanModel; CmdStanModel(stan_file='VIDRA.stan').compile(); CmdStanModel(stan_file='VIDRA_single_variant.stan').compile()"

WORKDIR /app

CMD ["bash"]
