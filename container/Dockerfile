FROM debian:bookworm-slim

LABEL org.opencontainers.image.authors="Johanna Girodolle"
LABEL org.opencontainers.image.version="0.4.0"
LABEL org.opencontainers.image.title="GeneModelTransfer dependencies"
LABEL org.opencontainers.image.description="Runtime dependencies for GeneModelTransfer"
LABEL org.opencontainers.image.source="https://github.com/ranwez/GeneModelTransfer"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
    && rm -rf /var/lib/apt/lists/*

# Install Pixi
RUN curl -fsSL https://pixi.sh/install.sh | bash
ENV PATH="/root/.pixi/bin:${PATH}"

# Install the locked dependency environment
WORKDIR /opt/pixi
COPY pixi.toml pixi.lock ./
RUN pixi install --locked

# Make the environment available without requiring `pixi run`
ENV PATH="/opt/pixi/.pixi/envs/default/bin:${PATH}"

# Basic runtime checks
RUN exonerate --version 2>&1 | grep -F "version 2.4.0" \
    && blastn -version \
    && python --version \
    && samtools --version \
    && bedtools --version \
    && mafft --version \
    && gawk --version \
    && bc --version
