# Image for the MEDICC2 tree-building environment shared by: medicc2,
# jitter_correct, foreground, preprocess_medicc2 (4 modules). It used to be
# duplicated as modules/fg_modules/*/medicc2.yaml (one byte-identical copy
# per module) -- it's now a single tracked spec here, medicc2.yaml, since no
# module needs its own local copy once every module runs from the built
# image instead of building its own conda env.
#
# Build (from the repo root, on a machine that has Docker):
#   docker build -t ghcr.io/molonc/mrca-foreground-medicc2:1.1.2 \
#       -f containers/medicc2/Dockerfile .
#
# Push:
#   docker push ghcr.io/molonc/mrca-foreground-medicc2:1.1.2
#
# If you push somewhere other than ghcr.io/molonc, update medicc2_container
# in nextflow.config to match.
#
# Run on this cluster via Nextflow's `singularity` profile (see
# containers/r44/Dockerfile for why -- no Docker daemon here, only
# Apptainer/Singularity, which pulls `docker://...` references directly):
#   nextflow run main.nf -c <config> -profile singularity

FROM mambaorg/micromamba:1.5.8

# procps (`ps`): Nextflow's task wrapper re-execs itself inside the
# container and hard-requires `ps` for its local-executor trace collection
# (see medicc2.def's comment) -- mambaorg/micromamba's base Debian image
# doesn't include it.
USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER containers/medicc2/medicc2.yaml /tmp/medicc2.yaml

RUN micromamba install -y -n base -f /tmp/medicc2.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Sanity-check the tools this image exists to provide.
RUN medicc2 --help > /dev/null && python -c "import scgenome"

# scgenome -> umap -> numba: numba JIT-caches compiled functions next to
# their source file by default. That works at build time (writable layer)
# but not at run time from a built (read-only) Apptainer .sif -- see
# medicc2.def's comment. Point the cache somewhere writable instead.
ENV NUMBA_CACHE_DIR=/tmp
