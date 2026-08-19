# "base" = the shared R 4.4 environment used by every R-based module in this
# pipeline: annotate_bps, filter_cells, get_bp_file, get_df, get_hdp_file,
# nnd, plot_heatmaps, punish_baf, signals, signals_prep (10 modules).
#
# This is the Docker-elsewhere alternative to the primary build path,
# containers/base/base.def (built directly on this cluster with Apptainer --
# see that file and containers/README.md for the full "why bioconductor_docker
# + pak, not conda" explanation). Keep the two in sync: same base image tag,
# same apt packages, same pak::pkg_install() package list.
#
# Build (from the repo root, on a machine that has Docker):
#   docker build -t ghcr.io/molonc/mrca-foreground-base:1.0 \
#       -f containers/base/Dockerfile .
#
# Push:
#   docker push ghcr.io/molonc/mrca-foreground-base:1.0
#
# If you push somewhere other than ghcr.io/molonc, update base_container in
# nextflow.config to match.
#
# Run on this cluster: no Docker daemon is installed here, only
# Apptainer/Singularity. Nextflow's existing `singularity` profile pulls
# `docker://...` image references directly and lets Apptainer convert/cache
# them, so nothing else needs to change to run this cluster-side:
#   nextflow run main.nf -c <config> -profile singularity

FROM bioconductor/bioconductor_docker:RELEASE_3_20-R-4.4.2

ENV DEBIAN_FRONTEND=noninteractive

# System libs R packages' ./configure scripts need at build time -- see
# base.def's %post comment for what each group is for.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        pkg-config \
        git \
        wget \
        libcairo2-dev \
        libxt-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libpng-dev \
        libjpeg-dev \
        libtiff5-dev \
        libgit2-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# pak: resolves CRAN + Bioconductor + GitHub dependency graphs (including
# cross-package version constraints) together in one shot -- see base.def.
RUN Rscript -e 'install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))'

# Everything every fg_module actually `library()`s, resolved and installed
# together in one pak call. Keep this list identical to base.def's.
#
# Version pins work around CRAN packages having moved on from what this
# image's Bioconductor-3.20-era packages (ggtree 3.14.0 et al) were
# built/tested against -- see base.def's comment for the full explanation:
# - Deriv@4.2.0: doBy (pulled in transitively) Imports Deriv; CRAN's current
#   Deriv (4.3.0) added Rcpp-backed code that doesn't compile here
#   (references `R_ClosureFormals`, not a real public R C API function).
#   4.2.0 predates that C++ code (pure R).
# - ggplot2@3.5.2: ggplot2 4.0.0's S3->S7 refactor removed/renamed internal
#   helpers ggtree needs (`check_linewidth`). 3.5.2 is the last pre-4.0
#   release.
# VGAM: `signals` itself needs this at runtime for its beta-binomial model,
# but it's only in signals' own Suggests, which pak doesn't install unless
# asked -- add it explicitly.
RUN Rscript -e 'pak::pkg_install(c("ape", "argparse", "bioc::ComplexHeatmap", "data.table", "Deriv@4.2.0", "dplyr", "fitdistrplus", "ggnewscale", "ggplot2@3.5.2", "bioc::ggtree", "googlesheets4", "gtools", "bioc::HMMcopy", "httr", "mclust", "optparse", "plyr", "purrr", "RColorBrewer", "stringr", "tibble", "tidyr", "tidyverse", "VGAM", "vroom", "shahcompbio/signals", "molonc/dlptools"), ask = FALSE)'

# Sanity-check every package this image exists to provide.
RUN Rscript -e 'invisible(lapply(c("ape", "argparse", "ComplexHeatmap", "data.table", "dplyr", "fitdistrplus", "ggnewscale", "ggplot2", "ggtree", "googlesheets4", "gtools", "HMMcopy", "httr", "mclust", "optparse", "plyr", "purrr", "RColorBrewer", "stringr", "tibble", "tidyr", "tidyverse", "vroom", "signals", "dlptools"), function(p) library(p, character.only = TRUE)))'

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
