# Containers

This pipeline uses exactly two environments across all 14 `fg_modules`:

| Image | Def file (primary) | Dockerfile (alternative) | Spec | Used by |
|---|---|---|---|---|
| `base` | [base/base.def](base/base.def) | [base/Dockerfile](base/Dockerfile) | package list inline in both (see below) | annotate_bps, filter_cells, foreground, get_bp_file, get_df, get_hdp_file, nnd, plot_heatmaps, punish_baf, signals, signals_prep (11 modules) |
| `medicc2` | [medicc2/medicc2.def](medicc2/medicc2.def) | [medicc2/Dockerfile](medicc2/Dockerfile) | [medicc2/medicc2.yaml](medicc2/medicc2.yaml) | medicc2, jitter_correct, preprocess_medicc2 (3 modules) |

Each module's `main.nf` carries only a `container` directive
(`params.base_container` / `params.medicc2_container`, set in
[nextflow.config](../nextflow.config)) -- there is no per-module `conda`
directive or env file anymore.

(The `base` image used to be called `r44`, after "R 4.4" -- renamed since
it's really just "the shared R environment," not something specific to R
4.4.)

## Why -- and why `base` isn't a conda/micromamba env anymore

`base` used to be a micromamba env built from a big pinned conda spec
(`base.yml`, ~470 packages, captured from a prior working environment) --
originally because `mamba env create` for the old `r44.yml` was getting
killed (SIGTERM / exit 143) building live on this pipeline's shared analysis
host under memory pressure, and containerizing meant that solve only had to
happen once, not on every run.

That got the image built, but the conda approach then fought itself over
every dependency `base` actually needs:
- Two packages this pipeline needs are GitHub-only --
  [signals](https://github.com/shahcompbio/signals)
  (`modules/fg_modules/signals/signals.R` calls `library(signals)`) and
  [dlptools](https://github.com/molonc/dlptools) (an internal Aparicio Lab
  package, `library()`'d by several other `fg_modules`) -- so neither could
  be conda-pinned at all; both had to be `remotes::install_github()`'d
  separately, with their *own* R-level Imports hand-added to the conda spec
  one at a time as each missing one surfaced (`ggtree` -> `treeio` ->
  `GenomicRanges`/`GenomeInfoDb`/`S4Vectors`/`ComplexHeatmap` -> `HMMcopy`,
  `Cairo`/`ggrastr`, and more).
- One of those transitive deps, `sigminer`, has no bioconda build for R 4.4
  on linux-64 at any point (its builds jump straight from an r43 build to
  r45) -- it needed installing from CRAN source directly, mid-`%post`.
- Bioconductor packages live on the separate `bioconda` channel, which
  conda's own solver doesn't reconcile against R-level version constraints:
  a pin that looked fine (`rlang=1.1.6`) still broke a from-source install
  that needed `rlang >= 1.1.7`, with no warning until that particular
  install step ran.

Each of those was a separate rebuild (each a 15-20+ minute conda solve) to
discover. `base.def`/`Dockerfile` now build on
[`bioconductor/bioconductor_docker`](https://hub.docker.com/r/bioconductor/bioconductor_docker)
(R 4.4 series) and use [`pak`](https://pak.r-lib.org) to install everything
in one call instead: pak resolves CRAN, Bioconductor (`bioc::pkg`), and
GitHub (`owner/repo`) package specs *together*, including cross-package
version requirements, in a single dependency-graph solve -- which is also
what signals' and dlptools' own upstream Dockerfiles use pak for. There's no
`base.yml` anymore; the package list lives directly in `base.def` (and is
mirrored in `Dockerfile` -- keep the two in sync).

## Primary path: build with Apptainer, directly on this cluster

This cluster has no Docker/Buildah/Podman (checked directly -- none
installed, and rootless container engines can't work here anyway since
unprivileged user namespaces are disabled at the kernel level:
`user.max_user_namespaces = 0`). It does have **Apptainer 1.3.1**, which
doesn't need either of those -- it builds via its own fakeroot emulation.

Build both images from the repo root:
```
apptainer build containers/base/base.sif containers/base/base.def
apptainer build containers/medicc2/medicc2.sif containers/medicc2/medicc2.def
```

No push, no registry, no auth -- `nextflow.config` already points
`base_container`/`medicc2_container` at these exact `.sif` paths, so once
they exist, `nextflow run main.nf -c <config> -profile singularity` (or no
`-profile` at all -- singularity is the default) just works.

`base`'s build no longer depends on solving a ~470-package conda env on this
shared, memory-contended host (see "Why" above) -- it's an apt-get install
plus one `pak::pkg_install()` call now, considerably lighter. `medicc2`
still is a micromamba env (`medicc2.yaml`); if that one ever gets killed
under memory pressure, fall back to the Docker path below for it.

`.sif` files are gitignored (multi-GB binaries) -- rebuild from the `.def`
files rather than committing them. Bump the filename/a version suffix if
you rebuild with a changed spec, so old runs stay reproducible against the
image they actually used.

## Alternative: build with Docker elsewhere, publish to a registry

If an Apptainer build gets killed under load, or you'd rather fully isolate
the build from this host's resource contention, the Dockerfiles are kept
as an equivalent alternative path:

```
docker build --platform linux/amd64 -t ghcr.io/molonc/mrca-foreground-base:1.0 -f containers/base/Dockerfile .
docker push ghcr.io/molonc/mrca-foreground-base:1.0

docker build --platform linux/amd64 -t ghcr.io/molonc/mrca-foreground-medicc2:1.1.2 -f containers/medicc2/Dockerfile .
docker push ghcr.io/molonc/mrca-foreground-medicc2:1.1.2
```
(`--platform linux/amd64` matters if your build machine isn't already
x86_64 Linux -- this cluster is, and both images' pinned/prebuilt package
binaries are Linux-x86_64-specific.)

Then switch `base_container`/`medicc2_container` in `nextflow.config` back
to the `docker://ghcr.io/...` references, and make sure the packages are
either public or that Apptainer has registry credentials configured on this
cluster (`apptainer remote login` / `$HOME/.singularity/docker-config.json`)
-- there's no interactive login flow otherwise.

## Caveats

- `base.def`/`Dockerfile` pin `bioconductor/bioconductor_docker:RELEASE_3_20-R-4.4.2`
  -- the newest Bioconductor-3.20 (R 4.4 series) tag on Docker Hub as of this
  writing. If that exact tag is ever removed, `3.20`/`3.19`/`3.19-R-4.4.1`
  are the other R-4.4 options; otherwise bump to a newer Bioconductor/R pair
  (and re-verify the package list still resolves).
- `base`'s package list (the `pak::pkg_install(...)` call) is intentionally
  a flat list of only what's directly `library()`'d by `fg_modules`/`helper`
  R scripts, not a full environment dump -- pak resolves and installs every
  transitive dependency (Bioconductor and CRAN and GitHub alike) on its own.
  Keep `base.def` and `Dockerfile`'s lists identical if you add another
  package.
- `medicc2` is unrelated to any of the above -- still its own micromamba env
  from `medicc2.yaml`, untouched.
- There is no conda fallback anymore -- `-profile conda` was removed from
  `nextflow.config` since no process declares a `conda` directive. If you
  ever need a non-container path again, you'd need to re-add both.
