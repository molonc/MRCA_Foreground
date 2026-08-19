# MEDICC2 Foreground Analysis Pipeline

A Nextflow (DSL2) pipeline for analyzing single-cell whole-genome sequencing (scWGS) data,
focusing on copy number analysis, phylogenetic reconstruction, and "foreground" (relative-to-
ancestor) copy number state identification using [MEDICC2](https://www.nature.com/articles/s41587-022-01340-z).

## Overview

This pipeline processes single-cell DNA sequencing data to:
- Identify and (optionally) filter divergent cells using nearest-neighbour-distance (NND) analysis
- Perform copy number segmentation and allele-specific analysis
- Reconstruct phylogenetic trees using MEDICC2
- Calculate foreground changes relative to ancestral states
- Generate visualizations of copy number profiles

Input hmmcopy/allele-count data must already be available locally (see below) — this pipeline
does not fetch data from any external system.

## Requirements

- [Nextflow](https://www.nextflow.io/) (`>=22.10`, DSL2)
- Conda or Mamba (every process supplies its own `conda` directive pointing at a tracked
  environment file — see `modules/fg_modules/*/{r44.yml,medicc2.yaml}`)

No manual R/Python package installation is required beyond that — Nextflow builds/caches each
module's environment automatically from its tracked `.yml` file the first time it runs.

## Required Input Files

Every sample needs three files: an hmmcopy reads file, an hmmcopy metrics file, and an
allele-count file. **The pipeline does not search for these itself** — no directory-naming
convention, no filename globbing. You tell it exactly where each file is, per sample, via a
required samplesheet.

Start a new sheet from [`input_files/samplesheets/TEMPLATE.csv`](input_files/samplesheets/TEMPLATE.csv).
It's a header CSV with one row per sample and four required columns:

| column | meaning |
|---|---|
| `sample_id` | This row's tracking id — see below |
| `reads` | Path to this sample's `*hmmcopy_reads.csv.gz` |
| `allele` | Path to this sample's `*allele_count.csv.gz` |
| `metrics` | Path to this sample's `*hmmcopy_metrics.csv.gz` |

`sample_id` is the pipeline's **only** tracking key: it's the tag on every process, and it names
that row's results subdirectory (`${out_dir}/${sample_id}/...`). It must be non-empty and unique
across the sheet — the pipeline checks both up front — but its structure is entirely up to you.
It can be a plain library ID, or a composite string if that's what's actually unique in your data;
e.g. this repo's own cohorts use `<library_id>_<sample_id>_<exp_con>` (matching their old
directory-naming convention) as the value in this one column.

`reads`/`allele`/`metrics` paths can be absolute, or relative to the samplesheet's own location
(so a samplesheet and its data can be shared as a self-contained unit without hardcoding a
machine-specific prefix). `params.input` pointing at a samplesheet is required — the pipeline
errors out immediately, with the exact missing column/file/duplicate-id named, if it isn't set or
a row is incomplete.

If your data already lives in the older `<library_id>_<sample_id>_<exp_con>/` per-folder layout
(one dir per library, each holding exactly one file matching each of the three suffixes above),
`helper_scripts/generate_samplesheet.py` will build/migrate a samplesheet for you from it — see
that script's docstring. It composes `sample_id` as `<library_id>_<sample_id>_<exp_con>` (so
existing result paths don't change) and flags (skips, rather than guesses) any library whose
folder has zero or more-than-one match for a given file type, so an ambiguous folder can't
silently produce a wrong run.

## Running the Pipeline

Cohort-specific parameters (`out_dir`, `phasing_object`, `input`, ...) are supplied via
a per-cohort config under `conf/main/` (production cohorts) or `conf/dev/` (fast/scratch runs) —
see those directories for examples, including `conf/main/paper_example.config`, a small
reviewer-facing example run.

```bash
nextflow run main.nf -c conf/main/<cohort>.config -profile conda -resume
```

Key params (set per-cohort in `conf/main/*.config` / `conf/dev/*.config`):

| param | meaning |
|---|---|
| `input` | **Required.** Path to the samplesheet (see above) |
| `out_dir` | Output directory |
| `phasing_object` | Path to a pre-computed SIGNALS phasing object (optional; a dummy empty file is used if unset) |
| `tcn` | `true` = total-copy-number-only MEDICC2 mode; `false` = allele-aware mode with WGD detection |
| `include_s_phase` | Whether to retain S-phase cells during QC filtering |
| `filter_divergent` | Whether to drop divergent cells identified by the NND/BAF steps |
| `num_cells_returned` | If set, keep only the top-N cells by quality after filtering |

> `data_dir` is no longer read by the pipeline (file locations now come entirely from the
> samplesheet) but is still set in some existing cohort configs; it's harmless to leave and safe
> to delete whenever those configs are next touched.

## Pipeline Steps

Each step is a Nextflow module under `modules/fg_modules/`; see each module's `meta.yml` for its
exact inputs/outputs.

1. **Data Preparation** (`signals_prep` — `hmmcopy_SIGNALS_prep.R`)
   Filters cells based on quality metrics; excludes S-phase and contaminated cells.
2. **Allele Assignment** (`signals` — `signals.R`)
   Assigns alleles to copy number bins and performs haplotype phasing (SIGNALS).
3. **Divergence Analysis** (`nnd` — `nnd_beta_divergence_punishment.R`)
   Calculates nearest-neighbour distances; identifies divergent cells via a fitted Beta distribution.
4. **BAF Quality Assessment** (`punish_baf` — `BAF_punishment_score_metric.R`)
   Evaluates B-allele-frequency integerness; calculates an error metric per cell.
5. **Cell Filtering** (`filter_cells` — `filter_cells.R`)
   Removes cells based on divergence and BAF scores (optional, `params.filter_divergent`).
6. **MEDICC2 Analysis** (`preprocess_medicc2`, `jitter_correct`, `medicc2`)
   Prepares MEDICC2 input, runs MSPCF-based jitter-corrected segmentation, reconstructs the
   phylogenetic tree, and calculates total and/or allele-specific copy numbers.
7. **Foreground Analysis** (`foreground` — `foreground_ancestor.R`)
   Computes copy number changes relative to ancestors; identifies gains/losses along branches.
8. **Breakpoint annotation & HDP export** (`get_df`, `annotate_bps`, `get_bp_file`, `get_hdp_file`)
   Assembles the full reads dataframe with foreground state, annotates segment edges relative to
   centromeres/telomeres, and exports the breakpoint and HDP files described below.
9. **Visualization** (`plot_heatmaps` — `plot_heatmaps.R`)
   Generates CNV/foreground/allele-specific heatmaps, phylogenetically ordered by the MEDICC2 tree.

## Output Files

Every output lands under a subdirectory named after that row's `sample_id`:

```
results/
├── ${sample_id}/
   ├── ${sample_id}_breakpoint_file.csv.gz  # See below on breakpoint_file format
   ├── ${sample_id}_hdp_segs.csv.gz         # See below on HDP file format
   ├── CNV_ABSOLUTE_${sample_id}.png        # Heatmap of absolute (background) CN
   ├── CNV_FOREGROUND_${sample_id}.png      # Heatmap of foreground CN
   ├── CNV_FOREGROUND_EDGES_${sample_id}.png  # Heatmap with breakpoint-edge annotation
```

### Breakpoint file

One of the main outputs of the pipeline is the breakpoint file.

The file is formatted kind of like a segment file. Each cell has the complete set of chromosomes,
whether there is a breakpoint there or not.

Here are the necessary columns and what they mean:

- `cell_id`: cell_id
- `chr`: chromosome
- `start`: start of segment
- `end`: end of segment
- `state`: FOREGROUND state of segment.
   - `-100` indicates a centromere
   - `-101` indicates a p-telomere (start telo)
   - `-102` indicates a q-telomere (end telo)
- `seg_width`: segment width
- `edge`:
   - `NA` indicates there is no breakpoint here (centromere/telomere regions, and regions with
     fg state == 0, are `NA`)
   - `TELO-BP`: segment is telomere-bound, upstream of a centromere (we start at the telomere and
     then there's a breakpoint)
   - `BP-TELO`: segment is telomere-bound, downstream of a centromere
   - `BP-CENTRO`: segment is centromere-bound, upstream of the centromere
   - `CENTRO-BP`: segment is centromere-bound, downstream of the centromere
   - `INNER`: segment is not centromere- or telomere-bound
- `breakpoint`: the exact location of the breakpoint; `NA` if the segment doesn't represent one
- `median_loh` / `LOH`: median loss-of-heterozygosity indicator across the segment's supporting
  reads (0 = retained heterozygosity, 1 = LOH), where allele-specific data is available
- `median_wgd` / `WGD`: median whole-genome-doubling indicator for the segment/cell, as reported
  by MEDICC2
- `background_CNV` / `background`: the background (non-foreground) copy-number state underlying
  the segment, i.e. the absolute copy number minus the foreground state

> The three fields above are populated from `helper/breakpoint_file_functions.R`'s
> `add_median_info`/`add_median_info_noallele` — double-check the exact semantics against that
> file if you're relying on them for a specific analysis, since this description was filled in
> from reading the code rather than from a domain write-up.

### HDP file

This file is the breakpoint file, reduced to only the true segment breakpoints, for input to HDP.

## Citation

If you use this pipeline, please cite:
- MEDICC2: [Kaufmann et al., 2022](https://www.nature.com/articles/s41587-022-01340-z)
- SIGNALS: Marc Williams, https://github.com/shahcompbio/signals
- Jitter-correction script: Matt Myers

See also `CITATION.cff`.

## License

This project is licensed under the MIT License — see the `LICENSE` file for details.

## Contact

For questions or issues, please open an issue on the
[molonc/medicc2_foreground](https://github.com/molonc/medicc2_foreground) repository, or contact
the Aparicio Lab.
