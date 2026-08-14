#!/usr/bin/env Rscript

#' DLP Data Processing Script
#' 
#' This script processes DLP (Direct Library Preparation) data by loading allele counts,
#' HMMcopy metrics, and read data, then filtering based on user-provided parameters.
#'
#' @param sample_id Sample identifier to filter data
#' @param library_id Library identifier to filter data
#' @param include_s_phase Logical, whether to include S phase cells (default: FALSE)
#' @param include_dead_cells Logical, whether to include dead cells (C2) along with alive cells (C1) (default: FALSE)
#' @param quality_threshold Numeric value between 0 and 1 for cell quality filtering (default: 0.75)
#' @param data_dir Base directory for data files (default: uses hardcoded path)
#' @param output_dir Directory to save output files (default: current working directory)
#' @param PDX Logical, whether to extract sample_id/library_id from cell_id for PDX samples (default: FALSE)

# Parse command line arguments
library(optparse)
library(vroom)
library(data.table)
library(dplyr)
library(stringr)

option_list <- list(
  make_option(c("--sample_id"), type="character", help="Sample identifier"),
  make_option(c("--library_id"), type="character", help="Library identifier"),
  make_option(c("--reads"), type = 'character'),
  make_option(c("--alleles"), type = 'character'),
  make_option(c("--metrics"), type = 'character'),
  make_option(c("--include_s_phase"), type="logical", default=FALSE, 
              help="Include S phase cells [default: %default]"),
  make_option(c("--s_phase_quality_threshold"), type="double", default=0.5,
              help="Quality threshold between 0 and 1 for s phase cells. if include sp hase cells is True [default: %default]"),
  make_option(c("--include_dead_cells"), type="logical", default=FALSE, 
              help="Include dead cells (C2) along with alive cells (C1) [default: %default]"),
  make_option(c("--include_low_quality_cells"), type="logical", default=FALSE, 
              help="Include Low quality cells identified by classifier [default: %default], uses quality_threshold"),
  make_option(c("--quality_threshold"), type="double", default=0.75,
              help="Quality threshold between 0 and 1 [default: %default]"),
  make_option(c("--PDX"), type="logical", default=FALSE,
              help="Extract sample_id/library_id from cell_id for PDX samples [default: %default]"),
  make_option(c("--num_cells_returned"), type="integer", default=NULL,
              help="If set, keep only the top N cells ranked by quality after all other filters [default: %default]")
)


opt_parser <- OptionParser(option_list=option_list, 
                           description="Process and filter DLP data based on specified parameters")
opt <- parse_args(opt_parser)
print(opt)
# Validate required arguments
if (is.null(opt$sample_id) || is.null(opt$library_id)) {
  stop("Required arguments missing. Please provide --scid, --sample_id, and --library_id")
}

# Allele file path
allele_file <- opt$alleles
# handle metrics
metrics_file <- opt$metrics
# handle reads
# For reads file, find any file ending with reads.csv.gz in the hmmcopy directory
reads_file <- opt$reads

print(paste("Allele File", allele_file))
print(paste("MEtrics File", metrics_file))
print(paste("reads File", reads_file))
# Load data with error handling
tryCatch({
  alleles <- vroom::vroom(allele_file)
  alleles <- alleles %>%
    distinct() %>%
    filter(!is.na(start))%>%
    mutate(
      start      = as.character(start),
      end        = as.character(end),
      hap_label  = as.character(hap_label),
      allele_id  = as.character(allele_id),
      readcount  = as.integer(readcount),
      chromosome = as.character(chromosome),
      cell_id    = as.character(cell_id)
    )

  if (nrow(alleles) ==0) {
    stop("ERROR: 'alleles' has 0 rows")
  }
  message("Successfully loaded allele counts data")
  message(paste(length(unique(alleles$cell_id)), "Number of cells found in alleles"))
}, error = function(e) {
  stop(paste("Error loading allele file:", e$message, allele_file))
})
tryCatch({
  mets <- vroom(metrics_file)
  mets <- mets %>%
    distinct() %>%
    mutate(
      is_s_phase      = as.logical(is_s_phase),
      is_control      = as.logical(is_control),
      is_contaminated = as.logical(is_contaminated),
      quality         = as.double(quality)
    )
  message("Successfully loaded metrics data")
  message(paste(length(unique(mets$cell_id)), "Number of cells found in metrics"))
  message(paste(length(unique(mets$library_id)), "Number of libraries found in metrics"))
}, error = function(e) {
  stop(paste("Error loading metrics file:", e$message))
})
tryCatch({
  reads <- vroom(reads_file)
  reads <- reads %>%
    distinct() %>%
    filter(!is.na(start))%>%
    mutate(
      cell_id           = as.character(cell_id),
      chr               = as.character(chr),
      start             = as.character(start),
      end               = as.character(end),
      reads             = as.double(reads),
      state             = as.integer(state),
      is_low_mappability = as.logical(is_low_mappability))
  message("Successfully loaded reads data")
  message(paste(length(unique(reads$cell_id)), "Number of cells found in reads"))
}, error = function(e) {
  stop(paste("Error loading reads file:", e$message))
})

message("Filtering metrics data...")
# Validate quality threshold is between 0 and 1
if (opt$quality_threshold < 0 || opt$quality_threshold > 1) {
  warning("Quality threshold should be between 0 and 1. Using default value of 0.75.")
  opt$quality_threshold <- 0.75
}

message("Allele Files has ", length(unique(alleles$cell_id)))
message("Reads Files has ", length(unique(reads$cell_id)))
message("Metrics Files has ", length(unique(mets$cell_id)))
# Apply filters to metrics data
message('column names in met')
mets %>% colnames() %>% print()
message("Quality sorted")
print(head(mets %>% arrange(desc(quality))), Inf)

n <- length(unique(mets$cell_id))
message(paste(n, "Cells were found in metrics file"))
# Remove Control Cells
filtered_cells <- mets %>%
  filter(is_control==FALSE) 
message(paste(length(unique(filtered_cells$cell_id)), "cells remain after control cell filtering"))
if (filtered_cells %>% nrow() == 0) {
  stop("ERROR: No cells found after filtering out control cells (is_control==FALSE)")
}

# Apply cell state filter (C1 = alive, C2 = dead)
if (!opt$include_dead_cells) {
  message("Including only alive cells (C1)")
  filtered_cells <- filtered_cells %>% filter(cell_call == "C1")
} else {
  message("Including both alive (C1) and dead (C2) cells")
  filtered_cells <- filtered_cells %>% filter(cell_call %in% c("C1", "C2"))
}
message(paste(length(unique(filtered_cells$cell_id)), "cells remain after dead cell filtering"))

# Always filter out contaminated cells
message("Filtering out contaminated cells")
# Check if is_contaminated column exists before filtering
if ("is_contaminated" %in% colnames(filtered_cells)) {
  filtered_cells <- filtered_cells %>% filter(!is_contaminated)
} else {
  message("Warning: is_contaminated column not found in metrics data. Skipping contamination filtering.")
}
message(paste(length(unique(filtered_cells$cell_id)), "cells remain after contaminated cell filtering"))

# Remove S phase cells
if (opt$include_s_phase == FALSE) {
  if ("is_s_phase" %in% colnames(filtered_cells)) {
    filtered_cells <- filtered_cells %>%
      filter(is_s_phase == opt$include_s_phase)}
    message(paste(length(unique(filtered_cells$cell_id)), "cells remain after s phase filtering"))
  if (filtered_cells %>% nrow() == 0) {
    stop("ERROR: No cells found after filtering out S phase cells")}
} else if (opt$include_s_phase == TRUE) {
  if ("is_s_phase" %in% colnames(filtered_cells)) {
    s_phase_cells <- filtered_cells %>%
      filter(is_s_phase == TRUE) %>%
      filter(quality >= opt$s_phase_quality_threshold)
  }
  filtered_cells <- filtered_cells %>%
    filter(is_s_phase == FALSE)
  message(paste(length(unique(filtered_cells$cell_id)), "cells remain after s phase filtering"))
}

# Remove low quality cells
if (opt$include_low_quality_cells == FALSE) {
  if ("quality" %in% colnames(filtered_cells)) {
    filtered_cells <- filtered_cells %>%
      filter(quality > opt$quality_threshold)
  }
  if (filtered_cells %>% nrow() == 0) {
    print(head(filtered_cells))
    message(opt$quality_threshold)
    stop("ERROR: No cells found after filtering out cells with low quality")
  }
  message(paste(length(unique(filtered_cells$cell_id)), "cells remain after quality filtering"))
}

if (opt$include_s_phase == TRUE) {
  filtered_cells <- bind_rows(filtered_cells, s_phase_cells)
  message(paste(length(unique(filtered_cells$cell_id)), "cells remain after adding back high-quality S phase cells"))
}

# Optionally keep only the top-N cells by quality
if (!is.null(opt$num_cells_returned)) {
  message(paste("Selecting top", opt$num_cells_returned, "cells by quality"))
  filtered_cells <- filtered_cells %>%
    arrange(desc(quality)) %>%
    slice_head(n = opt$num_cells_returned)
  message(paste(nrow(filtered_cells), "cells remain after top-N selection"))
}

# Select only cell_id for filtering
filtered_cells <- filtered_cells %>% select(cell_id)
write.table(filtered_cells, "filtered_cells.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
message(paste("Found", nrow(filtered_cells), "cells passing filters:"))
message(paste("- Sample ID:", opt$sample_id))
message(paste("- Library ID:", opt$library_id))
message(paste("- S phase inclusion:", opt$include_s_phase))
message(paste("- Include dead cells:", opt$include_dead_cells))
message(paste("- Include low quality cells:", opt$include_low_quality_cells))
message(paste("- Quality threshold:", opt$quality_threshold))
message(paste("- PDX mode (for ID extraction):", opt$PDX))

# Get top 10 unique cell IDs
# top10 <- filtered_cells %>%
#   distinct(cell_id) %>%
#   slice_head(n = 10) %>%
#   pull(cell_id)

# # Filter the data to only include those 10 cells
# filtered_cells <- filtered_cells %>%
#   filter(cell_id %in% top10)

# Filter reads data based on filtered cells
message("Filtering reads data based on filtered cells...")
filtered_reads <- reads %>%
  filter(cell_id %in% filtered_cells$cell_id)
filtered_alleles <- alleles %>%
  filter(cell_id %in% filtered_cells$cell_id)
filtered_mets <- mets %>%
  filter(cell_id %in% filtered_cells$cell_id)


## Filter out low mappability areas
filtered_reads <- filtered_reads %>%
   filter(is_low_mappability == FALSE)

# Output summary
message(paste("Processed data for SCID:", opt$scid))
message(paste("Sample ID:", opt$sample_id))
message(paste("Library ID:", opt$library_id))
message(paste("Final dataset contains", nrow(filtered_reads), "rows"))


# Create output filenames
reads_output <- file.path(paste0(opt$library_id, "_reads.csv.gz"))
alleles_output <- file.path(paste0(opt$library_id, "_alleles.csv.gz"))
metrics_output <- file.path(paste0(opt$library_id, "_metrics.csv.gz"))

print(head(filtered_reads))
print(head(filtered_alleles))
print(head(filtered_mets))
# Check for zero rows and stop with an error if found
if (nrow(filtered_reads) == 0) {
    stop("ERROR: 'filtered_reads' has 0 rows")
}
if (nrow(filtered_alleles) == 0) {
    stop("ERROR: 'filtered_alleles' has 0 rows")
}
if (nrow(filtered_mets) == 0) {
    stop("ERROR: 'filtered_mets' has 0 rows")
}



# Write filtered results to files
message("Writing output files...")
data.table::fwrite(filtered_reads, reads_output)
message(paste("Filtered reads written to:", reads_output))

data.table::fwrite(filtered_alleles, alleles_output)
message(paste("Filtered alleles written to:", alleles_output))

data.table::fwrite(filtered_mets, metrics_output) 
message(paste("Metrics written to:", metrics_output))