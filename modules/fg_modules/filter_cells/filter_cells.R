#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(vroom)
library(data.table)
library(argparse)

# Create argument parser
parser <- ArgumentParser(description='Filter cells based on BAF punishment and divergence')
parser$add_argument('--hscn', type='character', required=TRUE,
                    help='Path to hscn file')
parser$add_argument('--punishment_data', type='character', required=TRUE,
                    help='Path to punishment data file')
parser$add_argument('--divergence_data', type='character', required=TRUE,
                    help='Path to divergence data file')
parser$add_argument('--output', type='character', required=TRUE,
                    help='Output path for filtered hscn file')
parser$add_argument('--error_threshold', type='double', default=50,
                    help='Percent error threshold for bad allele calls (default: 50)')
parser$add_argument('--filter_divergent', type='character', default='true',
                    help='Whether to filter divergent cells (default: true)')

# Parse arguments
args <- parser$parse_args()

# Load in allele assigned copy-number data
cat("Loading hscn data...\n")
hscn <- vroom::vroom(args$hscn)

# Load punishment data
cat("Loading punishment data...\n")
pdat <- vroom::vroom(args$punishment_data)

# Identify bad allele calls
bad_allele_calls <- subset(pdat, percent_error > args$error_threshold)$cell_id
cat(sprintf("Found %d cells with percent_error > %g\n", 
            length(bad_allele_calls), args$error_threshold))

# Load divergent cell data
cat("Loading divergence data...\n")
ddat <- vroom::vroom(args$divergence_data)

filter_divergent <- tolower(args$filter_divergent) == "true"

# Identify divergent cells
if (filter_divergent) {
    divergent_cells <- subset(ddat, is_divergent == TRUE)$cell_id
    cat(sprintf("Found %d divergent cells\n", length(divergent_cells)))
} else {
    divergent_cells <- character(0)
    cat("Skipping divergent cell filtering (filter_divergent=false)\n")
}

# Combine cells to remove
combined_removal <- unique(c(bad_allele_calls, divergent_cells))
cat(sprintf("Total cells to remove: %d\n", length(combined_removal)))

# Filter hscn
initial_cells <- length(unique(hscn$cell_id))
hscn_filtered <- subset(hscn, !cell_id %in% combined_removal)
final_cells <- length(unique(hscn_filtered$cell_id))

cat(sprintf("Cells before filtering: %d\n", initial_cells))
cat(sprintf("Cells after filtering: %d\n", final_cells))
cat(sprintf("Cells removed: %d\n", initial_cells - final_cells))

# Write filtered data
cat(sprintf("Writing filtered data to %s\n", args$output))
data.table::fwrite(hscn_filtered, args$output)

cat("Filtering complete!\n")