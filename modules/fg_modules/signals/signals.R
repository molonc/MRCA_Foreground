#!/usr/bin/env Rscript
library(optparse)
library(vroom)
library(signals)
library(dplyr)
# Define command line options with mandatory arguments first
option_list <- list(
  # Mandatory arguments
  make_option("--reads", type = "character", default = NULL,
              help = "Path to reads/CNbins file (mandatory, e.g., test_CNbins.csv.gz)"),
    
  make_option("--allele_id", type = "character", default = NULL,
              help = "Path to allele_id file (mandatory, e.g., allele_id.csv.gz)"),
  
  # Optional arguments
  make_option("--phased_haplotypes", type = "character", default = NULL,
              help = "Path to phased_haplotypes file (premade phased object)"),
  
  make_option("--output", type = "character", default = "test_hscn.csv.gz",
              help = "Output file name [default= %default]"),
  
  make_option("--ncores", type = "integer", default = 4,
              help = "Number of cores to use [default= %default]"),
  
  make_option("--selftransitionprob", type = "double", default = 0.999,
              help = "Self-transition probability [default= %default]"),
  
  make_option("--include_chr_y", action = "store_true", default = FALSE,
              help = "Include chromosome Y in analysis [default= FALSE]"),

  # NOTE: type = "logical" (not action = "store_true") -- main.nf's SIGNALS
  # process always passes an explicit value, `--include_chr_x TRUE` or
  # `--include_chr_x FALSE` depending on task.attempt. A store_true switch
  # takes no argument, so that trailing TRUE/FALSE token was left over and
  # rejected by optparse ("'FALSE' is not a valid option") on every call.
  make_option("--include_chr_x", type = "logical", default = TRUE,
              help = "Include chromosome X in analysis [default= TRUE]")

)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list,
                          description = "Run signals haplotype-specific copy number analysis",
                          usage = "usage: %prog --reads READS --allele_id ALLELE_ID [options]")
# Parse the command line arguments
opt <- parse_args(opt_parser)
if (opt$phased_haplotypes == "NULL") {
  opt$phased_haplotypes <- NULL
}
if (basename(opt$phased_haplotypes) == "dummy_phasing.csv") {
  opt$phased_haplotypes <- NULL
  message("Dummy phasing file detected, proceeding without phased haplotypes.")
}
# Check required arguments
required_args <- c("reads", "allele_id")
missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]
if (length(missing_args) > 0) {
  cat("Error: Missing required arguments:", paste(missing_args, collapse = ", "), "\n")
  cat("Use --help for more information.\n")
  quit(status = 1)
}
# Print configuration
cat("Running with parameters:\n")
cat(paste("Reads/CNbins file:", opt$reads, "\n"))
cat(paste("Allele ID file:", opt$allele_id, "\n"))
cat(paste("Output file:", opt$output, "\n"))
cat(paste("Number of cores:", opt$ncores, "\n"))
cat(paste("Self-transition prob:", opt$selftransitionprob, "\n"))
cat(paste("Phased haplotypes:", ifelse(is.null(opt$phased_haplotypes), "Not provided", opt$phased_haplotypes), "\n"))
cat(paste("Include Chr Y:", opt$include_chr_y, "\n"))
cat(paste("Include Chr X:", opt$include_chr_x, "\n"))


# Load data
cat("Loading data...\n")
CNbins <- vroom(opt$reads, col_types = cols(
          chr = col_character(),
          start = col_double(),
          end = col_double())) %>%
  mutate(
    chr   = as.character(chr),
    start = as.double(start),
    end   = as.double(end)
  )

allele_id <- vroom(opt$allele_id) %>%
  mutate(
    chromosome   = as.character(chromosome),
    start = as.double(start),
    end   = as.double(end)
  )  %>%
  rename(chr = chromosome)

## Filter for low mappability
CNbins <- subset(CNbins, is_low_mappability == FALSE)

# Remove X and Y chr depending on params
if (!opt$include_chr_x) {
  CNbins <- subset(CNbins, chr != "X")
  allele_id <- subset(allele_id, chr != "X")
  print("Removed X chr from analysis")
}
if (!opt$include_chr_y) {
  CNbins <- subset(CNbins, chr != "Y")
  allele_id <- subset(allele_id, chr != "Y")
  print("Removed Y chr from analysis")
}

# Format object
print("Formatting haplotypes")
haplotypes <- format_haplotypes_dlp(allele_id, CNbins)
haplotypes <- haplotypes %>%
  mutate(chr = as.character(chr),
        start = as.double(start),
        end = as.double(end))


# Prepare the function arguments
args <- list(
  CNbins = CNbins,
  haplotypes = haplotypes,
  firstpassfiltering = FALSE,
  ncores = opt$ncores,
  selftransitionprob = opt$selftransitionprob
)

# Add phased_haplotypes if provided add it to the args for medic
if (!is.null(opt$phased_haplotypes) && !opt$phased_haplotypes=="NULL") {
  cat("Loading phased haplotypes...\n")
  print(typeof(opt$phased_haplotypes))
  phased_haplotypes <- vroom(opt$phased_haplotypes, col_types = cols(
    chr   = col_character(),
    start = col_double(),
    end   = col_double()
  ))
  phased_haplotypes <- phased_haplotypes %>%
      mutate(chr = as.character(chr),
            start = as.double(start),
            end   = as.double(end))
  args$phased_haplotypes <- phased_haplotypes
} else{
  print("No Phasing Object was used")
}

# Make object
cat("Calling haplotype-specific copy number...\n")
hscn <- do.call(callHaplotypeSpecificCN, args)
cat(paste("Number of cells pre- hsscn", length(unique(CNbins$cell_id))))
cat(paste("Number of cells in hsscn", length(unique(hscn$data$cell_id))))
# Write output
output_path <-  opt$output
cat(paste("Writing results to", output_path, "...\n"))
data.table::fwrite(hscn$data, output_path)
cat("hscn generated!\n")