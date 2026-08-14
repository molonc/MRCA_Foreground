#!/usr/bin/env Rscript

# BAF_punishment_score_metric.R
# Script to calculate BAF integerness and error metrics

suppressPackageStartupMessages({
  library(signals)
  library(dlptools)
  library(dplyr)
  library(tidyr)
  library(vroom)
  library(ggplot2)
  library(HMMcopy)
  library(tibble)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option("--hscn", type="character", help="Path to hscn file"),
  make_option("--metrics", type="character", help="Path to metrics file"),
  make_option("--tree", type="character", default=NULL, help="Path to tree file (optional)"),
  # make_option("--output_dir", type="character", default=".", help="Output directory"),
  make_option("--output", type="character", default=NULL, help="Output file name for results"),
  make_option("--epsilon_cutoff", type="numeric", default=0.10, help="Cutoff value for epsilon error"),
  make_option("--bin_size", type="numeric", default=5e5, help="Bin size for normalization")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Function to create segments with median BAF
readstoBaf <- function(reads, feature_to_seg, na.replace = -1) {
  
  reads$feature <- reads[[feature_to_seg]]
  
  longseg <- reads[, c("chr", "start", "end", "feature", "copy", "cell_id", "BAF")]
  longseg <- longseg[order(longseg$cell_id, longseg$chr, longseg$start), ]
  
  longseg$feature[is.na(longseg$feature)] <- na.replace
  longseg$copy[is.na(longseg$feature)] <- na.replace
  
  longseg$encode <- as.numeric(longseg$chr) * 100 + longseg$feature
  longrle <- rle(longseg$encode)
  longseg$rle <- rep(1:length(longrle$lengths), longrle$lengths)
  medseg <- group_by(.data=longseg, chr, feature, cell_id, rle) %>% 
    dplyr::summarise(
      start = min(start), 
      end = max(end), 
      median_copy = median(copy, na.rm = TRUE), 
      median_baf = median(BAF, na.rm = TRUE)
    )
  shortseg <- medseg[, c("chr", "start", "end", "feature", "median_copy", "cell_id", "median_baf")]
  shortseg$chr <- factor(shortseg$chr, c(1:22, "X", "Y"))
  shortseg <- shortseg[with(shortseg, order(cell_id, chr, start, end)), ]
  shortseg$width <- shortseg$end - shortseg$start + 1
  return(shortseg)
}

# Load colors
CNV_COLOURS <- dlptools::CNV_COLOURS

# Load data
message("Loading data...")
hscn <- vroom::vroom(opt$hscn)
metrics <- vroom::vroom(opt$metrics)

# Process cell ordering using tree if provided
if (!is.null(opt$tree) && file.exists(opt$tree)) {
  message("Reading tree file...")
  tree <- ape::read.tree(opt$tree)
  tree <- ape::drop.tip(tree, "diploid")
  order <- tree$tip.label
  hscn$cell_id <- factor(hscn$cell_id, levels = order)
}

# Factor chromosomes
hscn$chr <- factor(hscn$chr, levels = c(1:22, "X"))

# Make a segment with the median BAF for each segment
message("Creating BAF segments...")
baf_segs <- readstoBaf(hscn, "state")
baf_segs <- baf_segs %>% dplyr::rename(state = feature)

# Calculate integerness metrics
message("Calculating integerness metrics...")
baf_segs$non_integerness <- baf_segs$median_baf * baf_segs$state
baf_segs$round <- round(baf_segs$non_integerness)
baf_segs$epislon <- abs(baf_segs$round - baf_segs$non_integerness) # Difference between expected and observed

# Generate density plot for epsilon
p1 <- ggplot(baf_segs, aes(epislon)) + 
  geom_density() +
  ggtitle("Distribution of Epsilon (Non-integerness)") +
  xlab("Epsilon") + 
  ylab("Density")

# Save epsilon density plot
density_plot_path <-  "epsilon_density.png"
message("Saving epsilon density plot to ", density_plot_path)
ggsave(density_plot_path, p1, width = 8, height = 6)

# Calculate error metrics
message("Calculating error metrics...")
baf_segs$error <- ifelse(baf_segs$epislon > opt$epsilon_cutoff, TRUE, FALSE)
baf_segs$num_bins <- baf_segs$width/opt$bin_size

# Create percent error calculation
error <- baf_segs %>% group_by(cell_id, error) %>% summarize(error_factor = sum(num_bins, na.rm = T))
error_t <- subset(error, error == TRUE)
error_t <- error_t %>% dplyr::rename(error_factor_t = error_factor) %>% select(-error)
error_f <- subset(error, error == FALSE)
error_f <- error_f %>% dplyr::rename(error_factor_f = error_factor) %>% select(-error)

# Use full_join to capture all cells
overall_error <- full_join(error_t, error_f, by = "cell_id")

# Replace NA values with 0
overall_error <- overall_error %>%
  mutate(
    error_factor_t = replace_na(error_factor_t, 0),
    error_factor_f = replace_na(error_factor_f, 0),
    percent_error = (error_factor_t / (error_factor_t + error_factor_f)) * 100
  )

# Save results to file
results_output <- "baf_integerness_results.csv.gz"
# message("Saving results to ", results_output)
vroom::vroom_write(overall_error, results_output)

# Join with multiplier data if it exists in metrics
if ("multiplier" %in% colnames(metrics)) {
  message("Creating multiplier boxplot...")
  punishment <- left_join(overall_error, select(metrics, cell_id, multiplier))
  
  boxplot <- ggplot(punishment, aes(as.factor(multiplier), percent_error, col = as.factor(multiplier))) + 
    geom_boxplot() + 
    scale_color_manual(values = CNV_COLOURS, "Multiplier") +
    xlab("Multiplier") +
    ylab("Percent Error") +
    ggtitle("BAF Integerness Error by Multiplier")
  
  # Save boxplot
  boxplot_path <- "baf_punishment_boxplot.png"
  message("Saving boxplot to ", boxplot_path)
  ggsave(boxplot_path, boxplot, width = 10, height = 8)
  
  # Save punishment data
  punishment_path <-"punishment_data.csv.gz"
  message("Saving punishment data to ", punishment_path)
  vroom::vroom_write(punishment, punishment_path)
}

message("BAF punishment score calculation complete.")