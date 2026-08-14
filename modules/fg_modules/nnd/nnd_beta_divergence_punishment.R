#!/usr/bin/env Rscript

# ============================================================================
# NND Beta Distribution Divergence Analysis for scWGS Data
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(fitdistrplus)
  library(mclust)
  library(vroom)
})

# Create argument parser
parser <- ArgumentParser(description = "Identify divergent cells using NND and Beta distribution")

parser$add_argument("--hscn", 
                    type = "character", 
                    required = TRUE,
                    help = "Path to hscn file (e.g., AT28926_hscn.csv.gz)")

# parser$add_argument("--output_dir", 
#                     type = "character", 
#                     required = TRUE,
#                     help = "Output directory for results")

parser$add_argument("--sample_id",
                    type = "character",
                    required = FALSE,
                    default = NULL,
                    help = "Sample ID for output naming (extracted from hscn if not provided)")

parser$add_argument("--threshold_percentile",
                    type = "double",
                    default = 0.99,
                    help = "Percentile threshold for divergent cells (default: 0.99)")

parser$add_argument("--distance_method",
                    type = "character",
                    default = "manhattan",
                    choices = c("manhattan", "euclidean"),
                    help = "Distance metric to use (default: manhattan)")

parser$add_argument("--plot_heatmap",
                    action = "store_true",
                    default = FALSE,
                    help = "Generate copy number heatmap of divergent cells")

# Parse arguments
args <- parser$parse_args()
# Extract sample ID from filename if not provided
if (is.null(args$sample_id)) {
  args$sample_id <- gsub("_hscn\\.csv\\.gz$", "", basename(args$hscn))
}
print(args$sample_id)

cat("============================================================\n")
cat("NND Beta Distribution Divergence Analysis\n")
cat("============================================================\n")
cat(sprintf("Input file: %s\n", args$hscn))
# cat(sprintf("Output directory: %s\n", args$output_dir))
cat(sprintf("Sample ID: %s\n", args$sample_id))
cat(sprintf("Threshold percentile: %.1f%%\n", args$threshold_percentile * 100))
cat(sprintf("Distance method: %s\n", args$distance_method))
cat("============================================================\n\n")

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

cat("Loading data...\n")
hscn <- vroom::vroom(args$hscn, show_col_types = FALSE)

# Select relevant columns for distance calculation
min_hscn <- dplyr::select(hscn, cell_id, chr, start, end, state, A, B)

# Create a unique bin identifier
min_hscn <- min_hscn %>%
  dplyr::mutate(bin_id = paste(chr, start, end, sep = "_"))

# ============================================================================
# 2. CREATE COPY NUMBER MATRIX AND CALCULATE DISTANCES
# ============================================================================

cat("Creating copy number matrix...\n")

# Function to create copy number matrix
create_cn_matrix <- function(data) {
  # Create the matrix
  cn_matrix <- data %>%
    dplyr::select(cell_id, bin_id, state) %>%
    tidyr::pivot_wider(names_from = bin_id, 
                       values_from = state, 
                       values_fill = list(state = NA)) %>%
    tibble::column_to_rownames("cell_id") %>%
    as.matrix()
  
  return(cn_matrix)
}

# Function to calculate distance matrix
calculate_distances <- function(cn_matrix, method = "manhattan") {
  # Remove any bins with all NAs
  cn_matrix <- cn_matrix[, colSums(is.na(cn_matrix)) < nrow(cn_matrix)]
  
  # Impute NAs with column median
  for(i in 1:ncol(cn_matrix)) {
    if(any(is.na(cn_matrix[,i]))) {
      cn_matrix[is.na(cn_matrix[,i]), i] <- median(cn_matrix[,i], na.rm = TRUE)
    }
  }
  
  # Calculate distance matrix
  dist_matrix <- dist(cn_matrix, method = method)
  dist_matrix <- as.matrix(dist_matrix)
  
  return(dist_matrix)
}

# Create matrices and calculate distances
cn_matrix <- create_cn_matrix(min_hscn)

cat("Calculating distance matrix...\n")
dist_matrix <- calculate_distances(cn_matrix, method = args$distance_method)

# ============================================================================
# 3. CALCULATE NEAREST NEIGHBOR DISTANCES (NND)
# ============================================================================

cat("Calculating nearest neighbor distances...\n")

# Function to calculate NND
calculate_nnd <- function(dist_matrix) {
  nnd <- apply(dist_matrix, 1, function(x) {
    x_no_self <- x[x > 0]  # Remove self-distance
    if(length(x_no_self) > 0) {
      return(min(x_no_self))
    } else {
      return(NA)
    }
  })
  return(nnd)
}

# Calculate NND values
nnd_values <- calculate_nnd(dist_matrix)
nnd_clean <- nnd_values[!is.na(nnd_values)]

# Normalize NND to [0,1] range for Beta distribution fitting
nnd_normalized <- (nnd_clean - min(nnd_clean)) / (max(nnd_clean) - min(nnd_clean))

# ============================================================================
# 4. FIT BETA DISTRIBUTION AND IDENTIFY DIVERGENT CELLS
# ============================================================================

cat("Fitting Beta distribution...\n")

# Fit Beta distribution to normalized NND
beta_fit <- fitdist(nnd_normalized, "beta", method = "mle")
alpha_est <- beta_fit$estimate["shape1"]
beta_est <- beta_fit$estimate["shape2"]

cat("\nBeta distribution parameters:\n")
cat(sprintf("  Alpha (shape1) = %.3f\n", alpha_est))
cat(sprintf("  Beta (shape2) = %.3f\n", beta_est))

# Statistical threshold approach
threshold_value <- qbeta(args$threshold_percentile, alpha_est, beta_est)
is_divergent <- nnd_normalized > threshold_value

cat(sprintf("\nUsing %.1f%% threshold: %.3f\n", args$threshold_percentile * 100, threshold_value))

# Get divergent cell IDs
divergent_cells <- names(nnd_clean)[is_divergent]

# ============================================================================
# 5. CREATE VISUALIZATIONS
# ============================================================================

cat("\nGenerating plots...\n")

# Create output subdirectory for NND divergence analysis
nnd_output_dir <- ""
# dir.create(nnd_output_dir, showWarnings = FALSE, recursive = TRUE)

# Prepare QQ plot data
n <- length(nnd_normalized)
empirical_quantiles <- sort(nnd_normalized)
theoretical_quantiles <- qbeta(ppoints(n), alpha_est, beta_est)

qq_data <- data.frame(
  theoretical = theoretical_quantiles,
  empirical = empirical_quantiles,
  divergent = is_divergent[order(nnd_normalized)],
  cell_id = names(nnd_clean)[order(nnd_normalized)]
)

# Publication-style QQ plot
p_qq <- ggplot(qq_data, aes(x = theoretical, y = empirical, color = divergent)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "black", size = 1.5) +
  scale_color_manual(values = c("FALSE" = "#3366CC", "TRUE" = "#DC3912"),
                     labels = c("Non-divergent", "Divergent"),
                     name = "") +
  scale_x_continuous(limits = c(0, max(qq_data$theoretical) * 1.1), 
                     breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, max(qq_data$empirical) * 1.1), 
                     breaks = seq(0, 1, 0.2)) +
  labs(title = sprintf("Q-Q Plot of Data vs. Beta Distribution (%s)", args$sample_id),
       x = "Theoretical NND Quantiles (Beta)",
       y = "Empirical NND Quantiles") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        text = element_text(size = 14),
        legend.position = c(0.7, 0.3),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1))

# Save QQ plot
ggsave(filename = sprintf("%s_nnd_beta_qqplot.png", args$sample_id),
       plot = p_qq,
       width = 8,
       height = 7,
       dpi = 300)

# Density plot
# Check for after_stat availability (ggplot2 >= 3.3.0)
density_aes <- if(packageVersion("ggplot2") >= "3.3.0") {
  aes(y = after_stat(density))
} else {
  aes(y = ..density..)
}

p_density <- ggplot(data.frame(nnd = nnd_normalized, divergent = is_divergent), 
                    aes(x = nnd)) +
  geom_histogram(density_aes, bins = 30, alpha = 0.5, fill = "gray") +
  stat_function(fun = function(x) dbeta(x, alpha_est, beta_est), 
                color = "black", size = 1.2) +
  geom_density(color = "blue", size = 1) +
  geom_vline(xintercept = threshold_value, color = "red", linetype = "dashed") +
  annotate("text", x = threshold_value, y = 0, 
           label = sprintf("%.0f%% threshold", args$threshold_percentile * 100), 
           vjust = -0.5, hjust = -0.1, color = "red") +
  labs(title = sprintf("NND Distribution vs Fitted Beta (%s)", args$sample_id),
       x = "Normalized Nearest Neighbor Distance",
       y = "Density") +
  theme_minimal()

# Save density plot
ggsave(filename = sprintf("%s_nnd_density.png", args$sample_id),
       plot = p_density,
       width = 8,
       height = 6,
       dpi = 300)

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================

cat("\nSaving results...\n")

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat(sprintf("Total cells analyzed: %d\n", length(is_divergent)))
cat(sprintf("Number of divergent cells: %d\n", sum(is_divergent)))
cat(sprintf("Percentage divergent: %.1f%%\n", 100 * sum(is_divergent) / length(is_divergent)))

if (sum(is_divergent) > 0) {
  cat(sprintf("\nDivergent cells (top %d):\n", min(10, sum(is_divergent))))
  cat(paste(head(divergent_cells, 10), collapse = "\n"), "\n")
}

# Export divergent cell information
divergent_info <- data.frame(
  cell_id = divergent_cells, 
  nnd_raw = nnd_clean[is_divergent],
  nnd_normalized = nnd_normalized[is_divergent]
) %>%
  dplyr::arrange(desc(nnd_normalized))

# Save divergent cells list

write.csv(divergent_info, 
          file = sprintf("%s_divergent_cells_nnd.csv", args$sample_id), 
          row.names = FALSE)

# Save all NND values for reference
all_nnd_info <- data.frame(
  cell_id = names(nnd_clean),
  nnd_raw = nnd_clean,
  nnd_normalized = nnd_normalized,
  is_divergent = is_divergent
)

write.csv(all_nnd_info,
          file = sprintf("%s_all_cells_nnd.csv", args$sample_id),
          row.names = FALSE)

# Save beta distribution parameters
beta_params <- data.frame(
  sample_id = args$sample_id,
  alpha = alpha_est,
  beta = beta_est,
  threshold_percentile = args$threshold_percentile,
  threshold_value = threshold_value,
  n_cells = length(is_divergent),
  n_divergent = sum(is_divergent),
  pct_divergent = 100 * sum(is_divergent) / length(is_divergent)
)

write.csv(beta_params,
          file = sprintf("%s_beta_parameters.csv", args$sample_id),
          row.names = FALSE)

# ============================================================================
# 7. OPTIONAL: COPY NUMBER HEATMAP OF DIVERGENT CELLS
# ============================================================================

if (args$plot_heatmap && sum(is_divergent) > 0) {
  cat("\nGenerating copy number heatmap...\n")
  
  # Format chromosome levels
  hscn$chr <- factor(hscn$chr, levels = c(1:22, "X"))
  
  # Subset to divergent cells
  divergent_hscn <- hscn %>%
    dplyr::filter(cell_id %in% divergent_cells)
  
  # Order cells by mean copy number state
  mean_state_order <- divergent_hscn %>% 
    dplyr::group_by(cell_id) %>% 
    dplyr::summarize(mean_state = mean(state, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_state)
  
  cell_order <- mean_state_order$cell_id

  divergent_hscn$cell_id <- factor(divergent_hscn$cell_id, levels = cell_order)
  tile_plot <- dlptools::plot_read_bins_basic(divergent_hscn)

    # Save heatmap
  ggsave(filename = sprintf("%s_divergent_cells_heatmap.png", args$sample_id),
         plot = tile_plot,
         width = 12,
         height = max(4, 0.3 * length(divergent_cells)),
         dpi = 300)

  copy_plot <- ggplot(divergent_hscn, aes(start, copy, col = as.factor(state))) + 
                  geom_point(size = 0.5) + 
                  facet_grid(cell_id ~ chr, scales = "free", space = "free_x", switch = "x") + 
                  scale_x_continuous(expand = c(0, 0), breaks = NULL) + 
                  scale_colour_manual(values = dlptools::CNV_COLOURS, "CNV") + 
                  theme(panel.spacing = unit(0.1, "lines"))
  ggsave(filename = sprintf("%s_single_cell_copy_divergent_cells.png", args$sample_id),
       plot = copy_plot,
       width = 12,
       height = max(4, 0.3 * length(divergent_cells)),
       dpi = 300)

}

cat("\nAnalysis complete!\n")
cat(sprintf("Results saved to: %s\n", nnd_output_dir))