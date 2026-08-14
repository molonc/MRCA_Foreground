library(vroom)
library(dlptools)
library(dplyr)
library(tidyr)
library(data.table)
library(ComplexHeatmap)
ht_opt$message = FALSE
library(ape)
library(tibble)
library(ggplot2)
library(ggtree)
library(optparse)
library(purrr)

 
# use 
# Rscript plot_heatmaps.R --reads reads.csv.gz --metrics signals/test_metrics.csv.gz --tree medicc2_output_tCN/output_segments_all.tsv_final_tree.new --map test_map --allele signals/hscn.csv.gz
# Get the command line arguments

# Define command line options
# option_list <- list(
#   make_option("--reads", type = "character", default = NULL, 
#               help = "Path to reads file (e.g., reads.csv.gz)"),
#   make_option("--metrics", type = "character", default = NULL, 
#               help = "Path to metrics file (e.g., metrics.csv.gz)"),
#   make_option("--tree", type = "character", default = NULL, 
#               help = "Path to MEDICC2 tree file (e.g., medicc2_output_tCN/output_segments_all.tsv_final_tree.new)"),
#   make_option("--map", type = "character", default = NULL, 
#               help = "Path to map file"),
#   make_option("--allele", type = "character", default = FALSE,
#               help = "Path to allele file (e.g., signals/hscn.csv.gz) [default= %default]"),
#   make_option('--bp_file', type = 'character', help = "Breakpoint file")
# )

# # Create the parser object
# opt_parser <- OptionParser(option_list = option_list,
#                            description = "Plot heatmaps using reads, metrics, and MEDICC2 tree data",
#                            usage = "usage: %prog [options]")

# # Parse the command line arguments
# opt <- parse_args(opt_parser)

# # Check if required arguments are provided
# required_args <- c("reads", "metrics", "tree", "map")
# missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]

# if (length(missing_args) > 0) {
#   cat("Error: Missing required arguments:", paste(missing_args, collapse = ", "), "\n")
#   cat("Use --help for more information.\n")
#   quit(status = 1)
# # }


# arg1 <- opt$tree
# arg2 <- opt$reads
# arg3 <- opt$metrics
# arg4 <- opt$map
# ALLELE <- ifelse(is.character(opt$allele) && opt$allele != "FALSE", opt$allele, FALSE)

# # Print configuration for debugging
# cat("Running with parameters:\n")
# cat(paste("Tree file:", arg1, "\n"))
# cat(paste("Reads file:", arg2, "\n"))
# cat(paste("Metrics file:", arg3, "\n"))
# cat(paste("Map file:", arg4, "\n"))
# cat(paste("ALLELE file:", ifelse(ALLELE != FALSE, ALLELE, "Not specified"), "\n"))

reorder_reads_to_medicc_phylo <- function(reads, medicc_phylo) {
######
# I think this is designed to reorder or relabel a data frame of read data (reads) so that its cell_id column follows the order of samples (tips) as they appear in a phylogenetic tree (medicc_phylo)
# -alex
######

  # Get tip labels directly if it's an ape phylo object
  if (inherits(medicc_phylo, "phylo")) {
    tip_ord <- medicc_phylo$tip.label
  } else {
    # Try alternative methods to get the order
    # This depends on what class medic_phylo actually is
    tip_ord <- tryCatch({
      dlptools::cell_id_order_as_plotted(medicc_phylo)
    }, error = function(e) {
      # Fallback method
      warning("Failed to get cell order from tree, check tree format")
      return(NULL)
    })
  }
  
  if (!is.null(tip_ord) && length(tip_ord) > 0) {
    reads$cell_id <- factor(reads$cell_id, levels = rev(tip_ord))
  }
  return(reads)
}

# blackListReads <- function(reads, blfile, replace = NA, binsize = 5e5) {
#   if (missing(blfile)) {
#     stop("Needs segment format blacklist file")
#   }
  
#   if (!file.exists(blfile)) {
#     warning(paste("Blacklist file not found:", blfile))
#     warning("Proceeding without blacklisting")
#     reads$blacklist <- FALSE
#     return(reads)
#   }
  
#   bl <- read.delim(blfile)
#   bl$bins <- bl$width / binsize
#   blbig <- data.frame()
#   for (i in 1:nrow(bl)) {
#     row <- bl[i, ]
#     blbit <- data.frame(chr = row$seqnames, start = seq(row$start, row$end, by = 5e5))
#     blbig <- bind_rows(blbig, blbit)
#   }
#   blbig$blacklist <- TRUE
#   tmp_cell <- head(reads$cell_id, 1)
#   tmp <- subset(reads, cell_id == tmp_cell)
#   tmp <- merge(tmp, blbig, all = TRUE)
#   tmp$blacklist[is.na(tmp$blacklist)] <- FALSE
  
#   blbase <- tmp %>% arrange(cell_id, chr, start) %>% select(blacklist)
#   ucell_id <- unique(reads$cell_id)
#   reads <- reads %>% arrange(cell_id, chr, start)
  
#   if (nrow(reads) / nrow(tmp) == length(ucell_id)) {
#     reads$blacklist <- rep(blbase$blacklist, length(ucell_id))
#   } else {
#     stop("Reads input not of fixed width, function will not work")
#   }
  
#   if(!missing(replace)) {
#     reads$state[reads$blacklist] <- replace
#   }
  
#   return(reads)
# }
# CNV_COLOURS <- dlptools::CNV_COLOURS

# scCNAS_colors <- structure (
#   c("#3182BD","#9ECAE1","#CCCCCC","#666666","#FDCC8A","#FDCC8A","#FEE2BC",
#     "#FC8D59", "#FDC1A4","#FDC1A4", "#FB590E", "#E34A33", "#B30000", "#980043",
#     "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA", "#D4B9DA"),
#   names = c("0|0", "1|0","1|1","2|0","2|1","1|2", "3|0",
#             "2|2","3|1","1|3","4|0", "5", "6", "7",
#             "8", "9", "10", "11+", "11")
# )

# scCNphase_colors <- c(
#   `A-Hom` = "#56941E",
#   `B-Hom` = "#471871",
#   `A-Gained` = "#94C773",
#   `B-Gained` = "#7B52AE",
#   `Balanced` = "#d5d5d4"
# )

# CNV_COLOURS_FOREGROUND<-structure(
# c("#053061","#2166ac","#4393c3","#DDDDDD","#fb6a4a","#cb181d","#67001f", "#ffffff"),
# names=c(-3:3, NA)
# )


# # medicc phylo
# medic_phylo <- ape::read.tree(paste0(arg1))
# medic_phylo <- ape::drop.tip(medic_phylo, "diploid")
# medic_phylo <- as.phylo(medic_phylo)


# # reads df read and order
# reads <- vroom(arg2)
# reads <- reorder_reads_to_medicc_phylo(reads, medic_phylo)
# reads <- reads %>% dplyr::rename(state_CNV = state)
# message(colnames(reads))

# # metadata
# anno_df <- vroom(arg3)
# anno_df <- subset(anno_df, cell_id %in% reads$cell_id)
# anno_df$sample_id <- stringr::str_split(anno_df$cell_id, pattern = "-", simplify = TRUE)[,1]
# anno_df <- select(anno_df, cell_id, sample_id, multiplier, quality) ## cell_id required
# anno_df <- unique(anno_df)

# if(ALLELE != FALSE){
#   # Load the allele file when specified
#   # hscn <- vroom(ALLELE)
#   # reads <- left_join(reads, select(hscn, chr, start, end, cell_id, state_AS, state_phase))
#   cat(paste("Loaded allele file:", ALLELE, "\n"))

#   ## plotting
#   palettes <- list(CNV_COLOURS,
#                   scCNAS_colors, 
#                   scCNphase_colors,
#                   CNV_COLOURS_FOREGROUND)
#   names(palettes) <- c("CNV_COLOURS",
#                        "scCNAS_colors", 
#                        "scCNphase_colors",
#                        "CNV_COLOURS_FOREGROUND")

#   for(p in seq_along(palettes)){
#     q <- palettes[p]
#     print(p)
#     print(q)
#     if(names(q) == "CNV_COLOURS"){
#       reads$state <- reads$state_CNV
#       reads$state <- ifelse(reads$state > 11, 11, reads$state)
#     } else if(names(q) == "scCNAS_colors"){
      
#       reads$state <- reads$state_AS
#       reads$state <- ifelse(reads$state == 12, "11+", reads$state)

#     } else if(names(q) == "scCNphase_colors"){
#       reads$state <- reads$state_phase
#     } else if(names(q) == "CNV_COLOURS_FOREGROUND"){
#       reads$state <- reads$state_fg
#       reads$state <- ifelse(reads$state > 3, 3, ifelse(reads$state < -2, -3, reads$state))
#     }
    
#     cn_states_w <- convert_long_reads_to_wide(reads)

#     current_palette <- palettes[[p]]
    
#     output_file <- paste0(names(q),"_", arg4, ".png")
#     cat(paste("Saving plot to:", output_file, "\n"))
    
    # dlptools::plot_state_hm(
    #   states_df = reads, 
    #   anno_df = anno_df,
    #   phylogeny = medic_phylo,
    #   file_name = output_file,
    #   hm_discrete_colors = current_palette,
    #   state_col = "state"
    # )
#   }
# } else {
#   # Basic plotting without allele data
#   cat("No allele file specified, plotting CNV and foreground heatmaps\n")
  
#   # Plot with CNV_COLOURS
#   reads$state <- reads$state_CNV
#   cn_states_w <- convert_long_reads_to_wide(reads)
  
#   output_file_cnv <-paste0("CNV_COLOURS_", arg4, ".png")
#   cat(paste("Saving CNV plot to:", output_file_cnv, "\n"))
  
#   dlptools::plot_state_hm(
#     states_df = reads, 
#     anno_df = anno_df,
#     phylogeny = medic_phylo,
#     file_name = output_file_cnv,
#     hm_discrete_colors = CNV_COLOURS,
#     state_col = "state"
#   )
  
#   # Plot with FOREGROUND colors
#   if ("state_fg" %in% colnames(reads)) {
#     # Clone the data frame so we don't affect the original state values
#     reads$state <- reads$state_fg
#     reads$state <- ifelse(reads$state > 3, 3, ifelse(reads$state < -2, -3, reads$state))
    
#     output_file_fg <- paste0("CNV_COLOURS_FOREGROUND_", arg4, ".png")
#     cat(paste("Saving foreground plot to:", output_file_fg, "\n"))
    
#     dlptools::plot_state_hm(
#       states_df = reads, 
#       anno_df = anno_df,
#       phylogeny = medic_phylo,
#       file_name = output_file_fg,
#       hm_discrete_colors = CNV_COLOURS_FOREGROUND,
#       state_col = "state"
#     )
#   } else {
#     cat("Foreground state data not available. Skipping foreground plot.\n")
#   }
# }


########### EDGES

  # # Classify gains and losses
  # large_df$fg_type <- case_when(
  #   large_df$state_fg > 0 ~ "GAIN",
  #   large_df$state_fg < 0 ~ "LOSS",
  #   large_df$state_fg == 0 ~ "NEUTRAL"
  # )
  # large_df <- large_df %>%
  #   mutate(edge_simple =
  #     case_when(
  #       edge == "ARM" ~ "ARM",
  #       edge == "CENTRO-BP" ~ "CENTRO-PART",
  #       edge == "BP-CENTRO" ~ "CENTRO-PART",
  #       edge == "BP-TELO" ~ "TELO-PART",
  #       edge == "TELO-BP" ~ "TELO-PART",
  #       edge == "INTER" ~ "INTER",
  #       TRUE ~ NA_character_
  #     )) %>%
  # # Create combined edge and gain/loss type
  # large_df$edge_fg <- paste0(large_df$simple, "_", large_df$fg_type)
  # large_df$edge_fg <- ifelse(large_df$edge_fg == "NA_NEUTRAL", "NEUTRAL", large_df$edge_fg)
  
  # # Print table of edge-foreground state combinations
  # cat("Edge-foreground state combinations:\n")
  # print(table(large_df$edge_fg, large_df$state_fg))
  
  # # Define new color palette
  # new_palette <- structure(
  #   c("#a50026", "#f46d43", "#fdae61", "#fee090", "#f7f7f7", 
  #     "#FFFFFF", "#e0f3f8", "#abd9e9", "#74add1", "#313695"),
  #   names = c("ARM_GAIN", "CENTRO-PART_GAIN", "TELO-PART_GAIN", "INTER_GAIN", "NEUTRAL", 
  #            "NA_NA", "INTER_LOSS", "TELO-PART_LOSS", "CENTRO-PART_LOSS", "ARM_LOSS")
  # )
  
  # # Create final heatmap
  # cat("Creating final heatmap...\n")
  # dlptools::plot_state_hm(
  #   states_df = large_df, 
  #   anno_df = anno_df,
  #   phylogeny = medic_phylo,
  #   file_name = opt$output,
  #   hm_discrete_colors = new_palette,
  #   state_col = "edge_fg"
  # )

hdpsegs_to_reads <- function(hdpsegs, ref, state_col, default = 0) {
  cells  <- unique(hdpsegs$cell_id)
  n_bins <- nrow(ref)

  # Expand ref bins once per cell (cells × bins rows)
  ref_dt <- as.data.table(ref[, c("chr", "start", "end")])
  ref_dt[, chr := as.character(chr)]
  result <- ref_dt[rep(seq_len(n_bins), times = length(cells))]
  result[, cell_id    := rep(cells, each = n_bins)]
  result[, (state_col) := default]

  # Rename segment interval cols to avoid clashing with ref cols
  # chr is coerced to character on both sides: vroom infers its type per-file,
  # so a chromosome set with no "X"/"Y" entries gets guessed as numeric, which
  # breaks the join against ref's (always character) chr column
  seg_dt <- as.data.table(hdpsegs[, c("chr", "cell_id", "start", "end", state_col)])
  seg_dt[, chr := as.character(chr)]
  setnames(seg_dt, c("start", "end"), c("seg_start", "seg_end"))
  setnames(seg_dt, state_col, "seg_state")

  # Single vectorised non-equi join: bin overlaps segment when
  #   bin_start <= seg_end  AND  bin_end >= seg_start
  result[seg_dt,
         on = .(chr, cell_id, start <= seg_end, end >= seg_start),
         (state_col) := i.seg_state]

  setDF(result)
  result
}

make_missing_reads <- function(cell_id, ref) {
  missing_read <- ref
  missing_read$cell_id <- cell_id
  missing_read$state <- "Neutral"
  missing_read %>% select(cell_id, state, start, end, chr)
}



# Function to check if required arguments are provided
check_args <- function(opt) {
  missing_args <- c()
  if (is.null(opt$reads)) missing_args <- c(missing_args, "reads")
  if (is.null(opt$tree)) missing_args <- c(missing_args, "tree")
  if (is.null(opt$metrics)) missing_args <- c(missing_args, "metrics")
  if (is.null(opt$breakpoints)) missing_args <- c(missing_args, "breakpoints")
  if (is.null(opt$library_id)) missing_args <- c(missing_args, "library_id")
  if (is.null(opt$breakpoints)) missing_args <- c(missing_args, "breakpoints")
  if (is.null(opt$library_id)) missing_args <- c(missing_args, "annotated_segs")
  if (is.null(opt$sample_id)) missing_args <- c(missing_args, "sample_Id")

  if (length(missing_args) > 0) {
    stop(paste0("Missing required arguments: ", paste(missing_args, collapse=", ")))
  }
}

validate_inputs <- function(){
  #### VALIDATE INPUTS OF PLOTTING FUNCTIONS
  ### Returns input datframes
  # Define the command line arguments
  option_list <- list(
    make_option(c("-r", "--reads"), type="character", 
                help="Path to reads file that is FINAL (output of previous rule)"),
    # make_option(c("-s", "--segments-anno"), type="character", default=NULL, 
                # help="Path to segments file annotated with breakpoint locations"),
    make_option(c("-t", "--tree"), type="character", default=NULL, 
                help="Path to tree file (e.g. output_segments_all.tsv_final_tree.new)"),
    make_option(c("-m", "--metrics"), type="character", default=NULL, 
                help="Path to metrics file (AT*_metrics.csv.gz)"),
    make_option(c("-o", "--output_dir"), type="character", default=".", 
                help="Output dir for the heatmaps [default=%default]"),
    make_option(c("-b", "--breakpoints"), type="character"),
    make_option(c("-a", "--annotated_segs"), type='character'),
    make_option(c("-l", "--library_id"), type="character"),
    make_option(c("-s", "--sample_id"), type="character"),
    make_option(c("--hdp_segs_A"), type="character", default=NULL,
                help="Path to allele A HDP segs file (optional)"),
    make_option(c("--hdp_segs_B"), type="character", default=NULL,
                help="Path to allele B HDP segs file (optional)"),
    make_option(c("--segs_bp_anno_A"), type="character", default=NULL,
                help="Path to allele A annotated segments file (optional)"),
    make_option(c("--segs_bp_anno_B"), type="character", default=NULL,
                help="Path to allele B annotated segments file (optional)"),
    make_option(c("--ref_bins"), type="character", default=NULL,
                help="Path to reference genomic bins file")
  )

  # Parse the command line arguments
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  check_args(opt)
  return(opt)
}

make_background_CNV_plot <- function(reads, anno_df, medic_phylo, library_id){
  current_palette <- dlptools::CNV_COLOURS
#   current_palette <- c(
#   "0" = "#3182BD",
#   "1" = "#9ECAE1",
#   "2" = "#CCCCCC",
#   "3" = "#FDCC8A",
#   "4" = "#FC8D59",
#   "5" = "#E34A33",
#   "6" = "#B30000",
#   "7" = "#980043",
#   "8" = "#DD1C77",
#   "9" = "#DF65B0",
#   "10" = "#C994C7",
#   "11" = "#D4B9DA",
#   "11+" = "#D4B9DA"
#  )
  reads$state <- ifelse(reads$state > 11, 11, reads$state)
  # cn_states_w <- convert_long_reads_to_wide(reads)
  output_file <- paste0("CNV_ABSOLUTE_",library_id,".png")
  dlptools::plot_state_hm(
      states_df = reads, 
      anno_df = anno_df,
      phylogeny = medic_phylo,
      file_name = output_file,
      hm_discrete_colors = current_palette,
      state_col = "state",
      scale_y = TRUE,
      cell_height=15
    )

}
make_background_CNV_allele_ratio_plot <- function(reads, anno_df, medic_phylo, library_id){
  message("Plotting: Background allele-ratio CNV heatmap")
  current_palette <- structure (
    c("#3182BD","#9ECAE1","#CCCCCC","#666666","#FDCC8A","#FDCC8A","#FEE2BC",
      "#FC8D59", "#FDC1A4","#FDC1A4", "#FB590E", "#E34A33", "#B30000", "#980043",
      "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA", "#D4B9DA"),
    names = c("0|0", "1|0","1|1","2|0","2|1","1|2", "3|0",
              "2|2","3|1","1|3","4|0", "5", "6", "7",
              "8", "9", "10", "11+", "11")
    ) 
  # modify the state column to reflect what we want
  # state_AS is always an "A|B" pair regardless of total copy number, but the
  # palette above only distinguishes allele splits up to a total of 4 (as in
  # the old pipeline where higher states collapsed to a plain integer) - so
  # collapse any higher-total split back down to its total (capped at 11+)
  reads$state <- reads$state_AS
  ab_parts <- strsplit(as.character(reads$state), "\\|")
  ab_totals <- vapply(ab_parts, function(x) {
    if (length(x) != 2) return(NA_real_)
    sum(as.numeric(x))
  }, numeric(1))
  reads$state <- ifelse(
    !is.na(ab_totals) & ab_totals > 4,
    ifelse(ab_totals >= 12, "11+", as.character(ab_totals)),
    reads$state
  )

  output_file <- paste0("CNV_ABSOLUTE_ALLELE_RATIO_",library_id,".png")

  dlptools::plot_state_hm(
    states_df = reads, 
    anno_df = anno_df,
    phylogeny = medic_phylo,
    file_name = output_file,
    hm_discrete_colors = current_palette,
    state_col = "state",
    scale_y = TRUE,
    cell_height=15
  )

}


make_background_CNV_allele_phase_plot <- function(reads, anno_df, medic_phylo, library_id){
  message("Plotting: Background allele-phase CNV heatmap")
  current_palette <- c(
    `A-Hom` = "#56941E",
    `B-Hom` = "#471871",
    `A-Gained` = "#94C773",
    `B-Gained` = "#7B52AE",
    `Balanced` = "#d5d5d4"
    )
  # modify state to reflect what we want to pot
  reads$state <- reads$state_phase

  output_file <- paste0("CNV_ABSOLUTE_ALLELE_PHASE_",library_id,".png")

  dlptools::plot_state_hm(
    states_df = reads, 
    anno_df = anno_df,
    phylogeny = medic_phylo,
    file_name = output_file,
    hm_discrete_colors = current_palette,
    state_col = "state",
    scale_y = TRUE,
    cell_height=15
  )
}

make_fg_CNV_plot <- function(reads, anno_df, medic_phylo, library_id){
  message("Plotting: Foreground CNV heatmap")
  CNV_COLOURS_FOREGROUND<-structure(
  c("#053061","#2166ac","#4393c3","#DDDDDD","#fb6a4a","#cb181d","#67001f", "#ffffff"),
  names=c(-3:3, NA)
  )
  current_palette <- CNV_COLOURS_FOREGROUND
  reads$state <- reads$state_fg
  reads$state <- ifelse(reads$state > 3, 3, ifelse(reads$state < -2, -3, reads$state))

  output_file <- paste0("CNV_FOREGROUND_",library_id,".png")

  dlptools::plot_state_hm(
    states_df = reads,
    anno_df = anno_df,
    phylogeny = medic_phylo,
    file_name = output_file,
    hm_discrete_colors = current_palette,
    state_col = "state",
    scale_y = TRUE,
    cell_height=15
  )
}

make_fg_CNV_allele_A_plot <- function(reads, anno_df, medic_phylo, library_id){
  message("Plotting: Foreground allele-A CNV heatmap")
  CNV_COLOURS_FOREGROUND <- structure(
    c("#053061","#2166ac","#4393c3","#DDDDDD","#fb6a4a","#cb181d","#67001f", "#ffffff"),
    names = c(-3:3, NA)
  )
  reads$state <- reads$state_fg_A
  reads$state <- ifelse(reads$state > 3, 3, ifelse(reads$state < -2, -3, reads$state))

  output_file <- paste0("CNV_FOREGROUND_A_", library_id, ".png")

  dlptools::plot_state_hm(
    states_df = reads,
    anno_df = anno_df,
    phylogeny = medic_phylo,
    file_name = output_file,
    hm_discrete_colors = CNV_COLOURS_FOREGROUND,
    state_col = "state",
    scale_y = TRUE,
    cell_height = 15
  )
}

make_fg_CNV_allele_B_plot <- function(reads, anno_df, medic_phylo, library_id){
  message("Plotting: Foreground allele-B CNV heatmap")
  CNV_COLOURS_FOREGROUND <- structure(
    c("#053061","#2166ac","#4393c3","#DDDDDD","#fb6a4a","#cb181d","#67001f", "#ffffff"),
    names = c(-3:3, NA)
  )
  reads$state <- reads$state_fg_B
  reads$state <- ifelse(reads$state > 3, 3, ifelse(reads$state < -2, -3, reads$state))

  output_file <- paste0("CNV_FOREGROUND_B_", library_id, ".png")

  dlptools::plot_state_hm(
    states_df = reads,
    anno_df = anno_df,
    phylogeny = medic_phylo,
    file_name = output_file,
    hm_discrete_colors = CNV_COLOURS_FOREGROUND,
    state_col = "state",
    scale_y = TRUE,
    cell_height = 15
  )
}


make_fg_CNV_edge_plot_allele <- function(anno_segs_allele, anno_df, medic_phylo, library_id, allele_letter, ref_bins) {
  message("Plotting: Foreground edge-type heatmap (allele ", allele_letter, ")")
  if (is.null(anno_segs_allele) || nrow(anno_segs_allele) == 0) {
    message("Allele ", allele_letter, " annotated segs is empty, skipping edge plot")
    return(invisible(NULL))
  }

  current_palette <- structure(
    c("#a50026", "#fee090", "#56941E", "#313695", "#DDDDDD"),
    names = c("ARM", "CENTRO", "TELO", "INTER", "NA")
  )

  segs_w_edges <- anno_segs_allele %>%
    mutate(edge_simple = case_when(
      edge == "ARM"       ~ "ARM",
      edge == "CENTRO-BP" ~ "CENTRO",
      edge == "BP-CENTRO" ~ "CENTRO",
      edge == "BP-TELO"   ~ "TELO",
      edge == "TELO-BP"   ~ "TELO",
      edge == "INTER"     ~ "INTER",
      TRUE ~ "NA"
    ))
  segs_w_edges$edge <- segs_w_edges$edge_simple

  annotated_reads <- hdpsegs_to_reads(segs_w_edges, ref_bins, state_col = "edge", default = "NA") %>%
    select(chr, cell_id, start, end, edge) %>%
    filter(cell_id %in% anno_df$cell_id)

  output_file <- paste0("CNV_FOREGROUND_EDGES_", allele_letter, "_", library_id, ".png")

  dlptools::plot_state_hm(
    states_df = annotated_reads,
    anno_df = anno_df,
    phylogeny = medic_phylo,
    file_name = output_file,
    hm_discrete_colors = current_palette,
    state_col = "edge",
    scale_y = TRUE,
    cell_height = 15
  )
}

adjust_edges <- function(df) {
  df %>%
    mutate(edge = ifelse(is.na(edge), NA_character_, as.character(edge))) %>% 
    group_by(cell_id, chr) %>%
    arrange(start, .by_group = TRUE) %>%
    group_modify(~{
      dat <- .x
      n <- nrow(dat)
      
      # if only one row, just return it unchanged
      if (n == 1) return(dat)
      
      for (i in seq_len(n)) {
        if (!is.na(dat$edge[i]) && dat$edge[i] == "BREAKPOINT_LOC") {
          # previous row (if exists)
          if (i > 1) {
            dat$end[i - 1] <- min(dat$end[i - 1], dat$start[i] - 1)
          }
          # next row (if exists)
          if (i < n) {
            dat$start[i + 1] <- max(dat$start[i + 1], dat$end[i] + 1)
          }
        }
      }
      
      # drop any intervals that collapse (start > end)
      dat %>% filter(start <= end)
    }) %>%
    ungroup()
}

make_fg_breakpoint_plot <- function(
    reads, anno_df, medic_phylo, library_id, breakpoints, anno_segs, ref_bins,
    bp_padding = 200000, outdir = "breakpoint_cells"
) {
  message("Plotting: Per-cell breakpoint location heatmaps")
  # Color palette
  current_palette <- structure(
    c("#000000", "#fee090", "#56941E", "#313695", "#DDDDDD", "#a50026"),
    names = c("ARM", "CENTRO", "TELO", "INTER", "NA", "BREAKPOINT_LOC")
  )

  # Collapse edge types
  segs_w_edges <- anno_segs %>%
    mutate(edge_simple = case_when(
      edge == "ARM" ~ "ARM",
      edge %in% c("CENTRO-BP", "BP-CENTRO") ~ "CENTRO",
      edge %in% c("TELO-BP", "BP-TELO") ~ "TELO",
      edge == "INTER" ~ "INTER",
      TRUE ~ NA_character_
    )) %>%
    mutate(edge = edge_simple)

  # Take 5 cells for plotting
  five_cells <- head(unique(reads$cell_id), 5)
  message("Selecting 5 Cells")
  reads        <- filter(reads, cell_id %in% five_cells)
  breakpoints  <- filter(breakpoints, cell_id %in% five_cells)
  segs_w_edges <- filter(segs_w_edges, cell_id %in% five_cells)

  # Expand segs into reads
  message("Expanding Breakpoint Annotations from segs to reads")
  annotated_reads <- hdpsegs_to_reads(segs_w_edges, ref_bins, state_col = "edge", default = "NA") %>%
    select(chr, cell_id, start, end, edge)

  # Add artificial breakpoint bins
  message("Creating Breakpoint Bins (±", bp_padding, ")")
  breakpoint_segs <- breakpoints %>%
    mutate(
      start = breakpoint - bp_padding,
      end   = breakpoint + bp_padding-1,
      edge  = "BREAKPOINT_LOC",
      state = 0
    ) %>%
    select(chr, cell_id, start, end, state, edge)
  annotated_reads <- bind_rows(annotated_reads, breakpoint_segs)

  # Adjust bins to not overlap breakpoints
  message("Adjusting Bins to Have Breakpoints")
  annotated_reads <- adjust_edges(annotated_reads)
   # Ensure output dir exists
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  # Plot per cell
  for (cell in five_cells) {
    cell_df <- annotated_reads %>% 
      filter(cell_id == cell) %>%
      arrange(chr, start, end) %>%
      distinct(chr, start, end, .keep_all = TRUE)
    message("Plotting cell ", cell, " (", nrow(cell_df), " bins)")
    output_file <- file.path(outdir, paste0("CNV_FOREGROUND_EDGES_", cell, ".png"))
    plot_edge_heatmap(cell_df,current_palette, output_file)
  }
}
make_fg_CNV_edge_plot <- function(reads, anno_df, medic_phylo, library_id, breakpoints, anno_segs, ref_bins){
  message("Plotting: Foreground edge-type heatmap (all cells)")
  # Color palette
  current_palette <- structure(
    c("#a50026", "#fee090", "#56941E", "#313695","#DDDDDD"),
    names = c("ARM", "CENTRO", "TELO", "INTER", "NA")
  )


  # Extract only segements with edges from the breakpoints file to
  # get the boundedness of our CNVs
  segs_w_edges <- anno_segs %>%
    mutate(edge_simple = 
      case_when(
        edge == "ARM" ~ "ARM",
        edge == "CENTRO-BP" ~ "CENTRO",
        edge == "BP-CENTRO" ~ "CENTRO",
        edge == "BP-TELO" ~ "TELO",
        edge == "TELO-BP" ~ "TELO",
        edge == "INTER" ~ "INTER",
        TRUE ~ "NA"
      )) 
  segs_w_edges$edge <- segs_w_edges$edge_simple
  # Expand the breakpoint file
  print("Expanding Breakpoint Annotations from segs to reads")
  annotated_reads <- hdpsegs_to_reads(segs_w_edges, ref_bins, state_col = "edge", default = "NA") %>%
    select(chr, cell_id, start, end, edge) %>%
    filter(cell_id %in% anno_df$cell_id)
  # print("Arrging cells to be order of phylo")
  output_file <- paste0("CNV_FOREGROUND_EDGES_",library_id,".png")
  dlptools::plot_state_hm(
    states_df = annotated_reads, 
    anno_df = anno_df,
    phylogeny = medic_phylo,
    file_name = output_file,
    hm_discrete_colors = current_palette,
    state_col = "edge",
    scale_y = TRUE,
    cell_height=15
  )
}



make_fg_edge_status_plot <- function(hdp_segs_path, ref_bins, anno_df, medic_phylo, library_id, allele = NULL) {
  message("Plotting: Foreground edge-status (gain/loss × edge type) heatmap", if (!is.null(allele)) paste0(" (allele ", allele, ")") else "")
  current_palette <- c(
    "Arm-Gain"    = "#ff9693", "Telo-Gain"    = "#ba7aff",
    "Centro-Gain" = "#a4fffe", "Inter-Gain"   = "#c0ff80",
    "Arm-Loss"    = "#b10005", "Telo-Loss"    = "#43008E",
    "Centro-Loss" = "#009e9c", "Inter-Loss"   = "#498f00",
    "Neutral"     = "#C8C8C9"
  )

  segs <- vroom::vroom(hdp_segs_path, show_col_types = FALSE, progress = FALSE)
  segs$state_fg <- ifelse(segs$state_fg > 0, "Gain",
                     ifelse(segs$state_fg < 0, "Loss", "Neutral"))
  segs$edge <- ifelse(segs$edge == "ARM", "Arm",
                 ifelse(segs$edge == "TELO-PART", "Telo",
                   ifelse(segs$edge == "CENTRO-PART", "Centro", "Inter")))
  segs$state <- paste0(segs$edge, "-", segs$state_fg)
  segs <- segs %>% select(cell_id, state, start, end, chr) %>% distinct()

  reads <- hdpsegs_to_reads(segs, ref_bins, state_col = "state", default = "Neutral") %>%
    select(cell_id, state, start, end, chr) %>%
    distinct()

  lo_cell_ids <- c(medic_phylo$tip.label, medic_phylo$node.label)

  missing_cells <- lo_cell_ids[!lo_cell_ids %in% reads$cell_id]
  if (length(missing_cells) > 0) {
    message(paste("WARNING:", length(missing_cells), "cells were missing, adding diploid cells"))
    reads <- rbind(reads, do.call(rbind, lapply(missing_cells, make_missing_reads, ref = ref_bins)))
    missing_anno <- data.frame(
      cell_id    = missing_cells,
      sample_id  = NA_character_,
      multiplier = NA_real_,
      quality    = NA_real_
    )
    anno_df <- bind_rows(anno_df, missing_anno) %>% distinct(cell_id, .keep_all = TRUE)
  }

  reads <- reorder_reads_to_medicc_phylo(reads, medic_phylo) %>%
    filter(chr != "Y") %>%
    filter(as.character(cell_id) %in% medic_phylo$tip.label) %>%
    rename(Foreground_State = state)
  anno_df <- anno_df %>% filter(cell_id %in% medic_phylo$tip.label)

  suffix <- if (!is.null(allele)) paste0("_", allele) else ""
  output_file <- paste0("CNV_FOREGROUND_EDGE_STATUS", suffix, "_", library_id, ".png")

  if (nrow(anno_df) != length(unique(reads$cell_id))) {
    print(setdiff(anno_df$cell_id, unique(reads$cell_id)))
    stop("anno_df cell IDs do not match reads cell IDs")
  }
  dlptools::plot_state_hm(
    states_df = reads,
    anno_df   = anno_df,
    phylogeny = medic_phylo,
    file_name = output_file,
    hm_discrete_colors = current_palette,
    state_col  = "Foreground_State",
    scale_y    = TRUE,
    cell_height = 15
  )
}

plot_edge_heatmap <- function(df, current_palette, output_file) {
  # ensure edge is character with "NA"
  df <- df %>%
    group_by(chr) %>%
    arrange(start, .by_group = TRUE) %>%      # make sure it's ordered
    mutate(
      edge  = ifelse(is.na(edge), "NA", edge),
      start_rel = (start - min(start)) / (max(start) - min(start)),  # 0-1 relative
      y     = factor(1)
    ) %>%
    ungroup()
  df <- df %>%
    mutate(
      edge  = ifelse(is.na(edge), "NA", edge),
      chr   = factor(chr),
      start = factor(start),
      y     = factor(1)  # one strip per facet
    )

  # plot using ggplot with facets for each chromosome
  p <- ggplot(df, aes(x = start, y = y, fill = edge)) +
    geom_tile(height = 1) +   # makes it a horizontal strip
    scale_fill_manual(values = current_palette) +
    facet_wrap(~ chr, scales = "free_x", ncol = 1, strip.position = "right") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      panel.spacing = unit(-1, "lines"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(fill = "Edge", x = "Genomic bins", title = paste(unique(df$cell_id),": Breakpoint locations"))
  # ggtitle()
  # save as PNG
  ggsave(output_file, p, width = 12, height = 6, dpi = 300, bg='white')
  message("Faceted heatmap saved to: ", output_file)
}

main <-function(){
  print("Entering Main Exe")
  # Validate inputs
  opt <- validate_inputs()
  tree_path <- opt$tree
  reads_path <- opt$reads
  # segs_path <-opt$segments-anno
  metrics_path<- opt$metrics
  breakpoints_path <- opt$breakpoints
  annotated_segs_path <- opt$annotated_segs
  library_id <- opt$library_id
  sample_id <- opt$sample_id
  output_dir <- opt$output_dir
  hdp_segs_A_path <- opt$hdp_segs_A
  hdp_segs_B_path <- opt$hdp_segs_B
  segs_bp_anno_A_path <- opt$segs_bp_anno_A
  segs_bp_anno_B_path <- opt$segs_bp_anno_B


  # Load the tree
  medic_phylo <- ape::read.tree(tree_path)
  medic_phylo <- ape::drop.tip(medic_phylo, "diploid")
  medic_phylo <- as.phylo(medic_phylo)

  # reads df read and order
  reads <- vroom::vroom(reads_path,
    col_types = vroom::cols(
      chr = col_character(),
      start = col_double(),
      end = col_double(),
      width = col_double(),
      cell_id = col_character(),
      state = col_double(),
      state_fg = col_double(),
      state_fg_A = col_double(),
      state_fg_B = col_double(),
      state_AS = col_character(),
      state_phase = col_character(),
      LOH = col_character(),
      ptel = col_double(),
      cen_start = col_double(),
      cen = col_double(),
      cen_end = col_double(),
      qtel = col_double(),
      is_p_telo = col_logical(),
      is_q_telo = col_logical(),
      is_centro = col_logical(),
      blacklist = col_logical()
    )
  )
  # reads <- subset(reads, chr != "X")
  reads <- subset(reads, chr != "Y")
  reads <- subset(reads, !is.na(state))
  reads <- reorder_reads_to_medicc_phylo(reads, medic_phylo) # re-order cells to match tree
  # DEBUG: subsample to 100 cells for fast diagnostic runs
  # reads <- reads %>% filter(cell_id %in% head(unique(cell_id), 20))
  # reads <- reads %>% dplyr::rename(state_CNV = state)

  ### Load metrics (for quality)
  anno_df <- vroom(metrics_path)
  
  anno_df <- subset(anno_df, cell_id %in% reads$cell_id)
  anno_df$sample_id <- stringr::str_split(anno_df$cell_id, pattern = "-", simplify = TRUE)[,1]
  anno_df <- select(anno_df, cell_id, sample_id, multiplier, quality) ## cell_id required
  anno_df <- anno_df %>% distinct(cell_id, .keep_all = TRUE)

  # Pad anno_df with NA rows for any reads cells absent from metrics
  reads_cells <- as.character(unique(reads$cell_id))
  missing_from_metrics <- setdiff(reads_cells, anno_df$cell_id)
  if (length(missing_from_metrics) > 0) {
    message(sprintf("[DEBUG main] %d reads cells not in metrics — adding NA rows to anno_df", length(missing_from_metrics)))
    anno_df <- bind_rows(anno_df, data.frame(
      cell_id    = missing_from_metrics,
      sample_id  = NA_character_,
      multiplier = NA_real_,
      quality    = NA_real_
    ))
  }

  message(sprintf("[DEBUG main] reads unique cells: %d | anno_df rows: %d | tree tips: %d",
    length(unique(reads$cell_id)), nrow(anno_df), length(medic_phylo$tip.label)))
  print(head(anno_df))


  # Load reference bins (shared by all hdpsegs_to_reads calls)
  ref_bins <- NULL
  if (!is.null(opt$ref_bins)) {
    ref_bins <- vroom::vroom(opt$ref_bins, show_col_types = FALSE, progress = FALSE) %>%
      filter(!is_p_telo, !is_q_telo) %>%
      filter(!(chr %in% c("13","14","15","21","22") & is_centro)) %>%
      filter(!blacklist)
  }

  print("Plotting")
  # Make edge-foreground plots
  print("Loading breakpoints")
  breakpoints <- vroom::vroom(breakpoints_path, col_types = vroom::cols(
      chr = col_character(),
      cell_id = col_character(),
      start = col_double(),
      end = col_double(),
      state = col_double(),
      edge = col_character()
  ))
  print("Loaded annot segs")
  annotated_segs <- vroom::vroom(annotated_segs_path, 
                                col_types = 
                                  vroom::cols(
                                    edge=vroom::col_character(),
                                    breakpoint=vroom::col_double(),
                                    median_loh = vroom::col_double(),
                                    median_abs = vroom::col_double(),
                                    median_wgd = vroom::col_double(),
                                    median_bgCNV = vroom::col_double()
                                  )
  )
                                    
  # Plot the breakpoints with respect to the edge for select cells
  make_fg_breakpoint_plot(reads, anno_df, medic_phylo, library_id, breakpoints, annotated_segs, ref_bins)
  make_fg_CNV_edge_plot(reads, anno_df, medic_phylo, library_id, breakpoints, annotated_segs, ref_bins)

  # Make background plots
  make_background_CNV_plot(reads, anno_df, medic_phylo, library_id)


  # Make allele plots
  if ("state_AS" %in% colnames(reads)){
    make_background_CNV_allele_ratio_plot(reads, anno_df, medic_phylo, library_id)
    make_background_CNV_allele_phase_plot(reads, anno_df, medic_phylo, library_id)
  }

  # Make foreground plots
  make_fg_CNV_plot(reads, anno_df, medic_phylo, library_id)

  # Make allele-specific foreground plots when allele foreground columns are present
  if ("state_fg_A" %in% colnames(reads) && !all(is.na(reads$state_fg_A))) {
    make_fg_CNV_allele_A_plot(reads, anno_df, medic_phylo, library_id)
  }
  if ("state_fg_B" %in% colnames(reads) && !all(is.na(reads$state_fg_B))) {
    make_fg_CNV_allele_B_plot(reads, anno_df, medic_phylo, library_id)
  }

  # Make allele-specific edge plots using per-allele annotated segs (full genome coverage)
  load_anno_segs <- function(path) {
    if (is.null(path)) return(NULL)
    df <- vroom::vroom(path, col_types = vroom::cols(
      chr = col_character(), edge = col_character(), breakpoint = col_double()
    ))
    if (nrow(df) == 0) return(NULL)
    df
  }
  annotated_segs_A <- load_anno_segs(segs_bp_anno_A_path)
  annotated_segs_B <- load_anno_segs(segs_bp_anno_B_path)
  make_fg_CNV_edge_plot_allele(annotated_segs_A, anno_df, medic_phylo, library_id, "A", ref_bins)
  make_fg_CNV_edge_plot_allele(annotated_segs_B, anno_df, medic_phylo, library_id, "B", ref_bins)

  # Make fg edge status heatmap (edge type x gain/loss) for each allele
  if (!is.null(hdp_segs_A_path) && !is.null(ref_bins)) {
    make_fg_edge_status_plot(hdp_segs_A_path, ref_bins, anno_df, medic_phylo, library_id, allele = "A")
  }
  if (!is.null(hdp_segs_B_path) && !is.null(ref_bins)) {
    make_fg_edge_status_plot(hdp_segs_B_path, ref_bins, anno_df, medic_phylo, library_id, allele = "B")
  }
}


if (sys.nframe() == 0) {
  main()
}