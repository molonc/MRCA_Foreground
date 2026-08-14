# use
# Rscript foreground_ancestor.R output_segments_all.tsv_final_tree.new output_segments_all.tsv_final_cn_profiles.tsv NAME

library(vroom)
library(ape)
library(dplyr)
library(tidyr)
library(tibble)

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Arg1 = tree, Arg2 = ancestor and leaf file copy numbers, Arg3 = output file
arg1 <- args[1]
arg2 <- args[2]
arg3 <- args[3]


## functions
assembleSegsMedicc2 <- function(segs){
  if (all(c("cn_a", "cn_b") %in% colnames(segs))) {
    segs$state <- segs$cn_a + segs$cn_b
  } else {
    segs$state <- segs$state
  }
  segs <- segs %>% dplyr::rename(chr = chrom) %>% dplyr::rename(cell_id = sample_id)
  return(segs)
}

getmedicc2MRCA <- function(medicc2_tree, drop_diploid = TRUE){
  tree <- medicc2_tree
  distanceLeafNode <- dist.nodes(tree)
  colnames(distanceLeafNode) <- c(tree$tip.label, tree$node.label)
  rownames(distanceLeafNode) <- c(tree$tip.label, tree$node.label)
  distanceLeafNode <- data.frame(distanceLeafNode)

  colnames(distanceLeafNode) <- gsub("\\.", "-", colnames(distanceLeafNode))
  names <- c(tree$tip.label, tree$node.label[2:length(tree$node.label)])
  distanceLeafNodeValid <- distanceLeafNode[names, names]

  cell_ids <- tree$tip.label
  colnames <- cell_ids
  rownames <- tree$node.label[2:length(tree$node.label)]
  distanceLeafNodeFormat <- distanceLeafNodeValid[rownames, colnames]

  min_values <- apply(distanceLeafNodeFormat, 2, function(column) min(column))
  row_names_with_min_values <- apply(distanceLeafNodeFormat, 2, function(column) rownames(distanceLeafNodeFormat)[which.min(column)])

  result_df <- data.frame(
    cell_id = colnames(distanceLeafNodeFormat),
    parent = row_names_with_min_values,
    distance = min_values,
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- NULL

  if(drop_diploid == TRUE){
    result_df <- subset(result_df, !cell_id %in% "diploid")
  }
  return(result_df)
}

computeOneForeground <- function(medicc2_segs, medicc2_distances, value_col) {
  ## Generic helper: computes foreground (cell - ancestor) for any numeric column value_col.
  ## medicc2_segs must already have chr stripped of "chr" prefix and chr_desc column set.

  tmp <- subset(medicc2_segs, cell_id %in% medicc2_distances$cell_id) %>%
    select(cell_id, all_of(value_col), chr_desc)
  pivot <- pivot_wider(tmp, names_from = "cell_id", values_from = all_of(value_col))

  tmp2 <- subset(medicc2_segs, cell_id %in% medicc2_distances$parent) %>%
    select(cell_id, all_of(value_col), chr_desc)
  pivot2 <- pivot_wider(tmp2, names_from = "cell_id", values_from = all_of(value_col))

  pivot  <- column_to_rownames(pivot,  "chr_desc")
  pivot2 <- column_to_rownames(pivot2, "chr_desc")

  pivot_cells <- pivot[,  medicc2_distances$cell_id]
  pivot_nodes <- pivot2[, medicc2_distances$parent]

  fg_matrix <- pivot_cells - pivot_nodes

  fg_long <- rownames_to_column(fg_matrix, "chr_desc")
  segs_fg <- pivot_longer(fg_long, cols = -1, names_to = "cell_id", values_to = value_col)
  segs_fg <- segs_fg %>% separate(chr_desc, into = c("chr", "start", "end"), sep = "_")
  segs_fg$start <- as.numeric(segs_fg$start)
  segs_fg$end   <- as.numeric(segs_fg$end)
  return(segs_fg)
}

computeForegroundMRCA <- function(medicc2_segs, medicc2_distances){

  if (all(c("cn_a", "cn_b") %in% colnames(medicc2_segs))) {
    medicc2_segs$state <- medicc2_segs$cn_a + medicc2_segs$cn_b
  } else {
    medicc2_segs$state <- medicc2_segs$state
  }

  medicc2_segs$chr <- gsub("chr", "", medicc2_segs$chr)
  medicc2_segs$chr_desc <- paste0(medicc2_segs$chr, "_", medicc2_segs$start, "_", medicc2_segs$end)

  ## Total CN foreground
  segs_foreground <- computeOneForeground(medicc2_segs, medicc2_distances, "state")

  ## Allele-specific foreground when cn_a / cn_b are present (allele-aware medicc2 run)
  if (all(c("cn_a", "cn_b") %in% colnames(medicc2_segs))) {
    fg_a <- computeOneForeground(medicc2_segs, medicc2_distances, "cn_a")
    fg_b <- computeOneForeground(medicc2_segs, medicc2_distances, "cn_b")

    segs_foreground <- segs_foreground %>%
      left_join(
        fg_a %>% select(chr, start, end, cell_id, state_fg_A = cn_a),
        by = c("chr", "start", "end", "cell_id")
      ) %>%
      left_join(
        fg_b %>% select(chr, start, end, cell_id, state_fg_B = cn_b),
        by = c("chr", "start", "end", "cell_id")
      )
  }

  return(segs_foreground)
}


## calculations

## Tree
medic_phylo <- ape::read.tree(paste0(arg1))
medic_phylo <- ape::drop.tip(medic_phylo, "diploid")

segs <- vroom(arg2)
segs <- assembleSegsMedicc2(segs)

## Calculate Foreground with MRCA
mrcaDistance <- getmedicc2MRCA(medic_phylo, drop_diploid = TRUE)

segs_foreground <- computeForegroundMRCA(medicc2_segs = segs, medicc2_distances = mrcaDistance)

## annotate foreground state as "state_fg"
segs_foreground <- segs_foreground %>% dplyr::rename(state_fg = state)

segs$chr <- gsub("chr", "", segs$chr)

segs_foreground <- left_join(segs_foreground, segs, by = join_by(chr, start, end, cell_id))

data.table::fwrite(segs_foreground, arg3)
