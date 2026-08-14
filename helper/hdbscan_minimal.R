mixedsort <- gtools::mixedsort

CNV_COLOURS <- structure(
  c(
    "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA", "#D4B9DA"
  ),
  names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+", "11")
)

CNV_COLOURS_BLACK_HALPOID <- structure(
  c(
    "#3182BD", "#000000", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA", "#D4B9DA"
  ),
  names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+", "11")
)


# scCNAS_colors <- c(
#   `0|0` = "#3182BD",
#   `1|0` = "#9ECAE1",
#   `1|1` = "#CCCCCC",
#   `2|0` = "#666666",
#   `2|1` = "#FDCC8A",
#   `3|0` = "#FEE2BC",
#   `2|2` = "#FC8D59",
#   `3|1` = "#FDC1A4",
#   `4|0` = "#FB590E",
#   `5` = "#E34A33",
#   `6` = "#B30000",
#   `7` = "#980043",
#   `8` = "#DD1C77",
#   `9` = "#DF65B0",
#   `10` = "#C994C7",
#   `11+` = "#D4B9DA"
# )
# 
scCNAS_colors <- structure (
  c("#3182BD","#9ECAE1","#CCCCCC","#666666","#FDCC8A","#FDCC8A","#FEE2BC",
    "#FC8D59", "#FDC1A4","#FDC1A4", "#FB590E", "#E34A33", "#B30000", "#980043",
    "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA", "#D4B9DA"),
  names = c("0|0", "1|0","1|1","2|0","2|1","1|2", "3|0",
            "2|2","3|1","1|3","4|0", "5", "6", "7",
            "8", "9", "10", "11+", "11")
)

scCNphase_colors <- c(
  `A-Hom` = "#56941E",
  `B-Hom` = "#471871",
  `A-Gained` = "#94C773",
  `B-Gained` = "#7B52AE",
  `Balanced` = "#d5d5d4"
)



CNV_COLOURS_BLACK_FIX <- structure(
  c("#00FF00", "#000000",
    "#3182BD", "#9ECAE1","#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA", "#D4B9DA"
  ),
  names=c("-3", "-2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+", "11")
)

CNV_COLOURS_BLUE_FIX <- structure(
  c("#000000", "#FF0000", "#0000FF",
    "#3182BD", "#9ECAE1","#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA", "#D4B9DA"
  ),
  names=c("-4", "-3", "-2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+", "11")
)

CNV_COLOURS_BLUE_FOREGROUND_ONLY <- structure(
  c("#000000", "#FF0000", "#0000FF", 
    "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", 
    "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF"
  ),
  names=c("-4", "-3", "-2", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+", "11")
)

# CNV_COLOURS_FOREGROUND_DELTA <- structure(
#   c("#00006b","#1e1eff", "#6666FF", "#CCCCCC", "#ffa590", "#ff8164", "#ff6242", "#ff4122", "#fb3b1e", "#ed3419","#df2c14", "#c61a09", "#861a09" ),
#   names=c("-3", "-2","-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
# )

# CNV_COLOURS_FOREGROUND_DELTA <- structure(
#   c("#084594","#2171b5", "#6baed6", "#DDDDDD", "#fb6a4a", "#cb181d", "#99000d"),
#   names=c("-3", "-2","-1", "0", "1", "2", "3+")
# )

CNV_COLOURS_FOREGROUND_DELTA <- structure(
  c("#053061","#2166ac", "#4393c3", "#DDDDDD", "#fb6a4a", "#cb181d", "#67001f"),
  names=c("FG Loss -3", "FG Loss -2","FG Loss -1", "Median", "FG Gain 1", "FG Gain 2", "FG Gain 3+")
)

CNV_COLOURS_FOREGROUND<-structure(
c("#053061","#2166ac","#4393c3","#DDDDDD","#fb6a4a","#cb181d","#67001f"),
names=c("-3","-2","-1","0","1","2","3")
)

CNV_COLOURS_FOREGROUND_small <-structure(
  c("#053061","#DDDDDD","#cb181d"),
  names=c("loss","0","gain")
)

CNV_COLOURS_FOREGROUND_EDGE <- structure(
  c("#fa9fb5", "#c51b8a" , "#DDDDDD", "#6baed6", "#084594"),
  names=c("- Gain-INTER", "- Gain-PART", "median","- Loss-INTER", "- Loss-PART")
)


CNV_COLOURS_BLACK_SIGNAL_STATE <- structure(
  c("#000000",
    "#3182BD", "#9ECAE1","#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA", "#D4B9DA"
  ),
  names=c(".Signal", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+", "11")
)

make_corrupt_tree_heatmap <- function(tree_ggplot) {
  tree_annot_func = AnnotationFunction(
    fun=function(index) {
      pushViewport(viewport(height=1))
      grid.draw(ggplotGrob(tree_ggplot)$grobs[[5]])
      popViewport()
    },
    var_import=list(tree_ggplot=tree_ggplot),
    width=unit(4, "cm"),
    which="row"
  )
  tree_annot <- HeatmapAnnotation(
    tree=tree_annot_func, which="row", show_annotation_name=FALSE
  )
  
  n_cells <- sum(tree_ggplot$data$isTip)
  tree_hm <- Heatmap(matrix(nc=0, nr=n_cells), left_annotation=tree_annot)
  
  return(tree_hm)
}

format_clones <- function(clones, tree_plot_dat) {
  tree_cells <- get_ordered_cell_ids(tree_plot_dat)
  clones <- merge(clones, data.frame(cell_id=tree_cells), all=TRUE)
  clones[is.na(clones$clone_id), "clone_id"] <- "None"
  
  clone_counts <- clones %>% group_by(clone_id) %>% summarise(count=n())
  
  for(i in 1:nrow(clone_counts)) {
    clone <- unlist(clone_counts[i, "clone_id"], use.names=FALSE)
    clone_count <- unlist(clone_counts[i, "count"], use.names=FALSE)
    clone_label <- paste0(clone, " (", clone_count, ")")
    clones[clones$clone_id == clone, "clone_label"] <- clone_label
  }
  
  rownames(clones) <- clones$cell_id
  clones <- clones[tree_cells, ]
  
  return(clones)
}

get_ordered_cell_ids <- function(tree_plot_dat) {
  return(rev(arrange(tree_plot_dat[tree_plot_dat$isTip, ], y)$label))
}

get_clone_label_pos <- function(clones) {
  clone_label_pos <- list()
  for(clone in unique(clones$clone_id)) {
    if(!grepl("None", clone)) {
      clone_idx <- which(clones$clone_id == clone)
      clone_idx <- find_largest_contiguous_group(clone_idx)
      clone_label_pos[[as.character(clone)]] <-
        as.integer(round(mean(clone_idx)))
    }
  }
  return(clone_label_pos)
}

find_largest_contiguous_group <- function(x) {
  starts <- c(1, which(diff(x) != 1 & diff(x) != 0) + 1)
  ends <- c(starts[-1] - 1, length(x))
  largest <- which.max(ends - starts + 1)
  return(x[starts[largest]:ends[largest]])
}

format_copynumber_values <- function(copynumber) {
  copynumber[copynumber > 11] <- 11
  for(col in colnames(copynumber)) {
    values <- as.character(copynumber[, col])
    values[values == "11"] <- "11+"
    copynumber[, col] <- values
  }
  return(copynumber)
}


space_copynumber_columns <- function(copynumber, spacer_cols) {
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  spacer <- as.data.frame(matrix(
    data=NA, nrow=nrow(copynumber), ncol=spacer_cols
  ))
  chrom_copynumber_dfs <- list()
  for(chrom in mixedsort(unique(chroms))) {
    chrom_copynumber <- copynumber[, chroms == chrom, drop=FALSE]
    chrom_copynumber_dfs <- c(chrom_copynumber_dfs, list(chrom_copynumber))
    chrom_copynumber_dfs <- c(chrom_copynumber_dfs, list(spacer))
  }
  chrom_copynumber_dfs[length(chrom_copynumber_dfs)] <- NULL
  copynumber <- do.call(cbind, chrom_copynumber_dfs)
  
  return(copynumber)
}


format_copynumber <- function(copynumber, tree_plot_dat = NULL, spacer_cols=20) {
  if (!("chr" %in% colnames(copynumber))) {
    loci <- sapply(rownames(copynumber), strsplit, "_")
    copynumber$chr <- unname(sapply(loci, '[[', 1))
    copynumber$start <- as.numeric(unname(sapply(loci, '[[', 2)))
    copynumber$end <- as.numeric(unname(sapply(loci, '[[', 3)))
    copynumber$width <- (copynumber$end - copynumber$start + 1)
  }
  copynumber$chr <- gsub("chr", "", copynumber$chr)
  # copynumber <- arrange(copynumber, as.numeric(chr), chr, start)
  copynumber <- arrange(copynumber, as.character(chr), chr, start)
  
  rownames(copynumber) <- paste0(
    copynumber$chr, ":", copynumber$start, ":", copynumber$end
  )
  copynumber <- subset(copynumber, select=-c(chr, start, end, width))
  copynumber <- as.data.frame(t(copynumber))
  
  if (!is.null(tree_plot_dat)) {
  	  copynumber <- copynumber[get_ordered_cell_ids(tree_plot_dat), ]
  }
  
  copynumber <- format_copynumber_values(copynumber)
  copynumber <- space_copynumber_columns(copynumber, spacer_cols)
  
  return(copynumber)
}


dlp_long_to_wide <- function(raw) {
	raw$pos <- with(raw, paste0(chr, "_", start, "_", end))
	long <- raw[, c("pos", "state", "cell_id")]
	wide <- tidyr::pivot_wider(long, values_from = state, names_from = cell_id)
	wide <- tibble::column_to_rownames(wide, "pos")
	return(wide)
	
	# raw <- raw %>%
	#   dplyr::mutate(chr_desc=paste0(chr,'_', start, '_',end)) %>%
	#   dplyr::select(chr_desc, state, cell_id, clone_id) %>%
	#   tidyr::pivot_wider(values_from = 'state', names_from = 'cell_id') %>%
	#   tibble::column_to_rownames("chr_desc")
}

reads_to_raw <- function(raw) {
  raw <- raw %>%
    dplyr::mutate(chr_desc=paste0(chr,'_', start, '_',end)) %>%
    dplyr::select(chr_desc, state, cell_id, clone_id) %>%
    tidyr::pivot_wider(values_from = 'state', names_from = 'cell_id') %>%
    tibble::column_to_rownames("chr_desc")
  return(raw)
}

make_clone_palette <- function(levels) {
    clone_names <- sort(levels)
    n <- length(clone_names)
    # if (n < 3) {
    #     pal <- brewer.pal(3, "Set1")
    #     pal <- head(pal, n)
    # } else {
    #     pal <- brewer.pal(n, "Set1")
    #     if (n > 9) {
    #         pal <- brewer.pal(9, "Set1")
    #         pal <- c(pal, brewer.pal(n - 9, "Set2"))
    #         if (n > 17) {
    #             pal <- rainbow(n)
    #         }
    #     }
    # }

    if(n<=8){
      pal <- brewer.pal(n, "Dark2")
    }else{
      pal <- colorRampPalette(brewer.pal(8, "Set2"))(n)
    }
    pal <- as.character(pal)
    names(pal) <- clone_names

    pal <- pal[levels]
    return(pal)
}


make_clone_paletteNew <- function(levels) {
  clone_names <- sort(levels)
  n <- length(clone_names)
  
  if (n < 3) {
    pal <- brewer.pal(3, "Set1")
    pal <- head(pal, n)
  } else if (n <= 9) {
    pal <- brewer.pal(n, "Set1")
  } else if (n <= 17) {
    pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
  } else {
    pal <- rainbow(n)
  }
  
  pal <- as.character(pal)
  names(pal) <- clone_names
  
  return(pal)
}


make_annot_chr_labels <- function(cnv_matrix) {
  output <- c()
  chroms <- sapply(strsplit(colnames(cnv_matrix), ":"), function(x) x[[1]])
  uniq_chroms <- c(as.character(1:22), "X", "Y")
  for(chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    output[[chrom]] <- as.integer(round(min(chrom_idx)))
  }
  return(output)
}

