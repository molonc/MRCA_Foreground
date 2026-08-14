library(RColorBrewer)
library(ggtree)
library(signals)
library(dplyr)
library(ComplexHeatmap) # Bioconductor
library(ggplot2)
source("/projects/molonc/aparicio_lab/wdaniels/github/medicc2_foreground/helper/hdbscan_minimal.R")

plotDLPgeom_point <- function(reads) {
	g <- ggplot(reads, aes(start, copy, col = as.factor(state))) + geom_point(size = 0.5) + facet_grid(cell_id ~ chr, scales = "free", space = "free_x", switch = "x") + scale_x_continuous(expand = c(0, 0), breaks = NULL) + scale_colour_manual(values = CNV_COLOURS_BLUE_FIX, "CNV") + theme(panel.spacing = unit(0.1, "lines"))
	return(g)
}

plotDLPgeom_tile <- function(reads) {
	g <- ggplot(reads, aes(start, cell_id, fill = as.factor(state))) + geom_tile() + scale_fill_manual(values = CNV_COLOURS_BLUE_FIX, "CNV") + facet_grid(  ~chr, scales = "free", space = "free", switch = "x") + scale_x_continuous(expand = c(0, 0), breaks = NULL) + theme(panel.spacing = unit(0.1, "lines"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
	)
	return(g)
}

plotDLPgeom_tile_AS_phased <- function(reads) {
  g <- ggplot(reads, aes(start, cell_id, fill = as.factor(state))) + geom_tile() + scale_fill_manual(values = scCNAS_colors, "CNV") + facet_grid(Name  ~chr, scales = "free", space = "free", switch = "x") + scale_x_continuous(expand = c(0, 0), breaks = NULL) + theme(panel.spacing = unit(0.1, "lines"), axis.title.y=element_blank(), axis.ticks.y=element_blank()
  )
  return(g)
}

plotClusterHeatmap <- function(reads, clone_min_cells = 20, seed = 100, blacklist_file = "../data/blacklist_2023.07.17.txt") {

	reads_df <- reads
	reads_df$chr_desc <- paste(reads_df$chr, reads_df$start, reads_df$end, sep = "_")

	# CNV matrix
	copynumber <- dlp_long_to_wide(reads_df)

	bl <- read.delim(blacklist_file)
	bl$bins <- bl$width / 5e5
	blbig <- data.frame()
	for (i in 1:nrow(bl)) {
		row <- bl[i, ]
		blbit <- data.frame(chr = row$seqnames, start = seq(row$start, row$end, by = 5e5))
		blbig <- bind_rows(blbig, blbit)
	}
	blbig$blacklist <- TRUE
	tmp <- subset(reads_df, cell_id == head(reads_df$cell_id, 1))
	tmp <- merge(tmp, blbig, all = TRUE)

	exclude <- subset(tmp, blacklist)$chr_desc
	copynumber[rownames(copynumber) %in% exclude, ] <- NA

	reads_df <- reads_df %>%
		dplyr::filter(cell_id %in% colnames(copynumber) & 
		chr_desc %in% rownames(copynumber))

	# minimum number of cells to consider as a clone
	# pctcells = 0.01
	# ncells <- length(unique(reads_df$cell_id))
	# clone_min_cells <- max(round(pctcells * ncells), 2) 

	# TODO 3: SET YOUR SEED
	set.seed(seed)
# 	print(i)
	#?umap_clustering
	clusters <- signals::umap_clustering(reads_df, 
		minPts = clone_min_cells, 
		field = "copy")
	#### END TODO


	orig_tree <- clusters$tree
	tree <- orig_tree
	clones <- clusters$clustering
	print(unique(clones$clone_id))

	brlen <- NULL

	# make_cell_copynumber_tree_heatmap <- function(tree, copynumber, clones,
	#                                               brlen, grouping_file) {

	# FORMAT TREE, cleans it up, removes 
	# tree <- format_tree(tree, brlen)
	# locus_tips <- grep('locus', tree$tip.label, value=TRUE)
	# tree <- drop.tip(tree, locus_tips)

	# if (!is.null(brlen)) {
	# 	tree <- compute.brlen(tree, brlen)
	# }

	# print('tree_ggplot')
	# tree_ggplot <- make_tree_ggplot(tree, clones)
	tmp <- group_by(clones, clone_id) %>% summarise(clone_members = list(cell_id))
	names(tmp$clone_members) <- tmp$clone_id
	clone_members <- tmp$clone_members

	tree <- ggtree::groupOTU(tree, clone_members)

	clone_levels <- mixedsort(unique(clones$clone_id))
	clone_pal <- make_clone_palette(clone_levels)
	clone_pal["0"] <- "#1B1B1B"

	tree_aes <- aes(x, y, colour=group)

	p <- ggplot(tree, tree_aes) +
		geom_tree(size=0.25) +
		coord_cartesian(expand=FALSE) +
		ylim(0.5, length(tree$tip.label) + 0.5) +
		theme_void() + scale_colour_manual(values=clone_pal)

	tree_ggplot <- p

	old_clones <- clones

	# TREE ENDS HERE

	orig_copynumber <- copynumber
	
	# copynumber <- orig_copynumber
	copynumber <- format_copynumber(copynumber, tree_ggplot$data)

	
	
	clones <- format_clones(clones, tree_ggplot$data) #ADDS CLONE_LABEL HERE
	tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)

	# TEST THIS
	chrom_label_pos <- c()
	chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
	uniq_chroms <- c(as.character(1:22), "X", "Y")
	for(chrom in uniq_chroms) {
		chrom_idx <- which(chroms == chrom)
		chrom_label_pos[[chrom]] <- as.integer(round(min(chrom_idx)))
	}

	bottom_annot <- HeatmapAnnotation(
		chrom_labels =
		anno_mark(
			at = unlist(chrom_label_pos),
			labels = names(chrom_label_pos),
			side = "bottom",
			padding = 0.5, extend = 0.01
		),
		show_annotation_name=FALSE)

	# LEFT ANNOT
	# TODO 4 customize annotation/legend if required (DEFAULT SHOULD BE FAIRLY SAFE MOST OF THE TIME)
	annot_colours <- list()
	# annot_cols <- 3

	library_labels <- stringr::str_extract(rownames(copynumber), "A[0-9]+[A-D]?")
	library_levels <- mixedsort(unique(library_labels))

	# METRICS CHOSEN HERE
	tmp_exp <- subset(select(hq, cell_id, experimental_condition, quality, state_mode, Cell.Line, Drug.Treatment, Treatment..days., Conc.of.Drug, passage), cell_id %in% clones$cell_id)
	clones <- plyr::join(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
	# clones3 <- merge(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
	exp_levels <- mixedsort(unique(clones$Drug.Treatment)) # change this for fun

	n <- length(library_levels) + length(exp_levels)
	six_col <- make_clone_palette(c(library_levels, exp_levels))

	lib_col <- six_col[1:length(library_levels)]
	names(lib_col) <- library_levels
	annot_colours$Library <- lib_col

	clone_levels <- mixedsort(unique(clones$clone_label))
	clo_col <- make_clone_palette(clone_levels)
	annot_colours$Clone <- clo_col

	exp_col <- six_col[(length(library_levels) + 1):(length(library_levels) + length(exp_levels))]
	names(exp_col) <- exp_levels
	annot_colours$Drugs <- exp_col

	CLONE_LABEL_GENERATOR <- function(index) {
	  clone_label_pos <- get_clone_label_pos(clones)
	  y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
	  grid.text(
	    names(clone_label_pos), 0.5, y_pos,
	    just=c("centre", "centre")
	  )
	}

	left_annot <- HeatmapAnnotation(
		Clone=clones$clone_label, 
		clone_label=CLONE_LABEL_GENERATOR,
		Drugs = clones$Drug.Treatment,
		# fixation=fix_labels,
		# sample=sample_labels, 
		# drug=drug_labels, 
		# celltype=celltype_labels, 
		Library=library_labels,
		col=annot_colours, 
		# show_annotation_name=c(TRUE, FALSE, TRUE, TRUE, TRUE,TRUE),  #, TRUE
		show_annotation_name = c(TRUE, FALSE, TRUE, TRUE),
		which="row", 
		annotation_width=unit(rep(0.5, 4), "cm"), # MATCH NUMBER OF ANNOTS
		annotation_legend_param=list(
			# Clone=list(nrow=10),
			# fixation=list(nrow=library_legend_rows),
			# sample=list(nrow=library_legend_rows),
			# drug=list(nrow=library_legend_rows),
			# celltype=list(nrow=library_legend_rows),
			# Groupings=list(nrow=10)
		)
	)

	# TODO 4 END HERE
	# END LEFT ANNOT

	copynumber_hm <- Heatmap(
		name="Copy Number",
		as.matrix(copynumber),
		col=CNV_COLOURS,
		na_col="white",
		show_row_names=FALSE,
		cluster_rows=FALSE,
		cluster_columns=FALSE,
		show_column_names=FALSE,
		bottom_annotation=bottom_annot,
		left_annotation= left_annot,
		heatmap_legend_param=list(nrow=4),
		use_raster=FALSE,
		raster_quality=5
	)

	print('visualize tree')
	h <- tree_hm + copynumber_hm
	# h <- copynumber_hm
	print('Draw tree')

	return(list(h, clusters))

	# # TODO 5, give the proper title name + width (inches) + length (inches)
	# # You can also output PNG instead of PDF
	# pdf("hdbscan_clustering_all_cells2.pdf", 20, 20)
	# print(h)
	# dev.off()

}


plotClusterHeatmapBlack <- function(reads, clone_min_cells = 70, seed = 100, blacklist_file = "../data/blacklist_2023.07.17.txt") {
  
  reads_df <- reads
  reads_df$chr_desc <- paste(reads_df$chr, reads_df$start, reads_df$end, sep = "_")
  
  # CNV matrix
  copynumber <- dlp_long_to_wide(reads_df)
  
  bl <- read.delim(blacklist_file)
  bl$bins <- bl$width / 5e5
  blbig <- data.frame()
  for (i in 1:nrow(bl)) {
    row <- bl[i, ]
    blbit <- data.frame(chr = row$seqnames, start = seq(row$start, row$end, by = 5e5))
    blbig <- bind_rows(blbig, blbit)
  }
  blbig$blacklist <- TRUE
  tmp <- subset(reads_df, cell_id == head(reads_df$cell_id, 1))
  tmp <- merge(tmp, blbig, all = TRUE)
  
  exclude <- subset(tmp, blacklist)$chr_desc
  copynumber[rownames(copynumber) %in% exclude, ] <- NA
  
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% colnames(copynumber) & 
                    chr_desc %in% rownames(copynumber))
  
  # minimum number of cells to consider as a clone
  # pctcells = 0.01
  # ncells <- length(unique(reads_df$cell_id))
  # clone_min_cells <- max(round(pctcells * ncells), 2) 
  
  # TODO 3: SET YOUR SEED
  set.seed(seed)
  # 	print(i)
  #?umap_clustering
  clusters <- signals::umap_clustering(reads_df, 
                                       minPts = clone_min_cells, 
                                       field = "copy")
  #### END TODO
  
  orig_tree <- clusters$tree
  tree <- orig_tree
  clones <- clusters$clustering
  print(unique(clones$clone_id))
  
  brlen <- NULL
  
  # make_cell_copynumber_tree_heatmap <- function(tree, copynumber, clones,
  #                                               brlen, grouping_file) {
  
  # FORMAT TREE, cleans it up, removes 
  # tree <- format_tree(tree, brlen)
  # locus_tips <- grep('locus', tree$tip.label, value=TRUE)
  # tree <- drop.tip(tree, locus_tips)
  
  # if (!is.null(brlen)) {
  # 	tree <- compute.brlen(tree, brlen)
  # }
  
  # print('tree_ggplot')
  # tree_ggplot <- make_tree_ggplot(tree, clones)
  tmp <- group_by(clones, clone_id) %>% summarise(clone_members = list(cell_id))
  names(tmp$clone_members) <- tmp$clone_id
  clone_members <- tmp$clone_members
  
  tree <- ggtree::groupOTU(tree, clone_members)
  
  clone_levels <- mixedsort(unique(clones$clone_id))
  clone_pal <- make_clone_palette(clone_levels)
  clone_pal["0"] <- "#1B1B1B"
  
  tree_aes <- aes(x, y, colour=group)
  
  p <- ggplot(tree, tree_aes) +
    geom_tree(size=0.25) +
    coord_cartesian(expand=FALSE) +
    ylim(0.5, length(tree$tip.label) + 0.5) +
    theme_void() + scale_colour_manual(values=clone_pal)
  
  tree_ggplot <- p
  
  old_clones <- clones
  
  # TREE ENDS HERE
  
  orig_copynumber <- copynumber
  
  # copynumber <- orig_copynumber
  copynumber <- format_copynumber(copynumber, tree_ggplot$data)
  
  
  
  clones <- format_clones(clones, tree_ggplot$data) # ADDS CLONE_LABEL HERE
  tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
  
  # TEST THIS
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  uniq_chroms <- c(as.character(1:22), "X", "Y")
  for(chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    chrom_label_pos[[chrom]] <- as.integer(round(min(chrom_idx)))
  }
  
  bottom_annot <- HeatmapAnnotation(
    chrom_labels =
      anno_mark(
        at = unlist(chrom_label_pos),
        labels = names(chrom_label_pos),
        side = "bottom",
        padding = 0.5, extend = 0.01
      ),
    show_annotation_name=FALSE)
  
  # LEFT ANNOT
  # TODO 4 customize annotation/legend if required (DEFAULT SHOULD BE FAIRLY SAFE MOST OF THE TIME)
  annot_colours <- list()
  # annot_cols <- 3
  
  library_labels <- stringr::str_extract(rownames(copynumber), "A[0-9]+[A-D]?")
  library_levels <- mixedsort(unique(library_labels))
  
  # METRICS CHOSEN HERE
  tmp_exp <- subset(select(hq, cell_id, experimental_condition, quality, state_mode, Cell.Line, Drug.Treatment, Treatment..days., Conc.of.Drug, passage, Name, SCID), cell_id %in% clones$cell_id)
  clones <- plyr::join(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
  # clones3 <- merge(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
  exp_levels <- mixedsort(unique(clones$Drug.Treatment)) # change this for fun
  
  n <- length(library_levels) + length(exp_levels)
  six_col <- make_clone_palette(c(library_levels, exp_levels))
  
  lib_col <- six_col[1:length(library_levels)]
  names(lib_col) <- library_levels
  annot_colours$Library <- lib_col
  
  clone_levels <- mixedsort(unique(clones$clone_label))
  clo_col <- make_clone_palette(clone_levels)
  annot_colours$Clone <- clo_col
  
  exp_col <- six_col[(length(library_levels) + 1):(length(library_levels) + length(exp_levels))]
  names(exp_col) <- exp_levels
  annot_colours$Drugs <- exp_col
  
  CLONE_LABEL_GENERATOR <- function(index) {
    clone_label_pos <- get_clone_label_pos(clones)
    y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
    grid.text(
      names(clone_label_pos), 0.5, y_pos,
      just=c("centre", "centre")
    )
  }
  
  left_annot <- HeatmapAnnotation(
    Clone=clones$clone_label, 
    clone_label=CLONE_LABEL_GENERATOR,
    Drugs = clones$Drug.Treatment,
    # fixation=fix_labels,
    # sample=sample_labels, 
    # drug=drug_labels, 
    # celltype=celltype_labels, 
    Library=library_labels,
    col=annot_colours, 
    # show_annotation_name=c(TRUE, FALSE, TRUE, TRUE, TRUE,TRUE),  #, TRUE
    show_annotation_name = c(TRUE, FALSE, TRUE, TRUE),
    which="row", 
    annotation_width=unit(rep(0.5, 4), "cm"), # MATCH NUMBER OF ANNOTS
    annotation_legend_param=list(
      # Clone=list(nrow=10),
      # fixation=list(nrow=library_legend_rows),
      # sample=list(nrow=library_legend_rows),
      # drug=list(nrow=library_legend_rows),
      # celltype=list(nrow=library_legend_rows),
      # Groupings=list(nrow=10)
    )
  )
  
  # TODO 4 END HERE
  # END LEFT ANNOT
  
  copynumber_hm <- Heatmap(
    name="Copy Number",
    as.matrix(copynumber),
    col=CNV_COLOURS_BLACK_FIX,
    na_col="white",
    show_row_names=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=FALSE,
    bottom_annotation=bottom_annot,
    left_annotation= left_annot,
    heatmap_legend_param=list(nrow=4),
    use_raster=FALSE,
    raster_quality=10
  )
  
  print('visualize tree')
  h <- tree_hm + copynumber_hm
  # h <- copynumber_hm
  print('Draw tree')
  
  return(list(h, clusters))
  
  # # TODO 5, give the proper title name + width (inches) + length (inches)
  # # You can also output PNG instead of PDF
  # pdf("hdbscan_clustering_all_cells2.pdf", 20, 20)
  # print(h)
  # dev.off()
  
}


plotClusterHeatmapLossAmp <- function(reads, clone_min_cells = 70, seed = 100, blacklist_file = "../data/blacklist_2023.07.17.txt") {
  
  reads_df <- reads
  reads_df$chr_desc <- paste(reads_df$chr, reads_df$start, reads_df$end, sep = "_")
  
  # CNV matrix
  copynumber <- dlp_long_to_wide(reads_df)
  
  bl <- read.delim(blacklist_file)
  bl$bins <- bl$width / 5e5
  blbig <- data.frame()
  for (i in 1:nrow(bl)) {
    row <- bl[i, ]
    blbit <- data.frame(chr = row$seqnames, start = seq(row$start, row$end, by = 5e5))
    blbig <- bind_rows(blbig, blbit)
  }
  blbig$blacklist <- TRUE
  tmp <- subset(reads_df, cell_id == head(reads_df$cell_id, 1))
  tmp <- merge(tmp, blbig, all = TRUE)
  
  exclude <- subset(tmp, blacklist)$chr_desc
  copynumber[rownames(copynumber) %in% exclude, ] <- NA
  
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% colnames(copynumber) & 
                    chr_desc %in% rownames(copynumber))
  
  # minimum number of cells to consider as a clone
  # pctcells = 0.01
  # ncells <- length(unique(reads_df$cell_id))
  # clone_min_cells <- max(round(pctcells * ncells), 2) 
  
  # TODO 3: SET YOUR SEED
  set.seed(seed)
  # 	print(i)
  #?umap_clustering
  clusters <- signals::umap_clustering(reads_df, 
                                       minPts = clone_min_cells, 
                                       field = "copy")
  #### END TODO
  
  orig_tree <- clusters$tree
  tree <- orig_tree
  clones <- clusters$clustering
  print(unique(clones$clone_id))
  
  brlen <- NULL
  
  # make_cell_copynumber_tree_heatmap <- function(tree, copynumber, clones,
  #                                               brlen, grouping_file) {
  
  # FORMAT TREE, cleans it up, removes 
  # tree <- format_tree(tree, brlen)
  # locus_tips <- grep('locus', tree$tip.label, value=TRUE)
  # tree <- drop.tip(tree, locus_tips)
  
  # if (!is.null(brlen)) {
  # 	tree <- compute.brlen(tree, brlen)
  # }
  
  # print('tree_ggplot')
  # tree_ggplot <- make_tree_ggplot(tree, clones)
  tmp <- group_by(clones, clone_id) %>% summarise(clone_members = list(cell_id))
  names(tmp$clone_members) <- tmp$clone_id
  clone_members <- tmp$clone_members
  
  tree <- ggtree::groupOTU(tree, clone_members)
  
  clone_levels <- mixedsort(unique(clones$clone_id))
  clone_pal <- make_clone_palette(clone_levels)
  clone_pal["0"] <- "#1B1B1B"
  
  tree_aes <- aes(x, y, colour=group)
  
  p <- ggplot(tree, tree_aes) +
    geom_tree(size=0.25) +
    coord_cartesian(expand=FALSE) +
    ylim(0.5, length(tree$tip.label) + 0.5) +
    theme_void() + scale_colour_manual(values=clone_pal)
  
  tree_ggplot <- p
  
  old_clones <- clones
  
  # TREE ENDS HERE
  
  orig_copynumber <- copynumber
  
  # copynumber <- orig_copynumber
  copynumber <- format_copynumber(copynumber, tree_ggplot$data)
  
  
  
  clones <- format_clones(clones, tree_ggplot$data) # ADDS CLONE_LABEL HERE
  tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
  
  # TEST THIS
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  uniq_chroms <- c(as.character(1:22), "X", "Y")
  for(chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    chrom_label_pos[[chrom]] <- as.integer(round(min(chrom_idx)))
  }
  
  bottom_annot <- HeatmapAnnotation(
    chrom_labels =
      anno_mark(
        at = unlist(chrom_label_pos),
        labels = names(chrom_label_pos),
        side = "bottom",
        padding = 0.5, extend = 0.01
      ),
    show_annotation_name=FALSE)
  
  # LEFT ANNOT
  # TODO 4 customize annotation/legend if required (DEFAULT SHOULD BE FAIRLY SAFE MOST OF THE TIME)
  annot_colours <- list()
  # annot_cols <- 3
  
  library_labels <- stringr::str_extract(rownames(copynumber), "A[0-9]+[A-D]?")
  library_levels <- mixedsort(unique(library_labels))
  
  # METRICS CHOSEN HERE
  tmp_exp <- subset(select(hq, cell_id, experimental_condition, quality, state_mode, Cell.Line, Drug.Treatment, Treatment..days., Conc.of.Drug, passage, Name, SCID), cell_id %in% clones$cell_id)
  clones <- plyr::join(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
  # clones3 <- merge(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
  exp_levels <- mixedsort(unique(clones$Drug.Treatment)) # change this for fun
  
  n <- length(library_levels) + length(exp_levels)
  six_col <- make_clone_palette(c(library_levels, exp_levels))
  
  lib_col <- six_col[1:length(library_levels)]
  names(lib_col) <- library_levels
  annot_colours$Library <- lib_col
  
  clone_levels <- mixedsort(unique(clones$clone_label))
  clo_col <- make_clone_palette(clone_levels)
  annot_colours$Clone <- clo_col
  
  exp_col <- six_col[(length(library_levels) + 1):(length(library_levels) + length(exp_levels))]
  names(exp_col) <- exp_levels
  annot_colours$Drugs <- exp_col
  
  CLONE_LABEL_GENERATOR <- function(index) {
    clone_label_pos <- get_clone_label_pos(clones)
    y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
    grid.text(
      names(clone_label_pos), 0.5, y_pos,
      just=c("centre", "centre")
    )
  }
  
  left_annot <- HeatmapAnnotation(
    Clone=clones$clone_label, 
    clone_label=CLONE_LABEL_GENERATOR,
    Drugs = clones$Drug.Treatment,
    SCID = clones$SCID,
    Name = clones$Name,
    # fixation=fix_labels,
    # sample=sample_labels, 
    # drug=drug_labels, 
    # celltype=celltype_labels, 
    #Library=library_labels,
    #Passage=clones$passage,
    col=annot_colours, 
    # show_annotation_name=c(TRUE, FALSE, TRUE, TRUE, TRUE,TRUE),  #, TRUE
    show_annotation_name = c(TRUE, FALSE, TRUE, TRUE, TRUE),
    which="row", 
    annotation_width=unit(rep(0.5, 4), "cm"), # MATCH NUMBER OF ANNOTS
    annotation_legend_param=list(
      # Clone=list(nrow=10),
      # fixation=list(nrow=library_legend_rows),
      # sample=list(nrow=library_legend_rows),
      # drug=list(nrow=library_legend_rows),
      # celltype=list(nrow=library_legend_rows),
      # Groupings=list(nrow=10)
    )
  )
  
  # TODO 4 END HERE
  # END LEFT ANNOT
  
  copynumber_hm <- Heatmap(
    name="Copy Number",
    as.matrix(copynumber),
    col=CNV_COLOURS_BLACK_FIX,
    na_col="white",
    show_row_names=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=FALSE,
    bottom_annotation=bottom_annot,
    left_annotation= left_annot,
    heatmap_legend_param=list(nrow=4),
    use_raster=FALSE,
    raster_quality=10
  )
  
  print('visualize tree')
  h <- tree_hm + copynumber_hm
  # h <- copynumber_hm
  print('Draw tree')
  
  return(list(h, clusters))
  
  # # TODO 5, give the proper title name + width (inches) + length (inches)
  # # You can also output PNG instead of PDF
  # pdf("hdbscan_clustering_all_cells2.pdf", 20, 20)
  # print(h)
  # dev.off()
  
}


plotClusterHeatmapSignal <- function(reads, clone_min_cells = 70, seed = 100, blacklist_file = "../data/blacklist_2023.07.17.txt") {
  
  reads_df <- reads
  reads_df$chr_desc <- paste(reads_df$chr, reads_df$start, reads_df$end, sep = "_")
  
  # CNV matrix
  copynumber <- dlp_long_to_wide(reads_df)
  
  bl <- read.delim(blacklist_file)
  bl$bins <- bl$width / 5e5
  blbig <- data.frame()
  for (i in 1:nrow(bl)) {
    row <- bl[i, ]
    blbit <- data.frame(chr = row$seqnames, start = seq(row$start, row$end, by = 5e5))
    blbig <- bind_rows(blbig, blbit)
  }
  blbig$blacklist <- TRUE
  tmp <- subset(reads_df, cell_id == head(reads_df$cell_id, 1))
  tmp <- merge(tmp, blbig, all = TRUE)
  
  exclude <- subset(tmp, blacklist)$chr_desc
  copynumber[rownames(copynumber) %in% exclude, ] <- NA
  
  reads_df <- reads_df %>%
    dplyr::filter(cell_id %in% colnames(copynumber) & 
                    chr_desc %in% rownames(copynumber))
  
  # minimum number of cells to consider as a clone
  # pctcells = 0.01
  # ncells <- length(unique(reads_df$cell_id))
  # clone_min_cells <- max(round(pctcells * ncells), 2) 
  
  # TODO 3: SET YOUR SEED
  set.seed(seed)
  # 	print(i)
  #?umap_clustering
  clusters <- signals::umap_clustering(reads_df, 
                                       minPts = clone_min_cells, 
                                       field = "copy")
  #### END TODO
  
  
  orig_tree <- clusters$tree
  tree <- orig_tree
  ### REMOVE CLONE 0. Not as simple as removing. Can do another left join of ID remaining after 0 removed and print it that way but not sure
  clones <- clusters$clustering
  #clones <- subset(clones, clone_id != "0")
  print(unique(clones$clone_id))
  
  brlen <- NULL
  
  # make_cell_copynumber_tree_heatmap <- function(tree, copynumber, clones,
  #                                               brlen, grouping_file) {
  
  # FORMAT TREE, cleans it up, removes 
  # tree <- format_tree(tree, brlen)
  # locus_tips <- grep('locus', tree$tip.label, value=TRUE)
  # tree <- drop.tip(tree, locus_tips)
  
  # if (!is.null(brlen)) {
  # 	tree <- compute.brlen(tree, brlen)
  # }
  
  # print('tree_ggplot')
  # tree_ggplot <- make_tree_ggplot(tree, clones)
  tmp <- group_by(clones, clone_id) %>% summarise(clone_members = list(cell_id))
  names(tmp$clone_members) <- tmp$clone_id
  clone_members <- tmp$clone_members
  
  tree <- ggtree::groupOTU(tree, clone_members)
  
  clone_levels <- mixedsort(unique(clones$clone_id))
  clone_pal <- make_clone_palette(clone_levels)
  clone_pal["0"] <- "#1B1B1B"
  
  tree_aes <- aes(x, y, colour=group)
  
  p <- ggplot(tree, tree_aes) +
    geom_tree(size=0.25) +
    coord_cartesian(expand=FALSE) +
    ylim(0.5, length(tree$tip.label) + 0.5) +
    theme_void() + scale_colour_manual(values=clone_pal)
  
  tree_ggplot <- p
  
  old_clones <- clones
  
  # TREE ENDS HERE
  
  orig_copynumber <- copynumber
  
  # copynumber <- orig_copynumber
  copynumber <- format_copynumber(copynumber, tree_ggplot$data)
  
  
  
  clones <- format_clones(clones, tree_ggplot$data) # ADDS CLONE_LABEL HERE
  tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
  
  # TEST THIS
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  uniq_chroms <- c(as.character(1:22), "X", "Y")
  for(chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    chrom_label_pos[[chrom]] <- as.integer(round(min(chrom_idx)))
  }
  
  bottom_annot <- HeatmapAnnotation(
    chrom_labels =
      anno_mark(
        at = unlist(chrom_label_pos),
        labels = names(chrom_label_pos),
        side = "bottom",
        padding = 0.5, extend = 0.01
      ),
    show_annotation_name=FALSE)
  
  # LEFT ANNOT
  # TODO 4 customize annotation/legend if required (DEFAULT SHOULD BE FAIRLY SAFE MOST OF THE TIME)
  annot_colours <- list()
  # annot_cols <- 3
  
  library_labels <- stringr::str_extract(rownames(copynumber), "A[0-9]+[A-D]?")
  library_levels <- mixedsort(unique(library_labels))
  
  # METRICS CHOSEN HERE
  tmp_exp <- subset(select(hq, cell_id, experimental_condition, quality, state_mode, Cell.Line, Drug.Treatment, Treatment..days., Conc.of.Drug, passage, SCID, Name), cell_id %in% clones$cell_id)
  clones <- plyr::join(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
  # clones3 <- merge(clones, tmp_exp) # DO NOT REORDER THE DATA FRAME
  exp_levels <- mixedsort(unique(clones$Drug.Treatment)) # change this for fun
  scid_levels <- mixedsort(unique(clones$SCID)) # WILLIAM SCID
  passage_levels <- mixedsort(unique(clones$passage)) # WILLIAM SCID
  
  
  #plus 2 is required. weird. 
  n <- length(library_levels) + length(exp_levels) + length(scid_levels) + 2
  six_col <- make_clone_palette(c(library_levels, exp_levels, scid_levels))
  
  lib_col <- six_col[1:length(library_levels)]
  names(lib_col) <- library_levels
  annot_colours$Library <- lib_col
  
  clone_levels <- mixedsort(unique(clones$clone_label))
  clo_col <- make_clone_palette(clone_levels)
  annot_colours$Clone <- clo_col
  
  exp_col <- six_col[(length(library_levels) + 1):(length(library_levels) + length(exp_levels))]
  names(exp_col) <- exp_levels
  annot_colours$Drugs <- exp_col
  
  scid_col <- six_col[1:length(scid_levels)] # WILLIAM SCID
  names(scid_col) <- scid_levels
  annot_colours$SCID <- scid_col
  
  #passage_col <- six_col[1:length(passage_levels)] # WILLIAM SCID
  #names(passage_col) <- passage_levels
  #annot_colours$passage <- passage_col
  
  
  CLONE_LABEL_GENERATOR <- function(index) {
    clone_label_pos <- get_clone_label_pos(clones)
    y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
    grid.text(
      names(clone_label_pos), 0.5, y_pos,
      just=c("centre", "centre")
    )
  }
  
 left_annot <- HeatmapAnnotation(
    Clone=clones$clone_label, 
    clone_label=CLONE_LABEL_GENERATOR,
    Drugs = clones$Drug.Treatment,
    SCID = clones$SCID,
    Name = clones$Name,
    # fixation=fix_labels,
    # sample=sample_labels, 
    # drug=drug_labels, 
    # celltype=celltype_labels, 
    #Library=library_labels,
    #Passage=clones$passage,
    col=annot_colours, 
    # show_annotation_name=c(TRUE, FALSE, TRUE, TRUE, TRUE,TRUE),  #, TRUE
    show_annotation_name = c(TRUE, FALSE, TRUE, TRUE, TRUE),
    which="row", 
    annotation_width=unit(rep(0.5, 5), "cm"), # MATCH NUMBER OF ANNOTS
    annotation_legend_param=list(
      # Clone=list(nrow=10),
      # fixation=list(nrow=library_legend_rows),
      # sample=list(nrow=library_legend_rows),
      # drug=list(nrow=library_legend_rows),
      # celltype=list(nrow=library_legend_rows),
      # Groupings=list(nrow=10)
    )
  )
  
  # TODO 4 END HERE
  # END LEFT ANNOT
  
  copynumber_hm <- Heatmap(
    name="Copy Number",
    as.matrix(copynumber),
    col=CNV_COLOURS_BLACK_FIX,
    #col=CNV_COLOURS_BLACK_SIGNAL_STATE,
    na_col="white",
    show_row_names=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=FALSE,
    bottom_annotation=bottom_annot,
    left_annotation= left_annot,
    heatmap_legend_param=list(nrow=4),
    use_raster=FALSE,
    raster_quality=5
  )
  
  print('visualize tree')
  #h <- tree_hm + copynumber_hm # WILLIAM - IF WANT CLADE THIS IS ORIGINAL LINE TO ADD TO HEATMAP
  h <- copynumber_hm
  # h <- copynumber_hm
  print('Draw tree')
  
  return(list(h, clusters, n)) # ADDED CLONES, DATA RICH OBJECT
  print("returned h, clusters, n, copynumber")
  
  # # TODO 5, give the proper title name + width (inches) + length (inches)
  # # You can also output PNG instead of PDF
  # pdf("hdbscan_clustering_all_cells2.pdf", 20, 20)
  # print(h)
  # dev.off()
  
}
