# Converts Reads into Segments, by default will ignore low complexity regions
readsToSegs <- function(reads, na.replace = -1) {
  longseg <- reads[, c("chr", "start", "end", "state", "copy", "cell_id")]
  longseg <- longseg[order(longseg$cell_id, longseg$chr, longseg$start), ]
  
  longseg$state[is.na(longseg$state)] <- na.replace
  longseg$copy[is.na(longseg$state)] <- na.replace
  
  longseg$encode <- as.numeric(longseg$chr) * 100 + longseg$state
  longrle <- rle(longseg$encode)
  longseg$rle <- rep(1:length(longrle$lengths), longrle$lengths)
  medseg <- group_by(.data =longseg, chr, state, cell_id, rle) %>% summarise(start = min(start), end = max(end), median_copy = median(copy, na.rm = TRUE))
  shortseg <- medseg[, c("chr", "start", "end", "state", "median_copy", "cell_id")]
  shortseg$chr <- factor(shortseg$chr, c(1:22, "X", "Y"))
  shortseg <- shortseg[with(shortseg, order(cell_id, chr, start, end)), ]
  shortseg$width <- shortseg$end - shortseg$start + 1
  return(shortseg)
}

infer_missing_data <- function(reads, columns = c("state", "A", "B"), dummy_cell){
  cells <- unique(reads$cell_id)
  reads_with_blanks <- left_join(dummy_cell, reads)
  reads_with_blanks <- select(reads_with_blanks, chr, start, end, state, cell_id) %>% arrange(cell_id, chr, start, end) ## the dummy_cell of state thats real
  blank_bins <- subset(reads_with_blanks, is.na(state)) ## this has the NA bins that dont have value
  
  number_of_cells <- length(unique(reads$cell_id))
  
  # bins of all cells 
  boac <- bind_rows(rep(list(blank_bins), number_of_cells))
  boac$cell_id <- rep(cells, rep(nrow(blank_bins), length(cells)))
  boac <- select(boac, chr, start, end, cell_id)
  boac$enders <- boac$end[1:nrow(boac)] - c(boac$start[2:nrow(boac)], NA)+1
  boac$starter <- boac$start[1:nrow(boac)] - c(NA, boac$end[1:(nrow(boac)-1)])-1
  boac$starter[is.na(boac$starter)] = 1
  boac$enders[is.na(boac$enders)] = 1
  
  ender_rows <- subset(boac, enders != 0)
  start_rows <- subset(boac, starter != 0)
  
  blank_bins_segs = data.frame(chr = start_rows$chr, start = start_rows$start, end = ender_rows$end, cell_id = start_rows$cell_id)
  
  # ADD ALL COLUMNS YOU WANT HERE
  basecols <- c("cell_id", "chr", "start", "end")
  allcols <- c(basecols, columns)
  reads_with_data <- reads[, allcols]
  # reads_with_data <- select(reads, cell_id, chr, start, end, state, A, B, state_AS_phased, state_AS, LOH, phase, state_phase, state_BAF)
  reads_with_data <- rename(reads_with_data, bin_start = start, bin_end = end)
  
  # INFERENCE HERE.  WE ASSUME BLANK BINS TAKE VALUE OF PRIOR BINS
  blank_bins_segs$bin_end <- blank_bins_segs$start - 1
  segs_data_inferred <- left_join(blank_bins_segs, reads_with_data)
  good_data_inferred <- subset(segs_data_inferred, !is.na(state))
  
  # Fixing blank chromosal tips, with no prior bins, we instead of bins coming after
  blank_chr_start <- subset(segs_data_inferred, is.na(state)) %>% select(chr, start, end, cell_id)
  blank_chr_start$bin_start <- blank_chr_start$end + 1
  blank_chr_inferred <- left_join(blank_chr_start, reads_with_data)
  
  blank_all_inferred <- bind_rows(blank_chr_inferred, good_data_inferred)
  blank_all_inferred$bin_end <- NULL
  blank_all_inferred$bin_start <- NULL
  
  reads_all_inferred <- segsToReads2(blank_all_inferred)
  
  reads_return <- bind_rows(reads_all_inferred, reads)
  reads_return <- reads_return %>% arrange(cell_id, chr, start, end)
  return(reads_return)
}

segsToReads2 <- function(segments, binsize = 5e5) {
  collist <- colnames(segments)
  basecols <- c("chr", "start", "end")
  if(!all(basecols %in% collist)) {
    stop("Expecting columns `chr`, `start` and `end`")
  } else {
    datcols <- setdiff(collist, basecols)
  }
  
  expanded <- mapply(seq, from = segments$start, to = segments$end, by = binsize)
  seg_lengths <- sapply(expanded, length)
  chr <- rep(segments$chr, seg_lengths)
  start <- unlist(expanded)
  end <- unlist(expanded) + 5e5 - 1
  width <- binsize
  new_reads <- data.frame(chr, start, end, width)
  for (datcol in datcols) {
    tmp <- rep(segments[, datcol], seg_lengths)
    new_reads[, datcol] <- tmp
  }
  
  return(new_reads)
}

ReadsToSegsCopyless <- function(reads, na.replace = -1) {
  longseg <- reads[, c("chr", "start", "end", "state", "cell_id")]
  longseg <- longseg[order(longseg$cell_id, longseg$chr, longseg$start), ]
  
  longseg$state[is.na(longseg$state)] <- na.replace
  longseg$copy[is.na(longseg$state)] <- na.replace
  
  longseg$encode <- as.numeric(longseg$chr) * 100 + as.numeric(longseg$state) # can drop as.numeric() on the longseg$state
  longrle <- rle(longseg$encode)
  longseg$rle <- rep(1:length(longrle$lengths), longrle$lengths)
  medseg <- group_by(.data =longseg, chr, state, cell_id, rle) %>% summarise(start = min(start), end = max(end), median_state = median(state, na.rm = TRUE))
  shortseg <- medseg[, c("chr", "start", "end", "state", "median_state", "cell_id")]
  shortseg$chr <- factor(shortseg$chr, c(1:22, "X", "Y"))
  shortseg <- shortseg[with(shortseg, order(cell_id, chr, start, end)), ]
  shortseg$width <- shortseg$end - shortseg$start + 1
  
  return(shortseg)
}

ReadsToSegsCopylessTMP <- function(reads, na.replace = -1) {
  longseg <- reads[, c("chr", "start", "end", "state", "cell_id")]
  longseg <- longseg[order(longseg$cell_id, longseg$chr, longseg$start), ]
  
  longseg$state[is.na(longseg$state)] <- na.replace
  longseg$copy[is.na(longseg$state)] <- na.replace
  
  longseg$encode <- as.numeric(longseg$chr) * 100 + as.numeric(longseg$state) # can drop as.numeric() on the longseg$state
  longrle <- rle(longseg$encode)
  longseg$rle <- rep(1:length(longrle$lengths), longrle$lengths)
  medseg <- group_by(.data =longseg, chr, state, cell_id, rle) %>% summarise(start = min(start), end = max(end))
  shortseg <- medseg[, c("chr", "start", "end", "state", "cell_id")]
  shortseg$chr <- factor(shortseg$chr, c(1:22, "X", "Y"))
  shortseg <- shortseg[with(shortseg, order(cell_id, chr, start, end)), ]
  shortseg$width <- shortseg$end - shortseg$start + 1
  
  return(shortseg)
}

segsToReads <- function(segments, binsize = 5e5) {
  expanded <- mapply(seq, from = segments$start, to = segments$end, by = binsize)
  seg_lengths <- sapply(expanded, length)
  chr <- rep(segments$chr, seg_lengths)
  start <- unlist(expanded)
  end <- unlist(expanded) + 5e5 - 1
  width <- binsize
  state <- rep(segments$state, seg_lengths)
  cell_id <- rep(segments$cell_id, seg_lengths)
  # edge <- rep(segments$edge2, seg_lengths)
  new_reads <- data.frame(chr, start, end, width, state, cell_id)
  return(new_reads)
}

mediccSegsToReads <- function(segments, binsize = 5e5) {
  expanded <- mapply(seq, from = segments$start, to = segments$end, by = binsize)
  seg_lengths <- sapply(expanded, length)
  chr <- rep(segments$chr, seg_lengths)
  start <- unlist(expanded)
  end <- unlist(expanded) + 5e5 - 1
  width <- binsize
  state <- rep(segments$state, seg_lengths)
  cell_id <- rep(segments$cell_id, seg_lengths)
  cn_a <- rep(segments$cn_a, seg_lengths)
  cn_b <- rep(segments$cn_b, seg_lengths)
  is_loss <- rep(segments$is_loss, seg_lengths)
  is_gain <- rep(segments$is_gain, seg_lengths)
  is_wgd <- rep(segments$is_wgd, seg_lengths)
  is_normal <- rep(segments$is_normal, seg_lengths)
  is_clonal <- rep(segments$is_clonal, seg_lengths)
  new_reads <- data.frame(chr, start, end, width, state, cell_id, cn_a, cn_b, is_loss, is_gain, is_wgd, is_normal, is_clonal)  
  return(new_reads)
}

segsToReadsEdge <- function(segments, binsize = 5e5) {
  expanded <- mapply(seq, from = segments$start, to = segments$end, by = binsize)
  seg_lengths <- sapply(expanded, length)
  chr <- rep(segments$chr, seg_lengths)
  start <- unlist(expanded)
  end <- unlist(expanded) + 5e5 - 1
  width <- binsize
  state <- rep(segments$state, seg_lengths)
  cell_id <- rep(segments$cell_id, seg_lengths)
  edge <- rep(segments$edge2, seg_lengths)
  new_reads <- data.frame(chr, start, end, width, state, cell_id, edge)
  return(new_reads)
}

segsToReadsExpanded <- function(segments, binsize = 5e5) {
  expanded <- mapply(seq, from = segments$start, to = segments$end, by = binsize)
  seg_lengths <- sapply(expanded, length)
  chr <- rep(segments$chr, seg_lengths)
  start <- unlist(expanded)
  end <- unlist(expanded) + 5e5 - 1
  width <- binsize
  state <- rep(segments$state, seg_lengths)
  cell_id <- rep(segments$cell_id, seg_lengths)
  cn_a <- rep(segments$cn_a, seg_lengths)
  cn_b <- rep(segments$cn_b, seg_lengths)
  new_reads <- data.frame(chr, start, end, width, state, cell_id, cn_a, cn_b)
  return(new_reads)
}

#Annotates segments as 0 (institial), 1 (telomere-bounded), 2 (centromere-bounded), or 3 (whole arm)
annotSegs <- function(segs) {
  
  segs$centromere <- ifelse(segs$state == -1, TRUE, FALSE)
  segs$telostart <- !duplicated(segs[, c("cell_id", "chr")])
  segs$teloend <- !duplicated(segs[, c("cell_id", "chr")], fromLast = TRUE)
  
  segs$centroup <- FALSE
  segs$centroup[intersect(which(segs$centromere) - 1, which(!segs$teloend))] <- TRUE
  segs$centrodown <- FALSE
  segs$centrodown[intersect(which(segs$centromere) + 1, which(!segs$telostart))] <- TRUE
  
  segs$edge <- 0
  segs$edge <- ifelse(segs$telostart | segs$teloend, segs$edge + 1, segs$edge)
  segs$edge <- ifelse(segs$centroup | segs$centrodown, segs$edge + 2, segs$edge)
  
  segs$centromere <- NULL
  segs$telostart <- NULL
  segs$teloend <- NULL
  segs$centroup <- NULL
  segs$centrodown <- NULL
  
  return(segs)
}

annotSegsForeground <- function(segs, na.replace) {
  
  segs$centromere <- ifelse(segs$state == -10, TRUE, FALSE)
  segs$telostart <- !duplicated(segs[, c("cell_id", "chr")])
  segs$teloend <- !duplicated(segs[, c("cell_id", "chr")], fromLast = TRUE)
  
  segs$centroup <- FALSE
  segs$centroup[intersect(which(segs$centromere) - 1, which(!segs$teloend))] <- TRUE
  segs$centrodown <- FALSE
  segs$centrodown[intersect(which(segs$centromere) + 1, which(!segs$telostart))] <- TRUE
  
  segs$edge <- 0
  segs$edge <- ifelse(segs$telostart | segs$teloend, segs$edge + 1, segs$edge)
  segs$edge <- ifelse(segs$centroup | segs$centrodown, segs$edge + 2, segs$edge)
  
  segs$centromere <- NULL
  segs$telostart <- NULL
  segs$teloend <- NULL
  segs$centroup <- NULL
  segs$centrodown <- NULL
  
  return(segs)
}

blackListReads <- function(reads, blfile, replace = NA, binsize = 5e5) {
  
  if (missing(blfile)) {
    stop("Needs segment format blacklist file")
  }
  bl <- read.delim(blfile)
  bl$bins <- bl$width / binsize
  blbig <- data.frame()
  for (i in 1:nrow(bl)) {
    row <- bl[i, ]
    blbit <- data.frame(chr = row$seqnames, start = seq(row$start, row$end, by = 5e5))
    blbig <- bind_rows(blbig, blbit)
  }
  blbig$blacklist <- TRUE
  tmp_cell <- head(reads$cell_id, 1)
  tmp <- subset(reads, cell_id == tmp_cell)
  tmp <- merge(tmp, blbig, all = TRUE)
  tmp$blacklist[is.na(tmp$blacklist)] <- FALSE
  
  blbase <- tmp %>% arrange(cell_id, chr, start) %>% select(blacklist)
  ucell_id <- unique(reads$cell_id)
  reads <- reads %>% arrange(cell_id, chr, start)
  
  if (nrow(reads) / nrow(tmp) == length(ucell_id)) {
    reads$blacklist <- rep(blbase$blacklist, length(ucell_id))
  } else {
    stop("Reads input not of fixed width, function will not work")
  }
  
  if(!missing(replace)) {
    reads$state[reads$blacklist] <- replace
  }
  
  return(reads)
}

factorizeChr <- function(reads) {
  reads$chr <- factor(reads$chr, level = c(1:22, "X", "Y"))
  return(reads)
}

normalizePerChr <- function(reads) {
  chr_stats <- subset(reads, state != -1 & !blacklist) %>% group_by(cell_id, chr) %>% summarise(chr_state_median = median(state), length = max(end))
  reads <- inner_join(reads, select(chr_stats, cell_id, chr, chr_state_median))
  # reads$state <- ifelse(reads$chr_state_median > 0, round(2 * reads$state / reads$chr_state_median), reads$state)
  reads$state <- ifelse(reads$chr_state_median > 0, reads$state - reads$chr_state_median + 2, reads$state)
  return(reads)
}

slowSignal <- function(segs, sig_width = 10e6) {
  segs$small <- segs$width < sig_width
  segs$signal <- NA
  
  for (i in 2:(nrow(segs) - 1)) {
    
    print(paste("On line", i, "of", nrow(segs)))
    # Check each segment
    above <- segs[i - 1, ]
    current <- segs[i, ]
    below <- segs[i + 1, ]
    
    # Skip if blacklisted
    if (current$state == -1) {
      next
    }
    
    # Consider only if small
    if (current$small == TRUE) {
      
      # A signal will be flanked by two equal copy segments on the same chromosome
      if (above$chr == below$chr) {
        if (above$state == below$state) {
          
          # A signal will be one lower than it's adjacent segments
          if (current$state + 1 == above$state) {
            segs[i, "signal"] <- TRUE
          }
          
        } else{
          next
        }
      } else {
        next
      }
    } else {
      next
    }
  }
  return(segs)
}

fastSignal <- function(segs, sig_width = 12.5e6) {
  segs$forwardDelta <- c(0, head(segs$state, -1) - tail(segs$state, -1))
  segs$backwardDelta <- c(tail(segs$state, -1) - head(segs$state, -1), 0)
  segs$deltaSum <- segs$forwardDelta + segs$backwardDelta
  segs$signal <- segs$deltaSum == 2 & segs$width < sig_width & segs$edge == 0
  return(segs)	
}

fastSignalCooler <- function(segs, sig_width = 12.5e6) {
  segs$forwardDelta <- c(0, head(segs$state, -1) - tail(segs$state, -1))
  segs$backwardDelta <- c(tail(segs$state, -1) - head(segs$state, -1), 0)
  segs$deltaSum <- segs$forwardDelta + segs$backwardDelta
  segs$signalLoss <- segs$deltaSum == 2 & segs$width < sig_width & segs$edge == 0
  segs$signalAmp <- segs$deltaSum == -2 & segs$width < sig_width & segs$edge == 0  
  return(segs)	
}

get_chr_id <- function(chr_descs) {
  chrs <- sapply(strsplit(chr_descs, "_"), function(x) {
    return(x[1])
  })
  return(as.character(chrs))
}
get_start_id <- function(chr_descs) {
  chrs <- sapply(strsplit(chr_descs, "_"), function(x) {
    return(x[2])
  })
  return(as.character(chrs))
}

get_end_id <- function(chr_descs) {
  chrs <- sapply(strsplit(chr_descs, "_"), function(x) {
    return(x[3])
  })
  return(as.character(chrs))
}


get_library_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(as.character(labels))
}
get_sample_id <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[1])
  })
  return(as.character(labels))
}
# Turn this into a POINT + SEG plotting function
# # RANDOM
# mets <- subset(dat, SCID == sc)
# pcells <- subset(mets, quality > 0.75 & breakpoints > 40 & note == "cell" & mets$mean_copy < 3)$cell_id

# source("helper/plot_functions.R")
# pick <- subset(reads, cell_id %in% pcells)
# plotDLPgeom_tile(pick)

# psegs <- subset(segs, cell_id %in% pcells)

# EDGE_MAP <- c("INTER", "PART", "PART", "ARM")
# names(EDGE_MAP) <- 0:3
# psegs$edge2 <- EDGE_MAP[as.character(psegs$edge)] 

# ggplot(pick, aes(start, copy, col = as.factor(state))) + geom_point(size = 0.5) + facet_grid(cell_id ~ chr, scales = "free", space = "free_x", switch = "x") + scale_x_continuous(expand = c(0, 0), breaks = NULL) + scale_colour_manual(values = CNV_COLOURS, "CNV") + theme(panel.spacing = unit(0.1, "lines"))
# library(ggnewscale)
# ggplot(pick, aes(start, copy, col = as.factor(state))) + geom_point(size = 0.5) + facet_grid(cell_id ~ chr, scales = "free", space = "free_x", switch = "x") + scale_x_continuous(expand = c(0, 0), breaks = NULL) + scale_colour_manual(values = CNV_COLOURS, "CNV") + theme(panel.spacing = unit(0.1, "lines")) + new_scale_color() + geom_segment(data = subset(psegs, !state %in% c(-1, 2)), aes(start, xend = end, y = state, yend = state, colour = edge2), lwd = 2)

# ggplot(subset(psegs, !state %in% c(-1, 2)), aes(width, fill = edge2)) + geom_histogram(position = "dodge")

# ggplot(subset(psegs, !state %in% c(-1, 2)), aes(width / 5e5, fill = edge2)) + geom_histogram(position = "dodge") + scale_x_log10()