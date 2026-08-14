#!/usr/bin/env Rscript

# Install required packages if missing
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("vroom", quietly = TRUE)) {
  install.packages("vroom", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "https://cran.rstudio.com/")
}

library(optparse)
library(dplyr)
library(tidyr)
library(vroom)
library(data.table)

# Define command line options
option_list <- list(
  make_option("--segs", type = "character", default = NULL,
              help = "Path to segments file (e.g., test.csv.gz)"),
  
  make_option("--dummy", type = "character", default = NULL,
              help = "Path to dummy cell file (e.g., dummy_cell.csv)"),
  
  make_option("--profiles", type = "character", default = NULL,
              help = "Path to MEDICC2 profiles file"),
  
  make_option("--output", type = "character", default = NULL,
              help = "Output file prefix"),
  
  make_option("--tcn_bool", type = "logical", default = TRUE,
              help = "True if medicc2 allele is provided. Otherwise False"),

  make_option("--breakpoint_functions", type = "character", default = ".",
            help = "path to breakpoint functions"),

  make_option("--hscn_path", type = "character", default = NULL,
            help = "path to hscn"),

  make_option("--cen_info", type = "character", default = NULL,
            help = "Path to chromosome arm centromere info file (chr, ptel, cen_start, cen_end, qtel)")

)

# Create parser object
opt_parser <- OptionParser(option_list = option_list,
                          description = "Assemble full dataframe from segmentation data",
                          usage = "usage: %prog [options] or\n       %prog --segs SEGS_FILE --dummy DUMMY_FILE --profiles PROFILES_FILE --output OUTPUT_PREFIX")

opt <- parse_args(opt_parser)

required <- c("segs","output","cen_info")
missing <- required[sapply(required, function(x) is.null(opt[[x]]))]

if (length(missing) > 0) {
  cat("Error: Missing required arguments:", paste(missing, collapse=", "), "\n")
  print_help(opt_parser)
  quit(status=1)
}

segs      <- opt$segs
output    <- opt$output
tcn_bool    <- opt$tcn_bool
breakpoint_functions <- opt$breakpoint_functions
hscn_path            <- opt$hscn_path
cen_info_path        <- opt$cen_info


# Print configuration
cat("Running with parameters:\n")
cat(paste("Segments file:", opt$segs, "\n"))
cat(paste("Output prefix:", opt$output, "\n"))
cat(paste("Centromere info file:", opt$cen_info, "\n"))
if (!is.null(tcn_bool)) {
  cat(paste("Using TCN", tcn_bool, "\n"))
}

## FUNCTIONS

segsToReads <- function(segments, binsize = 5e5) {
  expanded <- mapply(seq, from = segments$start, to = segments$end, by = binsize)
  seg_lengths <- sapply(expanded, length)
  chr <- rep(segments$chr, seg_lengths)
  start <- unlist(expanded)
  end <- unlist(expanded) + 5e5 - 1
  width <- binsize
  state <- rep(segments$state, seg_lengths)
  state_fg <- rep(segments$state_fg, seg_lengths)
  state_fg_A <- if ("state_fg_A" %in% colnames(segments)) rep(segments$state_fg_A, seg_lengths) else rep(NA_real_, sum(seg_lengths))
  state_fg_B <- if ("state_fg_B" %in% colnames(segments)) rep(segments$state_fg_B, seg_lengths) else rep(NA_real_, sum(seg_lengths))
  state_AS <- rep(segments$state_AS, seg_lengths)
  state_phase <- rep(segments$state_phase, seg_lengths)
  LOH <- rep(segments$LOH, seg_lengths)
  is_wgd <- rep(segments$is_wgd, seg_lengths)

  cell_id <- rep(segments$cell_id, seg_lengths)
  new_reads <- data.frame(chr, start, end, width, state, cell_id, state_fg, state_fg_A, state_fg_B, state_AS, state_phase, LOH, is_wgd)
  return(new_reads)
}

segsToReads_tCN <- function(segments, binsize = 5e5) {
  expanded <- mapply(seq, from = segments$start, to = segments$end, by = binsize)
  seg_lengths <- sapply(expanded, length)
  chr <- rep(segments$chr, seg_lengths)
  start <- unlist(expanded)
  end <- unlist(expanded) + 5e5 - 1
  width <- binsize
  state <- rep(segments$state, seg_lengths)
  state_fg <- rep(segments$state_fg, seg_lengths)
  state_fg_A <- if ("state_fg_A" %in% colnames(segments)) rep(segments$state_fg_A, seg_lengths) else rep(NA_real_, sum(seg_lengths))
  state_fg_B <- if ("state_fg_B" %in% colnames(segments)) rep(segments$state_fg_B, seg_lengths) else rep(NA_real_, sum(seg_lengths))
  cell_id <- rep(segments$cell_id, seg_lengths)
  new_reads <- data.frame(chr, start, end, width, state, cell_id, state_fg, state_fg_A, state_fg_B)
  return(new_reads)
}

assembleSegsMedicc2 <- function(segs){
  if (all(c("cn_a", "cn_b") %in% colnames(segs))) {
    segs$state <- segs$cn_a + segs$cn_b
  } else {
    # If either cn_a or cn_b is missing, set state to NA or retain the existing state
    segs$state <- segs$state
  }

  segs <- segs %>% dplyr::rename(chr = chrom) %>% dplyr::rename(cell_id = sample_id)
  return(segs)
}

infer_missing_data <- function(reads, columns = c("state", "A", "B"), dummy_cell, allele){

  cells <- unique(reads$cell_id)
  reads$chr <- as.character(reads$chr)
  dummy_cell$chr <- as.character(dummy_cell$chr)
  reads_with_blanks <- left_join(dummy_cell, reads)

  reads_with_blanks <- select(reads_with_blanks, chr, start, end, state, cell_id) %>% arrange(cell_id, chr, start, end) ## the dummy_cell of state thats real
  
  blank_bins <- subset(reads_with_blanks, is.na(state)) ## this has the NA bins that dont have value # I feel this is the problem

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
  # print(head(reads_with_data))
  # print(dim(reads_with_data))

  # INFERENCE HERE.  WE ASSUME BLANK BINS TAKE VALUE OF PRIOR BINS
  blank_bins_segs$bin_end <- blank_bins_segs$start - 1

  blank_bins_segs$chr <- as.character(blank_bins_segs$chr)
  reads_with_data$chr <- as.character(reads_with_data$chr)

  segs_data_inferred <- left_join(blank_bins_segs, reads_with_data)
  good_data_inferred <- subset(segs_data_inferred, !is.na(state))
  
  # print(head(good_data_inferred))
  # print(dim(good_data_inferred))

  # Fixing blank chromosal tips, with no prior bins, we instead of bins coming after
  blank_chr_start <- subset(segs_data_inferred, is.na(state)) %>% select(chr, start, end, cell_id)
  blank_chr_start$bin_start <- blank_chr_start$end + 1
  blank_chr_inferred <- left_join(blank_chr_start, reads_with_data)
  
  # print(head(blank_chr_inferred))
  # print(dim(blank_chr_inferred))

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

blackListReads <- function(reads, blfile, replace = NA, binsize = 5e5) {
  if (missing(blfile)) {
    stop("Needs segment format blacklist file")
  }
  
  if (!file.exists(blfile)) {
    warning(paste("Blacklist file not found:", blfile))
    warning("Proceeding without blacklisting")
    reads$blacklist <- FALSE
    return(reads)
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


annotate_arm_boundaries <- function(reads) {
  reads <- as.data.table(reads)
  setorder(reads, cell_id, chr, arm, start)

  reads[, is_p_telo := FALSE]
  reads[, is_q_telo := FALSE]
  reads[, is_centro := FALSE]

  reads[arm == "p", is_p_telo := seq_len(.N) == 1L, by = .(cell_id, chr)]
  reads[arm == "p", is_centro := seq_len(.N) == .N,  by = .(cell_id, chr)]
  reads[arm == "q", is_centro := seq_len(.N) == 1L,  by = .(cell_id, chr)]
  reads[arm == "q", is_q_telo := seq_len(.N) == .N,  by = .(cell_id, chr)]

  return(reads)
}

assign_arm_to_bins <- function(segs, cen_info, columns, binsize = 5e5) {
  cen_info <- as.data.frame(cen_info)
  cen_info$chr     <- as.character(cen_info$chr)
  cen_info$is_acro <- cen_info$cen_start == 0
  cen_info$cen_mid <- ifelse(cen_info$is_acro,
                              as.numeric(cen_info$cen_end),
                              (as.numeric(cen_info$cen_start) + as.numeric(cen_info$cen_end)) / 2)
  cen_info <- cen_info[, c("chr", "cen_mid", "is_acro")]

  segs$chr <- as.character(segs$chr)
  segs <- left_join(segs, cen_info, by = "chr")

  expanded    <- mapply(seq, from = segs$start, to = segs$end, by = binsize, SIMPLIFY = FALSE)
  seg_lengths <- sapply(expanded, length)

  bin_start <- unlist(expanded)
  bin_end   <- bin_start + binsize - 1

  result <- data.frame(
    chr     = rep(segs$chr,     seg_lengths),
    start   = bin_start,
    end     = bin_end,
    width   = binsize,
    cell_id = rep(segs$cell_id, seg_lengths),
    cen_mid = rep(segs$cen_mid, seg_lengths),
    is_acro = rep(segs$is_acro, seg_lengths)
  )

  for (col in columns) {
    result[[col]] <- if (col %in% colnames(segs)) rep(segs[[col]], seg_lengths) else NA_real_
  }

  # Conservative: bin must lie entirely within one arm
  result$arm <- ifelse(result$end   < result$cen_mid, "p",
                ifelse(result$start > result$cen_mid, "q", NA_character_))

  # Acrocentric chromosomes have no true p-arm; drop bins before their centromere
  result <- result[!(result$arm == "p" & result$is_acro), ]

  result$cen_mid <- NULL
  result$is_acro <- NULL
  result <- result[!is.na(result$arm), ]
  return(result)
}


# main <- function() {

#   # Validate the inputs to the scrip
#   opt <- validate()



# }

# if (sys.nframe()==0){
#   main()
# }

# opt <-c()
# opt$segs <- "/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance2/A108790B/A108790B_seg.csv.gz"
# opt$dummy <- "/Users/ahamazaki/research/packages/medicc2_foreground-master/helper/dummy_cell.csv"
# opt$blacklist = "/Users/ahamazaki/research/packages/medicc2_foreground-master/helper/centro_telo_locs.csv"

# tcn_bool <- "/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance2/A108790B/A108790B_hscn_filtered.csv.gz"


cat("Loading data...\n")
cen_info <- vroom::vroom(cen_info_path)

if (tcn_bool == "true" | tcn_bool == "TRUE" | tcn_bool == TRUE) {
  cat("Reading segments file, not allele aware construction \n")
  segs_foreground <- vroom(opt$segs)

  reads_fg <- assign_arm_to_bins(
    segs_foreground, cen_info,
    columns = c("state", "state_fg", "state_fg_A", "state_fg_B")
  )
  print(head(reads_fg))

  # Merge allele info from separate hscn file
  hscn <- vroom(hscn_path) %>%
    select(c(cell_id, chr, start, end, A, B, state_AS, LOH, state_phase, state_BAF))
  reads <- merge(reads_fg, hscn, by = c("cell_id", "chr", "start", "end"), all.x = TRUE)
  print(head(reads))

} else {
  cat("Reading segments file, Allele Aware...\n")
  segs_foreground <- vroom(opt$segs)

  alleleCN <- segs_foreground %>% rename(A = cn_a, B = cn_b) %>% mutate(state = A + B) %>% as.data.table()
  alleleCN <- alleleCN %>%
    .[, state_AS_phased := paste0(A, "|", B)] %>%
    .[, state_AS := paste0(pmax(state - B, B), "|", pmin(state - B, B))] %>%
    .[, state_min := pmin(A, B)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, phase := c("Balanced", "A", "B")[1 +
      1 * ((B < A)) +
      2 * ((B > A))]] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-Hom", "B-Hom")[1 +
      1 * ((B < A) & (B != 0)) +
      2 * ((B > A) & (A != 0)) +
      3 * ((B < A) & (B == 0)) +
      4 * ((B > A) & (A == 0))]] %>%
    .[order(cell_id, chr, start)] %>%
    .[, state_BAF := round((B / state) / 0.1) * 0.1] %>%
    .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)]

  reads <- assign_arm_to_bins(
    alleleCN, cen_info,
    columns = c("state", "state_fg", "state_fg_A", "state_fg_B", "A", "B",
                "state_AS", "state_phase", "LOH", "is_wgd")
  )
}
print("Reads after inference")
print(head(reads))

reads <- subset(reads, !chr %in% "Y")
reads <- annotate_arm_boundaries(reads)

print("Reads right before save")
# print(head(reads))
# Save results
output_file <- opt$output
cat(paste("Saving results to", output_file, "...\n"))
data.table::fwrite(reads, output_file)

# if (is.logical(reads$state)) {
#   stop("PROBLEM")
# }
# if (output_file == "A108790B_reads_final.csv.gz")
# {stop()
# }

# if (any(reads$state == FALSE, na.rm = TRUE)) {
#   stop("HELP")
# }
# test <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance2/A108746B/A108746B_reads_final.csv.gz")
# any(test$state==FALSE, na.rm=TRUE)

# # stop("TERMINATED")

# cat("Processing complete!\n")


# head(reads)

# segs <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance2/A108790B/A108790B_seg.csv.gz")
# print(segs,Inf)

# test <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance2/A108790B/output_segments_all.tsv_final_cn_profiles.tsv")
# print(segs,Inf)

# segs[segs$cn_a != test$cn_a[1], ]

# # print(segs,Inf)
# # segsToReads_tCN(segs)

# # hsc<- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance2/A108790B/A108790B_hscn_filtered.csv.gz")
# # head(hsc)


# # hhh <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance2/A108790B/A108790B_reads_final.csv.gz")
# # head(hhh)

# test <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance2/A108790B/A108790B_reads_final.csv.gz")
# head(test)

# head(unique(test$cell_id))
# head(unique(reads$cell_id))


