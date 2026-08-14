###
# Helper functions for making breakpoint files. see get_bp_file.R
###
library(dplyr)
library(dlptools)
library(purrr)
library(tidyr)


mark_centromeres_and_telomeres <- function(reads, centro_telo_locs){
    ### This function takes in a reads file with all genomic bins. And marks the bins that are centromeres and telomeres with bools
    # If is_q_telo == TRUE, then the bin is a q-telomere
    # If is_p_telo == TRUE, then the bin is a p-telomere
    # If is_centro == TRUE, then the bin is a centromere

    # reads: You need to take a reads file and group it by cell_id and chromosome for this.
    stopifnot(length(unique(reads$cell_id)) == 1) # We only want to do this for one cell at a time
    stopifnot(length(unique(reads$chr)) == 1) # We only want to do this for one cell at a time

    ano_reads <- reads %>%
        left_join(centro_telo_locs, by = "chr") %>%
        mutate(
            is_p_telo = start < ptel,
            is_q_telo = end > qtel,
            is_centro = end > cen_start & start < cen_end,
            is_p_telo = is_p_telo & !is_centro,  # telocentric chroms: centromere at end, classify as centro only
            blacklist = is_p_telo | is_q_telo | is_centro
        )

    # Check the first row
    first_row <- ano_reads[1, ]
    if (!first_row$is_centro && !first_row$is_p_telo) {
        # Create the fake p-telo row
        p_telo_row <- first_row
        p_telo_row$start <- 0
        p_telo_row$end <- 0
        p_telo_row$is_p_telo <- TRUE
        p_telo_row$is_q_telo <- FALSE
        p_telo_row$is_centro <- FALSE
        p_telo_row$blacklist <- TRUE

        # Prepend to ano_reads
        ano_reads <- bind_rows(p_telo_row, ano_reads)
    }
    return(ano_reads)
}

update_state <-function(ano_reads){

    ### Given the p-bound, centro-bound, or t-bound (or inner) nature of the reads, we will update the state of the reads
    # with integers that are indicative of the state
    # ano_reads: reads dataframe annotated with if the segments are p-telo, centro, or q-telo bound

    # If you are a centromere, your state becomes -100
    ano_reads$state[ano_reads$is_centro==TRUE] <- -100
    # if you are atelomere, your state becomes -101 if p bound and 102 if q bound
    ano_reads$state[ano_reads$is_p_telo==TRUE] <- -101
    ano_reads$state[ano_reads$is_q_telo==TRUE] <- -102

    if (any(is.na(ano_reads$state))) {
    # print(subset(ano_reads, is.na(state)), Inf)
    test <- ano_reads %>%select(chr, start, end, state) %>% distinct()
    stop("Diagnostic: NA values found in 'state' column of ano_reads. Expected no missing values.")
    }
    return(ano_reads)
}


non_edge <- function(state){
    if (state %in% c(-100, -101, -102, 0)) {
        return(TRUE)  # Non-edge states
    } else {
        return(FALSE)  # Edge states
    }
}


real_edge <- function(state, prev_state, next_state) {
    if (prev_state == -101 && next_state == -100) {
        return("ARM")
    } else if (next_state == -102 && prev_state == -100) {
        return("ARM")
    } else if (prev_state == -101 && next_state == -102) {
        return("ARM")  # whole chrom single state, no centromere between telos
    } else if (prev_state == -101) {
        return("TELO-BP")
    } else if (next_state == -102) {
        return("BP-TELO")
    } else if (prev_state == -100) {
        return("CENTRO-BP")
    } else if (next_state == -100) {
        return("BP-CENTRO")
    } else {
        return("INTER")
    }
}

classify_edge <- function(state, prev_state, next_state){
    ### state is the state of the segment

    # If the state indicates we don't have an edge. Return None
    if (non_edge(state)) {
        return("NA")
    } else {# we have a real edge here
        edge <- real_edge(state,prev_state,next_state)
        return(edge)
    }
}


classify_edges <- function(chr_seg) {
    state <- as.character(chr_seg$state)
    n <- length(state)
    prev_state <- c("NA", state[-n])
    next_state <- c(state[-1], "NA")

    dplyr::case_when(
        state %in% c("-100", "-101", "-102", "0") ~ "NA",
        prev_state == "-101" & next_state == "-100" ~ "ARM",
        next_state == "-102" & prev_state == "-100" ~ "ARM",
        prev_state == "-101" & next_state == "-102" ~ "ARM",  # whole chrom single state, no centromere
        prev_state == "-101"                        ~ "TELO-BP",
        next_state == "-102"                        ~ "BP-TELO",
        prev_state == "-100"                        ~ "CENTRO-BP",
        next_state == "-100"                        ~ "BP-CENTRO",
        TRUE                                        ~ "INTER"
    )
}

extract_state_info <- function(chr_seg) {
    ### chr_seg: a data frame with a 'state' column for one chromosome from one cell
    state      <- chr_seg$state
    n          <- length(state)
    prev_state <- c("NA", state[-n])
    next_state <- c(state[-1], "NA")

    if (prev_state[1] == "NA" && state[1] == 0) {
        print(chr_seg)
        stop()
    }

    data.frame(state = state, prev_state = prev_state, next_state = next_state,
               stringsAsFactors = FALSE)
}


annotate_segs_with_reads_info <- function(segs, ano_reads){
    ### WHen you collapse ano_reads into segs with reads_to_segs you lose information
    ### This function recovers the median LOH and median background state over the collapsed segments
    required_cols_segs  <- c("chr", "cell_id", "start", "end", "state")
    required_cols_reads <- c("chr", "cell_id", "start", "end", "state", "state_abs")

    if (!all(required_cols_segs %in% colnames(segs)))
        stop("`segs` must contain columns: ", paste(required_cols_segs, collapse = ", "))
    if (!all(required_cols_reads %in% colnames(ano_reads)))
        stop("`ano_reads` must contain columns: ", paste(required_cols_reads, collapse = ", "))

    # Pre-process reads once: type coercions + LOH encoding
    reads_proc <- ano_reads %>%
        mutate(
            start     = as.numeric(start),
            end       = as.numeric(end),
            state_abs = as.numeric(state_abs),
            state     = as.numeric(state)
        )

    if (all(c("LOH", "is_wgd") %in% colnames(ano_reads))) {
        reads_proc <- reads_proc %>%
            mutate(
                LOH    = dplyr::case_when(LOH == "NO" ~ 0, LOH == "LOH" ~ 1),
                is_wgd = as.numeric(is_wgd)
            ) %>%
            select(chr, cell_id, start, end, state_abs, LOH, is_wgd, state)
    }

    segs <- segs %>% mutate(
        .seg_id = seq_len(n()),
        .start  = as.numeric(start),
        .end    = as.numeric(end)
    )

    medians <- segs %>%
        filter(!state %in% c("-101", "-100", "-102")) %>%
        select(.seg_id, cell_id, chr, .seg_start = .start, .seg_end = .end) %>%
        inner_join(reads_proc, by = c("cell_id", "chr")) %>%
        filter(start >= .seg_start, end <= .seg_end) %>%
        mutate(background_CNV = state_abs - state) %>%
        group_by(.seg_id) %>%
        summarise(
            median_loh   = median(LOH,            na.rm = TRUE),
            median_abs   = median(state_abs,      na.rm = TRUE),
            median_wgd   = median(is_wgd,         na.rm = TRUE),
            median_bgCNV = median(background_CNV, na.rm = TRUE),
            .groups = "drop"
        )

    segs %>%
        left_join(medians, by = ".seg_id") %>%
        select(-.seg_id, -.start, -.end)
}


annotate_segs_with_reads_info_noallele <- function(segs, ano_reads){
    ### WHen you collapse ano_reads into segs with reads_to_segs you lose information
    ### This function recovers the median LOH and median background state over the collapsed segments
    ### This function assumes we've added BACK the LOH information that was lost during medic tree construction
    ### Where we used TCN for lineage construction. So we don't have WGD info, but we can add back LOH from SIGNALS

    required_cols_segs  <- c("chr", "cell_id", "start", "end", "state")
    required_cols_reads <- c("chr", "cell_id", "start", "end", "state", "state_abs")

    if (!all(required_cols_segs %in% colnames(segs)))
        stop("`segs` must contain columns: ", paste(required_cols_segs, collapse = ", "))
    if (!all(required_cols_reads %in% colnames(ano_reads)))
        stop("`ano_reads` must contain columns: ", paste(required_cols_reads, collapse = ", "))

    reads_proc <- ano_reads %>%
        mutate(
            start     = as.numeric(start),
            end       = as.numeric(end),
            state_abs = as.numeric(state_abs),
            state     = as.numeric(state)
        )

    if ("LOH" %in% colnames(ano_reads)) {
        reads_proc <- reads_proc %>%
            mutate(LOH = dplyr::case_when(LOH == "NO" ~ 0, LOH == "LOH" ~ 1)) %>%
            select(chr, cell_id, start, end, state_abs, LOH, state)
    }

    segs <- segs %>% mutate(
        .seg_id = seq_len(n()),
        .start  = as.numeric(start),
        .end    = as.numeric(end)
    )

    medians <- segs %>%
        filter(!state %in% c("-101", "-100", "-102")) %>%
        select(.seg_id, cell_id, chr, .seg_start = .start, .seg_end = .end) %>%
        inner_join(reads_proc, by = c("cell_id", "chr")) %>%
        filter(start >= .seg_start, end <= .seg_end) %>%
        mutate(background_CNV = state_abs - state) %>%
        group_by(.seg_id) %>%
        summarise(
            median_abs   = median(state_abs,      na.rm = TRUE),
            median_bgCNV = median(background_CNV, na.rm = TRUE),
            median_loh   = median(LOH,            na.rm = TRUE),
            .groups = "drop"
        )

    segs %>%
        left_join(medians, by = ".seg_id") %>%
        select(-.seg_id, -.start, -.end)
}


#### EXTRACT THE BREAKPOINT LOCATIONS ONLY
get_breakpoint_inputs <- function(chr_seg) {
    ### chr_seg: a data frame with a 'state' column for one chromosome from one cell
    empty_row <- data.frame(
        cell_id   = "NA",
        chr       = "NA",
        start     = "NA",
        end       = "NA",
        state     = "NA",
        seg_width = "NA",
        edge      = "NA",
        stringsAsFactors = FALSE
    )
    n         <- nrow(chr_seg)
    prev_rows <- bind_rows(empty_row, chr_seg[-n, ])
    data.frame(row = chr_seg, prev_row = prev_rows)
}

extract_breakpoints <- function(chr_seg) {
    edge  <- chr_seg[["edge"]]
    start <- as.character(chr_seg$start)
    end_  <- as.character(chr_seg$end)

    result <- dplyr::case_when(
        edge %in% c("ARM", "NA")            ~ "NA",
        edge %in% c("BP-CENTRO", "BP-TELO") ~ start,
        edge %in% c("CENTRO-BP", "TELO-BP") ~ end_,
        edge == "INTER"                     ~ paste0(start, ",", end_),
        TRUE                                ~ NA_character_
    )

    if (any(is.na(result))) stop("Breakpoint location extraction failed")
    result
}

extract_breakpoint <- function(row_edge, prev_row_edge, row_start, row_end) {
    ###
    #This function extracts the breakpoint from a row of a chromosome segment.
    # row_edge: the edge of the current segment
    # prev_row_edge: the edge of the previous segment
    # row_start: the start of the current segment
    # row_end: the end of the current segment
    # Returns the breakpoint location as a numeric value or NA if no breakpoint exists.
    ###

    ### edge: the edge of the current segment
    ### next_edge: the edge of the next segment
    if (row_edge == 'ARM') {
        # Return the centromere breakpoint
        return(extract_ARM_bp(prev_row_edge, row_start, row_end))
    }
    if (row_edge == 'BP-CENTRO') {
        # Return the centromere breakpoint
        return(extract_BP_centro_bp(prev_row_edge, row_start, row_end))
    }
    if (row_edge == 'BP-TELO') {
        # Return the centromere breakpoint
        return(extract_BP_TELO_bp(prev_row_edge, row_start, row_end))
    }
    if (row_edge == 'CENTRO-BP') {
        # Return the centromere breakpoint
        return(extract_CENTRO_BP_bp(prev_row_edge, row_start, row_end))
    }
    if (row_edge == 'INTER') {
        # Return the centromere breakpoint
        return(extract_INTER_bp(prev_row_edge, row_start, row_end))
    }
    if (row_edge == 'TELO-BP') {
        # Return the centromere breakpoint
        return(extract_TELO_BP_bp(prev_row_edge, row_start, row_end))
    }
    if (row_edge == 'NA') {
        # Return the centromere breakpoint
        return(extract_NA_bp(prev_row_edge, row_start, row_end))
    }
    stop("Breakpoint location extraction failed")
    }

extract_ARM_bp <- function(prev_row_edge, row_start, row_end) "NA"

extract_BP_centro_bp <- function(prev_row_edge, row_start, row_end) row_start

extract_BP_TELO_bp <- function(prev_row_edge, row_start, row_end) row_start

extract_CENTRO_BP_bp <- function(prev_row_edge, row_start, row_end) row_end

extract_INTER_bp <- function(prev_row_edge, row_start, row_end) paste0(row_start, ",", row_end)

extract_TELO_BP_bp <- function(prev_row_edge, row_start, row_end) row_end

extract_NA_bp <- function(prev_row_edge, row_start, row_end) "NA"
####### Functions for processing the segment files annotated with breakpoints

process_breakpoints <- function(chr_seg){
    ### Take a chromosome from cell.
    # process the breakpoints such that we only have real breakpoints
    # so that we don't have duplicates
    bps <- chr_seg %>%
        select(cell_id, chr, breakpoint) %>%
        arrange(cell_id, chr, breakpoint) %>%         # ensure proper order
        separate_rows(breakpoint, sep = ",") %>%            # split rows on comma
        mutate(breakpoint = as.numeric(trimws(breakpoint))) %>%
        filter(is.na(lead(breakpoint)) | (lead(breakpoint) - breakpoint != 1))
    return(bps)
}
