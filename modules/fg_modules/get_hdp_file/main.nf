process GET_HDP_FILE {
    label 'copy_output'
    tag "${library_id}-${sample_id}-${exp_con}"
    label 'process_single'

    conda "${moduleDir}/r44.yml"

    input:
        tuple val(sample_id), val(library_id), val(exp_con), path(reads_2), path(alleles), path(metrics), path(hscn), path(tree), path(cnprofiles), path(segs), path(segs_bp_anno), path(segs_bp_anno_A), path(segs_bp_anno_B), path(breakpoint_functions)
    output:
        tuple val(sample_id), val(library_id), val(exp_con),
              path("${library_id}-${sample_id}_hdp_segs.csv.gz"),
              path("${library_id}-${sample_id}_hdp_segs_A.csv.gz"),
              path("${library_id}-${sample_id}_hdp_segs_B.csv.gz"),
              emit: master

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)

    ## ── Total-CN HDP segs (existing behaviour) ───────────────────────────────────
    df <- vroom::vroom("${segs_bp_anno}", delim = ",",
                       col_types = vroom::cols(edge=vroom::col_character(),
                                               breakpoint=vroom::col_double()))
    hdp <- subset(df, !is.na(edge))

    hdp <- hdp %>%
        mutate(edge = as.character(edge)) %>%
        rename(
            state_fg = state,
            width = seg_width,
            background_CNV = median_bgCNV)

    if ("median_loh" %in% colnames(hdp)) {
        hdp <-  hdp %>% rename(median_LOH = median_loh)
    }
    if ("median_wgd" %in% colnames(hdp)) {
        hdp <-  hdp %>% rename(is_wgd = median_wgd)
    }

    hdp <- hdp %>% mutate(Name = stringr::str_split_fixed(cell_id, "-", 2)[,1])
    hdp <- hdp %>% mutate(width = width + 1)

    if (all(c("median_LOH", "is_wgd") %in% colnames(hdp))) {
        hdp <- hdp %>% select(chr, start, end, state_fg, cell_id, Name,
                              median_LOH, background_CNV, is_wgd, width, edge)
    } else if ("median_LOH" %in% colnames(hdp)) {
        hdp <- hdp %>% select(chr, start, end, state_fg, cell_id, Name,
                              median_LOH, background_CNV, width, edge)
    } else if ("is_wgd" %in% colnames(hdp)) {
        hdp <- hdp %>% select(chr, start, end, state_fg, cell_id, Name,
                              is_wgd, background_CNV, width, edge)
    } else {
        hdp <- hdp %>% select(chr, start, end, state_fg, cell_id, Name,
                              background_CNV, width, edge)
    }

    hdp <- hdp %>% mutate(edge = case_when(
        edge == "TELO-BP"   ~ "TELO-PART",
        edge == "BP-TELO"   ~ "TELO-PART",
        edge == "BP-CENTRO" ~ "CENTRO-PART",
        edge == "CENTRO-BP" ~ "CENTRO-PART",
        TRUE ~ edge
    ))

    data.table::fwrite(hdp, file = "${library_id}-${sample_id}_hdp_segs.csv.gz",
                       sep = ",", row.names = FALSE, quote = FALSE, compress = "gzip")

    ## ── Allele-specific HDP segs (A and B) ───────────────────────────────────────
    reads_raw <- vroom::vroom("${reads_2}",
        col_types = vroom::cols(
            chr        = vroom::col_character(),
            state_fg_A = vroom::col_double(),
            state_fg_B = vroom::col_double(),
            state_fg   = vroom::col_double(),
            state      = vroom::col_double(),
            .default   = vroom::col_guess()
        ))

    if (all(c("state_fg_A", "state_fg_B") %in% colnames(reads_raw))) {
        source("${breakpoint_functions}")
        library(data.table)

        for (allele_letter in c("A", "B")) {
            fg_col <- paste0("state_fg_", allele_letter)

            ## Prepare reads: use allele-specific foreground as 'state'
            reads_allele <- reads_raw %>%
                mutate(
                    state_abs = state,
                    state     = .data[[fg_col]]
                )

            ## Mark centromere / telomere bins with sentinel values
            reads_allele <- update_state(reads_allele)

            ## Remove any remaining NA state bins
            reads_allele <- reads_allele %>% filter(!is.na(state))

            ## Segment by allele-specific foreground
            segs_allele <- dlptools::reads_to_segs(reads_allele)

            ## Classify edges
            edges_list <- segs_allele %>%
                group_by(cell_id, chr) %>%
                group_split() %>%
                lapply(classify_edges)
            segs_allele[["edge"]] <- unlist(edges_list)

            ## Extract breakpoints
            bps_list <- segs_allele %>%
                group_by(cell_id, chr) %>%
                group_split() %>%
                lapply(extract_breakpoints)
            segs_allele[["breakpoint"]] <- unlist(bps_list)

            ## Keep only edge-annotated segments (non-centromere/telo)
            hdp_allele <- subset(segs_allele, !is.na(edge) & edge != "NA")

            ## Compute per-segment median background (ancestor allele CN ≈ abs_allele - allele_fg)
            ## Use data.table foverlaps for efficient range join
            reads_bg <- reads_allele %>%
                mutate(abs_allele = state_abs,  ## total CN used as proxy for allele abs
                       bg_allele  = state_abs - .data[[fg_col]]) %>%
                select(cell_id, chr, start, end, bg_allele) %>%
                filter(!is.na(bg_allele))

            reads_dt <- as.data.table(reads_bg)
            segs_dt  <- as.data.table(hdp_allele)[, .(cell_id, chr, seg_start = start, seg_end = end)]
            setkey(reads_dt, cell_id, chr, start, end)
            setkey(segs_dt,  cell_id, chr, seg_start, seg_end)

            overlap <- foverlaps(reads_dt, segs_dt,
                                 by.x = c("cell_id", "chr", "start", "end"),
                                 by.y = c("cell_id", "chr", "seg_start", "seg_end"),
                                 type = "within", nomatch = 0L)

            seg_bg <- overlap[, .(background_CNV = median(bg_allele, na.rm = TRUE)),
                              by = .(cell_id, chr, seg_start, seg_end)]

            hdp_allele_dt <- as.data.table(hdp_allele)
            hdp_allele_dt <- merge(hdp_allele_dt, seg_bg,
                                   by.x = c("cell_id", "chr", "start", "end"),
                                   by.y = c("cell_id", "chr", "seg_start", "seg_end"),
                                   all.x = TRUE)

            hdp_allele_out <- hdp_allele_dt %>%
                mutate(
                    state_fg = state,
                    width    = seg_width + 1,
                    Name     = stringr::str_split_fixed(cell_id, "-", 2)[, 1],
                    edge     = case_when(
                        edge == "TELO-BP"   ~ "TELO-PART",
                        edge == "BP-TELO"   ~ "TELO-PART",
                        edge == "BP-CENTRO" ~ "CENTRO-PART",
                        edge == "CENTRO-BP" ~ "CENTRO-PART",
                        TRUE ~ edge
                    )
                ) %>%
                select(chr, start, end, state_fg, cell_id, Name, background_CNV, width, edge)

            out_name <- paste0("${library_id}-${sample_id}_hdp_segs_", allele_letter, ".csv.gz")
            data.table::fwrite(hdp_allele_out, file = out_name,
                               sep = ",", row.names = FALSE, quote = FALSE, compress = "gzip")
            message("Wrote allele ", allele_letter, " HDP segs to: ", out_name)
        }
    } else {
        message("state_fg_A / state_fg_B not found in reads — writing empty placeholder allele HDP segs")
        empty_hdp <- data.frame(chr=character(), start=numeric(), end=numeric(),
                                state_fg=numeric(), cell_id=character(), Name=character(),
                                background_CNV=numeric(), width=numeric(), edge=character())
        data.table::fwrite(empty_hdp, file="${library_id}-${sample_id}_hdp_segs_A.csv.gz",
                           sep=",", compress="gzip")
        data.table::fwrite(empty_hdp, file="${library_id}-${sample_id}_hdp_segs_B.csv.gz",
                           sep=",", compress="gzip")
    }
    """
}
