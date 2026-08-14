###
# This script takes in segs_hdp files and outputs a _bp file  where we can easily determine the exact location of each breakpoint
###
main <- function() {
  library(dplyr)
  library(dlptools)
  library(vroom)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Please provide both 'reads' and 'library_id' arguments.")
  }

  reads_file <- args[1]
  library_id <- args[2]
  tcn_bool<- args[3]
  breakpoint_functions <- args[4]

  source(breakpoint_functions)

  reads <- vroom::vroom(reads_file,
     col_types = vroom::cols(
      chr = col_character(),
      start = col_double(),
      end = col_double(),
      width = col_double(),
      cell_id = col_character(),
      state = col_character(),
      state_fg = col_character(),
      state_fg_A = col_double(),
      state_fg_B = col_double(),
      state_AS = col_character(),
      state_phase = col_character(),
      LOH = col_character(),
      is_wgd = col_logical(),
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

  # reads <- subset(reads, !(chr %in% c("X", "Y")))
  reads$state_abs <- reads$state
  reads$state <- reads$state_fg

  # ── Total foreground annotation ────────────────────────────────────────────
  ano_reads <- update_state(reads)
  print(head(ano_reads))
  segs <- dlptools::reads_to_segs(ano_reads)

  if (tcn_bool == 'false' || tcn_bool == "FALSE" || tcn_bool == FALSE){
    segs <- annotate_segs_with_reads_info(segs, ano_reads)
  }
  else if (tcn_bool == 'true' || tcn_bool == "TRUE" || tcn_bool == TRUE){
    print("Annotating segs with reads info no allele info used")
    segs <- annotate_segs_with_reads_info_noallele(segs, ano_reads)
  }
  print(head(segs))
  print("Extracting Edges")
  lo_edges <- segs %>%
    group_by(cell_id, chr) %>%
    group_split() %>%
    lapply(classify_edges)
  segs[['edge']] <- unlist(lo_edges)

  print("Extracting Breakpoints")
  lo_breakpoints <- segs %>%
    group_by(cell_id, chr) %>%
    group_split() %>%
    lapply(extract_breakpoints)
  segs[["breakpoint"]] <- unlist(lo_breakpoints)

  output_file <- paste0(library_id,"_segs_bp_anno.csv.gz")
  print(paste("Saving breakpoint file to:", output_file))
  write.csv(segs, output_file, row.names = FALSE, quote = TRUE, na = "NA")
  print("Done")

  # ── Allele-specific annotation (A and B) ───────────────────────────────────
  has_allele_fg <- all(c("state_fg_A", "state_fg_B") %in% colnames(reads))

  for (allele_letter in c("A", "B")) {
    output_allele <- paste0(library_id, "_segs_bp_anno_", allele_letter, ".csv.gz")

    if (!has_allele_fg) {
      message("state_fg_A/state_fg_B not found — writing empty placeholder for allele ", allele_letter)
      write.csv(data.frame(), output_allele, row.names = FALSE)
      next
    }

    fg_col <- paste0("state_fg_", allele_letter)
    reads_allele <- reads
    # Character cast matches the total-fg path (state_fg is col_character()), so
    # update_state sentinels and classify_edges comparisons behave identically.
    reads_allele$state <- as.character(reads_allele[[fg_col]])

    ano_reads_allele <- update_state(reads_allele)
    segs_allele <- dlptools::reads_to_segs(ano_reads_allele)
    segs_allele <- annotate_segs_with_reads_info(segs_allele, ano_reads_allele)

    lo_edges_allele <- segs_allele %>%
      group_by(cell_id, chr) %>%
      group_split() %>%
      lapply(classify_edges)
    segs_allele[['edge']] <- unlist(lo_edges_allele)

    lo_bps_allele <- segs_allele %>%
      group_by(cell_id, chr) %>%
      group_split() %>%
      lapply(extract_breakpoints)
    segs_allele[["breakpoint"]] <- unlist(lo_bps_allele)

    print(paste("Saving allele", allele_letter, "annotated segs to:", output_allele))
    write.csv(segs_allele, output_allele, row.names = FALSE, quote = TRUE, na = "NA")
  }
}

# This ensures the script only runs when executed directly
if (sys.nframe() == 0) {
  main()
}


### Debug stuff
# df <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance/A108790B/A108790B_reads_final.csv.gz")
# print(head(df),Inf)
# print(df,Inf)

# unique(df$cell_id)

# df2 <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance/A108790B/A108790B_seg.csv.gz")

# df2$state_fg
# df <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/persistance/A108790A/A108790A_reads_final.csv.gz")
# print(head(df),Inf)
# reads <- vroom::vroom("/Users/ahamazaki/research/packages/medicc2_foreground-master/work/b1/7e4808588a56b438d577455f5bdb01/A108746A_reads_final.csv.gz")
# reads <- subset(reads, chr != "X")
# reads <- subset(reads, chr != "Y")
# reads$state <- reads$state_fg
# ano_reads <- update_state(reads)

# subset(reads, chr=="X")
# reads
# head(df)



# df2 <- update_state(df)
# df2$state_fg
