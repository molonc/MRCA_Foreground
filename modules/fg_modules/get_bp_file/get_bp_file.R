###
# This script takes in segs_hdp files and outputs a _bp file  where we can easily determine the exact location of each breakpoint
###
main <- function() {
  library(dplyr)
  library(dlptools)
  library(vroom)
  
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Please provide both 'reads' and 'sample_id' arguments.")
  }

  annosegs_file <- args[1]
  sample_id <- args[2]
  breakpoint_functions <- args[3]
  allele_suffix <- if (length(args) >= 4 && nchar(args[4]) > 0) paste0("_", args[4]) else ""
  source(breakpoint_functions)

  annotated_segs <- vroom::vroom(annosegs_file,
     col_types = vroom::cols(
      chr = col_character(),
      start = col_double(),
      end = col_double(),
      width = col_double(),
      cell_id = col_character(),
      state = col_character(),
      edge = col_character(),
      breakpoint = col_character()

    )
  )
  print("Remove Fake brekapoints")
  lo_breakpoints <- annotated_segs %>%
    filter(!is.na(edge)) %>%
    group_by(cell_id, chr) %>%
    group_split() %>%
    lapply(process_breakpoints)
    
  breakpoints_df <- do.call(rbind, lo_breakpoints)

  output_file <- paste0(sample_id, "_breakpoint_file", allele_suffix, ".csv.gz")
  print(paste("Saving breakpoint file to:", output_file))
  write.csv(breakpoints_df, output_file, row.names = FALSE, quote = TRUE, na = "NA")
  print("Done")
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

### Got get breakpoint file. 
# process per chromosome per cell
# Order by start
# if previous row's end is the current row_s start -1, dont add a breakpoint
# then extract breakpoints accordingly.
# If any breakpoints are right next to eachother. Keep the 'start' column
