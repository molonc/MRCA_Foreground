# NOTE (moved from bin/tests/ during the publication cleanup pass): this suite was
# originally written against an older, loop-based implementation of
# breakpoint_file_functions.R (formerly bin/functions/, now deleted as dead code).
# The live implementation at helper/breakpoint_file_functions.R has since been
# rewritten to a vectorized form with different internal helpers (e.g. get_reads /
# get_reads_noallele no longer exist; classify_edges and extract_breakpoints are
# now vectorized). The source() path below has been updated to point at the live
# file, but the test expectations have NOT been re-verified against it — some
# test_that() blocks below may fail or may need updating to match current
# function signatures/behaviour before this suite can be trusted as a regression
# check. Treat this as a starting point to repair, not a passing suite.
library(testthat)
library(dplyr)
# Source the functions (replace with the correct path if needed)
source(here::here('helper/breakpoint_file_functions.R'))  # Only works if you refactor the logic into an R script without commandArgs code

########
### Generate some real test data
########

# # # Get base case where there is no segments in chr1
# no_segment <- reads %>%
#   filter(cell_id %in% {
#     reads %>%
#       filter(chr == "1") %>%
#       group_by(cell_id) %>%
#       summarise(all_state_fg_0 = all(state_fg == 0 | is.na(state_fg))) %>%
#       filter(all_state_fg_0) %>%
#       slice(1) %>%
#       pull(cell_id)
#   }) %>%
#   filter(chr == "1")

# write.csv(no_segment, file = here::here("bin/tests/chr1_no_segments.csv"), row.names = FALSE)

# cells_with_fg1 <- reads %>%
#   filter(chr == "1") %>%                        # only chromosome 1
#   group_by(cell_id) %>%
#   filter(any(state_fg == 1)) %>%  # keep cells where all state_fg == 1 or NA
#   ungroup() 
# write.csv(cells_with_fg1, file = here::here("bin/tests/chr1_fg1.csv"), row.names = FALSE)


# cells_with_fg2 <- reads %>%
#   filter(chr == "1") %>%                        # only chromosome 1
#   group_by(cell_id) %>%
#   filter(any(state_fg == 2)) %>%  # keep cells where all state_fg == 1 or NA
#   ungroup() 
# write.csv(cells_with_fg2, file = here::here("bin/tests/chr1_fg2.csv"), row.names = FALSE)

# cells_with_fgmin1 <- reads %>%
#   filter(chr == "1") %>%                        # only chromosome 1
#   group_by(cell_id) %>%
#   filter(any(state_fg == -1)) %>%  # keep cells where all state_fg == 1 or NA
#   ungroup()  
# write.csv(cells_with_fgmin1, file = here::here("bin/tests/chr1_fgmin1.csv"), row.names = FALSE)


############
#### TESTS FOR  mark_centromeres_and_telomeres
############

# This contains genomic bins of every potential is_p_telo, is _centro, and is _q_telo type that we can have for normal chromosomes
df_mark_centro_telo_1 <- rbind(
  data.frame(chr = "1", start = 50000,     end =  60000), # in p_telo. TRUE
  data.frame(chr = "1", start = 550000,     end = 750000), # on boundary of p-telo, TRUE
  data.frame(chr = "1", start = 740000,     end = 760000), # overlapping boundary of p-telo, TRUE
  data.frame(chr = "1", start = 750000,     end = 850000), # on boundary of p-telo, TRUE

  data.frame(chr = "1", start = 950000,     end = 860000), # inside. FALSE for everything
  data.frame(chr = "1", start =121260000,   end =121270001), # bounary of centromere# should , should be FALSE for is_centro
  data.frame(chr = "1", start =121270000,   end =121280000), # overlapping boundary of centromere should be TRUE for is_centro
  data.frame(chr = "1", start =121280001,   end =123025434),  # inside centromere should be TRUE for is_centro
  data.frame(chr = "1", start =140000000,   end =150000000), # ad end of centromere end # should be TRUE should be TRUE for is_centro
    data.frame(chr = "1", start =140000000,   end =160000000), # overlapping centromere endshould be TRUE should be TRUE for is_centro
  data.frame(chr = "1", start =150000001,   end =160000000), # at boundary of centroemre. should be false 
  data.frame(chr = "1", start = 160000000 ,   end =170000000), # inner
  data.frame(chr = "1", start = 249220000 ,   end =249220001), #  boundary of q-tel, should be FALSE for q-telo
  data.frame(chr = "1", start = 249220000 ,   end =249520000), #  overapping qtel, should be true
  data.frame(chr = "1", start =249220001,   end =249230001) # in q-tel, should be TRUE
)
test_that("mark_centromeres_and_telomeres annotates telomeres and centromeres correctly", {
  results <- mark_centromeres_and_telomeres(df_mark_centro_telo_1 , centro_telo_locs)
  expect_equal(results$is_p_telo, c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(results$is_q_telo, c(FALSE,FALSE,FALSE,FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,TRUE,TRUE))
    expect_equal(results$is_centro, c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,FALSE, FALSE, FALSE, FALSE, FALSE))
})

# This contains genomic bins of every potential is_p_telo, is _centro, and is _q_telo type that we can have for chromosomes that are telo-centromeric
df_mark_centro_telo_13 <- rbind(
  data.frame(chr = "13", start = 1,     end =  60000), # boundary. but  in p_telo. Yes, but FALSE because also in centromere , true for is_centro
  data.frame(chr = "13", start = 50000,     end =  60000), # in p_telo. Yes, but FALSE because also in centromere. true for is_centro
  data.frame(chr = "13", start = 19300000,     end = 19360000), # on boundary of centromere end. Should be false for p_tel and true for is_centro
    data.frame(chr = "13", start = 19300000,     end = 19400000), # overlapping boundary of centromere end. Should be false for p_tel and true for is_centro
    data.frame(chr = "13", start = 19360000 ,     end = 19400000 ), #  on boundary of centromere end. and inner Should be false for is_centro
    data.frame(chr = "13", start = 19400000 ,     end = 19500000 ), # inner. Should be false
    data.frame(chr = "13", start = 115100001,     end = 115110001 ), # boundary of q_tel. should be false for q-telo
    data.frame(chr = "13", start = 115000001,     end = 115210001 ), # overlapping qdtel. should be true fore q_telo
    data.frame(chr = "13", start = 115000001,     end = 115210001 ), # in q_telo should be true fore q_telo
    data.frame(chr = "13", start =115210001,     end = 115310001) # in q_telo should be true fore q_telo
)
test_that("mark_centromeres_and_telomeres annotates telomeres and centromeres correctly", {
  results <- mark_centromeres_and_telomeres(df_mark_centro_telo_13 , centro_telo_locs)
  results$is_centro
  expect_equal(results$is_p_telo, c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(results$is_q_telo, c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(results$is_centro, c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
})

# This contains genomic bins of every potential is_p_telo, is _centro, and is _q_telo type that we can have for chr17 as it donesn't have a marked p-telomere
df_mark_centro_telo_17 <- rbind(
    data.frame(chr = "17", start = 1,     end =  50000) # boundary. but  in p_telo. Yes, but FALSE because also in centromere , true for is_centro
)
test_that("mark_centromeres_and_telomeres annotates telomeres and centromeres correctly", {
  results <- mark_centromeres_and_telomeres(df_mark_centro_telo_17, centro_telo_locs)
  print(results)
  expect_equal(results$is_p_telo, c(TRUE, FALSE))
  expect_equal(results$is_q_telo, c(FALSE, FALSE))
  expect_equal(results$is_centro, c(FALSE, FALSE))
})


############
#### TESTS FOR  classify_edge
############
source(here::here('helper/breakpoint_file_functions.R'))  # Only works if you refactor the logic into an R script without commandArgs code

### Test non_edge() ###
test_that("non_edge returns correct TRUE/FALSE for states", {
  expect_true(non_edge(-100))
  expect_true(non_edge(-101))
  expect_true(non_edge(-102))
  expect_true(non_edge(0))
  expect_false(non_edge(1))
  expect_false(non_edge(5))
})

### Test real_edge() ###
test_that("real_edge classifies correct edge type", {
  expect_equal(real_edge(1, -101, 2), "TELO-BP")
  expect_equal(real_edge(1, 2, -102), "BP-TELO")
  expect_equal(real_edge(1, -100, 2), "CENTRO-BP")
  expect_equal(real_edge(1, 2, -100), "BP-CENTRO")
  expect_equal(real_edge(1, -101, -100), "ARM")
  expect_equal(real_edge(1, -100, -102), "ARM")
  expect_equal(real_edge(1, 3, 4), "INTER")
})

## Test classify_edges(). Just make sure it runs through ###
test_that("classify_edges failed to completely run though without any errors", {
    # first run on a bunch of examples to make sure I have all cases handeled
    reads <- vroom::vroom(here::here("bin/tests/chr1_no_segments.csv"))
    reads <- subset(reads, chr != "X")
    reads <- subset(reads, chr != "Y")
    reads$state <- reads$state_fg
    ano_reads <-  update_state(reads)
    segs <- dlptools::reads_to_segs(ano_reads)
    lo_edges <- classify_edges(segs %>% group_by(cell_id, chr))
    expect_equal(TRUE,TRUE)  # Test to make sure we get through. Everything is captured?
})

### Test classify_edge() ###
test_that("classify_edge handles all possible state combinations", {
  for (i in seq_len(nrow(distinct_states))) {
    row <- distinct_states[i, ]
    state <- row$state
    prev_state <- row$prev_state
    next_state <- row$next_state
    label <- as.character(row$label)
    
    edge <- classify_edge(state, prev_state, next_state)
    expect_equal(
      edge,
      label,
      info = paste0(
        "❌ Mismatch at row ", i, "\n",
        "  state       = ", state, "\n",
        "  prev_state  = ", prev_state, "\n",
        "  next_state  = ", next_state, "\n",
        "  expected    = ", label, "\n",
        "  actual edge = ", edge
      )
    )
}})

############
#### TESTS FOR  extract_breakpoint()
############
test_that("classify_edge handles all possible state combinations", {
  possible_bp <- vroom::vroom(here::here("bin/tests/possible_breakpoint_extractions.tsv"), na = character())
  possible_bp$returns <- as.character(possible_bp$returns)
  possible_bp$start <- as.character(possible_bp$start)
  possible_bp$end <- as.character(possible_bp$end)

  for (i in seq_len(nrow(possible_bp))) {
    row <- possible_bp[i, ]
    print(row)
    row_edge <- row$row
    prev_row_edge <- row$prev_row
    row_start <- row$start
    row_end <- row$end
    return <- row$returns

    result <- extract_breakpoint(row_edge, prev_row_edge, row_start, row_end)
    
    print(result)
    expect_equal(
      result,
      return,
      info = paste0(
        "❌ Mismatch at row ", i, "\n",
        "  row_edge  = ", row_edge, "\n",
        "  prevrow  = ", prev_row_edge, "\n",
        "  expected_ resu;t  = ", return, "\n",
        "  expected_ result type  = ", typeof(return), "\n",
        "  result       = ", result, "\n",
        "  result type = ", typeof(result))
    )
    }
})

############
#### TESTS FOR  extract_breakpoints()
############

test_that("classify_edges runs through without any errors when theres no segments",{
    # first run on a bunch of examples to make sure I have all cases handeled
    reads <- vroom::vroom(here::here("bin/tests/chr1_no_segments.csv"))
    reads <- subset(reads, chr != "X")
    reads <- subset(reads, chr != "Y")
    reads$state <- reads$state_fg
    ano_reads <-  update_state(reads)
    segs <- dlptools::reads_to_segs(ano_reads)
    lo_edges <- classify_edges(segs %>% group_by(cell_id, chr))
    segs[['edge']] <- lo_edges
    lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
    segs[["breakpoint"]] <- lo_breakpoints
    expect_equal(TRUE,TRUE)  # Test to make sure we get through. Everything is captured?
})

test_that("classify_edges runs through without any errors",{
    # first run on a bunch of examples to make sure I have all cases handeled
    reads <- vroom::vroom(here::here("bin/tests/chr1_fg1.csv"))
    reads <- subset(reads, chr != "X")
    reads <- subset(reads, chr != "Y")
    reads$state <- reads$state_fg
    ano_reads <-  update_state(reads)
    segs <- dlptools::reads_to_segs(ano_reads)
    lo_edges <- classify_edges(segs %>% group_by(cell_id, chr))
    segs[['edge']] <- lo_edges
    t <- segs %>% filter(edge == 'BP-TELO')
    lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
    segs[["breakpoint"]] <- lo_breakpoints
    expect_equal(TRUE,TRUE)  # Test to make sure we get through. Everything is captured?
})

test_that("classify_edges grabs correct bp for bp_telo edges",{
    # first run on a bunch of examples to make sure I have all cases handeled
    segs <- vroom::vroom(here::here("bin/tests/bp_examples/bp_telo.tsv"), col_types = vroom::cols(.default="c"),   na = character())   # don't treat anything as NA)
    lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
    segs[["breakpoint"]] <- lo_breakpoints
    print(segs)
    expect_equal(segs$breakpoint, segs$returns)  # Test to make sure we get through. Everything is captured?
})

test_that("classify_edges grabs correct bp for telo_bp edges",{
    # first run on a bunch of examples to make sure I have all cases handeled
    segs <- vroom::vroom(here::here("bin/tests/bp_examples/telo_bp.tsv"), col_types = vroom::cols(.default="c"),   na = character())   # don't treat anything as NA)
    lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
    segs[["breakpoint"]] <- lo_breakpoints
    print(segs)
    expect_equal(segs$breakpoint, segs$returns)  # Test to make sure we get through. Everything is captured?
})

test_that("classify_edges grabs correct bp for bp_centro edges",{
    # first run on a bunch of examples to make sure I have all cases handeled
    segs <- vroom::vroom(here::here("bin/tests/bp_examples/bp_centro.tsv"), col_types = vroom::cols(.default="c"),   na = character())   # don't treat anything as NA)
    lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
    segs[["breakpoint"]] <- lo_breakpoints
    print(segs)
    expect_equal(segs$breakpoint, segs$returns)  # Test to make sure we get through. Everything is captured?
})

test_that("classify_edges grabs correct bp for isioated inter segments",{
    # first run on a bunch of examples to make sure I have all cases handeled
    segs <- vroom::vroom(here::here("bin/tests/bp_examples/inter1.tsv"), col_types = vroom::cols(.default="c"),   na = character())   # don't treat anything as NA)
    lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
    segs[["breakpoint"]] <- lo_breakpoints
    print(segs)
    expect_equal(segs$breakpoint, segs$returns)  # Test to make sure we get through. Everything is captured?
})


test_that("classify_edges grabs correct bp for paired inter segments",{
    # first run on a bunch of examples to make sure I have all cases handeled
    segs <- vroom::vroom(here::here("bin/tests/bp_examples/inter2.tsv"), col_types = vroom::cols(.default="c"),   na = character())   # don't treat anything as NA)
    lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
    segs[["breakpoint"]] <- lo_breakpoints
    print(segs)
    expect_equal(segs$breakpoint, segs$returns)  # Test to make sure we get through. Everything is captured?
})

test_that("classify_edges grabs correct bp for inter with telo segments",{
    # first run on a bunch of examples to make sure I have all cases handeled
    segs <- vroom::vroom(here::here("bin/tests/bp_examples/intercombo.tsv"), col_types = vroom::cols(.default="c"),   na = character())   # don't treat anything as NA)
    lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
    segs[["breakpoint"]] <- lo_breakpoints
    print(segs)
    expect_equal(segs$breakpoint, segs$returns)  # Test to make sure we get through. Everything is captured?
})



#  TODO. GEt the real expected values from real data
# reads <- vroom::vroom(here::here("/Users/ahamazaki/research/packages/medicc2_foreground-master/results/A146279A/A146279A_readds_final.csv.gz"))
# reads <- subset(reads, chr != "X")
# reads <- subset(reads, chr != "Y")
# reads$state <- reads$state_fg
# print(subset(reads, chr=="17"),Inf)
# print(subset(reads, chr=="13"),Inf)
# ano_reads <- update_state(reads)
# segs <- dlptools::reads_to_segs(ano_reads)

# lo_edges <- classify_edges(segs %>% group_by(cell_id, chr))
# lo_edges
# segs[['edge']] <-lo_edges

# segs

# distinct_breakpoint_inputs  <- get_breakpoint_inputs(segs %>% group_by(cell_id, chr))
# test distinct_breakpoint_inputs[[1]]
# test <- bind_rows(distinct_breakpoint_inputs) %>%
#     select(row.edge, prev_row.edge) %>%
#     distinct() %>%
#     arrange(row.edge)
# test
# head(test)
# head(distinct_breakpoint_inputs)
# table(distinct_breakpoint_inputs)

# print("Extracting Breakpoints")
# lo_breakpoints <- extract_breakpoints(segs %>% group_by(cell_id, chr))
# lo_breakpoints

# bind_rows(lo_breakpoints)

# chr_groups <- segs %>%
#   group_by(chr, cell_id) %>%
#   group_split()
#   # Apply your function to each group
# results <-bind_rows(lapply(chr_groups, extract_state_info))
# results
# distinct_states <- results %>%
#     select(state, prev_state) %>%
#       mutate(
#             prev_state = ifelse(!(prev_state %in% valid_states) & !is.na(prev_state), 5, prev_state),
#             state = ifelse(!(state %in% valid_states) & !is.na(state), 3, state)
#       ) %>%
#     distinct()

# distinct_states
# head(results)
# # Optionally: bind all results back together
# valid_states <- c(-101, -102, -100, 0, NA)
# results
