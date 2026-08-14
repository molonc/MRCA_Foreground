## PDX_Helper.R



library(httr)
library(googlesheets4)
library(stringr)
library(vroom)
library(dplyr)

# query_colossus_analysis_by_library_id <- function(str) {
# 	query <- paste0("colossus.molonc.ca/api/analysis_information/?library__pool_id=", str)
# 	g <- GET(query, authenticate("dalai", "greenhelmetdonut"))
# 	json <- content(g)$results
# 	return(json)
# }

# lib <- "A98284A"

# results <- query_colossus_analysis_by_library_id(lib)

# if (length(results) != 1) {
# 	warning("More than 1 DLP+ results, odd but not unprecented")
# }

# for (result in results) {
# 	# Iterate through each analysis
# 	sc_id <- result$analysis_run$blob_path
# 	print(sc_id)
# }


# https://molonc.atlassian.net/browse/MAGIC-92
# GOOGLE SHEETS LINK
ss <- "https://docs.google.com/spreadsheets/d/1GKsLG1AF0gC7XPGaY1BrnYRi9Pr9wyAtRq4I12b_gFw/edit#gid=1550151387"


# WE REMOVED 9 duplicate mice
data <- read_sheet(ss = ss, sheet = "DLP_Metadata_PDX")
# libs <- data$jira_ticket_for_DLP

metrics <- vroom("metrics_pdx.csv")

# table(is.na(metrics$sample_id))
# FALSE  TRUE 
# 18472 95357

# SCID CHECK
datascid <- unlist(strsplit((subset(data, SCID_Broken != "NA")$SCID_Broken), split = " "))
metscid <- unique(metrics$SCID)

setdiff(datascid, metscid)
setdiff(metscid, datascid)
# SCID CHECK END

# THERE WILL BE PAIN for pooled libraries
metrics$library_id <- str_extract(metrics$cell_id, "A[0-9]{5,6}[A-Z]")
metrics$sample_id <- str_extract(metrics$cell_id, "^[^-]+")

lib_check <- metrics %>% group_by(library_id) %>% summarise(scid = unique(SCID))
print("WE HAVE A 1 to 1 mapping between LIBRARY to SCID")

# ROW 15 TEST
subset(lib_check, scid %in% c("SC-2878", "SC-3139"))
# PASSED

# WHY 1 MICE 2 LIBS??
test <- subset(metrics, SCID %in%  c("SC-2878", "SC-3139"))

# MOUSE_ID CHECK
datamouse <- data$mouse_id
metmouse <- metrics$sample_id

# THIS LOOKS GOOD :D
setdiff(datamouse, metmouse)

# THIS DOES NOT LOOK GOOD
setdiff(metmouse, datamouse)

# > setdiff(metmouse, datamouse)
# [1] "SA928"            "SA1015"           "SA039"            "SA1227"          
# [5] "SA1150"           "SA535X4XB05649"   "SA535X4XB09109"   "SA1142X2XB09526" 
# [9] "SA1050CX1XB01443"


# SA535X4XB05649 is not an expected library (from SC-3173)
# This was pulled in due to row 38 stating that A98232A was in the study
# However, we found that row 38 is actually for SA1035X8XB03425
# SA1035X8XB03425 is only in TWO chips: A98282A, A98283A
# The theory is that A98282A was mistyped as A98232A

# TODO, A98282A is https://colossus.molonc.ca/dlp/library/557 is SC-3045, GET IT!
# SA535X4XB09109 

# SA535X4XB05649 is an untreated sample from the metastasis project, REMOVE SC-3173

# SA535X4XB09109 is from row 47 library A98244A (A98168A, A98244A)
# HOWEVER, row 47 actually has mouse SA535X9XB03617
# SA535X9XB03617 is confirmed to be in both (A98168A, A98244A)
# We get SA535X4XB09109 CUZ A98244A has BOTH SA535X9XB03617 + SA535X4XB09109

# SA1142X2XB09526 is in library_? A108833B as a passenger
# lib <- subset(metrics, library_id == "A108833B")

# TODO:
# ADD SC-3045
# REMOVE SC-3173
# REMOVE SA535X4XB09109 from metrics with a subset
# REMOVE SA1142X2XB09526 from metrics with a subset
# REMOVE SA1050CX1XB01443 from metrics with a subset

metrics2 <- merge(data, metrics, by.x = "mouse_id", by.y = "sample_id")

name_check <- metrics2 %>% group_by(mouse_id) %>% summarise(scid = unique(Name))
