library(stringr)
library(vroom)
library(tibble)
library(dplyr)
library(ggplot2)

data <- read.csv("../orig_geninsta-data.csv")
data$SCID <- str_extract(data$Results, "SC-[0-9]+")
data$Cell.Line <- sub(" $", "", data$Cell.Line)

# Name
# Cell.Line
# passage
# Drug.Treatment
# Fixation
# Media
# Conc.of.Drug
# Treatment..days.

metrics_files <- Sys.glob("../data/*/*/annotation/*metrics.csv.gz")
names(metrics_files) <- str_extract(metrics_files, "SC-[0-9]+")
metrics_files <- subset(metrics_files, names(metrics_files) %in% data$SCID)

metrics_raw <- lapply(metrics_files, vroom)
metrics_raw <- mapply(add_column, .data = metrics_raw, file = metrics_files, SIMPLIFY = FALSE)
metrics <- bind_rows(metrics_raw)
metrics$SCID <- str_extract(metrics$file, "SC-[0-9]+")
datamin <- select(data, Name:Treatment..days., DLP.Seq.Lib, SCID)

# metrics2 <- merge(metrics, datamin, all.x = TRUE)

# ISSUE, merge changes row count of metrics, due to multiples in datamin

# > sort(table(datamin$SCID))

# SC-2610 SC-2619 SC-2963 SC-2964 SC-3439 SC-4334 SC-4337 SC-4349 SC-4370 SC-4392 
#       1       1       1       1       1       1       1       1       1       1 
# SC-4423 SC-4526 SC-7015 SC-7285 SC-7318 SC-7319 SC-7321 SC-7322 SC-7326 SC-7436 
#       1       1       1       1       1       1       1       1       1       1 
# SC-7438 SC-7464 SC-7468 SC-7471 SC-7550 SC-7758 SC-7759 SC-7787 SC-7788 SC-7789 
#       1       1       1       1       1       1       1       1       1       1 
# SC-7794 SC-7796 SC-7797 SC-7798 SC-7799 SC-7800 SC-7803 SC-7805 SC-7809 SC-7811 
#       1       1       1       1       1       1       1       1       1       1 
# SC-7812 SC-7814 SC-7815 SC-7816 SC-7837 SC-7838 SC-7840 SC-7841 SC-7842 SC-7848 
#       1       1       1       1       1       1       1       1       1       1 
# SC-7877 SC-7885 SC-7886 SC-7889 SC-7909 SC-7921 SC-7975 SC-7999 SC-8063 SC-8069 
#       1       1       1       1       1       1       1       1       1       1 
# SC-8086 SC-8122 SC-8244 SC-8245 SC-8259 SC-8260 SC-8261 SC-8271 SC-8274 SC-8326 
#       1       1       1       1       1       1       1       1       1       1 
# SC-8331 SC-8333 SC-8347 SC-8350 SC-8353 SC-8355 SC-8357 SC-8497 SC-8499 SC-8501 
#       1       1       1       1       1       1       1       1       1       1 
# SC-8505 SC-8507 SC-8509 SC-8525 SC-7917 SC-7971 
#       1       1       1       1       2       2 

# > subset(datamin, SCID == "SC-7917")
#                                                   Name
# 24 hTERT p53-/- BRCA2-/- 112.109-cisplatin-0.306 nM-7d
# 34      hTERT p53-/- BRCA1-/- 83.86-cisplatin-43 nM-7d
#                        Cell.Line passage Drug.Treatment Fixation    Media
# 24 hTERT p53-/- BRCA2-/- 112.109      56      cisplatin    Fresh DMEM+FBS
# 34   hTERT p53-/- BRCA1-/- 83.86      55      cisplatin    Fresh DMEM+FBS
#    Conc.of.Drug Treatment..days.     DLP.Seq.Lib    SCID
# 24     0.306 nM                7 AT13848 A98187B SC-7917 (EXP B)
# 34        43 nM                7 AT13847 A98187A SC-7917 (EXP A)

# > subset(datamin, SCID == "SC-7971")
#                                          Name                     Cell.Line
# 16 hTERT p53-/- BRCA2-/- 112.109-Untreated--d hTERT p53-/- BRCA2-/- 112.109
# 29   hTERT p53-/- BRCA1-/- 83.86-Untreated--d   hTERT p53-/- BRCA1-/- 83.86
#    passage Drug.Treatment Fixation              Media Conc.of.Drug
# 16      37      Untreated    Fresh MEBM defined media             
# 29      37      Untreated    Fresh MEBM defined media             
#    Treatment..days.      DLP.Seq.Lib    SCID
# 16                           AT10100 SC-7971 (EXP B)
# 29                  AT10099 A118413A SC-7971 (EXP A)


# FIX for libraries with two libraries... that are split by experimental_condition
idx <- which(metrics$SCID == "SC-7971" & grepl("^A", metrics$experimental_condition))
metrics$SCID[idx] <- paste0(metrics$SCID[idx], "-A")
idx <- which(metrics$SCID == "SC-7971" & grepl("^B", metrics$experimental_condition))
metrics$SCID[idx] <- paste0(metrics$SCID[idx], "-B")

idx <- which(metrics$SCID == "SC-7917" & grepl("^A", metrics$experimental_condition))
metrics$SCID[idx] <- paste0(metrics$SCID[idx], "-A")
idx <- which(metrics$SCID == "SC-7917" & grepl("^B", metrics$experimental_condition))
metrics$SCID[idx] <- paste0(metrics$SCID[idx], "-B")

idx <- which(datamin$SCID == "SC-7917" & datamin$DLP.Seq.Lib == "AT13847 A98187A")
datamin$SCID[idx] <- paste0(datamin$SCID[idx], "-A")
idx <- which(datamin$SCID == "SC-7917" & datamin$DLP.Seq.Lib == "AT13848 A98187B")
datamin$SCID[idx] <- paste0(datamin$SCID[idx], "-B")

idx <- which(datamin$SCID == "SC-7971" & datamin$DLP.Seq.Lib == "AT10099 A118413A")
datamin$SCID[idx] <- paste0(datamin$SCID[idx], "-A")
idx <- which(datamin$SCID == "SC-7971" & datamin$DLP.Seq.Lib == "AT10100")
datamin$SCID[idx] <- paste0(datamin$SCID[idx], "-B")

dat <- merge(metrics, datamin, all.x = TRUE)
dat$jira <- stringr::str_extract(dat$file, "SC-[0-9]+")

dat$note <- "negative"
dat$note[grep("hTERT", dat$experimental_condition)] <- "hTERT"
dat$note[dat$experimental_condition %in% c("A", "B", "p20b", "p30b")] <- "cell"

htert <- subset(dat, note == "hTERT")
htert_tmp <- htert %>% group_by(SCID) %>% summarise(avgq = mean(quality)) %>% arrange(desc(avgq))
htert$SCID <- factor(htert$SCID, levels = htert_tmp$SCID)
ggplot(htert, aes(SCID, quality, fill = Drug.Treatment)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cells <- subset(dat, note == "cell")
cells_tmp <- cells %>% group_by(SCID) %>% summarise(avgq = mean(quality)) %>% arrange(desc(avgq))
cells$SCID <- factor(cells$SCID, levels = cells_tmp$SCID)
ggplot(cells, aes(SCID, quality, fill = Drug.Treatment)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# > table(cells$Cell.Line, cells$Drug.Treatment)
                                 
#                                   AURKBi ZM447439 BMH-21 cisplatin CX5461 DMSO etoposide  PDS SB-124 Untreated voreloxin
#   184-hTERT-L9 (WT)                             0   1158      3141   8353    0         0 1113   2247      8331         0
#   HCT-116 BRCA2-/- B18                          0      0         0   1290    0         0    0      0      1194         0
#   HEK293                                     1114      0         0   1048 1150       936    0      0      1073      2216
#   hTERT p53-/- 95.22                            0   1171      2777   8020    0         0 1100   2178      4654         0
#   hTERT p53-/- BRCA1-/- 83.86                   0      0      2131   4259    0         0 1046   2104      4033         0
#   hTERT p53-/- BRCA1-/+ 85.14                   0      0         0      0    0         0    0      0       994         0
#   hTERT p53-/- BRCA2-/- 112.109                 0      0      2181   6842    0         0  915   1117      3917         0
#   hTERT p53-/- BRCA2-/+ 116.66                  0      0         0      0    0         0    0      0       496         0
#   RPE-1 cas9 p53 KO                             0      0         0      0    0         0    0      0      1166      1085
#   RPE-1 cas9 p53 KO ARID1A -/+ #5               0      0         0      0    0         0    0      0      1216         0
#   RPE-1 cas9 p53 WT                             0      0         0      0    0      1268    0      0      1148       925
#   RPE-1 cas9 p53 WT ARID1A -/- #7               0      0         0      0    0         0    0      0      1216         0

hq <- subset(cells, quality > 0.75)
ggplot(hq, aes(SCID, quality, fill = Drug.Treatment)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

table(hq$Cell.Line, hq$Drug.Treatment)

