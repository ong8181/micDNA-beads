####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No.3 STD sequences assignment
#### 2018.6.22 Ushio
#### 2019.2.20 revised, Ushio
#### R 3.5.2
####

# Load library
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(tidyverse); packageVersion("tidyverse") #v1.2.1

# Create output directory
dir.create("03_TaxaSTDcombineOut")

# Load Claident taxa assignment
claident.tax <- read.delim("01_ProSeqDADA2_Out/ProkASV_semiall_classigntax")

# Load Blastn taxa assignment for STD sequences
# Select sequences with < 3 mismatches and length > 250
blastn.std <- read.table("02_ident_STD_BLASTnOut/ProkSTD_out.txt") %>%
  filter(.$V5 < 3 & .$V6 < 1 & .$V4 > 250)

# Check claident taxa assignments of the potential std sequences
potential.std.id <- match(blastn.std$V1, claident.tax$query)

# Replace STD taxa names with claident taxa assigment
claident.tax$species <- as.character(claident.tax$species)
claident.tax[potential.std.id, "species"] <- as.character(blastn.std$V2)

# Output new claident tax table
write.csv(claident.tax, "03_TaxaSTDcombineOut/claident_tax_revise.csv", row.names = F)

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/03_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))