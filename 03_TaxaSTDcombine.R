####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### No.3 STD sequences assignment
#### 2018.6.22 Ushio
#### R 3.4.3
####

# Load library
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(tidyverse); packageVersion("tidyverse")

# Create output directory
dir.create("03_TaxaSTDcombineOut")

# Load Claident taxa assignment
claident.tax <- read.delim("01_ProSeqDADA2_r32Out/ProkASV_seqs_classigntax")

# Load Blastn taxa assignment for STD sequences
# Select sequences with < 3 mismatches and length > 250
blastn.std <- read.table("02_ident_STD_BLASTnOut/ProkSTD_out.txt") %>%
  filter(.$V5 < 3 & .$V6 < 1 & .$V4 > 250)

# Check claident taxa assignments of the potential std sequences
potential.std.id <- match(blastn.std$V1, claident.tax$query)
#claident.tax[potential.std.id, "family"] # Not close to any existing tax --> OK
#claident.tax[potential.std.id, "species"] # Not close to any existing tax --> OK

# Replace STD taxa names with claident taxa assigment
claident.tax$species <- as.character(claident.tax$species)
claident.tax[potential.std.id, "species"] <- as.character(blastn.std$V2)

# Output new claident tax table
write.csv(claident.tax, "03_TaxaSTDcombineOut/claident.tax.revise.csv", row.names = F)

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/03_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


#### Sanity check ####
# Quality check (Claident taxa assignment)
n.tot <- nrow(claident.tax)
n.bac <- sum(claident.tax[,"superkingdom"] == "Bacteria")
n.arc <- sum(claident.tax[,"superkingdom"] == "Archaea")
n.euk <- sum(claident.tax[,"superkingdom"] == "Eukaryota")
c(n.bac, n.arc, n.euk)/n.tot # bac = 86.27%, arc = 0.06%, euk = 0%
1 - sum(n.bac, n.arc, n.euk)/n.tot # Undetermined = 7.7%

(n.std.tax <- sum(substr(claident.tax[,"species"], 1, 7) == "STD_pro", na.rm = T))
# STD sequences = 6
# ----> all OK!!