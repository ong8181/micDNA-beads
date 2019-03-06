####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No.6.3 Field NC evaluations
#### 2018.6.22 Ushio
#### 2019.2.25 revised, Ushio
#### R 3.5.2
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("06_3_FieldNCEvaluationOut")

# Load library and functions
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(cowplot); packageVersion("cowplot")
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
source("functions/F1_HelperFunctions.R")

# Load workspace
load("05_PhyloseqPlotOut/05_PhyloseqPlotOut.RData")

# Extract subsets
ps.fnc.all <- tax_glom(ps.trim.tmp, taxrank = "phylum") %>%
  subset_samples(sample_nc == "field_nc") %>%
  merge_taxa(., taxa_names(.)[taxa_sums(.)/sum(taxa_sums(.)) < 0.001]) # merge rare taxa (<0.5%)
tax_table(ps.fnc.all)[is.na(tax_table(ps.fnc.all)[,"phylum"]), "phylum"] <- "Others"
sample_data(ps.fnc.all)$Sample_Name3 <- c(rep(c("Beads","NoBeads","PowerSoil"), 4))

f1 <- plot_bar(ps.fnc.all, x = "Method", y = "Abundance", fill = "phylum") +
  ylab("DNA (copies/ml water)") + facet_grid( ~ Site, scales = "free")

# Save workspace
save.image("06_3_FieldNCEvaluationOut/06_3_FieldNCEvaluationOut.RData")

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/06_3_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
