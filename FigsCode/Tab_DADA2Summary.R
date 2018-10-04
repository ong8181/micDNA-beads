####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Table: DADA2 summary
#### 2018.6.22 Ushio (Run on Mac OSX)
####

# Load library and functions
library("dada2"); packageVersion("dada2")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Collect information
track
sample.names
fnRs

summary.table <- data.frame(sample_name = sample.names,
                            original_path = fnRs,
                            track,
                            sum_std_seqs = rowSums(new.std.table),
                            sum_nonstd_seqs = track[,"nonchim"] - 
                              rowSums(new.std.table)
                            )
# Add sample names
write.csv(summary.table, "0_Table/CSV/Tab_DADA2_summary.csv", row.names = F)
