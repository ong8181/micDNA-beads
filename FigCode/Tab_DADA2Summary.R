####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Table: DADA2 summary
#### 2018.6.22 Ushio (Run on Mac OSX)
#### 2019.2.21 revised Ushio
####

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Collect information
track
sample.sheet

summary.table <- data.frame(track,
                            sum_std_seqs = rowSums(new.std.table),
                            sum_nonstd_seqs = track[,"nonchim"] - 
                              rowSums(new.std.table)
                            )
# Save DADA2 summary table and sample sheet
write.csv(summary.table, "0_Table/CSV/Tab_DADA2_summary.csv", row.names = F)
write.csv(sample.sheet, "0_Table/CSV/Tab_SampleSheet.csv", row.names = F)
