####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No.1.3 Merge sequence tables by DADA2
#### 2018.6.22 Ushio
#### 2019.2.19 Ushio
#### R 3.5.2
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)

# Load library and functions
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")

# Load DADA2 results
seqtab.nochim.1 <- readRDS("01_ProSeqDADA2_Out/objects/seqtab_nochim_1st.obj")
seqtab.nochim.2 <- readRDS("01_ProSeqDADA2_Out/objects/seqtab_nochim_2nd.obj")
seqtab.1 <- readRDS("01_ProSeqDADA2_Out/objects/seqtab_1st.obj")
seqtab.2 <- readRDS("01_ProSeqDADA2_Out/objects/seqtab_2nd.obj")
seqtab2.1 <- readRDS("01_ProSeqDADA2_Out/objects/seqtab2_1st.obj")
seqtab2.2 <- readRDS("01_ProSeqDADA2_Out/objects/seqtab2_2nd.obj")
track.1 <- readRDS("01_ProSeqDADA2_Out/objects/track_1st.obj")
track.2 <- readRDS("01_ProSeqDADA2_Out/objects/track_2nd.obj")
seqs.out.1 <- readRDS("01_ProSeqDADA2_Out/objects/seqs_out_1st.obj")
seqs.out.2 <- readRDS("01_ProSeqDADA2_Out/objects/seqs_out_2nd.obj")

# Merge sequence tables
seqtab.nochim <- dada2::mergeSequenceTables(seqtab.nochim.1, seqtab.nochim.2, repeats = "error", orderBy = "abundance")

# Save and output
seqs <- colnames(seqtab.nochim)
seqs.all <- matrix(rep(NA, 2*length(seqs)), ncol = 1)

# Output for Claident taxa assginment
seqs.out <- as.matrix(c(rbind(sprintf(">Taxa_%05d", 1:length(seqs)), seqs)), ncol=1)
write.table(seqs.out, "01_ProSeqDADA2_Out/ProkASV_seqs.fa", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Output track
write.csv(rbind(track.1, track.2), "01_ProSeqDADA2_Out/track.csv")

# Save workspace and session info
save.image("01_ProSeqDADA2_Out/01_ProSeqDADA2_Merged_Out.RData")
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/01_SessionInfo_Merged_%s.txt", substr(Sys.time(), 1, 10)))
