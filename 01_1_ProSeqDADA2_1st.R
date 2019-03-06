####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No.1.1 Sequence analysis by DADA2
#### 2018.6.22 Ushio
#### 2019.2.19 revised, Ushio
#### R 3.5.2
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("00_SessionInfo")
dir.create("01_ProSeqDADA2_Out")
dir.create("01_ProSeqDADA2_Out/objects")

# Load library and functions
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(dada2); packageVersion("dada2") # v1.10.1
library(ShortRead); packageVersion("ShortRead") # v1.40.0
library(ggplot2); packageVersion("ggplot2") # v3.1.0

# Load sequence reads
# (Generally follow DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html)
# To reproduce the results, fastq files should be downloaded from DDBJ DRA (Submission ID: DRA006959 for RMR-032)
# To reproduce the results, fastq files should be downloaded from DDBJ DRA (Submission ID: DRAXXXXXX for RMR-121)
sequence_id <- "RMR-032"
path <- sprintf("seqdata/%s", sequence_id)
fnFs <- sort(list.files(path, pattern=".forward.fastq", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern=".reverse.fastq", full.names = T)) # Reverse read files
# Get sample names
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 3) # This line should be changed depending on your computer environment so that the object "sample.names" represents your sample names.

# Visualize quality
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Performing filtering and trimming
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(215,160),
                     # Prokaryote primers; 515F = 19 bp, 806R = 20 bp
                     # Output sequences of Claident already trimmed Ns and primers
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=F,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Exclude 0 seq samples, rename filtFs and filtRs
if(length(sample.names[out[,2]<1 | out[,1]<1]) > 0){
  filtFs <- file.path(filt_path, paste0(sample.names[out[,2]>0 & out[,1]>0], "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names[out[,2]>0 & out[,1]>0], "_R_filt.fastq.gz"))
}

# Learn the error rates
min.nbases <- 5e+08
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min.nbases)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min.nbases)

# Visualize errors
#plotErrors(errF, nominalQ = T)
#plotErrors(errR, nominalQ = T)

# Dereplicatin
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[out[,2]>0 & out[,1]>0]
names(derepRs) <- sample.names[out[,2]>0 & out[,1]>0]

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Cutting unexpected length sequences
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(240,260)]
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads thourhg the pipeline
out2 <- out[out[,2]>0 & out[,1]>0,]
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2), rowSums(seqtab.nochim),  rowSums(seqtab.nochim)/out2[,1])
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "tabled2", "nonchim", "prop(last/first)")
rownames(track) <- sample.names[out[,2]>0 & out[,1]>0]
head(track)

# Reduce size
#rm(dadaFs)
rm(dadaRs)
rm(derepFs)
rm(derepRs)

# Save and output
seqs <- colnames(seqtab.nochim)
seqs.all <- matrix(rep(NA, 2*length(seqs)), ncol = 1)

# Output for Claident taxa assginment
seqs.out <- as.matrix(c(rbind(sprintf(">Taxa_seq1_%05d", 1:length(seqs)), seqs)), ncol=1)
write.table(seqs.out, "01_ProSeqDADA2_Out/ProkASV_1st_seqs.fa", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Save objects
saveRDS(seqtab, "01_ProSeqDADA2_Out/objects/seqtab_1st.obj")
saveRDS(seqtab2, "01_ProSeqDADA2_Out/objects/seqtab2_1st.obj")
saveRDS(seqtab.nochim, "01_ProSeqDADA2_Out/objects/seqtab_nochim_1st.obj")
saveRDS(track, "01_ProSeqDADA2_Out/objects/track_1st.obj")
saveRDS(seqs.out, "01_ProSeqDADA2_Out/objects/seqs_out_1st.obj")

# Save workspace and session info
save.image("01_ProSeqDADA2_Out/01_ProSeqDADA2_1st_seq_Out.RData")
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/01_SessionInfo_1st_seq_%s.txt", substr(Sys.time(), 1, 10)))
