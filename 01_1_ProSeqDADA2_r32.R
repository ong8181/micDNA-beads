####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### No.1 Sequence analysis by DADA2
#### 2018.6.22 Ushio
#### R 3.4.3
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("00_SessoinInfo")
dir.create("01_ProSeqDADA2_r32Out")
dir.create("DRA_RegisteredSeq")

# Load library and functions
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")

# Load sequence reads
# (Generally follow DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html)
# To reproduce the results, fastq files should be downloaded from DDBJ DRA (Submission ID: DRA006959)
path <- "DRA_RegisteredSeq" # Should be changed depending on your computer environment
fnFs <- sort(list.files(path, pattern=".forward.fastq", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern=".reverse.fastq", full.names = T)) # Reverse read files
# Get sample names
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 6) # This line should be changed depending on your computer environment so that the object "sample.names" represents your sample names.

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
                     #trimLeft = c(0, 0),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=F,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Exclude 0 seq samples, rename filtFs and filtRs
if(length(sample.names[out[,2]<1 | out[,1]<1]) > 0){
  filtFs <- file.path(filt_path, paste0(sample.names[out[,2]>0 & out[,1]>0], "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names[out[,2]>0 & out[,1]>0], "_R_filt.fastq.gz"))
}

# Learn the error rates
min.nreads <- 1e+06
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nreads = min.nreads)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nreads = min.nreads)

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
seqtab2 <- seqtab
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

# Output for Claident taxa assginment
seqs <- colnames(seqtab.nochim)
seqs.all <- matrix(rep(NA, 2*length(seqs)), ncol = 1)
seqs.out <- as.matrix(c(rbind(sprintf(">Taxa%05d", 1:length(seqs)), seqs)), ncol=1)
write.table(seqs.out, "01_ProSeqDADA2_r32Out/ProkASV_seqs.fa", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Reduce size
rm(dadaFs)
rm(dadaRs)
rm(derepFs)
rm(derepRs)

save.image("01_ProSeqDADA2_r32Out/01_ProSeqDADA2_r32Out.RData")

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/01_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
