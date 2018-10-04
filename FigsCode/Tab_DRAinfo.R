####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Table: DRA information
#### 2018.6.22 Ushio (Run on Mac OSX)
####

# Load library and workspace
library("dada2"); packageVersion("dada2") # 2.2.1, 2018.5.10
load("../01_ProSeqDADA2_r32Out/01_ProSeqDADA2_r32_AllOut.RData")

# Information for DRA submission
nominal.summary <- data.frame(NULL)

for(i in 1:length(mergers)){
  seq.length <- nchar(mergers[[i]]$sequence)
  seq.abundance <- mergers[[i]]$abundance
  
  all.accept <- all(mergers[[i]]$accept)
  seq.length.mean <- sum(seq.length * seq.abundance) / sum(seq.abundance)
  seq.length.sd <- sqrt(sum((seq.length - seq.length.mean)^2 * seq.abundance / sum(seq.abundance)))

  nominal.temp <- data.frame(sample_name = sample.names[i],
                             nominal.length = round(seq.length.mean, 0),
                             nominal.stdev = round(seq.length.sd, 3)
                             )
  nominal.summary <- rbind(nominal.summary, nominal.temp)
}