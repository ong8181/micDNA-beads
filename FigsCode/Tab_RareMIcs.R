####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Table: ASV detection rate table
#### 2018.6.22 Ushio (Run on Mac OSX)
####

# Load library
library(phyloseq)

# Load workspace
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Write cav files
# Lake samples
lsum.df <- as.data.frame(lsum)
rownames(lsum.df) <- c("Beads", "NoBeads", "PowerSoil")
lsum.df$site <- "Lake"

# River samples
rsum.df <- as.data.frame(rsum)
rownames(rsum.df) <- c("Beads", "NoBeads", "PowerSoil")
rsum.df$site <- "River"

sum.df <- rbind(lsum.df, rsum.df)
write.csv(sum.df, "0_Table/CSV/Tab_ProkASVsummary.csv")

# Calsulate abundance proportions
GetRA <- function(phyloseq_obj, n_rep){
  tax.reps <- taxa_names(phyloseq_obj)[colSums(otu_table(phyloseq_obj) > 0) == n_rep]
  return(sum(sample_sums(prune_taxa(tax.reps, phyloseq_obj))))
}

GetAllRA <- function(phyloseq_obj){
  RA0 <- GetRA(phyloseq_obj, 0)/sum(sample_sums(phyloseq_obj)) * 100
  RA1 <- GetRA(phyloseq_obj, 1)/sum(sample_sums(phyloseq_obj)) * 100
  RA2 <- GetRA(phyloseq_obj, 2)/sum(sample_sums(phyloseq_obj)) * 100
  RA3 <- GetRA(phyloseq_obj, 3)/sum(sample_sums(phyloseq_obj)) * 100
  RA4 <- GetRA(phyloseq_obj, 4)/sum(sample_sums(phyloseq_obj)) * 100
  RA5 <- GetRA(phyloseq_obj, 5)/sum(sample_sums(phyloseq_obj)) * 100
  return(data.frame(RA0 = RA0, RA1 = RA1, RA2 = RA2, RA3 = RA3, RA4 = RA4, RA5 = RA5))
}

lwb.ra <- GetAllRA(pl.wb)
lwo.ra <- GetAllRA(pl.wo)
lds.ra <- GetAllRA(pl.ds)
lsum.ra <- rbind(lwb.ra, lwo.ra, lds.ra)
rownames(lsum.ra) <- c("Beads", "NoBeads", "PowerSoil")
rwb.ra <- GetAllRA(pr.wb)
rwo.ra <- GetAllRA(pr.wo)
rds.ra <- GetAllRA(pr.ds)
rsum.ra <- rbind(rwb.ra, rwo.ra, rds.ra)
rownames(rsum.ra) <- c("Beads", "NoBeads", "PowerSoil")
# ASV proportion summary
lsum.ra; rsum.ra

# Lake samples
lsum.df.ra <- as.data.frame(lsum.ra)
rownames(lsum.df.ra) <- c("Beads", "NoBeads", "PowerSoil")
lsum.df.ra$site <- "Lake"

# River samples
rsum.df.ra <- as.data.frame(rsum.ra)
rownames(rsum.df.ra) <- c("Beads", "NoBeads", "PowerSoil")
rsum.df.ra$site <- "River"
sum.df.ra <- rbind(lsum.df.ra, rsum.df.ra)

write.csv(sum.df.ra, "0_Table/CSV/Tab_ProkASVPropsummary.csv")
