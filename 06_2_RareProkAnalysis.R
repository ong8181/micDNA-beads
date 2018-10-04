####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### No.6 Method-specific prokaryote (ver.2: quantitative criteria)
#### 2018.6.22 Ushio
#### R 3.4.3
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("06_2_RareProkAnalysisOut")

# Load library and functions
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(cowplot); packageVersion("cowplot")
library(reshape2); packageVersion("reshape2")
library(tidyverse); packageVersion("tidyverse")
library(ggsci); packageVersion("ggsci")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
source("functions/F1_HelperFunctions.R")

# Load workspace
load("05_PhyloseqPlotOut/05_PhyloseqPlotOut.RData")

# Extract subsets
ps.trim.lake <- subset_samples(ps.trim, Site == "lake" & sample_nc == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps.trim.river <- subset_samples(ps.trim, Site == "river" & sample_nc == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)

pl.wb <- subset_samples(ps.trim.lake, Method == "Beads")
pl.wo <- subset_samples(ps.trim.lake, Method == "NoBeads")
pl.ds <- subset_samples(ps.trim.lake, Method == "PowerSoil")
pr.wb <- subset_samples(ps.trim.river, Method == "Beads")
pr.wo <- subset_samples(ps.trim.river, Method == "NoBeads")
pr.ds <- subset_samples(ps.trim.river, Method == "PowerSoil")

# Calculate means
lwb.m <- colMeans(otu_table(pl.wb))
lwo.m <- colMeans(otu_table(pl.wo))
lds.m <- colMeans(otu_table(pl.ds))
rwb.m <- colMeans(otu_table(pr.wb))
rwo.m <- colMeans(otu_table(pr.wo))
rds.m <- colMeans(otu_table(pr.ds))

# Calculate S.D.
lwb.cv <- apply(otu_table(pl.wb), 2, sd)/lwb.m
lwo.cv <- apply(otu_table(pl.wo), 2, sd)/lwo.m
lds.cv <- apply(otu_table(pl.ds), 2, sd)/lds.m
rwb.cv <- apply(otu_table(pr.wb), 2, sd)/rwb.m
rwo.cv <- apply(otu_table(pr.wo), 2, sd)/rwo.m
rds.cv <- apply(otu_table(pr.ds), 2, sd)/rds.m

# Extract taxa
# Method-specific taxa criteria ver.2
m.cond <- 2
cv.cond <- 1

lwb.cv[is.na(lwb.cv)] <-
  lwo.cv[is.na(lwo.cv)] <- 
  lds.cv[is.na(lds.cv)] <-
  rwb.cv[is.na(rwb.cv)] <-
  rwo.cv[is.na(rwo.cv)] <-
  rds.cv[is.na(rds.cv)] <- cv.cond + 1

lwb.cv.ok <- lwb.cv < cv.cond
lwo.cv.ok <- lwo.cv < cv.cond
lds.cv.ok <- lds.cv < cv.cond
rwb.cv.ok <- rwb.cv < cv.cond
rwo.cv.ok <- rwo.cv < cv.cond
rds.cv.ok <- rds.cv < cv.cond

lwb.m.ok <- lwb.m - m.cond * apply(rbind(lwo.m, lds.m), 2, max) > 0
lwo.m.ok <- lwo.m - m.cond * apply(rbind(lwb.m, lds.m), 2, max) > 0
lds.m.ok <- lds.m - m.cond * apply(rbind(lwb.m, lwo.m), 2, max) > 0
rwb.m.ok <- rwb.m - m.cond * apply(rbind(rwo.m, rds.m), 2, max) > 0
rwo.m.ok <- rwo.m - m.cond * apply(rbind(rwb.m, rds.m), 2, max) > 0
rds.m.ok <- rds.m - m.cond * apply(rbind(rwb.m, rwo.m), 2, max) > 0

sum(lwb.m.ok & lwb.cv.ok)
sum(lwo.m.ok & lwo.cv.ok)
sum(lds.m.ok & lds.cv.ok)
lwb.sp.tax <- names(lwb.m.ok)[lwb.m.ok & lwb.cv.ok]
lwo.sp.tax <- names(lwo.m.ok)[lwo.m.ok & lwo.cv.ok]
lds.sp.tax <- names(lds.m.ok)[lds.m.ok & lds.cv.ok]

sum(rwb.m.ok & rwb.cv.ok)
sum(rwo.m.ok & rwo.cv.ok)
sum(rds.m.ok & rds.cv.ok)
rwb.sp.tax <- names(rwb.m.ok)[rwb.m.ok & rwb.cv.ok]
rwo.sp.tax <- names(rwo.m.ok)[rwo.m.ok & rwo.cv.ok]
rds.sp.tax <- names(rds.m.ok)[rds.m.ok & rds.cv.ok]

# Prune taxa
l.wb.ps <- prune_taxa(lwb.sp.tax, pl.wb)
l.wo.ps <- prune_taxa(lwo.sp.tax, pl.wo)
l.ds.ps <- prune_taxa(lds.sp.tax, pl.ds)
r.wb.ps <- prune_taxa(rwb.sp.tax, pr.wb)
r.wo.ps <- prune_taxa(rwo.sp.tax, pr.wo)
r.ds.ps <- prune_taxa(rds.sp.tax, pr.ds)

# Merge and compile phyloseq objects
rl.ps.merge <- merge_phyloseq(l.wb.ps, l.wo.ps, l.ds.ps, r.wb.ps, r.wo.ps, r.ds.ps)
tax_table(rl.ps.merge)[tax_table(rl.ps.merge)[,"phylum"] == "", "phylum"] <- "Undetermined"
rl.ps.merge <- tax_glom(rl.ps.merge, taxrank = "phylum")

rl1 <- plot_bar(rl.ps.merge, fill = "phylum") + scale_fill_ucscgb() +
  ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free")
(rl1.levels <- levels(rl1$data$phylum))
rl1$data$phylum <- factor(rl1$data$phylum, levels = c(rl1.levels[-match(c("Undetermined"), rl1.levels)], c("Undetermined")))
rl1 <- rl1 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                         axis.text = element_text(colour = "black"),
                         axis.title = element_text(colour = "black"),
                         panel.grid.minor = element_blank(),
                         axis.text.x = element_text(angle = -90, hjust = 1))

# Check figure
OpenDev(8, 5)
rl1
dev.off()

# Tax sequence and taxa name
dir.create("06_2_RareProkAnalysisOut/1_Lake")
dir.create("06_2_RareProkAnalysisOut/2_River")

taxa.seq <- colnames(seqtab.conv2) # taxa sequences
names(taxa.seq) <- rownames(taxa.wo.std2) # taxa IDs, phylogeny, and seq lengths

l.wb.spc.fa <- taxa.seq[match(lwb.sp.tax, names(taxa.seq))]
l.wo.spc.fa <- taxa.seq[match(lwo.sp.tax, names(taxa.seq))]
l.ds.spc.fa <- taxa.seq[match(lds.sp.tax, names(taxa.seq))]
# Output method-specific sequences
writeXStringSet(DNAStringSet(l.wb.spc.fa), "06_2_RareProkAnalysisOut/1_Lake/Beads-specific-seq.fa")
writeXStringSet(DNAStringSet(l.wo.spc.fa), "06_2_RareProkAnalysisOut/1_Lake/NoBeads-specific-seq.fa")
writeXStringSet(DNAStringSet(l.ds.spc.fa), "06_2_RareProkAnalysisOut/1_Lake/PowerSoil-specific-seq.fa")
# Output taxa information
write.csv(taxa.wo.std2[lwb.sp.tax,],  "06_2_RareProkAnalysisOut/1_Lake/Beads-specific-taxa.csv")
write.csv(taxa.wo.std2[lwo.sp.tax,],  "06_2_RareProkAnalysisOut/1_Lake/NoBeads-specific-taxa.csv")
write.csv(taxa.wo.std2[lds.sp.tax,],  "06_2_RareProkAnalysisOut/1_Lake/PowerSoil-specific-taxa.csv")

r.wb.spc.fa <- taxa.seq[match(rwb.sp.tax, names(taxa.seq))]
r.wo.spc.fa <- taxa.seq[match(rwo.sp.tax, names(taxa.seq))]
r.ds.spc.fa <- taxa.seq[match(rds.sp.tax, names(taxa.seq))]
# Output method-specific sequences
writeXStringSet(DNAStringSet(r.wb.spc.fa), "06_2_RareProkAnalysisOut/2_River/Beads-specific-seq.fa")
writeXStringSet(DNAStringSet(r.wo.spc.fa), "06_2_RareProkAnalysisOut/2_River/NoBeads-specific-seq.fa")
writeXStringSet(DNAStringSet(r.ds.spc.fa), "06_2_RareProkAnalysisOut/2_River/PowerSoil-specific-seq.fa")
# Output taxa information
write.csv(taxa.wo.std2[rwb.sp.tax,],  "06_2_RareProkAnalysisOut/2_River/Beads-specific-taxa.csv")
write.csv(taxa.wo.std2[rwo.sp.tax,],  "06_2_RareProkAnalysisOut/2_River/NoBeads-specific-taxa.csv")
write.csv(taxa.wo.std2[rds.sp.tax,],  "06_2_RareProkAnalysisOut/2_River/PowerSoil-specific-taxa.csv")

# Save workspace
save.image("06_2_RareProkAnalysisOut/06_2_RareProkAnalysisOut.RData")

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/06_2_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
