####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### No.6 Method-specific prokaryote (ver.1: qualitative criteria)
#### 2018.6.22 Ushio
#### R 3.4.3
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("06_1_RareProkAnalysisOut")

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

# Calculate sums
# Converted reads were used to detect rare phylotypes
lwb.cal <- table(colSums(otu_table(pl.wb) > 0))
lwo.cal <- table(colSums(otu_table(pl.wo) > 0))
lds.cal <- table(colSums(otu_table(pl.ds) > 0))
lsum <- rbind(lwb.cal, lwo.cal, lds.cal)
rownames(lsum) <- c("Beads", "NoBeads", "PowerSoil")
rwb.cal <- table(colSums(otu_table(pr.wb) > 0))
rwo.cal <- table(colSums(otu_table(pr.wo) > 0))
rds.cal <- table(colSums(otu_table(pr.ds) > 0))
rsum <- rbind(rwb.cal, rwo.cal, rds.cal)
rownames(rsum) <- c("Beads", "NoBeads", "PowerSoil")
# ASV summary
lsum; rsum

# Extract only repeatedly deteted taxa and merge
# (Exclude ASVs with only one or two times detection)
detect.th <- 2
pl.wb.trim <- subset_taxa(pl.wb, colSums(otu_table(pl.wb) > 0) > detect.th)
pl.wo.trim <- subset_taxa(pl.wo, colSums(otu_table(pl.wo) > 0) > detect.th)
pl.ds.trim <- subset_taxa(pl.ds, colSums(otu_table(pl.ds) > 0) > detect.th)
pr.wb.trim <- subset_taxa(pr.wb, colSums(otu_table(pr.wb) > 0) > detect.th)
pr.wo.trim <- subset_taxa(pr.wo, colSums(otu_table(pr.wo) > 0) > detect.th)
pr.ds.trim <- subset_taxa(pr.ds, colSums(otu_table(pr.ds) > 0) > detect.th)

# Extract method-specific taxa
# Lake data
l.wb.tax <- taxa_names(pl.wb.trim)[taxa_sums(pl.wb.trim) > 0]
l.wo.tax <- taxa_names(pl.wo.trim)[taxa_sums(pl.wo.trim) > 0]
l.ds.tax <- taxa_names(pl.ds.trim)[taxa_sums(pl.ds.trim) > 0]
l.wb.spc <- which(is.na(match(l.wb.tax, c(l.wo.tax, l.ds.tax))))
l.wo.spc <- which(is.na(match(l.wo.tax, c(l.wb.tax, l.ds.tax))))
l.ds.spc <- which(is.na(match(l.ds.tax, c(l.wb.tax, l.wo.tax))))
# The number of method-specific taxa of lake samples
length(l.wb.spc); length(l.wo.spc); length(l.ds.spc)

# River data
r.wb.tax <- names(taxa_sums(pr.wb.trim)[taxa_sums(pr.wb.trim) > 0])
r.wo.tax <- names(taxa_sums(pr.wo.trim)[taxa_sums(pr.wo.trim) > 0])
r.ds.tax <- names(taxa_sums(pr.ds.trim)[taxa_sums(pr.ds.trim) > 0])
r.wb.spc <- which(is.na(match(r.wb.tax, c(r.wo.tax, r.ds.tax))))
r.wo.spc <- which(is.na(match(r.wo.tax, c(r.wb.tax, r.ds.tax))))
r.ds.spc <- which(is.na(match(r.ds.tax, c(r.wb.tax, r.wo.tax))))
# The number of method-specific taxa of river samples
length(r.wb.spc); length(r.wo.spc); length(r.ds.spc)

# Extract taxa of lake samples as a phyloseq object
l.wb.ps0 <- prune_taxa(l.wb.tax[l.wb.spc], pl.wb.trim)
l.wo.ps0 <- prune_taxa(l.wo.tax[l.wo.spc], pl.wo.trim)
l.ds.ps0 <- prune_taxa(l.ds.tax[l.ds.spc], pl.ds.trim)
# Extract taxa of river samples as a phyloseq object
r.wb.ps0 <- prune_taxa(r.wb.tax[r.wb.spc], pr.wb.trim)
r.wo.ps0 <- prune_taxa(r.wo.tax[r.wo.spc], pr.wo.trim)
r.ds.ps0 <- prune_taxa(r.ds.tax[r.ds.spc], pr.ds.trim)

# Extract only repeatedly deteted taxa and merge (4 or 5 times detection)
detect.th2 <- 3
l.wb.ps <- subset_taxa(l.wb.ps0, colSums(otu_table(l.wb.ps0) > 0) > detect.th2)
l.wo.ps <- subset_taxa(l.wo.ps0, colSums(otu_table(l.wo.ps0) > 0) > detect.th2)
l.ds.ps <- subset_taxa(l.ds.ps0, colSums(otu_table(l.ds.ps0) > 0) > detect.th2)
r.wb.ps <- subset_taxa(r.wb.ps0, colSums(otu_table(r.wb.ps0) > 0) > detect.th2)
r.wo.ps <- subset_taxa(r.wo.ps0, colSums(otu_table(r.wo.ps0) > 0) > detect.th2)
r.ds.ps <- subset_taxa(r.ds.ps0, colSums(otu_table(r.ds.ps0) > 0) > detect.th2)

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

# NMDS using repeatedly detected ASVs (one or two detections were excluded)
detect.th3 <- 2
pl.wb.trim2 <- subset_taxa(pl.wb, colSums(otu_table(pl.wb) > 0) > detect.th3)
pl.wo.trim2 <- subset_taxa(pl.wo, colSums(otu_table(pl.wo) > 0) > detect.th3)
pl.ds.trim2 <- subset_taxa(pl.ds, colSums(otu_table(pl.ds) > 0) > detect.th3)
pr.wb.trim2 <- subset_taxa(pr.wb, colSums(otu_table(pr.wb) > 0) > detect.th3)
pr.wo.trim2 <- subset_taxa(pr.wo, colSums(otu_table(pr.wo) > 0) > detect.th3)
pr.ds.trim2 <- subset_taxa(pr.ds, colSums(otu_table(pr.ds) > 0) > detect.th3)
pl.merge <- merge_phyloseq(pl.wb.trim2, pl.wo.trim2, pl.ds.trim2)
pr.merge <- merge_phyloseq(pr.wb.trim2, pr.wo.trim2, pr.ds.trim2)

nmds.bray.lake2 <- ordinate(pl.merge, "NMDS", "bray")
nmds.bray.river2 <- ordinate(pr.merge, "NMDS", "bray")

p3 <- PlotStyle(plot_ordination(pl.merge, nmds.bray.lake2, color = "Method", shape = "Method", title="Lake samples, Bray NMDS"))
p3 <- p3 + geom_point(size = 2) + stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=Method))
p3 <- p3 + scale_color_igv() + scale_fill_igv()
p4 <- PlotStyle(plot_ordination(pr.merge, nmds.bray.river2, color = "Method", shape = "Method", title="River samples, Bray NMDS"))
p4 <- p4 + geom_point(size = 2) + stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=Method))
p4 <- p4 + scale_color_igv() + scale_fill_igv()

p.legend2 <- get_legend(p3)

OpenDev(9, 4.2)
plot_grid(p3 + theme(legend.position = "none"),
          p4 + theme(legend.position = "none"),
          p.legend2, rel_widths = c(1,1,0.3), nrow = 1, labels = c("a","b",NULL))
dev.off()

# Richness plot
rich.merge <- merge_phyloseq(pl.merge, pr.merge)
otu_table(rich.merge) <- round(otu_table(rich.merge), digits = 0)
b4 <- plot_richness(rich.merge, x = "Method", measures = c("Observed")) + 
  ggtitle("DNA extraction test") + facet_wrap("Site") +
  ylim(0, 250) + ylab("No. of ASVs")
b4$layers <- b4$layers[-1]
b4 <- b4 + geom_boxplot(width = 0.5, outlier.shape = NULL, outlier.size = 0, outlier.colour = "white")
b4 <- b4 + geom_jitter(size = 1.5, width = 0.1)


# Tax sequence and taxa name
dir.create("06_1_RareProkAnalysisOut/1_Lake")
dir.create("06_1_RareProkAnalysisOut/2_River")

taxa.seq <- colnames(seqtab.conv2) # taxa sequences
names(taxa.seq) <- rownames(taxa.wo.std2) # taxa IDs, phylogeny, and seq lengths

l.wb.spc.fa <- taxa.seq[match(taxa_names(l.wb.ps), names(taxa.seq))]
l.wo.spc.fa <- taxa.seq[match(taxa_names(l.wo.ps), names(taxa.seq))]
l.ds.spc.fa <- taxa.seq[match(taxa_names(l.ds.ps), names(taxa.seq))]
# Output method-specific sequences
writeXStringSet(DNAStringSet(l.wb.spc.fa), "06_1_RareProkAnalysisOut/1_Lake/Beads-specific-seq.fa")
writeXStringSet(DNAStringSet(l.wo.spc.fa), "06_1_RareProkAnalysisOut/1_Lake/NoBeads-specific-seq.fa")
writeXStringSet(DNAStringSet(l.ds.spc.fa), "06_1_RareProkAnalysisOut/1_Lake/PowerSoil-specific-seq.fa")
# Output taxa information
write.csv(taxa.wo.std2[taxa_names(l.wb.ps),],  "06_1_RareProkAnalysisOut/1_Lake/Beads-specific-taxa.csv")
write.csv(taxa.wo.std2[taxa_names(l.wo.ps),],  "06_1_RareProkAnalysisOut/1_Lake/NoBeads-specific-taxa.csv")
write.csv(taxa.wo.std2[taxa_names(l.ds.ps),],  "06_1_RareProkAnalysisOut/1_Lake/PowerSoil-specific-taxa.csv")

r.wb.spc.fa <- taxa.seq[match(taxa_names(r.wb.ps), names(taxa.seq))]
r.wo.spc.fa <- taxa.seq[match(taxa_names(r.wo.ps), names(taxa.seq))]
r.ds.spc.fa <- taxa.seq[match(taxa_names(r.ds.ps), names(taxa.seq))]
# Output method-specific sequences
writeXStringSet(DNAStringSet(r.wb.spc.fa), "06_1_RareProkAnalysisOut/2_River/Beads-specific-seq.fa")
writeXStringSet(DNAStringSet(r.wo.spc.fa), "06_1_RareProkAnalysisOut/2_River/NoBeads-specific-seq.fa")
writeXStringSet(DNAStringSet(r.ds.spc.fa), "06_1_RareProkAnalysisOut/2_River/PowerSoil-specific-seq.fa")
# Output taxa information
write.csv(taxa.wo.std2[taxa_names(r.wb.ps),],  "06_1_RareProkAnalysisOut/2_River/Beads-specific-taxa.csv")
write.csv(taxa.wo.std2[taxa_names(r.wo.ps),],  "06_1_RareProkAnalysisOut/2_River/NoBeads-specific-taxa.csv")
write.csv(taxa.wo.std2[taxa_names(r.ds.ps),],  "06_1_RareProkAnalysisOut/2_River/PowerSoil-specific-taxa.csv")

# Save workspace
save.image("06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/06_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
