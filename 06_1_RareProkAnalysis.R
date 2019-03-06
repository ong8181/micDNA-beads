####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No.6 Method-specific prokaryote (ver.1: qualitative criteria)
#### 2018.6.22 Ushio
#### 2019.2.20 revised, Ushio
#### R 3.5.2
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
ps.trim.pond <- subset_samples(ps.trim, Site == "pond" & sample_nc == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)
ps.trim.sea <- subset_samples(ps.trim, Site == "sea" & sample_nc == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)

pl.wb <- subset_samples(ps.trim.lake, Method == "Beads")
pl.wo <- subset_samples(ps.trim.lake, Method == "NoBeads")
pl.ds <- subset_samples(ps.trim.lake, Method == "PowerSoil")
pr.wb <- subset_samples(ps.trim.river, Method == "Beads")
pr.wo <- subset_samples(ps.trim.river, Method == "NoBeads")
pr.ds <- subset_samples(ps.trim.river, Method == "PowerSoil")
pp.wb <- subset_samples(ps.trim.pond, Method == "Beads")
pp.wo <- subset_samples(ps.trim.pond, Method == "NoBeads")
pp.ds <- subset_samples(ps.trim.pond, Method == "PowerSoil")
ps.wb <- subset_samples(ps.trim.sea, Method == "Beads")
ps.wo <- subset_samples(ps.trim.sea, Method == "NoBeads")
ps.ds <- subset_samples(ps.trim.sea, Method == "PowerSoil")

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
pwb.cal <- table(colSums(otu_table(pp.wb) > 0))
pwo.cal <- table(colSums(otu_table(pp.wo) > 0))
pds.cal <- table(colSums(otu_table(pp.ds) > 0))
psum <- rbind(pwb.cal, pwo.cal, pds.cal)
rownames(psum) <- c("Beads", "NoBeads", "PowerSoil")
swb.cal <- table(colSums(otu_table(ps.wb) > 0))
swo.cal <- table(colSums(otu_table(ps.wo) > 0))
sds.cal <- table(colSums(otu_table(ps.ds) > 0))
ssum <- rbind(swb.cal, swo.cal, sds.cal)
rownames(ssum) <- c("Beads", "NoBeads", "PowerSoil")
# ASV summary
lsum; rsum; psum; ssum

# Extract only repeatedly deteted taxa and merge
# (Exclude ASVs with only one or two times detection)
detect.th <- 2
pl.wb.trim <- subset_taxa(pl.wb, colSums(otu_table(pl.wb) > 0) > detect.th)
pl.wo.trim <- subset_taxa(pl.wo, colSums(otu_table(pl.wo) > 0) > detect.th)
pl.ds.trim <- subset_taxa(pl.ds, colSums(otu_table(pl.ds) > 0) > detect.th)
pr.wb.trim <- subset_taxa(pr.wb, colSums(otu_table(pr.wb) > 0) > detect.th)
pr.wo.trim <- subset_taxa(pr.wo, colSums(otu_table(pr.wo) > 0) > detect.th)
pr.ds.trim <- subset_taxa(pr.ds, colSums(otu_table(pr.ds) > 0) > detect.th)
pp.wb.trim <- subset_taxa(pp.wb, colSums(otu_table(pp.wb) > 0) > detect.th)
pp.wo.trim <- subset_taxa(pp.wo, colSums(otu_table(pp.wo) > 0) > detect.th)
pp.ds.trim <- subset_taxa(pp.ds, colSums(otu_table(pp.ds) > 0) > detect.th)
ps.wb.trim <- subset_taxa(ps.wb, colSums(otu_table(ps.wb) > 0) > detect.th)
ps.wo.trim <- subset_taxa(ps.wo, colSums(otu_table(ps.wo) > 0) > detect.th)
ps.ds.trim <- subset_taxa(ps.ds, colSums(otu_table(ps.ds) > 0) > detect.th)

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

# Pond data
p.wb.tax <- names(taxa_sums(pp.wb.trim)[taxa_sums(pp.wb.trim) > 0])
p.wo.tax <- names(taxa_sums(pp.wo.trim)[taxa_sums(pp.wo.trim) > 0])
p.ds.tax <- names(taxa_sums(pp.ds.trim)[taxa_sums(pp.ds.trim) > 0])
p.wb.spc <- which(is.na(match(p.wb.tax, c(p.wo.tax, p.ds.tax))))
p.wo.spc <- which(is.na(match(p.wo.tax, c(p.wb.tax, p.ds.tax))))
p.ds.spc <- which(is.na(match(p.ds.tax, c(p.wb.tax, p.wo.tax))))
# The number of method-specific taxa of pond samples
length(p.wb.spc); length(p.wo.spc); length(p.ds.spc)

# Sea data
s.wb.tax <- names(taxa_sums(ps.wb.trim)[taxa_sums(ps.wb.trim) > 0])
s.wo.tax <- names(taxa_sums(ps.wo.trim)[taxa_sums(ps.wo.trim) > 0])
s.ds.tax <- names(taxa_sums(ps.ds.trim)[taxa_sums(ps.ds.trim) > 0])
s.wb.spc <- which(is.na(match(s.wb.tax, c(s.wo.tax, s.ds.tax))))
s.wo.spc <- which(is.na(match(s.wo.tax, c(s.wb.tax, s.ds.tax))))
s.ds.spc <- which(is.na(match(s.ds.tax, c(s.wb.tax, s.wo.tax))))
# The number of method-specific taxa of sea samples
length(s.wb.spc); length(s.wo.spc); length(s.ds.spc)

# Extract taxa of lake samples as a phyloseq object
l.wb.ps0 <- prune_taxa(l.wb.tax[l.wb.spc], pl.wb.trim)
l.wo.ps0 <- prune_taxa(l.wo.tax[l.wo.spc], pl.wo.trim)
l.ds.ps0 <- prune_taxa(l.ds.tax[l.ds.spc], pl.ds.trim)
# Extract taxa of river samples as a phyloseq object
r.wb.ps0 <- prune_taxa(r.wb.tax[r.wb.spc], pr.wb.trim)
r.wo.ps0 <- prune_taxa(r.wo.tax[r.wo.spc], pr.wo.trim)
r.ds.ps0 <- prune_taxa(r.ds.tax[r.ds.spc], pr.ds.trim)
# Extract taxa of pond samples as a phyloseq object
p.wb.ps0 <- prune_taxa(p.wb.tax[p.wb.spc], pp.wb.trim)
p.wo.ps0 <- prune_taxa(p.wo.tax[p.wo.spc], pp.wo.trim)
p.ds.ps0 <- prune_taxa(p.ds.tax[p.ds.spc], pp.ds.trim)
# Extract taxa of sea samples as a phyloseq object
s.wb.ps0 <- prune_taxa(s.wb.tax[s.wb.spc], ps.wb.trim)
s.wo.ps0 <- prune_taxa(s.wo.tax[s.wo.spc], ps.wo.trim)
s.ds.ps0 <- prune_taxa(s.ds.tax[s.ds.spc], ps.ds.trim)


# Extract only repeatedly deteted taxa and merge (4 or 5 times detection)
detect.th2 <- 3
l.wb.ps <- subset_taxa(l.wb.ps0, colSums(otu_table(l.wb.ps0) > 0) > detect.th2)
l.wo.ps <- subset_taxa(l.wo.ps0, colSums(otu_table(l.wo.ps0) > 0) > detect.th2)
l.ds.ps <- subset_taxa(l.ds.ps0, colSums(otu_table(l.ds.ps0) > 0) > detect.th2)
r.wb.ps <- subset_taxa(r.wb.ps0, colSums(otu_table(r.wb.ps0) > 0) > detect.th2)
r.wo.ps <- subset_taxa(r.wo.ps0, colSums(otu_table(r.wo.ps0) > 0) > detect.th2)
r.ds.ps <- subset_taxa(r.ds.ps0, colSums(otu_table(r.ds.ps0) > 0) > detect.th2)
p.wb.ps <- subset_taxa(p.wb.ps0, colSums(otu_table(p.wb.ps0) > 0) > detect.th2)
p.wo.ps <- subset_taxa(p.wo.ps0, colSums(otu_table(p.wo.ps0) > 0) > detect.th2)
p.ds.ps <- subset_taxa(p.ds.ps0, colSums(otu_table(p.ds.ps0) > 0) > detect.th2)
s.wb.ps <- subset_taxa(s.wb.ps0, colSums(otu_table(s.wb.ps0) > 0) > detect.th2)
s.wo.ps <- subset_taxa(s.wo.ps0, colSums(otu_table(s.wo.ps0) > 0) > detect.th2)
s.ds.ps <- subset_taxa(s.ds.ps0, colSums(otu_table(s.ds.ps0) > 0) > detect.th2)

# Merge and compile phyloseq objects
rlps.ps.merge <- merge_phyloseq(l.wb.ps, l.wo.ps, l.ds.ps,
                                r.wb.ps, r.wo.ps, r.ds.ps,
                                p.wb.ps, p.wo.ps, p.ds.ps,
                                s.wb.ps, s.wo.ps, s.ds.ps)
tax_table(rlps.ps.merge)[tax_table(rlps.ps.merge)[,"phylum"] == "", "phylum"] <- "Undetermined"
rlps.ps.merge <- tax_glom(rlps.ps.merge, taxrank = "phylum")

rlps1 <- plot_bar(rlps.ps.merge, fill = "phylum") + scale_fill_ucscgb() +
  ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free")
(rlps1.levels <- levels(rlps1$data$phylum))
rlps1$data$phylum <- factor(rlps1$data$phylum, levels = c(rlps1.levels[-match(c("Undetermined"), rlps1.levels)], c("Undetermined")))
rlps1 <- rlps1 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                         axis.text = element_text(colour = "black"),
                         axis.title = element_text(colour = "black"),
                         panel.grid.minor = element_blank(),
                         axis.text.x = element_text(angle = -90, hjust = 1))

# NMDS using repeatedly detected ASVs (one or two detections were excluded)
detect.th3 <- 2
pl.wb.trim2 <- subset_taxa(pl.wb, colSums(otu_table(pl.wb) > 0) > detect.th3)
pl.wo.trim2 <- subset_taxa(pl.wo, colSums(otu_table(pl.wo) > 0) > detect.th3)
pl.ds.trim2 <- subset_taxa(pl.ds, colSums(otu_table(pl.ds) > 0) > detect.th3)
pr.wb.trim2 <- subset_taxa(pr.wb, colSums(otu_table(pr.wb) > 0) > detect.th3)
pr.wo.trim2 <- subset_taxa(pr.wo, colSums(otu_table(pr.wo) > 0) > detect.th3)
pr.ds.trim2 <- subset_taxa(pr.ds, colSums(otu_table(pr.ds) > 0) > detect.th3)
pp.wb.trim2 <- subset_taxa(pp.wb, colSums(otu_table(pp.wb) > 0) > detect.th3)
pp.wo.trim2 <- subset_taxa(pp.wo, colSums(otu_table(pp.wo) > 0) > detect.th3)
pp.ds.trim2 <- subset_taxa(pp.ds, colSums(otu_table(pp.ds) > 0) > detect.th3)
ps.wb.trim2 <- subset_taxa(ps.wb, colSums(otu_table(ps.wb) > 0) > detect.th3)
ps.wo.trim2 <- subset_taxa(ps.wo, colSums(otu_table(ps.wo) > 0) > detect.th3)
ps.ds.trim2 <- subset_taxa(ps.ds, colSums(otu_table(ps.ds) > 0) > detect.th3)
pl.merge <- merge_phyloseq(pl.wb.trim2, pl.wo.trim2, pl.ds.trim2)
pr.merge <- merge_phyloseq(pr.wb.trim2, pr.wo.trim2, pr.ds.trim2)
pp.merge <- merge_phyloseq(pp.wb.trim2, pp.wo.trim2, pp.ds.trim2)
ps.merge <- merge_phyloseq(ps.wb.trim2, ps.wo.trim2, ps.ds.trim2)

nmds.bray.lake2 <- ordinate(pl.merge, "NMDS", "bray")
nmds.bray.river2 <- ordinate(pr.merge, "NMDS", "bray")
nmds.bray.pond2 <- ordinate(pp.merge, "NMDS", "bray")
nmds.bray.sea2 <- ordinate(ps.merge, "NMDS", "bray")

q1 <- PlotStyle(plot_ordination(pl.merge, nmds.bray.lake2, color = "Method", shape = "Method", title="Lake samples, Bray NMDS"))
q1 <- q1 + geom_point(size = 2) + stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=Method))
q1 <- q1 + scale_color_igv() + scale_fill_igv()
q2 <- PlotStyle(plot_ordination(pr.merge, nmds.bray.river2, color = "Method", shape = "Method", title="River samples, Bray NMDS"))
q2 <- q2 + geom_point(size = 2) + stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=Method))
q2 <- q2 + scale_color_igv() + scale_fill_igv()
q3 <- PlotStyle(plot_ordination(pp.merge, nmds.bray.pond2, color = "Method", shape = "Method", title="Pond samples, Bray NMDS"))
q3 <- q3 + geom_point(size = 2) + stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=Method))
q3 <- q3 + scale_color_igv() + scale_fill_igv()
q4 <- PlotStyle(plot_ordination(ps.merge, nmds.bray.sea2, color = "Method", shape = "Method", title="Sea samples, Bray NMDS"))
q4 <- q4 + geom_point(size = 2) + stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=Method))
q4 <- q4 + scale_color_igv() + scale_fill_igv()

p.legend2 <- get_legend(q1)

# NMDS using all sites
all.merge <- merge_phyloseq(pl.merge, pr.merge, pp.merge, ps.merge)
nmds.bray.all2 <- ordinate(all.merge, "NMDS", "bray")

q5 <- PlotStyle(plot_ordination(all.merge, nmds.bray.all2, color="Method", shape = "Method", title="All samples, Bray NMDS"))
q5 <- q5 + geom_point(size = 2) + stat_ellipse(geom = "polygon", type="t", alpha=0, aes(fill=Method:Site))
q5 <- q5 + scale_color_igv() + scale_fill_igv()

# Richness plot
rich.merge <- merge_phyloseq(pl.merge, pr.merge, pp.merge, ps.merge)
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
dir.create("06_1_RareProkAnalysisOut/3_Pond")
dir.create("06_1_RareProkAnalysisOut/4_Sea")

taxa.seq <- colnames(seqtab.conv2) # taxa sequences
names(taxa.seq) <- rownames(taxa.wo.std2) # taxa IDs, phylogeny, and seq lengths

# Output lake sequendes
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

# Output river sequendes
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

# Output pond sequendes
p.wb.spc.fa <- taxa.seq[match(taxa_names(p.wb.ps), names(taxa.seq))]
p.wo.spc.fa <- taxa.seq[match(taxa_names(p.wo.ps), names(taxa.seq))]
p.ds.spc.fa <- taxa.seq[match(taxa_names(p.ds.ps), names(taxa.seq))]
# Output method-specific sequences
writeXStringSet(DNAStringSet(p.wb.spc.fa), "06_1_RareProkAnalysisOut/3_Pond/Beads-specific-seq.fa")
writeXStringSet(DNAStringSet(p.wo.spc.fa), "06_1_RareProkAnalysisOut/3_Pond/NoBeads-specific-seq.fa")
writeXStringSet(DNAStringSet(p.ds.spc.fa), "06_1_RareProkAnalysisOut/3_Pond/PowerSoil-specific-seq.fa")
# Output taxa information
write.csv(taxa.wo.std2[taxa_names(p.wb.ps),],  "06_1_RareProkAnalysisOut/3_Pond/Beads-specific-taxa.csv")
write.csv(taxa.wo.std2[taxa_names(p.wo.ps),],  "06_1_RareProkAnalysisOut/3_Pond/NoBeads-specific-taxa.csv")
write.csv(taxa.wo.std2[taxa_names(p.ds.ps),],  "06_1_RareProkAnalysisOut/3_Pond/PowerSoil-specific-taxa.csv")

# Output sea sequendes
s.wb.spc.fa <- taxa.seq[match(taxa_names(s.wb.ps), names(taxa.seq))]
s.wo.spc.fa <- taxa.seq[match(taxa_names(s.wo.ps), names(taxa.seq))]
s.ds.spc.fa <- taxa.seq[match(taxa_names(s.ds.ps), names(taxa.seq))]
# Output method-specific sequences
writeXStringSet(DNAStringSet(s.wb.spc.fa), "06_1_RareProkAnalysisOut/4_Sea/Beads-specific-seq.fa")
writeXStringSet(DNAStringSet(s.wo.spc.fa), "06_1_RareProkAnalysisOut/4_Sea/NoBeads-specific-seq.fa")
writeXStringSet(DNAStringSet(s.ds.spc.fa), "06_1_RareProkAnalysisOut/4_Sea/PowerSoil-specific-seq.fa")
# Output taxa information
write.csv(taxa.wo.std2[taxa_names(s.wb.ps),],  "06_1_RareProkAnalysisOut/4_Sea/Beads-specific-taxa.csv")
write.csv(taxa.wo.std2[taxa_names(s.wo.ps),],  "06_1_RareProkAnalysisOut/4_Sea/NoBeads-specific-taxa.csv")
write.csv(taxa.wo.std2[taxa_names(s.ds.ps),],  "06_1_RareProkAnalysisOut/4_Sea/PowerSoil-specific-taxa.csv")


# Save workspace
save.image("06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/06_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
