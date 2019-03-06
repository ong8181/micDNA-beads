####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No.5 Visualization using phyloseq
#### 2018.6.22 Ushio
#### 2019.2.20 revised, Ushio
#### R 3.5.2
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("05_PhyloseqPlotOut")

# Load library and functions
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(cowplot); packageVersion("cowplot")
library(reshape2); packageVersion("reshape2")
library(tidyverse); packageVersion("tidyverse")
library(ggsci); packageVersion("ggsci")
library(phyloseq); packageVersion("phyloseq")

# Load helper functions
source("functions/F1_HelperFunctions.R")

# Load workspace
load("04_STDCheckOut/04_STDCheck_Out.RData")

# DADA2 summary objects
sample.sheet # sample data

#seqtab.nochim # Not converted
sample.no.nc <- sample.sheet[sample.sheet$sample_nc != "pcr_nc",]
taxa <- tax.claident # Recreate taxa object
taxa.wo.std2 <- as.matrix(taxa.wo.std[,2:ncol(taxa.wo.std)])
rownames(taxa.wo.std2) <- taxa.wo.std[,1]
tax.seq.length <- nchar(colnames(seqtab.conv))
seq.all <- colnames(seqtab.conv)
# Re-sort seqtab.conv
#seqtab.conv2 <- seqtab.conv[sort(rownames(seqtab.conv), decreasing = F),]
seqtab.conv2 <- seqtab.conv # Not sort (the order is already correct)

# Dimension check
dim(seqtab.conv2); dim(sample.no.nc); dim(taxa.wo.std2)

# Convert seqtab.conv2 to DNA copies/ml water
# 100 µl DNA extract/ XX ml water
filt_ml_tmp <- sample.no.nc$water_vol
filt_ml_tmp[is.na(filt_ml_tmp)] <- 1
seqtab.conv3 <- seqtab.conv2*100/filt_ml_tmp

# Import to phyloseq
rownames(sample.no.nc) <- rownames(seqtab.conv3)
colnames(seqtab.conv3) <- rownames(taxa.wo.std2)
taxa.wo.std2 <- cbind(taxa.wo.std2, seq.all, tax.seq.length)
ps <- phyloseq(otu_table(seqtab.conv3, taxa_are_rows=FALSE),
                   sample_data(sample.no.nc),
                   tax_table(taxa.wo.std2))

# Label Methods
m.level <- as.character(sample_data(ps)$Method)
m.level[m.level == "w_beads"] <- "Beads"
m.level[m.level == "wo_beads"] <- "NoBeads"
m.level[m.level == "destruction"] <- "PowerSoil"
sample_data(ps)$Method <- factor(m.level, levels = c("Beads", "NoBeads", "PowerSoil"))

# Sequence length distribution
seq.l.info <- hist(as.numeric(tax_table(ps)[,"tax.seq.length"]))
rbind(seq.l.info$breaks[1:(length(seq.l.info$breaks)-1)], seq.l.info$counts)

# Remove short and long sequences
# Sequences were already trimmed at DADA2 processing, but further triming can be possible at this stage
ps.trim <- prune_taxa(as.numeric(tax_table(ps)[,"tax.seq.length"]) > 239 &
                 as.numeric(tax_table(ps)[,"tax.seq.length"]) < 260,
                 ps)
ps.trim <- subset_taxa(ps.trim, superkingdom == "Bacteria" | superkingdom == "Archaea")

# Visualize general pattern
ps.trim.tmp <- ps.trim
tax_table(ps.trim.tmp)[tax_table(ps.trim.tmp)[,"phylum"] == "", "phylum"] <- "Undetermined"
ps.phylum <- tax_glom(ps.trim.tmp, taxrank = "phylum") %>%
  subset_samples(sample_nc == "sample") %>%
  merge_taxa(., taxa_names(.)[taxa_sums(.)/sum(taxa_sums(.)) < 0.001]) # merge rare taxa (<0.5%)
tax_table(ps.phylum)[is.na(tax_table(ps.phylum)[,"phylum"]), "phylum"] <- "Others"

# Abundance plot
b1 <- plot_bar(ps.phylum, fill = "phylum") + ggtitle("DNA extraction test") +
               scale_fill_ucscgb(alpha = 1) +
               ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free") #+ scale_fill_brewer(palette = "Accent")
(b1.levels <- levels(b1$data$phylum))
b1$data$phylum <- factor(b1$data$phylum, levels = c(b1.levels[-match(c("Others", "Undetermined"), b1.levels)], c("Others", "Undetermined")))
b1 <- b1 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                                axis.text = element_text(colour = "black"),
                                axis.title = element_text(colour = "black"),
                                panel.grid.minor = element_blank(),
                                axis.text.x = element_text(angle = -90, hjust = 1))

# Richness plot
ps.rich <- subset_samples(ps.trim, sample_nc == "sample")
otu_table(ps.rich) <- round(otu_table(ps.rich), digits = 0)
b2 <- plot_richness(ps.rich, x = "Method", measures = c("Observed")) + 
  ggtitle("DNA extraction test") + facet_wrap("Site") +
  ylim(0, 1100) + ylab("No. of ASVs")
b2$layers <- b2$layers[-1]
b2 <- b2 + geom_boxplot(width = 0.5, outlier.shape = NULL, outlier.size = 0, outlier.colour = "white")
b2 <- b2 + geom_jitter(size = 1.5, width = 0.1)

ps.rich100 <- ps.rich
otu_table(ps.rich100)[otu_table(ps.rich100) < 100] <- 0
b3 <- plot_richness(ps.rich100, x = "Method", measures = c("Observed")) + 
  ggtitle("DNA extraction test") + facet_wrap("Site") +
  ylim(0, 500) + ylab("No. of ASVs (>100 copies/ml water)")
b3$layers <- b3$layers[-1]
b3 <- b3 + geom_boxplot(width = 0.5, outlier.shape = NULL, outlier.size = 0, outlier.colour = "white")
b3 <- b3 + geom_jitter(size = 1.5, width = 0.1)

# NMDS for each site
ps.lake <- prune_samples(sample_sums(ps.trim) > 0, ps.trim) %>%
  subset_samples(., sample_nc == "sample" & Site == "lake")
ps.river <- prune_samples(sample_sums(ps.trim) > 0, ps.trim) %>%
  subset_samples(., sample_nc == "sample" & Site == "river")
ps.pond <- prune_samples(sample_sums(ps.trim) > 0, ps.trim) %>%
  subset_samples(., sample_nc == "sample" & Site == "pond")
ps.sea <- prune_samples(sample_sums(ps.trim) > 0, ps.trim) %>%
  subset_samples(., sample_nc == "sample" & Site == "sea")
nmds.bray.lake <- ordinate(ps.lake, "NMDS", "bray")
nmds.bray.river <- ordinate(ps.river, "NMDS", "bray")
nmds.bray.pond <- ordinate(ps.pond, "NMDS", "bray")
nmds.bray.sea <- ordinate(ps.sea, "NMDS", "bray")

p1 <- PlotStyle(plot_ordination(ps.lake, nmds.bray.lake, color="Method", shape = "Method", title="Lake samples, Bray NMDS"))
p1 <- p1 + geom_point(size = 2) + stat_ellipse(geom = "polygon", type="t", alpha=0.1, aes(fill=Method))
p1 <- p1 + scale_color_igv() + scale_fill_igv()
p2 <- PlotStyle(plot_ordination(ps.river, nmds.bray.river, color="Method", shape = "Method", title="River samples, Bray NMDS"))
p2 <- p2 + geom_point(size = 2) + stat_ellipse(geom = "polygon", type="t", alpha=0.1, aes(fill=Method))
p2 <- p2 + scale_color_igv() + scale_fill_igv()
p3 <- PlotStyle(plot_ordination(ps.pond, nmds.bray.pond, color="Method", shape = "Method", title="Pond samples, Bray NMDS"))
p3 <- p3 + geom_point(size = 2) + stat_ellipse(geom = "polygon", type="t", alpha=0.1, aes(fill=Method))
p3 <- p3 + scale_color_igv() + scale_fill_igv()
p4 <- PlotStyle(plot_ordination(ps.sea, nmds.bray.sea, color="Method", shape = "Method", title="Sea samples, Bray NMDS"))
p4 <- p4 + geom_point(size = 2) + stat_ellipse(geom = "polygon", type="t", alpha=0.1, aes(fill=Method))
p4 <- p4 + scale_color_igv() + scale_fill_igv()

p.legend <- get_legend(p1)

# NMDS for all sites
ps.all <- prune_samples(sample_sums(ps.trim) > 0, ps.trim) %>% subset_samples(., sample_nc == "sample")
nmds.bray.all <- ordinate(ps.all, "NMDS", "bray")

p5 <- PlotStyle(plot_ordination(ps.all, nmds.bray.all, color="Method", shape = "Method", title="All samples, Bray NMDS"))
p5 <- p5 + geom_point(size = 2) + stat_ellipse(geom = "polygon", type="t", alpha=0, aes(fill=Method:Site))
p5 <- p5 + scale_color_igv() + scale_fill_igv()

# Remove unnecessary objects
rm(seqtab.valid)

# Save and output results
save.image("05_PhyloseqPlotOut/05_PhyloseqPlotOut.RData")

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/05_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))