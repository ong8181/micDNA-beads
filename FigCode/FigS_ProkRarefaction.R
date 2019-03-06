####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Supplementary Figure: Rarefaction curver
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
#### 2019.2.20 revised Ushio
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("ggsci"); packageVersion("ggsci")
library("phyloseq"); packageVersion("phyloseq")
library("reshape2"); packageVersion("reshape2")
library("ape"); packageVersion("ape")
library("gridExtra"); packageVersion("gridExtra")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# downloaded from https://github.com/mahendra-mariadassou/phyloseq-extended
set.seed(8181)

## Load helper functions
scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)
for (url in urls) source(url)

# Illustrate rarefaction curve
seqtab.wo.std <- seqtab.nochim[,-n.std.seq]
no.std.taxa <- tax.claident[-n.std.seq,]
dim(seqtab.wo.std); dim(no.std.taxa); dim(sample.sheet)

# Import to phyloseq
rownames(sample.sheet) <- rownames(seqtab.wo.std) <- sample.sheet$Sample_Name2
colnames(seqtab.wo.std) <- rownames(no.std.taxa) <- no.std.taxa$query

ps.rare0 <- phyloseq(otu_table(seqtab.wo.std, taxa_are_rows=FALSE),
                     sample_data(sample.sheet),
                     tax_table(as.matrix(no.std.taxa)))
ps.rare <- subset_samples(ps.rare0, sample_nc == "sample")
ps.rare <- subset_taxa(ps.rare, superkingdom == "Bacteria" | superkingdom == "Archaea")
sample_data(ps.rare)$Method <- rep(c(rep("NoBeads", 5), rep("Beads", 5), rep("PowerSoil", 5)), 4)
sample_data(ps.rare)$Method <- factor(sample_data(ps.rare)$Method, levels = c("Beads", "NoBeads", "PowerSoil"))

FigS.rarefy <- ggrare(ps.rare, step = 500, color = "Method", se = FALSE)
FigS.rarefy <- FigS.rarefy + facet_wrap(~Site, scales = "free") + scale_color_d3()
FigS.rarefy <- FigS.rarefy + xlab("Sequence reads") + ylab("No. of ASVs")
FigS.rarefy <- FigS.rarefy + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                              axis.text = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank())

quartz(type="pdf", file="0_RawFigs/PDF/FigS_Rarefy.pdf", width=7, height=5); FigS.rarefy; dev.off()
saveRDS(FigS.rarefy, "0_RawFigs/Robj/FigS_Rarefy.obj")
