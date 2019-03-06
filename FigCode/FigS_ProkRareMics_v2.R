####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Supplementary Figure: Bar plots and method-specific prokaryote
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
#### 2019.2.21 revised Ushio
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("scales"); packageVersion("scales")
library("phyloseq"); packageVersion("phyloseq")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_2_RareProkAnalysisOut/06_2_RareProkAnalysisOut.RData")
# Load my palette
source("0_FigFuncions/F2_HelperFunctions.R")

# Rare microbes barplot
sample_data(all.ps.merge)$Sample_Name3 <- c(rep(c(sprintf("Beads 0%s", 1:5),
                                                  sprintf("NoBeads 0%s", 1:5),
                                                  sprintf("PowerSoil 0%s", 1:5)), 4))
rl3 <- plot_bar(all.ps.merge, x = "Sample_Name3", fill = "phylum") + 
  ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free")
(rl3.levels <- levels(rl3$data$phylum))
rl3$data$phylum <- factor(rl3$data$phylum, levels = c(rl3.levels[-match(c("Undetermined"), rl3.levels)], c("Undetermined")))
rl3 <- rl3 + scale_fill_manual(values = ColorAssigner(levels(rl3$data$phylum)))
rl3 <- rl3 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                                axis.title.x = element_blank(),
                                axis.text = element_text(colour = "black"),
                                axis.title = element_text(colour = "black"),
                                panel.grid.minor = element_blank(),
                                panel.grid.major = element_blank(),
                                axis.text.x = element_text(angle = 90, hjust = 1))
rl3 <- rl3 + scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})

# Merge panels
Fig.Rare2 <- rl3

# Save and generate merged figure
quartz(type="pdf", file="0_RawFigs/PDF/FigS_Rare2.pdf", width=7, height=9); Fig.Rare2; dev.off()
saveRDS(Fig.Rare2, "0_RawFigs/Robj/FigS_Rare2.obj")
