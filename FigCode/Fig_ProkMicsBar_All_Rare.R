####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Figure: Bar plots
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
#### 2019.2.20 revised, Ushio
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("scales"); packageVersion("scales")
library("phyloseq"); packageVersion("phyloseq")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")
# Load my palette
source("0_FigFuncions/F2_HelperFunctions.R")

# Call rare microbes figures
ps.phylum.b5 <- ps.phylum
# Rename samples
sample_data(ps.phylum.b5)$Sample_Name3 <- c(rep(c(sprintf("NoBeads 0%s", 1:5),
                                            sprintf("Beads 0%s", 1:5),
                                            sprintf("PowerSoil 0%s", 1:5)), 4))
b5 <- plot_bar(ps.phylum.b5, x = "Sample_Name3", fill = "phylum") +
  ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free")
(b5.levels <- levels(b5$data$phylum))
b5$data$phylum <- factor(b5$data$phylum, levels = c(b5.levels[-match(c("Others", "Undetermined"), b5.levels)], c("Others", "Undetermined")))
b5 <- b5 + scale_fill_manual(values = ColorAssigner(levels(b5$data$phylum))) + scale_y_continuous(labels = scientific_format())
b5 <- b5 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                              axis.title.x = element_blank(),
                              axis.text = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 1))
b5 <- b5 + scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})


# Rare microbes barplot
rl.ps.merge2 <- rlps.ps.merge
sample_data(rl.ps.merge2)$Sample_Name3 <- c(rep(c(sprintf("Beads 0%s", 1:5),
                                                  sprintf("NoBeads 0%s", 1:5),
                                                  sprintf("PowerSoil 0%s", 1:5)), 4))
rl2 <- plot_bar(rl.ps.merge2, x = "Sample_Name3", fill = "phylum") +
  ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free")
(rl2.levels <- levels(rl2$data$phylum))
rl2$data$phylum <- factor(rl2$data$phylum, levels = c(rl2.levels[-match(c("Undetermined"), rl2.levels)], c("Undetermined")))
rl2 <- rl2 + scale_fill_manual(values = ColorAssigner(levels(rl2$data$phylum)))
rl2 <- rl2 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                                axis.title.x = element_blank(),
                                axis.text = element_text(colour = "black"),
                                axis.title = element_text(colour = "black"),
                                panel.grid.minor = element_blank(),
                                panel.grid.major = element_blank(),
                                axis.text.x = element_text(angle = 90, hjust = 1))

# Merge panelsf
Fig.Bar <- b5
Fig.Rare <- rl2

# Save and generate merged figure
quartz(type="pdf", file="0_RawFigs/PDF/Fig_Bar.pdf", width=7, height=9); Fig.Bar; dev.off()
saveRDS(Fig.Bar, "0_RawFigs/Robj/Fig_Bar.obj")

quartz(type="pdf", file="0_RawFigs/PDF/Fig_Rare.pdf", width=7, height=9); Fig.Rare; dev.off()
saveRDS(Fig.Rare, "0_RawFigs/Robj/Fig_Rare.obj")

write.csv(sample_data(ps.trim), "0_Table/CSV/TabS_SampleSheet.csv")
