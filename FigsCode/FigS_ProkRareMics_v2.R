####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Supplementary Figure: Bar plots and method-specific prokaryote
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library(phyloseq); packageVersion("phyloseq")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_2_RareProkAnalysisOut/06_2_RareProkAnalysisOut.RData")

# Make my palette
my.palette <- c("#FFC20A", #"#5050FF", # scale_fill_igv() in "ggsci" package
                "#CE3D32",
                "#F0E685", #"#749B58",
                "#1A0099", #"#F0E685",
                "#99CC00", #"#466983",
                
                "#BA6338",
                "#5DB1DD",
                "#990080", #"#802268",
                "#466983", #"#6BD76B",
                "#4775FF", #"#D595A7",
                
                "#924822",
                "#D58F5C", #"#837B8D",
                "#749B58", #"#C75127",
                "#837B8D" #"#D58F5C"
)

# Rare microbes barplot
my.palette3 <- c(my.palette[c(1:3,5:7)], "#CDDEB7", my.palette[c(8:10,12,14)])
rl.ps.merge3 <- rl.ps.merge
sample_data(rl.ps.merge3)$Sample_Name3 <- c(rep(c(sprintf("Beads S0%s", 1:5),
                                                  sprintf("NoBeads S0%s", 1:5),
                                                  sprintf("PowerSoil S0%s", 1:5)), 2))
rl3 <- plot_bar(rl.ps.merge3, x = "Sample_Name3", fill = "phylum") + scale_fill_manual(values = my.palette3) +
  ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free")
(rl3.levels <- levels(rl3$data$phylum))
rl3$data$phylum <- factor(rl3$data$phylum, levels = c(rl3.levels[-match(c("Undetermined"), rl3.levels)], c("Undetermined")))
rl3 <- rl3 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                                axis.title.x = element_blank(),
                                axis.text = element_text(colour = "black"),
                                axis.title = element_text(colour = "black"),
                                panel.grid.minor = element_blank(),
                                panel.grid.major = element_blank(),
                                axis.text.x = element_text(angle = 90, hjust = 1))

# Merge panels
Fig.Rare2 <- rl3

# Save and generate merged figure
save_plot("0_RawFigs/PNG/FigS_Rare2.png", Fig.Rare2, base_aspect_ratio = 2, ncol = 1, nrow = 2)
quartz(type="pdf", file="0_RawFigs/PDF/FigS_Rare2.pdf", width=7, height=7); Fig.Rare2; dev.off()