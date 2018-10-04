####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Figure: Bar plots
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library(phyloseq); packageVersion("phyloseq")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Make my palette
my.palette <- c("#FFC20A", # From scale_fill_igv() in "ggsci" package
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

# Call rare microbes figures
ps.phylum.b5 <- ps.phylum
# Rename samples
sample_data(ps.phylum.b5)$Sample_Name3 <- c(rep(c(sprintf("NoBeads S0%s", 1:5),
                                            sprintf("Beads S0%s", 1:5),
                                            sprintf("PowerSoil S0%s", 1:5)), 2))
b5 <- plot_bar(ps.phylum.b5, x = "Sample_Name3", fill = "phylum") +
  scale_fill_manual(values = my.palette) + 
  ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free")
(b5.levels <- levels(b5$data$phylum))
b5$data$phylum <- factor(b1$data$phylum, levels = c(b5.levels[-match(c("Others", "Undetermined"), b5.levels)], c("Others", "Undetermined")))
b5 <- b5 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                              axis.title.x = element_blank(),
                              axis.text = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 1))

# Rare microbes barplot
my.palette2 <- c(my.palette[c(3,5,7)], "#CDDEB7", my.palette[c(8:12,14)])
rl.ps.merge2 <- rl.ps.merge
sample_data(rl.ps.merge2)$Sample_Name3 <- c(rep(c(sprintf("Beads S0%s", 1:5),
                                                  sprintf("NoBeads S0%s", 1:5),
                                                  sprintf("PowerSoil S0%s", 1:5)), 2))
rl2 <- plot_bar(rl.ps.merge2, x = "Sample_Name3", fill = "phylum") + scale_fill_manual(values = my.palette2) +
  ylab("DNA (copies/ml water)") + facet_grid(Site ~ Method, scales = "free")
(rl2.levels <- levels(rl2$data$phylum))
rl2$data$phylum <- factor(rl2$data$phylum, levels = c(rl2.levels[-match(c("Undetermined"), rl2.levels)], c("Undetermined")))
rl2 <- rl2 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                                axis.title.x = element_blank(),
                                axis.text = element_text(colour = "black"),
                                axis.title = element_text(colour = "black"),
                                panel.grid.minor = element_blank(),
                                panel.grid.major = element_blank(),
                                axis.text.x = element_text(angle = 90, hjust = 1))

# Merge panels
Fig.Bar <- b5
Fig.Rare <- rl2

# Save and generate merged figure
save_plot("0_RawFigs/PNG/Fig_Bar.png", Fig.Bar, base_aspect_ratio = 2, ncol = 1, nrow = 2)
quartz(type="pdf", file="0_RawFigs/PDF/Fig_Bar.pdf", width=7, height=7); Fig.Bar; dev.off()

save_plot("0_RawFigs/PNG/Fig_Rare.png", Fig.Rare, base_aspect_ratio = 2, ncol = 1, nrow = 2)
quartz(type="pdf", file="0_RawFigs/PDF/Fig_Rare.pdf", width=7, height=7); Fig.Rare; dev.off()

write.csv(sample_data(ps.trim), "0_Table/CSV/TabS_SampleSheet.csv")
 