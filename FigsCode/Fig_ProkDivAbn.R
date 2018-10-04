####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Figure: Diversity and abundance
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("ggsci"); packageVersion("ggsci")
library("phyloseq"); packageVersion("phyloseq")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Adjust barplots
total.df <- data.frame(total_dna = sample_sums(ps.phylum),
                       Site = sample_data(ps.phylum)$Site,
                       Method = sample_data(ps.phylum)$Method)
# Statistical test
TukeyHSD(aov(total_dna ~ Method, data = subset(total.df, Site == "lake")))
TukeyHSD(aov(total_dna ~ Method, data = subset(total.df, Site == "river")))

# Visualize
p.t1 <- ggplot(total.df, aes(x = Method, y = total_dna))
p.t1 <- p.t1 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.1) + facet_grid(.~Site, scales = "free")
p.t1 <- p.t1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + xlab(NULL) + ylab("Total DNA\n(copies/ml water)")
# Add stat results
ann.text.pt1 <- data.frame(x = c(1.2, 2.2, 3.2),
                       y = c(220000, 150000, 120000),
                       lab = c("a","b","b"), Site = "lake")
p.t1 <- p.t1 + geom_text(data = ann.text.pt1, aes(x = x, y = y, label = lab), size = 5)
p.t1 <- p.t1 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                                  plot.title = element_blank(),
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.text = element_text(colour = "black"),
                                  axis.title = element_text(colour = "black"),
                                  panel.grid.minor = element_blank(),
                                  panel.grid.major = element_blank())

# ASV Diversity
# Statistical test
TukeyHSD(aov(value ~ Method, data = subset(b2$data, Site == "lake")))
TukeyHSD(aov(value ~ Method, data = subset(b2$data, Site == "river")))
TukeyHSD(aov(value ~ Method, data = subset(b4$data, Site == "lake")))
TukeyHSD(aov(value ~ Method, data = subset(b4$data, Site == "river")))

b2 <- b2 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                              plot.title = element_blank(),
                              axis.title.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.text = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank())

ann.text.b4.1 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(220, 180, 180),
                            lab = c("a","b","b"), Site = "lake")
ann.text.b4.2 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(180, 180, 140),
                            lab = c("a","a","b"), Site = "river")
b4 <- b4 + ylab("No. of repeatedly detected ASVs\n(detected at least 3 times)")
b4 <- b4 + geom_text(data = ann.text.b4.1, aes(x = x, y = y, label = lab), size = 5)
b4 <- b4 + geom_text(data = ann.text.b4.2, aes(x = x, y = y, label = lab), size = 5)
b4 <- b4 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                              plot.title = element_blank(),
                              axis.title.x = element_blank(),
                              axis.text = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Merge panels
Fig.Div <- plot_grid(p.t1, b2, b4,
                     ncol = 1,
                     labels = "auto",
                     align = "v",
                     axis = "l",
                     rel_heights = c(0.75,0.75,0.925))

# Save figures
save_plot("0_RawFigs/PNG/Fig_Div.png", Fig.Div, base_aspect_ratio = 1.5, ncol = 1, nrow = 3)
quartz(type="pdf", file="0_RawFigs/PDF/Fig_Div.pdf", width=5, height=9); Fig.Div; dev.off()