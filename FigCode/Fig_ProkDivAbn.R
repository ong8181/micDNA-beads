####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Figure: Diversity and abundance
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
#### 2019.2.20 Ushio (Figs generated on macOS)
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("scales"); packageVersion("scales")
library("ggsci"); packageVersion("ggsci")
library("phyloseq"); packageVersion("phyloseq")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")
source("0_FigFuncions/F1_HelperFunctions.R")
source("0_FigFuncions/F2_HelperFunctions.R")

# Adjust barplots
total.df <- data.frame(total_dna = sample_sums(ps.phylum),
                       Site = sample_data(ps.phylum)$Site,
                       Method = sample_data(ps.phylum)$Method)
# Statistical test
p.t0 <- ggplot(total.df, aes(x = Method, y = total_dna))
p.t0 <- p.t0 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.1)
p.t0 <- p.t0 + facet_wrap(~Site, scales = "free", ncol = 4)

# Lake site
TukeyHSD(aov(total_dna ~ Method, data = subset(total.df, Site == "lake")))
ann.text.pt1 <- data.frame(x = c(1.2, 2.2, 3.2),
                           y = c(150000, 110000, 100000),
                           lab = c("a","b","b"), Site = "lake")
p.t0 <- p.t0 + geom_text(data = ann.text.pt1, aes(x = x, y = y, label = lab), size = 5)

# River site
TukeyHSD(aov(total_dna ~ Method, data = subset(total.df, Site == "river")))

# Pond site
TukeyHSD(aov(total_dna ~ Method, data = subset(total.df, Site == "pond")))
ann.text.pt3 <- data.frame(x = c(1.4, 2.2, 3.2),
                           y = c(1500000, 900000, 700000),
                           lab = c("a","b","b"), Site = "pond")
p.t0 <- p.t0 + geom_text(data = ann.text.pt3, aes(x = x, y = y, label = lab), size = 5)

# Sea site
TukeyHSD(aov(total_dna ~ Method, data = subset(total.df, Site == "sea")))
ann.text.pt4 <- data.frame(x = c(1.4, 2.2, 3.2),
                           y = c(1000000, 550000, 350000),
                           lab = c("a","b","b"), Site = "sea")
p.t0 <- p.t0 + geom_text(data = ann.text.pt4, aes(x = x, y = y, label = lab), size = 5)
p.t0 <- p.t0 + scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})

# Define style
p.t0 <- PlotStyle4Pub(p.t0 + ylab("Total DNA\n(copies/ml water)"))

# ASV Diversity (All ASVs)
# Statistical test
TukeyHSD(aov(value ~ Method, data = subset(b2$data, Site == "lake")))
TukeyHSD(aov(value ~ Method, data = subset(b2$data, Site == "river")))
TukeyHSD(aov(value ~ Method, data = subset(b2$data, Site == "pond")))
TukeyHSD(aov(value ~ Method, data = subset(b2$data, Site == "sea")))

b2 <- PlotStyle4Pub(b2 + facet_grid(.~Site)) + ggtitle(NULL)

# ASV Diversity (ASVs > 100 copies/ml water)
# Statistical test
TukeyHSD(aov(value ~ Method, data = subset(b3$data, Site == "lake"))) # a, b, b
TukeyHSD(aov(value ~ Method, data = subset(b3$data, Site == "river"))) # N.S,
TukeyHSD(aov(value ~ Method, data = subset(b3$data, Site == "pond"))) # a, b, b
TukeyHSD(aov(value ~ Method, data = subset(b3$data, Site == "sea"))) # a, ab, b

b3 <- PlotStyle4Pub(b3) + facet_wrap(.~Site, ncol = 4) + ggtitle(NULL)

ann.text.b3.1 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(250, 190, 190),
                            lab = c("a","b","b"), Site = "lake")
ann.text.b3.3 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(355, 220, 240),
                            lab = c("a","b","b"), Site = "pond")
ann.text.b3.4 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(300, 230, 230),
                            lab = c("a","b","a"), Site = "sea")
b3 <- b3 + ylab("No. of ASVs\n (> 100 copies/ml water)")
b3 <- b3 + geom_text(data = ann.text.b3.1, aes(x = x, y = y, label = lab), size = 5)
b3 <- b3 + geom_text(data = ann.text.b3.3, aes(x = x, y = y, label = lab), size = 5)
b3 <- b3 + geom_text(data = ann.text.b3.4, aes(x = x, y = y, label = lab), size = 5)
b3 <- b3 + ylim(0,360)

# ASV Diversity (Repeatedly detected ASVs)
# Statistical test
TukeyHSD(aov(value ~ Method, data = subset(b4$data, Site == "lake")))
TukeyHSD(aov(value ~ Method, data = subset(b4$data, Site == "river")))
TukeyHSD(aov(value ~ Method, data = subset(b4$data, Site == "pond")))
TukeyHSD(aov(value ~ Method, data = subset(b4$data, Site == "sea")))

ann.text.b4.1 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(200, 160, 140),
                            lab = c("a","b","b"), Site = "lake")
ann.text.b4.2 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(160, 160, 120),
                            lab = c("a","a","b"), Site = "river")
ann.text.b4.3 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(230, 160, 180),
                            lab = c("a","b","b"), Site = "pond")
ann.text.b4.4 <- data.frame(x = c(1.2, 2.2, 3.2),
                            y = c(220, 160, 210),
                            lab = c("a","b","a"), Site = "sea")
b4 <- PlotStyle4Pub(b4 + facet_grid(.~Site)) + ggtitle(NULL)
b4 <- b4 + ylab("No. of repeatedly detected ASVs\n(detected at least 3 times)")
b4 <- b4 + geom_text(data = ann.text.b4.1, aes(x = x, y = y, label = lab), size = 5)
b4 <- b4 + geom_text(data = ann.text.b4.2, aes(x = x, y = y, label = lab), size = 5)
b4 <- b4 + geom_text(data = ann.text.b4.3, aes(x = x, y = y, label = lab), size = 5)
b4 <- b4 + geom_text(data = ann.text.b4.4, aes(x = x, y = y, label = lab), size = 5)
b4 <- b4 + ylim(0,270)

# Merge panels
Fig.Div <- plot_grid(p.t0, b2, b4,
                     ncol = 1,
                     labels = "auto",
                     align = "v",
                     axis = "l",
                     rel_heights = c(0.75,0.75,0.75))

# Save figures
quartz(type="pdf", file="0_RawFigs/PDF/Fig_Div.pdf", width=8, height=9); Fig.Div; dev.off()
saveRDS(Fig.Div, "0_RawFigs/Robj/Fig_Div.obj")
saveRDS(p.t0, "0_RawFigs/Robj/Fig_Div_a.obj")
saveRDS(b2, "0_RawFigs/Robj/Fig_Div_b.obj")
saveRDS(b4, "0_RawFigs/Robj/Fig_Div_c.obj")

quartz(type="pdf", file="0_RawFigs/PDF/Fig_Div100.pdf", width=8, height=3); b3; dev.off()
saveRDS(b3, "0_RawFigs/Robj/Fig_Div100.obj")
