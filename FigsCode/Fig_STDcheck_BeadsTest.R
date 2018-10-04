####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Figure: Standard seqeunce check
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("ggsci"); packageVersion("ggsci")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Call STD data
no.pcrnc.data <- subset(sample.sheet, sample_nc != "pcr_nc")
no.pcrnc.data$r2 <- r2.summary
no.pcrnc.data$slope <- coef.summary
no.pcrnc.data$sample_nc2 <- factor(c(rep(c(rep("Sample", 5), "Field NC"), 6),
                                     "PCR NC w/ STD", "PCR NC w/ STD"),
                                   levels = c("Sample", "Field NC", "PCR NC w/ STD"))

k1 <- ggplot(no.pcrnc.data, aes(x = sample_nc2, y = r2, colour = sample_nc2))
k1 <- k1 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2)
k1 <- k1 + scale_color_d3() + xlab(NULL) + ylab(expression(R^2))

k2 <- ggplot(no.pcrnc.data, aes(x = sample_nc2, y = slope, colour = sample_nc2))
k2 <- k2 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2)
k2 <- k2 + scale_color_d3() + xlab(NULL) + ylab("Regression slope")

# Sample data
sample.data <- subset(no.pcrnc.data, sample_nc == "sample")
sample.data$sample_nc3 <- factor(rep(c(rep("NoBeads", 5),
                                       rep("Beads", 5),
                                       rep("PowerSoil", 5)), 2),
                                 levels = c("Beads", "NoBeads", "PowerSoil"))

k3 <- ggplot(sample.data, aes(x = sample_nc3, y = r2, colour = sample_nc3))
k3 <- k3 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2)
k3 <- k3 + facet_grid(.~Site)
k3 <- k3 + scale_color_d3() + xlab(NULL) + ylab(expression(R^2))

# Generate combined figure
Fig.STD <- plot_grid(g1, g2, nrow = 1, labels = "auto")

save_plot("0_RawFigs/PNG/Fig_STD.png", Fig.STD, base_aspect_ratio = 1, ncol = 2, nrow = 1)
quartz(type="pdf", file="0_RawFigs/PDF/Fig_STD.pdf", width=6.8, height=3.5); Fig.STD; dev.off()
