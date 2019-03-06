####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Figure: Standard seqeunce check
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
#### 2019.2.20 revised, Ushio
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("scales"); packageVersion("scales")
library("ggsci"); packageVersion("ggsci")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")
source("0_FigFuncions/F2_HelperFunctions.R")

# Call STD data
no.pcrnc.data <- subset(sample.sheet, sample_nc != "pcr_nc")
no.pcrnc.data$r2 <- r2.summary
no.pcrnc.data$slope <- coef.summary
no.pcrnc.data$sample_nc2 <- factor(rep(c(rep(c(rep("Sample", 5), "Field NC"), 6),
                                     "PCR NC w/ STD", "PCR NC w/ STD"), 2),
                                   levels = c("Sample", "Field NC", "PCR NC w/ STD"))

k1 <- ggplot(no.pcrnc.data, aes(x = sample_nc2, y = r2, colour = sample_nc2))
k1 <- k1 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2)
k1 <- k1 + scale_color_d3() + xlab(NULL) + ylab(expression(R^2))

k2 <- ggplot(no.pcrnc.data, aes(x = sample_nc2, y = slope, colour = sample_nc2))
k2 <- k2 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2)
k2 <- k2 + scale_color_d3() + xlab(NULL) + ylab("Regression slope")

nonc.data <- subset(no.pcrnc.data, sample_nc2 == "Sample")
nonc.data$Method <- as.character(nonc.data$Method)
nonc.data$Method[nonc.data$Method == "w_beads"] <- "Beads"
nonc.data$Method[nonc.data$Method == "wo_beads"] <- "NoBeads"
nonc.data$Method[nonc.data$Method == "destruction"] <- "PowerSoil"
nonc.data$Method <- factor(nonc.data$Method, levels = c("Beads", "NoBeads", "PowerSoil"))
k3 <- ggplot(nonc.data, aes(x = Method, y = slope, colour = Method, facet = Site))
k3 <- k3 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2)
k3 <- k3 + facet_grid(~Site) + scale_color_d3() + xlab(NULL) + ylab("Regression slope")
k3 <- k3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Sample data
sample.data <- subset(no.pcrnc.data, sample_nc == "sample")
sample.data$sample_nc3 <- factor(rep(c(rep("NoBeads", 5),
                                       rep("Beads", 5),
                                       rep("PowerSoil", 5)), 4),
                                 levels = c("Beads", "NoBeads", "PowerSoil"))

k3 <- ggplot(sample.data, aes(x = sample_nc3, y = r2, colour = sample_nc3))
k3 <- k3 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2)
k3 <- k3 + facet_grid(.~Site)
k3 <- k3 + scale_color_d3() + xlab(NULL) + ylab(expression(R^2)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Change histogram colors
h1 <- g1 + scale_fill_manual(values = c("black", "gray70")) + xlim(0.6,1.05)
h2 <- g2 + scale_fill_manual(values = c("black", "gray70"))

h3 <- ggplot(slope.summary1, aes(x = copy, y = value, shape = variable, group = variable))
h3 <- h3 + geom_point(size = 2, color = "gray20")
h3 <- h3 + geom_smooth(method = "lm", size = 0.5, se = F, color = "gray50")
h3 <- h3 + xlab(expression(paste("Standard DNA copies (", mu, l^{-1}, ")"))) + ylab("Standard DNA reads")
h3 <- PlotStyle(h3) #+ theme(legend.position = c(0.25, 0.75), legend.text = element_text(size = 7), legend.title = element_text(size = 8))
h3 <- h3 + scale_x_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})
h3 <- h3 + scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})

# Generate combined figure
Fig.STD <- plot_grid(h3, h1, nrow = 1, labels = "auto", rel_widths = c(1, 0.8))

# Save figures
quartz(type="pdf", file="0_RawFigs/PDF/Fig_STD.pdf", width=9, height=3); Fig.STD; dev.off()
saveRDS(Fig.STD, "0_RawFigs/Robj/Fig_STD.obj")
saveRDS(h1, "0_RawFigs/Robj/Fig_STD_a.obj")
saveRDS(h3, "0_RawFigs/Robj/Fig_STD_b.obj")
saveRDS(h2, "0_RawFigs/Robj/Fig_STD_c.obj")

quartz(type="pdf", file="0_RawFigs/PDF/Fig_STDRegression.pdf", width=8, height=4)
k3; dev.off()
