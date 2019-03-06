####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Supplementaly Figure: Bead amount test
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
####

# Load library and functions
library("ggplot2")
library("cowplot")
source("../functions/F1_HelperFunctions.R")
par(family="Arial")

# Load data
d <- read.csv("0_nanodrop_data/beadstest_nanodrop.csv")
d$filter_um <- as.factor(d$filter_um)
d$dna_ng_ml_water <- d$dna_ng_ul*100/d$filtered_water

# ggplot
# Scatter plot with smoothing line
d_sample <- subset(d, sample_nc == "sample")
p1 <- ggplot(d_sample, aes(x=beads_g, y=dna_ng_ml_water))
p1 <- p1 + geom_point(size=4) + xlab("Bead amount (g)") + ylab("DNA (ng/µl)")
p1 <- p1 + ggtitle("Test for beads amount for DNA extraction")
p1 <- p1 + geom_smooth(method="loess", se=F)

d_sample$beads_g_fac <- as.factor(d_sample$beads_g)
p2 <- ggplot(d_sample, aes(x=beads_g_fac, y=dna_ng_ml_water))
p2 <- p2 + geom_boxplot(colour = "gray10") + xlab("Bead amount (g)") + ylab("DNA (ng/ml water)")

# Adjust plot setting
FigS.Beads <- p2 + theme_bw() +
  geom_jitter(colour = "gray30", width = 0.1) +
  theme(panel.grid = element_blank()) + ylim(0,45)

quartz(type="pdf", file="0_RawFigs/PDF/FigS_Beads.pdf", width=3, height=3.5); FigS.Beads; dev.off()
saveRDS(FigS.Beads, "0_RawFigs/Robj/FigS_Beads.obj")
