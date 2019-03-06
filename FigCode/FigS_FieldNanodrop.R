####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Supplementary Figure: Nanodrop results for DNA extraction test
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
#### 2019.2.20 revised Ushio
####

#### load data
d <- read.csv("0_nanodrop_data/field_nanodrop.csv")
d$site <- factor(d$site, levels = c("lake", "river", "pond", "sea"))

#### library load
library("ggplot2")
library("cowplot")
library("ggsci")

#### Select samples
d_sample <- subset(d, sample_nc=="sample")

#### Export to ggplot
p1 <- ggplot(d_sample, aes(x=method, y=DNA_ng_ml_water, color = method))
p1 <- p1 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.1, size=0.7)
p1 <- p1 + xlab(NULL) + ylab("DNA (ng/ml water)") + facet_wrap(.~site, scales = "free", ncol = 4)
p1 <- p1 + scale_color_d3()

# Statistical test
TukeyHSD(aov(DNA_ng_ml_water ~ method, data = subset(d_sample, site == "lake")))
TukeyHSD(aov(DNA_ng_ml_water ~ method, data = subset(d_sample, site == "river")))
TukeyHSD(aov(DNA_ng_ml_water ~ method, data = subset(d_sample, site == "pond")))
TukeyHSD(aov(DNA_ng_ml_water ~ method, data = subset(d_sample, site == "sea")))

ann.text.p1.1 <- data.frame(method = c("Beads", "NoBeads", "PowerSoil"),
                            DNA_ng_ml_water = c(35, 14, 12),
                            lab = c("a","b","b"), site = "lake")
ann.text.p1.2 <- data.frame(method = c("Beads", "NoBeads", "PowerSoil"),
                            DNA_ng_ml_water = c(9, 8, 5),
                            lab = c("a","a","b"), site = "river")
ann.text.p1.3 <- data.frame(method = c("Beads", "NoBeads", "PowerSoil"),
                            DNA_ng_ml_water = c(3.5, 3, 2.8),
                            lab = c("a","ab","b"), site = "pond")
ann.text.p1.4 <- data.frame(method = c("Beads", "NoBeads", "PowerSoil"),
                            DNA_ng_ml_water = c(2.3, 2.3, 2),
                            lab = c("ab","a","b"), site = "sea")

#p1 <- p1 + ylim(0,36)
p1 <- p1 + geom_text(data = ann.text.p1.1, aes(x = method, y = DNA_ng_ml_water, label = lab), size = 3, color = "black")
p1 <- p1 + geom_text(data = ann.text.p1.2, aes(x = method, y = DNA_ng_ml_water, label = lab), size = 3, color = "black")
p1 <- p1 + geom_text(data = ann.text.p1.3, aes(x = method, y = DNA_ng_ml_water, label = lab), size = 3, color = "black")
p1 <- p1 + geom_text(data = ann.text.p1.4, aes(x = method, y = DNA_ng_ml_water, label = lab), size = 3, color = "black")
p1 <- p1 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                              plot.title = element_blank(),
                              axis.title.x = element_blank(),
                              axis.text = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              #axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Output
quartz(type="pdf", file="0_RawFigs/PDF/FigS_NanoDrop.pdf", width=7, height=3); p1; dev.off()
saveRDS(p1, "0_RawFigs/Robj/FigS_NanoDrop.obj")
