####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Supplementary Figure: Nanodrop results for DNA extraction test
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
####

#### load data
d <- read.csv("0_nanodrop_data/field_nanodrop.csv")

#### library load
library("ggplot2")
library("cowplot")
library("ggsci")

#### define function to adjust figure details
style_plot <-  function(ggobject){
  return(ggobject + theme_bw() + theme(axis.text.x = element_text(angle=0),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.text = element_text(size=12),
                                       axis.title = element_text(size=12),
                                       panel.background=element_rect(colour="black", fill=NA, size=0.8)))
}

#### Select samples
d_sample <- subset(d, sample_nc=="sample")

#### Export to ggplot
p1 <- ggplot(d_sample, aes(x=method, y=DNA_pg_ml_water, color = method))
p1 <- p1 + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.1)
p1 <- p1 + xlab(NULL) + ylab("DNA (pg/ml water)") + facet_grid(.~site)
p1 <- p1 + scale_color_igv()

# Statistical test
TukeyHSD(aov(DNA_pg_ml_water ~ method, data = subset(d_sample, site == "lake")))
TukeyHSD(aov(DNA_pg_ml_water ~ method, data = subset(d_sample, site == "river")))

ann.text.p1.1 <- data.frame(method = c("Beads", "NoBeads", "PowerSoil"),
                            DNA_pg_ml_water = c(350, 140, 120),
                            lab = c("a","b","b"), site = "lake")
ann.text.p1.2 <- data.frame(method = c("Beads", "NoBeads", "PowerSoil"),
                            DNA_pg_ml_water = c(110, 110, 80),
                            lab = c("a","a","b"), site = "river")
p1 <- p1 + ylim(0,360)
p1 <- p1 + geom_text(data = ann.text.p1.1, aes(x = method, y = DNA_pg_ml_water, label = lab), size = 5, color = "black")
p1 <- p1 + geom_text(data = ann.text.p1.2, aes(x = method, y = DNA_pg_ml_water, label = lab), size = 5, color = "black")
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
save_plot("0_RawFigs/PNG/FigS_NanoDrop.png", p1, base_aspect_ratio = 1.5, ncol = 1, nrow = 1)
quartz(type="pdf", file="0_RawFigs/PDF/FigS_NanoDrop.pdf", width=5, height=3); p1; dev.off()
