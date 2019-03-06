####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### For adjusting and re-formatting all figures
#### 2019.2.27 Ushio (Figs generated on Mac OSX)
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("scales"); packageVersion("scales")
library("ggsci"); packageVersion("ggsci")
par(family="Arial")

# Creat output folder
dir.create("0_RawFigs/RefmtFigs")

# Load figure R objects
Fig.Div <- readRDS("0_RawFigs/Robj/Fig_Div.obj")
Fig.Div.a <- readRDS("0_RawFigs/Robj/Fig_Div_a.obj")
Fig.Div.b <- readRDS("0_RawFigs/Robj/Fig_Div_b.obj")
Fig.Div.c <- readRDS("0_RawFigs/Robj/Fig_Div_c.obj")
Fig.Div.100 <- readRDS("0_RawFigs/Robj/Fig_Div100.obj")

Fig.NMDS1 <- readRDS("0_RawFigs/Robj/Fig_NMDS1.obj")
Fig.NMDS1.a <- readRDS("0_RawFigs/Robj/Fig_NMDS1_a.obj")
Fig.NMDS1.b <- readRDS("0_RawFigs/Robj/Fig_NMDS1_b.obj")
Fig.NMDS1.c <- readRDS("0_RawFigs/Robj/Fig_NMDS1_c.obj")
Fig.NMDS1.d <- readRDS("0_RawFigs/Robj/Fig_NMDS1_d.obj")
Fig.NMDS.legend <- readRDS("0_RawFigs/Robj/Fig_NMDS_legend.obj")

Fig.NMDS2 <- readRDS("0_RawFigs/Robj/Fig_NMDS2.obj")
Fig.NMDS2.a <- readRDS("0_RawFigs/Robj/Fig_NMDS2_a.obj")
Fig.NMDS2.b <- readRDS("0_RawFigs/Robj/Fig_NMDS2_b.obj")

Fig.STD <- readRDS("0_RawFigs/Robj/Fig_STD.obj")
Fig.STD.a <- readRDS("0_RawFigs/Robj/Fig_STD_a.obj")
Fig.STD.b <- readRDS("0_RawFigs/Robj/Fig_STD_b.obj")

Fig.Bar <- readRDS("0_RawFigs/Robj/Fig_Bar.obj")
Fig.Rare <- readRDS("0_RawFigs/Robj/Fig_Rare.obj")
FigS.Beads <- readRDS("0_RawFigs/Robj/FigS_Beads.obj")
FigS.NanoDrop <- readRDS("0_RawFigs/Robj/FigS_NanoDrop.obj")
FigS.Rarefy <- readRDS("0_RawFigs/Robj/FigS_Rarefy.obj")
FigS.Rare2 <- readRDS("0_RawFigs/Robj/FigS_Rare2.obj")
FigS.NCBarProp <- readRDS("0_RawFigs/Robj/FigS_NCBarProp.obj")

#------------- Reformatting -------------#
# Figure 3
Fig.Div.a2 <- Fig.Div.a + theme_cowplot() + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
Fig.Div.b2 <- Fig.Div.b + theme_cowplot() + xlab(NULL) + ylab("No.\nof ASVs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
Fig.Div.c2 <- Fig.Div.c + theme_cowplot() + xlab(NULL) + ylab("No. of repeatedly\ndetected ASVs") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
Fig.Div.v2 <- plot_grid(Fig.Div.a2,
                        Fig.Div.b2,
                        Fig.Div.c2,
                        ncol = 1, scale = c(1,1,1), labels = "auto")
Fig.Div.100.v2 <- Fig.Div.100 + theme_cowplot() + xlab(NULL) + ylab("No. of ASVs > 100 copies/ml water") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

# Figure 4
Fig.NMDS1.a2 <- Fig.NMDS1.a + theme_cowplot() + theme(legend.position = "none") + ggtitle("Lake")
Fig.NMDS1.b2 <- Fig.NMDS1.b + theme_cowplot() + theme(legend.position = "none") + ggtitle("River")
Fig.NMDS1.c2 <- Fig.NMDS1.c + theme_cowplot() + theme(legend.position = "none") + ggtitle("Pond")
Fig.NMDS1.d2 <- Fig.NMDS1.d + theme_cowplot() + theme(legend.position = "none") + ggtitle("Sea")
Fig.NMDS1.v0 <- plot_grid(Fig.NMDS1.a2,
                          Fig.NMDS1.b2,
                          Fig.NMDS1.c2,
                          Fig.NMDS1.d2, ncol = 2,
                          align = "hv", labels = c("a", "b", "c", "d"),
                          scale = 1, label_x = 0.06)
Fig.NMDS1.v2 <- plot_grid(Fig.NMDS1.v0,
                          plot_grid(Fig.NMDS.legend, Fig.NMDS.legend, ncol = 1),
                          ncol = 2, rel_widths = c(2,0.5))

Fig.NMDS2.a2 <- Fig.NMDS2.a + theme_cowplot() + theme(legend.position = "none")
Fig.NMDS2.b2 <- Fig.NMDS2.b + theme_cowplot() + theme(legend.position = "none")
Fig.NMDS2.v0 <- plot_grid(Fig.NMDS2.a2,
                          Fig.NMDS2.b2,
                          ncol = 1, labels = c("a", "b"), align = "hv")
Fig.NMDS2.v2 <- plot_grid(Fig.NMDS2.v0,
                          plot_grid(Fig.NMDS.legend, Fig.NMDS.legend, ncol = 1),
                          ncol = 2,
                          rel_widths = c(1,0.5))

Fig.STD.a2 <- Fig.STD.a + theme_cowplot() + labs(fill="MiSeq run")
Fig.STD.b2 <- Fig.STD.b + theme_cowplot() + labs(shape="Regression\nslope") +
  scale_shape_manual(labels = c("max. slope", "med. slope", "min. slope"), values = c(16,17,15)) + geom_point(color="gray20")
Fig.STD.v2 <- plot_grid(Fig.STD.b2,
                        Fig.STD.a2, labels = "auto", align = "hv")

Fig.Bar.v2 <- Fig.Bar + theme_cowplot() + xlab(NULL) + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_bar(stat = "identity")
Fig.Rare.v2 <- Fig.Rare + theme_cowplot() + xlab(NULL) + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_bar(stat = "identity")
FigS.Beads.v2 <- FigS.Beads + theme_cowplot()
FigS.NanoDrop.v2 <- FigS.NanoDrop + theme_cowplot() + xlab(NULL) + scale_color_manual(values = c(1, 1, 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")
FigS.Rarefy.v2 <- FigS.Rarefy + theme_cowplot() + scale_color_igv() + theme(axis.text = element_text(size = 8))
FigS.Rare2.v2 <- FigS.Rare2 + theme_cowplot() + xlab(NULL) + geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_bar(stat = "identity")
FigS.NCBarProp.v2 <- FigS.NCBarProp + theme_cowplot() + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_bar(stat = "identity")

# Save reformat figures
quartz(type="pdf", file="0_RawFigs/RefmtFigs/Fig_Div.pdf", width=8, height=9)
Fig.Div.v2; dev.off()
quartz(type="pdf", file="0_RawFigs/RefmtFigs/Fig_Div100.pdf", width=8, height=4)
Fig.Div.100.v2; dev.off()

quartz(type="pdf", file="0_RawFigs/RefmtFigs/Fig_Bar.pdf", width=7, height=9)
Fig.Bar.v2; dev.off()
quartz(type="pdf", file="0_RawFigs/RefmtFigs/Fig_Rare.pdf", width=7, height=9)
Fig.Rare.v2; dev.off()

quartz(type="pdf", file="0_RawFigs/RefmtFigs/Fig_NMDS1.pdf", width=10, height=8.4)
Fig.NMDS1.v2; dev.off()
quartz(type="pdf", file="0_RawFigs/RefmtFigs/Fig_NMDS2.pdf", width=6, height=8)
Fig.NMDS2.v2; dev.off()

quartz(type="pdf", file="0_RawFigs/RefmtFigs/Fig_STD.pdf", width=9, height=3)
Fig.STD.v2; dev.off()
quartz(type="pdf", file="0_RawFigs/RefmtFigs/FigS_Beads.pdf", width=3, height=3.5)
FigS.Beads.v2; dev.off()
quartz(type="pdf", file="0_RawFigs/RefmtFigs/FigS_NanoDrop.pdf", width=7, height=3)
FigS.NanoDrop.v2; dev.off()
quartz(type="pdf", file="0_RawFigs/RefmtFigs/FigS_Rarefy.pdf", width=7, height=5)
FigS.Rarefy.v2; dev.off()
quartz(type="pdf", file="0_RawFigs/RefmtFigs/FigS_Rare2.pdf", width=7, height=9)
FigS.Rare2.v2; dev.off()
quartz(type="pdf", file="0_RawFigs/RefmtFigs/FigS_NCBarProp.pdf", width=8, height=5)
FigS.NCBarProp.v2; dev.off()

