####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Figure: NMDS
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("ggsci"); packageVersion("ggsci")
library("vegan"); packageVersion("vegan")
library("phyloseq"); packageVersion("phyloseq")
set.seed(8181)
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# PERMANOVA for NMDS data
# All data
permanova1 <- adonis(distance(ps.lake, method = "bray") ~ data.frame(sample_data(ps.lake))$Method, permutations = 9999)
# = 0.0112
permanova2 <- adonis(distance(ps.river, method = "bray") ~ data.frame(sample_data(ps.river))$Method, permutations = 9999)
# = 0.0111
# ASVs >= 3 data
permanova3 <- adonis(distance(pl.merge, method = "bray") ~ data.frame(sample_data(pl.merge))$Method, permutations = 9999)
# = 0.0042
permanova4 <- adonis(distance(pr.merge, method = "bray") ~ data.frame(sample_data(pr.merge))$Method, permutations = 9999)
# < 1e-04

# Add P value and stress values
strs1 <- round(nmds.bray.lake$stress, 4)
strs2 <- round(nmds.bray.lake2$stress, 4)
strs3 <- round(nmds.bray.river$stress, 4)
strs4 <- round(nmds.bray.river2$stress, 4)
pval1 <- round(permanova1$aov.tab[1,6], 4)
pval2 <- round(permanova2$aov.tab[1,6], 4)
pval3 <- round(permanova3$aov.tab[1,6], 4)
pval4 <- round(permanova4$aov.tab[1,6], 4)

p1 <- p1 +
  annotate("text", x = 1.3, y = 1.5, hjust = 1, label = sprintf("Stress = %s", strs1)) +
  annotate("text", x = 1.3, y = 1.35, hjust = 1, label = sprintf("italic(P) == %s", pval1), parse = T)
p2 <- p2 +
  annotate("text", x = 0.9, y = 0.35, hjust = 1, label = sprintf("Stress = %s", strs2)) +
  annotate("text", x = 0.9, y = 0.3, hjust = 1, label = sprintf("italic(P) == %s", pval2), parse = T)
p3 <- p3 +
  annotate("text", x = 0.75, y = 0.5, hjust = 1, label = sprintf("Stress = %s", strs3)) +
  annotate("text", x = 0.75, y = 0.45, hjust = 1, label = sprintf("italic(P) == %s", pval3), parse = T)
p4 <- p4 +
  annotate("text", x = 1.0, y = 0.55, hjust = 1, label = sprintf("Stress = %s", strs4)) +
  annotate("text", x = 1.0, y = 0.48, hjust = 1, label = sprintf("italic(P) == %s", pval4), parse = T)

p1 <- p1 + ggtitle("Lake (all ASVs)") + scale_shape_discrete() + geom_point(size = 2.5)
p2 <- p2 + ggtitle("River (all ASVs)") + geom_point(size = 2.5)
p3 <- p3 + ggtitle("Lake (repeatedly-detected ASVs)") + geom_point(size = 2.5)
p4 <- p4 + ggtitle("River (repeatedly-detected ASVs)") + geom_point(size = 2.5)

# Call plot
Fig.NMDS1 <- plot_grid(p1 + theme(legend.position = "none"),
                       p2 + theme(legend.position = "none"),
                       p.legend,
                       rel_widths = c(1,1,0.3),
                       nrow = 1,
                       labels = c("a","b",NULL))

Fig.NMDS2 <- plot_grid(p3 + theme(legend.position = "none"),
                       p4 + theme(legend.position = "none"),
                       p.legend2,
                       rel_widths = c(1,1,0.3),
                       nrow = 1,
                       labels = c("a","b",NULL))

Fig.NMDS.all0 <- plot_grid(p1 + theme(legend.position = "none"),
                       p2 + theme(legend.position = "none"),
                       p3 + theme(legend.position = "none"),
                       p4 + theme(legend.position = "none"),
                       nrow = 2,
                       labels = c("a","b", "c", "d"),
                       align = "hv")
legend.all <- plot_grid(p.legend, p.legend, ncol = 1)
Fig.NMDS.all <- plot_grid(Fig.NMDS.all0, legend.all, ncol = 2, rel_widths = c(1,0.2))

# Save and generate merged figure
save_plot("0_RawFigs/PNG/Fig_NMDS1.png", Fig.NMDS1, base_aspect_ratio = 1.1, ncol = 2, nrow = 1)
save_plot("0_RawFigs/PNG/Fig_NMDS2.png", Fig.NMDS2, base_aspect_ratio = 1.1, ncol = 2, nrow = 1)
save_plot("0_RawFigs/PNG/Fig_NMDS_all.png", Fig.NMDS.all, base_aspect_ratio = 1.1, ncol = 2, nrow = 2)
quartz(type="pdf", file="0_RawFigs/PDF/Fig_NMDS1.pdf", width=9.2, height=4); Fig.NMDS1; dev.off()
quartz(type="pdf", file="0_RawFigs/PDF/Fig_NMDS2.pdf", width=9.2, height=4); Fig.NMDS2; dev.off()
quartz(type="pdf", file="0_RawFigs/PDF/Fig_NMDS_all.pdf", width=9.2, height=8); Fig.NMDS.all; dev.off()
