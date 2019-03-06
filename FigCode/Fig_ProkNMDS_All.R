####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Figure: NMDS
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
#### 2019.2.20 revised, Ushio
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
(permanova1 <- adonis(distance(ps.lake, method = "bray") ~ data.frame(sample_data(ps.lake))$Method, permutations = 9999))
# = 0.0182
(permanova2 <- adonis(distance(ps.river, method = "bray") ~ data.frame(sample_data(ps.river))$Method, permutations = 9999))
# = 0.0621
(permanova3 <- adonis(distance(ps.pond, method = "bray") ~ data.frame(sample_data(ps.pond))$Method, permutations = 9999))
# < 2e-04
(permanova4 <- adonis(distance(ps.sea, method = "bray") ~ data.frame(sample_data(ps.sea))$Method, permutations = 9999))
# < 1e-04
# ASVs >= 3 data
(permanova5 <- adonis(distance(pl.merge, method = "bray") ~ data.frame(sample_data(pl.merge))$Method, permutations = 9999))
# = 0.0072
(permanova6 <- adonis(distance(pr.merge, method = "bray") ~ data.frame(sample_data(pr.merge))$Method, permutations = 9999))
# < 1e-04
(permanova7 <- adonis(distance(pp.merge, method = "bray") ~ data.frame(sample_data(pp.merge))$Method, permutations = 9999))
# < 1e-04
(permanova8 <- adonis(distance(ps.merge, method = "bray") ~ data.frame(sample_data(ps.merge))$Method, permutations = 9999))
# < 1e-04

# Add P value and stress values (only for repeatedly detected ASVs)
strs1 <- round(nmds.bray.lake2$stress, 4)
strs2 <- round(nmds.bray.river2$stress, 4)
strs3 <- round(nmds.bray.pond2$stress, 4)
strs4 <- round(nmds.bray.sea2$stress, 4)

pval1 <- round(permanova5$aov.tab[1,6], 4)
pval2 <- round(permanova6$aov.tab[1,6], 4)
pval3 <- round(permanova7$aov.tab[1,6], 4)
pval4 <- round(permanova8$aov.tab[1,6], 4)

q1 <- q1 +
  annotate("text", x = 0.9, y = 0.5, hjust = 1, label = sprintf("Stress = %s", strs1)) +
  annotate("text", x = 0.9, y = 0.45, hjust = 1, label = sprintf("italic(P) == %s", pval1), parse = T)
q2 <- q2 +
  annotate("text", x = 1.1, y = 0.52, hjust = 1, label = sprintf("Stress = %s", strs2)) +
  annotate("text", x = 1.1, y = 0.44, hjust = 1, label = sprintf("italic(P) == %s", pval2), parse = T)
q3 <- q3 +
  annotate("text", x = 1.3, y = 0.5, hjust = 1, label = sprintf("Stress = %s", strs3)) +
  annotate("text", x = 1.3, y = 0.45, hjust = 1, label = sprintf("italic(P) == %s", pval3), parse = T)
q4 <- q4 +
  annotate("text", x = 0.5, y = 0.6, hjust = 1, label = sprintf("Stress = %s", strs4)) +
  annotate("text", x = 0.5, y = 0.55, hjust = 1, label = sprintf("italic(P) == %s", pval4), parse = T)

q1 <- q1 + ggtitle("Lake (repeatedly-detected ASVs)") + scale_shape_discrete() + geom_point(size = 2.5)
q2 <- q2 + ggtitle("River (repeatedly-detected ASVs)") + geom_point(size = 2.5)
q3 <- q3 + ggtitle("Pond (repeatedly-detected ASVs)") + geom_point(size = 2.5)
q4 <- q4 + ggtitle("Sea (repeatedly-detected ASVs)") + geom_point(size = 2.5)

# Supplementary figures
r5 <- p5 + theme(legend.position = "none") + ggtitle("All ASVs") + 
  annotate("text", x = -0.2, y = 0.6, label = "Pond") +
  annotate("text", x = -0.7, y = 0.2, label = "River") +
  annotate("text", x = 0.5, y = 0.2, label = "Sea") +
  annotate("text", x = -0.1, y = -0.3, label = "Lake")
s5 <- q5 + theme(legend.position = "none") + ggtitle("Repeatedly detected ASVs") + 
  annotate("text", x = -0.5, y = 2, label = "Pond") +
  annotate("text", x = -2.5, y = 0.8, label = "River") +
  annotate("text", x = 2.5, y = 0.5, label = "Sea") +
  annotate("text", x = -0.5, y = -1, label = "Lake")

# Call plot
Fig.NMDS1 <- plot_grid(q1 + theme(legend.position = "none"),
                       q2 + theme(legend.position = "none"),
                       p.legend,
                       q3 + theme(legend.position = "none"),
                       q4 + theme(legend.position = "none"),
                       p.legend,
                       rel_widths = c(1,1,0.3),
                       nrow = 2,
                       labels = c("a","b","","c","d",""))

Fig.NMDS2 <- plot_grid(r5, s5,
                       p.legend,
                       rel_widths = c(1,1,0.3),
                       nrow = 1,
                       labels = c("a","b",NULL))

# Save and generate merged figure
quartz(type="pdf", file="0_RawFigs/PDF/Fig_NMDS1.pdf", width=9.2, height=8); Fig.NMDS1; dev.off()
quartz(type="pdf", file="0_RawFigs/PDF/Fig_NMDS2.pdf", width=9.2, height=4); Fig.NMDS2; dev.off()

saveRDS(Fig.NMDS1, "0_RawFigs/Robj/Fig_NMDS1.obj")
saveRDS(q1, "0_RawFigs/Robj/Fig_NMDS1_a.obj")
saveRDS(q2, "0_RawFigs/Robj/Fig_NMDS1_b.obj")
saveRDS(q3, "0_RawFigs/Robj/Fig_NMDS1_c.obj")
saveRDS(q4, "0_RawFigs/Robj/Fig_NMDS1_d.obj")
saveRDS(p.legend, "0_RawFigs/Robj/Fig_NMDS_legend.obj")

saveRDS(Fig.NMDS2, "0_RawFigs/Robj/Fig_NMDS2.obj")
saveRDS(r5, "0_RawFigs/Robj/Fig_NMDS2_a.obj")
saveRDS(s5, "0_RawFigs/Robj/Fig_NMDS2_b.obj")
