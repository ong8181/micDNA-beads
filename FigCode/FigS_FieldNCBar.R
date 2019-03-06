####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Figure: Bar plots for field NC samples
#### 2018.6.22 Ushio (Figs generated on Mac OSX)
#### 2019.2.25 revised, Ushio
####

# Load library and functions
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")
library("scales"); packageVersion("scales")
library(phyloseq); packageVersion("phyloseq")
par(family="Arial")

# Load workspace (Prokaryote data)
load("../06_3_FieldNCEvaluationOut/06_3_FieldNCEvaluationOut.RData")
# Load my palette
source("0_FigFuncions/F2_HelperFunctions.R")

# Barplot for field negative controls
(f1.levels <- levels(f1$data$phylum))
f1$data$phylum <- factor(f1$data$phylum, levels = c(f1.levels[-match(c("Others", "Undetermined"), f1.levels)], c("Others", "Undetermined")))
f1 <- f1 + scale_fill_manual(values = ColorAssigner(levels(f1$data$phylum))) + scale_y_continuous(labels = scientific_format())
f1 <- f1 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                              axis.title.x = element_blank(),
                              axis.text = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 1))
f1 <- f1 + scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})

# Calculate the proportion of NC DNA copies
mean_copy <- aggregate(sample_sums(ps.phylum),
                       by=list(sample_data(ps.phylum)$Site:sample_data(ps.phylum)$Method),
                       mean)
ps.fnc.all.prop <- ps.fnc.all
otu_table(ps.fnc.all.prop) <- otu_table(ps.fnc.all)/mean_copy[,2]

f2 <- plot_bar(ps.fnc.all.prop, x = "Method", fill = "phylum") +
  ylab("DNA % (NC copy numbers divided by\nmean copy numbers per positive sample)") + facet_grid( ~ Site, scales = "free")
(f2.levels <- levels(f2$data$phylum))
f2$data$phylum <- factor(f2$data$phylum, levels = c(f2.levels[-match(c("Others", "Undetermined"), f2.levels)], c("Others", "Undetermined")))
f2 <- f2 + scale_fill_manual(values = ColorAssigner(levels(f1$data$phylum)))
f2 <- f2 + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                              axis.title.x = element_blank(),
                              axis.text = element_text(colour = "black"),
                              axis.title = element_text(colour = "black"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 1))
f2 <- f2 + ylim(0, 0.1)

# Merge panelsf
Fig.NCBar <- f1
Fig.NCBarProp <- f2

# Save and generate merged figure
quartz(type="pdf", file="0_RawFigs/PDF/FigS_NCBar.pdf", width=7, height=5); Fig.NCBar; dev.off()
quartz(type="pdf", file="0_RawFigs/PDF/FigS_NCBarProp.pdf", width=7, height=5); Fig.NCBarProp; dev.off()
saveRDS(Fig.NCBar, "0_RawFigs/Robj/FigS_NCBar.obj")
saveRDS(Fig.NCBarProp, "0_RawFigs/Robj/FigS_NCBarProp.obj")
