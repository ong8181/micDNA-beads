####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### No.4 Check standard sequences and reads conversion
#### 2018.6.22 Ushio
#### R 3.4.3
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("04_STDCheckOut")

# Load library and functions
# For package versions used to analyze data,
# please see "00_SessionInfo_original"
library(cowplot); packageVersion("cowplot")
library(reshape2); packageVersion("reshape2")
library(tidyverse); packageVersion("tidyverse")
library(ggsci); packageVersion("ggsci")

# Load helper functions
source("functions/F1_HelperFunctions.R")

# Load workspace
load("01_ProSeqDADA2_r32Out/01_ProSeqDADA2_r32Out.RData")
std.copy.n <- c(26900, 12300, 6800, 3200, 540)
col.name.level <- "Family"
tax.claident <- read.csv("03_TaxaSTDcombineOut/claident.tax.revise.csv")
sample.sheet <- read.csv("sampledata/RMR-032_sampleP.csv")

# DADA2 summary objects
track

# Extract standard sequeces
detected.std.name <- unique(tax.claident[which(substr(tax.claident[,"species"],1,8) == "STD_prok"),"species"])
n.std.seq <- which(substr(tax.claident[,"species"],1,8) == "STD_prok")
std.table <- seqtab.nochim[,n.std.seq]
std.taxa <- tax.claident[n.std.seq, "species"]

# STD reads - copy number relationship
# Rename colnames
colnames(std.table) <- std.taxa
# Merge the same colnames
new.std.table <- data.frame(std_rank1 = MergeSTD(detected.std.name[1], std.data = std.table),
                            std_rank2 = MergeSTD(detected.std.name[2], std.data = std.table),
                            std_rank3 = MergeSTD(detected.std.name[3], std.data = std.table),
                            std_rank4 = MergeSTD(detected.std.name[4], std.data = std.table),
                            std_rank5 = MergeSTD(detected.std.name[5], std.data = std.table))
new.std.table2 <- new.std.table[sample.sheet$sample_nc != "pcr_nc",]

# Linear regression
adj.r.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$adj.r.squared
lm.coef.fun <- function(x) summary(lm(as.numeric(x) ~ std.copy.n + 0))$coefficients[1]
r2.summary <- apply(new.std.table2, 1, adj.r.fun)
coef.summary <- apply(new.std.table2, 1, lm.coef.fun)
new.seqtab <- as.data.frame(seqtab.nochim[,-n.std.seq]) # Make seq table without standard DNA
new.seqtab2 <- new.seqtab[sample.sheet$sample_nc != "pcr_nc",] # Make seq table without negative control samples

## Visualize regression results
# 1. R2 value distribution
g1 <- ggplot(data.frame(values = r2.summary), aes(x = values))
g1 <- g1 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2) + xlim(0.9, 1)
g1 <- g1 + xlab(expression(paste("R"^{2}, " values"))) + ylab("Count")
g1 <- PlotStyle(g1)

# 2. Slope distribution
g2 <- ggplot(data.frame(values = coef.summary), aes(x = values))
g2 <- g2 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2)
g2 <- g2 + xlab("Regression slope") + ylab("Count") + xlim(0, 1.5)
g2 <- PlotStyle(g2)

# 3. Regression examples
max.slope <- as.numeric(c(new.std.table2[which.max(coef.summary),], 0))
med.slope <- as.numeric(c(new.std.table2[which.min(abs(coef.summary - median(coef.summary))),], 0))
min.slope <- as.numeric(c(new.std.table2[which.min(coef.summary),], 0))
slope.summary <- melt(data.frame(copy =c(std.copy.n, 0),
                                 max_slope = max.slope,
                                 med_slope = med.slope,
                                 min_slope = min.slope), id.vars = "copy")
g3 <- ggplot(slope.summary, aes(x = copy, y = value, group = variable, colour = variable))
g3 <- g3 + geom_point(size = 2) + scale_color_d3(name = "Regression slope")
g3 <- g3 + geom_smooth(method = "lm", size = 0.5, se = F)
g3 <- g3 + xlab(expression(paste("Standard DNA copies (", mu, l^{-1}, ")"))) + ylab("Standard DNA reads")
g3 <- PlotStyle(g3) + theme(legend.position = c(0.25, 0.75), legend.text = element_text(size = 7), legend.title = element_text(size = 8))

# 4. Read visualization
new.seqtab2$sample <- sample.sheet$Sample_Name2[sample.sheet$sample_nc != "pcr_nc"]
seqtab.plot <- melt(new.seqtab2, id.vars = "sample") %>% filter(value > 0)
seqtab.plot$value2 <- seqtab.plot$value + 0.5
std.table.plot <- data.frame(sample = sample.sheet$Sample_Name2[sample.sheet$sample_nc != "pcr_nc"],
                             max = apply(new.std.table2, 1, max) + 0.5,
                             min = apply(new.std.table2, 1, min) + 0.5) # data.frame for overwrite standard DNA reads

# Box plot + jitter plot version
g4 <- ggplot(seqtab.plot, aes(x = sample, y = value2))
g4 <- g4 + geom_boxplot(shape = 16, alpha = 0.5)
g4 <- g4 + geom_jitter(shape = 16, size = 1.5, alpha = 0.5, position = position_jitter(0.2))
g4 <- g4 + scale_y_log10() + ylab("Sequence reads") + xlab("Sample ID")
g4 <- PlotStyle(g4) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Overwrite standard DNA reads
g4 <- g4 + geom_segment(data = std.table.plot, aes(x = sample, xend = sample, y = min, yend = max),
                        colour = "red3", size = 1.5, alpha = 0.5)
g4 <- g4 + ggtitle("RMR-032, DNA extraction prokaryote seqs")

# Summarize visualization
top.row <- plot_grid(g1, g2, g3, ncol = 3, align = "h", labels = c("a", "b", "c"), rel_widths = c(1,1,1.3))
Fig.std <- plot_grid(top.row, g4, ncol = 1, labels = c("", "d"))
OpenDev(8, 5.5)
Fig.std
dev.off()

# Raw data check
new.seqtab.print <- seqtab.nochim
colnames(new.seqtab.print) <- sprintf("Taxa%05d", 1:ncol(seqtab.nochim))
tax.claident2 <- as.data.frame(tax.claident)
rownames(tax.claident2) <- sprintf("Taxa%05d", 1:ncol(seqtab.nochim))

# Conversion of sample reads to calculated copy numbers
# Collect valid samples
valid.samples <- coef.summary[!(coef.summary == 0 | r2.summary < 0.5)]
seqtab.valid <- seqtab.nochim[match(names(valid.samples), rownames(seqtab.nochim)),]
seqtab.conv.w.std <- seqtab.valid/valid.samples # Conversion

taxa.w.std <- tax.claident2
dim(seqtab.conv.w.std)
dim(taxa.w.std)

# Exclude standard DNAs from converted tables
seqtab.conv <- seqtab.conv.w.std[,-which(substr(tax.claident[,"species"],1,8) == "STD_prok")]
taxa.wo.std <- tax.claident2[-which(substr(tax.claident[,"species"],1,8) == "STD_prok"),]
dim(seqtab.conv)
dim(taxa.wo.std)

# Save and output results
save.image("04_STDCheckOut/04_STDCheckOut.RData")

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/04_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))