####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No.4 Check standard sequences and reads conversion
#### 2018.6.22 Ushio
#### 2019.2.20 revised, Ushio
#### R 3.5.2
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
load("01_ProSeqDADA2_Out/01_ProSeqDADA2_Merged_Out.RData")
std.copy.n1 <- c(26900, 12300, 6800, 3200, 540) # For 1st run
std.copy.n2 <- c(40000, 20000, 10000, 4000, 2000) # For 2nd run
col.name.level <- "Family"
tax.claident <- read.csv("03_TaxaSTDcombineOut/claident_tax_revise.csv")
sample.sheet <- read.csv("sampledata/20190220_SampleData.csv")

# Adjust factor levels
sample.sheet$Sample_Name2 <- factor(sample.sheet$Sample_Name2, levels = c(sprintf("S%03d", 1:40), sprintf("R%03d", 1:40)))
sample.sheet$Site <- factor(sample.sheet$Site, levels = c("lake", "river", "pond", "sea"))

# DADA2 summary objects
track <- read.csv("01_ProSeqDADA2_Out/track.csv")

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
adj.r.fun1 <- function(x) summary(lm(as.numeric(x) ~ std.copy.n1 + 0))$adj.r.squared
adj.r.fun2 <- function(x) summary(lm(as.numeric(x) ~ std.copy.n2 + 0))$adj.r.squared
lm.coef.fun1 <- function(x) summary(lm(as.numeric(x) ~ std.copy.n1 + 0))$coefficients[1]
lm.coef.fun2 <- function(x) summary(lm(as.numeric(x) ~ std.copy.n2 + 0))$coefficients[1]
r2.summary <- c(apply(new.std.table2[1:38,], 1, adj.r.fun1), apply(new.std.table2[39:76,], 1, adj.r.fun2))
coef.summary <- c(apply(new.std.table2[1:38,], 1, lm.coef.fun1), apply(new.std.table2[39:76,], 1, lm.coef.fun2))
new.seqtab <- as.data.frame(seqtab.nochim[,-n.std.seq]) # Make seq table without standard DNA
new.seqtab2 <- new.seqtab[sample.sheet$sample_nc != "pcr_nc",] # Make seq table without negative control samples

## Visualize regression results
# 1. R2 value distribution
g1 <- ggplot(data.frame(values = r2.summary, run = c(rep("1st", 38), rep("2nd", 38))),
             aes(x = values, fill = run))
g1 <- g1 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2) + xlim(0.8, 1)
g1 <- g1 + xlab(expression(paste("R"^{2}, " values"))) + ylab("Count")
g1 <- PlotStyle(g1) + scale_fill_d3()

# 2. Slope distribution
g2 <- ggplot(data.frame(values = coef.summary, run = c(rep("1st", 38), rep("2nd", 38))),
             aes(x = values, fill = run))
g2 <- g2 + geom_histogram() + geom_hline(yintercept = 0, linetype = 2)
g2 <- g2 + xlab("Regression slope") + ylab("Count") + xlim(-0.05, 1.3)
g2 <- PlotStyle(g2) + scale_fill_d3()

# 3. Regression examples
max.slope1 <- as.numeric(c(new.std.table2[1:38,][which.max(coef.summary[1:38]),], 0))
med.slope1 <- as.numeric(c(new.std.table2[1:38,][which.min(abs(coef.summary[1:38] - median(coef.summary[1:38]))),], 0))
min.slope1 <- as.numeric(c(new.std.table2[1:38,][which.min(coef.summary[1:38]),], 0))
slope.summary1 <- melt(data.frame(copy =c(std.copy.n1, 0),
                                 max_slope = max.slope1,
                                 med_slope = med.slope1,
                                 min_slope = min.slope1), id.vars = "copy")
max.slope2 <- as.numeric(c(new.std.table2[39:76,][which.max(coef.summary[39:76]),], 0))
med.slope2 <- as.numeric(c(new.std.table2[39:76,][which.min(abs(coef.summary[39:76] - median(coef.summary[39:76]))),], 0))
min.slope2 <- as.numeric(c(new.std.table2[39:76,][which.min(coef.summary[39:76]),], 0))
slope.summary2 <- melt(data.frame(copy =c(std.copy.n2, 0),
                                  max_slope = max.slope2,
                                  med_slope = med.slope2,
                                  min_slope = min.slope2), id.vars = "copy")
slope.summary <- cbind(rbind(slope.summary1, slope.summary2), c(rep("1st", 18), rep("2nd", 18)))
colnames(slope.summary)[4] <- "run"
g3 <- ggplot(slope.summary, aes(x = copy, y = value, group = variable, colour = variable, facet = run))
g3 <- g3 + geom_point(size = 2) + scale_color_d3(name = "Regression slope")
g3 <- g3 + geom_smooth(method = "lm", size = 0.5, se = F) + facet_grid(.~run)
g3 <- g3 + xlab(expression(paste("Standard DNA copies (", mu, l^{-1}, ")"))) + ylab("Standard DNA reads")
g3 <- PlotStyle(g3) #+ theme(legend.position = c(0.25, 0.75), legend.text = element_text(size = 7), legend.title = element_text(size = 8))

# 4. Read visualization
new.seqtab2$sample <- sample.sheet$Sample_Name2[sample.sheet$sample_nc != "pcr_nc"]
seqtab.plot <- melt(new.seqtab2, id.vars = "sample") %>% filter(value > 0)
seqtab.plot$value2 <- seqtab.plot$value + 0.5
seqtab.plot$sample <- factor(seqtab.plot$sample, levels = c(sprintf("S%03d", 1:40), sprintf("R%03d", 1:40)))
std.table.plot <- data.frame(sample = sample.sheet$Sample_Name2[sample.sheet$sample_nc != "pcr_nc"],
                             max = apply(new.std.table2, 1, max) + 0.5,
                             min = apply(new.std.table2, 1, min) + 0.5) # data.frame for overwrite standard DNA reads
std.table.plot$sample <- factor(std.table.plot$sample, levels = c(sprintf("S%03d", 1:40), sprintf("R%03d", 1:40)))

# Box plot + jitter plot version
g4 <- ggplot(seqtab.plot, aes(x = sample, y = value2))
g4 <- g4 + geom_boxplot(shape = 16, alpha = 0.5)
g4 <- g4 + geom_jitter(shape = 16, size = 1.5, alpha = 0.5, position = position_jitter(0.2))
g4 <- g4 + scale_y_log10() + ylab("Sequence reads") + xlab("Sample ID")
g4 <- PlotStyle(g4) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Overwrite standard DNA reads
g4 <- g4 + geom_segment(data = std.table.plot, aes(x = sample, xend = sample, y = min, yend = max),
                        colour = "red3", size = 1.5, alpha = 0.5)
g4 <- g4 + ggtitle("DNA extraction prokaryote seqs")

# Summarize visualization
top.row <- plot_grid(g1 + theme(legend.position = "none"),
                     g2, g3,
                     ncol = 3, labels = c("a", "b", "c"), rel_widths = c(0.48,0.6,1.2))
Fig.std <- plot_grid(top.row, g4, ncol = 1, labels = c("", "d"))
ggsave("04_STDCheckOut/DNAextr_STD_summary.pdf", Fig.std, width = 14, height = 8)

# Raw data check
new.seqtab.print <- seqtab.nochim
colnames(new.seqtab.print) <- sprintf("Taxa%05d", 1:ncol(seqtab.nochim))
tax.claident2 <- as.data.frame(tax.claident)
rownames(tax.claident2) <- sprintf("Taxa%05d", 1:ncol(seqtab.nochim))

# Conversion of sample reads to calculated copy numbers
# Collect valid samples (all samples are valid!)
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
save.image("04_STDCheckOut/04_STDCheck_Out.RData")

# Save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/04_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))