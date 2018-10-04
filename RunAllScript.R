####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### 2018.6.22 Masayuki Ushio
#### R 3.4.3
####

# Impotant notes
# Sourcing this file will run the analyses and produce the figures for the paper.
# If you want to run all analyses including No. 00 (bcl2fastq and Claident demultiplex), you need raw MiSeq data files. Please contact ong8181@gmail.com.
# If you want to run analyses from No.01, you need fastq files registered in DRA. Please download sequence files (DRA Submission ID:DRA006959).
# Please not that Running scripts No.00, 01_1 and 01_2 may take time.

# Run all scripts
#system("./00_Demultiplex.sh") # Or run the commands in your terminal. You need to specify MiSeq run folder to use bcl2fastq and Claident commands.
#source("01_1_ProSeqDADA2_r32.R") # Running this script may take time depending on your computer.
#system("./01_2_QCauto.sh") # Or run the commands in your terminal.
#system("./02_ident_STD_BLASTn.sh") # Or run the commands in your terminal.
#source("03_TaxaSTDcombine.R")

# Standard sequence check and reads conversion
source("04_STDCheck.R")

# Statistical analyses using "phyloseq" package
source("05_PhyloseqPlot.R")
source("06_1_RareProkAnalysis.R")
source("06_2_RareProkAnalysis.R")

# Figure generation
setwd("FigsCode")
source("Fig_ProkDivAbn.R")
source("Fig_ProkMicsBar_All_Rare.R")
source("Fig_ProkNMDS_BeadsTest.R")
source("Fig_STDcheck_BeadsTest.R")
source("FigS_Beads.R")
source("FigS_FieldNanodrop.R")
source("FigS_ProkRarefaction.R")
source("FigS_ProkRareMics_v2.R")
source("Tab_DADA2Summary.R")
source("Tab_DRAinfo.R")
source("Tab_ProkRareMics.R")
source("Tab_RareMIcs.R")
