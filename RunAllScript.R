####
#### R script for Ushio (2019)
#### "Use of a filter cartridge combined with intra-cartridge bead-beating improves detection of microbial DNA from water samples"
#### 2018.6.22 Masayuki Ushio
#### 2019.2.27 revised, Masayuki Ushio
#### 2019.5.24 revised, Masayuki Ushio
#### R 3.5.2
####

# Impotant notes
# Sourcing this file will run the analyses and produce the figures for the paper.
# If you want to run all analyses including No. 00 (bcl2fastq and Claident demultiplex), you need raw MiSeq data files. Please contact ong8181@gmail.com.
# If you want to run analyses from No.01, you need fastq files registered in DRA. Please download sequence files (DRA BioProject ID: PRJDB7110).
# Please not that Running scripts No.00, 01_1 and 01_2 may take time.

# Run all scripts
#system("./00_Demultiplex.sh") # Or run the commands in your terminal. You need to specify MiSeq run folder to use bcl2fastq and Claident commands.
#source("01_1_ProSeqDADA2_1st.R") # Running this script may take time depending on your computer.
#source("01_2_ProSeqDADA2_2nd.R") # Running this script may take time depending on your computer.
#source("01_3_ProSeqDADA2_merge.R")
#01_4_QCauto.sh # Run the commands in your terminal.
#02_ident_STD_BLASTn.sh # Run the commands in your terminal.
#source("03_TaxaSTDcombine.R")

# Standard sequence check and reads conversion
source("04_STDCheck.R")

# Statistical analyses using "phyloseq" package
source("05_PhyloseqPlot.R")
source("06_1_RareProkAnalysis.R")
source("06_2_RareProkAnalysis.R")
source("06_3_FieldNCEvaluation.R")

# Figure generation
setwd("FigsCode")
source("Fig_ProkDivAbn.R")
source("Fig_ProkMicsBar_All_Rare.R")
source("Fig_ProkNMDS_BeadsTest.R")
source("Fig_STDcheck_BeadsTest.R")
source("FigS_Beads.R")
source("FigS_FieldNanodrop.R")
source("FigS_FieldNCBar.R")
source("FigS_ProkRarefaction.R")
source("FigS_ProkRareMics_v2.R")
source("Tab_DADA2Summary.R")
source("Tab_ProkRareMics.R")
source("Tab_RareMIcs.R")

# Re-formatting figures to adjust appearance
source("Reformat_AllFigs.R")
