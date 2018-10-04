####
#### Script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### No. 00 Shell commands to demultiplex fastq files
#### For information on Claident, please see https://www.fifthdimension.jp/documents/metabarcodingtextbook/metabarcodingtextbook.en.html
#### 2018.6.22 Ushio
####
#### bcl2fastq2 v2.18
#### Claident v0.2.2017.05.22
####

# Set parameters
RUN_FOLDER="" # Run folder name
FASTQ_OUT_FOLDER="fastq_out" # FASTQ output folder name
RUNNAME="RMR-032" # Run name
F_PRIMER="" # Foward primer text file name
R_PRIMER="" # Reverse primer text file name
I7_INDEX="" # i7 index primer text file name
I5_INDEX="" # i5 index primer text file name
DEMULTIPLEX_OUT="" # Demultiplex output foler name

# Convert Bcl to Fastq (bcl2fastq2 v2.18)
bcl2fastq --processing-threads 72 --use-bases-mask Y250n,I8,I8,Y250n --create-fastq-for-index-reads --runfolder-dir $RUN_FOLDER --output-dir $FASTQ_OUT_FOLDER

# Demultiplexing using the command implemented in Claident
clsplitseq --runname=$RUNNAME --index1file=$I7_INDEX --index2file=$I5_INDEX --primerfile=$F_PRIMER --reverseprimerfile=$R_PRIMER --minqualtag=30 --numthreads=72 --truncateN=enable *_R1_001.fastq.gz *_I1_001.fastq.gz *_I2_001.fastq.gz *_R2_001.fastq.gz $DEMULTIPLEX_OUT
