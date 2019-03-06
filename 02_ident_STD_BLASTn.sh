####
#### Script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No.2 Standard sequence identification using BLASTN
#### 2018.6.22 Ushio
#### 2019.2.20 revised, Ushio
####

# Identify Standard DNA of the 1st MiSeq run
DBPATH=~/DADA2_DB/STDseqs/ProkaryoteSTD/ProkSTD_515F # Path to the database
QUERYPATH=01_ProSeqDADA2_Out/ProkASV_seqs.fa # DADA2 ASV output
OUTPUT=02_ident_STD_BLASTnOut/ProkSTD_out.txt
EVALUE_SET=1e-50

mkdir 02_ident_STD_BLASTnOut
blastn -db $DBPATH -query $QUERYPATH -evalue $EVALUE_SET -outfmt 6 -out $OUTPUT
