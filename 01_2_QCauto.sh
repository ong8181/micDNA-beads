####
#### Script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### No. 01_1 Shell commands to indentify taxa
#### For information on Claident, please see https://www.fifthdimension.jp/documents/metabarcodingtextbook/metabarcodingtextbook.en.html
#### 2018.6.22 Ushio
####
#### Claident v0.2.2017.05.22
####

# Taxa assignment using claident (after picking ASV by DADA2)
# Commands to identify ASV taxa using QCauto method implemented in Claident
clidentseq --blastdb=prokaryota_16S_genus --numthreads=72 01_ProSeqDADA2_r32Out/ProkASV_seqs.fa 01_ProSeqDADA2_r32Out/ProkASV_seqs_clidentseq
classigntax --taxdb=prokaryota_16S_genus 01_ProSeqDADA2_r32Out/ProkASV_seqs_clidentseq 01_ProSeqDADA2_r32Out/ProkASV_seqs_classigntax

