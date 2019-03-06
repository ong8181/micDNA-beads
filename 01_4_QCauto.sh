####
#### Script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### No. 01_4 Shell commands to identify taxa
#### For information on Claident, please see https://www.fifthdimension.jp/documents/metabarcodingtextbook/metabarcodingtextbook.en.html
#### 2018.6.22 Ushio
#### 2019.2.20 revised, Ushio
####
#### Claident v0.2.2018.05.29
####

# Taxa assignment using claident (after picking ASV by DADA2)
# Commands to identify ASV taxa using QCauto method implemented in Claident
clmakecachedb --blastdb=semiall_genus --numthreads=72 ProkASV_seqs.fa ProkASV_semiall_cache
clidentseq --blastdb=ProkASV_semiall_cache --numthreads=72 ProkASV_seqs.fa ProkASV_semiall_clidentseq
classigntax --taxdb=semiall_genus ProkASV_semiall_clidentseq ProkASV_semiall_classigntax
