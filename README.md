# Shell and R scripts for Ushio (2019) _Methods in Ecology and Evolution_
The scripts are to reproduce results and figures in Ushio (2019) _Methods in Ecology and Evolution_ Citation information is as follows:

Published version: Ushio (2019) "Use of a filter cartridge combined with intra-cartridge bead-beating improves detection of microbial DNA from water samples" _Methods in Ecology and Evolution_ https://doi.org/10.1111/2041-210X.13204

bioRxiv version: Ushio (2018) "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples" _bioRxiv_ 10.1101/435305 https://doi.org/10.1101/435305

# Liscence
This code is released under the MIT License, see LICENSE.

# Impotant notes
If you want to run all analyses including No. 00 (bcl2fastq and Claident demultiplex), you need raw MiSeq data files. Please contact ong8181@gmail.com. If you want to run the analyses from No.01, you need fastq files registered in DRA. Please download the sequence files (see sequence data availability below). Please note that Running scripts No.00, 01_1 and 01_2 may take time.

# Software and package versions
### 1. For the sequence processing:
Claident 0.2.2018.05.29, Blast 2.8.1, R 3.5.2, Rcpp 1.0.0, dada2 1.10.1, ShortRead 1.40.1, tidyverse 1.2.1

### 2. For the statistical analyses:
phyloseq 1.26.1

### 3. For the visualization:
ggplot2 2.3.1, cowplot 0.9.4, reshape2 1.4.3, ggsci 2.9

For details, please see "00_SessionInfo_original" folder.

# Installations of core softwares
Claident: https://www.claident.org/

R: https://www.r-project.org/

DADA2: https://benjjneb.github.io/dada2/

phyloseq: https://joey711.github.io/phyloseq/

# Sequence data availability
Sequence data associated with this scripts are available in DDBJ Sequence Read Archive (DRA).

Accession numbers are as follows:

BioProject ID = PRJDB7110

BioSample ID = SAMD00127558-SAMD00127597 (1st run) and SAMD00163348-SAMD00163387 (2nd run)

Direct links to DNA sequences:

1st MiSeq run, https://ddbj.nig.ac.jp/DRASearch/submission?acc=DRA006959

2nd MiSeq run, https://ddbj.nig.ac.jp/DRASearch/submission?acc=DRA008110
