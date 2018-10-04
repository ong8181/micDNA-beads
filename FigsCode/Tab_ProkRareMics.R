####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with bead beating improves detection of microbial DNA from water samples"
#### Table: A list of method-specific microbial ASVs
#### 2018.6.22 Ushio (Run on Mac OSX)
####

# Load library
library(phyloseq); packageVersion("phyloseq")

# Load workspace (Prokaryote data)
load("../06_1_RareProkAnalysisOut/06_1_RareProkAnalysisOut.RData")

# Call method-specific microbe objects
l1 <- tax_table(l.wb.ps)
l2 <- tax_table(l.wo.ps)
l3 <- tax_table(l.ds.ps)
r1 <- tax_table(r.wb.ps)
r2 <- tax_table(r.wo.ps)
r3 <- tax_table(r.ds.ps)

all.list <- rbind(l1, l2, l3, r1, r2, r3)
tax.id.list <- rownames(all.list)

# Output 
write.csv(all.list, "0_Table/CSV/MethodSpecificMics.csv")
write.table(tax.id.list, "0_Table/CSV/MethodSpecificTaxIDs", row.names = F, col.names = F)
