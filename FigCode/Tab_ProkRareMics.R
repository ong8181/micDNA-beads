####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### Table: A list of method-specific microbial ASVs
#### 2018.6.22 Ushio (Run on Mac OSX)
#### 2019.2.21 revised Ushio
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
p1 <- tax_table(p.wb.ps)
p2 <- tax_table(p.wo.ps)
p3 <- tax_table(p.ds.ps)
s1 <- tax_table(s.wb.ps)
s2 <- tax_table(s.wo.ps)
s3 <- tax_table(s.ds.ps)

all.list <- rbind(l1, NA, l2, NA, l3, NA, NA,
                  r1, NA, r2, NA, r3, NA, NA,
                  p1, NA, p2, NA, p3, NA, NA,
                  s1, NA, s2, NA, s3)
tax.id.list <- rownames(all.list)

# Output 
write.csv(all.list, "0_Table/CSV/MethodSpecificMics.csv")
write.table(tax.id.list, "0_Table/CSV/MethodSpecificTaxIDs", row.names = F, col.names = F)
