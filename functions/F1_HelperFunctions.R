####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### F1. Collection of helper functions for DNA extration study
#### 2018.6.22 Ushio
####

#### ggplot function
PlotStyle <-  function(ggobject){
  return(ggobject + theme_bw() + theme(axis.text.x = element_text(angle=0),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.text = element_text(size=12),
                                       axis.title = element_text(size=12),
                                       panel.background=element_rect(colour="black", fill=NA, size=0.8)))
}

# Merge standard DNA sequences
MergeSTD <- function(std.i, std.data = std.table){
  index.std <- which(match(colnames(std.table), std.i) == 1)
  if(length(index.std) > 1){
    std.tmp <- rowSums(std.table[,index.std])
  }else{
    std.tmp <- std.table[,index.std]
  }
  return(std.tmp)
}

# Phylum sorting
PhylumSort <- function(phyloseq.object, select.n = 51){
  phylum.sort <- names(sort(taxa_sums(phyloseq.object), decreasing=TRUE))
  phylum.n <- length(phylum.sort)
  merged.phyloseq.object <- merge_taxa(phyloseq.object, phylum.sort[select.n:phylum.n])
  merged.phyloseq.object <- subset_taxa(merged.phyloseq.object, Phylum != "Not_Identified")
  return(merged.phyloseq.object)
}

# Make method-specific phyloseq object
ExtractMethodTaxa <- function(phyloseq.object){
  taxa.p.m1 <- names(sort(taxa_sums(phyloseq.object)[taxa_sums(phyloseq.object) > 0], decreasing=TRUE))
  return(taxa.p.m1)
}

ExtractSpecificTaxa <- function(comp.1, comp.2, comp.3){
  specific.taxa <- comp.1[is.na(match(comp.1, unique(c(comp.2, comp.3))))]
  return(specific.taxa)
}

# Identify OS and open figure device
OpenDev <- function(width = width, height = height){
  if(Sys.info()['sysname'] == "Linux"){
    x11(width = width, height = height)
  }else if(Sys.info()['sysname'] == "Darwin"){
    quartz(width = width, height = height)
  }else if(Sys.info()['sysname'] == "Windows"){
    windows(width = width, height = height)
  }
}
