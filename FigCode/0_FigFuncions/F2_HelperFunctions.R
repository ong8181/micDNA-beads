####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with  intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### F2. Color assigner for microbial taxa
#### 2019.2.21 Ushio
####

#### ColorAssigner function
ColorAssigner <- function(input_list){
  taxa.color.list <- matrix(c(
    "Acidobacteria", "#CE3D32", # From scale_fill_igv() in "ggsci" package
    "Actinobacteria", "#FFC20A", #"#749B58",
    "Bacteroidetes", "#1A0099", #"#F0E685",
    "Chloroflexi", "#0099CC",
    "Crenarchaeota","#FF0000", #"#466983",
    
    "Cyanobacteria", "#99CC00",
    "Euryarchaeota", "#996600",
    "Firmicutes", "#EE4C97", #"#802268", #"#D58F5C"
    "Fusobacteria", "#A6EEE6", #"#990080", #
    "Gemmatimonadetes", "#466983", #"#6BD76B",
    
    "Ignavibacteriae", "#00468B", #"#4775FF", #"#D595A7",
    "Nitrospirae", "#FF7F0E", #"#837B8D",
    "Planctomycetes", "#F0E685",
    "Proteobacteria", "#4775FF",
    "Thaumarchaeota", "#E4AF69",
    
    "Verrucomicrobia",  "#749B58", #"#C75127",
    "Others", "gray20", #"#749B58",
    "Undetermined", "#ADB6B6",
    "Deinococcus-Thermus", "#F0F8FF"
  ),
  ncol = 2, byrow = T
  )
  
  return(taxa.color.list[match(input_list, taxa.color.list[,1]),2])
}


# From https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

