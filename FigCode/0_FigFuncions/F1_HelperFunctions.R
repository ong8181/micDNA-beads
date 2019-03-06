####
#### R script for Ushio (2018)
#### "Use of a filter cartridge combined with intra-cartridge bead beating improves detection of microbial DNA from water samples"
#### F1. Collection of helper functions for DNA extration study
#### 2018.6.22 Ushio
#### 2019.2.20 Ushio
####

#### ggplot function
PlotStyle4Pub <-  function(ggobject){
  return(ggobject + theme_bw() + theme(panel.border = element_rect(color="black", fill=NA),
                                       axis.title.x = element_blank(),
                                       axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                                       axis.title = element_text(colour = "black"),
                                       panel.grid.minor = element_blank(),
                                       panel.grid.major = element_blank()))
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
