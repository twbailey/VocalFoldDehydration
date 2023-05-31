qPCR.OUTPUT <- function(DF, target, housekeeping="HPRT1", g1="CONT", g2="EXP", ...) {
  library(tidyverse)
  
  #Subset data into gene target groups
  tar = subset(DF, DF$Target== target)
  hou = subset(DF, DF$Target== housekeeping)
  
  #Calculate individual sample means, maintains a "Group" column.  
  tar.stat <- aggregate(tar$CT, list(Sample = tar$Sample, Group = tar$Group), mean, na.rm=TRUE)
  hou.stat <- aggregate(hou$CT, list(Sample = hou$Sample, Group = hou$Group), mean, na.rm=TRUE)
  
  #Calculate individual delCT: sample.target - sample.housekeeping
  tarDEL = tar.stat$x - hou.stat$x
  tar.stat <- mutate(tar.stat, tarDEL)
  
  #Calculate deldelCT: tarDEL - mean.of.control.tarDELS
  tar.stat <- mutate(tar.stat, tarDELDEL = tarDEL - mean(tar.stat$tarDEL[tar.stat$Group==g1], na.rm=TRUE))
  
  #Calculate aggregate mean: target(mean.CT.group1, mean.CT.group2); housekeeping(...)
  agg.mean <- tapply(tar$CT, tar$Group, mean, na.rm=TRUE) - tapply(hou$CT, hou$Group, mean, na.rm=TRUE)
  
  #Calculate figure mean: agg.mean - agg.mean.control
  fig.mean <- agg.mean - agg.mean[1]
  return(list(tar.stat, fig.mean))
}




