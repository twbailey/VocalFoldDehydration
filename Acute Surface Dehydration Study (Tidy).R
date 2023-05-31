# #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Code supporting: 
# Bailey, T. W., Dos Santos, A. P., do Nascimento, N. C., Xie, S., 
# Thimmapuram, J., Sivasankar, M. P., & Cox, A. (2020). RNA sequencing
# identifies transcriptional changes in the rabbit larynx in response to low
# humidity challenge. BMC Genomics, 21(1), 888–888.
# https://doi.org/10.1186/s12864-020-07301-7

# Taylor W. Bailey, PhD, MS, MPH
# Updated 19 Sept 2022
# #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Setup #########################################################################
library(readxl)
library(outliers)
library(tidyverse)
library(plotrix)

#Functions #####################################################################
delCTS <- function(df1,df2) {
### Calculates delta Ct from qPCR output
### Expects dataframe with columns: Sample, Target, Ct, Plate
  i<-1
  j<-1
  k<-1
  DCT <- vector()
  for(i in 1:nrow(df1)) {
    for( j in 1:nrow(df2)){
      if(df1$Sample[i] == df2$Sample[j] & df1$Plate[i] == df2$Plate[j]) 
      {DCT[k] <- df1$x[i]-df2$x[j]
      k<-k+1}
      
      j<- j+1}
    i<- i+1}
  return(DCT)
}

deldelCTS <- function(df1,df2) {
  i<-1
  j<-1
  k<-1
  DDCT <- vector()
  
  for(i in 1:nrow(df1)) {
    for( j in 1:nrow(df2)){
      
      if(df1$Target[i] == labels(df2)[[1]][j]) 
      {DDCT[k] <- df1$DCTMODS[i]-df2[j]  ###SWAP OUT "LOWS" AND "MODS" IN DF1
      k<-k+1}
      j<- j+1}
    i<- i+1}
  return(DDCT)
  print(k)
}

#Import QPCR data and formatting ###############################################
file="F:/Programming Things/Acute Surface Dehydration/Data_OutliersRemoved.xlsx"
REDO <- read_excel(file)

#Rename CACC as ECCP
for(i in 1:nrow(REDO)) {if (REDO$Target[i] == "CACC") {REDO$Target[i] <-"ECCP"}}
#SetColumnTypes
REDO$Sample <- as.factor(REDO$Sample)
REDO$Plate <- as.factor(REDO$Plate)
#Ignore valeus grater than 36 or "NA"
REDO <- subset(REDO, REDO$Ct < 36 & is.na(REDO$Ct)==FALSE)

# Initial data subsetting ######################################################
#Calculate Average for each Sample by Target and add group name
#HPRT1 is matched to plate
MEANS <- with(
  REDOsim, aggregate(Ct, list(Sample=Sample, Target=Target, Plate=Plate), mean))
MEANS <- mutate(MEANS, Group = ifelse(Sample %in% 34:39, "Moderate", "Low"))
LOWS <- subset(MEANS, Group == "Low")
MODS <- subset(MEANS, Group == "Moderate")
HPLOWS <- filter(LOWS, Target == "HPRT1")
HPMODS <- filter(MODS, Target== "HPRT1")

# Calculate Rq/FC ##############################################################
#Delta Cts (function above)
# Mean subtractions are matched by both animal and plate (genes were split 
#between plates)
# Individual's HPRT1 is subtracted from the individual's target genes
DCTLOWS <- delCTS(LOWS, HPLOWS)
DCTMODS <- delCTS(MODS, HPMODS)
DLOWS <- cbind(LOWS, DCTLOWS)
DMODS <- cbind(MODS, DCTMODS)
LOWAVGS <- tapply(DLOWS$DCTLOWS, DLOWS$Target, mean)
CONAVGS <- tapply(DMODS$DCTMODS, DMODS$Target, mean)

#Delta delta Cts (function above)
DDCTLOWS <- deldelCTS(DLOWS, CONAVGS)
DDCTMODS <- deldelCTS(DMODS, CONAVGS) #Have to switch name in ddCT function
DLOWS <- cbind(DLOWS, DDCTLOWS)
DMODS <- cbind(DMODS, DDCTMODS)

#Aggregate counts of each gene per group
COUNTL <- summary(as.factor(DLOWS$Target))
COUNTM <- summary(as.factor(DMODS$Target))

#Gene means
LMEANS <- tapply(DLOWS$DDCTLOWS, DLOWS$Target, mean)
LSEMS <- tapply(DLOWS$DDCTLOWS, DLOWS$Target, std.error)
MMEANS <- tapply(DMODS$DDCTMODS, DMODS$Target, mean)
LSEMS <- tapply(DMODS$DDCTMODS, DMODS$Target, std.error)

#Fold changes from mean values, SEM for graphing error bars
FLOWS <- 2**(-DLOWS$DDCTLOWS)
DLOWS <- cbind(DLOWS, FLOWS)
FLMEANS <- tapply(DLOWS$FLOWS, DLOWS$Target, mean)
FLSEMS <- tapply(DLOWS$FLOWS, DLOWS$Target, std.error)
FMODS <- 2**(-DMODS$DCTMODS)
DMODS <- cbind(DMODS, FMODS)
FMMEANS <- tapply(DMODS$FMODS, DMODS$Target, mean)
FMSEMS <- tapply(DMODS$FMODS, DMODS$Target, std.error)

Gene <- labels(COUNTL[c(1:9)])
GroupL <- rep("Low", 9)
GroupM <- rep("Moderate", 9)

#Group all together, means not for graph, SEM for graph
RQ <- data.frame(Gene=labels(c(FLMEANS, FMMEANS)), 
                 Group = c(rep("Low", 9), rep("Moderate",9)), 
                 Mean = c(FLMEANS, FMMEANS), 
                 SEM = c(FLSEMS, FMSEMS))

#Bargraph of the data for RQs with moderate humidity NOT standardized to 1
ggplot(data=RQ, 
       aes(x=Gene, y=Mean, ymin= Mean-SEM, ymax =Mean+SEM, fill=Group)) +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=50, vjust = 0.5, size=13)) +
  scale_fill_grey(start=0.2, end=0.6) +
  geom_errorbar(width=0.4,position = position_dodge(0.9)) +
  labs(x=NULL, y="RQ", fill = "Humidity Level")

ggplot(data=DMODS, aes(x=Target, y=FMODS)) + 
  geom_point(stat="identity") +
  theme(axis.text.x = element_text(angle=50, vjust = 0.5, size=13)) +
  scale_fill_grey(start=0.2, end=0.6)+
  labs(x=NULL, y="RQ", fill = "Humidity Level")

rm(CONAVGS, COUNTL, COUNTM, DCTLOWS, DCTMODS, DDCTLOWS, DDCTMODS, FLOWS,
   FMODS, GroupL, GroupM, LMEANS, LOWAVGS, LSEMS, MMEANS, HPLOWS, HPMODS, LOWS, 
   MEANS)

# Outlier Detection ############################################################

#Grubbs test on lows; one sided
GL1 <- tapply(DLOWS$FLOWS, DLOWS$Target, grubbs.test)
GL1OP <- tapply(DLOWS$FLOWS, DLOWS$Target, grubbs.test, opposite=TRUE)

#Need to remove ECCP(26,29) and MCP1(33)
DLOWS %>% filter(Target=="ECCP"|Target=="MCP1")
DLOWS.O1 <- DLOWS %>% filter(!(Sample == 26 & Target == "ECCP"), 
                             !(Sample == 33 & Target == "MCP1"))
tapply(DLOWS.O1$FLOWS, DLOWS.O1$Target, grubbs.test)
tapply(DLOWS.O1$FLOWS, DLOWS.O1$Target, grubbs.test, opposite=TRUE)
DLOWS.O1 %>% filter(Target=="ECCP")
DLOWS.O2 <- DLOWS.O1 %>% filter (!(Sample==29 & Target == "ECCP"))
tapply(DLOWS.O2$FLOWS, DLOWS.O2$Target, grubbs.test)

### Remove CTs for ECCP(26,29) and MCP1(33) {This is for LOWS only}
REDO.O <- REDO %>% filter(!(Sample == 26 & Target == "ECCP"), 
                          !(Sample==29 & Target == "ECCP"), 
                          !(Sample == 33 & Target == "MCP1"))

#Grubbs test on mods; one sided
tapply(DMODS$FMODS, DMODS$Target, grubbs.test)
tapply(DMODS$FMODS, DMODS$Target, grubbs.test, opposite=TRUE)
DMODS.O1 <- DMODS %>% filter(!(Sample==35&Target=="MMP12"),
                             !(Sample==39&Target=="MCP1"),
                             !(Sample==39&Target=="CDHR4"),
                             !(Sample==39&Target=="CDSM"),
                             !(Sample==39&Target=="ECCP"),
                             !(Sample==39&Target=="MUC21"))
tapply(DMODS.O1$FMODS, DMODS.O1$Target, grubbs.test)
DMODS.O1 %>% filter(Target=="ECCP")
DMODS.O2 <- DMODS.O1 %>% filter(!(Target=="ECCP"&Sample==38))
tapply(DMODS.O2$FMODS, DMODS.O2$Target, grubbs.test)

### Remove MMP12(35), MCP1-CDHR4-CDSM-MUC21(39), ECCP(38,39) {MODS only}
REDO.O2 <- REDO.O %>% filter(!(Sample==35&Target=="MMP12"),
                             !(Sample==39&Target=="MCP1"),
                             !(Sample==39&Target=="CDHR4"),
                             !(Sample==39&Target=="CDSM"),
                             !(Sample==39&Target=="ECCP"),
                             !(Sample==39&Target=="MUC21"),
                             !(Target=="ECCP"&Sample==38))

# Repeat entire Rq/FC calculation above with REDO.O2############################
#Verify outliers are taken care of
tapply(DLOWS$FLOWS, DLOWS$Target, grubbs.test, two.sided=TRUE)
tapply(DMODS$FMODS, DMODS$Target, grubbs.test)

tapply(DLOWS$FLOWS, DLOWS$Target, grubbs.test, opposite=TRUE)
tapply(DMODS$FMODS, DMODS$Target, grubbs.test, opposite=TRUE)

#Updated RQ Dataframe
RQ <- data.frame(Target=levels(REDO.O2$Target), 
                 Group = c(rep("Low", 9), rep("Moderate",9)), 
                 Mean = c(FLMEANS, FMMEANS), SEM = c(FLSEMS, FMSEMS))

# T-tests #######################################################################
#LOWS:
CDHR4L <- DLOWS %>% filter(Target=="CDHR4")
CDSML <- DLOWS %>% filter(Target=="CDSM")
ECCPL <- DLOWS %>% filter(Target=="ECCP")
MCP1L <- DLOWS %>% filter(Target=="MCP1")
MMP12L <- DLOWS %>% filter(Target=="MMP12")
MUC21L <- DLOWS %>% filter(Target=="MUC21")
SPBNL <- DLOWS %>% filter(Target=="SPBN")
ZACNL <- DLOWS %>% filter(Target=="ZACN")

#MODS:
CDHR4M <- DMODS %>% filter(Target=="CDHR4")
CDSMM <- DMODS %>% filter(Target=="CDSM")
ECCPM <- DMODS %>% filter(Target=="ECCP")
MCP1M <- DMODS %>% filter(Target=="MCP1")
MMP12M <- DMODS %>% filter(Target=="MMP12")
MUC21M <- DMODS %>% filter(Target=="MUC21")
SPBNM <- DMODS %>% filter(Target=="SPBN")
ZACNM <- DMODS %>% filter(Target=="ZACN")

t.test(CDHR4L$FLOWS, CDHR4M$FMODS, alternative = 'g')
t.test(CDSML$FLOWS, CDSMM$FMODS, alternative = 'l')
t.test(ECCPL$FLOWS, ECCPM$FMODS, alternative = 'l')
t.test(MCP1L$FLOWS, MCP1M$FMODS,  alternative = 'g')
t.test(MMP12L$FLOWS, MMP12M$FMOD, alternative = 'g')
t.test(MUC21L$FLOWS, MUC21M$FMODS, alternative = 'l')
t.test(SPBNL$FLOWS, SPBNM$FMODS, alternative = 'l')
t.test(ZACNL$FLOWS, ZACNM$FMODS, alternative = 'g')

# Figure with controls standardized to Rq 1 ####################################
#Format if necessary
REDO.O2$Sample <- as.factor(REDO.O2$Sample)
REDO.O2$Target <- as.factor(REDO.O2$Target)
REDO.O2$Plate <- as.factor(REDO.O2$Plate)
REDO.O2 <- mutate(REDO.O2, Group = ifelse(Sample %in% 34:39, "Moderate", "Low"))

#Collect HPRTs for same genes split across plates
HPRT13 <- REDO.O2 %>% filter(Target=="HPRT1" & 
                               (Plate == "Plate1" | Plate == "Plate3"))
tapply(HPRT13$Ct, HPRT13$Group, mean)

HPRT24 <- REDO.O2 %>% filter(Target=="HPRT1" & 
                               (Plate == "Plate2" | Plate == "Plate4"))
tapply(HPRT24$Ct, HPRT24$Group, mean)

HPRT5 <- REDO.O2 %>% filter(Target=="HPRT1" & (Plate == "Plate5"))
tapply(HPRT5$Ct, HPRT5$Group, mean)

#Gene means
ANALYSIS <- aggregate(REDO.O2$Ct, 
                      list(Target = REDO.O2$Target, Group = REDO.O2$Group),
                      mean)
#Avgerage HPRTs for respective gene
HAVGS <- c(24.2615,24.2615,24.2615,1,24.83451,24.9689,24.2615,24.83451,24.9689,
           24.32851,24.32851,24.32851,1,24.68314,24.95101,24.32851,24.68314,
           24.95101)
#Dataframe with all intermediate values, moderates set to FC1
ANALYSIS <- mutate(ANALYSIS, 
                   HPRT = HAVGS, dCT = x-HPRT , 
                   ddCT = dCT - c(dCT[10:18],dCT[10:18]),
                   RQ = 2**(-ddCT), FC = 1/RQ)
ANALYSIS$SEM <- RQ$SEM

#No error bars
ggplot(data=ANALYSIS, aes(x=Target, y=RQ, fill=Group)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=50, vjust = 0.5, size=13)) +
  scale_fill_grey(start=0.2, end=0.6) +
  labs(x=NULL, y="RQ", fill = "Humidity Level")

# FigStats <- read_excel("FigStats.xlsx", sheet = "Sheet2")
FIG <- ANALYSIS %>% filter(Target!="HPRT1")

# Final Figure #################################################################
ggplot(data=FIG, 
       aes(x=Target, y=RQ, ymin = RQ-SEM, ymax = RQ+SEM, fill=Group)) +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=50, vjust = 0.5, size=13)) +
  scale_fill_grey(start=0.2, end=0.6)+ 
  geom_errorbar(width=0.4,position = position_dodge(0.9)) +
  labs(x=NULL, y="RQ", fill = "Humidity Level") + 
  theme(text=element_text(size=14,  family="serif"))

######### USING NEW FUNCTIONS #######
library(openxlsx)
DF = read.xlsx("Data_OutliersRemoved.xlsx") %>% 
  mutate(Group = factor(ifelse(Sample %in% 34:39, "Moderate", "Low")),
         Plate = factor(Plate)) %>%
  rename(CT=Ct)
DF.PLATES = split(DF, DF$Plate)

targets = unique(DF$Target); targets = targets[targets != "HPRT1"]

GENES = sapply(DF.PLATES, function(x){
  sapply(intersect(targets, x$Target), function(y) {
    x %>% filter(Target %in% c(y, "HPRT1"))
  }, simplify=FALSE)},
  simplify=FALSE) %>% flatten()

COMBS = sapply(targets, function(x) do.call(bind_rows, GENES[names(GENES)==x]), simplify=FALSE)

PCR = sapply(targets, function(x) {
  CALC.PCR(COMBS[[x]], target=x, housekeeping="HPRT1", control="Moderate")
}, simplify=FALSE)

STATS = map(PCR, 1)
sapply(STATS, function(x) {
  sigma(lm(ddCT ~ Group, data=x, na.action=na.exclude))
}) %>% mean() #1.1896

FIG = do.call(bind_rows, map(PCR, 2)) %>% mutate(Gene = rep(names(PCR), each=2))
ggplot(data=FIG, aes(x=Gene, y=FC, fill=Group)) + geom_col(position="dodge")






CALC.PCR <- function(DF, target, housekeeping="HPRT1", control="Moderate", ...) {
  #Accept DATAFRAME with columns Sample, Target, and CT
  #Calculates qPCR parameters standardized to provided housekeeping gene and
  #control group designation
  
  tar = subset(DF, DF$Target == target)
  hou = subset(DF, DF$Target == housekeeping)
  
  tar = tar %>% filter(Sample %in% intersect(hou$Sample, tar$Sample))
  hou = hou %>% filter(Sample %in% intersect(hou$Sample, tar$Sample))
  
  tar.stat <- aggregate(CT ~ Sample + Group, data=tar, mean, na.rm=TRUE)
  hou.stat <- aggregate(CT ~ Sample + Group, data=hou, mean, na.rm=TRUE)
  
  tarDEL = tar.stat$CT - hou.stat$CT
  tar.stat <- mutate(tar.stat, tarDEL) 
  tar.stat <- tar.stat %>% 
    mutate(tarDELDEL = tarDEL - mean(tar.stat$tarDEL[tar.stat$Group==control], na.rm=TRUE),
           FC = 2^-tarDELDEL) %>% 
    rename(dCT = tarDEL, ddCT = tarDELDEL)
  
  agg.mean <- tapply(tar$CT, tar$Group, mean, na.rm=TRUE) - tapply(hou$CT, hou$Group, mean, na.rm=TRUE)
  fig.mean <- agg.mean - agg.mean[2]
  ddct.stderr <- with(tar.stat, tapply(ddCT, Group, function(x) sd(x)/sqrt(length(x))))
  fig.fc <- 2^-fig.mean
  fc.stderr <- with(tar.stat, tapply(FC, Group, function(x) sd(x)/sqrt(length(x))))

  fig = data.frame(Group = levels(DF$Group),
                   DDCT = fig.mean, DD.SEM = ddct.stderr,
                   FC = fig.fc, FC.SEM = fc.stderr, row.names=NULL)
 
    return(list(STATS = tar.stat, FIGURE = fig))
}