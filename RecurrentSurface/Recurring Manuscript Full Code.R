# ############################################################################# #
### R-Programming for Bailey, T. et al (2021):                                ###  
### Surface dehydration induced by recurring exposure to low humidity induces ###  
### transcriptional and protein level changes in the vocal folds of rabbits   ###
### Author: Taylor W. Bailey                                                  ###
### Updated 7.14.2021                                                         ###
# ############################################################################# # 

# ############################################################################# #
###                             Set up Environment                            ###
#################################################################################

##### Load libraries to be used #####
library(ggpubr);library(ggsignif);library(lme4);library(nlme);library(outliers);library(openxlsx);
library(plotrix);library(RColorBrewer);library(readxl);library(stringr);library(tidyverse);
library(VennDiagram)
library(plotly)
library(htmlwidgets)

##### Set working directory #####
datatables = "F:/Recurring Manuscript Revisions/Recurring Manuscript Revisions R Project/Revision Data Tables"
setwd(datatables)

# ############################################################################# #
##                           Environmental Conditions                          ##
# ############################################################################# #
##### Import Data #####------------------------------------------------------------
HUM = read_xlsx("RecurringFullHumidity.xlsx")
HUM = HUM %>% filter(!(COHORT %in% c("HIST 1", "HIST 2")))
##### Data Summaries #####------------------------------------------------------------
HUMSum <- with(HUM, aggregate(list(Summary=RH), list(GROUP=GROUP,COHORT=COHORT), summary)) %>% 
  do.call(data.frame,.) %>% arrange(GROUP)
colnames(HUMSum) <- c("Group","Cohort","Min","Q1","Median","Mean","Q3","Max")
HUMSum$Mean = round(HUMSum$Mean,1)

GROUPSum <-with(HUM, aggregate(list(Summary=RH), list(GROUP=GROUP), summary)) %>% do.call(data.frame,.)
colnames(GROUPSum) <- c("Group","Min","Q1","Median","Mean","Q3","Max")
GROUPSum$Mean = round(GROUPSum$Mean,1)
with(HUM, tapply(RH, GROUP, sd, na.rm=TRUE))
with(HUMSum,tapply(Mean, Group, mean, na.rm=TRUE))
with(HUM, aggregate(list(Summary=RH), list(GROUP=GROUP, COHORT=COHORT), sd)) %>% do.call(data.frame,.)

rm(HUM, HUMSum, GROUPSum)
# ############################################################################# #
###                                PCV Analysis                               ###
# ############################################################################# #
##                                    PCV                                      ##
#################################################################################
##### Import Data ##### --------------------------------------------------------------------------------
ALL = read_excel("Chronic Surface PCV Data.xlsx",col_types = c("numeric", "numeric", "numeric","text", "text")) %>% 
  mutate(Rabbit = as.factor(Rabbit), Day = as.factor(Day), Cohort = as.factor(Cohort), Group = as.factor(Group))
#Averages of technical replicates for each measurement
MEAN <- with(ALL, aggregate(list(PCV=PCV), list(Group=Group, Cohort=Cohort, Rabbit=Rabbit, Day=Day), mean, na.rm=TRUE))
PCV <- ALL %>% filter(Cohort != "HIST 1") %>% droplevels()

levels(PCV$Day) = c(1,8,15)
levels(PCV$Group) <- c("Moderate", "Low");PCV$Group = factor(PCV$Group, levels=c("Low","Moderate"))
levels(PCV$Cohort)=c("C","D","E","A","B");PCV$Cohort = factor(PCV$Cohort, levels=(c("A","B","C","D","E")))


##### Mixed Models ##### --------------------------------------------------------------------------------
PCV$Day <- as.numeric(PCV$Day)
PCV$CohortN <- as.numeric(PCV$Cohort)
#options(contrasts=c("contr.SAS", "contr.SAS"))
m1 = lmer(PCV ~Day*Group + (Day|Cohort/Rabbit), data=PCV, na.action=na.exclude)
m2 = lmer(PCV ~Day*Group + (1|Cohort/Rabbit), data=PCV, na.action=na.exclude)#2.2-16
m3 = lmer(PCV ~Day*Group + (1|Cohort), data=PCV, na.action=na.exclude) #<2.2-16
m4 = lme(PCV~Day*Group,random=list(Cohort=~Day, Rabbit=~Day), data=PCV, correlation=corAR1(), na.action=na.exclude)
m5 = lme(PCV~Day*Group,random=list(Cohort=~1, Rabbit=~Day), data=PCV, correlation=corAR1(), na.action=na.exclude)
m6 = lme(PCV~Day*Group,random=list(Cohort=~1, Rabbit=~1), data=PCV, correlation=corAR1(), na.action=na.exclude)
ml = lm(PCV ~Day*Group, data=PCV, na.action=na.exclude)

m1c = lmer(PCV~ Day*Group+Cohort+(Day|Cohort:Rabbit), data=PCV, na.action=na.exclude)
m1c2 = lmer(PCV~ Day*Group+Cohort+(1|Cohort:Rabbit), data=PCV, na.action=na.exclude)

m1ckr = update(m1c, REML=FALSE)
m1ckr2 = update(m1ckr, .~.-Cohort)
anova(m1c,m1c2)

mla = lm(PCV~Day*Group*Cohort, data=PCV, na.action=na.exclude)
mla2 = lm(PCV~(Day+Group+Cohort)^2, data=PCV, na.action=na.exclude)
mlc = lm(PCV~Day*Group + Cohort, data=PCV, na.action=na.exclude)

m1cN = lmer(PCV ~ Day*Group + CohortN + (Day|Cohort/Rabbit), data=PCV, na.action=na.exclude)

mgc = lmer(PCV ~ Group*Cohort+(Day|Cohort:Rabbit), data=PCV, na.action=na.exclude)
mgcf = lmer(PCV ~ Group*Cohort+(DayF|Cohort:Rabbit), data=PCV, na.action=na.exclude)


##### Percent Change ##### --------------------------------------------------------------------------------

PCV.W = pivot_wider(MEAN, names_from = Day, values_from = PCV) %>% mutate(PCT = (`3`-`1`)/(`1`)*100) %>% filter(Cohort != "HIST 1")
levels(PCV.W$Group) <- c("Moderate","Low"); PCV.W$Group = factor(PCV.W$Group, levels=c("Low", "Moderate"))

CONT = PCV.W %>% filter(Group=="CONT")
EXP = PCV.W %>% filter(Group=="EXP")
c(mean(EXP$PCT, na.rm=TRUE)-std.error(EXP$PCT, na.rm=TRUE)*1.96,mean(EXP$PCT, na.rm=TRUE)+std.error(EXP$PCT, na.rm=TRUE)*1.96)
c(mean(CONT$PCT, na.rm=TRUE)-std.error(CONT$PCT, na.rm=TRUE)*1.96,mean(CONT$PCT, na.rm=TRUE)+std.error(CONT$PCT, na.rm=TRUE)*1.96)
t.test(EXP$PCT, CONT$PCT)
t.test(EXP$PCT, CONT$PCT, alternative="g")

t.test(EXP$PCT)
t.test(CONT$PCT)

plot(PCV.W$PCT~PCV.W$Group)
write.xlsx(PCV.W, "PCVPercentChange.csv")

ggplot(data=PCV.W) + geom_point(mapping=aes(x=Group, y=PCT, color=Cohort), size=5)

### Model Comparison
anova(m3,m2,m1)

### New Simulations
#Create 4 parallel simulations using model 1 to assess subjective quality
sim = simulate(m1, 4)
sim1 = PCV; sim1$PCV = sim$sim_1; sim2 = PCV; sim2$PCV = sim$sim_2
sim3 = PCV; sim3$PCV = sim$sim_3; sim4 = PCV; sim4$PCV = sim$sim_4
s=list(sim1,sim2,sim3,sim4)

lapply(s, function(x) {ggplot()+stat_summary(data=na.exclude(x), fun=mean, geom="line", mapping=aes(x=Day, y=PCV, group=Rabbit, color=Group))+
    facet_wrap(~Cohort)})

### Refit Original Data to assess subjective quality
a1=refit(m2, sim1$PCV)
a2=refit(m2, sim2$PCV)
a3=refit(m2, sim3$PCV)
a4=refit(m2, sim4$PCV)
lapply(list(a1,a2,a4,a4,m3), function(x) summary(x)$coefficients)


### Account for correlated measures also
library(nlme)
m4 = lme(PCV~Day*Group,random=list(Cohort=~Day, Rabbit=~Day), data=PCV, correlation=corAR1(), na.action=na.exclude)
m5 = lme(PCV~Day*Group,random=list(Cohort=~1, Rabbit=~Day), data=PCV, correlation=corAR1(), na.action=na.exclude)


# ############################################################################# #
###                                 RT-qPCR                                   ###
#################################################################################
##### Functions #####------------------------------------------------------------
#This function will take as input a dataframe with qPCR data inlcuding
#Sammple names, target gene names, and a housekeeping gene, g1 and g2 are names of controls and experimental groups, respectively.
#Output will be a list with a dataframe with columns (1) Individual sample Ct means, (2)DCt's, and (3)DDCt's
#and a second dataframe with Aggregate means by group (with control set to 0).
qPCR.OUTPUT <- function(DF, target, housekeeping="HPRT1", g1="CONT", g2="EXP", ...) {
  library(tidyverse)
  
  tar = subset(DF, DF$Target== target)
  hou = subset(DF, DF$Target== housekeeping)
  
  tar.stat <- aggregate(tar$CT, list(Sample = tar$Sample, Group = tar$Group), mean, na.rm=TRUE)
  hou.stat <- aggregate(hou$CT, list(Sample = hou$Sample, Group = hou$Group), mean, na.rm=TRUE)
  tarDEL = tar.stat$x - hou.stat$x
  tar.stat <- mutate(tar.stat, tarDEL)
  tar.stat <- mutate(tar.stat, tarDELDEL = tarDEL - mean(tar.stat$tarDEL[tar.stat$Group==g1], na.rm=TRUE))
  agg.mean <- tapply(tar$CT, tar$Group, mean, na.rm=TRUE) - tapply(hou$CT, hou$Group, mean, na.rm=TRUE)
  fig.mean <- agg.mean - agg.mean[1]
  return(list(tar.stat, fig.mean))
}

##### Import Data #####----------------------------------------------------------
library(readxl)
ONE <- read_excel("Chronic PCR Results.xls", sheet = "AQP1,4", col_types = c("text","text", "numeric", "text"))
TWO <- read_excel("Chronic PCR Results.xls", sheet = "MMP1", col_types = c("text","text", "numeric", "text"))
THREE <- read_excel("Chronic PCR Results.xls", sheet = "BDK, SCNNA1", col_types = c("text","text", "numeric", "text"))
FOUR <- read_excel("Chronic PCR Results.xls", sheet = "MUC4", col_types = c("text","text", "numeric", "text"))
FIVE <- read_excel("Chronic PCR Results.xls", sheet = "MMP12, ZACN", col_types = c("text","text", "numeric", "text"))
SIX <- read_excel("Chronic PCR Results.xls", sheet = "SCNNA, SLC26", col_types = c("text","text", "numeric", "text"))
SEVEN <- read_excel("Chronic PCR Results.xls", sheet = "AQP5(2), MUC5AC(2)", col_types = c("text","text", "numeric", "text"))
EIGHT <- read_excel("Chronic PCR Results.xls", sheet = "CFTR(2)", col_types = c("text","text", "numeric", "text"))

##### Individual Gene Data Frames #####------------------------------------------
AQP1 <- qPCR.OUTPUT(ONE, "AQP1")
AQP4 <- qPCR.OUTPUT(ONE, "AQP4")
MMP1 <- qPCR.OUTPUT(TWO, "MMP1")
BDK <- qPCR.OUTPUT(THREE, "BDKR2B")
MUC4 <- qPCR.OUTPUT(FOUR, "MUC4")
MMP12 <- qPCR.OUTPUT(FIVE, "MMP12")
ZACN <- qPCR.OUTPUT(FIVE, "ZACN")
SCNNA <- qPCR.OUTPUT(SIX, "SCNNA")
SLC26 <- qPCR.OUTPUT(SIX, "SLC26")
AQP5 <- qPCR.OUTPUT(SEVEN, "AQP5")
MUC5AC <- qPCR.OUTPUT(SEVEN, "MUC5AC")
CFTR <- qPCR.OUTPUT(EIGHT, "CFTR(2)")

DATA <- list(AQP1, AQP4, AQP5, BDK, CFTR, MMP1, MMP12, MUC4, MUC5AC, SCNNA, SLC26, ZACN)
DELDELS <- data.frame(row.names=c("C2", "C3", "C4", "C5", "C6", "L1", "L2", "L3", "L4", "L5", "L6"))
for(i in 1:12){DELDELS[,i] <- DATA[[i]][[1]]$tarDELDEL}
colnames(DELDELS) <- c("AQP1", "AQP4", "AQP5", "BDK", "CFTR", "MMP1", "MMP12", "MUC4", "MUC5AC", "SCNNA", "SLC26", "ZACN")

##### Outlier Analysis #####-----------------------------------------------------

CONT <- DELDELS[1:5,]
LOW <- DELDELS[6:11,]

apply(CONT, 2, grubbs.test, two.sided=TRUE)
#MMP1, C5(0.481)
#MUC4 C3(-1.687)
apply(LOW, 2, grubbs.test, two.sided=TRUE)
#ZACN L1(3.00)

#Removal of identified outliers
TWO.R <- subset(TWO, TWO$Sample!="C5")
MMP1.R <- qPCR.OUTPUT(TWO.R, "MMP1")
grubbs.test(MMP1.R[[1]][[5]][1:4], two.sided=TRUE)

FOUR.R <- subset(FOUR, FOUR$Sample!="C3")
MUC4.R <- qPCR.OUTPUT(FOUR.R, "MUC4")
grubbs.test(MUC4.R[[1]][[5]][1:4], two.sided=TRUE)

FIVE.R <- subset(FIVE, FIVE$Sample!="L1")
ZACN.R <- qPCR.OUTPUT(FIVE.R, "ZACN")
grubbs.test(ZACN.R[[1]][[5]][7:11])

#Removal of identified outliers
DELDELS.R = DELDELS
DELDELS.R$MMP1 = c(MMP1.R[[1]][[5]][1:4],NaN, MMP1.R[[1]][[5]][5:10])
DELDELS.R$MUC4 = c(MUC4.R[[1]][[5]][1:4], NaN, MUC4.R[[1]][[5]][5:10])
DELDELS.R$ZACN = c(ZACN.R[[1]][[5]][1:5], NaN, ZACN.R[[1]][[5]][6:10])

##### Extract information from the Wilcoxon-tests--------------------------------
RQS.R = 2**(-DELDELS.R)
CONT.R = RQS.R[1:5,]
LOW.R = RQS.R[6:11,]
CONTDD.R = DELDELS.R[1:5,]
LOWDD.R = DELDELS.R[6:11,]

M <- list()
for(i in 1:12) {M[[i]] <- wilcox.test(CONTDD.R[,i], LOWDD.R[,i])}
P <- vector()
for(i in 1:12){P[i] <- M[[i]][3]$p.value}

##### Data for Graph###########

DATA.R = list(AQP1, AQP4, AQP5, BDK, CFTR, MMP1.R, MMP12, MUC4.R, MUC5AC, SCNNA, SLC26, ZACN.R)
GroupMeans <- data.frame(matrix(NA,nrow=12,ncol=2))

for(i in 1:12){GroupMeans[i,] <- DATA.R[[i]][[2]]}
GroupRQ <- 2**(-GroupMeans)
CONTSEM.R <- apply(CONT.R, 2, std.error, na.rm=TRUE)
LOWSEM.R <- apply(LOW.R, 2, std.error, na.rm=TRUE)
SEM.R <- c(CONTSEM.R, LOWSEM.R)
Genes = c("AQP1", "AQP4", "AQP5", "BDKR2B", "CFTR", "MMP1", "MMP12", "MUC4", "MUC5AC", "SCNNA1", "SLC26A9", "ZACN")
Group <- c(rep("Moderate",12),rep("Low",12))
GRAPH <- data.frame(Gene=Genes, Group=Group, Mean = c(GroupRQ[,1], GroupRQ[,2]), SEM=SEM.R)

rm(list=ls())
# ############################################################################# #
###                                Proteomics                                 ###
# ############################################################################# #
###                              Proteomics                                   ###
#################################################################################
##### Functions ##### --------------------------------------------------------------------------------
### Median centering
median_centering = function(df) {
  #Adapted from 
  #https://datascienceplus.com/proteomics-data-analysis-2-3-data-filtering-and-missing-value-imputation/
  # df = data frame containing LOG2 columns for normalization
  cols = grep("^C|^L", names(df), value = TRUE)
  
  df[, cols] = lapply(cols, 
                      function(x) {
                        LOG2 = df[[x]]
                        LOG2[!is.finite(LOG2)] = NA   # Exclude missing values from median calculation
                        gMedian = median(LOG2, na.rm = TRUE)
                        LOG2 - gMedian
                      }
  )
  return(df)
}

###Imputation from normal distribution
#Adapted from 
#https://datascienceplus.com/proteomics-data-analysis-2-3-data-filtering-and-missing-value-imputation/
imp <- function(data, shift=1.8, scale=0.3, seed=1988){
  #Will impute values randomly from a downshifted and scaled normal distribution column (rabbit)-wise
  cols = grep("^C|^L", names(data), value=TRUE) #returns names of C/L cols
  n.cols = grep("^C|^L", names(data)) #returns indices of C/L cols
  i.cols = sapply(cols, function(x) paste("i",x, sep='')) #creates "iCOL" variable
  
  #logical output for is.finite for column names from above
  data[i.cols] = sapply(n.cols, function(x) !is.finite(as.matrix(data[, x])))
  
  set.seed(seed)
  data[,cols] = lapply(cols,
                       function(x) {
                         temp = data[[x]]
                         temp[!is.finite(temp)] = NA #replace -Inf with NA
                         
                         temp.sd = scale * sd(temp, na.rm = TRUE)   # shrink sd width
                         temp.mean = mean(temp, na.rm = TRUE) - 
                           shift * sd(temp, na.rm = TRUE) # shift mean of imputed values
                         
                         n.missing = sum(is.na(temp))
                         temp[is.na(temp)] = rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                         return(temp)
                       })
  return(data)
}

### logP v Diff Plot
pdiff.plot <- function(df){
  ggplot(data=df,aes(x=D,y=logP,text=paste0(Name,"\n","D= ",D,"\n","LogP= ",logP))) + 
    geom_point(aes(col=PCA), size=1.5)+
    labs(x="Mean Difference", y="log10 P") +
    geom_hline(yintercept=-1, col="red", size=1.25)+
    geom_hline(yintercept=-1.3, col="blue", size=1.25)+
    geom_vline(xintercept=0.58,size=1.25)+ geom_vline(xintercept=-0.58, size=1.25)+
    scale_color_hue(c=75)+theme(legend.title=element_blank())
}

pdiff.plot2 <- function(df){
  ggplot(data=df,aes(x=D,y=logP,text=paste0(Name,"\n","D= ",D,"\n","LogP= ",logP))) + 
    geom_point(aes(col=C),size=1.5)+
    labs(x="Mean Difference", y="log10 P") +
    geom_hline(yintercept=-1, col="red", size=1.25)+
    geom_hline(yintercept=-1.3, col="blue", size=1.25)+
    geom_vline(xintercept=0.58,size=1.25)+ geom_vline(xintercept=-0.58, size=1.25)+
    scale_color_distiller(palette = "RdYlBu")+theme(legend.title=element_blank())
}


### PCA plot through ggplot
pca.plot <- function(df, cols, n=nrow(df), title){
  if (!require(stringr)){library(stringr)}
  pcomp = prcomp(t(df[1:n, cols]), scale=TRUE)
  per = pcomp$sdev^2/sum(pcomp$sdev^2) #percent variance
  val <- data.frame(pcomp$x) #PC individual coords
  val$name <- cols
  for(i in 1:length(cols)){ifelse(str_detect(cols[i],"^C"), val$Group[i] <-"Moderate", val$Group[i] <- "Low")}
  ggplot(data=val, aes_string(x="PC1", y="PC2", color="Group"), inherit.aes=FALSE)+
    geom_point(size=3) + geom_text(aes(label=name), hjust=0, vjust=-1)+
    xlab(paste0("PC1 (", round(per[1]*100,1),"%)"))+
    ylab(paste0("PC2 (", round(per[2]*100,1),"%)"))+
    ggtitle(title)
}

### Rename Uncharacterized Proteins with Unique Numbers
new.name <- function(df, label=NULL) {#label is optional additional label for protein names
  aa=1
  for(i in 1:nrow(df)){
    if (df$Name[i]=="Uncharacterized protein"){
      df$Name[i] <- paste0("Uncharacterized protein ",aa," ",label)
      aa=aa+1
    }
  }
  return(df)
}  

ttests <- function(df, exp, cont) {
  #"low" and "cont" are string vectors containing column names
  T = list(nrow(df))
  for(i in 1:nrow(df)){
    T[[i]] <- t.test(df[i,exp], df[i,cont])
    df$T[i] <- round(T[[i]]$statistic,2)
    df$P[i] <- round(T[[i]]$p.value,5)
    df$D[i] <- round(mean(unlist(df[i,exp])),3)-mean(unlist(df[i,cont]))
    df$CL[i] <- round(T[[i]]$conf.int[[1]],3)
    df$CU[i] <- round(T[[i]]$conf.int[[2]],3)
    df$logP[i] <- round(log10(T[[i]]$p.value),5)
    #DIFF = EXP - CONT --> (+) mean higher expression in EXP, (-) means lower expression in EXP
  }
  return(df)
}

### Calculate Variance Explained by First Principal Component for a Set of Protein Values
var.by.prot <- function(df, cols, start=5,by=5){
  s=seq(start,nrow(df),by)
  per1 = sapply(s, function(x) {
    pcomp=prcomp(t(df[1:x,rabbits.new]), scale=TRUE)
    100*pcomp$sdev[1]^2/sum(pcomp$sdev^2)
  })
  plot(x=s,y=per1, main=paste0(deparse(substitute(df)),": Variance Explained by PC1 (%)"), xlab="Proteins", ylab="%Var")
}

##### Import Data ##### --------------------------------------------------------------------------------
ALL <- read_excel("AllOldNewOnlyFiltered.xlsx",sheet = "AllOnlyFiltered") #2454
OLD <- read_excel("AllOldNewOnlyFiltered.xlsx",sheet = "OLDOnlyFiltered") #1002
NEW <- read_excel("AllOldNewOnlyFiltered.xlsx",sheet = "NewOnlyFiltered") #2466
for(i in 1:nrow(ALL)){ALL$Name[i]<- str_match(ALL$`Fasta headers`[i], "RABIT (.*?) OS")[2]
ALL$FID[i] = str_extract(ALL$`Protein IDs`[i], "^[a-zA-Z0-9]+")}
for(i in 1:nrow(OLD)){OLD$Name[i]<- str_match(OLD$`Fasta headers`[i], "RABIT (.*?) OS")[2]
OLD$FID[i] = str_extract(OLD$`Protein IDs`[i], "^[a-zA-Z0-9]+")}
for(i in 1:nrow(NEW)){NEW$Name[i]<- str_match(NEW$`Fasta headers`[i], "RABIT (.*?) OS")[2]
NEW$FID[i] = str_extract(NEW$`Protein IDs`[i], "^[a-zA-Z0-9]+")}

#Pull rabbit column names into a character vector
rabbits.all = grep("^C|^L", names(ALL), value=TRUE) 
mod.all = grep("^C", rabbits.all, value=TRUE)
low.all = grep("^L", rabbits.all, value=TRUE)
old.all = grep("7$|8$|9$", rabbits.all, value=TRUE)
rabbits.new = grep("^C|^L", names(NEW), value=TRUE)
mod.new = grep("^C", rabbits.new, value=TRUE)
low.new = grep("^L", rabbits.new, value=TRUE)

##### Counting Missing Values ##### ---------------------------------------------------------------------
### ALL: Cohorts Combined, 18 Rabbits ###
for(i in 1:nrow(ALL)){
  z=0
  for(j in 1:6){
    if(ALL[i,j]==0) {z=z+1}
  }
  ALL$N.MODM[i] <- z
}
for(i in 1:nrow(ALL)){
  z=0
  for(j in 7:9){
    if(ALL[i,j]==0) {z=z+1}
  }
  ALL$O.MODM[i] <- z
}
for(i in 1:nrow(ALL)){
  z=0
  for(j in 10:15){
    if(ALL[i,j]==0) {z=z+1}
  }
  ALL$N.LOWM[i] <- z
}
for(i in 1:nrow(ALL)){
  z=0
  for(j in 16:18){
    if(ALL[i,j]==0) {z=z+1}
  }
  ALL$O.LOWM[i] <- z
}
ALL$M.LOW = ALL$O.LOWM+ALL$N.LOWM
ALL$M.MOD = ALL$O.MODM+ALL$N.MODM
ALL$M.OLD = ALL$O.LOWM+ALL$O.MODM
ALL$M.NEW = ALL$N.LOWM+ALL$N.MODM

### NEW: Second Cohort, 12 Rabbits ###
for(i in 1:nrow(NEW)){
  z=0
  for(j in 1:6){
    if(NEW[i,j]==0) {z=z+1}
  }
  NEW$N.MODM[i] <- z
}
for(i in 1:nrow(NEW)){
  z=0
  for(j in 7:12){
    if(NEW[i,j]==0) {z=z+1}
  }
  NEW$N.LOWM[i] <- z
}

### OLD: First Cohort, 6 Rabbits ###
for(i in 1:nrow(OLD)){
  z=0
  for(j in 1:3){
    if(OLD[i,j]==0) {z=z+1}
  }
  OLD$O.MODM[i] <- z
}
for(i in 1:nrow(OLD)){
  z=0
  for(j in 4:6){
    if(OLD[i,j]==0) {z=z+1}
  }
  OLD$O.LOWM[i] <- z
}

ALL.M =ALL %>% filter(M.NEW+M.OLD != 18)
  NEW.M =NEW %>% filter(N.MODM+N.LOWM != 12)
  OLD.M = OLD %>% filter(O.MODM+O.LOWM != 6)
nrow(ALL %>% filter(M.NEW+M.OLD != 18)) #1696
nrow(NEW %>% filter(N.MODM+N.LOWM != 12)) #1685
nrow(OLD %>% filter(O.MODM+O.LOWM != 6)) #980
rm(i,j,z)

### Venn diagram for FASTA headers
# venn.diagram(x=list(OLD$`Fasta headers`,NEW$`Fasta headers`,ALL$`Fasta headers`),
#              filename="VennAllOldNewFASTA.tiff", category.names=c("Pilot","Exp 3","Combined"),
#              fill=brewer.pal(3,"Set1"),
#              resolution = 300, compression="lzw",
#              #height=, width=
#              cat.fontfamily="sans", cat.fontface="bold", cat.cex=3,
#              fontfamily="sans", fontface="bold", cex=2)

### Venn diagram for FASTA headers after removing All-Zero-Entries
venn.diagram(x=list(OLD.M$`Fasta headers`,NEW.M$`Fasta headers`,ALL.M$`Fasta headers`),
             filename="VennAllOldNewFASTA(AllZeroRemoved).tiff", category.names=c("Pilot","Exp 3","Combined"),
             fill=brewer.pal(3,"Set1"),
             resolution = 300, compression="lzw",
             #height=, width=
             cat.fontfamily="sans", cat.fontface="bold", cat.cex=3,
             fontfamily="sans", fontface="bold", cex=2)

##### TOG: Subset of ALL with at least 5 valid values (non-zero) for each humidity group AND at least 2 valid values in either group from OLD cohort #####

TOG <- ALL %>% filter(M.LOW<=4 & M.MOD<=4 & (O.LOWM <=1|O.MODM <=1)) %>% rename(FASTA=`Fasta headers`) %>%
  select(rabbits.all, Name, FID,`Protein IDs`,`Majority protein IDs`,id, FASTA, O.LOWM, O.MODM, N.LOWM, N.MODM, M.OLD, M.NEW, M.LOW, M.MOD) %>%
  new.name(., label="(ALL)")

### Centering and Imputation from downshifted normal
TOG[rabbits.all]=log2(TOG[rabbits.all]) #Log2 transforms LFQ values (results in -Inf)
TOG = TOG %>% median_centering() %>% imp()

### Create tables to validate imputed values
b = pivot_longer(TOG, rabbits.all, names_to="Rabbit", values_to="LFQ") #Single row per protein per rabbit
c = pivot_longer(TOG, iC20:iL9, names_to="iRabbit", values_to="Impute") %>% select(Impute) #Indication for if value was imputed
b = b %>% mutate(Impute=c$Impute)
#Create group variables
for(i in 1:nrow(b)){ifelse(b$Rabbit[i] %in% mod.all, b$HGroup[i]<-"MOD", b$HGroup[i]<-"LOW")}
for(i in 1:nrow(b)){ifelse(b$Rabbit[i] %in% old.all, b$EGroup[i]<-"OLD", b$EGroup[i]<-"NEW")}

#Tables for validating distributions of imputed values
ggplot()+geom_histogram(data=b, mapping=aes(x=LFQ, color=Impute), bins=50)+facet_wrap(~Rabbit, nrow=2)+ ggtitle("Imputations by Rabbit")
xtabs(~Impute+Rabbit, b)
ggplot()+geom_histogram(data=b, mapping=aes(x=LFQ, color=Impute), bins=50)+facet_wrap(~HGroup) +  ggtitle("Imputations by Humidity Group") 
xtabs(~Impute+HGroup,b)
ggplot()+geom_histogram(data=b, mapping=aes(x=LFQ, color=Impute), bins=50)+facet_wrap(~EGroup)+ ggtitle("Imputations by Experimental Cohort")
xtabs(~Impute+EGroup,b) 
rm(b,c)
##### TOG Unsupervised Learning: PCA based on set of proteins arranged by increasing P-value ##### -------
### T.test by row ###
TOG = TOG %>% ttests(exp=low.all, cont=mod.all) %>% arrange(P)

# PCA plots, function defined above
ggarrange(pca.plot(TOG, rabbits.all, title="TOG: Proteins by P-value"),
          pca.plot(TOG, rabbits.all, n=500, title="Top 500"),
          pca.plot(TOG, rabbits.all, n=250, title="Top 250"),
          pca.plot(TOG, rabbits.all, n=125, title="Top 125"),
          pca.plot(TOG, rabbits.all, n=75, title="Top 75"),
          pca.plot(TOG, rabbits.all, n=25, title="Top 25"),
          common.legend=TRUE, legend="bottom")

ggarrange(pca.plot(TOG, rabbits.all, n=110, title="TOG: Top 110 Proteins by P-value"),
          pca.plot(TOG, rabbits.all, n=105, title="Top 105"),
          pca.plot(TOG, rabbits.all, n=100, title="Top 100"),
          pca.plot(TOG, rabbits.all, n=95, title="Top 95"),
          pca.plot(TOG, rabbits.all, n=90, title="Top 90"),
          pca.plot(TOG, rabbits.all, n=85, title="Top 85"),
          common.legend=TRUE, legend="bottom")

pcomp.tog = prcomp(t(TOG[1:95, rabbits.all]), scale=TRUE)
km.tog2 = kmeans(pcomp.tog$x[,1:2], centers=2)
km.tog4 = kmeans(pcomp.tog$x[,1:2], centers=4)

# ggarrange(fviz_cluster(km.tog2, pcomp.tog$x[,1:2], main="TOG k=2"),
#          fviz_cluster(km.tog4, pcomp.tog$x[,1:2], main="TOG k=4"))

#Select top 95 proteins as arranged by ascending P-value and pull correlation to PC1		  
T.P95 <- head(TOG, n=95) %>% mutate(C = pcomp.tog$rotation[,1] * pcomp.tog$sdev[1]) %>%
  select(Name, P, D, CL, CU, C, N.MODM, N.LOWM, FID, `Protein IDs`, `Majority protein IDs`, id, FASTA, rabbits.all, logP)
rm(km.tog2, km.tog4, pcomp.tog)
###Heat map for top 25 proteins by mean difference
TOG.H <- TOG %>% mutate(aD = abs(D)) %>% arrange(desc(aD)) %>%
  select(Name,P,D,aD,,rabbits.all, FID) %>% head(n=50) %>% arrange(D)
TOG.H <- TOG.H %>% pivot_longer(rabbits.all, names_to="Rabbit", values_to="LFQ") %>% arrange(Rabbit)
TOG.H$Protein = factor(rep(1:50,18))
levels(TOG.H$Protein) <- paste0(TOG.H$FID[1:50]," (p = ",round(TOG.H$P[1:50],3),")")

ggplot(TOG.H, aes(x=Rabbit, y=Protein, fill=LFQ))+geom_raster()+scale_fill_distiller(palette = "RdYlBu")+
  labs(title="Top 50 Proteins TOG by Mean Difference between Groups (P-value)")

###Heat map for top 25 proteins by mean difference from T.P95
TOG.H2 <- T.P95 %>% mutate(aD = abs(D)) %>% arrange(desc(aD)) %>%
  select(Name, P, D, aD, rabbits.all, FID) %>% head(n=50) %>% arrange(D)
TOG.H2 <- TOG.H2 %>% pivot_longer(cols=rabbits.all, names_to="Rabbit", values_to="LFQ") %>% arrange(Rabbit)
TOG.H2$Protein = factor(rep(1:50,18))
levels(TOG.H2$Protein) <- paste0(TOG.H2$FID[1:50]," (p = ",round(TOG.H2$P[1:50],3),")")
for(i in 1:nrow(TOG.H2)){ifelse(TOG.H2$Rabbit[i] %in% rabbits.new, TOG.H2$Group[i] <- "NEW", TOG.H2$Group[i] <- "OLD")}
TOG.H2$Group = as.factor(unlist(TOG.H2$Group))

ggplot(TOG.H2, aes(x=Rabbit, y=Protein, fill=LFQ))+geom_raster()+scale_fill_distiller(palette = "RdYlBu")+
  labs(title="Top 50 Proteins T.P95 by Mean Difference between Groups (P-value)")

rm(TOG.H, TOG.H2)

##### Analysis of NEW: NEW 3 Subset Generation and Imputation, at least 3 in either group##### ---------------
NEW3 <- NEW %>% filter(N.MODM <=3 | N.LOWM <=3) %>% rename(FASTA=`Fasta headers`) %>%
  select(rabbits.new,Name,FID,`Protein IDs`,`Majority protein IDs`,id, FASTA, N.LOWM, N.MODM) %>%
  new.name(., label="(NEW)") # 1466

###Centering and Imputation from downshifted normal
NEW3[rabbits.new]=log2(NEW3[rabbits.new]) #Log2 transforms LFQ values (results in -Inf)
NEW3 = NEW3 %>% median_centering() %>% imp()

#Create tables to validate imputed values
b = pivot_longer(NEW3, rabbits.new, names_to="Rabbit", values_to="LFQ") #Single row per protein per rabbit
c = pivot_longer(NEW3, iC20:iL25, names_to="iRabbit", values_to="Impute") %>% select(Impute) #Indication for if value was imputed
b = b %>% mutate(Impute=c$Impute)
#Create group variables
for(i in 1:nrow(b)){ifelse(b$Rabbit[i] %in% mod.new, b$HGroup[i]<-"MOD", b$HGroup[i]<-"LOW")}

#Groups to validate distribution of imputed values
ggplot()+geom_histogram(data=b, mapping=aes(x=LFQ, color=Impute), bins=50)+facet_wrap(~Rabbit, nrow=2)+ ggtitle("NEW3 Imputations by Rabbit")
xtabs(~Impute+Rabbit, b)
ggplot()+geom_histogram(data=b, mapping=aes(x=LFQ, color=Impute), bins=50)+facet_wrap(~HGroup) +  ggtitle("NEW3 Imputations by Humidity Group") 
xtabs(~Impute+HGroup,b)

rm(i,b,c)

##### NEW3 Unsupervised Learning ##### ----------------------------------------------------------------------------
### T.test by row ###
NEW3 = NEW3 %>% ttests(exp=low.new, cont=mod.new) %>% arrange(P)
#ggarrange(pca.plot(NEW3, rabbits.new, title="NEW3: Proteins by P-value"),
#          pca.plot(NEW3, rabbits.new, n=500, title="Top 500"),
#          pca.plot(NEW3, rabbits.new, n=250, title="Top 250"),
#          pca.plot(NEW3, rabbits.new, n=125, title="Top 125"),
#          pca.plot(NEW3, rabbits.new, n=75, title="Top 75"),
#          pca.plot(NEW3, rabbits.new, n=25, title="Top 25"),
#          common.legend=TRUE, legend="bottom")

#var.by.prot(NEW3, rabbits.new)
#var.by.prot(NEW3, rabbits.new, start=5, stop=500, b=5)
#var.by.prot(NEW3, rabbits.new, start=200, stop=500, b=5)

s = seq(50,500,50)
prcomp.NEW3 = lapply(s, function(x) prcomp(t(NEW3[1:x, rabbits.new]), scale=TRUE))
kmeans.NEW3 = lapply(prcomp.NEW3, function(a) kmeans(a$x[,1:2],centers=2))
for(i in 1:length(kmeans.NEW3)){fviz_cluster(kmeans.NEW3[[i]], prcomp.NEW3[[i]]$x[,1:2], centers=2)}

n = 515 #Classifies into correct groups
pca.plot(NEW3, rabbits.new, n=n, title="Top 515: NEW3")
pcomp = prcomp(t(NEW3[1:n, rabbits.new]), scale=TRUE)
kmean = kmeans(pcomp$x[,1:2], centers=2)
#fviz_cluster(kmean,pcomp$x[,1:2], main=paste0(n," Proteins"))
100*pcomp$sdev[1]^2/sum(pcomp$sdev^2) #33.3% of variance explained
100*pcomp$sdev[2]^2/sum(pcomp$sdev^2) #14.4
100*pcomp$sdev[3]^2/sum(pcomp$sdev^2) #8.2

kmean2 = kmeans(pcomp$x[,1], center=2)
pcomp$rotation[,1] #first PC

N.P515 = head(NEW3, n=515) %>% mutate(C = pcomp$rotation[,1] * pcomp$sdev[1]) %>%
select(Name, P, D, CL, CU, C,N.MODM, N.LOWM, FID, `Protein IDs`, `Majority protein IDs`, id, FASTA, logP)
#because data was centered and scaled, variances are 1 and these are perfectly correlated with PCA weights
pos = which(N.P515$C > 0) 
neg = which(N.P515$C <= 0)
N.P515POS = N.P515[pos,] %>% filter(P <= 0.1) %>% arrange(desc(C))
N.P515NEG = N.P515[neg,] %>% filter(P <= 0.1) %>% arrange(C)

rm(kmean, kmean2, kmeans.NEW3, pcomp, prcomp.NEW3, s, n)

##### NEW3 Heatmap ##### -------------------------------------------------------------------------------------------
###Heat map for top 50 proteins by mean difference
#NEW3, find and arranged by absolute differences, take top 50 and arrange by difference.
NEW3.H = NEW3 %>% mutate(aD=abs(D)) %>% arrange(desc(aD)) %>% head(n=50) %>% 
  select(Name,P,D,aD,T,starts_with("C"),starts_with("L"), FID) %>% arrange(D)
NEW3.H <- NEW3.H %>% pivot_longer(rabbits.new, names_to="Rabbit", values_to="LFQ") %>% arrange(Rabbit)
NEW3.H$Protein = factor(rep(1:50,12))
levels(NEW3.H$Protein) <- paste0(NEW3.H$FID[1:50]," (p = ",round(NEW3.H$P[1:50],3),")")

ggplot(NEW3.H, aes(x=Rabbit, y=Protein, fill=LFQ))+geom_tile()+scale_fill_distiller(palette = "RdYlBu")+
  labs(title="Top 25 Proteins NEW3 by Mean Difference between Groups (P-value)")

### Top 25 Proteins by mean difference within the p <= 0.1 subset
NEW3.H2 = NEW3 %>% filter(P <= 0.1) %>% mutate(aD=abs(D)) %>% arrange(desc(aD)) %>% head(n=50) %>%
  select(Name,P,D,aD,T,starts_with("C"),starts_with("L"), FID) %>% arrange(D)
NEW3.H2 <- NEW3.H2 %>% pivot_longer(rabbits.new, names_to="Rabbit", values_to="LFQ") %>% arrange(Rabbit)
NEW3.H2$Protein = factor(rep(1:50,12))
levels(NEW3.H2$Protein) <- paste0(NEW3.H2$FID[1:50]," (p = ",round(NEW3.H2$P[1:50],3),")")

ggplot(NEW3.H2, aes(x=Rabbit, y=Protein, fill=LFQ))+geom_tile()+scale_fill_distiller(palette = "RdYlBu")+
  labs(title="NEW3: Top 25 Proteins by Mean Difference from Top 250 by P (P-value)")+
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text(angle = 45, vjust=0.5))+
  labs(y=NULL)

rm(NEW3.H, NEW3.H2)

##### Joint Groups Analyses on Significant Subsets ##### ------------------------------------------------------------
overlap = na.omit(match(T.P95$FASTA, N.P515$FASTA)) #45
TN.P = T.P95 %>% filter(FASTA %in% N.P515$FASTA) %>%
  mutate(N.P = N.P515$P[overlap],N.D=N.P515$D[overlap],N.CL=N.P515$CL[overlap],N.CU=N.P515$CU[overlap],N.C=N.P515$C[overlap])

nrow(TN.P %>% filter(P <= 0.05)) #7
nrow(T.P95 %>% filter(P <= 0.05)) #7
nrow(N.P515 %>% filter(P <= 0.05)) #124
#Of these, 98 have a |PC1 correlation| greater than 0.044
#72 have higher expression in LOW and 26 have higher expression in MOD 
#Will analyze downstream from the N.P515POS and N.P515NEG sets.


# ############################################################################# #
##                           Gene Enrichment                                   ##
#################################################################################
##### Functions #####------------------------------------------------------------
#Will expand cluster variables into their unique values
intersections <- function(df, df2=NAMES){
  GO=sapply(df$GO, function(x) {str_extract_all(x, between.commas.GO)})
  GO.U = unique(unlist(GO))
  CAT=sapply(df$Category, function(x) {str_extract_all(x, between.commas)})
  CAT.U = unique(trimws(unlist(CAT)))
  HIT=sapply(df$Hits, function(x) {str_extract_all(x, between.bars)})
  HIT.U = unique(unlist(HIT))
  NAM.U = df2$NAME[match(HIT.U, df2$GENE)]
  FID.U = df2$PROT[match(HIT.U, df2$GENE)]
  l = list(GO = GO.U, CAT = CAT.U, NAME = NAM.U, HITS = HIT.U, FID = FID.U)
  return(l)
}
 
##### Import Data ##### --------------------------------------------------------------------------------
POS.GO <- read.csv("POS Final GO.csv", stringsAsFactors=TRUE)
NEG.GO <- read.csv("NEG Final GO.csv", stringsAsFactors=TRUE)
NAMES <- read_excel("UniProtGeneNames.xlsx")
#PROTEIN- Uniprot ID, GENE- Gene mapped to by UniProt, GROUP- POS or NEG correlated to PC1

##### POSITIVE LIST ##### ------------------------------------------------------------
table(POS.GO$Category)
length(unique(POS.GO$Hits)) #313

#Identify the set of the largest unique hits
between.bars = "(|)[a-zA-Z0-9]+(|)"

#Test if each full hit is a subset of any other or if any other hit is nested within given hit
for(i in 1:nrow(POS.GO)){
  z=numeric()
  for(j in 1:nrow(POS.GO)){
    if(all(unlist(str_extract_all(POS.GO$Hits[i], between.bars)) %in% unlist(str_extract_all(POS.GO$Hits[j], between.bars)))){
      z=append(z,j)
    }
  }
  ifelse(length(z)==1, POS.GO$SUB[i] <- 0, POS.GO$SUB[i] <- 1)
  m=which.max(sapply(z, function(x) length(unlist(str_extract_all(POS.GO$Hits[x], between.bars)))))
  ifelse(length(z)==1, POS.GO$MAX[i] <- i, POS.GO$MAX[i] <- z[m])
};rm(z)
#Pull out the indices of the largest substrings
U.MAX = unique(POS.GO$MAX)
length(U.MAX) #131

#Pull out descriptions for each substring of the largest strings
a=sapply(U.MAX, function(x) POS.GO$Description[POS.GO$MAX == x])
d=sapply(U.MAX, function(x) POS.GO$GO[POS.GO$MAX == x])
g=sapply(U.MAX, function(x) POS.GO$Category[POS.GO$MAX == x])
max.desc = data.frame(MAX=U.MAX)
for(z in 1:length(U.MAX)){
  b = unlist(a[[z]])
  c = paste(b, collapse=", ")
  e = unlist(d[[z]])
  f = paste(e, collapse=", ")
  h = unique(unlist(g[[z]]))
  i = paste(h, collapse=", ")
  max.desc$ALLS[z] <- c
  max.desc$TERMS[z] <- f
  max.desc$CATS[z] <- i
}

POS.U = POS.GO[U.MAX,] %>% mutate(Description = max.desc$ALLS, GO = max.desc$TERMS, Category = max.desc$CATS)
rm(z,a,b,c,d,e,f,g,h,i,U.MAX, max.desc)

##### NEGATIVE LIST ##### ------------------------------------------------------------
table(NEG.GO$Category)
length(unique(NEG.GO$Hits)) #136

#Test if each full hit is a subset of any other or if any other hit is nested within given hit
for(i in 1:nrow(NEG.GO)){
  z=numeric()
  for(j in 1:nrow(NEG.GO)){
    if(all(unlist(str_extract_all(NEG.GO$Hits[i], between.bars)) %in% unlist(str_extract_all(NEG.GO$Hits[j], between.bars)))){
      z=append(z,j)
    }
  }
  ifelse(length(z)==1, NEG.GO$SUB[i] <- 0, NEG.GO$SUB[i] <- 1)
  m=which.max(sapply(z, function(x) length(unlist(str_extract_all(NEG.GO$Hits[x], between.bars)))))
  ifelse(length(z)==1, NEG.GO$MAX[i] <- i, NEG.GO$MAX[i] <- z[m])
};rm(z)
#Pull out the indices of the largest substrings
U.MAX = unique(NEG.GO$MAX)
length(U.MAX) #131

#Pull out descriptions for each substring of the largest strings
a=sapply(U.MAX, function(x) NEG.GO$Description[NEG.GO$MAX == x])
d=sapply(U.MAX, function(x) NEG.GO$GO[NEG.GO$MAX == x])
g=sapply(U.MAX, function(x) NEG.GO$Category[NEG.GO$MAX == x])
max.desc = data.frame(MAX=U.MAX)
for(z in 1:length(U.MAX)){
  b = unlist(a[[z]])
  c = paste(b, collapse=", ")
  e = unlist(d[[z]])
  f = paste(e, collapse=", ")
  h = unique(unlist(g[[z]]))
  i = paste(h, collapse=", ")
  max.desc$ALLS[z] <- c
  max.desc$TERMS[z] <- f
  max.desc$CATS[z] <- i
}

NEG.U = NEG.GO[U.MAX,] %>% mutate(Description = max.desc$ALLS, GO = max.desc$TERMS, Category = max.desc$CATS)
rm(z,a,b,c,d,e,f,g,h,i, U.MAX)

##### Add Protein Names ##### ------------------------------------------------------------
#Requires N.P515 from PROTEOMICS above.
NAMES$NAME = N.P515$Name[na.omit(match(NAMES$PROT,N.P515$FID))]
for(i in 1:nrow(POS.U)){
  a=str_extract_all(POS.U$Hits[i],between.bars)
  b=match(unlist(a), NAMES$GENE)
  c=NAMES$NAME[b]
  d=paste(c, collapse=", ")
  POS.U$Names[i] = d
  e=NAMES$PROT[b]
  f=paste(e, collapse=", ")
  POS.U$FID[i] = f
}
rm(a,b,d,c,e,f,i)
for(i in 1:nrow(NEG.U)){
  a=str_extract_all(NEG.U$Hits[i],between.bars)
  b=match(unlist(a), NAMES$GENE)
  c=NAMES$NAME[b]
  d=paste(c, collapse=", ")
  NEG.U$Names[i] = d
  e=NAMES$PROT[b]
  f=paste(e, collapse=", ")
  NEG.U$FID[i] = f
}
rm(a,b,d,c,e,f,i)

##### Cluster List ##### ------------------------------------------------------------ 
CLUST = read_excel("CombinedUniqueClusters.xlsx")
between.commas.GO = "[^,\\s]+" #not commas or white space
between.commas = "[^,]+"

CLUSTS = CLUST %>% mutate(Cluster = as.factor(Cluster))
clusters = levels(CLUSTS$Cluster)
CLUSTU = lapply(clusters, function(x) {intersections(CLUSTS %>% filter(Cluster == x))})
names(CLUSTU) <- clusters

CLUSTU.df = data.frame(CLUSTER = unlist(sapply(1:8, function(x) {rep(clusters[x], length(unlist(CLUSTU[[x]][3])))})),
                       NAME = unlist(sapply(1:8, function(x) {append(NULL, CLUSTU[[x]][3])})),
                       HITS = unlist(sapply(1:8, function(x) {append(NULL, CLUSTU[[x]][4])})),
                       FID = unlist(sapply(1:8, function(x) {append(NULL, CLUSTU[[x]][5])})))
                       
COMB = merge(CLUSTU.df, N.P515) %>% select(CLUSTER, FID, Name,P, D, CL,CU,D,C) %>% arrange(CLUSTER, P)
##### Combined Cluster and NEW Statistics ##### ------------------------------------------------------------
COMB = merge(CLUSTU.df, N.P515) %>% select(-Name, -starts_with("i")) %>% mutate(CLUSTER = as.factor(CLUSTER))
xtabs(~CLUSTER, COMB)
xtabs(~CLUSTER[which(COMB$P <= 0.05)], COMB)

pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[1]),rabbits.new]), cols=rabbits.new, title=levels(COMB$CLUSTER)[1])
pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[2]),rabbits.new]), cols=rabbits.new, title=levels(COMB$CLUSTER)[2])
pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[3]),rabbits.new]), cols=rabbits.new, title=levels(COMB$CLUSTER)[3])
pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[4]),rabbits.new]), cols=rabbits.new, title=levels(COMB$CLUSTER)[4])
pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[5]),rabbits.new]), cols=rabbits.new, title=levels(COMB$CLUSTER)[5])
pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[6]),rabbits.new]), cols=rabbits.new, title=levels(COMB$CLUSTER)[6])
pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[7]),rabbits.new]), cols=rabbits.new, title=levels(COMB$CLUSTER)[7])
pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[8]),rabbits.new]), cols=rabbits.new, title=levels(COMB$CLUSTER)[8])

kmeanplot = function(i){
  pcomp = prcomp(t(COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[i]),rabbits.new]), scale=TRUE)
  k=kmeans(pcomp$x[,1:2], centers=2)
  fviz_cluster(k,pcomp$x[,1:2], main=levels(COMB$CLUSTER)[i])} 
kmeanplot(1); kmeanplot(2); kmeanplot(3); kmeanplot(4); kmeanplot(5); kmeanplot(6); kmeanplot(7); kmeanplot(8)

COMB2 = COMB %>% select(FID, NAME, P, D, CL, CU, C, CLUSTER) %>% arrange(CLUSTER, P)
#write.xlsx(COMB2, "F:/Recurring Manuscript Revisions/Proteomics/Revision Analysis/Final Analysis/Cluster Summaries/COMB2.xlsx")

rm(CLUST, CLUSTS, CLUSTU, CLUSTU.df, COMB)

# ############################################################################# #
###                          Manuscript Figures                               ###
#################################################################################

#Figure 2: Aggregate Cohort Environmental Conditions
a <- HUM %>%  na.omit() %>% mutate(GROUP = as.factor(GROUP), COHORT = as.factor(COHORT), DAY = as.factor(DAY))
levels(a$GROUP) = c("Low", "Moderate"); levels(a$COHORT)=c("C","D","E","A","B")
a$COHORT = factor(a$COHORT, levels=(c("A","B","C","D","E")))


tiff("Recurrent Humidity Aggregate Data (for Revision).tiff", units="in", height=4, width=8, res=300, compression="lzw")
ggplot(data=a, mapping=aes(x=COHORT,y=RH)) + geom_boxplot() + facet_wrap(~GROUP) +
  labs(x="Experimental Cohort", y="Relative Humidity (%)") +theme(text=element_text(size=14,  family="sans")) +
  theme(axis.text.x = element_text(size=13))
dev.off()

### Supplementary Figure 1: 15-day Environmental Conditions by Cohort
tiff("Recurrent Humidity Exposures by Day (for Revision).tiff", units="in", height=10, width=8, res=300, compression="lzw")
ggplot(data=a, mapping=aes(x=DAY, y=RH))+geom_boxplot()+facet_grid(COHORT~GROUP)
dev.off()


### Figure 4: PCR
#Figure that has comparison bars and stars
tiff("Recurrent PCR Data.tiff", units="in", height=4, width=8, res=300, compression="lzw")

ggplot(data=GRAPH, aes(x=Gene, y=Mean, ymin = Mean-SEM, ymax = Mean+SEM, fill=Group)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=50, vjust = 0.5, size=13)) +
  #scale_fill_grey(start=0.2, end=0.6)+
  geom_errorbar(width=0.4,position = position_dodge(0.9)) + 
  labs(x=NULL, y="Relative Quantification", fill = "Humidity Group") +
  theme(text=element_text(size=14,  family="sans")) + 
  
  geom_signif(stat="identity", 
              aes(x=c(NA, NA, NA, NA, NA, NA, NA, NA ,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 7.55,NA, NA, 10.55, NA),
                  xend=c(NA, NA, NA, NA, NA, NA, NA, NA ,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 8.45,NA, NA, 11.45, NA), 
                  y=c(NA, NA, NA, NA, NA, NA, NA, NA ,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, Mean[20]+SEM[20]+0.5,NA, NA, Mean[23]+SEM[23]+0.5, NA), 
                  yend=c(NA, NA, NA, NA, NA, NA, NA, NA ,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, Mean[20]+SEM[20]+0.5,NA, NA, Mean[23]+SEM[23]+0.5, NA), 
                  annotation=c(NA, NA, NA, NA, NA, NA, NA, NA ,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "",NA, NA, "", NA)),
              tip_length=0) +
  geom_point(aes(x=8, y=Mean[20]+SEM[20]+1), size=4, shape=8, show.legend=FALSE) +
  geom_point(aes(x=11,y=Mean[23]+SEM[23]+1), size=4, shape=8, show.legend=FALSE)

dev.off()

### Figure 5: Principal Component Analysis
tiff("Fig4. Principal Component Analyses.tiff", units="in", res=300, width=8, height=6, compression="lzw")
ggarrange(pca.plot(TOG, rabbits.all, n=95, title="a"),
          pca.plot(NEW3, rabbits.new, n=515, title="b"),
          pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[3]),rabbits.new]), cols=rabbits.new, title="c"),
          pca.plot((COMB[which(COMB$CLUSTER==levels(COMB$CLUSTER)[5]),rabbits.new]), cols=rabbits.new, title="d"),
          common.legend=TRUE, legend="bottom")
dev.off()

### Figure 6: Heat Maps
tiff("NEW3HeatMaps50Prot.tiff", units="in", res=300, heigh=8, width=8, compression="lzw")
ggarrange(ggplot(NEW3.H, aes(x=Rabbit, y=Protein, fill=LFQ))+geom_tile()+scale_fill_distiller(palette = "RdYlBu")+labs(y=NULL)+
            labs(title="a ")+theme(axis.text.x = element_text(angle = 45, vjust=0.5), legend.position="none", plot.title.position = "plot"),
          ggplot(NEW3.H2, aes(x=Rabbit, y=Protein, fill=LFQ))+geom_tile()+scale_fill_distiller(palette = "RdYlBu")+labs(y=NULL)+
            labs(title="b ")+scale_y_discrete(position = "right")+theme(axis.text.x = element_text(angle = 45, vjust=0.5), legend.position="none"))
dev.off()

### Supplementary Figure 2: PCV
tiff("Revision PCV Supp 2.tiff", units="in", compression="lzw", width=6, height=4, res=300)
a=ggplot()+stat_summary(data=na.exclude(PCV), fun=mean, geom="line", aes(x=Day, y=PCV, group=Rabbit, color=Group), size=1.25)+facet_wrap(~Cohort)+
  labs(x="Day of Measurement", y="PCV (%)", color="Humidity Group") + ggtitle("a ")+theme(legend.position="bottom")
b=ggplot(data=PCV.W) + geom_boxplot(mapping=aes(x=Group, y=PCT)) + ggtitle("b ") + labs(y="%Change in PCV", x=NULL)
ggarrange(a,b, widths=c(1, 0.5))
dev.off()

### Supplementary Figure 3: Proteomics Work-flow

for(i in 1:nrow(TOG)){
  ifelse(TOG$FID[i] %in% T.P95$FID,TOG$PCA[i]<-"in PCA",TOG$PCA[i]<- "Not in PCA")
  next}
for(i in 1:nrow(NEW3)){
  ifelse(NEW3$FID[i] %in% N.P515$FID,NEW3$PCA[i]<-"in PCA",NEW3$PCA[i]<- "Not in PCA")
  next}
for(i in 1:nrow(N.P515)){
  ifelse(NEW3$FID[i] %in% N.P515$FID,NEW3$PCA[i]<-"in PCA",NEW3$PCA[i]<- "Not in PCA")
  next}


ggarrange(pdiff.plot(TOG), pdiff.plot(NEW3), pdiff.plot(T.P95), pdiff.plot(N.P515))
tiff("Revision pdiff TOG.tiff", units="in", compression="lzw", width=4, height=3, res=300)
pdiff.plot(TOG)
dev.off()

tiff("Revision pdiff NEW3.tiff", units="in", compression="lzw", width=4, height=3, res=300)
pdiff.plot(NEW3)
dev.off()

tiff("Revision pdiff NP515.tiff", units="in", compression="lzw", width=4, height=3, res=300)
pdiff.plot2(N.P515)
dev.off()

#####
NEW3
for(i in 1:nrow(NEW3)){
  NEW3$SDPool[i] = sd_pooled(NEW3[i,low.new], NEW3[i,mod.new])
  NEW3$
}

e.newsdlow = ecdf(NEW3$SD.LOW)
e.newsdmod = ecdf(NEW3$SD.MOD)
plot(e.newsdlow, main="eCDF of NEW3$SD.LOW", xlab="Standard Deviation")
plot(e.newsdmod, main="eCDF of NEW3$SD.MOD", xlab="Standard Deviation")
hist(NEW3$SD.LOW);hist(NEW3$SD.MOD)
summary(NEW3$SD.LOW); summary(NEW3$SD.MOD)

NEW3 = NEW3 %>% mutate(CV =)

for(i in 1:nrow(NEW3)){
  NEW3$POWER[i] = pwr.2p2n.test(h=)
}
pwr.2p.test()