# #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Code supporting: 
# Bailey, T. W., dos Santos, A. P., do Nascimento, N. C., Xie, J., 
# Sivasankar, M. P., & Cox, A. (2021). Recurring exposure to low humidity
# induces transcriptional and protein level changes in the vocal folds of
# rabbits. Scientific Report. (2021)11:24180.
# https://doi.org/10.1038/s41598-021-03489-0

# Taylor W. Bailey, PhD, MS, MPH
# Updated 2 Nov 2022
# #  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(tidyverse)
library(ggplotify)
library(patchwork)
library(factoextra)
library(readxl)
library(openxlsx)
library(VennDiagram)
library(lme4)
library(emmeans)

file = "WaterRestrictRehydrationProteinsRaw.xlsx"

##### Functions ################################################################

### Median centering
median_centering = function(df) {
  #Adapted from 
  #https://datascienceplus.com/proteomics-data-analysis-2-3-data-filtering-and-
  # missing-value-imputation/
  # df = data frame containing LOG2 columns for normalization
  cols = grep("^C|^D|^R", names(df), value = TRUE)
  
  df[, cols] = lapply(cols, 
                      function(x) {
                        LOG2 = df[[x]]
                        LOG2[!is.finite(LOG2)] = NA   # Exclude missing values
                        gMedian = median(LOG2, na.rm = TRUE)
                        LOG2 - gMedian
                      }
  )
  return(df)
}

###Imputation from normal distribution
#Adapted from 
#https://datascienceplus.com/proteomics-data-analysis-2-3-data-filtering-and-
# missing-value-imputation/

imp <- function(data, shift=1.8, scale=0.3, seed=1988){
  #Will impute values randomly from a downshifted and scaled normal distribution
  # column (rabbit)-wise
  #Variance is scaled and the mean is shifted in units of the standard deviation
  
  cols = grep("^C|^D|^R", names(data), value=TRUE) #returns names of C/D/R cols
  n.cols = grep("^C|^D|^R", names(data)) #returns indices of C/D/R cols
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

### Produces Volcano Plot
pdiff.plot <- function(df, D, P){
  #Volcano plot
  ggplot(data=df, aes(x=unlist(df[D]),y=-log10(unlist(df[P])))) + 
    geom_point()+
    geom_hline(yintercept=1, col="red", size=1)+
    geom_hline(yintercept=1.3, col="blue", size=1)+
    geom_vline(xintercept=0.58,size=1, color="gray35")+ geom_vline(xintercept=-0.58, size=1, color="gray35")+
    labs(x="Mean Difference (log2)", y="-log10 P") +
    scale_color_hue(c=75)+theme(legend.title=element_blank())
}

### PCA plot through ggplot, default color by "Group"
pca.plot <- function(df, cols, n=nrow(df), title=NULL, colors=NULL){
  ## Input as variables in rows and samples in columns  
  pcomp <- prcomp(t(df[1:n,cols]), scale=TRUE)
  val <- as.data.frame(pcomp$x) #PCA vectors
  per <- pcomp$sdev^2/sum(pcomp$sdev^2) #percent variance
  grps = factor(str_extract(cols, "^[A-Z]"))
  grps = gsub("C", "Control", grps);
  grps = gsub("D", "Dehydrated", grps)
  grps = gsub("R", "Rehydrated", grps)
  val$Group <- grps
  val$name <- rownames(val)
  ggplot(data=val, aes_string(x="PC1", y="PC2", color="Group"))+
    geom_point(size=3) + geom_text(aes(label=name), hjust=0, vjust=-1)+
    xlab(paste0("PC1 (", round(per[1]*100,1),"%)"))+
    ylab(paste0("PC2 (", round(per[2]*100,1),"%)"))+
    ggtitle(title) + colors
}

### Rename Uncharacterized Proteins with Unique Numbers
new.name <- function(df, label=NULL) {#label is optional additional label for protein names
  if(is.null(label)) {add <- NULL}
  else {add <- paste0(" ",label)}
  aa=1
  for(i in 1:nrow(df)){
    if (df$NAME[i]=="Uncharacterized protein"){
      df$NAME[i] <- paste0("Uncharacterized protein ",aa,add)
      aa=aa+1
    }
  }
  return(df)
} 


### Calculate Variance Explained by First Principal Component for a Set of 
# with N's defined within the function
var.by.prot <- function(df, cols, start=5,by=5){
  s=seq(start,nrow(df),by)
  per1 = sapply(s, function(x) {
    pcomp=prcomp(t(df[1:x,cols]), scale=TRUE)
    100*pcomp$sdev[1]^2/sum(pcomp$sdev^2)})
  per2 = sapply(s, function(x) {
    pcomp=prcomp(t(df[1:x,cols]), scale=TRUE)
    100*pcomp$sdev[2]^2/sum(pcomp$sdev^2)})
  par(mfrow=c(2,2))
  plot(x=s,y=per1, main=paste0(deparse(substitute(df)),": Variance Explained by PC1 (%)"), xlab="Proteins", ylab="%Var")
  plot(x=s,y=per2, main=paste0(deparse(substitute(df)),": Variance Explained by PC2 (%)"), xlab="Proteins", ylab="%Var")
  plot(unlist(sapply(1:length(per1)-1, function(x) {per1[x]-per1[x-1]})), type="l", main="Incremental Change VE by PC1", ylab="")
  plot(unlist(sapply(1:length(per2)-1, function(x) {per2[x]-per2[x-1]})), type="l",  main="Incremental Change VE by PC2", ylab="")
  par(mfrow=c(1,1))
  
  lst <- list(s=s, per1=per1, per2=per2)
  return(lst)
}

### Add names, used in enrichment analysis
addNamesFID <- function(x, names = NAMES){
  between.bars = "(|)[a-zA-Z0-9]+(|)"
  for(i in 1:nrow(x)){
    #Extract each Hits into a single list
    a=str_extract_all(x$Hits[i],between.bars)
    #Find indices of associated Genes in NAMES
    b=match(unlist(a), names$GENE)
    #Get the names at these indices
    c=names$NAME[b]
    #Collapse them together into a single string
    d=paste(c, collapse="|")
    #Assign values to column named Names
    x$Names[i] = d
    
    #Get the FIDs at the indices above.
    e=names$FID[b]
    #Collapse them together into a single string
    f=paste(e, collapse="|")
    x$FIDs[i] = f
  }
  
  return(x)}

collapsed.groups <- function(df){
  #Function will return results for collapsed of HITS into largest unique groups.
  #Will account for GO terms, descriptions, Categories, and Groups#
  
  between.bars = "(|)[a-zA-Z0-9]+(|)"
  #Test if each full hit is a subset of any other or if any other hit is nested within given hit  
  
  for(i in 1:nrow(df)){
    z=numeric();y=numeric()
    for(j in 1:nrow(df)){
      if(all(unlist(str_extract_all(df$Hits[i], between.bars)) %in% unlist(str_extract_all(df$Hits[j], between.bars)))){
        z=append(z,j)
      }
    }
    ifelse(length(z)==1, df$SUB[i] <- 0, df$SUB[i] <- 1)
    m=which.max(sapply(z, function(x) length(unlist(str_extract_all(df$Hits[x], between.bars)))))
    ifelse(length(z)==1, df$MAX[i] <- i, df$MAX[i] <- z[m])
  }
  
  for(i in 1:nrow(df)){
    df$Cluster[i] = paste(unique(df$Cluster[which(df$MAX==df$MAX[i])]), collapse=",")
  }
  
  #Pull out the indices of the largest substrings
  U.MAX = unique(df$MAX)
  
  #Pull out descriptions for each substring of the largest strings
  a=sapply(U.MAX, function(x) df$Description[df$MAX == x])
  d=sapply(U.MAX, function(x) df$GO[df$MAX == x])
  g=sapply(U.MAX, function(x) df$Category[df$MAX == x])
  j=sapply(U.MAX, function(x) df$Cluster[df$MAX == x])
  max.desc = data.frame(MAX=U.MAX)
  for(z in 1:length(U.MAX)){
    b = sort(unlist(a[[z]]))
    c = paste(b, collapse="|")
    e = sort(unlist(d[[z]]))
    f = paste(e, collapse="|")
    h = sort(unique(unlist(g[[z]])))
    i = paste(h, collapse="|")
    max.desc$ALLS[z] <- c
    max.desc$TERMS[z] <- f
    max.desc$CATS[z] <- i
    max.desc$GROUPS[z] <- unique(j[[z]])
  }
  
  #Returns rows with the largest unique set of hits, having collapsed nested results upward.
  df.new = df[U.MAX,] %>% mutate(Descriptions = max.desc$ALLS, GOs = max.desc$TERMS, Categories = max.desc$CATS, Groups = max.desc$GROUPS)
  return(df.new)
}


intersections <- function(df, df2=NAMES){
  between.bars = "[^|]+"
  between.commas = "[^,]+"
  GO=sapply(df$GOs, function(x) {str_extract_all(x, between.commas)})
  GO.U = unique(unlist(GO), incomparables = "NA")
  CAT=sapply(df$Categories, function(x) {str_extract_all(x, between.commas)})
  CAT.U = unique(trimws(unlist(CAT)), incomparables = "NA")
  HIT=sapply(df$Hits, function(x) {str_extract_all(x, between.bars)})
  HIT.U = unique(unlist(HIT), incomparables = "NA")
  FID = sapply(df$FIDs, function(x) {str_extract_all(x, between.bars)})
  FID.U = unique(unlist(FID), incomparables = "NA")
  NAM = sapply(df$Names, function(x) {str_extract_all(x, between.bars)})
  NAM.U = unique(unlist(NAM), incomparables = "NA")
  l = list(GO = trimws(GO.U), CAT = trimws(CAT.U), NAME = trimws(NAM.U), HITS = trimws(HIT.U), FID = trim(FID.U))
  return(l)
}

separate.terms = function(df, how.many=5){
  clusts = unique(df$Cluster)
  a = lapply(clusts, function(x) {b = df %>% filter(Cluster == x)
  return(intersections(b))})}


##### Colors #####
colors1 = c("#00CC00","#CC0000") #Green, Red
names(colors1) = c("Control", "Dehydrated")
colorset1 = scale_color_manual(name="Group", values=colors1)

colors2 = c("#00CC00","#0000CC") #Green, Blue
names(colors2) = c("Control", "Rehydrated")
colorset2 = scale_color_manual(name="Group", values=colors2)

colors3 = c("#CC0000","#CC0000","#0000CC") #Green, Red, Blue
names(colors3) = c("Control", "Dehydrated", "Rehydrated")
colorset3 = scale_color_manual(name="Group", values=colors3)

colors4 = c("#CC0000","#0000CC","#6600CC","#00CC00") #Red, Blue, Purple, Green
names(colors4) = c("Dehydrated", "Rehydrated", "Both", "Neither")
colorset4 = scale_color_manual(name="SET", values=colors4)


##### Import Data ##############################################################
RAW <- read_excel(file)
PROT <- RAW #So an original unprocessed dataframe is retained.

### Column (rabbit) names for groups
names(PROT) = gsub("B","C", names(PROT)); controls = grep("^C", names(PROT), value=TRUE)
dehyds = grep("^D", names(PROT), value=TRUE)
rehyds = grep("^R", names(PROT), value=TRUE)
alls = c(controls, dehyds, rehyds)

### Pull names out of FASTA Header and pull first Uniprot ID and add universal 
# #'s to "Uncharacterized Proteins [new.name()]" 
for(i in 1:nrow(PROT))
  {PROT$NAME[i]<- str_match(PROT$FASTA[i], "RABIT (.*?) OS")[2]
   PROT$FID[i] = str_extract(PROT$ProtID[i], "^[a-zA-Z0-9]+")}
PROT$FID[1127] = "P11974-2" #Isoform not caught by the regex

which(duplicated(PROT$FID))
# [1] 1127 2482 2627 2719 2808 2929
contam = grep("^CON",PROT$ProtID) #Labeled as contaminants (CON prefix)
# [1]  119 2627 2808 2929
for(i in contam) {
  PROT$FID[i] <- paste0(str_extract_all(PROT$ProtID[i], "[a-zA-Z0-9]+")[[1]][2],"(CON)")}
PROT$FID[c(2482, 2719)] <- c("P12798-3", "P02602-2")

PROT = PROT %>% new.name()

###Count missing values in each group
for(i in 1:nrow(PROT)){
  PROT$M.C[i] = sum(PROT[i,controls]==0)
  PROT$M.D[i] = sum(PROT[i,dehyds]==0)
  PROT$M.R[i] = sum(PROT[i,rehyds]==0)
}


##### Analytical Subset  #######################################################
PROT3 = PROT %>% filter(M.D < 5 | M.R < 5 | M.C < 5) #1827

#Log transform and imputation
PROT3[alls] = log2(PROT3[alls])
PROT3 = PROT3 %>% median_centering() %>% imp()
PROT3.L <- PROT3 %>% 
  pivot_longer(all_of(alls), names_to="Rabbit", values_to="LFQ")
PROT3.L$Group <- factor(sapply(PROT3.L$Rabbit, function(x) 
  str_extract(x, "^[A-Z]{1}"), USE.NAMES=FALSE))
levels(PROT3.L$Group) <- c("Control", "Dehydration", "Rehydration")


##### Different Linear Models with Different Contrasts UNCORRECTED #####
contrALL <- list(D_C= c(-1, 1, 0),
                 R_C= c(-1, 0, 1),
                 D_R= c(0, 1, -1),
                 DR_C= c(-2, 1, 1))
PROT3.LM <- lmList(LFQ~Group|FID, data=PROT3.L, na.action=na.exclude)
mod3.order <- names(PROT3.LM)

# Takes a while  
PROT3.ALL = lapply(mod3.order, function(x) {
  as.data.frame(contrast(emmeans(PROT3.LM[[x]], ~Group), contrALL))})

names(PROT3.ALL) <- mod3.order
ALLS <- do.call(rbind, PROT3.ALL)
ALLS <- ALLS %>% mutate(FID = rep(mod3.order, each=4)) %>%
  select(-SE, -df, -t.ratio) %>% 
  pivot_wider(id_cols="FID", names_from="contrast", 
              values_from=c("estimate", "p.value")) %>%
  rename(EST_CD = estimate_D_C, 
         EST_CR = estimate_R_C, 
         EST_DR = estimate_D_R, 
         EST_CDR = estimate_DR_C, 
         P_CD = p.value_D_C, 
         P_CR = p.value_R_C, 
         P_DR = p.value_D_R, 
         P_CDR = p.value_DR_C) %>% 
  mutate(Q_CD = p.adjust(P_CD, method="fdr"), 
         Q_CR = p.adjust(P_CR, method="fdr"), 
         Q_DR = p.adjust(P_DR, method="fdr"), 
         Q_CDR = p.adjust(P_CDR, method="fdr"))

ALLS <- ALLS %>% merge(., PROT3, by="FID")
CD <- ALLS %>% select(FID, NAME, EST_CD, P_CD, Q_CD) %>% 
  filter(P_CD < 0.05) %>% merge(., PROT3, by=c("FID")) %>% arrange(P_CD)
CR <- ALLS %>% select(FID, NAME, EST_CD, P_CD, Q_CD, EST_CR, P_CR, Q_CR) %>% 
  filter(P_CR < 0.05) %>% merge(., PROT3, by="FID") %>% arrange(P_CR)
CDR <- ALLS %>% select(FID, NAME, EST_CDR, P_CDR, Q_CDR) %>% 
  filter(P_CDR < 0.05) %>% merge(., PROT3, by="FID") %>% arrange(P_CDR)
CDRC <- ALLS %>% select(FID, NAME, EST_CD, P_CD, Q_CD, EST_CDR, P_CDR, Q_CDR)%>%  
  filter(P_CD < 0.05 | P_CDR < 0.05)  %>% merge(., PROT3, by="FID")
DR <- ALLS %>% select(FID, NAME, EST_DR, P_DR, Q_DR) %>% 
  filter(P_DR < 0.05) %>% merge(., PROT3, by="FID")


##### Comparison Control to Dehydration ########################################

cd.dist <- get_dist(t(ALLS[cols1]), method="spearman")
cd.tree <- hclust(cd.dist, method="average")

# tiff("F:/Rabbit Water Restriction Rehydration Spring 2021/Proteomics R Project/Figures/ContToDehyComparison.tiff", units="in", heigh=5, width=5, res=300, compression="lzw")

pdiff.plot(ALLS %>% filter(-log10(P_CD)<5), "EST_CD", "P_CD") / 
  #left out one upper outlying value
(pca.plot(PROT3, cols1) + colorset1 + guides(color="none") + ylim(c(-40,32.5))
  + xlim(c(-40, 30)) |
  pca.plot(CD, cols1) + colorset1 + guides(color="none") + ylim(c(-15,15)) + 
    xlim(c(-20,17.5)) ) /
as.ggplot(~plot(cd.tree, ann=FALSE)) +
  plot_annotation(tag_levels = 'A')+plot_layout(heights = c(1,2,2))

#dev.off()

##### Comparison of all Three Groups #####

# tiff("F:/Rabbit Water Restriction Rehydration Spring 2021/Proteomics R Project/Figures/AllGroupComparison.tiff", units="in", heigh=5, width=5, res=300, compression="lzw")

pca.plot(ALLS, cols2) + colorset2 + guides(color="none") + 
  ylim(c(-30,30)) + xlim(c(-40,25)) +
  pca.plot(ALLS, alls) + colorset3 + guides(color="none") + 
  ylim(c(-40,30)) + xlim(c(-40,30)) +
  pca.plot(CDR, alls) + colorset3 + guides(color="none") + 
  ylim(c(-15,15)) + xlim(c(-30,20)) +
  pca.plot(CR, alls) + colorset3 + guides(color="none") + 
  ylim(c(-15,15)) + xlim(c(-30,20)) +
  plot_annotation(tag_levels = 'A')

# dev.off()

##### Other Hierachical Clustering #####

cr.dist <- get_dist(t(ALLS[cols2]), method="spearman")
cr.tree <- hclust(cr.dist, method="average")

cdr.dist <- get_dist(t(ALLS[alls]), method="spearman")
cdr.tree <- hclust(cdr.dist, method="average")

cdrp.dist <- get_dist(t(CDR[alls]), method="spearman")
cdrp.tree <- hclust(cdrp.dist, method="average")

# tiff("F:/Rabbit Water Restriction Rehydration Spring 2021/Proteomics R Project/Figures/RehyClustering.tiff", units="in", heigh=5, width=5, res=300, compression="lzw")

as.ggplot(~plot(cr.tree, ann=FALSE)) /
  as.ggplot(~plot(cdr.tree, ann=FALSE)) /
  as.ggplot(~plot(cdrp.tree, ann=FALSE)) +
  plot_annotation(tag_levels = 'A')

# dev.off()

dr.dist <- get_dist(t(ALLS[c(dehyds, rehyds)]), method="spearman")
dr.cdr.dist <- get_dist(t(CDR[c(dehyds, rehyds)]), method="spearman")
dr.p.dist <- 
  plot(hclust(get_dist(
    t(head(ALLS%>%arrange(P_DR) %>% select(dehyds,rehyds),n=500)), 
    method="spearman"), method="average"), ann=FALSE)

plot(hclust(get_dist(
  t(head(ALLS%>%arrange(P_DR)%>%select(dehyds,rehyds),n=500)),
  method="spearman"), method="average"), ann=FALSE)

dr.tree <- hclust(dr.dist, method="average")
dr.cdr.tree <- hclust(dr.cdr.dist, method="average")


##### Comparison of C and R #####

# tiff("F:/Rabbit Water Restriction Rehydration Spring 2021/Proteomics R Project/Figures/DehytoREhyComparison.tiff", units="in", heigh=5, width=5, res=300, compression="lzw")

(pca.plot(ALLS, c(rehyds,dehyds)) + colorset3 + guides(color="none") +
  xlim(c(-25,35)) + ylim(c(-20, 25)) +
pca.plot(CDR, c(rehyds,dehyds)) + colorset3 + guides(color="none") + 
  xlim(c(-15, 20)) + ylim(c(-10,15)) +
pca.plot(CR, c(rehyds,dehyds)) + colorset3 + guides(color="none") + 
  xlim(c(-12.5, 20)) + ylim(c(-15,15)) +
ggplot(ALLS) + 
  geom_point(aes(x=EST_CD, y=EST_CR, color=-log10(P_DR)), size=2) + 
  geom_abline(slope=1, intercept=0, linetype="dotted", size=1) + 
  geom_hline(yintercept=0, size=1) + geom_vline(xintercept=0, size=1)+
  geom_abline(slope=1, intercept=-0.58, size=1, linetype="dashed") + 
  geom_abline(slope=1, intercept=0.58, size=1, linetype="dashed") + 
  scale_color_distiller(palette = "RdYlBu")  + 
  theme(panel.background = element_rect(fill = "#AAAAAA")) +
  ylim(c(-5,6)) + xlim(c(-4,5)) + guides(color="none") + 
  xlab("Difference Control-Dehy") + ylab("Difference Control-Rehy")) +

  plot_annotation(tag_levels = 'A')

# dev.off()

##### Tables #####
CD <- CD %>% arrange(P_CD)
p=prcomp(t(CD[cols1]), scale=TRUE)
CD <- CD %>% mutate(C.1=p$rotation[,1] * p$sdev[1])
CD.topP <- CD %>% head(n=15)
CD.botE <- CD %>% arrange(desc(EST_CD)) %>% head(n=15)

CD.Table <- rbind(CD.topP, CD.botE) %>% 
  select(FID, NAME.x, P_CD, Q_CD, EST_CD, C.1) %>% 
  rename(NAME=NAME.x, P=P_CD, Q=Q_CD, D=EST_CD, C=C.1) %>%
  mutate(P=log10(P), Q=log10(Q))


CDR <- CDR %>% arrange(P_CDR)
p=prcomp(t(CDR[alls]), scale=TRUE)
CDR <- CDR %>% mutate(C.1=p$rotation[,1] * p$sdev[1])
CDR.topP <- CDR %>% head(n=15)
CDR.topE <- CDR %>% arrange(desc(EST_CDR)) %>% head(n=15)

CDR.Table <- rbind(CDR.topP, CDR.topE) %>% 
  select(FID, NAME.x, P_CDR, Q_CDR, EST_CDR, C.1) %>% 
  rename(NAME=NAME.x, P=P_CDR, Q=Q_CDR, D=EST_CDR, C=C.1) %>%
  mutate(P=log10(P), Q=log10(Q))

# write.xlsx(list(CD.Table=CD.Table, CDR.Table=CDR.Table), "NewTopTablesAfterLabMeeting.xlsx")

##### Enrichement #####

ONTOLOGY <- read_excel("F:/Rabbit Water Restriction Rehydration Spring 2021/Proteomics Metascape/NewAfterLabMeeting/NewALMEnrichment.xlsx", sheet="GOS")
DISGENET <- read_excel("F:/Rabbit Water Restriction Rehydration Spring 2021/Proteomics Metascape/NewAfterLabMeeting/NewALMEnrichment.xlsx", sheet="DisGeNET")
TF <- read_excel("F:/Rabbit Water Restriction Rehydration Spring 2021/Proteomics Metascape/NewAfterLabMeeting/NewALMEnrichment.xlsx", sheet="TFT")

ONTOLOGY.P <- ONTOLOGY %>% group_by(COMP) %>% arrange(LogP, .by_group=TRUE) %>% 
  slice(match(1:10, Cluster))
ONTOLOGY.E <- ONTOLOGY %>% group_by(COMP) %>% 
  arrange(desc(Enrichment), .by_group=TRUE) %>% 
  slice(match(unique(Cluster)[1:10], Cluster))
# write.xlsx(rbind(ONTOLOGY.P, ONTOLOGY.E), "NEWOntologyTopTable.xlsx")

NAMES <- read_xlsx("NEWTablesAfterLabMeeting.xlsx", sheet="FID_GENE")
NAMES <- left_join(NAMES, ALLS %>% select(FID, NAME))

CD.SEPS <- separate.terms(ONTOLOGY.CD)\CD.SEPS[[2]][58]<-""
CDR.SEPS <- separate.terms(ONTOLOGY.CDR)

DEHY <- ONTOLOGY %>% filter(COMP == "DEHY")
dehy.list <- lapply(1:5, function(x) {
  a = DEHY %>% filter(Cluster == x)
  b = unlist(a$GO)
  d = unlist(a$Description)
  e = rep(x, length(b))
  list(CLUSTER = e, GO = b, DESC = d)
})

DEHY.GOS <- do.call(rbind, lapply(1:5, 
  function(x) as.data.frame(dehy.list[[x]])))
DEHY.GOS.R <- DEHY %>% collapsed.groups() %>% filter(Cluster %in% 1:5) %>% 
  group_by(Cluster) %>% arrange(LogP) %>% slice(1:3, by.group=TRUE)

CDR.go <- ONTOLOGY %>% filter(COMP == "CDR")
cdr.list <- lapply(1:5, function(x) {
  a = CDR.go %>% filter(Cluster == x)
  b = unlist(a$GO)
  d = unlist(a$Description)
  e = rep(x, length(b))
  list(CLUSTER = e, GO = b, DESC = d)
})

CDR.GOS <- do.call(rbind, lapply(1:5, function(x) as.data.frame(cdr.list[[x]])))
CDR.GOS.R <- CDR.go %>% collapsed.groups() %>% filter(Cluster %in% 1:5) %>% 
  group_by(Cluster) %>% arrange(LogP) %>% slice(1:3, by.group=TRUE)

ONTOLOGY.names <- addNamesFID(ONTOLOGY)

BOTH.GOS.R <- rbind(DEHY.GOS.R, CDR.GOS.R)
BOTH.GOS.R <- addNamesFID(BOTH.GOS.R)

DEHYEnrich <- DEHY %>% addNamesFID()
CDREnrich <- CDR.go %>% addNamesFID()

# venn.diagram(list(DEHY = (ONTOLOGY %>% filter(COMP=="DEHY") %>% select(GO) %>% unlist()), 
#                   CDR = (ONTOLOGY %>% filter(COMP=="CDR") %>% select(GO) %>% unlist())), "test.venn.tiff")


##### PPI #####

PPIData <- read_excel("F:/Rabbit Water Restriction Rehydration Spring 2021/Proteomics Metascape/NewAfterLabMeeting/PPIData.xlsx") %>%
  mutate(Name = as.factor(Name))
with(PPIData, tapply(Name, COMP, table))
with(PPIData, tapply(Name, COMP, tail))


PPI.fil <- PPIData %>% filter(LogP <= -3 & Enrichment >= 5 & `#GeneInGOAndHitList` >= 6) %>% filter(!(Name %in% c("MyList","MyList_MCODE_ALL")))
PPI.fil5 <- PPIData %>% filter(LogP <= -3 & Enrichment >= 5 & `#GeneInGOAndHitList` >= 5) %>% filter(!(Name %in% c("MyList","MyList_MCODE_ALL")))
PPI.d  <- PPI.fil %>% filter(COMP=="DEHY")
PPI.cdr <- PPI.fil %>% filter(COMP=="CDR")


ppi.d.1 <- PPI.d %>% filter(Name == "MyList_SUB1_MCODE_1")
ppi.d.2 <- PPI.d %>% filter(Name == "MyList_SUB1_MCODE_2") %>% slice(1:3)
ppi.d.3 <- PPI.d %>% filter(Name == "MyList_SUB1_MCODE_4") %>% slice(1:3)

ppi.cdr.1 <- PPI.cdr %>% filter(Name == "MyList_SUB1_MCODE_1") %>% slice(1:3)
ppi.cdr.2 <- PPI.cdr %>% filter(Name == "MyList_SUB1_MCODE_2") %>% slice(1:3)
ppi.cdr.3 <- PPI.cdr %>% filter(Name == "MyList_SUB1_MCODE_3") %>% slice(1:3)
ppi.cdr.4 <- PPI.cdr %>% filter(Name == "MyList_SUB1_MCODE_4") %>% slice(1:3)
ppi.cdr.5 <- PPI.cdr %>% filter(Name == "MyList_SUB1_MCODE_5") %>% slice(1:3)

ppi.table <- rbind(ppi.d.1, ppi.d.2, ppi.d.3, ppi.cdr.1, ppi.cdr.2, ppi.cdr.3, ppi.cdr.4, ppi.cdr.5)
ppi.table <- addNamesFID(ppi.table)