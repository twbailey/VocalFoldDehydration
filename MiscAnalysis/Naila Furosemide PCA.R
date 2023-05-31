library(RColorBrewer)
library(tidyverse)

colors2 = c("#88BACD","#FFD966")
names(colors2) = c("Control", "Dehydrated")
colorset2 = scale_color_manual(name="Group", values=colors2)

PROT <- read_excel("Furosemide Proteomics_Table S1_Naila.xlsx")
rabbits <- grep("^C|^D",names(PROT), value=TRUE)

pca.plot <- function(df, cols, n=nrow(df), title=NULL, colors=NULL){
  ## Input as variables in rows and samples in columns  
  pcomp <- prcomp(t(df[1:n,cols]), scale=TRUE)
  val <- as.data.frame(pcomp$x) #PCA vectors
  per <- pcomp$sdev^2/sum(pcomp$sdev^2) #percent variance
  grps = factor(str_extract(cols, "^[A-Z]"))
  grps = gsub("C", "Control", grps);
  grps = gsub("D", "Dehydrated", grps)
  val$Group <- grps
  val$name <- rownames(val)
  ggplot(data=val, aes_string(x="PC1", y="PC2", color="Group"))+
    geom_point(size=5) + geom_text(aes(label=name,fontface = "bold"), hjust=0, vjust=-1)+
    xlab(paste0("PC1 (", round(per[1]*100,1),"%)"))+
    ylab(paste0("PC2 (", round(per[2]*100,1),"%)"))+
    ggtitle(title) + colors
}


tiff("Furosemide All Proteins.tiff", units="in", width=4, height=4, compression="lzw", res=300)
pca.plot(PROT, rabbits, colors=colorset2)+xlim(c(-30,60))+ylim(c(-30,15)) #all proteins
dev.off()

tiff("Furosemide 500 Proteins.tiff", units="in", width=4, height=4, compression="lzw", res=300)
pca.plot(head(PROT, n=500), rabbits, colors=colorset2)+xlim(c(-20,30))+ylim(c(-12.5,15)) #top 500 proteins
dev.off()

tiff("Furosemide 250 Proteins.tiff", units="in", width=4, height=4, compression="lzw", res=300)
pca.plot(head(PROT, n=250), rabbits, colors=colorset2)+xlim(c(-15,20))+ylim(c(-10,10)) #top 250 proteins
dev.off()

tiff("Furosemide 100 Proteins.tiff", units="in", width=4, height=4, compression="lzw", res=300)
pca.plot(head(PROT, n=100), rabbits, colors=colorset2)+xlim(c(-7.5,12.5))+ylim(c(-7,5)) #top 100 proteins
dev.off()

tiff("Furosemide Significant P Proteins.tiff", units="in", width=4, height=4, compression="lzw", res=300)
pca.plot(PROT %>% filter(P < 0.05), rabbits, colors=colorset2)+xlim(c(-6,6))+ylim(c(-2,4))
dev.off()