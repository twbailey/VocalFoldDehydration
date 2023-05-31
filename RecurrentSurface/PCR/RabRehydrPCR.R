library(readxl)
# working <- read_xlsx("Export Data/ReducedPCR.xlsx")

library(tidyverse)
# working <- working %>% mutate(PLATE = rep(1:4, each=144))

library(openxlsx)
# write.xlsx(list(PCR = working), "Export Data/ReducedPCR.xlsx")


library(plotrix)
PCR <- read_xlsx("Export Data/ReducedPCR.xlsx")

#Remove individual technical reps that do not fit in the trio.
PCR <- PCR %>% filter(is.na(FLAG)) %>% mutate(SAMPLE = as.factor(SAMPLE), TARGET = as.factor(TARGET), PLATE = as.factor(PLATE), 
                                              GROUP = factor(str_extract(PCR$SAMPLE, "^[A-Z]")),
                                              CT = round(as.numeric(CT),2)) %>% filter(!is.na(CT))

genes = levels(PCT$TARGET)
targets = sort(unique(PCR$TARGET[PCR$TARGET!="HPRT1"]) %>% droplevels())
levels(PCR$PLATE) <- targets
samples = levels(PCR$SAMPLE)
groups = levels(PCR$GROUP)
controls = grep("^C", samples, value=TRUE)
dehyds = grep("^D", samples, value=TRUE)
rehyds = grep("^R", samples, value=TRUE)


with(PCR, aggregate(list(mean.CT = CT), list(PLATE=PLATE, TARGET=TARGET, SAMPLE=SAMPLE, GROUP=GROUP), mean, na.rm=TRUE))


pcr.outs <- lapply(targets, function(x){
  
  plate <- PCR %>% filter(PLATE == x)
  means <- with(plate, aggregate(list(mean.CT = CT), list(TARGET = TARGET, SAMPLE=SAMPLE, GROUP=GROUP), mean, na.rm=TRUE))
  house <- means[which(means$TARGET == "HPRT1"),] %>% arrange(SAMPLE)
  
  target <- means[which(means$TARGET != "HPRT1"),] %>% arrange(SAMPLE)
  target <- target %>% mutate( d.CT = target$mean.CT - house$mean.CT,
                              dd.CT = d.CT - mean(d.CT[target$GROUP=="C"]),
                              RQ = 2^(-dd.CT) )
  agg.dCT <- with(target, tapply(mean.CT, GROUP, mean, na.rm=TRUE)) - with(house, tapply(mean.CT, GROUP, mean, na.rm=TRUE))
  agg.ddCT <- agg.dCT - agg.dCT["C"]
  agg.RQ <- 2^(-agg.ddCT)
  agg.stdev <- with(target, tapply(RQ, GROUP, sd, na.rm=TRUE))
  agg.stderr <- with(target, tapply(RQ, GROUP, std.error, na.rm=TRUE))
  
  agg.df <- data.frame(DCT=agg.dCT, DDCT=agg.ddCT, RQ=agg.RQ, STDEV=agg.stdev, STDER=agg.stderr)
  
  return(list(TARGET=target, AGGREGATE=agg.df))
})
names(pcr.outs) <- targets

pcr.outs <- do.call(Map, c(f=rbind, pcr.outs))
TARGET <- pcr.outs$TARGET %>% mutate(target=TARGET)
AGG <- pcr.outs$AGGREGATE
AGG <- AGG %>% mutate(GROUP = str_extract(rownames(AGG),"[A-Z]$"), TARGET = str_extract(rownames(AGG),"[A-Z0-9]+"))



ggplot(data=TARGET) + geom_point(mapping=aes(x=target, y=RQ, color=GROUP, group=GROUP), position=position_dodge2(width=0.5))

library(plotly)
ggplot(TARGET, aes(x = TARGET, fill = GROUP, y = RQ)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", binwidth=0.1) +
    geom_label(aes(label=SAMPLE), position=position_jitterdodge())


ggplotly(ggplot(data=AGG) + geom_col(mapping=aes(x=TARGET, y=RQ, fill=GROUP), position="dodge"))


ggplot(data=TARGET %>% filter(PLATE == "ABG1")) + geom_point(aes(x=SAMPLE, y=CT, color=TARGET))+facet_wrap(~GROUP, scales="free")
ggplot(data=PCR %>% filter(PLATE == "ATP1")) + geom_point(aes(x=SAMPLE, y=CT, color=TARGET))+facet_wrap(~GROUP, scales="free")
ggplot(data=PCR %>% filter(PLATE == "COL2")) + geom_point(aes(x=SAMPLE, y=CT, color=TARGET))+facet_wrap(~GROUP, scales="free")
ggplot(data=PCR %>% filter(PLATE == "IL1RAP")) + geom_point(aes(x=SAMPLE, y=CT, color=TARGET))+facet_wrap(~GROUP, scales="free")

