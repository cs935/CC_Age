library(limma)
library(ggplot2)
library(tidyverse)

target <- "Young" #"Young"/"Middle"/"Old"
Non_target <- "Middle" #"Young"/"Middle"/"Old"
exprSet <- read.csv("GSE26511mergeGSE52903_expData.csv",header = T,row.names = 1,check.names = F)
pdata <- read.csv("GSE26511mergeGSE52903_clinicalData.csv",header = T,row.names = 1,check.names = F)
list.group_target <- filter(pdata,cluster==target) 
group_target <- intersect(rownames(list.group_target),colnames(exprSet))
target_number <- length(group_target) 
list.group_Nontarget <- filter(pdata,cluster == Non_target) 
group_Nontarget <- intersect(rownames(list.group_Nontarget),colnames(exprSet))
Nontarget_number <- length(group_Nontarget) 
data.target <- exprSet[,group_target]
data.Nontarget<- exprSet[,group_Nontarget]
all(row.names(data.target)==row.names(data.Nontarget))
data <- cbind(data.target,data.Nontarget)
group_list<-pdata[colnames(data),]
group<-group_list[,2]
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(data)
contrast.matrix<-makeContrasts(paste0(target,"-",Non_target),levels=design) 
fit <- lmFit(data,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
write.csv(nrDEG,"GEO.limma.csv")