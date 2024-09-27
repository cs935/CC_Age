library(limma)
library(sva)

eset_1 <- read.csv("GSE26511_expData.csv",header = T,row.names = 1,check.names = F)
eset_2 <- read.csv("GSE52903_expData.csv",header = T,row.names = 1,check.names = F)
same_genes <- intersect(rownames(eset_1),rownames(eset_2))
eset_merge<- cbind(eset_1[same_genes,],eset_2[same_genes,])
col <- c(rep("#54D9ED",length(colnames(eset_1))),rep("#E93639",length(colnames(eset_2))))
boxplot(eset_merge,col=col,las=2)
group1 <- read.csv("GSE26511_clinicalData.csv",header = T,row.names = 1,check.names = F)
group1 <- group1[,c("age at diagnosis:ch1"),drop = F]
colnames(group1) <- "age"
group2 <- read.csv("GSE52903_clinicalData.csv",header = T,row.names = 1,check.names = F)
group2 <- group2[,c("characteristics_ch1.1"),drop = F]
colnames(group2) <- "age"
group2$age <- sub(".*: ", "", group2$age)
group <- rbind(group1,group2)
group$cluster <- ifelse(group$age >=30,ifelse(group$age >=65,'Old','Middle'),'Young') 
GSE<- c(rep('GSE26511',39),rep('GSE52903',55)) 
group_list<- group$cluster 
data <- eset_merge
batch <- c(rep('GSE26511',39),rep('GSE52903',55))
design <- model.matrix(~group_list)
expr_limma <- removeBatchEffect(data,batch = batch,design = design)
boxplot(expr_limma,col=col,las=2)