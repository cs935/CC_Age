library(ggplot2) 
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA)
library(dplyr)

target = "c6"
geneSetFile = "c6.all.v2023.1.Hs.symbols.gmt"
gsva_data <- read.csv("expdataTumor_TPM.csv",row.names = 1,check.names = F)
gene_set <- getGmt(geneSetFile)
gsva_result <- gsva(as.matrix(gsva_data), gene_set, method = "gsva", min.sz=1,
                    max.sz=Inf, kcdf="Gaussian", parallel.sz=1L) 
scores <- as.data.frame(t(gsva_result))
cluster <- read.table("cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
cluster$group <- factor(cluster$group,levels= c('young',"middle",'old'))
sameSample <- intersect(row.names(scores), row.names(cluster))
rt <- cbind(scores[sameSample,,drop=F], cluster[sameSample,,drop=F])
rt2 <- data.frame(score=colnames(scores),Pvalue=NA)
rownames(rt2) <- rt2[,1]
rt2 <- rt2[,-1,drop=F]
for(i in rownames(rt2)){
  result1 <- rt 
  result1 <- rename(result1, marker = i)
  test=kruskal.test(marker ~ group, data=result1) 
  diffPvalue=test$p.value
  rt2[i,1]=diffPvalue
}
positive_rt <- rt[,rownames(filter(rt2,Pvalue<0.05))]
positive_rt1 <- data.frame(matrix(ncol = length(colnames(positive_rt)), nrow = 3))
colnames(positive_rt1) <- colnames(positive_rt)
rownames(positive_rt1) <- c('Young',"Middle",'Old')
positive_rt1[1,] <- colMeans(positive_rt[1:22,])  
positive_rt1[2,] <- colMeans(positive_rt[(22+1):(22+243),])
positive_rt1[3,] <- colMeans(positive_rt[(22+243+1):(22+243+39),])
positive_rt1 <- t(scale((positive_rt1)))
ht <- Heatmap(positive_rt1, name = "z-scores", cluster_rows = F, cluster_columns = F,
        show_column_names = T, show_row_dend = FALSE, show_row_names = T)
draw(ht)