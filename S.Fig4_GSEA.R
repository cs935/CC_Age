library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggridges)
library(ggplot2)
library(ggsci)

lim <- read.csv("GEO.limma.csv", header = T, row.names = 1, check.names = F)
symbol <- rownames(lim)
entrez <- bitr(symbol,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
genelist <- lim$logFC
names(genelist) <- rownames(lim)
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
genelist <- sort(genelist, decreasing = T)
KEGG_ges <- gseKEGG(geneList = genelist,organism = "hsa",minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.05,
  pAdjustMethod = "BH",verbose = FALSE,eps = 0)
KEGG_ges <- setReadable(KEGG_ges,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
gseaplot2(KEGG_ges,geneSetID = items,pvalue_table = TRUE,ES_geom = "line")