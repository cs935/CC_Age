library(DESeq2)
library(ggplot2)
library(tidyverse)
library(ggrepel)

target <- "Young" # "Young"/"Middle"/"Old"
Non_target <- "Old" # "Young"/"Middle"/"Old"
volcano.data <- res.cpm
rownames(volcano.data) <- volcano.data$name 
volcano.data <- volcano.data[complete.cases(volcano.data$padj),] 
up <- volcano.data[volcano.data$log2FoldChange>1 & volcano.data$padj<0.05,]
down <- volcano.data[volcano.data$log2FoldChange< -1 & volcano.data$padj<0.05,]
up.top <- rownames(up[order(up$padj,-up$log2FoldChange),])[1:5]
down.top <- rownames(down[order(down$padj,down$log2FoldChange),])[1:5]
top_genes <- c(up.top,down.top)
top_genes_data <- volcano.data[top_genes,]
valcano <- ggplot(data=volcano.data, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  ggtitle(paste0(target,"VS",Non_target)) + 
  scale_color_manual(name="", values=c("red", "blue", "black"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-log2(2), log2(2)), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.05), lty=2, col="gray", lwd=0.5)+
  geom_text_repel(data = top_genes_data, aes(label = name), size = 2.5,color = "black")