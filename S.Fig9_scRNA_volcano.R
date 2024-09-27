library(Seurat) 
library(tidyverse) 

Cells.sub<- subset(sce@meta.data,celltype=="Epithelial") 
scRNAsub <- subset(sce,cells=row.names(Cells.sub))
scRNAsub <- FindVariableFeatures(scRNAsub,selection.method = "vst",nfeatures=2000) 
scale.genes <- rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub,features = scale.genes)
scRNAsub <- RunPCA(scRNAsub,features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub,ndims=20,reduction="pca")
pc.num=1:10
scRNAsub <- FindNeighbors(scRNAsub,dims = pc.num) 
scRNAsub <- FindClusters(scRNAsub,resolution =0.15) 
cell_cluster <- data.frame(cell_ID=rownames(metadata),cluster_ID=metadata$seurat_clusters) 
scRNAsub <- RunUMAP(scRNAsub,dims = pc.num) 
DimPlot(scRNAsub,reduction = "umap",label = T) 
markers <- FindAllMarkers(object = scRNAsub,test.use="wilcox",
                          only.pos = F, 
                          logfc.threshold =0.25) 
all.markers =markers %>% dplyr::select(gene,everything()) %>% subset(p_val<0.05) 
jjVolcano(diffData = all.markers,
          log2FC.cutoff = 0.25,
          size =3.5,
          fontface = 'italic',
          topGeneN =3
)