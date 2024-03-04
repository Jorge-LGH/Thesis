load("4_Intermediate/rna.RData")

Idents(rna) <- rna$RNA_snn_res.0.7

levels(factor(rna$RNA_snn_res.0.7))
clusters <- c(0:13)

rna <- rna[,rna$RNA_snn_res.0.7 %in% clusters]
rna <- NormalizeData(rna)
rna <- ScaleData(rna,features = rownames(rna))
markers <- FindAllMarkers(rna,only.pos =T)
markers <- markers[markers$p_val_adj <= 0.01,]

markers.top <- markers %>% group_by(cluster) %>% top_n(n=30)

# Downsample cells from each cluster
rna.sub <- subset(rna,downsample =457)
rna.sub <- NormalizeData(rna.sub)
rna.sub <- ScaleData(rna.sub,features = rownames(rna.sub))
rna.sub <- FindVariableFeatures(rna.sub,nfeatures = 2000)


rna.pseudobulk <- data.frame(rownames(rna.sub))
for (i in levels(factor(rna.sub$RNA_snn_res.0.7))){
  cells <- rownames(dplyr::filter(as.data.frame(rna.sub@meta.data),RNA_snn_res.0.7 == i))
  
  rna.sub.new <- rna.sub@assays$RNA@data[,colnames(rna.sub@assays$RNA@data) %in% cells]
  
  rna.bulk <- Matrix::rowSums(rna.sub.new)
  
  rna.pseudobulk$i <- rna.bulk
  colnames(rna.pseudobulk)[dim(rna.pseudobulk)[2]] <- i
  
  print("Iteration complete")
}
mat <- rna.pseudobulk
rownames(mat) <- mat$rownames.rna.sub.
mat <- mat[,-1]
head(rownames(mat))

mat <- mat[rownames(mat) %in% markers.top$gene,]
#mat <- mat[VariableFeatures(rna.sub),]

cluster_anno<-levels(factor(rna.sub$RNA_snn_res.0.7))
library(viridis)
library(ComplexHeatmap)
col_fun = circlize::colorRamp2(c(-2, 0, 2),viridis(n = 3))

heatmapRNA <- Heatmap(t(scale(t(mat))), name = "Expression",
                      cluster_rows = T,
                      show_row_dend = F,
                      clustering_method_columns = "ward.D2",
                      clustering_method_rows="ward.D2",
                      col = col_fun)
heatmapRNA

