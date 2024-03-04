#--------------------Libraries---------------------
# Load necessary libraries for analysis
library(dplyr) 
library(Seurat) 
library(ggplot2) 
library(tidyverse)
library(psych)
library(SingleR)

#--------------------Load previously created seurat object---------------------
readRDS("1_Data/Created_data/rna_postdoublet_SkipChecks.rds")

#--------------------Further cell type annotation---------------------
# SingleR predicted cell types 
rna$SingleR <- rna$SingleR.endo
rna$Sample <- "3533EL"

# Adding module scores top Seurat object
rna <- AddModuleScore(rna,features = list(panglaodb$`B cells`,
                                          panglaodb$`Plasma cells`,
                                          panglaodb$`Mast cells`,
                                          panglaodb$Macrophages,
                                          panglaodb$`Dendritic cells`,
                                          panglaodb$`T cells`,
                                          panglaodb$`NK cells`,
                                          panglaodb$`Endothelial cells`,
                                          panglaodb$Fibroblasts,
                                          panglaodb$`Epithelial cells`,
                                          panglaodb$`Smooth muscle cells`,
                                          c("TPSB2","TPSAB1","KIT")),#Three gene Mast signature
                      name = c("B.","Plasma.","Mast.","Macrophage.","DC.",
                               "T.","NK.","Endothelial.","Fibroblast.","Epithelial.","Smooth_muscle.","Mast_3_gene."),search = T)

#### Continuation of cell annotation #

# Assess Mast cell enrichment to potentially rename clusters
mast1 <- StackedVlnPlot(rna,features = c("B.1","Plasma.2","Mast.3","Macrophage.4",
                                         "DC.5","T.6","NK.7","Endothelial.8","Fibroblast.9",
                                         "Epithelial.10","Smooth_muscle.11"))
mast1
mast2 <- StackedVlnPlot(rna,features = c("TPSB2","TPSAB1","KIT"))
mast2
vln.df <- VlnPlot(rna,features = "Mast_3_gene.12")                                 # Expression of specific genes pertaining to mast cells
vln.df
data.mast <- describeBy(vln.df$data$Mast_3_gene.12, vln.df$data$ident, mat = TRUE) # Describe grouping by variable in clusters
data.mast <- dplyr::filter(data.mast,median > 0.225)                               # Threshold for mast cell assignation
data.mast                                                                          # No cell complies with threshold

# Assess B cell enrichment to potentially rename clusters
vln.df <- VlnPlot(rna,features = "B.1")                                            # Violin plot for expression of B cells' genes
vln.df
data.B <- describeBy(vln.df$data$B.1, vln.df$data$ident, mat = TRUE)
data.B <- dplyr::filter(data.B,median > 0.225)
data.B                                                                             # No cell complies with threshold

# Annotate mast/b cells
rna$mast.cell <- ifelse(rna$RNA_snn_res.0.7 %in% as.character(data.mast$group1),TRUE,FALSE) # No cells were annotated as mast cell
rna$B.cell <- ifelse(rna$RNA_snn_res.0.7 %in% as.character(data.B$group1),TRUE,FALSE)       # No cells were annotated as B cell

# Append SingleR annotations to cluster labels:
# The most common SingleR label in each cluster becomes the cluster label 
cells <- rna@meta.data %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(SingleR) # cluster| Cell type| number of cells
cluster.ids <- rep("fill",length(levels(Idents(rna))))                                # Make a vector of 14 "fills"
names(cluster.ids) <- levels(Idents(rna))                                             # Assign cluster to each element in vector

#### Aggregate all cells pertaining to their specific clusters #

# Group all cell types per cluster and number of cells for each cell type
# Group the data per cluster per cell type to each vector in clusters
for ( i in factor(cells$RNA_snn_res.0.7)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.7 ==i) %>% arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$SingleR[1]
}
cluster.ids

# Rename cluster if median enrichment score is greater than 0.1  
if(nrow(data.mast) > 0){
  for (i in 1:nrow(data.mast)){
    cluster.ids[[data.mast$group1[i]]] <- "Mast cell"           # Marker Mast cell cluster 
  }
}else{cluster.ids <- cluster.ids} 
cluster.ids                                                     # No apparent changes in cell annotations

# Rename cluster if median enrichment score is greater than 0.1  
if (nrow(data.B) >0 ){
  for (i in 1:nrow(data.B)){
    cluster.ids[[data.B$group1[i]]] <- "B cell"                 # Marker B cell cluster 
  }
}else{cluster.ids <- cluster.ids}
cluster.ids                                                     # No apparent changes in cell annotations

cluster.ids <- as.data.frame(cluster.ids)                       # Convert character object into a dataframe
levels(Idents(rna)) <- cluster.ids$cluster.ids                  # Assign cell types to each cluster previously created
rna$cell.type <- Idents(rna)                                    # New metadata containing cell types
rna$cell.type <- paste0(rna$RNA_snn_res.0.7,"-",rna$cell.type)  # Cell type + cluster
Idents(rna) <- rna$cell.type                                    # Assign cell types and clusters in metadata
DimPlot(rna, label = T, repel = T,group.by = "cell.type")

Idents(rna) <- "cell.type"
Wilcox.markers <- FindAllMarkers(object = rna, min.pct = 0.25,only.pos = T, test.use = "wilcox") # Find differentially expressed genes

saveRDS(rna, "1_Data/Created_data/rna_postinfercnv.rds")
save(Wilcox.markers, file = "1_Data/Created_data/wilcox_markers.RData")