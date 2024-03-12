#--------------------Libraries---------------------
# Load necessary libraries for analysis
library(dplyr) 
library(Seurat) 
library(ggplot2) 
library(tidyverse)
library(psych)
library(SingleR)

#--------------------Load previously created seurat object---------------------
rna <- readRDS("1_Data/Created_data/endo_36186L_scRNA_processed.rds")

#--------------------Further cell type annotation---------------------
# SingleR predicted cell types 
rna$SingleR <- rna$SingleR.endo  # Previously annotated cell types
rna$Sample <- "3533EL"           # Sample ID 

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
                                          c("TPSB2","TPSAB1","KIT")),     # Three gene Mast signature
                      name = c("B.","Plasma.","Mast.","Macrophage.","DC.",
                               "T.","NK.","Endothelial.","Fibroblast.","Epithelial.","Smooth_muscle.","Mast_3_gene."),search = T)

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

#### Aggregate all cells pertaining to their specific clusters 

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

cell_umap <- DimPlot(rna, label = T, repel = T,group.by = "cell.type") + xlab("UMAP 1") +ylab("UMAP 2")
cell_umap

Idents(rna) <- "cell.type"                                      # Make identification == cell type
Wilcox.markers <- FindAllMarkers(object = rna, min.pct = 0.25,only.pos = T, test.use = "wilcox") # Find differentially expressed genes

saveRDS(rna, "1_Data/Created_data/rna_postinfercnv.rds")
save(Wilcox.markers, file = "1_Data/Created_data/wilcox_markers.RData")

#--------------------UMAP plots and reads per cluster plot-------------------
# Establish levels for existing clusters, only from 0-13
my_levels <- c(11,3,9,10,0,6,8,12,7,1,2,4,5,13)

# Create maneable dataframe
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings) # Coordinates for each cell in UMAP visualization
rna.df$cluster <- rna$RNA_snn_res.0.7                        # Assign cluster to each cell in new data frame object
rna.df$cell.type <- rna$cell.type                            # Assign cell type to each cell
head(rna.df)

rna.df$cluster <- as.factor(rna.df$cluster)                  # Change from character into factor
rna.df$cell.type <- as.factor(rna.df$cell.type)              # Change from character into factor
levels(rna.df$cluster)                                       # Levels from 0 to 13
levels(rna.df$cell.type)                                     # Levels from 0 to 13 with cell types

# Relevel object@ident
rna.df$cluster.new <- factor(x = rna.df$cluster, levels = my_levels)

epithelial.cols <- colorRampPalette(c("#A0E989", "#245719")) # Set colors for epithelial cells
epithelial.cols <- epithelial.cols(14)
fibro.cols <-colorRampPalette(c("#FABFD2", "#B1339E"))       # Set color for fibroblasts
fibro.cols <- fibro.cols(10)
smooth.cols <- c("#b47fe5","#d8b7e8")                        # Set colors for smooth muscle cells
endo.cols <- c("#93CEFF","#4A99FF")                          # Set colors for endothelial cells
t.cols <- c("gray60","gray20","gray40")                      # Set colors for T-cells
macro.cols <- c("#ff6600","#ff9d5c")                         # Set colors for macrophages
mast.cols <- "gold3"                                         # Set color for mast cells
b.cols <- c("#B22222","#CD5C5C")                             # Set colors for B-cells

cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols) # Concatenate color tags

plot_umap_cols <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))

umap_no_labels <- LabelClusters(plot_umap_cols,id="cluster.new",color="black",repel = T,size=8)   # UMAP plot without cell id's
umap_with_labels <- LabelClusters(plot_umap_cols,id="cell.type",color="black",repel = T,size=4)   # UMAP plot with cell id's 

# Reads per cluster plot based on cell types 
meta <- rna@meta.data
df <- meta %>% group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c(11,3,9,10,0,6,8,12,7,1,2,4,5,13))) %>% 
  ggplot(aes(fill=Cluster, y=Cells, x= fct_rev(cell.type))) + geom_bar(stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("Cells per cluster")

#--------------------Cancer Cell detection-------------------
meta <- rna@meta.data
meta$cluster <- factor(meta$RNA_snn_res.0.7, levels = rev(my_levels))

#### Violin plots
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))

# MUC16 expression
muc_p1 <- VlnPlot(rna,features = "MUC16",pt.size = 0)+coord_flip()+NoLegend()
muc_p1

# CA125 protein expression by MUC16 gene
p2 <- ggplot(muc_p1$data,aes(y=ident,x=MUC16))+
  geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CA125 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")
p2

# WFDC2 gene expression
wf_p1 <- VlnPlot(rna,features = "WFDC2",pt.size = 0)+coord_flip()+NoLegend()
wf_p1

# HE4 protein expression by WFDC2
p3 <- ggplot(wf_p1$data,aes(y=ident,x=WFDC2))+
  geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("HE4 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")
p3

# Make violin plots for KIT expression
c117_p1 <- VlnPlot(rna,features = "KIT",pt.size = 0)+coord_flip()+NoLegend()
c117_p1

# CD117 protein 
p4 <- ggplot(c117_p1$data,aes(y=ident,x=KIT))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CD117 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")
p4

# Previous plots combined
p.CA125 <- FeaturePlot(rna,features = "MUC16")
p.HE4 <- FeaturePlot(rna,features = "WFDC2")
p.CD117 <- FeaturePlot(rna,features = "KIT") # Kinda weird it is so low
                                             # Highest reads goes to 5 and only present in 63 cells

exp_plot <- CombinePlots(list(p4,p3,p2),ncol=3)
expression_combined <- CombinePlots(list(umap_no_labels, p.CA125, p.HE4, p.CD117),ncol=2)

#### Detect cancer clusters
# Clusters selected based on visual analysis of marker genes' expression
malignant <- c(0,1,2,3,4,7,8,13)
remaining <- c(5,6,9,10,11,12)

# Get marker gene information on malignant clusters
for ( i in malignant){
  tryCatch(
    markers <- FindMarkers(rna,ident.1=i,
                           ident.2=remaining,
                           features=c("MUC16","WFDC2","KIT"),
                           min.pct = 0,only.pos = T,
                           logfc.threshold = 0.5))
  if(nrow(markers) >0){
    print(paste0("Cancer markers present for ",i,":"))
    print(markers)
    print(paste0("End markers for ",i))
  }
}

# Make empty column in the object's metadata
rna$cancer <- " "

# Assign label corresponding to cancer or maintain cell type annotation
for(i in 1:length(rna$seurat_clusters)){
  if(sum(rna$seurat_clusters[i] == malignant) != 0){
    rna$cancer[i] <- paste0("Endometrial cancer")
  }else{
    rna$cancer[i] <- rna$cell.type[i]
  }
}

# Plot with endometrial cancer clusters assigned as well as other cell types
p6 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cancer))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()+ 
  annotate(geom="text", x=1, y=-3.5, label="Cáncer endometrial",color="magenta", size = 6)+
  annotate(geom="text", x=12, y=14, label="Linfocitos",color="red", size = 6)+
  annotate(geom="text", x=13, y=2, label="Endotelio",color="green", size = 6)+
  annotate(geom="text", x=-3, y=-12, label="Fibroblastos",color="deepskyblue3", size = 6)+
  annotate(geom="text", x=-3, y=-13, label="estromales",color="deepskyblue3", size = 6)+
  annotate(geom="text", x=-8, y=15, label="Macrófagos",color="brown", size = 6)+
  annotate(geom="text", x=6, y=-9, label="Células de",color="blue", size = 6)+
  annotate(geom="text", x=6, y=-10, label="músculo liso",color="blue", size = 6)

final_umap <- LabelClusters(p6, id="cluster.new",color="black",repel = T,size=8)
final_umap

saveRDS(rna, "1_Data/Created_data/pre_network_rna.rds")
