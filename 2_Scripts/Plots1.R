#############
# Libraries #
#############
library(ggplot2)
library(Seurat)
library(forcats)
library(RColorBrewer)

#########################################
# Plotting with cell types and clusters #
#########################################
load("4_Intermediate/rna.RData")

sampleColors <- RColorBrewer::brewer.pal(11,"Paired")                                                               # Assign color codes
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],sampleColors[10],                 # Order color codes
                  sampleColors[9],sampleColors[4],sampleColors[3],sampleColors[2],sampleColors[11],sampleColors[1])
sampleColors[11] <- "#8c8b8b"                                                                                       # Reassign 11th color code

rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)          # Cells' coordinates at UMAP reduction
length(which(rownames(rna.df)==rownames(rna@meta.data)))              # 5124 cells and their UMAP's reduction coordinates
levels(factor(rna$RNA_snn_res.0.7))                                   # 14 levels 0-13

rna.df$cluster <- rna$RNA_snn_res.0.7                                 # Assign cluster to each cell in new data frame object
rna.df$cell.type <- rna$cell.type                                     # Assign cell type to each cell
head(rna.df)

# Manually annotate 23-cluster as smooth muscle
# It will replace a specific cell type to Smooth Muscle Cells
rna.df$cell.type <- str_replace_all(rna.df$cell.type,"23-Stromal fibroblast","23-Smooth muscle cells")
rna.df$cluster <- as.factor(rna.df$cluster)                           # Change from character into factor
rna.df$cell.type <- as.factor(rna.df$cell.type)                       # Change from character into factor
levels(rna.df$cluster)                                                # Levels from 0 to 13
levels(rna.df$cell.type)                                              # Levels from 0 to 13 with cell types

##################################
# Help sort the cluster numbers: #
##################################

epi <- grep("pitheli",levels(rna.df$cell.type))            # Identify clusters with Uniciliated epithelia cell type
epi.new <- grep("-Ciliated",levels(rna.df$cell.type))      # Identify clusters with Ciliated epithelia cell type
epi <- c(epi,epi.new)                                      # Combined

fibro <- grep("ibro",levels(rna.df$cell.type))             # Identify clusters with Stromal fibroblasts cell type
smooth <- grep("mooth",levels(rna.df$cell.type))           # Identify clusters with Smooth Muscle cells cell type
endo <- grep("dothel",levels(rna.df$cell.type))            # Identify clusters with Endothelia cell type
t.nk <- grep("T cell",levels(rna.df$cell.type))            # Identify clusters with t cells and nk cells cell types. NO CELLS DETECTED
t.nk.new <- grep("Lym",levels(rna.df$cell.type))           # Identify clusters with Lymphocytes cell type
t.nk <- c(t.nk,t.nk.new)

mac <- grep("acrophage",levels(rna.df$cell.type))          # Identify clusters with Macrophages cell type
mast <- grep("Mast",levels(rna.df$cell.type))              # Identify clusters with Mast cells cell types. NO CELLS DETECTED
b <- grep("B cell",levels(rna.df$cell.type))               # Identify clusters with B cells cell type. NO CELLS DETECTED
cell.types.idx <- c(epi,fibro,smooth,endo,t.nk,mac,mast,b) # Store clusters by cell type
store <- numeric(0)

for(i in 1:length(cell.types.idx)){
  name <- levels(rna.df$cell.type)[cell.types.idx[i]]
  print(gsub("-.*","",name))
  new.name <- gsub("-.*","",name)
  new.num <- as.numeric(new.name)
  store[i] <- new.num
}
print(store)

#####################################################
# Establish levels for existing clusters, only from 0-13
my_levels <- c(11,
               3,
               9,10,
               0,
               6,8,12,
               7,
               1,
               2,4,
               5,13)

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

p1 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)

p2 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
LabelClusters(p2,id="cell.type",color="black",repel = T,size=4)

meta <- rna@meta.data
df <- meta %>% group_by(RNA_snn_res.0.7) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c(11,
                                                      3,
                                                      9,10,
                                                      0,
                                                      6,8,12,
                                                      7,
                                                      1,
                                                      2,4,
                                                      5,13))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + geom_bar(stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+scale_fill_manual(values = sampleColors)

total.rna <- as.data.frame(table(rna$cell.type))
colnames(total.rna) <- c("Cluster","RNA cells")

