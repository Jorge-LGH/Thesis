library(Seurat)
library(RColorBrewer)
library(ggplot2)

load("4_Intermediate/rna.RData")

cancer.clusters <- c("0-Unciliated epithelia 1","1-Ciliated","2-Unciliated epithelia 1",
                     "3-Unciliated epithelia 1","4-Unciliated epithelia 1","6-Ciliated",
                     "8-Unciliated epithelia 1","13-Unciliated epithelia 1")


# Color RNA UMAP according to main Figure 1:

# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cluster <- rna$RNA_snn_res.0.7
rna.df$cell.type <- rna$cell.type
#Manually annotate 23-cluster as smooth muscle
rna.df$cell.type <- str_replace_all(rna.df$cell.type,"23-Stromal fibroblast","23-Smooth muscle cells")

rna.df$cluster <- as.factor(rna.df$cluster)
rna.df$cell.type <- as.factor(rna.df$cell.type)
levels(rna.df$cluster)
levels(rna.df$cell.type)


# Help sort the cluster numbers:
###############################
epi <- grep("pitheli",levels(rna.df$cell.type))
epi.new <- grep("-Ciliated",levels(rna.df$cell.type))
epi <- c(epi,epi.new)

fibro <- grep("ibro",levels(rna.df$cell.type))

smooth <- grep("mooth",levels(rna.df$cell.type))

endo <- grep("dothel",levels(rna.df$cell.type))

t.nk <- grep("T cell",levels(rna.df$cell.type))
t.nk.new <- grep("Lym",levels(rna.df$cell.type))
t.nk <- c(t.nk,t.nk.new)

mac <- grep("acrophage",levels(rna.df$cell.type))

mast <- grep("Mast",levels(rna.df$cell.type))

b <- grep("B cell",levels(rna.df$cell.type))


cell.types.idx <- c(epi,fibro,smooth,endo,t.nk,mac,mast,b)

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

my_levels <- c(0:13)

# Relevel object@ident
rna.df$cluster.new <- factor(x = rna.df$cluster, levels = my_levels)


epithelial.cols <- colorRampPalette(c("#A0E989", "#245719"))
epithelial.cols <- epithelial.cols(14)
fibro.cols <-colorRampPalette(c("#FABFD2", "#B1339E"))
fibro.cols <- fibro.cols(10)
smooth.cols <- c("#b47fe5","#d8b7e8")
endo.cols <- c("#93CEFF","#4A99FF")
t.cols <- c("gray60","gray20","gray40")
macro.cols <- c("#ff6600","#ff9d5c")
mast.cols <- "gold3"
b.cols <- c("#B22222","#CD5C5C")


cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols)


cells <- colnames(rna[,rna$cell.type %in% cancer.clusters])

p1 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
p1 <- LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)
p1
