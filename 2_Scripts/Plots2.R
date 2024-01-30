#############
# Libraries #
#############
library(ggplot2)
library(Seurat)
library(forcats)
library(RColorBrewer)
library(stringr)

#########################################
#  #
#########################################
load("4_Intermediate/rna.RData")

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


my_levels <- c(11,
               3,
               9,10,
               0,
               6,8,12,
               7,
               1,
               2,4,
               5,13)

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

p1 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))
p1

p2 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)

meta <- rna@meta.data
meta$cluster <- factor(meta$RNA_snn_res.0.7,levels = rev(c(11,
                                                           3,
                                                           9,10,
                                                           0,
                                                           6,8,12,
                                                           7,
                                                           1,
                                                           2,4,
                                                           5,13)))

p0<-ggplot(meta,aes(x=cluster,y=Total_CNVs,fill=cluster))+geom_boxplot(lwd=0.45,outlier.size = 0.55,fatten=0.95)+coord_flip()+
  theme_classic()+scale_fill_manual(values = rev(cols))+NoLegend()+
  scale_y_continuous(position = "right")
p0

# Make violin plots for CA125 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "MUC16",pt.size = 0)+coord_flip()+NoLegend()
p1

p2 <- ggplot(p1$data,aes(y=ident,x=MUC16))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CA125 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")
p2

# Make violin plots for CA125 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "WFDC2",pt.size = 0)+coord_flip()+NoLegend()
p1

p3 <- ggplot(p1$data,aes(y=ident,x=WFDC2))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("HE4 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")
p3

# Make violin plots for CD117 expression
# Relevel object@ident
rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna,features = "KIT",pt.size = 0)+coord_flip()+NoLegend()
p1

p4 <- ggplot(p1$data,aes(y=ident,x=KIT))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.55,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("CD117 expression")+
  scale_fill_manual(values=rev(cols))+
  scale_x_continuous(position = "top")
p4

sampleColors <- RColorBrewer::brewer.pal(11,"Paired")
sampleColors[11] <- "#8c8b8b"
pie(rep(1,11), col=sampleColors) 


# Color patient tumors to resemble the cancer ribbon color 
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],sampleColors[10],sampleColors[9],sampleColors[4],sampleColors[3],sampleColors[2],sampleColors[11],sampleColors[1])

# Patient proportion per subcluster in RNA:
meta <- rna@meta.data

df <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

CombinePlots(list(p4,p3,p2),ncol=3)

p.CNV <- FeaturePlot(rna,features = "Total_CNVs")
p.CA125 <- FeaturePlot(rna,features = "MUC16")
p.HE4 <- FeaturePlot(rna,features = "WFDC2")
p.CD117 <- FeaturePlot(rna,features = "KIT")

CombinePlots(list(p.CA125,p.HE4,p.CD117),ncol=1)

malignant <- c(0,1,2,3,4,6,8,13)
remaining <- c(5,7,9,10,11,12)

# find markers avg logFC >= 0.5 and adj pval <=0.01
for ( i in malignant){
  tryCatch(
    markers <- FindMarkers(rna,ident.1=i,
                           ident.2=remaining,
                           features=c("MUC16",
                                      "WFDC2",
                                      "KIT"),
                           min.pct = 0,only.pos = T,
                           logfc.threshold = 0.5)
  )
  if(nrow(markers) >0){
    print(paste0("Cancer markers present for ",i,":"))
    print(markers)
    print(paste0("End markers for ",i))
  }else{
  }
}


for(i in 0:(nrow(cluster.ids)-1)){
  mark <- FindMarkers(rna,ident.1=i,
              ident.2=remaining,
              features=c("MUC16",
                         "WFDC2",
                         "KIT"),
              min.pct = 0,only.pos = T,
              logfc.threshold = 0.5)
  print(paste0("cluster",i))
  print(mark)
  print("####################")
}

rna$cancer <- " "

for(i in 1:length(rna$seurat_clusters)){
  if(sum(rna$seurat_clusters[i] == malignant) != 0){
    rna$cancer[i] <- paste0("Endometrial cancer")
  }else{
    rna$cancer[i] <- rna$cell.type[i]
  }
}

rna.df$cancer <- rna$cancer

p6 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cancer))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()+ 
  annotate(geom="text", x=3, y=-3.5, label="Cáncer endometrial",color="magenta", size = 6)+
  annotate(geom="text", x=6, y=16, label="Linfocitos",color="red", size = 6)+
  annotate(geom="text", x=-13, y=-3.5, label="Endotelio",color="green", size = 6)+
  annotate(geom="text", x=1, y=-13, label="Fibroblastos",color="deepskyblue3", size = 6)+
  annotate(geom="text", x=1, y=-14, label="estromales",color="deepskyblue3", size = 6)+
  annotate(geom="text", x=-8, y=15, label="Macrófagos",color="brown", size = 6)+
  annotate(geom="text", x=-12, y=-12, label="Células de",color="blue", size = 6)+
  annotate(geom="text", x=-12, y=-13.5, label="músculo liso",color="blue", size = 6)
  
LabelClusters(p6, id="cluster.new",color="black",repel = T,size=8)

