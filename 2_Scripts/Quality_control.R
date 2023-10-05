#----------Load necessary libraries----------
library(Seurat)
library(ggplot2)
library(ggpubr)
library(scater)

#----------Load objects-----------
load("1_Data/seu_ob.RData")

#----------Mitochondrial genes----------
# Add a column with the percentage of reads that correspond to mitochondrial genes
seu_ob[["percent_mt"]] <- PercentageFeatureSet(seu_ob, pattern = "^MT-")
head(seu_ob)

#----------Quality plots--------------
# nFeature_RNA = total number of features that have at least 1 read per cell
# nCount_RNA = total number of reads per cell
# percent_mt = percent of features that are mitochondrial genes

vln_pre <- VlnPlot(seu_ob, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"))
vln_pre # Violin plot shows the features, rna counts, and mitochondrial gene percentage

count_x_feature <- FeatureScatter(object = seu_ob, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
count_x_feature # Scatterplot reads_per_cell**features_per_cell

mt_x_nfeat <- FeatureScatter(object = seu_ob, feature1 = "percent_mt", feature2 = "nFeature_RNA")
mt_x_ncount <- FeatureScatter(object = seu_ob, feature1 = "percent_mt", feature2 = "nCount_RNA")
ncount_x_mfeat <- FeatureScatter(object = seu_ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


