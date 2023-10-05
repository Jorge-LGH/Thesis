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
# nFeature_RNA = total number of features that have at least 1 read per cell
# nCount_RNA = total number of reads per cell
# percent_mt = percent of features that are mitochondrial genes

#----------Quality plots--------------
# Violin plot shows the features, rna counts, and mitochondrial gene percentage
vln_pre <- VlnPlot(seu_ob, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"))

# Scatterplot reads_per_cell**features_per_cell
count_x_feature <- FeatureScatter(object = seu_ob, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

# Scatterplot reads_per_cell**mitochondiral_genes_percent
count_x_mt <- FeatureScatter(object = seu_ob, feature1 = "nCount_RNA", feature2 = "percent_mt")

plts1 <- ggarrange(vln_pre, count_x_feature, count_x_mt, ncol = 1, nrow = 3)
plts1
ggsave("3_Figures/pre_quality_plots.png")

#----------Distribution plots----------------
# Reads distribution
reads_hist <- ggplot(as.data.frame(log(seu_ob@meta.data$nCount_RNA)), 
                     aes(x = log(seu_ob@meta.data$nCount_RNA))) + 
  geom_histogram(bins = 30, binwidth=.1, color="darkblue", fill="lightblue") +
  xlab("log Reads per cell") + ylab("Frequency") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, linetype = "solid",
                                 colour = "black"))

# Mitochondrial
mito_hist <- ggplot(as.data.frame(log1p(seu_ob@meta.data$percent_mt)), 
                    aes(x = log1p(seu_ob@meta.data$percent_mt))) + 
  geom_histogram(bins = 30, binwidth=.1, color="darkblue", fill="lightblue") +
  xlab("log1p Mitochondrial percent") + ylab("Frequency") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))

# Features
feats_hist <- ggplot(as.data.frame(log(seu_ob@meta.data$nFeature_RNA)), 
                     aes(x = log(seu_ob@meta.data$nFeature_RNA))) + 
  geom_histogram(bins = 30, binwidth=.1, color="darkblue", fill="lightblue") +
  xlab("log Features per cell") + ylab("Frequency") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))

plts2 <- ggarrange(reads_hist, feats_hist, mito_hist, ncol = 1, nrow = 3)
plts2
ggsave("3_Figures/pre_quality_dist.png")

#----------Filtering---------
# A threshold is established in order to determine if the data is an outlier or not. In this case, 2 nmads
# (median absolute deviation) are used. Number or nmads away from the median required for a value
# to be called an outlier
seu_ob@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(seu_ob@meta.data$nCount_RNA),
                                                      log = F,type = "lower",nmads = 2)

seu_ob@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(seu_ob@meta.data$nFeature_RNA),
                                                        log = F,type = "lower",nmads = 2)

seu_ob@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(seu_ob@meta.data$percent_mt),
                                                      log = F,type = "higher",nmads = 2)

PreQC <- length(colnames(seu_ob)) # Total number of cells before quality control

seu_ob <- subset(seu_ob, subset = nCount_RNA_outlier_2mad == "FALSE" &
                     nFeature_RNA_outlier_2mad == 'FALSE' & 
                     percent_mt_outlier_2mad == "FALSE")

AfterQC <- length(colnames(seu_ob))
PreQC - AfterQC # 929 cells were trimmed based on the selection threshold

#----------Quality plots after filtering-----------
# Violin plot shows the features, rna counts, and mitochondrial gene percentage
vln_pre2 <- VlnPlot(seu_ob, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"))

# Scatterplot reads_per_cell**features_per_cell
count_x_feature2 <- FeatureScatter(object = seu_ob, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")

# Scatterplot reads_per_cell**mitochondiral_genes_percent
count_x_mt2 <- FeatureScatter(object = seu_ob, feature1 = "nCount_RNA", feature2 = "percent_mt")

plts3 <- ggarrange(vln_pre, count_x_feature, count_x_mt, ncol = 1, nrow = 3)
plts3
ggsave("3_Figures/post_quality_plots.png")

#----------Normalization----------
# Normalization methods. The default normalization method is the log normalization, however, since the data is so biased
# towards zero it won't do a great job by default. Another normalization method available is a 
# centered log ratio transformation. We can normalize based on cells or genes, so both methods can be tested
seu_ob <- NormalizeData(seu_ob, normalization.method = "CLR", margin = 2)
seu_ob2<- NormalizeData(seu_ob, normalization.method = "LogNormalize", scale.factor = 10000)
seu_ob3 <- NormalizeData(seu_ob, normalization.method = "CLR", margin = 1)

save(seu_ob, file = "1_Data/seu_ob.RData")
save(seu_ob2, file = "1_Data/seu_ob2.RData")
save(seu_ob3, file = "1_Data/seu_ob3.RData")
