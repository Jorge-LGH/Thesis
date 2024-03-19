#--------------------Libraries---------------------
# Load necessary libraries for analysis
library(scater) 
library(dplyr) 
library(Seurat) 
library(ggplot2) 
library(tidyverse) 

#--------------------Datasets--------------------
# Load necessaary data for analysis
# These datasets are all previously downloaded from different sources

#### Panglao dataset
# The panglao datasets contains genes for cell type identification
tsv <- gzfile("1_Data/PanglaoDB_markers_27_Mar_2020.tsv.gz") # This is a database with cell type markers
                                                             # it will be needed in order to assign the 
                                                             #putative cell types to the cells

panglaodb <- read.csv(tsv, header = T, sep = "\t")               # make tsv file for managing the data
panglaodb <- dplyr::filter(panglaodb, 
                           species == "Hs" | species == "Mm Hs") # subsetting only for human cell type markers

panglaodb <- dplyr::filter(panglaodb,                            # selecting only the desired cell type markers
                             organ == "Connective tissue" |        
                             organ == "Epithelium" |
                             organ == "Immune system" |
                             organ == "Reproductive"|
                             organ == "Vasculature" |
                             organ == "Smooth muscle")
panglaodb <- split(as.character(panglaodb$official.gene.symbol), # only keep the official gene
                   panglaodb$cell.type)                          # symbols and to which cell type they belong to 
                                                                 # in order to reduce the data


#### Loading Barcodes
# barcodes serving as identifiers for each cell
barcodes_tsv <- read.table("1_Data/GSM5276935_barcodes-36186L.tsv")

#### Loading matrix RNA counts
# Counts for each gene per cell gene x cell
# 33538 features (genes) and 6054 cells
matrix_tsv <- Matrix::readMM("1_Data/GSM5276935_matrix-36186L.mtx")

#### Loading features database
# This is the table with all the genes and their identifiers
features_tsv <- read.table("1_Data/GSM5276935_features-36186L.tsv")

#### Combine datasets
# This is done for the creation of the Seurat object that will be manipulated down the line
colnames(matrix_tsv) <- barcodes_tsv[[1]] # add column names, they represent the cell barcode
rownames(matrix_tsv) <- features_tsv[ ,2] # add row names, they represent the features or the
                                          # gene identifiers

#### Load estimate signatures
ESTIMATE.signatures <- "1_Data/ESTIMATE_signatures.csv" # CSV with genes associated to stromal or
                                                        # immune cell types

#### Load homo sapiens gene annotation
GRCH38.annotations <- "1_Data/Homo_sapiens.GRCh38.86.txt" # Dataset with gene ID's and chromosome
                                                          # coordinates

doublet.rate <- 0.0460     # Doublet rate appearance estimation

SAMPLE.ID <- "endo_36186L" # Identify sample, this case endometrial tumor

#--------------------scRNA-seq processing before doublet removal--------------------
# Basic analysis perfomred previously to doublet removal
# The steps include seurat object creation, quality control, feature selection,
# dimension reduction, immune signatures, indetifying clusters, and another reduction

#### Create Seurat object
rna <- CreateSeuratObject(matrix_tsv,         # reads matrix
                          project = "Endo_3", # Name of the project, could have more data
                          assay = "RNA",      # Name of the initial assay
                          min.cells = 3,      # min of cells that a feature must be present in
                          min.features = 1)   # min of features a cell must have to be kept

PreQCNumCells <- length(colnames(rna))        # 6054 cells

#### Quality control
# Add metadata with percentage of all the counts corresponding to mitochondrial genes
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

# Add metadata with percentage of all the counts corresponding to ribosomal genes
rna[["percent.Ribosomal"]] <- PercentageFeatureSet(rna, pattern = "^RP[LS]")

#### Quality Plots
# Violin plot with read count, gene count per cell, mitochondrial gene reads, and ribosomal gene counts
vln_pre_quality <- VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", 
                                             "percent.mt", "percent.Ribosomal"), ncol = 4)

# Scatter plot read count per cell by feature count per cell
# 0.88 correlation between both axes 
scatter_pre_quality_1 <- FeatureScatter(rna, "nCount_RNA", "nFeature_RNA") + theme(legend.position="none") +
  xlab("Reads per cell") + ylab("Features per cell") 

# Scatter plot read count per cell by mitochondrial gene reads per cell
# -0.24 correlation between both axes 
scatter_pre_quality_2 <- FeatureScatter(rna, "nCount_RNA", "percent.mt") + theme(legend.position="none") + 
  xlab("Reads per cell") + ylab("% of mitochondrial genes")

scatters_pre_quality <- scatter_pre_quality_1 + scatter_pre_quality_2

# Histograms to see percent of mitochondrial and ribosomal genes
qc.metrics <- as_tibble(rna[[]],rownames="Cell.Barcode") # Create tibble for easier management

mito_percent_distri <- qc.metrics %>% ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  xlab("Mitochondrial genes percent") + ylab("Cells")+ 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"))

ribo_percent_distri <- qc.metrics %>% ggplot(aes(percent.Ribosomal)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  xlab("Ribosomal genes percent") + ylab("Cells") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"))

mito_ribo_pre_distri <- mito_percent_distri + ribo_percent_distri

# Gene expression distribution pre-normalization and quality thresholding
pre_normalization_expression <- as_tibble(rna@assays$RNA@data[,1:50]) %>% 
  pivot_longer(cols=everything(), names_to="cell", values_to="expression") %>%
  ggplot(aes(x=cell, y=expression)) + xlab("Cell ID") + ylab("Gene expression") +
  stat_summary(geom="crossbar", fun.data=mean_sdl) + 
  ggtitle("Pre-Normalization Gene Expression Distribution") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 0.5,linetype = "solid", colour = "black"))

#### Quality control thresholding
# Will use Median absolute deviations (MADs) for thresholding 
# nCount_RNA, nFeature_RNA, and percent mitochrondial counts will be thresholded
# Outliers are >2 MADs

rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna@meta.data$nCount_RNA),
                                                   log = F,type = "lower",nmads = 2)

rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna@meta.data$nFeature_RNA),
                                                     log = F,type = "lower",nmads = 2)

rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna@meta.data$percent.mt),
                                                   log = F,type = "higher",nmads = 2)

# Subsetting
rna <- subset(rna, subset = nCount_RNA_outlier_2mad == "FALSE" &
                nFeature_RNA_outlier_2mad == 'FALSE' & 
                percent_mt_outlier_2mad == "FALSE")

PostQCNumCells <- length(colnames(rna)) # 5125 cells remaining
PreQCNumCells - PostQCNumCells          # 929 cells are lost

#### Post quality thresholding plots
vln_post_quality <- VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", 
                                              "percent.mt", "percent.Ribosomal"), ncol = 4)

scatter_post_quality_1 <- FeatureScatter(rna, "nCount_RNA", "nFeature_RNA") + theme(legend.position="none") +
  xlab("Reads per cell") + ylab("Features per cell") 

scatter_post_quality_2 <- FeatureScatter(rna, "nCount_RNA", "percent.mt") + theme(legend.position="none") + 
  xlab("Reads per cell") + ylab("% of mitochondrial genes")

scatters_post_quality <- scatter_post_quality_1 + scatter_post_quality_2

qc.metrics <- as_tibble(rna[[]],rownames="Cell.Barcode") 

mito_percent_distri_post_quality <- qc.metrics %>% ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  xlab("Mitochondrial genes percent") + ylab("Cells")+ 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"))

ribo_percent_distri_post_quality  <- qc.metrics %>% ggplot(aes(percent.Ribosomal)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  xlab("Ribosomal genes percent") + ylab("Cells") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"))

mito_ribo_pre_distri <- mito_percent_distri + ribo_percent_distri

#### Data normalization 
# Log normalization, still need to get why this method and what does it do specifically
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

# Gene expression distribution post-normalization and quality thresholding
post_normalization_expression <- as_tibble(rna@assays$RNA@data[,1:50]) %>% 
  pivot_longer(cols=everything(), names_to="cell", values_to="expression") %>%
  ggplot(aes(x=cell, y=expression)) + xlab("Cell ID") + ylab("Gene expression") +
  stat_summary(geom="crossbar", fun.data=mean_sdl) + 
  ggtitle("Pre-Normalization Gene Expression Distribution") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 0.5,linetype = "solid", colour = "black"))

#### Feature selection
# Find most variable genes based on the "vst" selection method
# Need more information about selection method
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

# Plot variable genes
top10 <- head(VariableFeatures(rna), 10)
var_genes_plt <- VariableFeaturePlot(rna)
varg_genes_plot_pre <- LabelPoints(var_genes_plt, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

#### Scaling data and PCA dimensionality reduction
# Scale data:  each feature will be centered to have a mean of 0 
# and scaled by the standard deviation of each feature
# PCA: Principal component analysis dimensionality reduction

all.genes <- rownames(rna)                  # Character of all gene ID's
rna <- ScaleData(rna, features = all.genes) # Scales and centers features in the dataset
rna <- RunPCA(rna)

#### Score cells for immune/stromal/fibroblast/endothelial signatures
# Read immune and stromal cell type signatures
immune.stromal <- read.csv(ESTIMATE.signatures,header = F) 

# signal subset for immune cells
stromal <- immune.stromal$V1[1:141]
immune <- immune.stromal$V1[142:282]
fibroblast <- panglaodb$Fibroblasts
endothelial <- panglaodb$`Endothelial cells`
epithelial <- panglaodb$`Epithelial cells`
smooth <- panglaodb$`Smooth muscle cells`
plasma <- panglaodb$`Plasma cells`

# Combine all cell types' genes into a list
feats <- list(stromal,immune,fibroblast,endothelial,epithelial,smooth,plasma)

# Adding module score:  It is the relative expression, the higher the gene expression in a cell compared to 
# the others and control feature stes, the higher the score
rna <- AddModuleScore(rna, features = feats, name = c("stromal.","immune.","fibroblast.","endothelial.",
                                                      "epithelial.","smooth.","plasma."), search = T)

#### PC significance and analysis
# PC1 variance explained for each cell
rna@meta.data$PC1 <- rna@reductions$pca@cell.embeddings[,1]

# Spearman correlation between PC1 and the reads per gene per cell
count_cor_PC1 <- cor(rna$PC1, rna$nCount_RNA, method = "spearman") # 0.251827

# Spearman correlation between 
stromal.cor     <- cor(rna$PC1, rna$stromal.1, method = "spearman")     # -0.7573302
immune.cor      <- cor(rna$PC1, rna$immune.2, method = "spearman")      # -0.6883634
fibroblast.cor  <- cor(rna$PC1, rna$fibroblast.3, method = "spearman")  # -0.7558407
endothelial.cor <- cor(rna$PC1, rna$endothelial.4, method = "spearman") # -0.5271504
epithelial.cor  <- cor(rna$PC1, rna$epithelial.5, method = "spearman")  # 0.3915784
smooth.cor      <- cor(rna$PC1, rna$smooth.6, method = "spearman")      # -0.5014221
plasma.cor      <- cor(rna$PC1, rna$plasma.7, method = "spearman")      # 0.02972904

# Jackstraw: Identification of statistically significant genomic variables associated
# with subset or linear combination of PCs
rna <- JackStraw(rna, num.replicate = 100, dims = 50) # 100 random subset sampling and 50 PCs
rna <- ScoreJackStraw(rna, dims = 1:50)
pre_doublet_jackstraw <- JackStrawPlot(rna, dims = 1:50)

# Construct a KNN graph based on the euclidean distance in PCA space
# Calculate the neighborhood overlap (Jaccard index) between every cell and its k param nearest neighbors
rna <- FindNeighbors(rna, dims = 1:50, reduction = "pca")

# Identify clusters by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
# First calculate k-nearest neighbors and construct the SNN graph
# Clusters are stores in metadata (13 clusters appear at this moment)
rna <- FindClusters(rna, resolution = 0.5)

# Run UMAP to PCs previously performed
rna <- RunUMAP(rna, dims = 1:50)

# Change rna identifiers
Idents(rna) <- "RNA_snn_res.0.7" 

# Save object
saveRDS(rna,"1_Data/Created_data/rna_predoublet_SkipChecks.rds")

