#--------------------Libraries---------------------
# Load necessary libraries for analysis
library(dplyr) 
library(Seurat) 
library(ggplot2) 
library(tidyverse) 
library(infercnv) 
library(psych) 
library(celldex) 
library(SingleR) 
library(SingleCellExperiment) 
library(EnsDb.Hsapiens.v86) 

#--------------------Load previously created seurat object---------------------
readRDS("1_Data/Created_data/rna_postdoublet_SkipChecks.rds")

#--------------------Pre-InferCNV analysis---------------------
rna <- rna[,cells]
rna <- rna[,!(colnames(rna) %in% doublets)]
dim(rna)

#Normalize
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

#Scaling
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna)

feats <- list(stromal,immune,fibroblast,endothelial,epithelial,smooth,plasma)

rna <- AddModuleScore(rna,features = feats,name = c("stromal.","immune.","fibroblast.","endothelial.",
                                                    "epithelial.","smooth.","plasma."),search = T)

# Add PC1 to metadata
rna@meta.data$PC1 <- rna@reductions$pca@cell.embeddings[,1]

count_cor_PC1 <- cor(rna$PC1,rna$nCount_RNA,method = "spearman")

stromal.cor <- cor(rna$PC1,rna$stromal.1,method = "spearman")
immune.cor <- cor(rna$PC1,rna$immune.2,method = "spearman")
fibroblast.cor <- cor(rna$PC1,rna$fibroblast.3,method = "spearman")
endothelial.cor <- cor(rna$PC1,rna$endothelial.4,method = "spearman")
epithelial.cor <- cor(rna$PC1,rna$epithelial.5,method = "spearman")
smooth.cor <- cor(rna$PC1,rna$smooth.6,method = "spearman")
plasma.cor <- cor(rna$PC1,rna$plasma.7,method = "spearman")

# Make JackStraw Plot
rna <- JackStraw(rna, num.replicate = 100,dims = 50)
rna <- ScoreJackStraw(rna,dims = 1:50)
post_doublet_jackstraw <- JackStrawPlot(rna, dims = 1:50)

# Construct a KNN graph based on the euclidean distance in PCA space
# Calculate the neighborhood overlap (Jaccard index) between every cell and its k param nearest neighbors
rna <- FindNeighbors(rna,dims = 1:50)
rna <- FindClusters(rna,resolution = 0.7)
rna <- RunUMAP(rna,dims = 1:50)
Idents(rna) <- "RNA_snn_res.0.7"

# Proceed with CNV
counts_matrix = GetAssayData(rna, layer="counts") # Get expression matrix

# Identify immune clusters
# Find immune cells by relative enrichment of ESTIMATE immune signature
# Plot violinplot with only features for immune cells signals
test <- VlnPlot(rna, features = "immune.2")
# Describe data stored in plot with basic summary statistics per cluster
data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
# Filter which cluster may correspond to immune cells based on median
data.immune <- dplyr::filter(data, median > 0.1)

# Plot violinplot with only features for plasma cells signals
test <- VlnPlot(rna,features = "plasma.7")
# Describe data stored in plot with basic summary statistics per cluster
data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
# Filter which cluster may correspond to immune cells based on median
data.plasma <- dplyr::filter(data, median > 0.1)

# Determine clusters corresponding to plasma and immune cells
immune.clusters <- intersect(data.immune$group1, levels(Idents(rna)))
plasma.clusters <- intersect(data.plasma$group1, levels(Idents(rna)))

# Select only non-repeated clusters 
immune.clusters <- unique(append(immune.clusters, plasma.clusters))

# Assign cluster identification
for (i in 1:length(immune.clusters)){
  j <- which(levels(Idents(rna)) == immune.clusters[i])
  levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
}
rna@meta.data$postdoublet.idents <- Idents(rna)

# Dataframe with cell IDs and cluster IDs
idents <- data.frame(rownames(rna@meta.data), rna@meta.data$postdoublet.idents)

# Only 3 immune clusters were found to be assigned a signature
colnames(idents) <- c("V1","V2")          # Rename columns
saveRDS(rna,"1_Data/Created_data/rna_postdoublet_preinferCNV.rds")

# Make inferCNV inputs
rownames(idents) <- NULL
colnames(idents) <- NULL
# Create table with idents
write.table(idents,"1_Data/Created_data/sample_annotation_file_inferCNV.txt", sep = "\t", row.names = FALSE)
idents <- read.delim("1_Data/Created_data/sample_annotation_file_inferCNV.txt", header = F)

#### Convert original annotation data into table with symbol and chromosome position
gtf <- read.delim(GRCH38.annotations, header = F)
convert.symbol = function(Data){
  ensembls <- Data$V1
  ensembls <- gsub("\\.[0-9]*$", "", ensembls)
  geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembls, keytype = "GENEID", columns = "SYMBOL")
  Data <- rowr::cbind.fill(Data, geneIDs1, fill = NA)
  Data <- na.omit(Data)
  Data$feature <- Data$SYMBOL
  Data.new <- data.frame(Data$SYMBOL,Data$V2,Data$V3,Data$V4)
  Data.new$Data.V2 <- paste("chr",Data.new$Data.V2,sep = "")
  Data.new$Data.SYMBOL <- make.unique(Data.new$Data.SYMBOL)
  return(Data.new)
}
gtf <- convert.symbol(gtf)
head(gtf)

write.table(gtf, "1_Data/Created_data/Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",
            row.names = FALSE,col.names = FALSE)

# create the infercnv object

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                      annotations_file="1_Data/Created_data/sample_annotation_file_inferCNV.txt",
                                      gene_order_file="1_Data/Created_data/Homo_sapiens.GRCh38.86.symbol.txt",
                                      ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                        paste0("immune.",immune.clusters[2]),
                                                        paste0("immune.",immune.clusters[3])))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="1_Data/Created_data/output_dir_CNV_postdoublet_SkipChecks", 
                             cluster_by_groups=T,   # cluster
                             denoise=T,scale_data = T,
                             HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                             BayesMaxPNormal = 0.4,
                             num_threads = 8)

regions <- read.delim("1_Data/Created_data/output_dir_CNV_postdoublet_SkipChecks/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
probs <- read.delim("1_Data/Created_data/output_dir_CNV_postdoublet_SkipChecks/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")

probs <- as.data.frame(t(probs[3,]))
colnames(probs) <- "Prob.Normal"
probs <- dplyr::filter(probs,Prob.Normal < 0.05)
cnvs <- rownames(probs)
cnvs <- gsub("\\.","-",cnvs)

regions <- regions[regions$cnv_name %in% cnvs, ]

cnv.groups <- sub("\\..*", "", regions$cell_group_name)

length(which(rownames(rna@reductions$pca@cell.embeddings) == rownames(rna@meta.data)))
rna$PC1.loading <- rna@reductions$pca@cell.embeddings[,1]
rna$cell.barcode <- rownames(rna@meta.data)
rna$CNV.Pos <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.groups,1,0)

cnv.freq <- data.frame(table(regions$cell_group_name))
cnv.freq$Var1 <- sub("\\..*", "", cnv.freq$Var1)

rna$Total_CNVs <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.freq$Var1,cnv.freq$Freq,0)

boxplot.cnv <- ggplot(rna@meta.data,aes(x= postdoublet.idents,y=PC1.loading,color = as.factor(CNV.Pos)))+geom_boxplot()
boxplot.cnv
ggsave(filename = "Postdoublet_CNV_PC1_boxplot.png", plot = last_plot(), path = "3_Figures/")

saveRDS(rna,"1_Data/Created_data/rna_postdoublet_SkipChecks.rds")

<<<<<<< HEAD
#--------------------Post-InferCN cell type annotation--------------------

=======
#--------------------Post-InferCN analysis--------------------
>>>>>>> 0303a3dc43490856db6df8faa7243731efa8307b
Idents(rna)<- "RNA_snn_res.0.7"

# DEG analysis with Wilcox
# DEG analysis with Wilcoxon
Wilcox.markers <- FindAllMarkers(object =rna, min.pct = 0.25,only.pos = F, test.use = "wilcox")

# Part 4: SingleR cell typing 
# SingleR labeling of celltypes
rna.sce <- as.SingleCellExperiment(rna)

#### Cell type annotation
# Wang et al. GSE111976 single cell, cell type annotation as reference 
# Set reference paths:
ref.data.counts <- readRDS("1_Data/GSE111976_ct_endo_10x.rds")      # Read previously downloaded object
meta <- read.csv("1_Data/GSE111976_summary_10x_day_donor_ctype.csv")# CSV with cell types 
rownames(meta) <- meta$X                                            # Add cell types to metadata

length(which(rownames(meta) == colnames(ref.data.counts)))
ref.data.endo <- CreateSeuratObject(ref.data.counts,meta.data = meta)
Idents(ref.data.endo) <- "cell_type"
ref.data.endo <- NormalizeData(ref.data.endo)

#### SingleR annotation with different datasets
# 1) Wang et al. GSE111976
ref.data.endo <- as.SingleCellExperiment(ref.data.endo) # Preivously created 

# 2) Human Primary Cell Atlas Data (microarray)
ref.data.HPCA <- celldex::HumanPrimaryCellAtlasData()

# # 3) BluePrint Encode (bulk RNA-seq)
ref.data.BED <- celldex::BlueprintEncodeData()

#### Actual cell type annotations
# 2) Single-cell level label transfer: 
predictions.HPCA.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts",
                               ref=ref.data.HPCA, labels=ref.data.HPCA$label.main)
predictions.BED.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts",
                              ref=ref.data.BED, labels=ref.data.BED$label.main)
# Use de.method wilcox with scRNA-seq reference b/c the reference data is more sparse
predictions.endo.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts",
                               ref=ref.data.endo, labels=ref.data.endo$cell_type,de.method = "wilcox")

rna$SingleR.HPCA <- predictions.HPCA.sc$pruned.labels
rna$SingleR.BED <- predictions.BED.sc$pruned.labels
rna$SingleR.endo <- predictions.endo.sc$pruned.labels

saveRDS(rna,paste0("1_Data/Created_data/",SAMPLE.ID,"_scRNA_processed.rds"))

#### Data summary and saving
# # Starting cells, PostQC cells, doublets, Post doublet/QC cells, Cluster #
output.meta <- data.frame(StartingNumCells=length(colnames(matrix_tsv)),
                          nMADLogCounts =2,
                          nMADLogFeatures = 2,
                          nMADLog1pMito =2,
                          PostQCNumCells=PostQCNumCells,
                          ExpectedDoubletFraction=doublet.rate,
                          ObservedDoubletFraction=length(doublets)/length(colnames(matrix_tsv)),
                          PostDoubletNumCells=length(colnames(rna)),
                          NumClusters=length(levels(Idents(rna))),
                          DoubletFinderpK = pK.1,
                          MinNumCounts=min(rna$nCount_RNA),
                          MaxNumCounts= max(rna$nCount_RNA),
                          MedianNumbCounts = median(rna$nCount_RNA),
                          MinNumFeats=min(rna$nFeature_RNA),
                          MaxNumFeats= max(rna$nFeature_RNA),
                          MedianNumbFeats = median(rna$nFeature_RNA),
                          stringsAsFactors=FALSE)
output <- as.data.frame(t(output.meta))
colnames(output) <- SAMPLE.ID
xlsx::write.xlsx(output, "scRNA_pipeline_summary.xlsx",
                 row.names = T, col.names = TRUE)

