# Libraries 
library(scater)
library(dplyr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(infercnv)
library(stringr)
library(psych)
library(celldex)
library(SingleR)
library(scran)
library(patchwork)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(msigdbr)
library(fgsea)
library(tibble)
library(DoubletFinder)
library(Signac)

#----------------------------------------------------#
#-------------------Datasets-------------------------
#----------------------------------------------------#

tsv <- gzfile("1_Data/PanglaoDB_markers_27_Mar_2020.tsv.gz") # This is a database with cell type markers
# it will be needed in order to assign the putative cell types to the cells
panglaodb <- read.csv(tsv, header = T, sep = "\t") # make tsv file for managing the data
panglaodb <- dplyr::filter(panglaodb, species == "Hs" | species == "Mm Hs") # subsetting only for human cell type markers
panglaodb <- dplyr::filter(panglaodb, organ == "Connective tissue" | # selecting only the desired cell type markers
                             organ == "Epithelium" |
                             organ == "Immune system" |
                             organ == "Reproductive"|
                             organ == "Vasculature" |
                             organ == "Smooth muscle")
panglaodb <- split(as.character(panglaodb$official.gene.symbol), panglaodb$cell.type) # only keep the official gene
# gene symbols and to which cell type they belong to in order to reduce the data


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
# identifiers for the cells in this case
rownames(matrix_tsv) <- features_tsv[ ,2] # add row names, they represent the features or the
# gene identifiers

#### Load estimate signatures
ESTIMATE.signatures <- "1_Data/ESTIMATE_signatures.csv"

#### Load homo sapiens gene annotation
GRCH38.annotations <- "1_Data/Homo_sapiens.GRCh38.86.txt"

doublet.rate <- 0.0460

SAMPLE.ID <- "endo_36186L"

#--------------------------------------------------------------------------#
#-------------------Part 1: scRNA-seq processing before doublet removal----
#--------------------------------------------------------------------------#

#### Basic structure
# Create the Seurat object
# The idea for making this object is it's easier and more straightforward manipulation. It will contain the reads,
# the cell barcode identities and the genes all the samples have. It is created based on the integration of the previous
# data that was downloaded and set into the environment. However, it doesn't contemplate the gene identifiers for the cell
# types that may exist in the samples
rna <- CreateSeuratObject(matrix_tsv,
                          project = "Endo_3", # Name of the project, could have more data
                          assay = "RNA",      # Name of the initial assay
                          min.cells = 3,      # min of cells that a feature must be present in
                          min.features = 1)   # min of features a cell must have to be kept

PreQCNumCells <- length(colnames(rna)) # 6054 cells

#### Quality control
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
# QC metrics: nCount_RNA, nFeature_RNA, and percent mitochrondial counts
# Outliers are >2 MADs
rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna@meta.data$nCount_RNA),log = F,type = "lower",nmads = 2)
rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna@meta.data$nFeature_RNA),log = F,type = "lower",nmads = 2)
rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna@meta.data$percent.mt),log = F,type = "higher",nmads = 2)

rna <- subset(rna, subset = nCount_RNA_outlier_2mad == "FALSE" &
                nFeature_RNA_outlier_2mad == 'FALSE' & 
                percent_mt_outlier_2mad == "FALSE")

PostQCNumCells <- length(colnames(rna))

PreQCNumCells - PostQCNumCells # 929 cells are lost, total of 5125 cells remaining

# Default Seurat processing
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

#Scaling & PCA
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna)

# Score cells for immune/stromal/fibroblast/endothelial signatures
############################################################################################

# Read immune signatures available at panglaodb 
immune.stromal <- read.csv(ESTIMATE.signatures,header = F)

# signal subset for immune cells
stromal <- immune.stromal$V1[1:141]
immune <- immune.stromal$V1[142:282]
fibroblast <- panglaodb$Fibroblasts
endothelial <- panglaodb$`Endothelial cells`
epithelial <- panglaodb$`Epithelial cells`
smooth <- panglaodb$`Smooth muscle cells`
plasma <- panglaodb$`Plasma cells`

feats <- list(stromal,immune,fibroblast,endothelial,epithelial,smooth,plasma)

# add to immune signatures to cells
rna <- AddModuleScore(rna, features = feats, name = c("stromal.","immune.","fibroblast.","endothelial.",
                                                      "epithelial.","smooth.","plasma."), search = T)
#########################################################################

#######################################################################

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
JackStrawPlot(rna, dims = 1:50)

# If PC1 is correlated with read depth, check to see if biological variation is correlated to PC1
if (round(abs(count_cor_PC1),2) > 0.5){
  if( round(abs(stromal.cor),2) >= 0.5 |
      round(abs(immune.cor),2) >= 0.5 |
      round(abs(fibroblast.cor),2) >= 0.5 |
      round(abs(endothelial.cor),2) >= 0.5 |
      round(abs(epithelial.cor),2) >= 0.5 |
      round(abs(smooth.cor),2) >= 0.5 |
      round(abs(plasma.cor),2) >= 0.5){
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    dents(rna) <- "RNA_snn_res.0.7"
    
    # Verify with inferCNV: is PC1 correlated with CNV events/Malignancy?
    #########################################################################
    # inferCNV: is PC1 also correlated with CNV/malignancy status?
    counts_matrix <- GetAssayData(rna, layer="counts")
    
    # Identify immune clusters
    #######################################################
    # Find immune cells by relative enrichment of ESTIMATE immune signature
    test <- VlnPlot(rna,features = "immune.2")
    data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
    data.immune <- dplyr::filter(data,median > 0.1)
    
    test <- VlnPlot(rna,features = "plasma.7")
    data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
    data.plasma <- dplyr::filter(data,median > 0.1)
    
    immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
    plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
    
    immune.clusters <- unique(append(immune.clusters,plasma.clusters))
    
    for (i in 1:length(immune.clusters)){
      j <- which(levels(Idents(rna)) == immune.clusters[i])
      levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
    }
    rna@meta.data$predoublet.idents <- Idents(rna)
    idents <- data.frame(rownames(rna@meta.data),rna@meta.data$predoublet.idents)
    
    colnames(idents) <- c("V1","V2")
    
    saveRDS(rna,"1_Data/rna_predoublet_preinferCNV.rds")
    # Make inferCNV input files
    rownames(idents) <- NULL
    colnames(idents) <- NULL
    write.table(idents,"1_Data/sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
    idents <- read.delim("1_Data/sample_annotation_file_inferCNV.txt",header = F)
    
    
    gtf <- read.delim(GRCH38.annotations,header = F)
    library(EnsDb.Hsapiens.v86)
    convert.symbol = function(Data){
      ensembls <- Data$V1
      ensembls <- gsub("\\.[0-9]*$", "", ensembls)
      geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembls, keytype = "GENEID", columns = "SYMBOL")
      Data <- cbind.fill(Data, geneIDs1, fill = NA)
      Data <- na.omit(Data)
      Data$feature <- Data$SYMBOL
      Data.new <- data.frame(Data$SYMBOL,Data$V2,Data$V3,Data$V4)
      Data.new$Data.V2 <- paste("chr",Data.new$Data.V2,sep = "")
      Data.new$Data.SYMBOL <- make.unique(Data.new$Data.SYMBOL)
      return(Data.new)
    }
    gtf <- convert.symbol(gtf)
    head(gtf)
    
    write.table(gtf,"1_Data/Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
    
    
    num.immune.clusters = length(immune.clusters)
    # create the infercnv object
    if ( num.immune.clusters == 1) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=paste0("immune.",immune.clusters[1]))
      
    } else if (num.immune.clusters == 2) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2])))
      
    } else if ( num.immune.clusters == 3 ) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2]),
                                                             paste0("immune.",immune.clusters[3])))
    } else if (num.immune.clusters == 4) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2]),
                                                             paste0("immune.",immune.clusters[3]),
                                                             paste0("immune.",immune.clusters[4])))
    } else if (num.immune.clusters == 5) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2]),
                                                             paste0("immune.",immune.clusters[3]),
                                                             paste0("immune.",immune.clusters[4]),
                                                             paste0("immune.",immune.clusters[5])))
    } else if (num.immune.clusters == 6) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2]),
                                                             paste0("immune.",immune.clusters[3]),
                                                             paste0("immune.",immune.clusters[4]),
                                                             paste0("immune.",immune.clusters[5]),
                                                             paste0("immune.",immune.clusters[6])))
    }else if (num.immune.clusters == 7) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2]),
                                                             paste0("immune.",immune.clusters[3]),
                                                             paste0("immune.",immune.clusters[4]),
                                                             paste0("immune.",immune.clusters[5]),
                                                             paste0("immune.",immune.clusters[6]),
                                                             paste0("immune.",immune.clusters[7])))
    }else if (num.immune.clusters == 8) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2]),
                                                             paste0("immune.",immune.clusters[3]),
                                                             paste0("immune.",immune.clusters[4]),
                                                             paste0("immune.",immune.clusters[5]),
                                                             paste0("immune.",immune.clusters[6]),
                                                             paste0("immune.",immune.clusters[7]),
                                                             paste0("immune.",immune.clusters[8])))
    }else if (num.immune.clusters == 9) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2]),
                                                             paste0("immune.",immune.clusters[3]),
                                                             paste0("immune.",immune.clusters[4]),
                                                             paste0("immune.",immune.clusters[5]),
                                                             paste0("immune.",immune.clusters[6]),
                                                             paste0("immune.",immune.clusters[7]),
                                                             paste0("immune.",immune.clusters[8]),
                                                             paste0("immune.",immune.clusters[9])))
    }else if (num.immune.clusters == 10) {
      infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                           annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                           gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                           ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                             paste0("immune.",immune.clusters[2]),
                                                             paste0("immune.",immune.clusters[3]),
                                                             paste0("immune.",immune.clusters[4]),
                                                             paste0("immune.",immune.clusters[5]),
                                                             paste0("immune.",immune.clusters[6]),
                                                             paste0("immune.",immune.clusters[7]),
                                                             paste0("immune.",immune.clusters[8]),
                                                             paste0("immune.",immune.clusters[9]),
                                                             paste0("immune.",immune.clusters[10])))
    }
    
    # perform infercnv operations to reveal cnv signal
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                  out_dir="1_Data/output_dir_CNV_predoublet",  # dir is auto-created for storing outputs
                                  cluster_by_groups=T,   # cluster
                                  denoise=T,scale_data = T,
                                  HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                                  BayesMaxPNormal = 0.4, num_threads = 8
                                  
    )
    
    
    
    regions <- read.delim("1_Data/output_dir_CNV_predoublet/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
    probs <- read.delim("1_Data/output_dir_CNV_predoublet/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
    
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
    rna$CNV.Pos <- ifelse(as.character(rna$predoublet.idents) %in% cnv.groups,1,0)
    
    
    cnv.freq <- data.frame(table(regions$cell_group_name))
    cnv.freq$Var1 <- sub("\\..*", "", cnv.freq$Var1)
    
    rna$Total_CNVs <- ifelse(as.character(rna$predoublet.idents) %in% cnv.freq$Var1,cnv.freq$Freq,0)
    
    boxplot.cnv <- ggplot(rna@meta.data,aes(x= predoublet.idents,y=PC1.loading,color = as.factor(CNV.Pos)))+geom_boxplot()
    boxplot.cnv
    ggsave(filename = "Predoublet_CNV_PC1_boxplot.png", plot = last_plot(), path = "1_Data/")
    
    data <- describeBy(boxplot.cnv$data$PC1.loading, boxplot.cnv$data$predoublet.idents, mat = TRUE)
    data$CNV <- ifelse(data$group1 %in% cnv.groups,1,0)
    
    wilcox <- wilcox.test(data = rna@meta.data,PC1.loading~CNV.Pos)
    
    if (wilcox$p.value < 0.05){
      rna <- rna
      library(stringr)
      levels(Idents(rna)) <- str_remove(levels(Idents(rna)),"immune.")
      saveRDS(rna,"1_Data/rna_predoublet_PassedPC1Checks.rds")
    }else{
      all.genes <- rownames(rna)
      rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
      rna <- FindNeighbors(rna,dims = 1:50)
      rna <- FindClusters(rna,resolution = 0.7)
      rna <- RunUMAP(rna,dims = 1:50)
      Idents(rna) <- "RNA_snn_res.0.7"
      
      saveRDS(rna,"1_Data/rna_predoublet_FailedCNVTest.rds")
    }
    
  }else{
    all.genes <- rownames(rna)
    rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    saveRDS(rna,"1_Data/rna_predoublet_FailedCorTest.rds")
  }
}else{
  rna <- FindNeighbors(rna,dims = 1:50)
  rna <- FindClusters(rna,resolution = 0.7)
  rna <- RunUMAP(rna,dims = 1:50)
  Idents(rna) <- "RNA_snn_res.0.7"
  saveRDS(rna,"1_Data/rna_predoublet_SkipChecks.rds")
}

###########################################################
# Part 2: Doublet detection and removal
###########################################################

#' # Doublet detection: 1) DoubletDecon 2) DoubletFinder 3) Take intersection of calls
######################################################################################

# # Doublet Decon
# # Add if else statement to regress out nCount RNA if needed before running DoubletDecon
library(DoubletDecon)
idents.length <- length(levels(Idents(rna)))
for (i in levels(Idents(rna))){
  i <- as.numeric(i)
  levels(Idents(rna))[i] <- i -1
}

#Improved_Seurat_Pre_Process()
#Idents(rna) <- as.factor(Idents(rna))
seuratObject <- rna
newFiles <- Improved_Seurat_Pre_Process(seuratObject, num_genes=50, write_files=F)

# make newFiles$newGroupsFile into dataframe, otherwise it will not work
newFiles$newGroupsFile <- as.data.frame.array(newFiles$newGroupsFile)

#dev.off()
results <- Main_Doublet_Decon(rawDataFile= newFiles$newExpressionFile,
                              groupsFile = newFiles$newGroupsFile,
                              filename= SAMPLE.ID,
                              location= "1_Data/DoubletDecon",
                              fullDataFile= NULL,
                              removeCC= FALSE,
                              species= "hsa",
                              rhop= 1.1,
                              write= F,
                              PMF = TRUE,
                              useFull= FALSE,
                              heatmap= FALSE,
                              centroids= TRUE,
                              num_doubs= 100,
                              only50= FALSE,
                              min_uniq= 4,
                              nCores= 4)

decon.doublets <- rownames(results$Final_doublets_groups)
decon.doublets <- gsub("\\.","-",decon.doublets)

# DoubletFinder
# Add modfieid Parameter sweep function to regress out nCOunt RNA if needed
# Add if else statement to regress out nCount RNA if needed before running DoubletFinder
#########################################################################################
library(DoubletFinder)
# Estimate Doublets
#################################################################################################
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep(rna, PCs = 1:50, sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res.list , GT = FALSE)

bcmvn<- find.pK(sweep.stats)

pK.1 <- as.numeric(unfactor(dplyr::arrange(bcmvn,desc(BCmetric))$pK[1]))


nExp_poi <- round(doublet.rate*length(colnames(counts)))  ## Assuming 4.6% doublet formation rate - tailor for your dataset

# Run doubletFinder with pk.1
rna <- doubletFinder(rna, PCs = 1:50, pN = 0.25, pK = pK.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
meta <- rna@meta.data


doublet.column <- paste0("DF.classifications_0.25_",pK.1,"_",nExp_poi)


doublet.calls <- rna[[doublet.column]]
colnames(doublet.calls) <- "Call"

rna.dub <- dplyr::filter(doublet.calls, Call == "Doublet")
rna.singlet <- dplyr::filter(doublet.calls, Call == "Singlet")

DF.doublets <- rownames(rna.dub)

# Intersect doublet calls
###########################################################################

head(DF.doublets)
head(decon.doublets)

doublets <- intersect(decon.doublets,DF.doublets)
rna@meta.data$Doublet.Call <- ifelse(rownames(rna@meta.data) %in% doublets,"TRUE","FALSE")
FeatureScatter(rna,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "Doublet.Call")
DimPlot(rna,group.by = "Doublet.Call")
cells <- colnames(rna)
##################################################################################################################

###########################################################
# Part 3: scRNA-seq processing after doublet removal
###########################################################
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
JackStrawPlot(rna, dims = 1:50)

# If PC1 is correalted with read depth, check to see if biological variation is corralted to PC1
if (round(abs(count_cor_PC1),2) > 0.5){
  
  if( round(abs(stromal.cor),2) >= 0.5 |
      round(abs(immune.cor),2) >= 0.5 |
      round(abs(fibroblast.cor),2) >= 0.5 |
      round(abs(endothelial.cor),2) >= 0.5 |
      round(abs(epithelial.cor),2) >= 0.5 |
      round(abs(smooth.cor),2) >= 0.5 |
      round(abs(plasma.cor),2) >= 0.5){
    
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    
    
    # Verify with inferCNV: is PC1 correlated with CNV events/Malignancy?
    #########################################################################
    # inferCNV: does PC1 also correlated with CNV/malignancy status?
    library(infercnv)
    library(stringr)
    library(Seurat)
    counts_matrix = GetAssayData(rna, layer="counts")
    
    # Identify immune clusters
    #######################################################
    # Find immune cells by relative enrichment of ESTIMATE immune signature
    library(psych)
    test <- VlnPlot(rna,features = "immune.2")
    data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
    data.immune <- dplyr::filter(data,median > 0.1)
    
    test <- VlnPlot(rna,features = "plasma.7")
    data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
    data.plasma <- dplyr::filter(data,median > 0.1)
    
    immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
    plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
    
    immune.clusters <- unique(append(immune.clusters,plasma.clusters))
    
    for (i in 1:length(immune.clusters)){
      j <- which(levels(Idents(rna)) == immune.clusters[i])
      levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
    }
    rna@meta.data$postdoublet.idents <- Idents(rna)
    idents <- data.frame(rownames(rna@meta.data),rna@meta.data$postdoublet.idents)
    
    
    colnames(idents) <- c("V1","V2")
    saveRDS(rna,"1_Data/rna_postdoublet_preinferCNV.rds")
    # Make inferCNV inputs
    
    rownames(idents) <- NULL
    colnames(idents) <- NULL
    write.table(idents,"1_Data/sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
    idents <- read.delim("1_Data/sample_annotation_file_inferCNV.txt",header = F)
    
    
    
    gtf <- read.delim(GRCH38.annotations,header = F)
    library(EnsDb.Hsapiens.v86)
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
    
    write.table(gtf,"1_Data/Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
    
    
    num.immune.clusters = length(immune.clusters)
    # create the infercnv object
    if ( num.immune.clusters == 1) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=paste0("immune.",immune.clusters[1]))
      
    } else if (num.immune.clusters == 2) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2])))
      
    } else if ( num.immune.clusters == 3 ) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3])))
    } else if (num.immune.clusters == 4) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4])))
    } else if (num.immune.clusters == 5) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5])))
    } else if (num.immune.clusters == 6) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6])))
    }else if (num.immune.clusters == 7) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7])))
    }else if (num.immune.clusters == 8) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8])))
    }else if (num.immune.clusters == 9) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9])))
    }else if (num.immune.clusters == 10) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9]),
                                                            paste0("immune.",immune.clusters[10])))
    }
    # perform infercnv operations to reveal cnv signal
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir="1_Data/output_dir_CNV_postdoublet_PassedPC1Checks",  # dir is auto-created for storing outputs
                                 cluster_by_groups=T,   # cluster
                                 denoise=T,scale_data = T,
                                 HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                                 BayesMaxPNormal = 0.4, num_threads = 8
                                 
    )
    
    
    
    regions <- read.delim("1_Data/output_dir_CNV_postdoublet_PassedPC1Checks/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
    probs <- read.delim("1_Data/output_dir_CNV_postdoublet_PassedPC1Checks/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
    
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
    ggsave(filename = "postdoublet_CNV_PC1_boxplot.png", plot = last_plot(), path = "1_Data/")
    
    data <- describeBy(boxplot.cnv$data$PC1.loading, boxplot.cnv$data$postdoublet.idents, mat = TRUE)
    data$CNV <- ifelse(data$group1 %in% cnv.groups,1,0)
    
    wilcox <- wilcox.test(data = rna@meta.data,PC1.loading~CNV.Pos)
    
    if (wilcox$p.value < 0.05){
      rna <- rna
      library(stringr)
      levels(Idents(rna)) <- str_remove(levels(Idents(rna)),"immune.")
      saveRDS(rna,"1_Data/rna_postdoublet_PassedPC1Checks.rds")
    }else{
      all.genes <- rownames(rna)
      rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
      rna <- FindNeighbors(rna,dims = 1:50)
      rna <- FindClusters(rna,resolution = 0.7)
      rna <- RunUMAP(rna,dims = 1:50)
      Idents(rna) <- "RNA_snn_res.0.7"
      
      library(infercnv)
      library(stringr)
      library(Seurat)
      counts_matrix = GetAssayData(rna, layer="counts")
      
      
      # Identify immune clusters
      #######################################################
      # Find immune cells by relative enrichment of ESTIMATE immune signature
      library(psych)
      test <- VlnPlot(rna,features = "immune.2")
      data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
      data.immune <- dplyr::filter(data,median > 0.1)
      
      test <- VlnPlot(rna,features = "plasma.7")
      data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
      data.plasma <- dplyr::filter(data,median > 0.1)
      
      immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
      plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
      
      immune.clusters <- unique(append(immune.clusters,plasma.clusters))
      
      for (i in 1:length(immune.clusters)){
        j <- which(levels(Idents(rna)) == immune.clusters[i])
        levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
      }
      rna@meta.data$postdoublet.idents <- Idents(rna)
      idents <- data.frame(rownames(rna@meta.data),rna@meta.data$postdoublet.idents)
      saveRDS(rna,"1_Data/rna_postdoublet_preinferCNV.rds")
      #MAke inferCNV input
      
      colnames(idents) <- c("V1","V2")
      
      
      rownames(idents) <- NULL
      colnames(idents) <- NULL
      write.table(idents,"1_Data/sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
      idents <- read.delim("1_Data/sample_annotation_file_inferCNV.txt",header = F)
      
      
      
      gtf <- read.delim(GRCH38.annotations,header = F)
      library(EnsDb.Hsapiens.v86)
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
      
      write.table(gtf,"1_Data/Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
      
      
      num.immune.clusters = length(immune.clusters)
      # create the infercnv object
      if ( num.immune.clusters == 1) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=paste0("immune.",immune.clusters[1]))
        
      } else if (num.immune.clusters == 2) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2])))
        
      } else if ( num.immune.clusters == 3 ) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3])))
      } else if (num.immune.clusters == 4) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4])))
      } else if (num.immune.clusters == 5) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5])))
      } else if (num.immune.clusters == 6) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6])))
      }else if (num.immune.clusters == 7) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6]),
                                                              paste0("immune.",immune.clusters[7])))
      }else if (num.immune.clusters == 8) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6]),
                                                              paste0("immune.",immune.clusters[7]),
                                                              paste0("immune.",immune.clusters[8])))
      }else if (num.immune.clusters == 9) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6]),
                                                              paste0("immune.",immune.clusters[7]),
                                                              paste0("immune.",immune.clusters[8]),
                                                              paste0("immune.",immune.clusters[9])))
      }else if (num.immune.clusters == 10) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                            gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6]),
                                                              paste0("immune.",immune.clusters[7]),
                                                              paste0("immune.",immune.clusters[8]),
                                                              paste0("immune.",immune.clusters[9]),
                                                              paste0("immune.",immune.clusters[10])))
      }
      # perform infercnv operations to reveal cnv signal
      infercnv_obj = infercnv::run(infercnv_obj,
                                   cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                   out_dir="1_Data/output_dir_CNV_postdoublet_FailedCNVTest",  # dir is auto-created for storing outputs
                                   cluster_by_groups=T,   # cluster
                                   denoise=T,scale_data = T,
                                   HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                                   BayesMaxPNormal = 0.4, num_threads = 8
                                   
      )
      
      
      
      regions <- read.delim("1_Data/output_dir_CNV_postdoublet_FailedCNVTest/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
      probs <- read.delim("1_Data/output_dir_CNV_postdoublet_FailedCNVTest/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
      
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
      
      saveRDS(rna,"1_Data/rna_postdoublet_FailedCNVTest.rds")
    }
    
  }else{
    all.genes <- rownames(rna)
    rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    
    # Proceed with CNV
    library(infercnv)
    library(stringr)
    library(Seurat)
    counts_matrix = GetAssayData(rna, layer="counts")
    
    # Identify immune clusters
    #######################################################
    # Find immune cells by relative enrichment of ESTIMATE immune signature
    library(psych)
    test <- VlnPlot(rna,features = "immune.2")
    data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
    data.immune <- dplyr::filter(data,median > 0.1)
    
    test <- VlnPlot(rna,features = "plasma.7")
    data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
    data.plasma <- dplyr::filter(data,median > 0.1)
    
    immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
    plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
    
    immune.clusters <- unique(append(immune.clusters,plasma.clusters))
    
    for (i in 1:length(immune.clusters)){
      j <- which(levels(Idents(rna)) == immune.clusters[i])
      levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
    }
    rna@meta.data$postdoublet.idents <- Idents(rna)
    idents <- data.frame(rownames(rna@meta.data),rna@meta.data$postdoublet.idents)
    
    
    
    colnames(idents) <- c("V1","V2")
    
    saveRDS(rna,"1_Data/rna_postdoublet_preinferCNV.rds")
    # Make inferCNV inputs
    
    rownames(idents) <- NULL
    colnames(idents) <- NULL
    write.table(idents,"1_Data/sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
    idents <- read.delim("1_Data/sample_annotation_file_inferCNV.txt",header = F)
    
    
    
    gtf <- read.delim(GRCH38.annotations,header = F)
    library(EnsDb.Hsapiens.v86)
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
    
    write.table(gtf,"1_Data/Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
    
    
    num.immune.clusters = length(immune.clusters)
    # create the infercnv object
    if ( num.immune.clusters == 1) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=paste0("immune.",immune.clusters[1]))
      
    } else if (num.immune.clusters == 2) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2])))
      
    } else if ( num.immune.clusters == 3 ) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3])))
    } else if (num.immune.clusters == 4) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4])))
    } else if (num.immune.clusters == 5) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5])))
    } else if (num.immune.clusters == 6) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6])))
    }else if (num.immune.clusters == 7) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7])))
    }else if (num.immune.clusters == 8) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8])))
    }else if (num.immune.clusters == 9) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9])))
    }else if (num.immune.clusters == 10) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                          gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9]),
                                                            paste0("immune.",immune.clusters[10])))
    }
    # perform infercnv operations to reveal cnv signal
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir="1_Data/output_dir_CNV_postdoublet_FailedCorTest",  # dir is auto-created for storing outputs
                                 cluster_by_groups=T,   # cluster
                                 denoise=T,scale_data = T,
                                 HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                                 BayesMaxPNormal = 0.4, num_threads = 8
                                 
    )
    
    
    
    regions <- read.delim("1_Data/output_dir_CNV_postdoublet_FailedCorTest/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
    probs <- read.delim("1_Data/output_dir_CNV_postdoublet_FailedCorTest/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
    
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
    saveRDS(rna,"1_Data/rna_postdoublet_FailedCorTest.rds")
  }
}else{
  
  rna <- FindNeighbors(rna,dims = 1:50)
  rna <- FindClusters(rna,resolution = 0.7)
  rna <- RunUMAP(rna,dims = 1:50)
  Idents(rna) <- "RNA_snn_res.0.7"
  
  # Proceed with CNV
  library(infercnv)
  library(stringr)
  library(Seurat)
  counts_matrix = GetAssayData(rna, layer="counts")
  
  
  # Identify immune clusters
  #######################################################
  # Find immune cells by relative enrichment of ESTIMATE immune signature
  library(psych)
  test <- VlnPlot(rna,features = "immune.2")
  data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
  data.immune <- dplyr::filter(data,median > 0.1)
  
  test <- VlnPlot(rna,features = "plasma.7")
  data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
  data.plasma <- dplyr::filter(data,median > 0.1)
  
  immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
  plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
  
  immune.clusters <- unique(append(immune.clusters,plasma.clusters))
  
  for (i in 1:length(immune.clusters)){
    j <- which(levels(Idents(rna)) == immune.clusters[i])
    levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
  }
  rna@meta.data$postdoublet.idents <- Idents(rna)
  idents <- data.frame(rownames(rna@meta.data),rna@meta.data$postdoublet.idents)
  
  
  colnames(idents) <- c("V1","V2")
  saveRDS(rna,"1_Data/rna_postdoublet_preinferCNV.rds")
  # Make inferCNV inputs
  rownames(idents) <- NULL
  colnames(idents) <- NULL
  write.table(idents,"1_Data/sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
  idents <- read.delim("1_Data/sample_annotation_file_inferCNV.txt",header = F)
  
  
  
  gtf <- read.delim(GRCH38.annotations,header = F)
  library(EnsDb.Hsapiens.v86)
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
  
  write.table(gtf,"1_Data/Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
  
  
  num.immune.clusters = length(immune.clusters)
  # create the infercnv object
  # create the infercnv object
  if ( num.immune.clusters == 1) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=paste0("immune.",immune.clusters[1]))
    
  } else if (num.immune.clusters == 2) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2])))
    
  } else if ( num.immune.clusters == 3 ) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3])))
  } else if (num.immune.clusters == 4) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4])))
  } else if (num.immune.clusters == 5) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5])))
  } else if (num.immune.clusters == 6) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6])))
  }else if (num.immune.clusters == 7) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6]),
                                                          paste0("immune.",immune.clusters[7])))
  }else if (num.immune.clusters == 8) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6]),
                                                          paste0("immune.",immune.clusters[7]),
                                                          paste0("immune.",immune.clusters[8])))
  }else if (num.immune.clusters == 9) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6]),
                                                          paste0("immune.",immune.clusters[7]),
                                                          paste0("immune.",immune.clusters[8]),
                                                          paste0("immune.",immune.clusters[9])))
  }else if (num.immune.clusters == 10) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="1_Data/sample_annotation_file_inferCNV.txt",
                                        gene_order_file="1_Data/Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6]),
                                                          paste0("immune.",immune.clusters[7]),
                                                          paste0("immune.",immune.clusters[8]),
                                                          paste0("immune.",immune.clusters[9]),
                                                          paste0("immune.",immune.clusters[10])))
  }
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir="1_Data/output_dir_CNV_postdoublet_SkipChecks",  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               denoise=T,scale_data = T,
                               HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                               BayesMaxPNormal = 0.4, num_threads = 8
                               
  )
  
  
  
  regions <- read.delim("1_Data/output_dir_CNV_postdoublet_SkipChecks/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
  probs <- read.delim("1_Data/output_dir_CNV_postdoublet_SkipChecks/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
  
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
  ggsave(filename = "Postdoublet_CNV_PC1_boxplot.png", plot = last_plot(), path = "1_Data/")
  saveRDS(rna,"1_Data/rna_postdoublet_SkipChecks.rds")
}

Idents(rna)<- "RNA_snn_res.0.7"

# DEG analysis with Wilcox
##########################################################################

# Wilcox
Wilcox.markers <- FindAllMarkers(object =rna, min.pct = 0.25,only.pos = F,
                                 test.use = "wilcox")
saveRDS(Wilcox.markers,"1_Data/wilcox_DEGs.rds")


###########################################################
# Part 4: SingleR cell typing 
###########################################################

# SingleR labeling of celltypes
##########################################################################
rna.sce <- as.SingleCellExperiment(rna)

# Set reference paths:
# Wang et al. GSE111976
ref.data.counts <- readRDS("1_Data/GSE111976_ct_endo_10x.rds")
meta <- read.csv("1_Data/GSE111976_summary_10x_day_donor_ctype.csv")
rownames(meta) <- meta$X

length(which(rownames(meta) == colnames(ref.data.counts)))
ref.data.endo <- CreateSeuratObject(ref.data.counts,meta.data = meta)
Idents(ref.data.endo) <- "cell_type"
ref.data.endo <- NormalizeData(ref.data.endo)

# SingleR annotation
#####################################################

# Read in reference datasets for SingleR annotation 

# 1) Slyper et al. Nat. Medicine 2020 scRNA-seq ovarian tumor 
ref.data.endo <- as.SingleCellExperiment(ref.data.endo)

# 2) Human Primary Cell Atlas Data (microarray)
# HPCA_celldex <- celldex::HumanPrimaryCellAtlasData()
# save(HPCA_celldex, file = "1_Data/HPCA_celldex.rds")
load("1_Data/HPCA_celldex.rds")
ref.data.HPCA <- HPCA_celldex
 
# # 3) BluePrint Encode (bulk RNA-seq)
# BlueprintEncode_celldex <- celldex::BlueprintEncodeData()
# save(BlueprintEncode_celldex, file = "1_Data/BlueprintEncode_celldex.rds")
load("1_Data/BlueprintEncode_celldex.rds")
ref.data.BED <- BlueprintEncode_celldex

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

# Save Seurat object 
#date <- Sys.Date()
saveRDS(rna,paste0(SAMPLE.ID,"_scRNA_processed.rds"))


# ###########################################################################################################
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
#save(rna, file = "4_Intermediate/rna.RData")

###################################
# Continuation of cell annotation #
###################################

# Assess Mast cell enrichment to potentially rename clusters
mast1 <- StackedVlnPlot(rna,features = c("B.1","Plasma.2","Mast.3","Macrophage.4",
                                         "DC.5","T.6","NK.7","Endothelial.8","Fibroblast.9",
                                         "Epithelial.10","Smooth_muscle.11"))
mast1
mast2 <- StackedVlnPlot(rna,features = c("TPSB2","TPSAB1","KIT"))
mast2
vln.df <- VlnPlot(rna,features = "Mast_3_gene.12")                                 # Expression of specific genes pertaning to mast cells
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

#############################################################
# Aggregate all cells pertaining to their specific clusters #
#############################################################

# Group all cell types per cluster and number of cells for each cell type
# Group te data per cluster per cell type to each vector in clusters
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

# save(rna, file = "4_Intermediate/rna.RData")
# save(Wilcox.markers, file = "4_Intermediate/wilcox_markers.RData")
