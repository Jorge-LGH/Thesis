#--------------------Libraries---------------------
# Load necessary libraries for analysis
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
library(msigdbr)
library(fgsea)
library(tibble)
library(DoubletDecon)
library(DoubletFinder)
library(Signac)

#--------------------Load previously created seurat object---------------------
rna <- readRDS("1_Data/Created_data/rna_predoublet_SkipChecks.rds")

#--------------------Doublet removal---------------------
# Doublet removal consists of data deconvolution, doublet finder use and intersection
# Data deconvolution is performed by DoubletDecon package
# DoubletFinder performs identification of putative existing doublets
# Intersection of deconvoluted data and doublet finder signals "real" doublets

#### DoubletDecon
# Add if else statement to regress out nCount RNA if needed before running DoubletDecon
Idents(rna) <- rna$seurat_clusters
idents.length <- length(levels(Idents(rna)))
for (i in levels(Idents(rna))){
  i <- as.numeric(i)
  levels(Idents(rna))[i] <- i -1
}

# Assign object to new variable
seuratObject <- rna

# Takes normalized expression matrix, discriminating genes
newFiles <- Improved_Seurat_Pre_Process(seuratObject, num_genes=50, write_files=F)

# make newFiles$newGroupsFile into dataframe, otherwise it will not work
newFiles$newGroupsFile <- as.data.frame.array(newFiles$newGroupsFile)

# Function identifies clusters of doublets via deconvolution analysis and unique
# gene expression and individual doublet cells with deconvolution analysis
results <- Main_Doublet_Decon(rawDataFile= newFiles$newExpressionFile,
                              groupsFile = newFiles$newGroupsFile,
                              filename= SAMPLE.ID,
                              location= "1_Data/Created_data/DoubletDecon",
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

# Doublets according to deconvolution analysis
decon.doublets <- rownames(results$Final_doublets_groups)
decon.doublets <- gsub("\\.","-",decon.doublets)

#### Doublet Finder
## pK Identification (no ground-truth)
sweep.res.list <- paramSweep(rna, PCs = 1:50, sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res.list , GT = FALSE)

bcmvn<- find.pK(sweep.stats) # Makes plot

pK.1 <- as.numeric(unfactor(dplyr::arrange(bcmvn,desc(BCmetric))$pK[1]))


nExp_poi <- round(doublet.rate*length(colnames(counts)))  ## Assuming 4.6% doublet formation rate-tailor for your dataset

# Generates artificial doublets from pre-processed Seurat object
rna <- doubletFinder(rna, PCs = 1:50, pN = 0.25, pK = pK.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
meta <- rna@meta.data
doublet.column <- paste0("DF.classifications_0.25_",pK.1,"_",nExp_poi)

# Doublet or singleton identification
doublet.calls <- rna[[doublet.column]]
colnames(doublet.calls) <- "Call"

# Split data based on singletons or doublts
rna.dub <- dplyr::filter(doublet.calls, Call == "Doublet")
rna.singlet <- dplyr::filter(doublet.calls, Call == "Singlet")

DF.doublets <- rownames(rna.dub) # Only own doublet detected 

#### Intersect doublet calls
head(DF.doublets)
head(decon.doublets)

# Intersect doublets
doublets <- intersect(decon.doublets,DF.doublets)

# Metadata with doublet or singleton calls
rna@meta.data$Doublet.Call <- ifelse(rownames(rna@meta.data) %in% doublets,"TRUE","FALSE")

# Doublets in plots
FeatureScatter(rna,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "Doublet.Call")
DimPlot(rna,group.by = "Doublet.Call")
cells <- colnames(rna)

saveRDS(rna,"1_Data/Created_data/rna_postdoublet_SkipChecks.rds")
