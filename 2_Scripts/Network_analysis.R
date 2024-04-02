#--------------------Libraries--------------------
library(WGCNA)
library(Seurat)
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(DESeq2)
library(gridExtra)
library(CorLevelPlot)
library(igraph)
library(enrichR)
library(RCy3)

#--------------------Load data--------------------
rna <- readRDS("1_Data/Created_data/pre_network_rna.rds")

#--------------------Processing--------------------
# Firstly, the counts per gene have to be collapsed into cell types. The data should contain
# reads per gene by each cell type as a column. The normalized data will be the reads' source.
# Additionally, meta data will be selected as features can be selected during th WGCNA analysis

Idents(rna) <- rna$cell.type # ID for clusters is their cell type

#### Separate the pre-existing seurat object based on cell type 
for(i in levels(factor(Idents(rna)))){
  clus_num <- strsplit(i,"-")[[1]][1]               # extract cluster number
  new_name_cluster <- paste0("Cluster_",clus_num)   # create new cluster name
  assign(new_name_cluster, subset(rna, idents = i)) # assign the cluster seurat data to the new object
}

#### Collapse clusters' reads and make one single matrix 
# The following function extracts the total reads per gene per cluster
mat_count <- function(cluster, experiment = "RNA"){
  count_per_gene <- rowSums(cluster[[as.character(experiment)]]$counts)
  count_per_gene <- as.data.frame(count_per_gene)
  return(count_per_gene)
}

# The following function extracts the total reads per gene per cluster. Afterwards, it creates a dataframe
# combining the total reads of all clusters
combined_mat <- function(cluster_list){
  empty_vec <- c()
  col_names <- c()
  num <- 1
  for(i in cluster_list){
    new_col <- mat_count(i)
    empty_vec[num] <- new_col
    num <- num +1
  }
  empty_vec <- as.data.frame(empty_vec)
  for(x in 0:length(cluster_list)){
    col_names[x] <- sprintf("Cluster_%d", x-1)
  }
  colnames(empty_vec) <- col_names
  return(empty_vec)
}

cl_list <- c(Cluster_0,Cluster_1,Cluster_2,Cluster_3,Cluster_4,Cluster_5,Cluster_6,
             Cluster_7,Cluster_8,Cluster_9,Cluster_10,Cluster_11,Cluster_12,Cluster_13)

# Final dataframe with all reads per gene per cluster
clusters_total_counts <- combined_mat(cl_list)
rownames(clusters_total_counts) <- rownames(rna)

#### Create metadata dataframe for the samples
f <- 1
cell_type <- c()
for(i in cl_list){
  fast <- levels(factor(i$cell.type))
  fast <- strsplit(fast, "-")[[1]][2]
  fast <- gsub(" ", "_", fast)
  cell_type[f] <- fast
  f <- f + 1
}

f <- 1
cancer <- c()
for(i in cl_list){
  if(levels(factor(i$cancer)) == "Endometrial cancer"){
    cancer[f] <- 1
  }else{
    cancer[f] <- 0
  }
  f <- f + 1
}

data_fr <- cbind(cell_type,cancer)
colData <- data.frame(data_fr)
rownames(colData) <- colnames(clusters_total_counts)

#--------------------Quality control and outlier detection--------------------
#### Basic outlier detection
gsg <- goodSamplesGenes(t(clusters_total_counts))     # Requires the rows to be the samples, so transpose is necessary
gsg$allOK                                             # Value == T, therefore no apparent outlier sample and/or gene

#### Detect outlier samples (hierarchical clustering)
htree <- hclust(dist(t(clusters_total_counts)), method = "average") # Hierarchical clustering based on distances between samples
plot(htree, xlab = "Distance", sub ="")                             # Seems as Cluster_0 may be an outlier

#### Detect outlier samples based on PCA
pca <- prcomp(t(clusters_total_counts))                         # Performs PCA on given data matrix
pca.dat <- pca$x                                                # PCA information stored
pca.var <- pca$sdev^2                                           # Calculate variance explained by each PC
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)  # Calculate percentage of variance explained by each PC

# Plot PCA. Cluster_0 seems to differ from the rest, therefore will be considered outlier
cluster_pca_outlier <- ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point()+
  geom_text(label = rownames(pca.dat), nudge_y = 10000) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))+
  theme_classic()
cluster_pca_outlier

#### Remove outlier sample
sample_exclude <- "Cluster_0"                                                           # Excluded samples
clusters_total_counts <- clusters_total_counts[,!(colnames(clusters_total_counts) %in%  # Subset data
                                                    sample_exclude)] 
colData <- colData %>% filter(!row.names(.) %in% sample_exclude)                        # Remove outliers from meta data too

#----------Data normalization----------
#### Normalize data based on VST method
# Create dataset used for the DESeq2 normalization method
dds <- DESeqDataSetFromMatrix(countData = clusters_total_counts,  # Read matrix for each gene per sample
                              colData = colData,                  # Meta data associated with reads and samples
                              design = ~ 1)                       # not specifying model 

# Remove all genes with <10 reads in more than 90% of samples (13*0.9=12) suggested by WGCNA on RNAseq FAQ
dds75 <- dds[rowSums(counts(dds) >= 10) >= 12,]  # Subset based on parameters established above
dim(dds75)                                       # 5,535 genes remain across 13 samples

# Perform variance stabilization
# Calculates a variance stabilizing transformation 
# Yields a matrix of values which are now approximately homoskedastic (having constant variance along the range of mean values)
dds_norm <- vst(dds75)                 # vst normalization
head(assay(dds_norm))                  # Visualize normalized data
norm.counts <- assay(dds_norm) %>% t() # Store normalized counts and transpose for further manipulation


#----------Network Construction----------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,            # Expresison matrix (sample x genes)
                         powerVector = power,    # vector of soft thresholding powers for which the scale free topology fit 
                         # indices are to be calculated.
                         networkType = "signed", # network type
                         verbose = 5)            # integer level of verbosity. Zero means silent, higher values make the 
# output progressively more and more verbose.
sft.data <- sft$fitIndices                       # Store data

# Visualization to pick power
# Ideally, the best power  value has the highest R² value and the lowest mean connectivity (mean.k)

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)                   # Selected power will be 26, but could also be 28

# Convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric) # Assign numerical and normalized gene counts
soft_power <- 26                                 # Previously assigned soft power parameter
temp_cor <- cor                                  # Assign native correlation function to temporal location
cor <- WGCNA::cor                                # We need to mask the correlation function nartive to R with the WGCNA's function

# Memory estimate w.r.t blocksize
# The function performs an automatic construction of the nwtworks. It also detects modules (clusters) 
# on the expression datasets in a block-wise manner
bwnet <- blockwiseModules(norm.counts,           # Data matrix (genes x samples)
                          maxBlockSize = 14000,  # Maximum block size for module detection. It is dependant on RAM. The bigger the better
                          TOMType = "signed",    # Type of adjacency network 
                          power = soft_power,    # soft-thresholding value previously established
                          mergeCutHeight = 0.25, # Dendrogram cut height for module merging
                          numericLabels = F,     # should the returned modules be labeled by colors (FALSE), or by numbers (TRUE)
                          randomSeed = 1234,     # Seed for the random number generator before the function starts
                          verbose = 3)           # How much verbose would we like 
cor <- temp_cor                                  # Return original correlation function to avoid masking

#----------Module Eigengenes----------
module_eigengenes <- bwnet$MEs # Dataframe with modules' eigen genes (a summary of the gene expression patterns within a module)
                               # Representing genes that explain most variance in performed PCA

# Print out a preview
head(module_eigengenes)        # Each sample, the module's eigen genes have been calculated

# Get number of genes for each module (color)
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

#----------Relate modules to specific traits----------
# My traits will be cell type and if said cell types are cancerous
traits <- as.numeric(colData[,2]) # Binarize cancer state or not for each cell cluster/cell type

# Take as factor each possible cell type level
colData$cell_type <- factor(colData$cell_type, levels = c("Ciliated", "Endothelia", "Lymphocytes", "Macrophages", 
                                                          "Smooth_muscle_cells", "Stromal_fibroblasts", "Unciliated_epithelia_1"))

# Binarization for each cell type
ce_type_out <- binarizeCategoricalColumns(as.data.frame(colData$cell_type),
                                          includePairwise = F,
                                          includeLevelVsAll = T,
                                          minCount = 1)

# Apparently, the first cell type is not being detected, so it will be manually added
ce_type_out$`colData$cell_type.Ciliated.vs.all` <- c(1,0,0,0,0,0,1,0,0,0,0,0,0)
ce_type_out <- ce_type_out[,c(7,1,2,3,4,5,6)]      # Reorder columns

# Combine binarized traits
traits <- cbind(traits, ce_type_out)             # Combine by column the cancer trait and the cell type
rownames(traits) <- rownames(module_eigengenes)  # Assign row names

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)                    # 13 clusters remaining
nGenes <- ncol(norm.counts)                      # 5,535 genes reamining

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')           # Correlation between the module eigen genes and the selected traits, in this case,
                                                                         # the traits are the cell types as well as cancerous or not
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples) # Determine statistical significance of the modules' influence on determining if the
                                                                         # clusters can explain the cell type and if they are cancerous or not

# Change the column names for easier plot understanding
colnames(traits) <- c("Cancerous", "Ciliated", "Endothelia", "Lymphocytes", "Macrophages",
                      "Smooth_muscle", "Stromal_fiborblasts", "Unciliated")

# Visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

# Significant correlations will be the ones with "*"
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:25],
             y = names(heatmap.data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# Selected modules are: turquoise for cancerous cells, and lightcyan for ciliated cell type
module.gene.mapping <- as.data.frame(bwnet$colors)                                               # Extract which genes correspond to which module
turq_genes <- module.gene.mapping %>% filter(`bwnet$colors` == 'turquoise') %>% rownames()       # Genes associated to turquoise module (Cancerous)
lightcyan_genes <- module.gene.mapping %>% filter(`bwnet$colors` == 'lightcyan') %>% rownames()  # Genes associated to light cyan module (Ciliated)

#----------Gene enrichment analysis----------
# Set databases for enrichment analysis
dbs <- listEnrichrDbs()
View(dbs)
dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023")

#### Turquoise genes
enriched_t <- enrichr(turq_genes, dbs)
plotEnrich(enriched_t[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", xlab = "") 

#### Light cyan genes
enriched_c <- enrichr(lightcyan_genes, dbs)
plotEnrich(enriched_c[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", xlab = "")

#----------Edge list creation for network----------
module_edge_df <- data.frame(gene_id = names(bwnet$colors), colors = bwnet$colors)   # Dataframe with modules' genes and colors
rownames(module_edge_df) <- c()                                                      # Remove row names
modules_of_interest <- c("lightcyan", "turquoise")                                   # Set modules of interest

#### CYAN
module_genes_cyan <- module_edge_df %>% subset(colors %in% modules_of_interest[1])   # Subset only genes associated with modules of interest
gene_expr_cyan <- norm.counts[,module_genes_cyan$gene_id]                            # Normalized expression matrix only for modules' of interest genes
gene_expr_cyan <- t(gene_expr_cyan)                                                  # Matrix transpose

TOM_cyan <- TOMsimilarityFromExpr(t(gene_expr_cyan), power = soft_power)             # Calculate topological overlap matrix from matrix
row.names(TOM_cyan) <- row.names(gene_expr_cyan)                                     # Assign genes' id as row names
colnames(TOM_cyan) = row.names(gene_expr_cyan)                                       # Assign genes' id as column names
diag(TOM_cyan) <- 0                                                                  # Make all diagonal values turn to 0

#### TURQUOISE
module_genes_turq <- module_edge_df %>% subset(colors %in% modules_of_interest[2])   # Subset only genes associated with modules of interest
gene_expr_turq <- norm.counts[,module_genes_turq$gene_id]                            # Normalized expression matrix only for modules' of interest genes
gene_expr_turq <- t(gene_expr_turq)                                                  # Matrix transpose

TOM_turq <- TOMsimilarityFromExpr(t(gene_expr_turq), power = soft_power)             # Calculate topological overlap matrix from matrix
row.names(TOM_turq) <- row.names(gene_expr_turq)                                     # Assign genes' id as row names
colnames(TOM_turq) = row.names(gene_expr_turq)                                       # Assign genes' id as column names
diag(TOM_turq) <- 0                                                                  # Make all diagonal values turn to 0

#---------Reduce size for bigger module----------
top_30 <- 30                                                     # Top 30 
norm_counts_turq <- norm.counts[,turq_genes]                     # Get normalized counts only for the genes associated with turquoise module
IMConn <- softConnectivity(norm_counts_turq, power = soft_power) # Construct adjacency matrix for each node (gene) and calculates connectivity
                                                                 # which is the sum of the adjacency to the other nodes
top_genes <- (rank(-IMConn) <= top_30)                           # Select which genes are the top 30
TOM_turq_30 <- TOM_turq[top_genes, top_genes]                    # Subset data 

#---------Network Visualization and analysis----------
#### Create igraph network
cyan_net <- graph_from_adjacency_matrix(TOM_cyan, mode = "undirected", weighted = T, diag = F)
turq_net <- graph_from_adjacency_matrix(TOM_turq_30, mode = "undirected", weighted = T, diag = F)

#### Thresholding based on adjacency
summary(c(TOM_turq_30)) # Find where only 25% of connections remain: 0.1556
summary(c(TOM_cyan))    # Find where only 25% of connections remain: 0.031478

turq_edges_to_trim <- E(turq_net)[E(turq_net)$weight < 0.1575]   # Select edges to delete based on threshold
turq_net_trimmed <- delete_edges(turq_net, turq_edges_to_trim)   # Remove selected edges

cyan_edges_to_trim <- E(cyan_net)[E(cyan_net)$weight < 0.031478] # Select edges to delete based on threshold
cyan_net_trimmed <- delete_edges(cyan_net, cyan_edges_to_trim)   # Remove selected edges

# Set edges' weight into variables
turq_weight <- E(turq_net_trimmed)$weight
cyan_weight <- E(cyan_net_trimmed)$weight

# Plot
plot.igraph(cyan_net_trimmed, vertex.size = 12, edge.width=E(cyan_net_trimmed)$weight*40) + title("Light cyan module network")
plot.igraph(turq_net_trimmed, vertex.size = 12, edge.width=E(turq_net_trimmed)$weight*20) + title("Turquoise module network")

#### Clusterizatrion and centralities
# Cluster optimal
cl_opt_turq <- cluster_optimal(turq_net_trimmed)
cl_opt_cyan <- cluster_optimal(cyan_net_trimmed)

# Edge betweennes
bet_turq <- betweenness(turq_net_trimmed, directed = F, weights = turq_weight)
edge_b_turq <- cluster_edge_betweenness(turq_net_trimmed, weights = turq_weight, directed = F)
bet_cyan <- betweenness(cyan_net_trimmed, directed = F, weights = cyan_weight)
edge_b_cyan <- cluster_edge_betweenness(cyan_net_trimmed, weights = cyan_weight, directed = F)

# Eigen centrality
ei_turq <- eigen_centrality(turq_net_trimmed, weights = turq_weight)$vector
eigen_turq <- cluster_leading_eigen(turq_net_trimmed, weights = turq_weight)
ei_cyan <- eigen_centrality(cyan_net_trimmed, weights = cyan_weight)$vector
eigen_cyan <- cluster_leading_eigen(cyan_net_trimmed, weights = cyan_weight)

# Fast and Greedy
fast_turq <- cluster_fast_greedy(turq_net_trimmed)
fast_cyan <- cluster_fast_greedy(cyan_net_trimmed)

# Optimal
opti_turq <- cluster_optimal(turq_net_trimmed, weights = turq_weight)
opti_cyan <- cluster_optimal(cyan_net_trimmed, weights = cyan_weight)

# Walk
walk_turq <- cluster_walktrap(turq_net_trimmed, weights = turq_weight)
walk_cyan <- cluster_walktrap(cyan_net_trimmed, weights = cyan_weight)

# Closeness
close_turq <- as.data.frame(closeness(turq_net_trimmed, weights = turq_weight))[,1]
close_cyan <- as.data.frame(closeness(cyan_net_trimmed, weights = cyan_weight))[,1]

#### Turq clusterization plots
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
plot(cl_opt_turq, turq_net_trimmed)
plot(edge_b_turq, turq_net_trimmed)
plot(eigen_turq, turq_net_trimmed)
plot(fast_turq, turq_net_trimmed)
plot(walk_turq, turq_net_trimmed)
plot(opti_turq, turq_net_trimmed)

#### Cyan clusterization plots
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
plot(cl_opt_cyan, cyan_net_trimmed)
plot(edge_b_cyan, cyan_net_trimmed)
plot(eigen_cyan, cyan_net_trimmed)
plot(fast_cyan, cyan_net_trimmed)
plot(walk_cyan, cyan_net_trimmed)
plot(opti_cyan, cyan_net_trimmed)

#### Comparison between communities turquoise
Cluster_optimal_turq <- c(compare(cl_opt_turq,cl_opt_turq, method = "rand.adjusted"),compare(cl_opt_turq,edge_b_turq, method = "rand.adjusted"),
                          compare(cl_opt_turq,eigen_turq, method = "rand.adjusted"),compare(cl_opt_turq,fast_turq, method = "rand.adjusted"),
                          compare(cl_opt_turq,walk_turq, method = "rand.adjusted"),compare(cl_opt_turq,opti_turq, method = "rand.adjusted"))

Edge_bet_turq <- c(compare(edge_b_turq,cl_opt_turq, method = "rand.adjusted"),compare(edge_b_turq,edge_b_turq, method = "rand.adjusted"),
                   compare(edge_b_turq,eigen_turq, method = "rand.adjusted"),compare(edge_b_turq,fast_turq, method = "rand.adjusted"),
                   compare(edge_b_turq,walk_turq, method = "rand.adjusted"),compare(edge_b_turq,opti_turq, method = "rand.adjusted"))

Eigen_cent_turq <- c(compare(eigen_turq,cl_opt_turq, method = "rand.adjusted"),compare(eigen_turq,edge_b_turq, method = "rand.adjusted"),
                     compare(eigen_turq,eigen_turq, method = "rand.adjusted"),compare(eigen_turq,fast_turq, method = "rand.adjusted"),
                     compare(eigen_turq,walk_turq, method = "rand.adjusted"),compare(eigen_turq,opti_turq, method = "rand.adjusted"))

Fast_cent_turq <- c(compare(fast_turq,cl_opt_turq, method = "rand.adjusted"),compare(fast_turq,edge_b_turq, method = "rand.adjusted"),
                    compare(fast_turq,eigen_turq, method = "rand.adjusted"),compare(fast_turq,fast_turq, method = "rand.adjusted"),
                    compare(fast_turq,walk_turq, method = "rand.adjusted"),compare(fast_turq,opti_turq, method = "rand.adjusted"))

Walk_cent_turq <- c(compare(walk_turq,cl_opt_turq, method = "rand.adjusted"),compare(walk_turq,edge_b_turq, method = "rand.adjusted"),
                    compare(walk_turq,eigen_turq, method = "rand.adjusted"),compare(walk_turq,fast_turq, method = "rand.adjusted"),
                    compare(walk_turq,walk_turq, method = "rand.adjusted"),compare(walk_turq,opti_turq, method = "rand.adjusted"))

Opti_cent_turq <- c(compare(opti_turq,cl_opt_turq, method = "rand.adjusted"),compare(opti_turq,edge_b_turq, method = "rand.adjusted"),
                    compare(opti_turq,eigen_turq, method = "rand.adjusted"),compare(opti_turq,fast_turq, method = "rand.adjusted"),
                    compare(opti_turq,walk_turq, method = "rand.adjusted"),compare(opti_turq,opti_turq, method = "rand.adjusted"))

compare_clust_turq <- matrix(c(Cluster_optimal_turq,Edge_bet_turq,Eigen_cent_turq,Fast_cent_turq,Walk_cent_turq,Opti_cent_turq), ncol = 6)
compare_clust_turq <- as.data.frame(compare_clust_turq)
names_mat <- c("Cl_optimal", "Edge_between", "Eigen_cluster", "Fast_cluster", "Walk_cluster","Opti_cluster")
rownames(compare_clust_turq) <- names_mat
colnames(compare_clust_turq) <- names_mat
compare_clust_turq

#### Comparison between communities cyan
Cluster_optimal_cyan <- c(compare(cl_opt_cyan,cl_opt_cyan, method = "rand.adjusted"),compare(cl_opt_cyan,edge_b_cyan, method = "rand.adjusted"),
                          compare(cl_opt_cyan,eigen_cyan, method = "rand.adjusted"),compare(cl_opt_cyan,fast_cyan, method = "rand.adjusted"),
                          compare(cl_opt_cyan,walk_cyan, method = "rand.adjusted"),compare(cl_opt_cyan,opti_cyan, method = "rand.adjusted"))

Edge_bet_cyan <- c(compare(edge_b_cyan,cl_opt_cyan, method = "rand.adjusted"),compare(edge_b_cyan,edge_b_cyan, method = "rand.adjusted"),
                   compare(edge_b_cyan,eigen_cyan, method = "rand.adjusted"),compare(edge_b_cyan,fast_cyan, method = "rand.adjusted"),
                   compare(edge_b_cyan,walk_cyan, method = "rand.adjusted"),compare(edge_b_cyan,opti_cyan, method = "rand.adjusted"))

Eigen_cent_cyan <- c(compare(eigen_cyan,cl_opt_cyan, method = "rand.adjusted"),compare(eigen_cyan,edge_b_cyan, method = "rand.adjusted"),
                     compare(eigen_cyan,eigen_cyan, method = "rand.adjusted"),compare(eigen_cyan,fast_cyan, method = "rand.adjusted"),
                     compare(eigen_cyan,walk_cyan, method = "rand.adjusted"),compare(eigen_cyan,opti_cyan, method = "rand.adjusted"))

Fast_cent_cyan <- c(compare(fast_cyan,cl_opt_cyan, method = "rand.adjusted"),compare(fast_cyan,edge_b_cyan, method = "rand.adjusted"),
                    compare(fast_cyan,eigen_cyan, method = "rand.adjusted"),compare(fast_cyan,fast_cyan, method = "rand.adjusted"),
                    compare(fast_cyan,walk_cyan, method = "rand.adjusted"),compare(fast_cyan,opti_cyan, method = "rand.adjusted"))

Walk_cent_cyan <- c(compare(walk_cyan,cl_opt_cyan, method = "rand.adjusted"),compare(walk_cyan,edge_b_cyan, method = "rand.adjusted"),
                    compare(walk_cyan,eigen_cyan, method = "rand.adjusted"),compare(walk_cyan,fast_cyan, method = "rand.adjusted"),
                    compare(walk_cyan,walk_cyan, method = "rand.adjusted"),compare(walk_cyan,opti_cyan, method = "rand.adjusted"))

Opti_cent_cyan <- c(compare(opti_cyan,cl_opt_cyan, method = "rand.adjusted"),compare(opti_cyan,edge_b_cyan, method = "rand.adjusted"),
                    compare(opti_cyan,eigen_cyan, method = "rand.adjusted"),compare(opti_cyan,fast_cyan, method = "rand.adjusted"),
                    compare(opti_cyan,walk_cyan, method = "rand.adjusted"),compare(opti_cyan,opti_cyan, method = "rand.adjusted"))

compare_clust_cyan <- matrix(c(Cluster_optimal_cyan,Edge_bet_cyan,Eigen_cent_cyan,Fast_cent_cyan,Walk_cent_cyan,Opti_cent_cyan), ncol = 6)
compare_clust_cyan <- as.data.frame(compare_clust_cyan)
names_mat <- c("Cl_optimal", "Edge_between", "Eigen_cluster", "Fast_cluster", "Walk_cluster","Opti_cluster")
rownames(compare_clust_cyan) <- names_mat
colnames(compare_clust_cyan) <- names_mat
compare_clust_cyan

#-----Cytoscape
createNetworkFromIgraph(turq_net_trimmed, title = "Trimmed turquoise module network")
createNetworkFromIgraph(cyan_net_trimmed, title = "Trimmed cyan module network")

ff <- cbind(bet_turq, ei_turq, close_turq)
cor(ff, use = "complete.obs")

