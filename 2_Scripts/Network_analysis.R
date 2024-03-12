#--------------------Libraries--------------------
library(WGCNA)
library(Seurat)

#--------------------Load data--------------------
rna <- readRDS("1_Data/Created_data/pre_network_rna.rds")

#--------------------Processing--------------------
# Firstly, the counts per gene have to be collapsed into cell types. The data should contain
# reads per gene by each cell type as a column. The normalized data will be the reads' source.
# Additionally, meta data will be selected as features can be selected during th WGCNA analysis

Idents(rna) <- rna$cell.type # ID for clusters is their cell type

# Separate the pre-existing seurat object based on cell type 
for(i in levels(factor(Idents(rna)))){
  clus_num <- strsplit(i,"-")[[1]][1]               # extract cluster number
  new_name_cluster <- paste0("Cluster_",clus_num)   # create new cluster name
  assign(new_name_cluster, subset(rna, idents = i)) # assign the cluster seurat data to the new object
}

# Collapse clusters' reads and make one single matrix 
mat_count <- function(cluster, experiment = "RNA"){
  count_per_gene <- rowSums(cluster[[as.character(experiment)]]$data)
  count_per_gene <- as.data.frame(count_per_gene)
  return(count_per_gene)
}

exp_1 <- mat_count(Cluster_1)
exp_2 <- mat_count(Cluster_2)
exp_3 <- mat_count(Cluster_3)
exp_4 <- mat_count(Cluster_4)
exp_5 <- mat_count(Cluster_5)
exp_6 <- mat_count(Cluster_6)
exp_7 <- mat_count(Cluster_7)
exp_8 <- mat_count(Cluster_8)
exp_9 <- mat_count(Cluster_9)
exp_10 <- mat_count(Cluster_10)
exp_11 <- mat_count(Cluster_11)
exp_12 <- mat_count(Cluster_12)
exp_13 <- mat_count(Cluster_13)

df_list <- c(exp_1, exp_2, exp_3, exp_4, exp_5, exp_6, exp_7,
             exp_8, exp_9, exp_10, exp_11, exp_12, exp_13)

total_df_norm

