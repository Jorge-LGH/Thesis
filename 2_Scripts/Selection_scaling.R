#-----------Load libraries---------
library(Seurat)
library(tidyverse)

#-----------Load data--------------
load("1_Data/seu_ob.RData")


#-----------Preview----------------
# Will check if the normalization method was correctly executed 
# Average of most expressed genes across all cells
exp_genes<- apply(seu_ob@assays$RNA@data, 1, mean)
exp_genes<- sort(exp_genes, decreasing = T)
head(exp_genes, 50) # Most are ribosomal proteins, check it

express_distri <- as.data.frame(seu_ob@assays$RNA@data[,1:100]) %>% # Convert into dataframe and keep reads per cell
  pivot_longer(cols = everything(),                                 # increase number of rows and diminish columns 
                names_to = "cell",                                  # makes a column of cells, they repeat 27259 times, one per gene
                values_to = "expression") %>%                       # expression per cell per gene
  ggplot(aes(x=expression, group=cell)) +                           # plot
  geom_density() +
  coord_cartesian(ylim=c(0,1), xlim=c(0,3)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                  colour = "black"))

express_distri
ggsave("3_Figures/express_distri.png")

#---------Feature selection-------------
# Now we'll be selecting the features that actually have some differences and can give more information
# The features are the genes, obviously, as they are the ones changing expression between each cell type and will be of use

# Select variable features
seu_ob <- FindVariableFeatures(seu_ob, selection.method = "vst", nfeatures = 2000) # features that are outliers on a 'mean variability plot'.
# I must check the other feature selection methods so I can actually give a reason of choosing one over the other

variance_data <- as_tibble(HVFInfo(seu_ob), rownames = "Gene") %>% # make a dataframe of highly variable features
  mutate(hypervariable = Gene %in% VariableFeatures(seu_ob))       # select features that are considered variable features
head(variance_data, 20)

var_genes_plot <- variance_data %>% 
  ggplot(aes(log(mean),log(variance),color=hypervariable)) + 
  geom_point() + 
  scale_color_manual(values=c("black","red")) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))

var_genes_plot
ggsave("3_Figures/var_genes_plot.png")

#-----------Scaling data----------------
all_genes <- rownames(seu_ob) # all gene names
seu_ob <- ScaleData(seu_ob, features = all_genes) # check for possible models like linear and whatnot

# Running dimensional reductions, they might be PCA or TSNE
# Must check the pros and cons for each method, I think the standard is PCA but it might be a little bit skewed
# Perhaps even UMAP might be the answer
seu_ob <- RunPCA(seu_ob)
seu_ob <- RunTSNE(seu_ob)

pca_plot <- DimPlot(seu_ob, reduction = "pca")
tsne_plot <- DimPlot(seu_ob, reduction = "tsne")
pca_plot
tsne_plot
save(pca_plot, file = "3_Figures/pca_plot.RData")
# save(tsne_plot, file = "3_Figures/tsne_plot.RData")

