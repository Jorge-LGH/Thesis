# Load necessary libraries
library(Seurat)
library(ggplot2)

# Load the previous objects
load("1_Data/matrix_tsv.RData")

#----------Basic structure----------

# Create the Seurat object
# The idea for making this object is it's easier and more straightforward manipulation. It will contain the reads,
# the cell barcode identities and the genes all the samples have. It is created based on the integration of the previous
# data that was downloaded and set into the environment. However, it doesn't contemplate the gene identifiers for the cell
# types that may exist in the samples
seu_ob <- CreateSeuratObject(matrix_tsv,
                             project = "Endo_3", # Name of the project, could have more data
                             assay = "RNA",      # Name of the initial assay
                             min.cells = 1,      # min of cells that a feature must be present in
                             min.features = 1)   # min of features a cell must have to be kept

dim(seu_ob)  # 27259 genes * 6054 cells. The original matrix had 33538 genes, but a warning did say there were non-unique
             # features (rownames) and had to make them unique. Checked and only 24 gene id's were not unique. 
head(seu_ob) # barcode id, orig.ident, nCount_RNA, nFeature_RNA
             # nCount_RNA = reads per cell, nFeature_RNA = present genes per cell
str(seu_ob)  # complete seurat object's structure and stored data

#----------Expression observation----------

# Reads per feature. This section shows the total read count per feature. It also shows the average read 
# count of each gene (gene_reads/number_cells)
mst_exp <- sort(Matrix::rowSums(seu_ob@assays$RNA@counts), decreasing = T)
ave_exp <- sort(Matrix::rowMeans(seu_ob@assays$RNA@counts), decreasing = T)
head(mst_exp, 30) # Total expression/reads for the highest expressed genes
head(ave_exp, 30) # Average reads for most expressed genes

# The next plot will show the distribution of the average reads per gene
average_feats_plot <- ggplot() + aes(as.vector(ave_exp)) + 
  geom_histogram(binwidth = 10, color="darkblue", fill="lightblue") + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black")) +
  xlab("Average reads per feature") + ylab("Features")
average_feats_plot
ggsave("3_Figures/average_feats_plot.png")

# Reads per cell. This section shows the total read count per cell. It also shows the average read
# count of each feature (cell_reads/number_features). This means how many reads EACH feature has.
# Example: average of 5 means each gene of the cell has 5 reads in average
mst_reads <- sort(colSums(seu_ob@assays$RNA@counts), decreasing = T)
ave_reads <- sort(colMeans(seu_ob@assays$RNA@counts), decreasing = T)
head(mst_reads, 30) # Total gene reads in the top 30 cells with most reads
head(ave_reads, 30) # Average reads in the top 30 cells

# This plot shows the average of feature read per cell distribution
average_reads_plot <- ggplot() + aes(as.vector(ave_reads)) + 
  geom_histogram(bins = 50, color = "darkblue", fill = "lightblue") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black")) + 
  xlab("Average feature reads per cell") + ylab("Number of cells")
average_reads_plot
ggsave("3_Figures/average_reads_plot.png")

#----------Save objects for later use----------
save(seu_ob, file = "1_Data/seu_ob.RData")
# ggsave("3_Figures/average_feats_plot.png") below the ggplot
# ggsave("3_Figures/average_reads_plot.png") below the ggplot

