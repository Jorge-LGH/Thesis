# Load all necessary libraries
library(dplyr)


# Loading necessary data and explain a little 

#---------panglao data base----------------
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

# 