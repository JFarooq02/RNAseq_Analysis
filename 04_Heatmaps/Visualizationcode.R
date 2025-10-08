# =========================================
# RNA-seq Heatmap Visualization with DESeq2
# =========================================

# Set working directory
setwd("C:/Users/Javeria/Desktop/RNA_SEQ")

# ----------------------------
# Install / load required libraries
# ----------------------------


library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(matrixStats)

# ----------------------------
# Load matrices
# ----------------------------
gene_matrix <- read.csv("nuclear_gene_count_matrix.csv", row.names = 1)
transcript_matrix <- read.csv("transcript_count_matrix.csv", row.names = 1)

organelle_matrix <- read_excel("Organelle_gene_count_matrix.xlsx")
organelle_matrix <- as.data.frame(organelle_matrix)
rownames(organelle_matrix) <- organelle_matrix[,1]
organelle_matrix <- organelle_matrix[,-1]

# ----------------------------
# DESeq2 normalization
# ----------------------------
normalize_counts <- function(count_matrix){
  col_data <- data.frame(row.names = colnames(count_matrix))
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = col_data,
                                design = ~ 1)
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized=TRUE)
  return(log2(norm_counts + 1))
}

gene_norm <- normalize_counts(gene_matrix)
transcript_norm <- normalize_counts(transcript_matrix)
organelle_norm <- normalize_counts(organelle_matrix)

# ----------------------------
# Filter zero/constant rows
# ----------------------------
filter_matrix <- function(mat, min_var = 0.01){
  mat[rowVars(as.matrix(mat)) > min_var, ]
}

gene_norm <- filter_matrix(gene_norm)
transcript_norm <- filter_matrix(transcript_norm)
organelle_norm <- filter_matrix(organelle_norm)

# ----------------------------
# Function: select top N most variable genes
# ----------------------------
select_top_variable <- function(mat, top_n = 200){
  vars <- rowVars(as.matrix(mat))
  top_genes <- order(vars, decreasing = TRUE)[1:min(top_n, nrow(mat))]
  mat[top_genes, ]
}

# ----------------------------
# Subset top variable genes (adjusted per matrix)
# ----------------------------
gene_norm <- select_top_variable(gene_norm, top_n = min(200, nrow(gene_norm)))
transcript_norm <- select_top_variable(transcript_norm, top_n = min(200, nrow(transcript_norm)))
organelle_norm <- select_top_variable(organelle_norm, top_n = min(50, nrow(organelle_norm)))

# ----------------------------
# Function: Plot heatmap
# ----------------------------
plot_heatmap <- function(mat, title, filename=NULL){
  if(!is.null(filename)) png(filename, width=1200, height=800)
  pheatmap(mat,
           scale="row",
           clustering_distance_rows="euclidean",
           clustering_distance_cols="euclidean",
           clustering_method="complete",
           color=colorRampPalette(rev(brewer.pal(7,"RdYlBu")))(100),
           main=title,
           show_rownames=TRUE,
           show_colnames=TRUE,
           cluster_cols = ifelse(ncol(mat) > 2, TRUE, FALSE))  # avoid useless clustering for 2 samples
  if(!is.null(filename)) dev.off()
}

# ----------------------------
# Generate heatmaps
# ----------------------------
plot_heatmap(gene_norm, "Nuclear Gene Count Matrix", "Nuclear_Gene_Heatmap.png")
plot_heatmap(transcript_norm, "Transcript Count Matrix", "Transcript_Heatmap.png")
plot_heatmap(organelle_norm, "Organelle Gene Matrix", "Organelle_Heatmap.png")
