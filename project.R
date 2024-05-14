# annotation
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE, force = TRUE)
library(organism, character.only = TRUE)

install.packages("devtools", force = TRUE)
devtools::install_github("stephenturner/annotables")
library(annotables)

# Load Packages
library(rstudioapi)
library(matrixStats)
library(readxl)
library(dplyr)
library(ggplot2)
library(factoextra)
library(ashr)
library(plotrix)
library(reshape2)
library(magick)
library(writexl)
library(ggrepel)
library(lubridate)
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(topGO)
library(limma)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(GSVA)
library(GSEABase)
library(dendsort)
library(EnhancedVolcano)
library(RColorBrewer)
library(tibble)
library(export)
library(devEMF)
library(svglite)
library(scales)
library(enrichR)
library(readr)
library(fgsea)
library(rentrez)
library(GenomicFeatures)
library(biomaRt)
library(esquisse)
library(officer)
library(rvg)
library(dplyr)
library(DOSE)
library(stats)

# Set working directory
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)

#################################
########### MSF-6 ###############
#################################
data <- as.data.frame(read_excel(paste(dir, "/data/GSE36854/MSF-6.xlsx", sep = "")), stringsAsFactors = FALSE)
colnames(data)
colnames(data) <- c("ID", "Mock1", "Mock2", "MPXV1", "MPXV2")

# annotation
annot_data <- as.data.frame(read_excel(paste(dir, "/data/GSE36854/GPL4133_old_annotations.xlsx", sep = "")), stringsAsFactors = FALSE)

# join annotation
data <- inner_join(data, annot_data, by = 'ID')
data <- data[!is.na(data$"GENE_SYMBOL"),]
data <- data[!duplicated(data$"GENE_SYMBOL"),]
rownames(data) <- data$"GENE_SYMBOL"

data_1 <- data[, c("Mock1", "Mock2", "MPXV1", "MPXV2")]

for (i in colnames(data_1)){
  data_1[, i] <- as.numeric(data_1[, i])
  data_1[, i] <- round(1000 * data_1[, i])
}

colnames(data_1) <- c(rep("Mock", 2), rep("MSF6", 2))
condition_1 <- c(rep("Mock", 2), rep("MSF6", 2))
df_meta_1 <- data.frame(condition_1)


# Create DESeqDataSet object
data_mat_1 <- DESeqDataSetFromMatrix(
  countData = data_1,
  colData = df_meta_1,
  design = ~ condition_1
)

# counts_ij: count for gene j in sample i
data_mat_1 <- estimateSizeFactors(data_mat_1)

# Normalize using size factor
data_normalized_1 <- counts(data_mat_1, normalized = TRUE)

# Perform differential expression analysis
# DESeq: provides statistical measures such as log2 fold changes, p-values
data_DE_1 <- DESeq(data_mat_1)
contrast_1 <- c("condition_1", "MSF6", "Mock")

# alpha: the significance cutoff used for optimizing the independent filtering
#        adjusted p-value cutoff (FDR, False Discovery Rate)
data_result_1 <- results(data_DE_1, contrast = contrast_1, alpha = 0.05)

# Adds shrunk log2 fold changes
# Shrunk LFC: adjusted or "shrunken" towards a common value by borrowing
#             information from the overall distribution of FC across genes.
#             Helpful when dealing small sample sizes or genes with low counts.
data_result_1 <- lfcShrink(data_DE_1, contrast = contrast_1, res = data_result_1, type = 'ashr')

df_normalized_1 <- as.data.frame(data_normalized_1)

result_list <- list(data_result_1)
result_final <- list()
go_analysis <- list()


# %>%: pipe operator, pass the result of LHS function
#      directly as the first argument to the RHS function
j <- result_list[[1]] %>% data.frame()
j['symbol'] <- rownames(j)

# sort j in ascending order based on padj
j <- arrange(j, padj)
j <- left_join(j, annotables::grch38[,c("symbol", "entrez", "description")], by="symbol")
result_final[[1]] <- j

###################################################
###################### GSEA #######################
###################################################

result_df_1 <- data.frame(data_result_1)
gsea_1 <- result_df_1$log2FoldChange
names(gsea_1) <- rownames(data_result_1)
gsea_1 <- na.omit(gsea_1)
gsea_1 <- sort(gsea_1, decreasing = TRUE)

gse <- gseGO(geneList=gsea_1,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 5,
             maxGSSize = 800,
             pvalueCutoff = 0.99,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none")
require(DOSE)
dotplot(gse, showCategory=10, split=".sign", font.size = 9, title = "MSF-6") + facet_grid(.~.sign)


#################################
############ ZAI ################
#################################

data2 <- as.data.frame(read_excel(paste(dir, "/data/GSE24125/ZAI.xlsx", sep = "")), stringsAsFactors = FALSE)
colnames(data2)
colnames(data2) <- c("ID", "Mock1", "Mock2", "ZAI1", "ZAI2")

# annotation
annot_data2 <- as.data.frame(read_excel(paste(dir, "/data/GSE24125/GPL10913.xlsx", sep = "")), stringsAsFactors = FALSE)

# join annotation
data2 <- inner_join(data2, annot_data2, by = 'ID')
data2 <- data2[!is.na(data2$"GENE_SYMBOL"),]
data2 <- data2[!duplicated(data2$"GENE_SYMBOL"),]
rownames(data2) <- data2$"GENE_SYMBOL"

data_2 <- data2[, c("Mock1", "Mock2", "ZAI1", "ZAI2")]

for (i in colnames(data_2)){
  data_2[, i] <- as.numeric(data_2[, i])
  data_2[, i] <- 2^data_2[, i]
  data_2[, i] <- round(1000 * data_2[, i])
}

data_2 <- na.omit(data_2)

colnames(data_2) <- c(rep("Mock", 2), rep("ZAI", 2))
condition_2 <- c(rep("ZAI", 2), rep("Mock", 2))
df_meta_2 <- data.frame(condition_2)


# Create DESeqDataSet object
data_mat_2 <- DESeqDataSetFromMatrix(
  countData = data_2,
  colData = df_meta_2,
  design = ~ condition_2
)

# counts_ij: count for gene j in sample i
data_mat_2 <- estimateSizeFactors(data_mat_2)

# Normalize using size factor
data_normalized_2 <- counts(data_mat_2, normalized = TRUE)

# Perform differential expression analysis
# DESeq: provides statistical measures such as log2 fold changes, p-values
data_DE_2 <- DESeq(data_mat_2)
contrast_2 <- c("condition_2","Mock", "ZAI")

# alpha: the significance cutoff used for optimizing the independent filtering
#        adjusted p-value cutoff (FDR, False Discovery Rate)
data_result_2 <- results(data_DE_2, contrast = contrast_2, alpha = 0.05)

# Adds shrunk log2 fold changes
# Shrunk LFC: adjusted or "shrunken" towards a common value by borrowing
#             information from the overall distribution of FC across genes.
#             Helpful when dealing small sample sizes or genes with low counts.
data_result_2 <- lfcShrink(data_DE_2, contrast = contrast_2, res = data_result_2, type = 'ashr')

df_normalized_2 <- as.data.frame(data_normalized_2)

result_list_2 <- list(data_result_2)
result_final_2 <- list()
go_analysis_2 <- list()


# %>%: pipe operator, pass the result of LHS function
#      directly as the first argument to the RHS function
j <- result_list_2[[1]] %>% data.frame()
j['symbol'] <- rownames(j)

# sort j in ascending order based on padj
j <- arrange(j, padj)
j <- left_join(j, annotables::grch38[,c("symbol", "entrez", "description")], by="symbol")
result_final_2[[1]] <- j

###################################################
###################### GSEA #######################
###################################################

result_df_2 <- data.frame(data_result_2)
gsea_2 <- result_df_2$log2FoldChange
names(gsea_2) <- rownames(data_result_2)
gsea_2 <- na.omit(gsea_2)
gsea_2 <- sort(gsea_2, decreasing = TRUE)

gse2 <- gseGO(geneList=gsea_2,
             ont ="BP",
             keyType = "SYMBOL",
             nPerm = 10000,
             minGSSize = 5,
             maxGSSize = 800,
             pvalueCutoff = 0.99,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none")
require(DOSE)
dotplot(gse2, showCategory=10, split=".sign", font.size = 9) + facet_grid(.~.sign)