# this project is on Cancer biomarker discovery and Machine Learning of Gliomas

# Project: Biomarker Discovery and Machine Learning of Gliomas
# Cancer Subtype: Low-grade gliomas (LGGs) and Glioblastoma multiforme (GBM)

# setting working directory
getwd()
setwd("C:/Users/nigus/HackBio/Raw_Data")
getwd()

# enable libraries
library("TCGAbiolinks")
library(SummarizedExperiment) #install alongside with TCGAbiolinks
library(biomaRt)
library("edgeR")
library("limma")
library("EDASeq")
library(gplots)
library(ggplot2)
library(ggrepel)

# project information 
getProjectSummary("TCGA-LGG")
getProjectSummary("TCGA-GBM")

# download LGG and GBM datasets from TCGA
# download LGG datasets
tcga_lgg <- GDCquery(project = "TCGA-LGG",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification")

GDCdownload(tcga_lgg) 
lgg_data <- GDCprepare(tcga_lgg) 

# download GBM datasets
tcga_gbm <- GDCquery(project = "TCGA-GBM",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification")
GDCdownload(tcga_gbm) 
gbm_data <- GDCprepare(tcga_gbm) 

# get raw counts from lgg_data
lgg_rawdata <- assays(lgg_data) 
dim(lgg_rawdata$unstranded) 

# get raw counts from gbm_data
gbm_rawdata <- assays(gbm_data) 
dim(gbm_rawdata$unstranded) 

# prepare metadata for lgg
lgg_data$barcode
lgg_data$paper_IDH.status
table(lgg_data$paper_IDH.status) 

meta_lgg <- data.frame("barcode" = lgg_data$barcode,
                       "IDH status" = lgg_data$paper_IDH.status)
View(meta_lgg)
table(is.na(meta_lgg)) 
meta_lgg <- na.omit(meta_lgg)

# prepare metadata for gbm 
gbm_data$barcode
gbm_data$paper_IDH.status
table(gbm_data$paper_IDH.status) 

meta_gbm <- data.frame("barcode" = gbm_data$barcode,
                       "IDH status" = gbm_data$paper_IDH.status)
View(meta_gbm)
table(is.na(meta_gbm))
meta_gbm <- na.omit(meta_gbm)

# final dataset for lgg (513 samples)
lgg <- lgg_rawdata$unstranded
samples_lgg <- meta_lgg$barcode
lgg <- lgg[, colnames(lgg) %in% samples_lgg]
all(colnames(lgg) %in% samples_lgg) # TRUE

# final dataset for gbm (154 samples)
gbm <- gbm_rawdata$unstranded
samples_gbm <- meta_gbm$barcode
gbm <- gbm[, colnames(gbm) %in% samples_gbm]
all(colnames(gbm) %in% samples_gbm)

# merge the two datasets (both counts and metadata)
glioma <- cbind(lgg, gbm)
meta_glioma <- rbind(meta_lgg, meta_gbm)

# export the final dfs
write.csv(glioma, "glioma_rawcounts.csv", row.names = TRUE)
write.csv(meta_glioma, "glioma_metadata.csv", row.names = FALSE)

# clean and preprocess the data
table(is.na(glioma))
glioma <- TCGAanalyze_Normalization(tabDF = glioma, geneInfo = geneInfoHT, method = "geneLength")
glioma <- TCGAanalyze_Normalization(tabDF = glioma, geneInfo = geneInfoHT, method = "gcContent")
glioma <- TCGAanalyze_Filtering(tabDF = glioma,
                                method = "quantile",
                                qnt.cut = 0.25)

# differential expression analysis 

group.wt <- meta_glioma$barcode[meta_glioma$IDH.status == "WT"]
group.mut <- meta_glioma$barcode[meta_glioma$IDH.status == "Mutant"]

DEA <- TCGAanalyze_DEA(mat1 = glioma[ , group.wt], 
                       mat2 = glioma[ , group.mut],
                       Cond1type = "WT",
                       Cond2type = "Mutant",
                       pipeline = "edgeR")
DEA.Level <- 
  TCGAanalyze_LevelTab(DEA, "WT", "Mutant",
                       glioma[ , group.wt],
                       glioma[ , group.mut])

# annotation
gene_names <- rownames(DEA.Level)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = gene_names,
  mart = mart)
gene_with_symbols <- merge(DEA.Level, annot, 
                           by.x = "mRNA", 
                           by.y = "ensembl_gene_id", 
                           all.x = TRUE)

# volcano plot
genes_to_label <- subset(gene_with_symbols, abs(logFC) > 5 & FDR < 0.05)

ggplot(DEA.Level, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(FDR < 0.05 & abs(logFC) > 1,
                                ifelse(logFC > 1, "Upregulated", "Downregulated"), "Not significant")), size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c("Upregulated" = "green", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 FDR", title = "Volcano Plot of WT vs Mut", color = "Gene Regulation") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top")
  #ggrepel::geom_label_repel(data = genes_to_label, aes(label = hgnc_symbol))


# heatmap of DEGs
# top differentially expressed genes 
DEGs <- subset(DEA.Level, abs(logFC) > 3 & FDR < 0.01) # 905 DEGs   
heat.DEGs <- glioma[rownames(DEGs), ]

# color code
ccodes <- ifelse(colnames(heat.DEGs) %in% group.wt, "green", "blue") # wt in red, mut in blue

# heatmap
red_green_palette <- colorRampPalette(c("green", "black"))(10)
heatmap.2(x = as.matrix(heat.DEGs),
          col = red_green_palette,
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.9, cexCol = 0.7,
          main = "Heatmap of WT vs Mut",
          na.color = 'black',
          ColSideColors = ccodes,
          margins = c(11,10))
legend("left",                       
       legend = c("WT", "Mut"),     
       col = c("green", "blue"),           
       lty = 1,                          
       lwd = 4,                         
       cex = 0.7,
       xpd = TRUE,
       inset = c(-0.1, 0))



########## Machine Learning ####################

# load necessary Libraries
library(caret)
library(DALEX)
library(pROC)
library(iml)
set.seed(1234)

# Load the main and meta data
glioma.data <- read.csv(file = "glioma_rawcounts.csv", header = TRUE)
glioma.meta <- read.csv(file = "glioma_metadata.csv", header = TRUE)
dim(glioma.data)
rownames(glioma.data) <- glioma.data$X

glioma.data$X <- NULL

colnames(glioma.data) <- gsub("\\.", "-", colnames(glioma.data))
glioma.data <- log10(glioma.data + 1)

glioma.data <- data.frame(t(glioma.data))


SDs = apply(glioma.data, 2, sd)
# This line calculates the standard deviation for each column (each gene) in lgg.data and stores the results in the vector SDs. 
# A higher standard deviation indicates greater variability in the gene expression levels.

topPredicts = order(SDs, decreasing = T)[1:3000]
#After calculating the standard deviation, this lines ensures that the top 3000 genes with the highest SD are selected, and in descending order.


glioma.data = glioma.data[, topPredicts]
#The glioma.data is now reduced to only include the top 3000 genes that have the highest variability, 
# which are often more informative for subsequent analyses like classification

# To remove near zero variation
zer <- preProcess(glioma.data, method = "nzv", uniqueCut = 15)
glioma.data <- predict(zer, glioma.data)

# center. Centering helps to stabilize numerical computations 
#and can be particularly important for algorithms sensitive to the scale of the data (like PCA, KNN, etc.).
cent <- preProcess(glioma.data, method = "center")
glioma.data <- predict(cent, glioma.data)

# to remove highly correlated values.
#Reduce Multicollinearity: High correlations between features can lead to issues in model training, 
#such as inflated standard errors and difficulties in interpreting coefficients. 
cor <- preProcess(glioma.data, method = "corr", cutoff = 0.5)
glioma.data <- predict(cor, glioma.data)


glioma.data <- scale(glioma.data)  # Scaling the data

# Perform K-Means clustering
set.seed(123)
kmeans_result <- kmeans(glioma.data, centers = 4, nstart = 25)

fviz_cluster(kmeans_result, data = glioma.data)



# Create a dataframe from scaled data
scaled_data_df <- as.data.frame(glioma.data)

# Add cluster assignments to the dataframe
scaled_data_df$Cluster <- kmeans_result$cluster

# View the updated dataframe
head(scaled_data_df)


## Process the metadata
rownames(glioma.meta) <- glioma.meta$barcode

glioma.meta$barcode <- NULL

merge.data <- merge(scaled_data_df, glioma.meta, by = "row.names")

rownames(merge.data) <- merge.data$Row.names


fviz_cluster(kmeans_result, 
             data = glioma.data,  # Use the original scaled data
             geom = "point",       # You can also use "text" to label points
             pointsize = 2,        # Adjust point size for better visibility
             ggtheme = theme_minimal()) +
  # Add colors based on IDH status from the combined_data
  geom_point(aes(color = combined_data$IDH.status)) + 
  labs(color = "IDH Status") +  # Label for the legend
  scale_color_manual(values = c("Mutant" = "blue", "WT" = "red"))