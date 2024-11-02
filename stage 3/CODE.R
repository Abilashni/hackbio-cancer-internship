# loading required libraries 

library("TCGAbiolinks")
library(SummarizedExperiment)
library(edgeR)
library(gplots)
library(ggplot2)
library(biomaRt)

# project information 
getProjectSummary("TCGA-STAD")

# downloading the dataset
tcga_stad <- GDCquery(project = "TCGA-STAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification")
GDCdownload(tcga_stad) 
stad_data <- GDCprepare(tcga_stad) 
head(stad_data) 
View(stad_data)

# explore metadata information
stad_data$barcode
table(stad_data$barcode)
stad_data$race 
table(stad_data$race)
stad_data$tumor_descriptor 
table(stad_data$tumor_descriptor)
stad_data$ajcc_pathologic_stage 
table(stad_data$ajcc_pathologic_stage)
stad_data$ajcc_pathologic_m 
table(stad_data$ajcc_pathologic_m)
stad_data$gender
table(stad_data$gender)
stad_data$primary_site
table(stad_data$primary_site)
stad_data$progression_or_recurrence
table(stad_data$progression_or_recurrence)

# creating a simple metadata 
metadata_df <- data.frame("barcode" = stad_data$barcode,
                          "race" = stad_data$race,
                          'tumor_type' = stad_data$tumor_descriptor,
                          'stage' = stad_data$ajcc_pathologic_stage,
                          'metastasis_status' = stad_data$ajcc_pathologic_m,
                          'gender' = stad_data$gender)
View(metadata_df)

# categorize the metadata as female vs male 
stad_rawdata <- assays(stad_data) 
dim(stad_rawdata$unstranded) 
View(stad_rawdata$unstranded)

# reducing the data to 20 by 20  
samples <- c(subset(metadata_df, gender == "female")$barcode[c(1:20)],
             subset(metadata_df, gender == "male")$barcode[c(1:20)])
final_df <- stad_rawdata$unstranded[ , c(samples)]
dim(final_df)
View(final_df)

# clean and preprocess the data, handling missing values, normalizing gene expression data ####
table(is.na(final_df)) 
# no NAs
norm_data <- TCGAanalyze_Normalization(tabDF = final_df, geneInfo = geneInfoHT, method = "geneLength")

filt_data <- TCGAanalyze_Filtering(tabDF = norm_data,
                                   method = "quantile",
                                   qnt.cut = 0.25)

# DIFFERENTIAL EXPRESSION ANALYSIS (DEA) 
DEA <- TCGAanalyze_DEA(mat1 = filt_data[ , c(samples)[1:20]], 
                       mat2 = filt_data[ , c(samples)[21:40]],
                       Cond1type = "female",
                       Cond2type = "male",
                       pipeline = "edgeR",
                       fdr.cut = 0.05,
                       logFC.cut = 1)

DEA.Level <- 
  TCGAanalyze_LevelTab(DEA, "female", "male",
                       filt_data[ , c(samples)[1:20]],
                       filt_data[ , c(samples)[21:40]])
View(DEA.Level)

# visualization of top DEGs with a heatmap color coded based on the samples (female = green, male = red)
heat.DEGs <- filt_data[rownames(DEA.Level), ]

# color code
gender <- c(rep("female", 20), rep("male", 20))
ccodes <- c()
for (i in gender) {
  if ( i == "female") {
    ccodes <- c(ccodes, "green") 
  } else {
    ccodes <- c(ccodes, "red") 
  }
}

# generating heatmap 
heatmap.2(x = as.matrix(heat.DEGs),
          col = hcl.colors(10, palette = 'Blue-Red 2'),
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'grey',
          trace = "none",
          key = TRUE,
          cexRow = 0.9, cexCol = 0.7,
          main = "Heatmap (Females vs Males)",
          na.color = 'grey',
          ColSideColors = ccodes,
          margins = c(11,10))
# generating heatmap with dendogram
heatmap.2(x = as.matrix(heat.DEGs),
          col = hcl.colors(10, palette = 'Blue-Red 2'),
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'grey',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.9, cexCol = 0.7,
          main = "Heatmap (Females vs Males)",
          na.color = 'grey',
          ColSideColors = ccodes,
          margins = c(11,10))
legend("left",                       
       legend = c("Females", "Males"),     
       col = c("green", "red"),           
       lty = 1,                          
       lwd = 4,                         
       cex = 0.7,
       xpd = TRUE,
       inset = c(-0.2, 2))

# volcano plot
ggplot(DEA.Level, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(FDR < 0.05 & abs(logFC) > 1,
                                ifelse(logFC > 1, "Upregulated", "Downregulated"), "Not significant")), size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c("Upregulated" = "green", "Downregulated" = "red", "Not significant" = "purple")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 FDR", title = "Volcano Plot (Females vs Males)", color = "Gene Regulation") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top")

# FUNCTIONAL ENRICHMENT ANALYSIS

# selecting up- and down-regulated genes from the DEA 
upreg.genes <- rownames(subset(DEA.Level, logFC > 1 & FDR < 0.05))
dnreg.genes <- rownames(subset(DEA.Level, logFC < -1 & FDR < 0.05))

# converting ensemble IDs to gene IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg.genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = upreg.genes,
                     mart = mart)$hgnc_symbol
View(upreg.genes)
dnreg.genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = dnreg.genes,
                     mart = mart)$hgnc_symbol
View(dnreg.genes)
# EA for up- and down-regulated genes 
up.EA <- TCGAanalyze_EAcomplete(TFname = "Upregulated", upreg.genes) 
dn.EA <- TCGAanalyze_EAcomplete(TFname = "Downregulated", dnreg.genes)

# barplot for enriched pathways in up-regulated genes 
TCGAvisualize_EAbarplot(tf = rownames(up.EA$ResBP), 
                        GOBPTab = up.EA$ResBP, 
                        GOCCTab = up.EA$ResCC,
                        GOMFTab = up.EA$ResMF,
                        PathTab = up.EA$ResPat, 
                        nRGTab = upreg.genes, 
                        nBar = 5, 
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)
# barplot for enriched pathways in down-regulated genes
TCGAvisualize_EAbarplot(tf = rownames(dn.EA$ResBP), 
                        GOBPTab = dn.EA$ResBP, 
                        GOCCTab = dn.EA$ResCC,
                        GOMFTab = dn.EA$ResMF,
                        PathTab = dn.EA$ResPat, 
                        nRGTab = dnreg.genes, 
                        nBar = 5, 
                        text.size = 2, 
                        fig.width = 30,
                        fig.height = 15)

# MACHINE LEARNING MODEL 

# load required libraries and set seed (for random number generation). This ensures that the sequence of random numbers generated will be consistent each time the code is run. It ensures reproducibility of results

install.packages("caret")
install.packages("DALEX")
install.packages("DALEXtra")
install.packages("pROC")
install.packages(tidyr)

library(caret)
library(DALEX)
library(DALEXtra)
library(pROC)
library(dplyr) 
library(tidyr)

set.seed(123)
# Preparing the data for ML
stcarcinoma_data <- final_df

boxplot(stcarcinoma_data, col = "orange")
# log transformation to normalize the data
stcarcinoma_data <- log2(stcarcinoma_data + 1) 
boxplot(stcarcinoma_data, col = "orange")

# Editing the main data so it have rownames of the meta data
colnames(stcarcinoma_data) <- gsub("\\.", "-", colnames(stcarcinoma_data)) 

# Transpose main data
stcarcinoma_data <- data.frame(t(stcarcinoma_data))
# Calculates the standard deviation (SD) for each column (each gene) in stcarcinoma_data and the results are saved in the vector SDs. Higher SD shows great variability in the gene expression levels.
SDs = apply(stcarcinoma_data, 2, sd)
# Selecting the top 3000 genes with the highest SD in descending order, to give more information for further analyses.
topPredicts = order(SDs, decreasing = T)[1:3000]
stcarcinoma_data = stcarcinoma_data[, topPredicts]

# Preparing the meta data
stcarcinoma_meta <-  metadata_df

anyNA(stcarcinoma_meta)
sum(is.na(stcarcinoma_meta))
stcarcinoma_meta <- stcarcinoma_meta %>% drop_na()
anyNA(stcarcinoma_meta)
#changed the rownames from numerical "1, 2, 3, 4, 5...." to the Barcode
rownames(stcarcinoma_meta) <- stcarcinoma_meta$barcode 
#to remove the row name duplicate
stcarcinoma_meta$barcode <- NULL 

# Merging both main and meta data
stcarcinoma_merged_data <- merge(stcarcinoma_data, stcarcinoma_meta, by = "row.names")
dim(stcarcinoma_merged_data)
View(stcarcinoma_merged_data)
# make row names the samples
rownames(stcarcinoma_merged_data) <- stcarcinoma_merged_data$Row.names 
# remove duplicate columns of row names
stcarcinoma_merged_data$Row.names <- NULL 

# remove near zero variation
all_zero <- preProcess(stcarcinoma_merged_data, method = "nzv", uniqueCut = 15)
stcarcinoma_merged_data <- predict(all_zero, stcarcinoma_merged_data)
dim(stcarcinoma_merged_data)

# centering merged data to stabilize numerical computations, this is important for algorithms sensitive to the scale of the data (like KNN).
all_center <- preProcess(stcarcinoma_merged_data, method = "center")
stcarcinoma_merged_data <- predict(all_center, stcarcinoma_merged_data)
dim(stcarcinoma_merged_data)

# remove highly correlated values between features to prevent high standard errors and difficulties in interpreting coefficients, in model training and improve the efficiency and performance of certain machine learning algorithms.
all_corr <-preProcess(stcarcinoma_merged_data, method = "corr", cutoff = 0.5)
stcarcinoma_merged_data <- predict(all_corr, stcarcinoma_merged_data)
dim(stcarcinoma_merged_data)

# Preparing the data for the ML model
# Assuming data is loaded in the variable stcarcinoma_merged_data
stdata <- stcarcinoma_merged_data
# One-hot encoding the gender (convert gender into 0 and 1)
stdata$gender <- ifelse(stdata$gender == "female", 0, 1)
# Select only the gene expression columns (assuming columns starting with 'ENSG' are gene expression data)
gene_expression_columns <- grep("^ENSG", colnames(stdata), value = TRUE)

X <- stdata[, gene_expression_columns]
y <- stdata$gender
# Split the data into training and testing sets
set.seed(123)  
trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Conducting KNN model
knn_model <- train(X_train, as.factor(y_train), method = "knn", tuneLength = 5)
# Predict on the test set
y_pred <- predict(knn_model, X_test)
x_pred <- predict(knn_model, X_train)

# Evaluating model using confusion matrix
conf_matrix <- confusionMatrix(y_pred, as.factor(y_test))
conf_matrix1 <- confusionMatrix(x_pred, as.factor(y_train))
print(conf_matrix)
print(conf_matrix1)

# Perform permutation importance
set.seed(123) 
importance <- varImp(knn_model, scale = FALSE)
print(importance)
