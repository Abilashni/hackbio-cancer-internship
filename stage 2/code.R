# install gplots & RColorBrewer
options(repos = c(CRAN = "https://cloud.r-project.org"))
#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")
install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)

# load dataset
gl<-read.csv("https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv", row.names = 1)
View(gl)

gl_num<-as.matrix(gl)
rownames(gl_num)<-rownames(gl)
gl_norm<-t(scale(t(log2(gl_num+1))))
rownames(gl_norm)<-rownames(gl)

#heatmaps
# Sequential Color Palette Heatmap
jpeg("plots/glioma_heatmap_sequential.jpeg", width=12, height=12, units="in", res=300, quality=100)

heatmap.2(gl_norm,
          col = brewer.pal(9, "Blues"),  
          scale = "none",
	  na.rm=TRUE,
          key = TRUE,         
          density.info = "none", 
          trace = "none",     
          margins = c(5, 10),
	  main = "Heatmap with sequential palette")

dev.off()  

# Diverging Color Palette Heatmap
jpeg("plots/glioma_heatmap_diverging.jpeg", width=12, height=12, units="in", res=300, quality=100)

heatmap.2(gl_norm,
          col = brewer.pal(9, "RdYlBu"),  
          scale = "none",
	  na.rm=TRUE,
          key = TRUE,         
          density.info = "none",  
          trace = "none",     
          margins = c(5, 10),
	  main = "Heatmap with diverging palette")

dev.off()  

# Clustering with sequentila & diverging colors

# sequential colors
# rows
jpeg("plots/seq_row_cluster.jpeg", width=12, height=12, units="in", res=300, quality=100)

heatmap.2(gl_norm,
          Rowv = TRUE,         
          Colv = FALSE,        
          col = brewer.pal(9, "Blues"),
          scale = "none",
	  na.rm=TRUE,
          key = TRUE,          
          density.info = "none",
	  dendrogram = "row",
          trace = "none",     
          margins = c(10, 10), 
	  main = "Clustering of genes")

dev.off()

# cols
jpeg("plots/seq_col_cluster.jpeg", width=12, height=12, units="in", res=300, quality=100)

heatmap.2(gl_norm,
          Rowv = FALSE,         
          Colv = TRUE,        
          col = brewer.pal(9, "Blues"),   
          scale = "none",
	  na.rm=TRUE,
          key = TRUE,          
          density.info = "none",
	  dendrogram = "column",
          trace = "none",     
          margins = c(5, 10), 
	  main = "Clustering of samples")

dev.off()

# rows & cols
jpeg("plots/seq_both_cluster.jpeg", width=12, height=12, units="in", res=300, quality=100)

heatmap.2(gl_norm,
          Rowv = TRUE,         
          Colv = TRUE,        
          col = brewer.pal(9, "Blues"),
          scale = "none",
	  na.rm=TRUE,
          key = TRUE,          
          density.info = "none",
          trace = "none",     
          margins = c(5, 10), 
	  main = "Clustering of both genes & samples")

dev.off()

# diverging colors 
# rows
jpeg("plots/div_row_cluster.jpeg", width=12, height=12, units="in", res=300, quality=100)

heatmap.2(gl_norm,
          Rowv = TRUE,         
          Colv = FALSE,        
          col = brewer.pal(9, "RdYlBu"),
          scale = "none",
	  na.rm=TRUE,
          key = TRUE,          
          density.info = "none",
	  dendrogram = "row",
          trace = "none",     
          margins = c(5, 10), 
	  main = "Clustering of genes")

dev.off()

# col
jpeg("plots/div_col_cluster.jpeg", width=12, height=12, units="in", res=300, quality=100)

heatmap.2(gl_norm,
          Rowv = FALSE,         
          Colv = TRUE,        
          col = brewer.pal(9, "RdYlBu"),   
          scale = "none",
	  na.rm=TRUE,
          key = TRUE,          
          density.info = "none",
	  dendrogram = "column",
          trace = "none",     
          margins = c(5, 10), 
	  main = "Clustering of samples")

dev.off()

# rows & cols
jpeg("plots/div_both_cluster.jpeg", width=12, height=12, units="in", res=300, quality=100)

heatmap.2(gl_norm,
          Rowv = TRUE,         
          Colv = TRUE,        
          col = brewer.pal(9, "RdYlBu"),
          scale = "none",
	  na.rm=TRUE,
          key = TRUE,          
          density.info = "none",
          trace = "none",     
          margins = c(5, 10), 
	  main = "Clustering of both genes & samples")

dev.off()

# metadata
gl_fil<-gl[rowSums(gl >= 3) >=3, ] # filter genes ::: 582-->461 genes
gl_fil<-as.matrix(gl_fil)

sample_names<-colnames(gl_fil)
condition<-factor(rep(c("A", "B"), each = length(sample_names)/2))
colData<-data.frame(row.names = sample_names, condition = condition)

dds<-DESeqDataSetFromMatrix(countData = gl_fil,
                              colData = colData,
                              design = ~ condition)
dds<-DESeq(dds)
res<-results(dds)
summary(res)

up.genes<-subset(res, log2FoldChange > 1 & padj < 0.05)
write.csv(up.genes, "upreg_genes.csv")
down.genes<-subset(res, log2FoldChange < -1 & padj < 0.05)
write.csv(down.genes, "downreg_genes.csv")

jpeg("plots/volcanoplot.jpeg", width=8, height=6, units="in", res=300, quality=100)

plot(res$log2FoldChange, -log10(res$padj),
     pch=20, 
     main="Volcano Plot",
     xlab="Log2 Fold Change", 
     ylab="-Log10 Adjusted P-value",
     col=ifelse(res$padj < 0.05 & res$log2FoldChange > 1, "red",
                 ifelse(res$padj < 0.05 & res$log2FoldChange < -1,"blue","grey")))

abline(h=-log10(0.05), col="black", lty=2)

abline(v=c(-1, 1), col="black", lty=2)

legend("topright", legend = c("Upregulated", "Downregulated", "Not Significant"),
       col = c("red", "blue", "grey"), pch = 20)

dev.off()
