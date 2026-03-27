library(DESeq2)
library(ggplot2)
library(dplyr)

counts <- read.table(
  "GSE152418_p20047_Study1_RawCounts.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

head(counts)
dim(counts)

#metadata
condition <- c(
  rep("COVID",17),
  rep("Healthy",17)
)

coldata <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

coldata

#create deseq dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)
dds<-dds[rowSums(counts(dds))>10,]

#deg analysis
dds <- DESeq(dds)

res <- results(dds)

summary(res)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

#upregulated genes
upregulated <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1)

dim(upregulated)
head(upregulated)

write.csv(
  upregulated,
  "Upregulated_genes.csv",
  row.names = FALSE
)
nrow(upregulated)

#downregulated genes
downregulated <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < -1)

dim(downregulated)
head(downregulated)

write.csv(
  downregulated,
  "Downregulated_genes.csv",
  row.names = FALSE
)
nrow(downregulated)

#deg results
write.csv(
  res_df,
  "All_DEG_results.csv",
  row.names = FALSE
)
nrow(res_df[res_df$padj<0.05&abs(res_df$log2FoldChange)>1,])
nrow(res_df[!res_df$padj<0.05&abs(res_df$log2FoldChange)>1,])
table(res_df$gene_type)


#volcano plots
library(ggplot2)

res_df <- as.data.frame(res)

res_df$gene_type <- "Non-significant"
res_df$gene_type[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$gene_type[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

volcano_plot <- ggplot(res_df,
                       aes(log2FoldChange, -log10(padj), color=gene_type)) +
  geom_point(alpha=0.8, size=1.5) +
  scale_color_manual(values=c(
    "Upregulated"="red",
    "Downregulated"="blue",
    "Non-significant"="grey"
  )) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_minimal() +
  labs(title="Volcano Plot")



# save image
ggsave("Volcano_plot.png", volcano_plot, width=8, height=6)

#ma plot
library(ggplot2)

# Convert results to dataframe
res_df <- as.data.frame(res)

# Remove NA values
res_df <- res_df[!is.na(res_df$padj), ]

# Create gene categories
res_df$gene_type <- "Non-significant"

res_df$gene_type[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$gene_type[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# MA plot
ma_plot <- ggplot(res_df,
                  aes(x=baseMean,
                      y=log2FoldChange,
                      color=gene_type)) +
  
  geom_point(alpha=0.7, size=1.5) +
  
  scale_x_log10() +
  
  scale_color_manual(values=c(
    "Upregulated"="red",
    "Downregulated"="blue",
    "Non-significant"="grey"
  )) +
  
  geom_hline(yintercept=c(-1,1), linetype="dashed") +
  
  theme_minimal() +
  
  labs(
    title="MA Plot (Healthy vs Infected)",
    x="Mean Expression (log scale)",
    y="Log2 Fold Change",
    color="Gene Category"
  )



# Save image
ggsave("MA_plot_colored.png", ma_plot, width=8, height=6)

#heatmaps
library(pheatmap)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Order genes by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Remove NA padj
resOrdered <- resOrdered[!is.na(resOrdered$padj), ]

# Select top 20 genes
topgenes <- head(rownames(resOrdered), 20)

# Extract expression matrix
mat <- assay(vsd)[topgenes, ]

# Remove zero variance genes
mat <- mat[apply(mat,1,sd) != 0,]

# Scale rows
mat <- t(scale(t(mat)))

# Remove NA rows if any
mat <- mat[complete.cases(mat),]

# Create annotation for sample condition
annotation_col <- data.frame(
  Condition = colData(dds)$condition
)

rownames(annotation_col) <- colnames(mat)

# Heatmap
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  main = "Top 20 Differentially Expressed Genes (Healthy vs Infected)"
)

#pca plots
library(ggplot2)

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  theme_minimal() +
  labs(
    title="PCA Plot (Healthy vs Infected)",
    x=paste0("PC1: ", percentVar[1], "% variance"),
    y=paste0("PC2: ", percentVar[2], "% variance")
  )

# show plot
pca_plot

# save image
ggsave("PCA_plot.png", pca_plot, width=8, height=6)


up <- read.csv("Upregulated_genes.csv")
down <- read.csv("Downregulated_genes.csv")

# Assuming gene symbols are in column "gene"
write.table(up$gene, "up_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(down$gene, "down_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)