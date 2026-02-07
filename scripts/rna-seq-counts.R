# Install necessary libraries
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)

# Utility functions
rowMax <- function(x) apply(x,1,max)

# Read data
fileList <- list.files(pattern = "\\.cts.txt$")
metadata <- read.delim("keratoconus_design.txt", header = TRUE)

# Combine all available counts to single dataframe
countList <- lapply(fileList, function(file) {
  
  currentCount <- read.delim(file, header = TRUE, comment = "#")
  currentCount <- currentCount[, c(1, ncol(currentCount))]
  colnames(currentCount) <- sub("^.*/", "", sub("\\..*", "", colnames(currentCount)))
  
  return(currentCount)
  
})

counts <- Reduce(function(x, y) merge(x, y, by = "Geneid"), countList)
counts <- as.data.frame(lapply(counts, function(x) if (is.numeric(x)) round(x, digits = 0) else x))

# Prepare data for DESeq
rownames(counts) <- counts$Geneid
genes <- counts[, c("Geneid"), drop = FALSE]
counts <- counts[colnames(counts) != "Geneid"]

# Prepare metadata for DESeq
rownames(metadata) <- metadata$SampleID
metadata <- metadata[colnames(counts),]
metadata <- metadata[colnames(metadata) != "SampleID"]
metadata$SampleGroup <- as.factor(metadata$SampleGroup)

# DESeq
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~SampleGroup)

# Ignore the following:
# Genes with no sample over threshold reads
# Samples with low reads
dds <- dds[rowMax(counts(dds)) > 30,]
dds <- dds[colSums(counts(dds)) > 1000000]

if (nrow(counts(dds)) < 1) {
  print(paste("Samples are filtered if there is < 1M reads.  There are less than no remaining sample(s) after this filter",sep=' '))
  q()
}

# Run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast = c("SampleGroup", levels(metadata$SampleGroup)))


if (nrow(metadata) > 50) {
  
  vsd <- vst(dds, blind = FALSE)
  
} else {
  
  vsd <- rlog(dds, blind = FALSE)
  
}

topGenes <- head(row.names(res[order(res$padj),]), 20)

# Visualization

# Volcano plot
png(file = "volcano.png", height = 768, width = 1024)
EnhancedVolcano(res, lab = rownames(res), x = "log2FoldChange", y = 'pvalue')
dev.off()

# Compare samples using PCA
png(file = "pca.png", height = 768, width = 1024)
plotPCA(vsd, intgroup = "SampleGroup")
dev.off()

# Heatmap of samples and genes
x <- colnames(assay(vsd))
y <- rownames(assay(vsd))
data <- expand.grid(X=x, Y=y)
data$Z <- as.vector(t(assay(vsd)))
data$group <- rep(metadata$SampleGroup, times = length(rownames(assay(vsd))))

data <- data %>% 
  mutate(text = paste0("Sample: ", x, "\n", "Gene: ", y, "\n",
                                      "Normalized Gene Count: ", round(Z,3), "\n",
                                      "Sample Group: ", group)) %>%
  filter(Y %in% topGenes)

png(filename = "heatmap.png", width = 768, height = 1024, res = 100)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()
dev.off()

# Pairwise comparison
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- vsd$SampleGroup

png(filename = "distance_heatmap.png", width = 768, height = 1024, res = 100)
pheatmap(sampleDistMatrix, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

# GO Analysis
ensemblIds <- sub("\\..*$", "", rownames(res))
entrezIds <- mapIds(org.Hs.eg.db,
                    keys = ensemblIds,
                    column = "ENTREZID",
                    keytype = "ENSEMBL",
                    multiVals = "first")

deGenesEnsembl <- (res$padj < 0.1) & (abs(res$log2FoldChange) > 1)
deGenesEntrez <- unname(entrezIds[deGenesEnsembl])
deGenesEntrez <- deGenesEntrez[ ! is.na(deGenesEntrez) ]

geneUniverse <- unname(entrezIds[ensemblIds])

go <- goana(de = deGenesEntrez,
            universe = geneUniverse,
            species = "Hs")

# Multiple testing correction with BH
go <- go[order(go$P.DE), ]
go$FDR <- p.adjust(go$P.DE, method = "fdr")

go <- cbind(IDs = rownames(go), go)
write.csv(go, file = "GO_analysis.csv", row.names = FALSE)
