library(DESeq2)
count <- read.csv("Normalize_count.csv", header = TRUE, row.names = 1)
condition <- factor(c("A","A","A","B","B","B"))
dds <- DESeqDataSetFromMatrix(round(count), DataFrame(condition), ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res
dds <- lfcShrink(dds, coef=2)
library(BiocParallel)
register(MulticoreParam(4))
resOrder <- dds[order(dds$pvalue), ] 
summary(resOrder)
sum(dds$padj < 0.1, na.rm = TRUE)
plotMA(dds)
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')
