## DESeq2 script for RNA-Seq analysis

library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
#library(ggfortify)
#library(devtools)
#library(limma)
#library(FactoMineR)

args = commandArgs(trailingOnly = TRUE)

compare_col <- args[1]
compare_char1 <- args[2]
compare_char2 <- args[3]

## Variables
# Directory with htseq results
directory <- "../../07-featureCounts/"
# File names present in directory
sampleFiles <- grep("_gene.featureCounts.txt.summary",list.files(directory,recursive=TRUE,include.dirs=FALSE),value=TRUE)
# Get only sampleNames
sampleNames <- sub("(.*)/(.*)_gene.featureCounts.txt.summary","\\2",sampleFiles)
# Metadata obtained from user table.
# state
pasCts <- "../../07-featureCounts/merged_gene_counts.txt"

pasAnno <- "../../clinical_data.txt"

cts <- as.matrix(read.table(pasCts,sep="\t", header = T))
header_with_period <- colnames(cts)
header_with_underscore <- gsub(x = header_with_period, pattern = "\\.", replacement = "_")
colnames(cts) <- header_with_underscore
gene_ids <- cts[,1]
rownames(cts) <- gene_ids

coldata <- read.table(pasAnno, sep = "\t", header = T, row.names = 1)

#all(rownames(coldata) %in% colnames(cts)) #test if the have same name
cts <- cts[, rownames(coldata)] #order the colums of cts to be the same as in coldata
#all(rownames(coldata) == colnames(cts)) #test if previous step worked
cts <- as.matrix(cts)
mode(cts) <- "integer" #convert values to integer


#####DIFERENTIAL EXPRESSION###################
dds_full <- DESeqDataSetFromMatrix(countData = cts,
                                   colData = coldata,
                                   design = formula(paste("~",compare_col)))
featureData <- data.frame(gene=rownames(cts))
mcols(dds_full) <- DataFrame(mcols(dds_full), featureData)
dds <- dds_full[ rowSums(counts(dds_full)) >= 1, ]
dds <- DESeq(dds)
res <- results(dds,contrast = c(compare_col,compare_char1,compare_char2))
mcols(res, use.names=TRUE)

DE_results <- as.data.frame(res)
DE_results$GeneID <- row.names(DE_results)
DE_results <- DE_results[,c(ncol(DE_results), 1:ncol(DE_results)-1)]
write.table(x = DE_results, file = "DE_matrix.txt", sep = "\t", quote = F, col.names = T, row.names = F)

ntd <- normTransform(dds)
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

## Plots Quality control
#MA-plotThe  MA-plot  shows  the  log2  fold  changes  from  the  treatment  over  the  meanof  normalized  counts.
#The  average  of  counts  normalized  by  size  factor.
pdf(file="maPlot_all.pdf")
plotMA( res, ylim = c(-1, 1) )
dev.off()


###########SAMPLE DISTANCE##############
sampleDists <- dist( t( assay(rld) ) )

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- names(rld$sizeFactor)
colnames(sampleDistMatrix) <- names(rld$sizeFactor)
colours = colorRampPalette(rev(brewer.pal(9, "Blues"))) (255)
pdf(file="heatmap_sample_to_sample.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colours)
dev.off()


#############PCA PLOTS################
pcaData <- plotPCA(rld, intgroup=c(compare_col), returnData=TRUE)
pcaData_2 <- plotPCA(vsd, intgroup=c(compare_col), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file="plotPCA.pdf")
ggplot(pcaData, aes_string("PC1", "PC2", color=compare_col)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label = name), color = "black", size=2, position = position_nudge(y = 0.8)) +
  labs(title="PCA: rlog") +
  coord_fixed()
ggplot(pcaData_2, aes_string("PC1", "PC2", color=compare_col)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label = name), color = "black", size=2, position = position_nudge(y = 0.8)) +
  labs(title="PCA: vsd") +
  coord_fixed()
dev.off()

#############BOX PLOTS################
pdf(file="boxplot.pdf")
boxplot(assay(ntd), col="blue", las =2)
title(main="Boxplot: normalized counts")
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
title(main="Boxplot see outliers: cooks distance")
dev.off()

#############DISPERSION PLOTS################
pdf(file="plotDispersions.pdf")
plotDispEsts(dds)
title(main="Dispersion calc")
hist( res$pvalue, breaks=20, col="grey", main = "pvalues test for differential expression")
dev.off()


#############DESVIATION PLOT################
pdf(file="plotSD.pdf")
meanSdPlot(assay(ntd))
dev.off()


##############HCLUST###################
assay_ntd <- assay(ntd)
pdf(file="cluster_dendrogram.pdf")
plot(hclust(dist(t(assay_ntd)),method="average"))
dev.off()


##############PHEATMAP##############
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c(compare_col)])
colnames(df) <- c(compare_col)
rownames(df) <- colnames(ntd)
pdf(file="heatmapCount_top20.pdf")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, main="Normalized counts top 20 genes")
dev.off()



########################################################################
######FULL PHEATMAP#################
pdf(file="heatmapCount_all_genes.pdf")
pheatmap(assay(ntd), cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE,main="Normalized counts", annotation_col=df)
dev.off()

save.image()
if(FALSE) {
########NEW PCA info###############
ntop=500
returnData=FALSE
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(assay(rld)[select,]))

###Manually Plot#####
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
intgroup = "samples"
intgroup.df <- as.data.frame(colData(rld)[, intgroup, drop=FALSE])
group <- if (length(intgroup) > 1) {
  factor(apply( intgroup.df, 1, paste, collapse=":"))
} else {
  colData(rld)[[intgroup]]
}
d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(rld))
if (returnData) {
  attr(d, "percentVar") <- percentVar[1:2]
  return(d)
}

ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  labs(title="PCA: rlog") +
  geom_text(aes_string(x = "PC1", y = "PC2", label = "samples"), color = "black", size=2, position = position_nudge(y = 2)) +
  theme(legend.position = "none") +
  coord_fixed()

###END Plot#####

axes <- (1:length(pca$sdev))
proba = 0.05
new_res <- NULL
new_res$ind$coord <- as.data.frame(pca$x)
new_res$call$X <- t(assay(rld)[select,])
new_res$call$row.w.init <- rep(1, length(pca$sdev))
result = structure(vector(mode = "list", length = length(axes)), names = colnames(new_res$ind$coord)[axes])
for (k in 1:length(axes)) {
  tableau = cbind.data.frame(new_res$ind$coord[, axes[k], drop = FALSE], new_res$call$X)
  result[[k]] <- condes(tableau, 1, proba = proba, weights = new_res$call$row.w.init)
}

gene_parser <- read.table(file = "../gene_parser.txt", header = T, sep = "\t")

dimension_1 <- na.omit(result$PC1[["quanti"]])
genes_D1 <- rownames(dimension_1)
genes_D1 <- as.data.frame(genes_D1)
genes_D1 <- merge(genes_D1, gene_parser, by.x="genes_D1", by.y="Geneid", all.x=T, all.y=F)
genes_D1 <- as.vector(genes_D1$gene_name)
genes_D1 <- unique(genes_D1)
write(x = genes_D1, file = "PC1_genes.txt", sep = "\n")

dimension_2 <- na.omit(result$PC2[["quanti"]])
genes_D2 <- rownames(dimension_2)
genes_D2 <- as.data.frame(genes_D2)
genes_D2 <- merge(genes_D2, gene_parser, by.x="genes_D2", by.y="Geneid", all.x=T, all.y=F)
genes_D2 <- as.vector(genes_D2$gene_name)
genes_D2 <- unique(genes_D2)
write(x = genes_D2, file = "PC2_genes.txt", sep = "\n")

save.image()
}
