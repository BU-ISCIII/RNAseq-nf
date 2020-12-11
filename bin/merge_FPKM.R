#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

require(ballgown)
my.data = ballgown(dataDir="./", samplePattern='*_ballgown', meas='all')
expression <- gexpr(my.data)
expression <- as.data.frame(expression)
expression$gene <- rownames(expression)
expression <- expression[, c(ncol(expression), 1:c(ncol(expression)-1))]
write.table(x = expression, file = "merged_FPKM.txt", col.names = T, row.names = F, sep = "\t", quote = F)
