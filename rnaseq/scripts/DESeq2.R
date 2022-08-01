library(BiocParallel)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
register(MulticoreParam(16))

args <- commandArgs(trailingOnly = TRUE)

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]

deseq_function <- function(counts_file, design_file, threshold, out_prefix){
  counts = read.csv(counts_file, sep="\t", header = TRUE,
                    row.names = 1, check.names = FALSE)
  design = read.csv(design_file, header=TRUE, sep="\t", row.names=1)
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = design,
                                design = ~ condition)
  
  dds <- dds[rowSums(counts(dds)) > threshold,]
  
  # Performing DESeq2 analysis
  dds <- DESeq(dds, parallel=TRUE)
  saveRDS(dds, file=paste(out_prefix, "dds.rds", collapse="", sep=""))
  rld <- rlog(dds)

  pcaPlotData <- plotPCA(rld, intgroup=c("condition"), returnData=T)  
  percentVar <- round(100 * attr(pcaPlotData, "percentVar"))
  pcaPlot <- ggplot(pcaPlotData, aes(PC1, PC2)) + 
    geom_point(aes(colour=condition), alpha=0.9) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) 
  
  theme_set(theme_bw(base_size=15, base_family="Source Sans"))
  pcaPlot
  ggsave(paste(out_prefix, "PCA.png", collapse="", sep=""), width=5,
         height=3.5, dpi=300, units="in")
  
  write.table(assay(rld), file=paste(c(out_prefix, "rld.tsv"), collapse="", sep=""), sep="\t")

  ko_vs_wt <- results(dds, c("condition", "KO", "WT"), independentFiltering = TRUE)

  write.table(as.data.frame(ko_vs_wt),
      file=paste( out_prefix, "ko_vs_wt.tsv", collapse = "", sep=""),
      quote=F, col.names=NA, sep="\t")
  
  print(paste(c(counts_file, "finished")))
}

deseq_function(args[1], args[2], as.integer(args[3]), args[4])