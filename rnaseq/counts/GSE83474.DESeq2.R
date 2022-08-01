library(BiocParallel)
library(DESeq2)
register(MulticoreParam(4))

deseq_function <- function(counts_file, design_file, threshold, out_prefix){
    
    counts = read.csv(counts_file, sep="\t", header = TRUE, row.names = 1, check.names = FALSE)
    design = read.csv(design_file, header=TRUE, sep="\t", row.names=1)
    
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~ condition)
    dds <- dds[rowSums(counts(dds)) > threshold,]

    # Performing DESeq2 analysis
    dds <- DESeq(dds, parallel=TRUE)    
    saveRDS(dds, file=paste(out_prefix, "dds.rds", collapse="", sep=""))
    
    rld <- rlog(dds)
    r106w_vs_wt <- results(dds, c("condition", "R106W", "WT"), independentFiltering = TRUE)
    write.table(as.data.frame(r106w_vs_wt), file=paste( out_prefix, "R106W_vs_WT.tsv", collapse = "", sep=""), quote=F, col.names=NA, sep="\t")
    
    t158m_vs_wt <- results(dds, c("condition", "T158M", "WT"), independentFiltering = TRUE)
    write.table(as.data.frame(t158m_vs_wt), file=paste( out_prefix, "T158M_vs_WT.tsv", collapse = "", sep=""), quote=F, col.names=NA, sep="\t")
    print(paste(c(counts_file, "finished")))
}

deseq_function("GSE83474.tsv", "GSE83474.design.tsv", 20, "GSE83474.DESeq2.results/")