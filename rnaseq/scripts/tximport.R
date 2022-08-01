library(tximport)

args <- commandArgs(trailingOnly = TRUE)

gse <- paste("kallisto", args[1], sep="/")
design <- read.table(args[2], sep="\t", header=TRUE)

files <- file.path(gse, design$name, "abundance.h5")
names(files) <- design$name

tx2gene <- read.table("annotation/mus_musculus/transcripts_to_genes.txt", sep="\t", header=FALSE)
colnames(tx2gene) <- c("TXNAME", "GENEID", "GENESYMBOL")

txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)

# I am giving up trying to convert a data frame to integer from float. FUCK YOU R
write.table(txi.kallisto$counts, args[3], sep="\t", quote=F)