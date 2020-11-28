library(tximport)
library(readr)
library(dplyr)
library(jsonlite)
library(EnsDb.Hsapiens.v86)

#Annotation tables can be found on the gencode website
t10b <- read.table("../index/gencodev27.annotation.txt", stringsAsFactors=FALSE, col.names ="Name")
TranscriptID <- t10b[,1]
Transcripts <- sapply(strsplit(TranscriptID,'|',fixed=TRUE), "[[", 1)
GeneID <- sapply(strsplit(TranscriptID,'|',fixed=TRUE), "[[", 6)
RNAClass <- sapply(strsplit(TranscriptID,'|',fixed=TRUE), "[[", 8)
tx2gene <- data.frame(Transcripts, GeneID, RNAClass)

# Read which samples to include
samples <- list.files(pattern = "_quant")
salmon.files = file.path(fsep="", samples, "/quant.sf")
names(salmon.files) = samples

# Check all exist
salmon.files[!file.exists(salmon.files)] #name character(0) means yes

# import
salmon.tx = tximport(salmon.files, type="salmon", tx2gene=tx2gene, countsFromAbundance="no",ignoreAfterBar=TRUE)
genes_all = rownames(salmon.tx$abundance)
salmon.counts = salmon.tx$abundance
write.table(salmon.counts, file="your.salmon.counts.tsv",sep='\t',quote=FALSE,row.names=TRUE)

###save it as a file
save(salmon.tx, file = "your.salmon.counts.RData")
