#required libraries
library(DESeq2)
library(tximport)
library(ggplot2)
library(vsn)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library("genefilter")
#importing Gene ID's for later
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
edb <- EnsDb.Hsapiens.v79

#########################################################
#########################################################

### PART 1 TRANSCRIPT LVL ANALYSIS ###

# optional portion of the script for collapsing transcript reads to a gene level matrix. If tximport script is used, start at part 2

metadataTC <- read.table('metadata/exp_cond_hour.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)
files <- file.path("quants_hours", metadataTC$sample, "abundance.tsv")
names(files) <- paste0(metadataTC$sample)
txi.kallisto.tsvTC <- tximport(files, type = "kallisto", txOut = T, ignoreAfterBar = T)
head(txi.kallisto.tsvTC$counts)

ddsTC <- DESeqDataSetFromTximport(txi.kallisto.tsvTC,
                                  colData = metadataTC,
                                  design = ~line+hour)

#converting hour/line column to factor
ddsTC$hour <- factor(ddsTC$hour, levels=c("0","4","12","18"))
ddsTC$line <- factor(ddsTC$line, levels=c("1","2","3"))

#throw out low count samples
keep <- rowSums(counts(ddsTC)) >= 1
ddsTC <- ddsTC[keep,]
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~line)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)

#importing gene ID's
seqlevelsStyle(edb) <- "UCSC"

ens.str <- substr(rownames(resTC),1,15)
resTC$symbol <- mapIds(edb,
                       keys=ens.str,
                       column="SYMBOL",
                       keytype="TXID",
                       multiVals="first")
resTC$entrez <- mapIds(EnsDb.Hsapiens.v79,
                       keys=ens.str,
                       column="ENTREZID",
                       keytype="TXID",
                       multiVals="first")


#set rownames to gene symbols for resTC
resTCgenes<-resTC[complete.cases(resTC[ , "symbol"]),]
rownames(resTCgenes) <- resTCgenes[["symbol"]]
#set rownames to gene symbols for ddsTC for heatmap
ddsTCgenes<-coef(ddsTC)
ddsTCgenes<-coef(ddsTC[complete.cases(resTC[ , "symbol"]),])
rownames(ddsTCgenes)<-resTCgenes[["symbol"]]

#Timecourse single gene plot
fiss <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup = c("hour","line"), returnData = TRUE)
fiss$hour <- as.numeric(as.character(fiss$hour))
ggplot(fiss,
       aes(x = hour, y = count, color = line, group = line)) +
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()

#GOI list
goi <- c("ATF4", "XBP1", "ERN1", "NFKB1", "DDIT3", "RIPK1", "TNFRSF10B", "HSPA5")
gois <-which(rownames(resTCgenes) %in% goi)

#TC heatmap
betas <- ddsTCgenes
colnames(betas)
topGenes <- gois
mat <- betas[topGenes,c(2:6)]
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)


#PCA plot for biological replicates
plotPCA(rld, intgroup = c("hour", "line"))

#generate "distances" between samples
sampleDists <- dist(t(assay(rld)))

#plot sample distance matrix
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition, rld$sample, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#export count matrix into wd for gene level analysis
csv.write(txi.kallisto.tsvTC, "tx_lvl_counts_TC.csv")

#for use in PART 2, you must aggregate by row for gene names.
#I did this in google sheets using a pivot table adding values and summarizing by gene, but it can be done in R as well.


#################################################################
#################################################################

### PART 2 GENE LVL ANALYSIS ###

#Input should be a raw count matrix as produced by the tximport script, or by step 1

#import gene_lvl counts table
setwd('/Users/tsears/code/NereaZappaProj/')
cts <- as.matrix(read.csv("gene_experiment_aggregated_xbp1.csv", sep=",",header=T))
anno <- read.table('exp_cond_hour_genelvl.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)

#prepare counts matrix
rownames(cts)<-cts[,1]
cts<-cts[,2:13]
cts<-na.omit(cts)
cts<-as.matrix(cts)
mode(cts)<-"integer"

#keep only highly expressed genes
keep <- rowMeans(cts) >= 30
cts <- cts[keep,]

#construct dds object. Design = whichever variables are included in your metadata table,
#as well as combinations of them if that applies
dds_genelvl <- DESeqDataSetFromMatrix(countData = cts,
                                      colData = anno,
                                      design = ~ line+hour)

#convert metadata integers to factors--this helps to resolve the warning message issued above
dds_genelvl$hour <- factor(dds_genelvl$hour, levels=c("0","4","12","18"))
dds_genelvl$line <- factor(dds_genelvl$line, levels=c("1","2","3"))


#run LRT on reduced version of test to see how HOUR impacts counts. Specify ~line in the reduced model.
dds_genelvl <- DESeq(dds_genelvl, test="LRT", reduced = ~line)
res_g <- results(dds_genelvl)
res_g$symbol <- mcols(dds_genelvl)$symbol
head(res_g[order(res_g$padj),], 4)

#export normalized counts table for GSEA. This is used for seperate ontology analysis in GSEA
write.table(counts(dds_genelvl, normalized=TRUE), "normalized_counts.txt",sep="\t")

#show a plot of all counts highlighting statistically significant samples
plotMA(res_g, alpha=0.1)

#filter by padj value and see which genes in genelist are significant
res0.05 <- subset(res12, (res12$padj <= (0.05))) 
res0.05 <- subset(res0.05, abs(res0.05$log2FoldChange) > 1.2)
res0.05 <- res0.05[order(abs(res0.05$log2FoldChange),decreasing = T),]
top300 <- res0.05[1:1000,]
top300 <- top300[order(top300$log2FoldChange,decreasing =T),]
top300list <- top300$log2FoldChange
names(top300list)<-rownames(top300)
top300list<-as.data.frame(top300list)
write.table(top300list,"0v12_top1000_byL2FC.rnk.txt",sep = "\t") #Quick output for inspection

sig_genes<-(which(rownames(res0.05)%in%gois))
sig_gene_names<-(res0.05[sig_genes,])
sig_genespvals<-res0.05[sig_gene_names,]
sig_genespvals<-sig_genespvals[order(sig_genespvals$padj,decreasing = F),]


#List of genes of interest for heatmap plotting
gois<-scan(text="DNAJB9 TMED2 SERP1 VIMP 
ERLEC1 ARMCX3 TIMM17A SEC61A1 SEC61B PPIB 
TMED9 C19orf10 SSR2 NANS OSTC SSR3 SSR1 MTDH YIF1A TMEM165 SLC25A3 CD59",what="")

#TimeCourse heatmap
betas <- coef(dds_genelvl)
colnames(betas)
mat <- betas[intersect(rownames(betas),gois),c(4:6)]
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, fontsize_row = 5,main="UPR All Genes",labels_col = c("4 vs 0","12 vs 0","18 vs 0"))

#graph of one GOI over time
fiss <- plotCounts(dds_genelvl, which(rownames(res_g) %in% c("TFEB")),
                   intgroup = c("hour"), returnData = TRUE)
fiss$hour <- as.numeric(as.character(fiss$hour))

ggplot(fiss,
       aes(x = hour, y = count)) +
  #geom_point() 
  stat_summary(fun=mean, geom="line") +
  scale_y_log10()

ggsave("TFEB.png")

#export statistically significant head to head LFC results
resultsNames(dds_genelvl)
res4 <- results(dds_genelvl, name="hour_4_vs_0", test="LRT")
res4 <- res4[res4$padj<0.05,]
res12 <- results(dds_genelvl, name="hour_12_vs_0", test="LRT")
res12 <- res12[res12$padj<0.05,]
res18 <- results(dds_genelvl, name="hour_18_vs_0", test="LRT")
res18 <- res18[res18$padj<0.05,]

write.csv(res4,"Secretome_4v0_gene_lvl_FC_1.5.csv")
write.csv(res12,"Secretome_12v0_gene_lvl_FC_1.5.csv")
write.csv(res18,"Secretome_18v0_gene_lvl_FC_1.5.csv")


#exporting filtered head to head GOI lists. Uncomment portions below for threshold filtering
res4_h2h <- subset(res4,
                   #(res4$log2FoldChange <= (-1.5))|(res4$log2FoldChange >= (1.5))&
                     (res4$padj <= (0.05)))
res12_h2h <- subset(res12,
                   # (res12$log2FoldChange <= (-1.5))|(res12$log2FoldChange >= (1.5))&
                      (res12$padj <= (0.05)))
res18_h2h <- subset(res18,
                   # (res18$log2FoldChange <= (-1.75))|(res18$log2FoldChange >= (1.75))&
                      (res18$padj <= (0.05)))

write.csv(res4_h2h,"Secretome_0vs4_genelist.csv")
write.csv(res12_h2h,"Secretome_0vs12_genelist.csv")
write.csv(res18_h2h,"Secretome_0vs18_genelist.csv")

####################################
####################################

### PART 3 GENE ONTOLOGY PLOTTING ###
#This portion requires that GSEA is run on the previous normalized count matrix output in step 2
library(forcats)

#Below are bubble plots as described in the README

#0V4
d1 = read.table("/Users/tsears/Downloads/Hallmark_0v4/excel_files/gsea_report_for_04_1587425425860.xls",sep="\t",header=T)
d1$NES<-abs(d1$NES)
d1$NAME<-substr(d1$NAME,10,nchar(d1$NAME))
ggplot(d1, aes(x = NES, y = fct_reorder(NAME, NES))) + 
  geom_point(aes(size = SIZE, color = NOM.p.val)) +
  theme_bw(base_size = 15) +
  scale_colour_gradient(limits=c(0, 0.75), low="red") +
  ylab(NULL) + 
  scale_size(name = "Count", limits = c(10,200)) +
  ggtitle("0 vs 4")

ggsave("/Users/tsears/code/NereaZappaProj/0v4_bubble.pdf")

#0V12
d2 = read.csv("/Users/tsears/Downloads/Hallmark_0v12/excel_files/gsea_report_for_12_1588020441104.xls",sep="\t",header=T)
d2$NES<-abs(d2$NES)
d2$NES<-abs(d2$NES)
d2$NAME<-substr(d2$NAME,10,nchar(d2$NAME))

ggplot(d2, aes(x = NES, y = fct_reorder(NAME, NES))) + 
  geom_point(aes(size = SIZE, color = NOM.p.val)) +
  theme_bw(base_size = 15) + 
  scale_colour_gradient(limits=c(0, 0.4), low="red") +
  scale_size(name = "Count", limits = c(10,200)) +
  ylab(NULL) +
  ggtitle("0 vs 12")

ggsave("0v12_bubble.pdf")

#0V18
d3 = read.csv("/Users/tsears/Downloads/Hallmark_0v18/excel_files/gsea_report_for_18_1588020680952.xls",sep="\t",header=T)
d3$NES<-abs(d3$NES)
d3$NAME<-substr(d3$NAME,10,nchar(d3$NAME))

ggplot(d3, aes(x = NES, y = fct_reorder(NAME, NES))) + 
  geom_point(aes(size = SIZE, color = NOM.p.val)) +
  theme_bw(base_size = 15) +
  scale_colour_gradient(limits=c(0, 0.4), low="red") +
  scale_size(name = "Count", limits = c(10,200)) +
  ylab(NULL) +
  ggtitle("0 vs 18")

ggsave("0v18_bubble.pdf")

#### Plotting from external papers ####
#GSEA was run on public data for comparison with our data

d4 = read.csv("gsea0v0.5_polyic.csv")

ggplot(d4, aes(x = NES, y = fct_reorder(NAME, NES))) + 
  geom_point(aes(size = SIZE, color = pval)) +
  theme_bw(base_size = 15) +
  scale_colour_gradient(limits=c(0, 1), low="red") +
  scale_size(name = "Count", limits = c(25,200)) +
  ylab(NULL) +
  ggtitle("GO Pathway Enrichment Hour 0 vs 0.5 (rath2019 PolyIC)")

ggsave("0v0.5_bubble.png")

d5 = read.csv("gsea0v1_polyic.csv")

ggplot(d5, aes(x = NES, y = fct_reorder(NAME, NES))) + 
  geom_point(aes(size = SIZE, color = pval)) +
  theme_bw(base_size = 15) +
  scale_colour_gradient(limits=c(0, 1), low="red") +
  scale_size(name = "Count", limits = c(25,200)) +
  ylab(NULL) +
  ggtitle("GO Pathway Enrichment Hour 0 vs 1 (rath2019 PolyIC)")

ggsave("0v1_bubble.png")

d6 = read.csv("gsea0v2_polyic.csv")

ggplot(d6, aes(x = NES, y = fct_reorder(NAME, NES))) + 
  geom_point(aes(size = SIZE, color = pval)) +
  theme_bw(base_size = 15) +
  scale_colour_gradient(limits=c(0, 1), low="red") +
  scale_size(name = "Count", limits = c(25,200)) +
  ylab(NULL) +
  ggtitle("GO Pathway Enrichment Hour 0 vs 2 (rath2019 PolyIC)")

ggsave("0v2_bubble.png")

d7 = read.csv("gsea0v4_polyic.csv")

ggplot(d7, aes(x = NES, y = fct_reorder(NAME, NES))) + 
  geom_point(aes(size = SIZE, color = pval)) +
  theme_bw(base_size = 15) +
  scale_colour_gradient(limits=c(0, 0.5), low="red") +
  scale_size(name = "Count", limits = c(25,200)) +
  ylab(NULL) +
  ggtitle("GO Pathway Enrichment Hour 0 vs 4 (rath2019 PolyIC)")

ggsave("0v4_bubble.png")

logmatrix <- data.frame("0"<-(matrix(0, length(res4$log2FoldChange), 1)), "04"<-res4$log2FoldChange, "12"<-res12$log2FoldChange, "18"<-res18$log2FoldChange)
rownames(logmatrix)<-rownames(res4)
goi2 <-which(rownames(logmatrix) %in% c("PARP16", "ERN1", "XBP1U", "XBP1", "BLOC1S1"))
logmatrix<-logmatrix[goi2,]
colnames(logmatrix)<-c("0 Hours","04 Hours","12 Hours","18 Hours")

logmatrix_1<-read.csv("logmatrix.csv")

ggplot(logmatrix_1, aes(x=Time, y=Log2FoldChange, group=Gene, color=Gene)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE, name="") +
  ggtitle("Selected Genes Log2FoldChange")

countmatrix <- data.frame("0"<-(matrix(0, length(res4$baseMean), 1)), "04"<-res4$baseMean, "12"<-res12$baseMean, "18"<-res18$baseMean)
rownames(countmatrix)<-rownames(res4)
goi3 <-which(rownames(countmatrix) %in% c("PARP16", "ERN1", "XBP1U", "XBP1", "BLOC1S1"))
countmatrix<-countmatrix[goi3,]
colnames(countmatrix)<-c("0 Hours","04 Hours","12 Hours","18 Hours")

write.csv(countmatrix,"countmatrix.csv")  
countmatrix_1 <- read.csv("countmatrix.csv")

ggplot(countmatrix_1, aes(x=Time, y=Raw_Count, group=ATF5, color=ATF5)) +
  geom_line() +
  scale_color_viridis(discrete = TRUE, name="") +
  ggtitle("Selected Genes Raw Count")

## heatmaps for top 50 in each GSEA group for each condition. This must be tailored to respective GSEA output
library(stringr)
file.names <- dir("/Users/tsears/code/NereaZappaProj/GSEA/Hallmark_0v18/excel_files", pattern ="HALLMARK")
file.names2<-file.names[c(3,15,18,21,22,25,26)]
for (i in file.names2) {
 
    gois<-read.table(paste("/Users/tsears/code/NereaZappaProj/GSEA/Hallmark_0v18/excel_files/",i,sep="")
                     ,sep="\t",header=T,fill=T,row.names = NULL)
    genes<-gois$GENE.SYMBOL
    
    #TC heatmap
    betas <- coef(dds_genelvl)
    colnames(betas)
    mat <- betas[intersect(rownames(betas),genes),c(4:6)]
    absmat<-abs(mat)
    mat<-mat[(rowMax(absmat))>0.9,] #Filter by magnitude of betas
    thr <- 3
    mat[mat < -thr] <- -thr
    mat[mat > thr] <- thr
    mat<-mat[order(mat[,3],decreasing = T),]
    p<-pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
           cluster_col=FALSE, cluster_rows = F,
            labels_col = c("4 vs 0","12 vs 0","18 vs 0"), fontsize_row = round(500/length(rownames(mat)))
           ,filename = paste("/Users/tsears/code/NereaZappaProj/Plots/full_pathway_all_genes/Test/All_Genes_",i,".pdf",sep=""),main = (str_sub(i,10,-5)))
    
  }

## top 50 genes in waterfall plot 

top50<-res18[order(abs(res18$log2FoldChange),decreasing = T),][1:50,]
bot50<-res18[order((res18$log2FoldChange),decreasing = F),][1:50,]
comb50<-rbind(top50,bot50)
top50<-top50[1:50,]
top50<-as.data.frame(top50)
top50<-top50[order(top50$log2FoldChange,decreasing = T),]
top50<-res18[intersect(rownames(res18),gois),] #For custom Waterfall Plot
df<-data.frame(names=c(rownames(top50)),log2FoldChange=c(top50$log2FoldChange))
ggplot(data=df, aes(y=reorder(names,log2FoldChange), x=log2FoldChange)) + ylab("Gene ID") +
  theme(axis.text.y = element_text(size = 5), axis.ticks.y = element_blank()) + scale_x_continuous(breaks=seq(-25,10,2.5)) +
  geom_bar(stat="identity", width=0.7, position = position_dodge(width=0.4),color="white", fill="steelblue") + 
  ggtitle("Top and Bottom 50 Genes Hour 18")

## top 300 genes in waterfall plot 

top50<-res_g[order(abs(res_g$log2FoldChange),decreasing = T),]
top50<-top50[1:300,]
top50<-as.data.frame(top50)
top50<-top50[order(top50$log2FoldChange,decreasing = T),]
df<-data.frame(names=c(rownames(top50)),log2FoldChange=c(top50$log2FoldChange))
ggplot(data=df, aes(y=reorder(names,log2FoldChange), x=log2FoldChange)) + ylab("Gene ID") +
  theme(axis.text.y = element_text(size = 1.3), axis.ticks.y = element_blank()) + scale_x_continuous(breaks=seq(-10,18,5)) +
  geom_bar(stat="identity", width=0.5, position = position_dodge(width=0.1), fill="steelblue") + 
  ggtitle("Top 300 Genes by Abs LFC")

## Check to see if a list of genes overlaps with the secretome

#import secretome data
setwd("/Users/tsears/Downloads/")
uniprot_signal<-read.csv("Uniprot_signal_peptide_Hsapiens.csv")
uniprot_transmembrane<-read.csv("Uniprot_transmembrane_Hsapiens.csv")
genes<-rbind(uniprot_signal,uniprot_transmembrane)
deduped.data <- unique(as.data.frame(genes)[,2] )
deduped.data<-as.data.frame(deduped.data)

#Repeat GO loop from earlier but calculate secretome percentages. This will output a histogram of secretome overlap 

file.names <- dir("/Users/tsears/code/NereaZappaProj/GSEA/Hallmark_0v18/excel_files", pattern ="HALLMARK")
file.names2<-file.names[c(15,18,21,25,26)]
for (i in file.names2) {
  holder<-c()
  counter<-1
  for (j in c("0v4","0v12","0v18")) {
    gois<-read.table(paste("/Users/tsears/code/NereaZappaProj/GSEA/Hallmark_", j, "/excel_files/",i,sep="")
                     ,sep="\t",header=T,fill=T,row.names = NULL)
    gois$RANK.METRIC.SCORE<-as.numeric(gois$RANK.METRIC.SCORE)
    gois1<-gois[order(abs(gois$RANK.METRIC.SCORE),decreasing = T),]
    gois1<-gois1[1:70,]
    holder[counter]<-(length((which(gois1$GENE.SYMBOL %in% deduped.data$deduped.data)))/length(gois1$GENE.SYMBOL))
    names(holder)[counter]<-j
    counter<-counter+1
  }
  holderdf<-as.data.frame(holder)
  holderdf<-cbind(holderdf,factor(rownames(holderdf),levels = c("0v4","0v12","0v18")))
  colnames(holderdf)<-c("Percentage","Timepoint")
  ggplot(holderdf, aes(x=Timepoint, y=Percentage)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Secretome Proportion Change Top 70"),i)
  ggsave(paste("SecPropChangeTop70",i,".pdf"))
}

## Plot multi-line graph of time couse LFC data

#Read in GOI's, split into multiple groups for clarity. Depends on respective GOIs
genes<-scan(text="TXNIP ATF4 XBP1 PPP1R15A ATF3 TRIB3 DDIT3 ATF5 CEBPA CEBPB IBTK SLAC35A4 CAT1 EPRS", what = "")
genes<-scan(text="TNFRSF1A TNFRSF1B FAS TNFRSF10A TNFRSF10B TNFRSF10C TNFRSF10D TNFSF10 FASLG TNF",what="")
genes<-scan(text="SQSTM1  BECLIN1 ULK1 MAP1LC3A MAP1LC3B NBR1 AMBRA1 VPS34 TFEB",what="")
genes<-scan(text="ATG1 ATG2A ATG2B ATG3 ATG4A ATG4B ATG4C ATG4D ATG5 ATG6 ATG7 ATG8a 
            ATG8b ATG9A ATG10 ATG11 ATG12 ATG18a ATG18b ATG13 ATG14 ATG15 ATG16L1 ATG16L2 ATG17 ATG101 TFEB",what="")

#Prepare dataframe for plotting
library(reshape2)
library(ggrepel)
library(viridis)
library(RColorBrewer)

betas <- coef(dds_genelvl)
res0.05 <- subset(res_g, (res_g$padj <= (0.05))) 
betas <- betas[rownames(res0.05),]
colnames(betas)
mat <- betas[intersect(rownames(betas),genes),c(4:6)]
absmat<-abs(mat)
mat<-mat[(rowMax(absmat))>0.585,]
zero_fill<-rep(0, each = nrow(mat))
plot_mat<-cbind(zero_fill,mat)
colnames(plot_mat)<-c(0,4,12,18)

Plot_df<-melt(plot_mat,id=c(rownames(plot_mat)))

Plot_df %>%
  mutate(label = if_else(Var2 == max(Var2), as.character(Var1), NA_character_)) %>%
  ggplot(aes(x = Var2, y = value, group = Var1, colour = Var1)) + 
  geom_line(size=1) + 
  geom_label_repel(aes(label = label),
                   nudge_x = 1,
                   na.rm = TRUE,
                   segment.colour = "grey64",
                   segment.size = 0.2,
                   point.padding = NA,
                   xlim=c(18,NA),) + 
  xlab("Hour") + ylab("Log2FoldChange") + ggtitle("Genes by LFC") +
  scale_x_continuous(breaks = c(4,12,18)) +
  coord_cartesian(clip = 'off',xlim = c(NA,21)) +
  theme(legend.position = 'none') + 
  scale_color_brewer(type="qual",palette = "Paired")
   
ggsave("/Users/tsears/code/NereaZappaProj/Plots/LineGraphs/Group3.pdf")




       