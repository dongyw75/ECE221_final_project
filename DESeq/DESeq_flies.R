#http://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.R

library(DESeq)
library("lattice")
library(biomaRt)
library(tximport)
library(readr)
library("RColorBrewer")
setwd("~/Documents/UCDavis/ECE221/ECE221_final_project")
flymatrix<-read.csv("flies_space_counts_salmon.csv")
rownames(flymatrix)<-flymatrix$Name
flymatrix<-flymatrix[,c(2:9)]
colnames(flymatrix)
head(flymatrix)
# SRR3390478.quant = G1R1
# SRR3390479.quant = G1R2
# SRR3390480.quant = G1R3
# SRR3390481.quant = G1R4
# SRR3390482.quant = G3R1
# SRR3390483.quant = G3R2
# SRR3390484.quant = G3R3 = G1R1
# SRR3390485.quant = G3R4

# import salmon
dir<-"/Users/cohenl06/Documents/UCDavis/ECE221/ECE221_final_project"
list.files("~/Documents/UCDavis/ECE221/ECE221_final_project")
files <- file.path(dir, "salmon",c("SRR3390478.quant","SRR3390479.quant","SRR3390480.quant","SRR3390481.quant","SRR3390482.quant","SRR3390483.quant","SRR3390484.quant","SRR3390485.quant"), "quant.sf")
#names(files) <- c("G1R1","G1R2","G1R3","G1R4","G3R1","G3R2","G3R3","G3R4")
names(files) <- c("G3R3","G1R2","G1R3","G1R4","G3R1","G3R2","G1R1","G3R4")
files
# get annotations
ensembl=useMart("ensembl")
ensembl_id<-rownames(flymatrix)
ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
query<-getBM(attributes=c('flybase_transcript_id','ensembl_gene_id','external_gene_name','description'), filters = 'flybase_transcript_id', values = ensembl_id, mart=ensembl)
tx2gene<-query[,c(1,2)]
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
head(txi.salmon$counts)
dim(txi.salmon$counts)
# convert to integers
salmon_counts_int<-round(txi.salmon$counts,digits=0)
head(salmon_counts_int)
#write.csv(txi.salmon$counts,"salmon_counts_est_genes.csv")
colnames(salmon_counts_int)
condition = factor(c("G3","G1","G1","G1","G3","G3","G1","G3"))
ExpDesign <- data.frame(row.names=colnames(salmon_counts_int), condition = c("G3","G1","G1","G1","G3","G3","G1","G3"))
ExpDesign
cds = newCountDataSet( salmon_counts_int, condition )
cds = estimateSizeFactors( cds )
notAllZero = (rowSums(counts(cds))>0)
cdsFullBlind = estimateDispersions( cds, method = "blind" )
vsdFull = varianceStabilizingTransformation( cdsFullBlind )
plotPCA(vsdFull, intgroup="condition")

sizeFactors( cds )

head( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds )
str( fitInfo(cds) )
plotDispEsts( cds )
all(table(conditions(cds))==2)
head( fData(cds) )
res = nbinomTest( cds, "G1", "G3" )
head(res)
stopifnot(identical(colnames(res), c("id", "baseMean", "baseMeanA", "baseMeanB", "foldChange",
                                     "log2FoldChange", "pval", "padj")))
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")


# get counts
counts_table = counts( cds, normalized=TRUE )
filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1, ]
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
dim(filtered_norm_counts)
# get annotations
ensembl=useMart("ensembl")
ensembl_id<-rownames(flymatrix)
ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
query<-getBM(attributes=c('flybase_transcript_id','ensembl_gene_id','external_gene_name','description'), filters = 'flybase_transcript_id', values = ensembl_id, mart=ensembl)
#
colnames(query)
col.names<-c("flybase_transcript_id","GeneID","external_gene_id","description")
colnames(query)<-col.names
query_gene<-query[,c(2,3,4)]
dim(query_gene)
merge_biomart_res_counts <- merge(filtered_norm_counts,query_gene,by=unique("GeneID"))
merge_biomart_res_counts<-merge_biomart_res_counts[!duplicated(merge_biomart_res_counts[,c(1)]),]
dim(merge_biomart_res_counts)
head(merge_biomart_res_counts)
head(res)
colnames(res)
colnames(res)<-c("GeneID","baseMean","baseMeanA","baseMeanB","foldChange","log2FoldChange","pval","padj")
merged<-merge(merge_biomart_res_counts,res,by="GeneID")
merged<-merged[order(merged[,18]), ]
write.csv(merged, file="DESeq/hypergravity_space_flies_DESeq_all.csv" )
resSig = merged[ merged$padj < 0.05, ]
#write.csv(resSig,file="DESeq/hypergravity_space_flies_DESeq_padj0.05.csv")
#
plot(log2(merged$baseMean), merged$log2FoldChange, col=ifelse(merged$padj < 0.05, "red","gray67"),main="(DESeq2) G3 vs. G1 (padj<0.05, log2FC = ±1)",xlim=c(1,20),pch=20,cex=1,ylim=c(-12,12))
abline(h=c(-1,1), col="blue")
# heatmap
up_down<-resSig
dim(up_down)
up_down_1FC<-subset(up_down,up_down$log2FoldChange>1 | up_down$log2FoldChange< -1)
write.csv(up_down_1FC,file="DESeq/hypergravity_space_flies_DESeq_padj0.05.csv")
dim(up_down_1FC)
d<-up_down_1FC
d<-na.omit(d)
dim(d)
d<-up_down_1FC[,c(2:9)]
d<-as.matrix(d)
d<-as.data.frame(d)
d<-as.matrix(d)
rownames(d) <- up_down_1FC[,10]
d<-type.convert(d)
head(d)
colnames(d)
d<-d[,c(7,2,3,4,5,6,1,8)]
colnames(d)

hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="G3 vs. G1 (padj<0.05, log2FC = ±1)", 
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)

