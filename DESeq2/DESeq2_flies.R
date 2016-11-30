#detach("package:DESeq", unload=TRUE)
library(DESeq2)
library("lattice")
library(biomaRt)
library(tximport)
library(readr)
library(gplots)
library(RColorBrewer)
setwd("~/Documents/UCDavis/ECE221/ECE221_final_project")
source('~/Documents/UCDavis/ECE221/ECE221_final_project/extra_scripts/plotPCAWithSampleNames.R')
flymatrix<-read.csv("flies_space_counts_salmon.csv")
head(flymatrix)
# SRR3390478.quant = G1R1
# SRR3390479.quant = G1R2
# SRR3390480.quant = G1R3
# SRR3390481.quant = G1R4
# SRR3390482.quant = G3R1
# SRR3390483.quant = G3R2
# SRR3390484.quant = G3R3 = G1R1 (fixed)
# SRR3390485.quant = G3R4
# this is correct:
colnames(flymatrix)<-c("Name","G3R3","G1R2","G1R3","G1R4","G3R1","G3R2","G1R1","G3R4")
# this is wrong:
#colnames(flymatrix)<-c("Name","G1R1","G1R2","G1R3","G1R4","G3R1","G3R2","G3R3","G3R4")
head(flymatrix)
rownames(flymatrix)<-flymatrix$Name
flymatrix<-flymatrix[,c(2:9)]
colnames(flymatrix)
head(flymatrix)
# get annotations
ensembl=useMart("ensembl")
ensembl_id<-rownames(flymatrix)
ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
query<-getBM(attributes=c('flybase_transcript_id','ensembl_gene_id','external_gene_name','description'), filters = 'flybase_transcript_id', values = ensembl_id, mart=ensembl)
tx2gene<-query[,c(1,2)]
colnames(tx2gene)<-c("TXNAME","GENEID")
####
# import salmon
dir<-"/Users/cohenl06/Documents/UCDavis/ECE221/ECE221_final_project"
list.files("~/Documents/UCDavis/ECE221/ECE221_final_project")
files <- file.path(dir, "salmon",c("SRR3390478.quant","SRR3390479.quant","SRR3390480.quant","SRR3390481.quant","SRR3390482.quant","SRR3390483.quant","SRR3390484.quant","SRR3390485.quant"), "quant.sf")
# this is incorrect:
#names(files) <- c("G1R1","G1R2","G1R3","G1R4","G3R1","G3R2","G3R3","G3R4")
# this is correct:
names(files) <- c("G3R3","G1R2","G1R3","G1R4","G3R1","G3R2","G1R1","G3R4")
files
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
head(txi.salmon$counts)
dim(txi.salmon$counts)
# do not need to do this each time:
#write.csv(txi.salmon$counts,"salmon_counts_est_genes.csv")
# this is for the data from flymatrix (transcript-level counts)
#ExpDesign <- data.frame(row.names=colnames(flymatrix), condition = c("G1","G1","G1","G1","G3","G3","G3","G3"))
# this is for the estimated gene-level counts from tximport - use this:
ExpDesign <- data.frame(row.names=colnames(txi.salmon$counts), condition = c("G3","G1","G1","G1","G3","G3","G1","G3"))
dds <- DESeqDataSetFromTximport(txi.salmon, ExpDesign, ~condition)
#dds<-DESeqDataSetFromMatrix(countData=flymatrix, colData=ExpDesign,design=~condition)
dds$condition <- relevel(dds$condition, "G1")
dds<-DESeq(dds,betaPrior=FALSE)
norm_counts<-counts(dds,normalized=TRUE)
norm_counts_data<-as.data.frame(norm_counts)
dim(norm_counts_data)
filtered_norm_counts<-norm_counts_data[!rowSums(norm_counts_data==0)>=1, ]
dim(filtered_norm_counts)
head(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
dim(filtered_norm_counts)
colnames(query)
col.names<-c("flybase_transcript_id","GeneID","external_gene_id","description")
colnames(query)<-col.names
query_gene<-query[,c(2,3,4)]
dim(query_gene)
merge_biomart_res_counts <- merge(filtered_norm_counts,query_gene,by=unique("GeneID"))
counts<-merge_biomart_res_counts[!duplicated(merge_biomart_res_counts$GeneID), ]
head(counts)
dim(counts)
###
plotDispEsts(dds)
log_dds<-rlog(dds)
plotPCAWithSampleNames(log_dds, intgroup="condition", ntop=40000)
###
res<-results(dds,contrast=c("condition","G3","G1"))
res_ordered<-res[order(res$padj),]
GeneID<-rownames(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_ordered<-cbind(res_ordered,GeneID)
merge_biomart_res_counts <- merge(counts,res_ordered,by="GeneID")
merge_biomart_res_all<-subset(merge_biomart_res_counts,merge_biomart_res_counts$padj!="NA")
merge_biomart_res_all<-merge_biomart_res_all[order(merge_biomart_res_all$padj),]
dim(merge_biomart_res_all)
head(merge_biomart_res_all)
###
res_merged_cutoff<-subset(merge_biomart_res_all,merge_biomart_res_all$padj<0.05)
dim(res_merged_cutoff)
plot(log2(res$baseMean), res$log2FoldChange, col=ifelse(res$padj < 0.05, "red","gray67"),main="(DESeq2) G3 vs. G1 (padj<0.05, log2FC = ±1)",xlim=c(1,20),pch=20,cex=1,ylim=c(-12,12))
abline(h=c(-1,1), col="blue")
###
# make a table of differentially expressed genes
##
# don't have to do this each time

res_merged_cutoff<-subset(merge_biomart_res_all,merge_biomart_res_all$padj<0.05)
# don't have to do this each time
#write.csv(res_merged_cutoff,"hypergravity_space_flies_DESeq2_padj0.05.csv")
# heatmap
###
up_down<-res_merged_cutoff
dim(up_down)
up_down_1FC<-subset(res_merged_cutoff,res_merged_cutoff$log2FoldChange>1 | res_merged_cutoff
                    $log2FoldChange< -1)
write.csv(up_down_1FC,"DESeq2/hypergravity_space_flies_DESeq2_padj0.05.csv")
dim(up_down_1FC)
d<-up_down_1FC
d<-as.matrix(up_down_1FC[,c(2:9)])
rownames(d) <- up_down_1FC[,10]
d<-na.omit(d)
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
