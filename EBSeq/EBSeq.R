#http://www.bioconductor.org/packages/devel/bioc/vignettes/EBSeq/inst/doc/EBSeq_Vignette.R
library(EBSeq)
library(biomaRt)
library(tximport)
library(readr)
library("RColorBrewer")
setwd("~/Documents/UCDavis/ECE221/ECE221_final_project")
source("extra_scripts/intensityFilter.R")
dir<-"/Users/cohenl06/Documents/UCDavis/ECE221/ECE221_final_project"
flymatrix<-read.csv("flies_space_counts_salmon.csv")
# get annotations
ensembl=useMart("ensembl")
ensembl_id<-rownames(flymatrix)
ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
query<-getBM(attributes=c('flybase_transcript_id','ensembl_gene_id','external_gene_name','description'), filters = 'flybase_transcript_id', values = ensembl_id, mart=ensembl)
tx2gene<-query[,c(1,2)]
list.files("~/Documents/UCDavis/ECE221/ECE221_final_project")
files <- file.path(dir, "salmon",c("SRR3390478.quant","SRR3390479.quant","SRR3390480.quant","SRR3390481.quant","SRR3390482.quant","SRR3390483.quant","SRR3390484.quant","SRR3390485.quant"), "quant.sf")
#names(files) <- c("G1R1","G1R2","G1R3","G1R4","G3R1","G3R2","G3R3","G3R4")
names(files) <- c("G3R3","G1R2","G1R3","G1R4","G3R1","G3R2","G1R1","G3R4")
files
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
head(txi.salmon$counts)
dim(txi.salmon$counts)
# EBSeq
Sizes=MedianNorm(txi.salmon$counts)
#ExpDesign <- data.frame(row.names=colnames(txi.salmon$counts), condition = c("G3","G1","G1","G1","G3","G3","G1","G3"))
#ExpDesign
x<-as.factor(c("G3","G1","G1","G1","G3","G3","G1","G3"))
EBOut=EBTest(Data=txi.salmon$counts, 
             Conditions=x,sizeFactors=Sizes, maxround=5)
EBDERes=GetDEResults(EBOut, FDR=0.05)
str(EBDERes$DEfound)
head(EBDERes$PPMat)
str(EBDERes$Status)
dim(EBDERes$PPMat)
head(EBDERes$PPMat)
GeneFC=PostFC(EBOut)
PlotPostVsRawFC(EBOut,GeneFC)
write.csv(EBDERes$DEfound,"EBSeq/DEfound.csv")
table<-as.data.frame(EBDERes$DEfound)
results<-as.data.frame(EBDERes$PPMat)
head(results)
colnames(table)<-c("GeneID")
head(table)

## Simulate some data
set.seed(1234)

## Expression level
logCPM<- sort(rnorm(n= 10000, mean= 5, sd= 1))

## Log fold change
set.seed(12345)
logFC<- sapply(length(logCPM):1, function(i) rnorm(n= 1, mean= 0, sd= i^(1/1.5)))

smoothScatter(x= logCPM, y= logFC, nrpoints= 1000)

## z-score
z<- localZ(logCPM, logFC, nbins= 20)

## Highlight genes on the periphery of the cloud
points(logCPM, logFC, col= ifelse(abs(z) > 1.5, 'red', NA), pch= '.', cex= 0.5)