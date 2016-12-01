source('extra_scripts/overLapper_original.R')
deseq1<-read.csv("DESeq/hypergravity_space_flies_DESeq_padj0.05.csv")
deseq2<-read.csv("DESeq2/hypergravity_space_flies_DESeq2_padj0.05.csv")
edger<-read.csv("edgeR/EdgeR_sig_gene.csv")
dim(deseq2)
dim(deseq1)
dim(edger)
colnames(deseq2)
colnames(deseq1)
# venn diagram
# DESeq1
m<-deseq1$GeneID
length(m)
# DESeq2
n<-deseq2$GeneID
length(n)
# edgeR

o<-edger$GeneID
length(o)
p<-EBDERes$DEfound
length(p)
setlist <- list(DESeq1=as.vector(m),DESeq2=as.vector(n),edgeR=as.vector(o),EBSeq=as.vector(p))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)

