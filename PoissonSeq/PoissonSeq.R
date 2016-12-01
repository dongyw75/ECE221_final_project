#install.packages("PoissonSeq")
library("PoissonSeq")
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
class(salmon_counts_int)
# PoissonSeq
y <- c(2, 1, 1, 1, 2, 2, 1, 2)
type <- "twoclass"
pair <- FALSE
gname=rownames(salmon_counts_int)
dat <- list(n=txi.salmon$counts, y=y, type=type, pair=pair, gname=gname)
para <- list(npermu=200, seed=100)
res <- PS.Main(dat=dat, para=para)
write.table(res, file="PoissonSeq_all.txt", quote=F, row.names=F, col.names=T, sep="\t")
head(res)
plot(res$nc, res$fdr, xlab="Number called", ylab="estimated FDR", type="l")
depth <- PS.Est.Depth(salmon_counts_int)
baseMean_df<-as.data.frame(rowMeans(dat$n))
baseMean_df<-cbind(baseMean_df,rownames(baseMean_df))
colnames(baseMean_df)<-c("rowMeans","GeneID")
res_df<-cbind(res,rownames(res))
colnames(res_df)<-c("nc","gname","tt","pval","fdr","log.fc","GeneID")
res_df<-merge(res_df,baseMean_df,by="GeneID")
filtered_norm_counts<-as.data.frame(dat$n[!rowSums(dat$n==0)>=1, ])
filtered_norm_counts<-cbind(filtered_norm_counts,rownames(filtered_norm_counts))
colnames(filtered_norm_counts)<-c("G3R3","G1R2","G1R3","G1R4","G3R1","G3R2","G1R1","G3R4","GeneID")
res_df<-merge(filtered_norm_counts,res_df,by="GeneID")
res_df<-na.omit(res_df)
res_df<-res_df[order(res_df[,14]), ]
write.csv(res_df,"PoissonSeq/PoissonSeq_all.csv")
plot(log2(res_df$rowMeans), res_df$log.fc, col=ifelse(res$fdr < 0.05, "red","gray67"),main="(PoissonSeq) G3 vs. G1 (fdr<0.05, log2FC = ±1)",xlim=c(1,20),pch=20,cex=1,ylim=c(-12,12))
abline(h=c(-1,1), col="blue")
resSig = res_df[ res_df$fdr < 0.05, ]
dim(resSig)
up_down_1FC<-subset(resSig,resSig$log.fc>1 | resSig$log.fc< -1)
dim(up_down_1FC)
write.csv(up_down_1FC,"PoissonSeq/PoissonSeq_sig.csv")
