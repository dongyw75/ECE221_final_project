setwd("C:/Users/JiyingLi/Desktop/project221limma")
library(limma)
x<-read.csv("space_fly_gene_counts.csv")
rownames(x)=x[,1]
x<-subset(x, select=-X)
x<-as.matrix(x)
pca = prcomp(t(x))
names = colnames(x)
fac= factor(c("G3","G1","G1","G1","G3","G3","G1","G3"))
colours = c("red","blue")

install.packages( pkgs= "gplots" )
library(gplots)
library(lattice)
xyplot(
PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
panel=function(x, y, ...) {
panel.xyplot(x, y, ...);
ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=1)
},
aspect = "fill", col=colours
#main = draw.key(key = list(rect = list(col = list(col=colours), text = list(levels(fac)), rep = FALSE)))
)

library("edgeR")
dge <- DGEList(counts=x)
dge <- calcNormFactors(dge)


design <- model.matrix(~ 0+factor(c(2,1,1,1,2,2,1,2)))
colnames(design) <- c("G1 ",  "G3")


v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef=ncol(design))
write.csv( topTable(fit, coef=ncol(design), number=15000), "limma_siggenes.csv")


