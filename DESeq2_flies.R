library(DESeq2)
setwd("~/Documents/UCDavis/ECE221/ECE221_final_project")
flymatrix<-read.csv("flies_space_counts.csv")
# SRR3390478.quant = G1R1
# SRR3390479.quant = G1R2
# SRR3390480.quant = G1R3
# SRR3390481.quant = G1R4
# SRR3390482.quant = G3R1
# SRR3390483.quant = G3R2
# SRR3390484.quant = G3R3
# SRR3390485.quant = G3R4
colnames(flymatrix)<-c("Name","G1R1","G1R2","G1R3","G1R4","G3R1","G3R2","G3R3","G3R4")
head(flymatrix)
group <- factor(c(1,1,1,1,2,2,2,2))
design <- model.matrix(~group)
