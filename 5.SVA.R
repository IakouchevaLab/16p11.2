#SVA
library(sva)

dge.cpm=cpm(calcNormFactors(DGEList(counts = counts),method = 'TMM'))
mod1 = model.matrix(~Genotype, data=info) 
mod0 = mod1[,!grepl("Genotype", colnames(mod1))] 

sva = svaseq(dge.cpm,mod1,mod0,n.sv=5)
sv=sva$sv
colnames(sv)=paste("sv",seq(1,5),sep = "")
info=cbind(info,sv)
