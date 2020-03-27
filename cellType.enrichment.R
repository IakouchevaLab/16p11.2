rm(list = ls());gc()
options(stringsAsFactors = F)
library(WGCNA)
library(EWCE)

#load("F:/projects/ipsc/finallme/wgcna/sv1/ds0.mod100/ipsc.lme.network.recut.RData")
load("D:/projects/organoid/rna1month/finallme2/wgcna/ds2.mod50/organoid.lme.network.recut.RData")
#load("D:/projects/organoid/rna3month/sv1.6/ds2.mod70.RData")

prefix="1month_"

geneTree = networks$tree
datExpr=networks$datExpr
merged = networks$merged
modules = merged$colors
MEs = networks$MEs
kMEtable = networks$kMEtable

annot=read.delim("D:/references/gencode.v19.gene.name.txt",header = T,sep = "\t")
annot=annot[match(rownames(datExpr),annot$geneid),]


#Load cell-type expression signatures for enrichment
reps=10000
level=1
load("CellTypeData_DamonNeuralFetalOnly.rda")
ctd_fetal = ctd
rm(ctd)
bg=annot$genename
background_fetal = intersect(bg,rownames(ctd_fetal[[level]]$specificity))

# loop through each module except gray
fetal_all=data.frame()
for(i in 2:ncol(MEs$eigengenes)) {
  
  me = MEs$eigengenes[,i]
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  moduleGenes = annot$genename[modules==moduleNumber]

  input_fetal = intersect(moduleGenes,rownames(ctd_fetal[[level]]$specificity))
  
  
  #Cell-Type enrichment
  fetal_cells=bootstrap.enrichment.test(sct_data=ctd_fetal,hits=input_fetal,bg=background_fetal,reps=reps,annotLevel=level,genelistSpecies="human",sctSpecies="human")
  
  fetal_results=fetal_cells$results
  fetal_results$fdr=p.adjust(fetal_results$p)
  fetal_results$module=paste(prefix,moduleNumber,moduleColor,sep = "")
  
  fetal_all=rbind(fetal_all,fetal_results)
  
}

save(fetal_all,file = "ewce.1moOrg.RData")



# plot fetal

flag=unique(fetal_all$module[fetal_all$fdr < 0.05])
dat=fetal_all[fetal_all$module %in% flag,]
#dat$realp=2*pnorm(-abs(dat$sd_from_mean))
dat$sd_from_mean[dat$sd_from_mean < 0]=0
dat$fdr[dat$fdr == 0]=min(dat$fdr[dat$fdr>0]) * 0.1
dat$sd_from_mean=round(dat$sd_from_mean,1)
dat$sd_from_mean[dat$fdr > 0.05]=""
dat$moduleNumber=as.numeric(str_sub(str_extract(dat$module,"_[0-9]*"),2,-1))
dat$module=reorder(dat$module,dat$moduleNumber)
ggplot(dat,aes(CellType,module,label=sd_from_mean))+
  geom_tile(aes(fill=-log10(fdr)))+
  geom_text(size=3, color="black")+
  scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = 1.3)
