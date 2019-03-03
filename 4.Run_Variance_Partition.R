library(variancePartition)
library(limma)
library(edgeR)

# Variance Partition
gExpr <- DGEList(counts=counts)
gExpr <- calcNormFactors(gExpr)
design <- model.matrix( ~ Genotype + Run +Lab, info)
vobjGenes <- voom(gExpr, design )
geneExpres <-vobjGenes$E

form <- ~ (1|Genotype) + (1|Individual) + (1|Run) + (1|Lab) 
varPart <- fitExtractVarPartModel( vobjGenes, form, info )
vp <- sortCols(varPart )
plotVarPart(vp)

write.table(vp, file="1month_Org_VariancePartition_variability.txt", sep="\t")

# Top most variable genes
sortedLab <- (varPart[order(varPart$Lab, decreasing=TRUE),])
i <- sortedLab$Lab > 0.5
topLab <- sortedLab[i, ]
write.table(topLab, file="Genes_mostVariable_byLab_0.5.txt", sep="\t")

# SECOND FILTERING STEP - REMOVE GENES WITH HIGH VARIANCE BETWEEN LABS
hvLab <- rownames(topLab)
hvLab <-as.vector(hvLab)

tpm <- tpm[! rownames(tpm) %in% hvLab, ]
counts <- counts[!rownames(counts) %in% hvLab, ]


