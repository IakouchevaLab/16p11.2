options(stringsAsFactors = TRUE)
library(limma) 
library(edgeR)
library(DESeq2)
library(calibrate)

# DGE using limma-voom + dupCorr + lmFit
genes = DGEList(counts)
genes = calcNormFactors( genes )
design = model.matrix( ~ Genotype+sv1+sv2+sv3+sv4+sv5, info)
colnames(design)[1]="Intercept"
vobj_tmp = voom( genes, design, plot=FALSE)

#run duplicate correlation by Individual
dupcor <- duplicateCorrelation(vobj_tmp,design,block=info$Individual)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights

vobj = voom( genes, design, plot=FALSE, block=info$Individual, correlation=dupcor$consensus)
fitDupCor <- lmFit(vobj, design, block=info$Individual, correlation=dupcor$consensus.correlation)

# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes( fitDupCor )

# DGE tables
summary(decideTests(fitDupCor))
delTable=topTable(fitDupCor,coef = "GenotypeDEL", number = Inf, sort.by = "p")
dupTable=topTable(fitDupCor,coef = "GenotypeDUP", number = Inf, sort.by = "p")
contrast=makeContrasts(DUPDEL=GenotypeDUP-GenotypeDEL,levels = design)
contraFit=contrasts.fit(fitDupCor,contrast)
contraFit=eBayes(contraFit)
summary(decideTests(contraFit))
dupdelTable=topTable(contraFit,coef = "DUPDEL", number = Inf, sort.by = "p")

# annotation 
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") 
attributes = c("ensembl_gene_id","gene_biotype","description","external_gene_id")
annotation = getBM(attributes, mart=ensembl)

# 1. DEL vs CTL
delTable$GeneSymbol <- annotation$external_gene_id[match(row.names(delTable),annotation$ensembl_gene_id)]
delTable$Gene_type <- annotation$gene_biotype[match(row.names(delTable),annotation$ensembl_gene_id)]
delTable$Description <- annotation$description[match(row.names(delTable),annotation$ensembl_gene_id)]

# 2. DUP vs CTL
dupTable$GeneSymbol <- annotation$external_gene_id[match(row.names(dupTable),annotation$ensembl_gene_id)]
dupTable$Gene_type <- annotation$gene_biotype[match(row.names(dupTable),annotation$ensembl_gene_id)]
dupTable$Description <- annotation$description[match(row.names(dupTable),annotation$ensembl_gene_id)]

#3. DUP vs DEL
dupdelTable$GeneSymbol <- annotation$external_gene_id[match(row.names(dupdelTable),annotation$ensembl_gene_id)]
dupdelTable$Gene_type <- annotation$gene_biotype[match(row.names(dupdelTable),annotation$ensembl_gene_id)]
dupdelTable$Description <- annotation$description[match(row.names(dupdelTable),annotation$ensembl_gene_id)]

write.table(delTable,file="1m_limma_voom_dupCorr_lmFit_2rounds_voom_DELvsCTL_SV1toSV5_annotated.txt",sep="\t")
write.table(dupTable,file="1m_limma_voom_dupCorr_lmFit_2rounds_voom_DUPvsCTL_SV1toSV5_annotated.txt",sep="\t")
write.table(dupdelTable,file="1m_limma_voom_dupCorr_lmFit_2rounds_voom_DUPvsDEL_SV1toSV5_annotated.txt",sep="\t")

# Volcano Plot
#list of genes within the 16p11.2 CNV
genes_16p <- read.table(file="16p_genes.txt",header=FALSE,sep='\n')

with(delTable, plot(logFC, adj.P.Val, pch=20, main="Volcano plot"))

# A. Volcano Plot DEL vs CTL
# Add colored points:
# blue if padj<0.1 and logFC < -0.58 (i.e. FC=1.5) 
# orange if padj<0.1 and logFC > 0.58 (i.e. FC=1.5)
with(subset(delTable, padj< 0.1 & logFC < -0.58), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
with(subset(delTable, padj< 0.1 & logFC > 0.58), points(logFC, -log10(adj.P.Val), pch=20, col="orange"))
with(subset(delTable, genes_16p), points(logFC, -log10(adj.P.Val), pch=20, col="deeppink"))

# Label genes
with(subset(delTable, padj< 0.1 & logFC < -0.58), textxy(logFC, -log10(adj.P.Val), labs=GeneSymbol, cex=.8))
with(subset(delTable, padj< 0.1 & logFC > 0.58), textxy(logFC, -log10(adj.P.Val), labs=GeneSymbol, cex=.8))

# B. Volcano Plot DUP vs CTL
# Add colored points:
# blue if padj<0.1 and logFC < -0.58 (i.e. FC=1.5) 
# orange if padj<0.1 and logFC > 0.58 (i.e. FC=1.5)
with(subset(dupTable, padj< 0.1 & logFC < -0.58), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
with(subset(dupTable, padj< 0.1 & logFC > 0.58), points(logFC, -log10(adj.P.Val), pch=20, col="orange"))
with(subset(dupTable, genes_16p), points(logFC, -log10(adj.P.Val), pch=20, col="deeppink"))

# Label genes
with(subset(dupTable, padj< 0.1 & logFC < -0.58), textxy(logFC, -log10(adj.P.Val), labs=GeneSymbol, cex=.8))
with(subset(dupTable, padj< 0.1 & logFC > 0.58), textxy(logFC, -log10(adj.P.Val), labs=GeneSymbol, cex=.8))

# C. Volcano Plot DUP vs DEL
# Add colored points:
# blue if padj<0.1 and logFC < -0.58 (i.e. FC=1.5) 
# orange if padj<0.1 and logFC > 0.58 (i.e. FC=1.5)
with(subset(dupdelTable, padj< 0.1 & logFC < -0.58), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
with(subset(dupdelTable, padj< 0.1 & logFC > 0.58), points(logFC, -log10(adj.P.Val), pch=20, col="orange"))
with(subset(dupdelTable, genes_16p), points(logFC, -log10(adj.P.Val), pch=20, col="deeppink"))

# Label genes
with(subset(dupdelTable, padj< 0.1 & logFC < -0.58), textxy(logFC, -log10(adj.P.Val), labs=GeneSymbol, cex=.8))
with(subset(dupdelTable, padj< 0.1 & logFC > 0.58), textxy(logFC, -log10(adj.P.Val), labs=GeneSymbol, cex=.8))

