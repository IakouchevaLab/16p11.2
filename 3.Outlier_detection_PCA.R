library(limma) 
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)
library(WGCNA)
library(DESeq2)
library(ggplot2)
library(magrittr)
library(Rmisc) 
library(ggrepel)
library(dplyr)
library(WGCNA); 
library(SummarizedExperiment)


# load counts and tpm files
load("./RSEM_Quant.genes.counts.RData")
load("./RSEM_Quant.genes.tpm.RData")

# metadata file
info <- read.table(file="./metadata.txt",header=TRUE,sep='\t',row.names=1)

# Filter expressed genes
ctl <- rownames(info[info$Genotype== "CTL", ])
dup <- rownames(info[info$Genotype== "DUP", ])
del <- rownames(info[info$Genotype== "DEL", ])

ctl_tpm <- tpm [,colnames(tpm) %in% ctl, drop=F]
dup_tpm <- tpm [,colnames(tpm) %in% dup, drop=F]
del_tpm <- tpm [,colnames(tpm) %in% del, drop=F]

# CTL : filter genes with TPM 0.5 in 80% of CTL samples
pres = apply(ctl_tpm>0.5,1,sum) 
Expressed_ctl = (pres > 0.80*ncol(ctl_tpm))
ctl_tpm = ctl_tpm[Expressed_ctl,]
gene_names_ctl <- rownames(ctl_tpm)

# DUP : filter genes with TPM 0.5 in 80% of DUP samples
pres_dup = apply(dup_tpm>0.5,1,sum) 
Expressed_dup = (pres_dup > 0.80*ncol(dup_tpm))
dup_tpm = dup_tpm[Expressed_dup,]
gene_names_dup <- rownames(dup_tpm)

# DEL : filter genes with TPM 0.5 in 80% of DEL samples
pres_del = apply(del_tpm>0.5,1,sum) 
Expressed_del = (pres_del > 0.80*ncol(del_tpm))
del_tpm = del_tpm[Expressed_del,]
gene_names_del <- rownames(del_tpm)

# Expressed genes
ctl_dup <- union(gene_names_ctl,gene_names_dup)
Expressed_genes <- union (ctl_dup,gene_names_del )
write.table(Expressed_genes,file="Expressed_genes_TPM0.5_in80percSamples.txt")

# Expressed genes file
Expressed<- read.table(file="./Expressed_genes_TPM0.5_in80percSamples.txt",header=FALSE,sep="\n")
Expressed <- as.vector(t(Expressed))
tpm <- tpm[rownames(tpm) %in% Expressed, ]
counts <- counts[rownames(counts) %in% Expressed, ]

#outlier detection based on gene expression connectivity
all_outliers = character()
par(mfrow=c(2,4), mar=c(5,4,2,2))
info$ConnectivityZscore = NA

idx = match(colnames(tpm), rownames(info))
normadj <- (0.5+0.5*bicor(tpm))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
info$ConnectivityZscore[idx] = z.ku

#Plot outliers
g1 = ggplot(info, aes(x=info$Genotype, y=ConnectivityZscore, color=Genotype))+ theme_classic() + geom_point(size = 5) + geom_hline(yintercept = -2, lty=2) +  
  geom_hline(yintercept = 2, lty=2) +theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + xlab("") + ggtitle("Outlier Sample Removal") +
  ylab("Network Connectivity (Z score)")   + geom_text_repel(aes(label=Sample),size=6) +  theme(text = element_text(size=26))

# Remove outliers
dont.want <- c("1m_101C1_a")
tpm<- tpm[, ! colnames(tpm) %in% dont.want, drop = F]
counts<-counts[, ! colnames(counts) %in% dont.want, drop = F]
info<-info[ ! rownames(info) %in% dont.want, ]

#PCA
counts <- round(counts)
dds<-DESeqDataSetFromMatrix(countData = counts,
                            colData = info,
                            design = ~ 1)

rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup = c("Genotype"))
z <- plotPCA(rld, intgroup = c("Genotype"))
z + geom_text_repel(aes(label = name))
