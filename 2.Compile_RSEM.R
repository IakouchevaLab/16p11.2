#Patricia Moran 02/27/2019
#compile RSEM gene matrix

options(stringsAsFactors=F)

setwd("./pmlosada/16p11.2-RNASeq/RSEM_files")
files = dir(pattern="*.RSEM_Quant.genes.results",recursive=T)

inFile=read.delim(files[[1]],row.names=1)
genes = row.names(inFile)
counts=inFile$expected_count
tpm=inFile$TPM

for(i in 2:length(files)) {
  print(i)
  inFile=read.delim(files[[i]],row.names=1)
  inFile = inFile[match(genes,rownames(inFile)),]
  counts = cbind(counts, inFile$expected_count)
  tpm = cbind(tpm, inFile$TPM)
}


rownames(counts) = rownames(tpm) = substr(genes,1,15)
colnames(counts) = colnames(tpm) = gsub("/","_", gsub(".RSEM_Quant.genes.results","",files))

save(file="./RSEM_Quant.genes.counts.RData",counts)
save(file="./RSEM_Quant.genes.tpm.RData",tpm)
