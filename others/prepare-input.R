rm(list = ls())
library(data.table)
setwd('~/Documents/GitHub/INFIMA/data/')

load('~/Dropbox/EM/run090319/input_data_filtered.RData')
load('~/Dropbox/RNA-seq/countmat_13568_protein_coding_UQUA.RData')
# filter out different chr cases
indicator <- do.eqtl$chr == do.eqtl$qtl_chr # 10396 rows left
YY <- YY[indicator,] # 10396*8
do.eqtl <- do.eqtl[indicator]
index.founder <- do.eqtl$index_founder

#### INPUT DATASET 1: DO-eQTL data
# A=A/J, B=B6, C=129, D=NOD, E=NZO, F=Cast, G=PWK, H=WSB
Y <- do.eqtl[,A:H]
Y <- Y[,c(3,1,2,6,4,5,7,8)]
colnames(Y) <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')
do.eqtl <- cbind(do.eqtl[,-(A:H)], Y)
do.eqtl$cis <- NULL
do.eqtl$index_founder <- NULL
dt1 <- do.eqtl

#### INPUT DATASET 2: RNA-seq data
rna.seq <- as.data.table(countmat.n)[index.founder]
samples <- sapply(colnames(rna.seq), function(x){
  strsplit(x, '-')[[1]][1]
})
names(samples) <- NULL
colnames(rna.seq) <- samples
dt2 <- rna.seq


#### INPUT DATASET 3: ATAC-seq data and the local-ATAC-QTL data
colnames(snpData)
atac.seq <- as.data.table(AA)
colnames(atac.seq) <- paste0(colnames(atac.seq),'-ATAC-seq')
genotype <- as.data.table(EE)
colnames(genotype) <- paste0(colnames(genotype),'-genotype')
snp.info <- as.data.table(snpData)[,c('Name','Chr','SNP','REF','ALT')]
colnames(snp.info) <- c('snp_id','chr','snp_pos','ref','alt')
dt3 <- cbind(snp.info, genotype, atac.seq)
dt3$footprint <- FF


save(dt1, dt2, dt3, file = 'inputs.RData')
write.table(dt1[chr == 'chr1'], file = 'input1-chr1.csv', sep = ',',
            quote = F, row.names = F, col.names = T)
write.table(dt2[dt1$chr == 'chr1'], file = 'input2-chr1.csv', sep = ',',
            quote = F, row.names = F, col.names = T)
write.table(dt3[chr == 'chr1'], file = 'input3-chr1.csv', sep = ',',
            quote = F, row.names = F, col.names = T)
