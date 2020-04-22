### input_data.RData to input_data_filtered.RData ###
### input_data_trinary.RData to input_data_trinary_filtered.RData ###
library(data.table)
library(GenomicAlignments)
library(parallel)

rm(list = ls())
outdir <- '/p/keles/Collab_2014/volumeK/model/run090319/'
setwd(outdir)
load('input_data.RData')
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/DO_eQTL/DO_eQTL.RData')
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/ATAC_QTL/ATAC-QTL.RData')
do.eqtl <- as.data.table(do.eqtl)

snps <- GRanges(seqnames = snpData$Chr,
                ranges = IRanges(start = snpData$SNP, end = snpData$SNP),
                strand = '*')
snp.names <- snpData$Name
snp.is.cisEQTL <- snp.names %in% unique(cisEQTL$SNP)

# load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/peaks/IDR_master_peak_list_SNP_with_counts.RData')
# widths <- width(peaks)
# widths <- widths[snpData$`Which Peak`]
# snp.is.in.narrowPeak <- widths <= 3000
# 
# indicator <- snp.is.cisEQTL & snp.is.in.narrowPeak


snps <- snps[snp.is.cisEQTL]
AA <- AA[snp.is.cisEQTL,] # 34711 ATAC-QTLs are significant cis-EQTLs
FF <- FF[snp.is.cisEQTL] # 966/34711
EE <- EE[snp.is.cisEQTL,]
snpData <- snpData[snp.is.cisEQTL,]
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/RNA-seq/countmat_13568_protein_coding_UQUA.RData')
genotype <- fread('/p/keles/Collab_2014/volumeK/eQTL/data/SNP.txt')
genotype <- genotype[snp.is.cisEQTL]

### make sure that we only use the DO data with support from founder RNA-seq
do.genes <- as.character(do.eqtl$ensembl)
index.founder <- mclapply(1:length(do.genes), function(x){
  print(x)
  res = grep(do.genes[x], rownames(countmat))
  ifelse(identical(res, integer(0)), 0, res)
}, mc.cores = 20)

# table(unlist(lapply(index.founder, length))) good

index.founder <- unlist(index.founder)
# > sum(index.founder != 0)
# [1] 27090
do.eqtl <- do.eqtl[index.founder != 0,]
YY <- YY[index.founder != 0,]
index.founder <- index.founder[index.founder != 0]
do.eqtl$index_founder <- index.founder

# get countmat.n.s, 13568*8
countmat.n.s <- matrix(0, nrow = 13568, ncol = 8)
colnames(countmat.n.s) <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')
strains <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')
colnames.s <- unlist(lapply(colnames(countmat.n), function(x){strsplit(x,'-')[[1]][1]}))

for(i in 1:8){
  countmat.n.s[,i] <- apply(countmat.n[,colnames.s == strains[i]], 1, mean)
}


save(YY,AA,FF,EE,cisEQTL, genotype, do.eqtl, snps, snpData, countmat.n.s, file = 'input_data_filtered.RData')





##################
rm(list = ls())
outdir <- '/p/keles/Collab_2014/volumeK/model/run090319/'
setwd(outdir)
load('input_data_trinary.RData')
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/DO_eQTL/DO_eQTL.RData')
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/ATAC_QTL/ATAC-QTL.RData')
do.eqtl <- as.data.table(do.eqtl)

snps <- GRanges(seqnames = snpData$Chr,
                ranges = IRanges(start = snpData$SNP, end = snpData$SNP),
                strand = '*')
snp.names <- snpData$Name
snp.is.cisEQTL <- snp.names %in% unique(cisEQTL$SNP)
snps <- snps[snp.is.cisEQTL]
AA <- AA[snp.is.cisEQTL,] # 34711 ATAC-QTLs are significant cis-EQTLs
FF <- FF[snp.is.cisEQTL] # 6744/34711
EE <- EE[snp.is.cisEQTL,]
snpData <- snpData[snp.is.cisEQTL,]
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/RNA-seq/countmat_13568_protein_coding_UQUA.RData')
genotype <- fread('/p/keles/Collab_2014/volumeK/eQTL/data/SNP.txt')
genotype <- genotype[snp.is.cisEQTL]

### make sure that we only use the DO data with support from founder RNA-seq
do.genes <- as.character(do.eqtl$ensembl)
index.founder <- mclapply(1:length(do.genes), function(x){
  print(x)
  res = grep(do.genes[x], rownames(countmat))
  ifelse(identical(res, integer(0)), 0, res)
}, mc.cores = 20)

# table(unlist(lapply(index.founder, length))) good

index.founder <- unlist(index.founder)
# > sum(index.founder != 0)
# [1] 27090
do.eqtl <- do.eqtl[index.founder != 0,]
YY <- YY[index.founder != 0,]
index.founder <- index.founder[index.founder != 0]
do.eqtl$index_founder <- index.founder

# get countmat.n.s, 13568*8
countmat.n.s <- matrix(0, nrow = 13568, ncol = 8)
colnames(countmat.n.s) <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')
strains <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')
colnames.s <- unlist(lapply(colnames(countmat.n), function(x){strsplit(x,'-')[[1]][1]}))

for(i in 1:8){
  countmat.n.s[,i] <- apply(countmat.n[,colnames.s == strains[i]], 1, mean)
}


save(YY,AA,FF,EE,cisEQTL, genotype, do.eqtl, snps, snpData, countmat.n.s, file = 'input_data_trinary_filtered.RData')


