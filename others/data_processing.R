# Build a data table of all candidate ATAC-QTLs
# for each gene g, do the following
# columns:
# 1. 


library(data.table)
library(GenomicAlignments)

outdir <- '/p/keles/Collab_2014/volumeK/model/run032719/'
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
index.founder <- lapply(1:length(do.genes), function(x){
  print(x)
  res = grep(do.genes[x], rownames(countmat))
  ifelse(identical(res, integer(0)), 0, res)
})

# table(unlist(lapply(index.founder, length))) good

index.founder <- unlist(index.founder)
# > sum(index.founder != 0)
# [1] 27090
do.eqtl <- do.eqtl[index.founder != 0,]
index.founder <- index.founder[index.founder != 0]

YY <- do.eqtl[,A:H]
YY <- YY[,c(3,1,2,6,4,5,7,8)]
colnames(YY) <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')
YY <- as.matrix(YY) # 27090*8


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







### focus on the first row of do.eqtl data
#load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/RNA-seq/countmat_13568_protein_coding_UQUA.RData')
load('input_data_filtered.RData')


row.no <- 1

dt <- NULL
no_overlap <- 0
for(row.no in 1:nrow(do.eqtl)){
print(row.no)
Y.g <- YY[row.no,]
# find associated SNPs
window <- 1000000
snp.loc <- GRanges(seqnames = do.eqtl$qtl_chr[row.no],
                   ranges = IRanges(start = do.eqtl$qtl_pos[row.no] - window, 
                                    end = do.eqtl$qtl_pos[row.no] + window),
                   strand = '*')
index <- subjectHits(findOverlaps(snp.loc, snps))

if(!identical(index, integer(0))){
  A.g <- AA[index,]
  F.g <- FF[index]
  E.g <- EE[index,] # genotype.g
  # g <- do.eqtl$ensembl[row.no] # current gene
  exp.g <- countmat.n.s[do.eqtl$index_founder[row.no],] # gene expression, 8*1
  
  if(is.null(dim(E.g))){ # only one SNP found
    E.g <- t(as.matrix(E.g))
    A.g <- t(as.matrix(A.g))
  }
  # founder RNA-seq vs ATAC-seq
  c1 <- apply(A.g, 1, function(x) cor(x,exp.g))
  # ATAC-seq vs SNP
  c2 <- unlist(lapply(1:nrow(A.g), function(x) cor(A.g[x,], E.g[x,])))
  # DO effect size vs ATAC-seq
  c3 <- apply(A.g, 1, function(x) cor(x, Y.g))
  # DO effect size vs SNP
  c4 <- apply(E.g, 1, function(x) cor(x, Y.g))
  
  dt.g <- data.table(atSNP = F.g,
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    c4 = c4,
                    row = row.no)
  
  dt <- rbind(dt, dt.g)
}
else{
  no_overlap <- no_overlap + 1
}
}

save(dt, file = 'dt_all.RData')


load('~/Dropbox/EM/dt_all.RData')
library(Rtsne)
tsne <- Rtsne(dt[,2:5], check_duplicates = FALSE, pca = FALSE, perplexity=30, theta=0.5, dims=2)

t.sne <- data.table(x = tsne$Y[,1], y = tsne$Y[,2], atSNP = dt$atSNP)

g <- ggplot(t.sne, aes(x = x, y = y, color = atSNP)) + geom_point(size = 2) + scale_color_manual(values = colors) 
g <- g + labs(x = 't-SNE1', y = 't-SNE2', color = 'atSNP')

pdf('./tSNE.pdf', height = 5, width = 5)
g + theme_classic()
dev.off()

library(ggplot2)
pca <- prcomp(dt[,2:5])
dt.pca <- data.table(PC2 = pca$x[,2], PC3 = pca$x[,3], PC1 = ifelse(pca$x[,1] > 0, 1, -1))
g <- ggplot(dt.pca, aes(x = PC2, y = PC3, color = PC1)) + geom_point(size = 0.2)
g <- g + labs(x = 'PC2', y = 'PC3', color = 'PC1')
png('~/Dropbox/EM/PC2_PC3.png', width = 600, height = 600)
g + theme_classic()
dev.off()


dt.pca <- data.table(PC2 = pca$x[,2], PC4 = pca$x[,4], PC1 = ifelse(pca$x[,1] > 0, 1, -1))
g <- ggplot(dt.pca, aes(x = PC2, y = PC4, color = PC1)) + geom_point(size = 0.2)
g <- g + labs(x = 'PC2', y = 'PC4', color = 'PC1')
png('~/Dropbox/EM/PC2_PC4.png', width = 600, height = 600)
g + theme_classic()
dev.off()

dt.pca <- data.table(PC3 = pca$x[,3], PC4 = pca$x[,4], PC1 = ifelse(pca$x[,1] > 0, 1, -1))
g <- ggplot(dt.pca, aes(x = PC3, y = PC4, color = PC1)) + geom_point(size = 0.2)
g <- g + labs(x = 'PC3', y = 'PC4', color = 'PC1')
png('~/Dropbox/EM/PC3_PC4.png', width = 600, height = 600)
g + theme_classic()
dev.off()
