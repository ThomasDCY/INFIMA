#### load the data ####
library(data.table)
library(GenomicAlignments)
rm(list = ls())
outdir <- '/p/keles/Collab_2014/volumeK/model/run090319/'
setwd(outdir)
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/RNA-seq/countmat_13568_protein_coding_UQUA.RData')
load('input_data_trinary_filtered.RData')

# filter out different chr cases
indicator <- do.eqtl$chr == do.eqtl$qtl_chr # 10396 rows left
YY <- YY[indicator,] # 10396*8
do.eqtl <- do.eqtl[indicator]


#### filter out those do.eqtl rows that do not have overlapping SNPs ####
G <- nrow(do.eqtl) # 10396
E.g <- vector('list', G)
Y.g <- vector('list', G)
A.g <- vector('list', G)
F.g <- vector('list', G)
D.g <- vector('list', G)
B.g <- vector('list', G)
snp_index <- vector('list', G)

for(g in 1:G){
  print(g)
  Y.g[[g]] <- YY[g,]
  B.g[[g]] <- countmat.n.s[do.eqtl$index_founder[g],]
  # find associated SNPs
  window <- 1000000
  snp.loc <- GRanges(seqnames = do.eqtl$qtl_chr[g],
                     ranges = IRanges(start = do.eqtl$qtl_pos[g] - window, 
                                      end = do.eqtl$qtl_pos[g] + window),
                     strand = '*')
  index <- subjectHits(findOverlaps(snp.loc, snps))
  snp_index[[g]] <- index
  
  if(!identical(index, integer(0))){ # there exist at least one overlap
    A.g[[g]] <- AA[index,]
    F.g[[g]] <- FF[index]
    E.g[[g]] <- EE[index,] # genotype.g
    
    if(is.null(dim(E.g[[g]]))){ # only one SNP found, correct the dimension
      E.g[[g]] <- t(as.matrix(E.g[[g]]))
      A.g[[g]] <- t(as.matrix(A.g[[g]]))
    }
    
    # get effect size E.g
    exp.g <- countmat.n[do.eqtl$index_founder[g],] # gene expression
    atac.qtls.g <- genotype[index] # regressors
    atac.qtls.g <- t(as.matrix(atac.qtls.g[,2:92]))
    
    # run eQTL and trinarize E.g
    for(i in 1:length(index)){
      lm <- lm(exp.g ~ atac.qtls.g[,i])
      estimate <- summary(lm)$coefficients[2,1]
      pval <- summary(lm)$coefficients[2,4]
      M <- 0
      if(pval < 0.05){
        M <- ifelse(estimate > 0, 1, -1)
      }
      E.g[[g]][i,] <- E.g[[g]][i,]/2*M
    }
    
    # compute the absolute distance
    D.g[[g]] <- matrix(0, nrow(E.g[[g]]), ncol(E.g[[g]]))
    for(i in 1:nrow(E.g[[g]])){
      D.g[[g]][i,] <- abs(E.g[[g]][i,] - Y.g[[g]])
    }
    
  }
}

p.g <- rep(0, G)
for(g in 1:G){
  if(!identical(snp_index[[g]], integer(0))){
    p.g[g] <- length(snp_index[[g]])
  }
}

# > sum(p.g == 0)
# [1] 205

save(E.g, Y.g, A.g, F.g, D.g, B.g, snp_index, p.g, file = './EM_input_main.RData')
###############################################


