#### compute pseudo counts EM ####
#### load the data ####
library(data.table)
library(GenomicAlignments)

rm(list = ls())
outdir <- '/p/keles/Collab_2014/volumeK/model/run090319/'
setwd(outdir)
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/RNA-seq/countmat_13568_protein_coding_UQUA.RData')
load('input_data_filtered.RData') # do not use the trinarized data for priors
# YY,AA,FF,EE,cisEQTL, genotype, do.eqtl, snps, snpData, countmat.n.s
load('EM_input_main.RData') # E.g, Y.g, A.g, F.g, D.g, snp_index, p.g

# filter out different chr cases
indicator <- do.eqtl$chr == do.eqtl$qtl_chr # 10936 rows left
YY <- YY[indicator,] # 10396*8
do.eqtl <- do.eqtl[indicator]
G <- nrow(do.eqtl) # 10396

# initialize the pseudo counts
PI <- vector('list', G)
# compute pseudo_count components
dist <- vector('list', G)
cor.A.E <- vector('list', G)
cor.A.R <- vector('list', G)

for(g in 1:G){
  print(g)
  if(p.g[g] > 0){
    PI[[g]] <- rep(0, p.g[g])
    dist[[g]] <- rep(0, p.g[g])
    cor.A.E[[g]] <- rep(0, p.g[g])
    cor.A.R[[g]] <- rep(0, p.g[g])
    
    # Modification: distance should be changed to the distance between TSS and ATAC-QTL
    if(do.eqtl$strand[g] == '+'){
      tss <- c(do.eqtl$start[g] - 5000, do.eqtl$start[g] + 500)
    }
    else{
      tss <- c(do.eqtl$end[g] - 500, do.eqtl$end[g] + 5000)
    }
    
    for(k in 1:p.g[g]){
      # distance score
      # <= 250K: score = 0.5
      # > 250K & <= 1M: score from 0.5 to 0.25
      snp_pos <- start(snps[snp_index[[g]][k]])
      d <- min(abs(tss[1] - snp_pos), abs(tss[2] - snp_pos))/1e5
      if(d <= 2.5){
        dist[[g]][k] <- 0.5
      }
      else{
        dist[[g]][k] <- 5/12/(d - 5/3)
      }
      # correlation of founder RNA-seq vs ATAC-seq score
      cor.A.R[[g]][k] <- cor(countmat.n.s[do.eqtl$index_founder[g],], AA[snp_index[[g]][k],])
      # correlation of ATAC-seq and SNP score
      cor.A.E[[g]][k] <- cor(AA[snp_index[[g]][k],], EE[snp_index[[g]][k],])
      
      # distance prior
      PI[[g]][k] <- dist[[g]][k]
      # footprint analysis score
      PI[[g]][k] <- PI[[g]][k] + F.g[[g]][k]
      # correlation of founder RNA-seq vs ATAC-seq score
      PI[[g]][k] <- PI[[g]][k] + abs(cor.A.R[[g]][k])
      # correlation of ATAC-seq and SNP score
      PI[[g]][k] <- PI[[g]][k] + abs(cor.A.E[[g]][k])
    }
  }
}

save(PI, file = 'pseudo_counts.RData')
save(dist, F.g, cor.A.R, cor.A.E, file = 'pseudo_counts_components.RData')

