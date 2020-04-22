##### distance summary #####
# summary of distances between ATAC-QTLs and genes
library(data.table)
library(GenomicAlignments)

outdir <- '/p/keles/Collab_2014/volumeK/model/run032719/'
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

# distance to gene
dist.gene <- vector('list', G)
# distance to eqtl marker
dist.mark <- vector('list', G)
for(g in 1:G){
  print(g)
  if(p.g[g] > 0){
    dist.gene[[g]] <- rep(0, p.g[g])
    dist.mark[[g]] <- rep(0, p.g[g])
    for(k in 1:p.g[g]){
      gene_start <- do.eqtl$start[g]
      gene_end <- do.eqtl$end[g]
      qtl_pos <- do.eqtl$qtl_pos[g]
      snp_pos <- start(snps[snp_index[[g]][k]])
      # if snp is outside of the gene
      if(snp_pos < gene_start | snp_pos > gene_end){
        dist.gene[[g]][k] <- min(abs(gene_start-snp_pos), abs(gene_end-snp_pos))
      }
      dist.mark[[g]][k] <- abs(qtl_pos - snp_pos)
    }
  }
}

save(dist.gene, dist.mark, file = './dist_summary.RData')

indicator <- unlist(lapply(dist, is.null))
dist2 <- dist[!indicator]
dist.min <- unlist(lapply(dist2, min))
summary(dist.min)
dist.max <- unlist(lapply(dist2, max))
summary(dist.max)

# > summary(dist.min)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 1      1132     14439   1864540     82252 170072559
# > summary(dist.max)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 64    891899   1027821   2897784   1332322 171630154

