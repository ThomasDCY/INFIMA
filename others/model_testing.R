library(data.table)
library(GenomicAlignments)

outdir <- '/p/keles/Collab_2014/volumeK/model/run032719/'
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

save(YY,AA,FF,EE,cisEQTL, genotype, do.eqtl, snps, snpData, file = 'input_data_trinary_filtered.RData')



### focus on the first row of do.eqtl data
row.no <- 1
Y.g <- YY[row.no,]
# find associated SNPs
window <- 1000000
snp.loc <- GRanges(seqnames = do.eqtl$qtl_chr[row.no],
                   ranges = IRanges(start = do.eqtl$qtl_pos[row.no] - window, 
                                    end = do.eqtl$qtl_pos[row.no] + window),
                   strand = '*')
index <- subjectHits(findOverlaps(snp.loc, snps))
A.g <- AA[index,]
F.g <- FF[index]
E.g <- EE[index,]
g <- do.eqtl$ensembl[row.no] # current gene
exp.g <- countmat.n[grep(g, rownames(countmat)),] # gene expression
atac.qtls.g <- genotype[index] # regressors
# get effect size
atac.qtls.g <- t(as.matrix(atac.qtls.g[,2:92]))

# trinarize E.g and run eQTL
for(i in 1:length(index)){
  lm <- lm(exp.g ~ atac.qtls.g[,i] - 1)
  estimate <- summary(lm)$coefficients[1]
  pval <- summary(lm)$coefficients[4]
  M <- 0
  if(pval < 0.05){
    M <- ifelse(estimate > 0, 1, -1)
  }
  E.g[i,] <- E.g[i,]/2*M
}

### Now we have Y.g, A.g, F.g, E.g ready to run EM.
# set parameters and simulate Y.g given A.g, F.g, E.g
theta.true <- c(0, 0, 0, 0, 0, 0.9, 0.1, 0, 0, 0)
alpha.true <- c(0.7,0.3,0)
beta.true <- c(0.6,0.4,0)
p <- length(theta.true)
set.seed(2019)

Y.g.simulate <- function(R.g, alpha.true, beta.true){
  p1 <- sum(R.g == -1)/length(R.g)
  p2 <- sum(R.g == 0)/length(R.g)
  p3 <- sum(R.g == 1)/length(R.g)
  Y.g.sim <- rep(0,8)
  for(i in 1:8){
    if(i %in% c(1,2,3,5,6)){
      params = alpha.true
    }
    else{
      params = beta.true
    }
    
    if(R.g[i] == -1){
      probs = c(1-params[3]/(p1+p3)-params[2]/(p2+1),
                params[2]/(p2+1),
                params[3]/(p1+p3))
      Y.g.sim[i] <- which(as.numeric(rmultinom(1,1,probs)) == 1) - 2
    }
    if(R.g[i] == 0){
      probs = c(params[2]/(p2+1),
                1-2*params[2]/(p2+1),
                params[2]/(p2+1))
      Y.g.sim[i] <- which(as.numeric(rmultinom(1,1,probs)) == 1) - 2
    }
    if(R.g[i] == 1){
      probs = c(params[3]/(p1+p3),
                params[2]/(p2+1),
                1-params[3]/(p1+p3)-params[2]/(p2+1))
      Y.g.sim[i] <- which(as.numeric(rmultinom(1,1,probs)) == 1) - 2
    }
  }
  return(Y.g.sim)
}

for(iter in 1:1000){
  Z.g.true <- as.vector(rmultinom(1,1,theta.true))
  R.g <- E.g[as.logical(Z.g.true),]
  Y.g <- Y.g.simulate(R.g, alpha.true, beta.true)

  ### run EM ###
  ### try to recover theta, alpha and beta ###
  ### initial values
  alpha <- c(0.65,0.35,0)
  beta <- c(0.55,0.45,0)
  theta <- rep(1/p, p)
  Z.g <- theta
  
  D.g <- matrix(0, nrow(E.g), ncol(E.g))
  for(i in 1:nrow(E.g)){
    D.g[i,] <- abs(E.g[i,] - Y.g)
  }
  D1.g <- D.g[,c(1,2,3,5,6)]
  D2.g <- D.g[,c(4,7,8)]
  n10 <- apply(D1.g, 1, function(x){sum(x==0)})
  n11 <- apply(D1.g, 1, function(x){sum(x==1)})
  n12 <- apply(D1.g, 1, function(x){sum(x==2)})
  n20 <- apply(D2.g, 1, function(x){sum(x==0)})
  n21 <- apply(D2.g, 1, function(x){sum(x==1)})
  n22 <- apply(D2.g, 1, function(x){sum(x==2)})
  
  for(em.iter in 1:100){
    # E-step
    for(k in 1:p){
      Z.g[k] <- theta[k]*prod(c(alpha,beta)^c(n10[k],n11[k],n12[k],n20[k],n21[k],n22[k]))
    }
    Z.g <- Z.g/sum(Z.g)
    
    # M-step
    PI <- rep(0,p)
    for(k in 1:p){
      PI[k] <- F.g[k] + abs(cor(Y.g, A.g[k,]))
    }
    theta <- PI + Z.g
    theta <- theta/sum(theta)
    nn10 <- sum(Z.g*n10)
    nn11 <- sum(Z.g*n11)
    nn12 <- sum(Z.g*n12)
    nn20 <- sum(Z.g*n20)
    nn21 <- sum(Z.g*n21)
    nn22 <- sum(Z.g*n22)
    alpha <- c(nn10 + p, nn11, nn12)
    beta <- c(nn20 + p, nn21, nn22)
    alpha <- alpha/sum(alpha)
    beta <- beta/sum(beta)
  }

  
}
