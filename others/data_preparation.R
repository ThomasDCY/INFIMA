# data preparation for EM model
# A=A/J, B=B6, C=129, D=NOD, E=NZO, F=Cast, G=PWK, H=WSB
library(data.table)
setwd('~/Dropbox/AlanAttie4/MethodA/')
load('DO_eQTL.RData')
load('ATAC-QTL.RData')

outdir <- '/p/keles/Collab_2014/volumeK/model/run090319/'
do.eqtl <- as.data.table(do.eqtl)

#### YY:39958*8, DO-RNA-seq effect sizes ####
Y <- do.eqtl[,A:H]
Y <- Y[,c(3,1,2,6,4,5,7,8)]
colnames(Y) <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')
Y <- as.matrix(Y) # 39958*8
YY <- Y

#### AA:47062*8, ATAC-seq peak signal ####
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/peaks/IDR_master_peak_list_SNP_with_counts.RData')
countMatScaled <- countMatScaled[snpData$`Which Peak`,] # 47062*16
A <- matrix(0, nrow = 47062, ncol = 8)
for(i in 1:8){
  A[,i] <- (countMatScaled[,i*2-1] + countMatScaled[,i*2])/2
}
colnames(A) <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')
A <- as.matrix(A) # 47062*8
AA <- A

#### FF:47062*1, involved in significant TF binding ####
load('/p/keles/Collab_2014/volumeK/AlanAttie4/atSNP/run080519/atsnp.final.RData')
snps <- unique(atsnp.final$snpid) # 1350 SNPs
FF <- rep(0, 47062)
FF[match(snps, snpData$Name)] <- 1
# sum(FF) 1350

#### EE:ATAC-QTL genotype ####
setwd('/p/keles/Collab_2014/volumeK/eQTL/data')
cisEQTL <- fread('/p/keles/Collab_2014/volumeK/eQTL/output/cisEQTL.txt')
snpData <- as.data.table(snpData)
EE <- as.matrix(snpData[,8:15])
colnames(EE) <- c('129','AJ','B6','Cast','NOD','NZO','PWK','WSB')

outdir <- '/p/keles/Collab_2014/volumeK/model/run090319/'
setwd(outdir)
save(YY,AA,FF,EE,cisEQTL, file = 'input_data.RData')


##################################
rm(list = ls())
library(data.table)
# trinarize data #
outdir <- '/p/keles/Collab_2014/volumeK/model/run090319/'
setwd(outdir)
load('input_data.RData')
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/DO_eQTL/DO_eQTL.RData')
load('/p/keles/Collab_2014/volumeK/AlanAttie4/RData/ATAC_QTL/ATAC-QTL.RData')
do.eqtl <- as.data.table(do.eqtl)

# trinarize YY #
# standardize to [0,1] 
row.min <- apply(YY,1,min)
row.max <- apply(YY,1,max)
for(i in 1:8){
  YY[,i] <- (YY[,i] - row.min)/(row.max - row.min)
}
# set B6 to be zero
B6 <- YY[,3]
for(i in 1:8){
  YY[,i] <- YY[,i] - B6
}
# value < -0.2, becomes -1
# value > 0.2, becomes 1
# else, becomes 0
YY[YY < -0.2] <- -1
YY[YY > 0.2] <- 1
YY[YY >= -0.2 & YY <= 0.2] <- 0


### apply the same technique to normalize AA ###

row.min <- apply(AA,1,min)
row.max <- apply(AA,1,max)
for(i in 1:8){
  AA[,i] <- (AA[,i] - row.min)/(row.max - row.min)
}
# set B6 to be zero
B6 <- AA[,3]
for(i in 1:8){
  AA[,i] <- AA[,i] - B6
}
# value < -0.2, becomes -1
# value > 0.2, becomes 1
# else, becomes 0
AA[AA < -0.2] <- -1
AA[AA > 0.2] <- 1
AA[AA >= -0.2 & AA <= 0.2] <- 0


#### prepare EE ####
### if effect size > 0 then strains with alternative allele has a 1
### otherwise -1

#effectsize <- ifelse(cisEQTL$beta > 0, 1, -1)
#for(i in 1:8){
#  EE[,i] <- EE[,i]/2*effectsize
#}

save(YY,AA,FF,EE,cisEQTL, file = 'input_data_trinary.RData')

