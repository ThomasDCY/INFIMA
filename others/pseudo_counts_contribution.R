### quantify the contribution of each piece of information in the prior ###
library(data.table)
library(GenomicAlignments)
library(ggplot2)
setwd('~/Dropbox/EM/run090319/')

load('validation/EM_results_090419.RData')
load('pseudo_counts_components.RData')
cor.A.E <- lapply(cor.A.E, function(x) {
  if(length(x) > 0) return(abs(x))
  else return(x)
})
cor.A.R <- lapply(cor.A.R, function(x) {
  if(length(x) > 0) return(abs(x))
  else return(x)
})

##### get the rank of each piece of information for the predicted causal SNP
which.max(Z.g[[1]])
len <- sapply(Z.g, length)
Z.g[[366]]
which.max(NULL)

# rank(): the higher the rank, the larger the true value
getRank <- function(i, quantity){
  if(length(Z.g[[i]]) > 0){
    ind <- which.max(Z.g[[i]])
    return(rank(quantity[[i]])[ind]/length(quantity[[i]]))
  }
  else{
    return(NA)
  }
}

cor.A.E.rank <- sapply(1:length(Z.g), function(x) getRank(x, cor.A.E))
summary(cor.A.E.rank)
cor.A.R.rank <- sapply(1:length(Z.g), function(x) getRank(x, cor.A.R))
summary(cor.A.R.rank)
F.g.rank <- sapply(1:length(Z.g), function(x) getRank(x, F.g))
summary(F.g.rank)
dist.rank <- sapply(1:length(Z.g), function(x) getRank(x, dist))
summary(dist.rank)

dt <- data.table(cor.A.E.rank = cor.A.E.rank,
                 cor.A.R.rank = cor.A.R.rank,
                 F.g.rank = F.g.rank,
                 dist.rank = dist.rank)
dt <- na.omit(dt)
dt

# for each DO gene, rank the contribution from each piece of information
dt$id <- 1:nrow(dt)
dt2 <- melt(dt, id.vars = 'id')
dt3 <- dt2[, variable[which.max(value)], by = id]
table(dt3$V1)/nrow(dt3)
# cor.A.E.rank cor.A.R.rank     F.g.rank    dist.rank 
# 0.2013541    0.3507016    0.2444314    0.2035129

dt4 <- dt2[, max(value), by = id]
dt4$variable <- dt3$V1
colnames(dt4) <- c('id', 'Score', 'Variable')
ggplot(dt4, aes(x = Variable, y = Score, color = Variable)) +
  geom_boxplot()

V.g

