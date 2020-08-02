#' @name compute_prior
#' @title Compute Dirichlet prior
#' @description For the candidate SNPs of each row of DO-eQTL data, 
#' keep track of their the distance scores, footprinting annotations,
#' correlations between local ATAC-seq signal and effect size,
#' correlations between local ATAC-seq signal and gene expression.
#' @param raw_data The output from \code{raw_input_data}
#' @param model_data The output from  \code{model_input_data}
#' @param n_cores The number of cores for parallel computing.
#' @param verbose Print messages or not. 
#' Default = \code{detectCores() - 2}.
#' 
#' 
#' @return A prior data object.
#' 
#' \tabular{ll}{
#' \code{PI} \tab The overall prior term. \cr
#' \code{dist} \tab The distance score. \cr
#' \code{F.g} \tab The footprint annotation. \cr
#' \code{cor.A.B} \tab The Correlation betwwen ATAC-seq signal and founder effect size. \cr
#' \code{cor.A.E} \tab The Correlation between ATAC-seq signal and gene expression. \cr
#' }
#' @examples
#' data('example-10-genes')
#' raw_data <- raw_input_data(dt1, dt2, dt3)
#' model_data <- model_input_data(raw_data)
#' prior <- compute_prior(raw_data, model_data)
#' @seealso \code{\link{raw_input_data}}, \code{\link{model_input_data}}.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import parallel
#' @import doParallel
#' @import foreach
#' @rawNamespace import(data.table, except = shift)
#' @export
compute_prior <- function(raw_data = NULL,
                                model_data = NULL,
                                n_cores = detectCores() - 2,
                                verbose = TRUE) {
  if (is.null(raw_data)) {
    stop('Error: please input the raw_data.')
  }
  
  if (class(raw_data) != 'raw_data') {
    stop('Error: please input the S3 class raw_data.')
  }
  
  if (is.null(model_data)) {
    stop('Error: please input the model_data.')
  }
  
  if (class(model_data) != 'model_data') {
    stop('Error: please input the S3 class model_data.')
  }
  
  time.start <- proc.time()
  
  AA <- raw_data$AA
  EE <- raw_data$EE
  BB <- raw_data$BB
  do.eqtl <- raw_data$do.eqtl
  snpData <- raw_data$snpData
  
  F.g <- model_data$F.g
  snp_index <- model_data$snp_index
  p.g <- model_data$p.g
  window <- model_data$window
  
  stopifnot(exprs = {
    ! is.null(AA)
    ! is.null(EE)
    ! is.null(BB)
    ! is.null(do.eqtl)
    ! is.null(snpData)
    ! is.null(F.g)
    ! is.null(snp_index)
    ! is.null(p.g)
    ! is.null(window)
  })
  
  
  G <- length(p.g)
  
  # # initialize the pseudo counts
  # PI <- vector('list', G)
  # # compute pseudo_count components
  # dist <- vector('list', G)
  # cor.A.E <- vector('list', G)
  # cor.A.B <- vector('list', G)
  
  registerDoParallel(cores = n_cores)
  
  if(verbose){
    message(paste('n_cores =', n_cores))
  }
  
  r <- foreach(g = 1:G) %dopar% {
    # print(g)
    
    # if no local-ATAC-QTL, then return NULL
    PI <- NULL
    dist <- NULL
    cor.A.E <- NULL
    cor.A.B <- NULL
    
    if (p.g[g] > 0) {
      PI <- rep(0, p.g[g])
      dist <- rep(0, p.g[g])
      cor.A.E <- rep(0, p.g[g])
      cor.A.B <- rep(0, p.g[g])
      
      # Distance between TSS and ATAC-QTL
      if (do.eqtl$strand[g] == '+') {
        tss <- c(do.eqtl$start[g] - 5000, do.eqtl$start[g] + 500)
      }
      else{
        tss <- c(do.eqtl$end[g] - 500, do.eqtl$end[g] + 5000)
      }
      
      for (k in 1:p.g[g]) {
        # distance score when window = 1M
        # <= 250K: score = 0.5
        # > 250K & <= 1M: score from 0.5 to 0.25
        # snp_pos <- start(snps[snp_index[[g]][k]])
        snp_pos <- snpData[snp_index[[g]][k]]$snp_pos
        d <- min(abs(tss[1] - snp_pos), abs(tss[2] - snp_pos)) / 1e5
        d <- d * window / 1000000
        if (d <= 2.5) {
          dist[k] <- 0.5
        }
        else{
          dist[k] <- 5 / 12 / (d - 5 / 3)
        }
        
        # correlation of founder RNA-seq vs ATAC-seq score
        cor.A.B[k] <- cor(BB[g, ], AA[snp_index[[g]][k], ])
        # correlation of ATAC-seq and allelic effect of founder gene
        cor.A.E[k] <-
          cor(AA[snp_index[[g]][k], ], EE[snp_index[[g]][k], ])
        
        # distance prior
        PI[k] <- dist[k]
        # footprint analysis score
        PI[k] <- PI[k] + F.g[[g]][k]
        # correlation of founder RNA-seq vs ATAC-seq score
        PI[k] <- PI[k] + abs(cor.A.B[k])
        # correlation of ATAC-seq and SNP genotype
        PI[k] <- PI[k] + abs(cor.A.E[k])
      }
    }
    
    return(list(
      PI = PI,
      dist = dist,
      cor.A.B = cor.A.B,
      cor.A.E = cor.A.E
    ))
  }
  
  PI <- lapply(r, function(x) x$PI)
  dist <- lapply(r, function(x) x$dist)
  cor.A.B <- lapply(r, function(x) x$cor.A.B)
  cor.A.E <- lapply(r, function(x) x$cor.A.E)
  
  prior <- list(
    PI = PI,
    dist = dist,
    F.g = F.g,
    cor.A.B = cor.A.B,
    cor.A.E = cor.A.E
  )
  
  if(verbose){
    cat("Time taken", proc.time()[3] - time.start[3])
  }
  
  class(prior) <- 'prior'
  
  return(prior)
}





