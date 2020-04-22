#' Compute pseudocount prior
#' @param raw_data The output from @seealso raw_input_data
#' @param model_data The output from @seealso model_input_data
#' @param window The window around DO-eQTL marker contains candidate SNPs.
#' Default = 1 Mbp
#' @return A list of data objects.
#' PI: The pseudocount term.
#' dist: The distance score.
#' F.g: The footprint information.
#' cor.A.R: correlation betwwen ATAC-seq signal and effect size.
#' cor.A.E: correlation between ATAC-seq signal and gene expression.
#' @example
#' pseudocount <- compute_pseudocount(raw_data, model_data)
#' @import data.table
#' @export
compute_pseudocount <- function(raw_data = NULL,
                                model_data = NULL,
                                window = 1000000) {
  if (is.null(raw_data)) {
    stop('Error: please input the raw_data.')
  }
  
  if (is.null(model_data)) {
    stop('Error: please input the model_data.')
  }
  
  AA <- raw_data$AA
  EE <- raw_data$EE
  BB <- raw_data$BB
  do.eqtl <- raw_data$do.eqtl
  snpData <- raw_data$snpData
  
  F.g <- model_data$F.g
  snp_index <- model_data$snp_index
  p.g <- model_data$p.g
  
  
  stopifnot(exprs = {
    !is.null(AA)
    ! is.null(EE)
    ! is.null(BB)
    ! is.null(do.eqtl)
    ! is.null(snpData)
    ! is.null(F.g)
    ! is.null(snp_index)
    ! is.null(p.g)
  })
  
  
  G <- length(p.g)
  
  
  
  
  # initialize the pseudo counts
  PI <- vector('list', G)
  # compute pseudo_count components
  dist <- vector('list', G)
  cor.A.E <- vector('list', G)
  cor.A.R <- vector('list', G)
  
  for (g in 1:G) {
    # print(g)
    if (p.g[g] > 0) {
      PI[[g]] <- rep(0, p.g[g])
      dist[[g]] <- rep(0, p.g[g])
      cor.A.E[[g]] <- rep(0, p.g[g])
      cor.A.R[[g]] <- rep(0, p.g[g])
      
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
          dist[[g]][k] <- 0.5
        }
        else{
          dist[[g]][k] <- 5 / 12 / (d - 5 / 3)
        }
        
        # correlation of founder RNA-seq vs ATAC-seq score
        cor.A.R[[g]][k] <- cor(BB[g, ], AA[snp_index[[g]][k], ])
        # correlation of ATAC-seq and SNP genotype
        cor.A.E[[g]][k] <-
          cor(AA[snp_index[[g]][k], ], EE[snp_index[[g]][k], ])
        
        # distance prior
        PI[[g]][k] <- dist[[g]][k]
        # footprint analysis score
        PI[[g]][k] <- PI[[g]][k] + F.g[[g]][k]
        # correlation of founder RNA-seq vs ATAC-seq score
        PI[[g]][k] <- PI[[g]][k] + abs(cor.A.R[[g]][k])
        # correlation of ATAC-seq and SNP genotype
        PI[[g]][k] <- PI[[g]][k] + abs(cor.A.E[[g]][k])
      }
    }
  }
  
  return(list(
    PI = PI,
    dist = dist,
    F.g = F.g,
    cor.A.R = cor.A.R,
    cor.A.E = cor.A.E
  ))
}





