#' @name snp_link_gene
#' @title Link SNPs to effector genes
#'
#' @description Analyze the outputs from INFIMA and link local-ATAC-QTLs to effector genes
#'
#' @param infima The output from \code{model_fitting}
#' @param model_data The output from \code{\link{model_input_data}}.
#' @param fdr The False Discovery Rate cutoff with direct posterior probability
#' approach (Newton, 2004). Default 0.1.
#' @param cum.pprob The maximum cumulative posterior probability for
#' the target credible set. Default 0.8. Take values between 0 and 1.
#' @param cred.set The maximum proportion of candidate SNPs allowed
#' in the credible set. Default 0.5. Take values between 0 and 1.
#' @return A data table with the following columns.
#' \tabular{ll}{
#' \code{snp_id} \tab The SNP ID. \cr
#' \code{chr} \tab The chromosome of the SNP. \cr
#' \code{snp_pos} \tab The locus of the SNP. \cr
#' \code{ref} \tab The reference allele of the SNP. \cr
#' \code{alt} \tab The alternative allele of the SNP. \cr
#' \code{ensembl} \tab The ensemble ID of the effector gene. \cr
#' \code{symbol} \tab The gene symbol of the effector gene. \cr
#' }
#' @examples
#' data('example-10-genes')
#' raw_data <- raw_input_data(dt1, dt2, dt3)
#' model_data <- model_input_data(raw_data)
#' pseudocount <- compute_pseudocount(raw_data, model_data)
#' infima <- model_fitting(model_data, pseudocount)
#' effector_genes <- snp_link_gene(infima, model_data, fdr = 0.05, cum.pprob = 0.8, cred.set = 0.5)
#' @seealso \code{\link{model_fitting}}, \code{\link{model_input_data}}.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @rawNamespace import(data.table, except = shift)
#' @export
snp_link_gene <-
  function(infima,
           model_data,
           fdr = 0.1,
           cum.pprob = 0.8,
           cred.set = 0.5) {
    V.g <- infima$V.g
    Z.g <- infima$Z.g
    snp_index <- model_data$snp_index
    # consider V.g >= tau as causal
    direct.pprob <- function(p, fdr = 0.1) {
      err <- 1 - p # the error rate
      err <- sort(err) # ascending order
      tau <-
        1 - err[max(which((cumsum(err) / 1:length(err) <= fdr) == 1))]
      return(tau)
    }
    
    tau <- direct.pprob(V.g)
    #### direct pprob approach cutoff
    
    cum.pprob <- 0.85 # cumulative pprob
    cred.set <-
      0.5 # the maximum proportion of candidates in credible sets
    
    # for each local-ATAC-QTL, we record all the genes it is affecting based on INFIMA.
    # We are keeping track of the do.eqtl gene
    
    res <- NULL
    # snps.link.gene <- vector('list', nrow(dt3)) # number of SNPs
    for (x in 1:length(Z.g)) {
      print(x)
      if (V.g[x] >= tau) {
        # remain top k snps as causal
        k <- max(2, floor(length(Z.g[[x]]) * cred.set) + 1)
        k <- min(k, length(Z.g[[x]]))
        # check cumulative pprobs
        if (k >= 2) {
          pprobs <- sort(Z.g[[x]], decreasing = T)
          if (sum(Z.g[[x]]) > cum.pprob) {
            k <- min(k, which(cumsum(pprobs) >= cum.pprob)[1])
          }
        }
        ranks <- rank(Z.g[[x]], ties.method = 'random')
        cutoff <- Z.g[[x]][ranks == length(Z.g[[x]]) - k + 1]
        causal <- snp_index[[x]][Z.g[[x]] >= cutoff]
        for (i in 1:length(causal)) {
          # snps.link.gene[[causal[i]]] <- c(snps.link.gene[[causal[i]]], x)
          snps.link.gene <- dt3[causal[i], snp_id:alt]
          snps.link.gene$ensembl <- dt1$ensembl[x]
          snps.link.gene$symbol <- dt1$symbol[x]
          res <- rbind(res, snps.link.gene)
        }
      }
    }
    return(res)
  }
