#' @name snp_link_gene
#' @title Link SNPs to effector genes
#'
#' @description Analyze the outputs from INFIMA and link local-ATAC-QTLs to effector genes.
#' Including the posterior probabilities, target genes, each piece of input information,
#' and the relative ranks of each piece of information.
#'
#' @param infima The output from \code{model_fitting}
#' @param raw_data The output from \code{\link{raw_input_data}}.
#' @param model_data The output from \code{\link{model_input_data}}.
#' @param fdr The False Discovery Rate cutoff with direct posterior probability
#' approach (Newton, 2004). Default 0.1.
#' @param cum.pprob The maximum cumulative posterior probability for
#' the target credible set. Default 0.8. Take values between 0 and 1.
#' @param cred.set The maximum proportion of candidate SNPs allowed
#' in the credible set. Default 0.5. Take values between 0 and 1.
#' @param n_cores The number of cores for parallel computing.
#' @param verbose Print messages or not.
#' Default = \code{detectCores() - 2}.
#'
#'
#' @return An \code{infima_results} object. It is a list with length equal to the total number
#' of instances in the DO-eQTL data. \code{infima_results[[i]]} is \code{NULL} if the i-th DO gene does
#' not pass the FDR cutoff. Otherwise, \code{infima_results[[i]]} contains two lists \code{input_data}
#' and \code{output_data}.
#'
#' \code{input_data} contains:
#' \tabular{ll}{
#' \code{Y} \tab raw DO-eQTL signal \cr
#' \code{Y.t} \tab trinarized DO-eQTL signal \cr
#' \code{A} \tab raw ATAC-seq signal \cr
#' \code{A.t} \tab trinarized ATAC-seq signal \cr
#' \code{B} \tab raw RNA-seq data (all samples) \cr
#' \code{B.avg} \tab raw RNA-seq data (averaged) \cr
#' \code{B.t} \tab trinarized RNA-seq data \cr
#' \code{do.eqtl} \tab DO-eQTL data: DO gene & marker information \cr
#' \code{snpData} \tab SNP information of causal local-ATAC-QTLs \cr
#' \code{D} \tab the edit distance \cr
#' \code{footprint} \tab the footprint information \cr
#' \code{E.t} \tab the trinarized allelic effects \cr
#' \code{dist} \tab the distance score component in the prior \cr
#' \code{cor.A.E} \tab the Correlation between ATAC-seq signal and gene expression \cr
#' \code{cor.A.B} \tab the Correlation betwwen ATAC-seq signal and founder effect size \cr
#' \code{p} \tab the total number of candidates \cr
#' \code{k} \tab the number of candidates in the credible set \cr
#' }
#'
#' \code{output_data} contains:
#' \tabular{ll}{
#' \code{footprint.rs} \tab the footprint information rank score \cr
#' \code{dist.rs} \tab the distance rank score \cr
#' \code{cor.A.E.rs} \tab the cor.A.E rank score \cr
#' \code{cor.A.B.rs} \tab the cor.A.B rank score \cr
#' \code{Z} \tab the posterior probabilities of the output \cr
#' \code{Z.rs} \tab the rank score of the posterior probabilities \cr
#' }
#'
#' @examples
#' data('example-10-genes')
#' raw_data <- raw_input_data(dt1, dt2, dt3)
#' model_data <- model_input_data(raw_data)
#' prior <- compute_prior(raw_data, model_data)
#' infima <- model_fitting(model_data, prior)
#' infima_results <- snp_link_gene(infima, raw_data, model_data, prior, fdr = 0.05, cum.pprob = 0.8, cred.set = 0.5)
#'
#'
#' @seealso \code{\link{model_fitting}}, \code{\link{model_input_data}}.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import parallel
#' @import doParallel
#' @import foreach
#' @rawNamespace import(data.table, except = shift)
#' @export
snp_link_gene <-
  function(infima,
           raw_data,
           model_data,
           prior,
           fdr = 0.1,
           cum.pprob = 0.8,
           cred.set = 0.5,
           n_cores = detectCores() - 2,
           verbose = TRUE) {
    if (is.null(infima)) {
      stop('Error: please input the infima.')
    }
    
    if (class(infima) != 'infima') {
      stop('Error: please input the S3 class infima.')
    }
    
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
    
    if (is.null(prior)) {
      stop('Error: please input the prior.')
    }
    
    if (class(prior) != 'prior') {
      stop('Error: please input the S3 class prior.')
    }
    
    stopifnot(exprs = {
      fdr > 0 && fdr < 1
      cum.pprob > 0 && cum.pprob < 1
      cred.set > 0 && cred.set < 1
    })
    
    time.start <- proc.time()
    
    ###### data objects to be extracted ######
    V.g <- infima$V.g
    Z.g <- infima$Z.g
    
    YY <- raw_data$YY # DO-eQTL signal
    AA <- raw_data$AA # ATAC-seq signal
    EE <- raw_data$EE # local-ATAC-QTL genotype
    rna.seq <- raw_data$rna.seq # RNA-seq data (all samples)
    BB <- raw_data$BB # RNA-seq data (averaged)
    FF <- raw_data$FF # footprint information
    do.eqtl <- raw_data$do.eqtl # DO-eQTL data
    snpData <- raw_data$snpData # the SNP information
    
    E.g <- model_data$E.g # effect size, trinarized
    Y.g <- model_data$Y.g # trinarized DO-eQTL signal
    A.g <- model_data$A.g # trinarized ATAC-seq signal
    F.g <- model_data$F.g # footprint information
    D.g <- model_data$D.g # edit distance
    B.g <-
      model_data$B.g # trinarized founder RNA-seq gene expression
    snp_index <-
      model_data$snp_index # the index of local-ATAC-QTLs in the original list
    p.g <- model_data$p.g # the total number of candidates
    
    dist <- prior$dist # dist score
    cor.A.E <-
      prior$cor.A.E # correlation between ATAC-seq signal and founder allelic effects
    cor.A.B <-
      prior$cor.A.B # correlation between ATAC-seq signal and founder gene expression
    
    
    # consider V.g >= tau as causal
    tau <- .direct_pprob(V.g, fdr)
    # direct pprob approach cutoff
    
    
    # for each local-ATAC-QTL, we record all the genes it is affecting based on INFIMA.
    # We are keeping track of the do.eqtl gene
    
    registerDoParallel(cores = n_cores)
    
    if (verbose) {
      message(paste('n_cores =', n_cores))
    }
    
    res <- NULL
    
    r <- foreach(x = 1:length(Z.g)) %dopar% {
      # print(x)
      if (V.g[x] >= tau) {
        # treat top k snps as causal under constraints
        
        # credible set constraint
        k <- max(2, floor(length(Z.g[[x]]) * cred.set) + 1)
        k <- min(k, length(Z.g[[x]]))
        
        # cumulative pprobs constraint
        pprobs <- sort(Z.g[[x]], decreasing = T)
        if (k >= 2) {
          if (sum(Z.g[[x]]) > cum.pprob) {
            k <- min(k, which(cumsum(pprobs) >= cum.pprob)[1])
          }
        }
        
        cutoff <- pprobs[k]
        # indicators of causal local-ATAC-QTLs among candidates
        subset <- Z.g[[x]] >= cutoff
        
        # the SNP indices in the original list
        causal <- snp_index[[x]][Z.g[[x]] >= cutoff]
        
        ###### organize each piece of input data as a list ######
        input_data <- list(
          Y = YY[x,],
          # raw DO-eQTL signal
          Y.t = Y.g[[x]],
          # trinarized DO-eQTL signal
          A = AA[causal,],
          # raw ATAC-seq signal
          A.t = A.g[[x]][subset,],
          # trinarized ATAC-seq signal
          B = rna.seq[x,],
          # raw RNA-seq data (all samples)
          B.avg = BB[x,],
          # raw RNA-seq data (averaged)
          B.t = B.g[[x]],
          # trinarized RNA-seq data
          do.eqtl = do.eqtl[x],
          # DO-eQTL data: DO gene & marker information
          snpData = snpData[causal],
          # SNP information of causal local-ATAC-QTLs
          D = D.g[[x]][subset,],
          # the edit distance
          footprint = FF[causal],
          # the footprint information
          E.t = E.g[[x]][subset,],
          # the trinarized allelic effects
          dist = dist[[x]][subset],
          # the distance score component in the prior
          cor.A.E = cor.A.E[[x]][subset],
          # the cor.A.E component
          cor.A.B = cor.A.B[[x]][subset],
          # the cor.A.B component
          p = p.g[[x]],
          # the total number of candidates
          k = k # the number of candidates in the credible set
        )
        
        ###### organize the output data ######
        output_data <- list(
          footprint.rs = rank(F.g[[x]])[subset] / p.g[[x]],
          # the footprint information rank score
          dist.rs = rank(dist[[x]])[subset] / p.g[[x]],
          # the distance rank score
          cor.A.E.rs = rank(cor.A.E[[x]])[subset] / p.g[[x]],
          # the cor.A.E rank score
          cor.A.B.rs = rank(cor.A.B[[x]])[subset] / p.g[[x]],
          # the cor.A.B rank score
          Z = Z.g[[x]][subset],
          # the posterior probabilities of the output
          Z.rs = rank(Z.g[[x]])[subset] / p.g[[x]] # the rank score of the posterior probabilities
        )
        
        
        return(list(input_data = input_data,
                    output_data = output_data))
      }
    }
    
    infima$tau <- tau
    infima_results <- list(infima = infima,
                           results = r)
    class(infima_results) <- 'infima_results'
    
    if (verbose) {
      cat("Time taken", proc.time()[3] - time.start[3])
    }
    
    return(infima_results)
  }



# Direct posterior probability approach
# Micheal A. Newton, 2014
.direct_pprob <- function(p, fdr) {
  err <- 1 - p # the error rate
  err <- sort(err) # ascending order
  tau <-
    1 - err[max(which((cumsum(err) / 1:length(err) <= fdr) == 1))]
  return(tau)
}



#' @name query_input_data
#' @title Given query snp_id and ensembl id, return the multi-omics input data.
#' @description The S3 method for the class \code{infima_results}
#'
#' @param x The S3 class \code{infima_results}
#' @param snp_id The query snp_id of local-ATAC-QTL
#' @param ensembl The query ensembl column value of DO-eQTL data
#' @param qtl_marker The query qtl_marker column value of DO-eQTL data
#' @return an input_query class, which is a list of input data corresponding to the query
#' 
#' @examples 
#' input_query <- query_input_data(infima_results, snp_id = 'rs51076312', ensembl = 'ENSMUSG00000037995', qtl_marker = '1_172713578')
#' 
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' 
#' @rawNamespace import(data.table, except = shift)
#' @export
query_input_data <- function(x,
                             snp_id = NULL,
                             ensembl = NULL,
                             qtl_marker = NULL,
                             ...) {
  stopifnot(class(x) == 'infima_results')
  
  stopifnot(exprs = {
    !is.null(snp_id)
    ! is.null(ensembl)
    ! is.null(qtl_marker)
  })
  
  infima <- x$infima
  results <- x$results
  results <- results[sapply(results, function(x)
    ! is.null(x))]
  
  # find the indices
  ensembls <-
    sapply(results, function(x)
      x$input_data$do.eqtl$ensembl)
  qtl_markers <-
    sapply(results, function(x)
      x$input_data$do.eqtl$qtl_marker)
  
  if (!ensembl %in% ensembls) {
    message('ensembl id not found!')
    return(NULL)
  }
  
  if (!qtl_marker %in% qtl_markers) {
    message('qtl_marker id not found!')
    return(NULL)
  }
  
  if (!(ensembl %in% ensembls && qtl_marker %in% qtl_markers)) {
    message('ensembl and qtl marker combination not found!')
    return(NULL)
  }
  
  # This is uniquely defined
  current <-
    results[ensembls == ensembl & qtl_markers == qtl_marker][[1]]
  input_data <- current$input_data
  
  snps <- input_data$snpData$snp_id
  
  if (!snp_id %in% snps) {
    message('snp id not found!')
    return(NULL)
  }
  
  ind <- snps == snp_id
  
  res <- list(
    Y = input_data$Y,
    Y.t = input_data$Y.t,
    A = input_data$A[ind, ],
    A.t = input_data$A.t[ind, ],
    B = input_data$B,
    B.avg = input_data$B.avg,
    B.t = input_data$B.t,
    D = input_data$D[ind, ],
    E.t = input_data$E.t[ind, ],
    do.eqtl = input_data$do.eqtl,
    snpData = input_data$snpData[ind, ],
    snp_id = snp_id,
    ensembl = ensembl,
    qtl_marker = qtl_marker
  )
  
  class(res) <- 'input_query'
  return(res)
}



summary <- function(x) {
  UseMethod('summary', x)
}

#' @name summary
#' @title summary of infima_results
#' @method summary infima_results
#'
#'
#' @description S3 method for the class \code{infima_results}
#'
#' @param x The \code{infima_results} S3 class
#' @param ... Additional arguments
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @export
summary.infima_results <- function(x, ...) {
  stopifnot(class(x) == 'infima_results')
  
  infima <- x$infima
  results <- x$results
  G <- length(infima$V.g)
  cat(paste('Total # of DO-eQTL genes:', G, '\n'))
  cat(paste(
    'Total # of DO-eQTL genes containing causal SNP:',
    sum(infima$V.g >= infima$tau),
    '\n'
  ))
  # cat('The contribution of each component of prior')
}


#' @name summary
#' @title summary of infima
#' @method summary infima
#'
#' @description S3 method for the class \code{infima}
#'
#' @param x The \code{infima} S3 class
#' @param ... Additional arguments
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @export
summary.infima <- function(x, ...) {
  stopifnot(class(x) == 'infima')
  
  cat(paste('INFIMA converges after', x$em.iter, 'iterations!\n'))
  cat('Null generative model parameters:\n')
  cat(paste0('a0 = (', paste(
    format(infima$alpha.0, digits = 4), collapse = ', '
  ), ')\n'))
  cat(paste0('b0 = (', paste(
    format(infima$beta.0, digits = 4), collapse = ', '
  ), ')\n'))
  cat('Causal generative model parameters:\n')
  cat(paste0('a1 = (', paste(
    format(infima$alpha, digits = 4), collapse = ', '
  ), ')\n'))
  cat(paste0('b1 = (', paste(
    format(infima$beta, digits = 4), collapse = ', '
  ), ')\n'))
  cat(paste('gamma =', format(infima$gamma, digits = 4)))
  
  summary <- list(
    a0 = infima$alpha.0,
    b0 = infima$beta.0,
    a1 = infima$alpha,
    b1 = infima$beta,
    gamma = infima$gamma,
    n.iter = infima$em.iter
  )
  invisible(summary)
}



#' @name plot
#' @title plot of input_query
#' @method plot input_query
#'
#' @description S3 method for the class \code{input_query}
#'
#' @param x The \code{input_query} S3 class
#' @param option Which piece of data to plot
#' 1: DO allele effect data
#' 2: ATAC-seq data
#' 3: founder gene expression data
#' 4: edit distance and founder allele efect
#' @param ... Additional arguments
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import ggplot2
#' @export
plot.input_query <- function(x, option = NULL){
  stopifnot(class(x) == 'input_query')
  
  if(is.null(option)){
    stop('Please input the option argument!')
  }
  
  if(! option %in% c(1,2,3,4)){
    stop('Invalid value for the option argument!')
  }
  
  colors <- c(rgb(240, 128, 128, maxColorValue = 255, alpha = 255),
              rgb(218, 165, 32, maxColorValue = 255, alpha = 255),
              rgb(128, 128, 128, maxColorValue = 255, alpha = 255),
              rgb(0, 160, 0, maxColorValue = 255, alpha = 255),
              rgb(16, 16, 240, maxColorValue = 255, alpha = 255),
              rgb(0, 160, 240, maxColorValue = 255, alpha = 255),
              rgb(240, 0, 0, maxColorValue = 255, alpha = 255),
              rgb(144, 0, 224, maxColorValue = 255, alpha = 255))

  if(option == 1){
    grid <- par(mfrow = c(1,2), 
                mar = c(1,1,1,1) + 3)
    plot(x$Y, main = '', xaxt = 'n', 
         xlab = 'Strain', ylab = 'DO allele effect',
         col = colors, pch = 19)
    axis(1, at=1:8, labels = names(x$Y))
    
    plot(x$Y.t, main = '', xaxt = 'n', 
         xlab = 'Strain', ylab = 'DO allele effect (trinarized)',
         col = colors, pch = 19)
    axis(1, at=1:8, labels = names(x$Y.t))
    par(grid)
  }
  
  if(option == 2){
    plot(x$A, main = '', xaxt = 'n', 
         xlab = 'Strain', ylab = 'Local ATAC-seq signal',
         col = colors, pch = 19)
    axis(1, at=1:8, labels = names(x$A))
    
    plot(x$A.t, main = '', xaxt = 'n', 
         xlab = 'Strain', ylab = 'Local ATAC-seq signal (trinarized)',
         col = colors, pch = 19)
    axis(1, at=1:8, labels = names(x$A.t))
    par(grid)
  }
  
  if(option == 3){
    plot(x$B.avg, main = '', xaxt = 'n', 
         xlab = 'Strain', ylab = 'Founder gene expression',
         col = colors, pch = 19)
    axis(1, at=1:8, labels = names(x$B.avg))
    
    plot(x$B.t, main = '', xaxt = 'n', 
         xlab = 'Strain', ylab = 'Founder gene expression (trinarized)',
         col = colors, pch = 19)
    axis(1, at=1:8, labels = names(x$B.t))
    par(grid)
  }
  
  if(option == 4){
    plot(x$E.t, main = '', xaxt = 'n', 
         xlab = 'Strain', ylab = 'Founder allele effect (trinarized)',
         col = colors, pch = 19)
    axis(1, at=1:8, labels = names(x$E.t))
    
    plot(x$D, main = '', xaxt = 'n', 
         xlab = 'Strain', ylab = 'Edit distance',
         col = colors, pch = 19)
    axis(1, at=1:8, labels = names(x$E.t))
    
    par(grid)
  }
}



as.data.frame <- function(x) {
  UseMethod('as.data.frame', x)
}


#' @name as.data.frame.infima_results
#' @title as.data.frame for infima_results
#' @description The S3 method for the class \code{infima_results}
#'
#' @param x The S3 class \code{infima_results}
#' @param n_cores The number of cores
#' @param ... Additional arguments
#' @return a data frame containing the following colums:
#' \tabular{ll}{
#' \code{snp_id ... alt} \tab SNP information of the causal local-ATAC-QTLs \cr
#' \code{ensembl ... lod} \tab DO-eQTL data: DO gene & marker information \cr
#' \code{footprint} \tab the footprint information \cr
#' \code{footprint.rs} \tab the footprint information rank score \cr
#' \code{dist} \tab the distance score component in the prior \cr
#' \code{dist.rs} \tab the distance rank score \cr
#' \code{cor.A.E} \tab the Correlation between ATAC-seq signal and gene expression \cr
#' \code{cor.A.E.rs} \tab the cor.A.E rank score \cr
#' \code{cor.A.B} \tab the Correlation between ATAC-seq signal and founder effect size \cr
#' \code{cor.A.B.rs} \tab the cor.A.B rank score \cr
#' \code{p} \tab the total number of candidates \cr
#' \code{k} \tab the number of candidates in the credible set \cr
#' \code{Z} \tab the posterior probabilities of the output \cr
#' \code{Z.rs} \tab the rank score of the posterior probabilities \cr
#' }
#' 
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' 
#' @rawNamespace import(data.table, except = shift)
#' @import GenomicRanges
#' @import parallel
#' @import doParallel
#' @import foreach
#' @export
as.data.frame.infima_results <-
  function(x, n_cores = detectCores() - 2, ...) {
    stopifnot(class(x) == 'infima_results')
    
    registerDoParallel(cores = n_cores)
    
    infima <- x$infima
    results <- x$results
    results <- results[sapply(results, function(x)
      ! is.null(x))]
    
    res <- foreach(i = 1:length(results)) %dopar% {
      current <- results[[i]]
      input_data <- current$input_data
      output_data <- current$output_data
      
      
      dt1 <- input_data$snpData[, snp_id:alt]
      dt2 <- input_data$do.eqtl
      colnames(dt2)[colnames(dt2) == 'chr'] <- 'gene_chr'
      dt3 <- data.table(
        cor.A.E = input_data$cor.A.E,
        cor.A.E.rs = output_data$cor.A.E.rs,
        cor.A.B = input_data$cor.A.B,
        cor.A.B.rs = output_data$cor.A.B.rs,
        footprint = input_data$footprint,
        footprint.rs = output_data$footprint.rs,
        dist = input_data$dist,
        dist.rs = output_data$dist.rs,
        p = input_data$p,
        k = input_data$k,
        Z = output_data$Z,
        Z.rs = output_data$Z.rs
      )
      return(cbind(dt1, dt2, dt3))
    }
    
    # the combined data table for INFIMA fine-mapping results
    return(do.call('rbind', res))
  }



