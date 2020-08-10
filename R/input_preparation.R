#' @name raw_input_data
#' @title Load INFIMA inputs from files.
#' @description This function reads in the processed DO-eQTL, ATAC-seq
#' and RNA-seq data as inputs for INFIMA.
#'
#'
#'
#' @param filename1 DO-eQTL data, a data frame object or
#' the path to a comma separated .csv file.
#' Must contain the following columns with exactly the following column names in order:
#' \tabular{ll}{
#' ensembl \tab The ensembl ID of DO-eQTL gene.\cr
#' symbol \tab The symbol of DO-eQTL gene.\cr
#' chr \tab The chromosome of DO-eQTL gene.\cr
#' start \tab The start position of DO-eQTL gene.\cr
#' end \tab The end position of DO-eQTL gene.\cr
#' strand \tab The strand of DO-eQTL gene.\cr
#' qtl_marker \tab The ID of DO-eQTL marker.\cr
#' qtl_chr \tab The chromosome of DO-eQTL marker.\cr
#' qtl_pos \tab The position of DO-eQTL marker.\cr
#' lod \tab The LOD score of DO-eQTL marker.\cr
#' }
#' Plus 8 columns named \code{129, AJ, B6, CAST, NOD, NZO, PWK, WSB}
#' recording the DO-eQTL allelic dependence measures for the 8 strains.
#'
#'
#' @param filename2 Founder RNA-seq data, a matrix, a data frame or
#' the path to a comma separated .csv file.
#' DO gene (row) by sample (column) matrix with normalized gene expressions.
#' The column names must be strain names \code{129, AJ, B6, CAST, NOD, NZO, PWK, WSB}
#' and there must be at least one sample for each strain.
#' One caveat is that the rows of founder RNA-seq data should
#' correspond to the rows of DO-eQTL data (the same DO gene), i.e. the first two input files
#' have the same number of rows.
#'
#'
#' @param filename3 SNP and local ATAC-seq signal data, a data frame object
#' or the path to a comma separated .csv file.
#' Must contain the following columns
#' with exactly the following column names in order:
#' \tabular{ll}{
#' snp_id \tab The SNP ID. \cr
#' chr \tab The chromosome of the SNP. \cr
#' snp_pos \tab The position of the SNP. \cr
#' ref \tab The reference nucleotide. \cr
#' alt \tab The alternative nucleotide. \cr
#' footprint \tab The user-defined footprint annotation.
#' The recommended range is between 0 and 1. The larger the value,
#' the more confidence the SNP is affecting a footprint. Put 0 if
#' the footprint analysis is not applicable. \cr
#' }
#' Plus
#'
#' (a) 8 columns named as \code{strain-genotype} where \code{strain}
#' is \code{129, AJ, B6, CAST, NOD, NZO, PWK, WSB} in order. 0 means reference
#' nucleotide while 2 means alternative nucleotide.
#'
#' (b) 8 columns named as \code{strain-ATAC-seq} where \code{strain}
#' is \code{129, AJ, B6, CAST, NOD, NZO, PWK, WSB} in order. The normalized
#' local ATAC-seq signal of the SNP.
#'
#'
#' @return A raw_data object including:
#' \tabular{ll}{
#' \code{YY} \tab DO-eQTL allelic dependence matrix. \cr
#' \code{AA} \tab Local ATAC-seq signal matrix. \cr
#' \code{EE} \tab Local-ATAC-QTL genotypes. \cr
#' \code{rna.seq} \tab Full founder RNA-seq data matrix. \cr
#' \code{BB} \tab Averaged founder RNA-seq data matrix. \cr
#' \code{FF} \tab Footprint annnotation vector. \cr
#' \code{do.eqtl} \tab DO-eQTL marker and target DO gene information. \cr
#' \code{snpData} \tab SNP information. \cr
#' }
#' @examples
#' data('example-10-genes')
#' raw_data <- raw_input_data(dt1, dt2, dt3)
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @rawNamespace import(data.table, except = shift)
#' @export
raw_input_data <- function(filename1 = NULL,
                           filename2 = NULL,
                           filename3 = NULL,
                           ...) {
  if (is.null(filename1) ||
      is.null(filename2) || is.null(filename3)) {
    stop('Error: input files are missing.')
  }
  
  strains <-
    c('129', 'AJ', 'B6', 'CAST', 'NOD', 'NZO', 'PWK', 'WSB')
  
  if (is.data.frame(filename1)) {
    dt1 <- as.data.table(filename1)
  }
  else{
    dt1 <- fread(filename1)
  }
  refcols1 <- c(
    'ensembl',
    'symbol',
    'chr',
    'start',
    'end',
    'strand',
    'qtl_marker',
    'qtl_chr',
    'qtl_pos',
    'lod',
    strains
  )
  if (any((refcols1 %in% colnames(dt1)) == FALSE)) {
    stop(paste(
      'Error: column(s)',
      paste(refcols1[which(!refcols1 %in% colnames(dt1))], collapse = ', '),
      'are missing in the first input file!'
    ))
  }
  
  if (is.matrix(filename2)) {
    dt2 <- as.data.table(as.numeric(filename2))
  }
  else if (is.data.frame(filename2)) {
    dt2 <- as.data.table(filename2)
  }
  else{
    dt2 <- fread(filename2)
  }
  refcols2 <- strains
  if (any((refcols2 %in% colnames(dt2)) == FALSE)) {
    stop(paste(
      'Error: column(s)',
      paste(refcols2[which(!refcols2 %in% colnames(dt2))], collapse = ', '),
      'are missing in the second input file!'
    ))
  }
  if (any((colnames(dt2) %in% refcols2) == FALSE)) {
    stop(
      paste(
        'Error: the second input file cannot contain column names other than',
        paste(strains, collapse = ', '),
        '!'
      )
    )
  }
  
  if (nrow(dt1) != nrow(dt2)) {
    stop('Error: unequal number of rows in file1 and file2.')
  }
  
  
  if (is.data.frame(filename3)) {
    dt3 <- as.data.table(filename3)
  }
  else{
    dt3 <- fread(filename3)
  }
  refcols3 <- c(
    'snp_id',
    'chr',
    'snp_pos',
    'ref',
    'alt',
    'footprint',
    paste0(strains, '-genotype'),
    paste0(strains, '-ATAC-seq')
  )
  if (any((refcols3 %in% colnames(dt3)) == FALSE)) {
    stop(paste(
      'Error: column(s)',
      paste(refcols3[which(!refcols3 %in% colnames(dt3))], collapse = ', '),
      'are missing in the third input file!'
    ))
  }
  
  DiffChr <- dt1[chr != qtl_chr, .N]
  if (DiffChr > 0) {
    dt2 <- dt2[dt1$chr == dt1$qtl_chr]
    dt1 <- dt1[chr == qtl_chr]
    message(
      paste(
        'Warning:',
        'remove',
        DiffChr,
        'rows that eQTL and gene reside in different chromosomes.'
      )
    )
  }
  
  # prepare the input data
  # YY: the DO-eQTL effect size
  YY <- dt1[, .SD, .SDcols = strains]
  colnames(YY) <- strains
  YY <- as.matrix(YY)
  
  # AA: the ATAC-seq signal
  AA <- dt3[, .SD, .SDcols = paste0(strains, '-ATAC-seq')]
  colnames(AA) <- strains
  AA <- as.matrix(AA)
  
  # EE: the genotype of the SNPs
  EE <- dt3[, .SD, .SDcols = paste0(strains, '-genotype')]
  colnames(EE) <- strains
  EE <- as.matrix(EE)
  
  # rna.seq: the full founder RNA-seq data, to compute effect size later
  rna.seq <- as.matrix(dt2)
  
  # BB: the average founder RNA-seq data
  BB <- matrix(0, nrow = nrow(rna.seq), ncol = 8)
  colnames(BB) <- strains
  for (i in 1:8) {
    BB[, i] <-
      apply(rna.seq[, colnames(rna.seq) == strains[i]], 1, mean)
  }
  
  # FF: the footprint analysis mark for SNPs
  FF <- dt3$footprint
  
  # do.eqtl: the DO-eQTL information
  do.eqtl <- dt1[, ensembl:lod]
  
  # snpData: the SNP information
  snpData <- dt3[, `snp_id`:`WSB-genotype`]
  
  
  raw_data <- list(
    YY = YY,
    AA = AA,
    EE = EE,
    rna.seq = rna.seq,
    BB = BB,
    FF = FF,
    do.eqtl = do.eqtl,
    snpData = snpData
  )
  
  class(raw_data) <- 'raw_data'
  return(raw_data)
}



#' @name data_trinarize
#' @title Data matrix trinarization across 8 founder strains.
#' @description This function trinarizes numeric matrix rowwisely
#' to \{-1,0,1\}. The matrix should have 8 columns corresponding to
#' 8 founder strains \code{129, AJ, B6, CAST, NOD, NZO, PWK, WSB}.
#'
#' First, each row of the matrix is standardized to [0,1].
#'
#' Then, all values subtracted by the \code{B6} value in order to
#' set \code{B6} as the reference.
#'
#' Finally, anything larger than \code{cutoff} are transformed to 1.
#' Anything smaller than -\code{cutoff} are transformed to -1.
#' Anything between are transformed to 0.
#' @param mat The numeric matrix with 8 columns.
#' @param cutoff A value between 0 and 0.5. Default 0.2.
#' @return A numeric matrix after trinarization.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @export
data_trinarize <- function(mat, cutoff = 0.2) {
  if (ncol(mat) != 8) {
    stop('data_trinarize: The number of columns does not equal to 8!')
  }
  
  if (cutoff <= 0 || cutoff >= 0.5) {
    stop('data_trinarize: cutoff parameter out of range (0, 0.5)!')
  }
  
  row.min <- apply(mat, 1, min)
  row.max <- apply(mat, 1, max)
  for (i in 1:8) {
    mat[, i] <- (mat[, i] - row.min) / (row.max - row.min)
  }
  # set B6 to be zero
  ref <- mat[, 3]
  for (i in 1:8) {
    mat[, i] <- mat[, i] - ref
  }
  
  mat[mat < -cutoff] <- -1
  mat[mat > cutoff] <- 1
  mat[mat >= -cutoff & mat <= cutoff] <- 0
  return(mat)
}



#' @name model_input_data
#' @title Process the INFIMA inputs for each DO gene
#' @description This function converts the raw input data into the data
#' for each row of the DO-eQTL data. For each eQTL marker, candidate SNPs
#' are searched from raw inputs within a user defined window size.
#'
#' @param raw_data The output from \code{\link{raw_input_data}}.
#' @param window The window size (Mbp) around DO-eQTL marker contains candidate SNPs.
#' Default = 1 Mbp.
#' @param cutoff Trinarization cutoff. A value between 0 and 0.5.
#' Default 0.2. See \code{\link{data_trinarize}}.
#' @param n_cores The number of cores for parallel computing.
#' Default = \code{detectCores() - 2}.
#' @param verbose Print messages or not. 
#' Default \code{verbose = TRUE}.
#' @return A model_data object: the processed data objects are broken down
#' into small elements for each row of the DO-eQTL data.
#'
#' \tabular{ll}{
#' \code{E.g} \tab Trinarized effect size of candidate SNPs with respect to DO gene. \cr
#' \code{Y.g} \tab Trinarized DO-eQTL allelic dependence. \cr
#' \code{A.g} \tab Trinarized local ATAC-seq signal for candidate SNPs. \cr
#' \code{F.g} \tab Footprint annotations for candidate SNPs. \cr
#' \code{D.g} \tab Edited distances between candidate SNP and DO-eQTL signal. \cr
#' \code{B.g} \tab Trinarized gene expression in founder RNA-seq data. \cr
#' \code{snp_index} \tab The indices of candidate SNPs in the original list of SNPs. \cr
#' \code{p.g} \tab The number of candidate SNPs for each DO gene g. \cr
#' \code{window} \tab The window size in Mb. \cr
#' }
#' @examples
#' data('example-10-genes')
#' raw_data <- raw_input_data(dt1, dt2, dt3)
#' model_data <- model_input_data(raw_data)
#' @seealso \code{\link{raw_input_data}}, \code{\link{data_trinarize}}.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @rawNamespace import(data.table, except = shift)
#' @import GenomicRanges
#' @import parallel
#' @import doParallel
#' @import foreach
#' @rawNamespace importFrom(IRanges,IRanges)
#' @rawNamespace importFrom(S4Vectors,queryHits)
#' @rawNamespace importFrom(S4Vectors,subjectHits)
#' @importClassesFrom GenomicRanges GRanges
#' @export
model_input_data <- function(raw_data = NULL,
                             window = 1,
                             cutoff = 0.2,
                             n_cores = detectCores() - 2,
                             verbose = TRUE,
                             ...) {
  if (is.null(raw_data)) {
    stop('Error: please input the raw_data.')
  }
  
  if (class(raw_data) != 'raw_data') {
    stop('Error: please input the S3 class raw_data.')
  }
  
  
  if (window <= 0) {
    stop('Error: the window around eQTL marker must be positive!')
  }
  
  time.start <- proc.time()
  
  window <- window * 1000000
  
  YY <- raw_data$YY
  AA <- raw_data$AA
  EE <- raw_data$EE
  rna.seq <- raw_data$rna.seq
  BB <- raw_data$BB
  FF <- raw_data$FF
  do.eqtl <- raw_data$do.eqtl
  snpData <- raw_data$snpData
  
  stopifnot(exprs = {
    ! is.null(YY)
    ! is.null(AA)
    ! is.null(EE)
    ! is.null(rna.seq)
    ! is.null(BB)
    ! is.null(FF)
    ! is.null(do.eqtl)
    ! is.null(snpData)
  })
  
  #### filter out those do.eqtl rows that do not have overlapping SNPs ####
  G <- nrow(do.eqtl)
  
  YY.t <- data_trinarize(YY, cutoff)
  AA.t <- data_trinarize(AA, cutoff)
  BB.t <- data_trinarize(BB, cutoff)
  
  
  cnt <- as.vector(table(colnames(rna.seq)))
  genotype <- apply(EE, 1, function(x) {
    rep(x, cnt)
  })
  
  snps <- GRanges(
    seqnames = snpData$chr,
    ranges = IRanges(start = snpData$snp_pos, end = snpData$snp_pos),
    strand = '*'
  )
  
  registerDoParallel(cores = n_cores)
  
  if(verbose){
    message(paste('n_cores =', n_cores))
  }
  
  r <- foreach(g = 1:G) %dopar% {
    # if (verbose) {
    #   message(paste('Current row in DO-eQTL data:', g))
    # }
    Y.g <- YY.t[g, ]
    B.g <- BB.t[g, ]
    # find associated SNPs
    snp.loc <- GRanges(
      seqnames = do.eqtl$qtl_chr[g],
      ranges = IRanges(
        start = do.eqtl$qtl_pos[g] - window,
        end = do.eqtl$qtl_pos[g] + window
      ),
      strand = '*'
    )
    index <- subjectHits(findOverlaps(snp.loc, snps))
    snp_index <- index
    p.g <- length(snp_index)
    
    # if no local-ATAC-QTL, then return NULL
    E.g <- NULL
    A.g <- NULL
    F.g <- NULL
    D.g <- NULL
    
    
    if (!identical(index, integer(0))) {
      # there exist at least one overlap
      A.g <- AA.t[index, ]
      F.g <- FF[index]
      E.g <- EE[index, ] # genotype.g
      
      # get effect size E.g
      exp.g <- rna.seq[g, ] # gene expression
      atac.qtls.g <- genotype[, index] # regressors
      
      if (is.null(dim(E.g))) {
        # only one SNP found, correct the dimension
        E.g <- t(as.matrix(E.g))
        A.g <- t(as.matrix(A.g))
        atac.qtls.g <- t(t(atac.qtls.g))
      }
      
      
      # run eQTL and trinarize E.g
      for (i in 1:length(index)) {
        lm <- lm(exp.g ~ atac.qtls.g[, i])
        estimate <- summary(lm)$coefficients[2, 1]
        pval <- summary(lm)$coefficients[2, 4]
        M <- 0
        if (pval < 0.05) {
          M <- ifelse(estimate > 0, 1,-1)
        }
        E.g[i, ] <- E.g[i, ] / 2 * M
      }
      
      # compute the absolute distance
      D.g <- matrix(0, nrow(E.g), ncol(E.g))
      for (i in 1:nrow(E.g)) {
        D.g[i, ] <- abs(E.g[i, ] - Y.g)
      }
      
    }
    
    return(list(
      E.g = E.g,
      Y.g = Y.g,
      A.g = A.g,
      F.g = F.g,
      D.g = D.g,
      B.g = B.g,
      snp_index = snp_index,
      p.g = p.g
    ))
  }
  
  
  E.g <- lapply(r, function(x) x$E.g)
  Y.g <- lapply(r, function(x) x$Y.g)
  A.g <- lapply(r, function(x) x$A.g)
  F.g <- lapply(r, function(x) x$F.g)
  D.g <- lapply(r, function(x) x$D.g)
  B.g <- lapply(r, function(x) x$B.g)
  p.g <- sapply(r, function(x) x$p.g)
  snp_index <- lapply(r, function(x) x$snp_index)
  
  if(verbose){
    cat("Time taken", proc.time()[3] - time.start[3])
  }
  
  model_data <- list(
    E.g = E.g,
    Y.g = Y.g,
    A.g = A.g,
    F.g = F.g,
    D.g = D.g,
    B.g = B.g,
    snp_index = snp_index,
    p.g = p.g,
    window = window
  )
  class(model_data) <- 'model_data'
  return(model_data)
}







