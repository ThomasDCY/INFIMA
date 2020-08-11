#' @name raw_data
#' @title raw_data data object
#' @description The organized raw input data.
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
#' @docType data
#' @format A list object.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name model_data
#' @title model_data data object
#' @description The processed input data.
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
#' @docType data
#' @format A list object.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name prior
#' @title prior data object
#' @description The Dirichlet prior information
#'
#' \tabular{ll}{
#' \code{PI} \tab The overall prior term. \cr
#' \code{dist} \tab The distance score. \cr
#' \code{F.g} \tab The footprint annotation. \cr
#' \code{cor.A.B} \tab The Correlation betwwen ATAC-seq signal and founder effect size. \cr
#' \code{cor.A.E} \tab The Correlation between ATAC-seq signal and gene expression. \cr
#' }
#'
#' @docType data
#' @format A list object.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name infima
#' @title infima data object
#' @description An infima object containing the results after model fitting.
#'
#' Parameters (see the paper for details):
#' \tabular{ll}{
#' \code{alpha} \tab The estimated \code{a1}. \cr
#' \code{beta} \tab The estimated \code{b1}. \cr
#' \code{alpha.0} \tab The estimated \code{a0}. \cr
#' \code{beta.0} \tab The estimated \code{b0}. \cr
#' \code{gamma} \tab The estimated \code{gamma}. \cr
#' }
#' The results of INFIMA model:
#' \tabular{ll}{
#' \code{em.iter} \tab The number of iterations used. \cr
#' \code{V.g} \tab The posterior probability for each DO gene to be causal. \cr
#' \code{Z.g} \tab The posterior probability of candidate SNPs for each DO gene. \cr
#' }
#'
#' @docType data
#' @format A list object.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name infima_results
#' @title infima_results data object
#' @description
#' It is a list with length equal to the total number
#' of instances in the DO-eQTL data. \code{infima_results[[i]]} is \code{NULL} if the i-th DO gene does
#' not pass the FDR cutoff. Otherwise, \code{infima_results[[i]]} contains two lists \code{input_data}
#' and \code{output_data}.
#'
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
#' \code{cor.A.B} \tab the Correlation between ATAC-seq signal and founder effect size \cr
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
#'
#' @docType data
#' @format A list object.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name input_query
#' @title input_query data object
#' @description A list of multi-omics input data correspond to the query (snp_id, ensembl, qtl_marker).
#' \tabular{ll}{
#'   \code{Y} \tab DO allelic effects \cr
#'   \code{Y.t} \tab DO allelic effects (trinarized, \code{B6} as reference) \cr
#'   \code{A} \tab local ATAC-seq signal \cr
#'   \code{A.t} \tab local ATAC-seq signal (trinarized, \code{B6} as reference) \cr
#'   \code{B} \tab founder gene expressions from all samples \cr
#'   \code{B.avg} \tab average founder gene expressions from all strains \cr
#'   \code{B.t} \tab average founder gene expressions from all strains (trinarized, \code{B6} as reference) \cr
#'   \code{D} \tab edit distance, the difference between Y.t and E.t \cr
#'   \code{E.t} \tab founder allelic effects (trinarized, \code{B6} as reference) \cr
#'   \code{snpData} \tab the information of the query local-ATAC-QTL \cr
#'   \code{snp_id} \tab the SNP id of the query local-ATAC-QTL \cr
#'   \code{ensembl} \tab the ensembl id of the query DO gene \cr
#'   \code{qtl_marker} \tab the name of the QTL marker of the query DO gene, a DO gene may correspond to multiple markers \cr
#' }
#' @docType data
#' @format A list object.
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name Y
#' @title Y
#' @description raw DO-eQTL signal
#' @docType data
#' @format 8-dimension vector
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name Y.t
#' @title Y.t
#' @description trinarized DO-eQTL signal
#' @docType data
#' @format 8-dimension vector
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name A
#' @title A
#' @description local ATAC-seq signal
#' @docType data
#' @format 8-column matrix, each row correspond to a local-ATAC-QTL
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name A.t
#' @title A.t
#' @description trinarized local ATAC-seq signal
#' @docType data
#' @format 8-column matrix, each row correspond to a local-ATAC-QTL
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name B
#' @title B
#' @description founder gene expressions (all samples)
#' @docType data
#' @format a vector
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name B.avg
#' @title B.avg
#' @description founder gene expressions (averaged over strains)
#' @docType data
#' @format a vector
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name B.t
#' @title B.t
#' @description trinarized founder gene expressions (averaged over strains)
#' @docType data
#' @format a vector
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name do.eqtl
#' @title do.eqtl
#' @description DO-eQTL data: DO gene & marker information
#' @docType data
#' @format a data table
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name snpData
#' @title snpData
#' @description SNP information of causal local-ATAC-QTLs
#' @docType data
#' @format a data table
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name D
#' @title D
#' @description the edit distance, absolute difference between DO and founder allelic effects
#' @docType data
#' @format 8-column matrix, each row correspond to a local-ATAC-QTL
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name E.t
#' @title E.t
#' @description the trinarized founder allelic effects
#' @docType data
#' @format 8-column matrix, each row correspond to a local-ATAC-QTL
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name footprint
#' @title footprint
#' @description the footprint vector
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name footprint.rs
#' @title footprint.rs
#' @description the footprint rank score vector
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL


#' @name dist
#' @title dist
#' @description the distance score component in the prior
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name dist.rs
#' @title dist.rs
#' @description the distance rank score component in the prior
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name cor.A.E
#' @title cor.A.E
#' @description the correlation between ATAC-seq signal and gene expression
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name cor.A.E.rs
#' @title cor.A.E.rs
#' @description the rank score of the correlation between ATAC-seq signal and gene expression
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name cor.A.R
#' @title cor.A.R
#' @description the correlation between ATAC-seq signal and founder effect size
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name cor.A.R.rs
#' @title cor.A.R.rs
#' @description the rank score of the correlation between ATAC-seq signal and founder effect size
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name Z
#' @title Z
#' @description the posterior probabilities of the output
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL



#' @name Z.rs
#' @title Z.rs
#' @description the rank score of the posterior probabilities of the output
#' @docType data
#' @format a vector, the dimension of which equals to the number of local-ATAC-QTLs
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
NULL
