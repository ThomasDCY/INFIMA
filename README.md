## INFIMA
**INFIMA** is an R package for the **In**tegrative **Fi**ne-**Ma**pping with Model Organism Multi-Omics Data. **INFIMA** utilizes the diversity outbred (DO) mice population as the model organsim. The major usage of the **INFIMA** package is to fine-map the eQTL markers in DO studies (DO-eQTL). **INFIMA** implements an empirical Bayes model which quantifies how well each non-coding SNP explains the observed DO allelic patern through consistency of founder mice RNA-seq data, founder mice ATAC-seq data (including the existence and consistency of a footprint) with the observed allelic pattern. **INFIMA** assumes 0 or 1 causal SNPs for a DO-eQTL gene. Given a library of mouse SNPs, the corresponding local ATAC-seq signals, DO-eQTL data as well as the gene expression in founder mice, **INFIMA** provides the following functionalities:

* Fine-mapping DO-eQTLs and estimate SNP-level posterior probabilities.
* Linking SNPs or the local ATAC-seq peaks to effector genes.

**INFIMA** results from DO studies can be mapped to human, which provides putative effector genes of human GWAS SNPs.

## Installation

**INFIMA** will be submitted to Bioconductor. Currently, **INFIMA** can be downloaded and installed in R by: 

```r
# install.packages("devtools")
devtools::install_github("ThomasDCY/INFIMA")
```

**INFIMA** depends on the following R packages:

* `data.table` is used for formatting results that are easy for users to query.
* `GenomicRanges` is used for operating genomic intervals.
* `doParallel` is used for parallel computing.

## A quick start

```r
library(INFIMA)

## Use the example data
data('example-10-genes')

## dt1: DO mouse eQTL results
## dt2: founder RNA-seq count matrix
## dt3: SNP information (with footprint annotation) and the local ATAC-seq signal
## For more information about the input data format, please refer to the vignette or ?raw_input_data
raw_data <- raw_input_data(dt1, dt2, dt3)
model_data <- model_input_data(raw_data)

## Dirichlet prior for the INFIMA model
prior <- compute_prior(raw_data, model_data)

## Model fitting
infima <- model_fitting(model_data, prior)

## After fitting the INFIMA model, we focus on the DO genes containing a causal SNP
infima_results <- snp_link_gene(infima, raw_data, model_data, prior, fdr = 0.05, cum.pprob = 0.8, cred.set = 0.5)

## Obtain the fine-mapping results
## Link mouse SNPs to target genes
results <- as.data.frame(infima_results)

## An example input data query
input_query <- query_input_data(infima_results, snp_id = 'rs51076312', ensembl = 'ENSMUSG00000037995', qtl_marker = '1_172713578')

## Plot the query input data
plot_input(input_query, option = 1) ## DO allele effect data
plot_input(input_query, option = 2) ## ATAC-seq data
plot_input(input_query, option = 3) ## RNA-seq data
plot_input(input_query, option = 4) ## edit distance and founder allele effect


```

See the vignette for more information!

### Reference

**C. Dong**, S. Simonett, S. Shin, D. Stapleton, K. Schueler, G. Churchill, L. Lu, X. Liu, F. Jin, Y. Li, A. Attie, M. Keller, and **S. Keleş**, "INFIMA leverages multi-omics model organism data to identify effector genes of human GWAS variants", *Manuscript*.


