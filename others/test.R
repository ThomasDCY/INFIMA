rm(list = ls())
library(INFIMA)
# filename1 <- '~/Documents/GitHub/INFIMA/data/input1-chr1.csv'
# filename2 <- '~/Documents/GitHub/INFIMA/data/input2-chr1.csv'
# filename3 <- '~/Documents/GitHub/INFIMA/data/input3-chr1.csv'
# dt1 <- fread(filename1)
# dt2 <- fread(filename2)
# dt3 <- fread(filename3)
# dt1 <- dt1[1:10]
# dt2 <- dt2[1:10]
# save(dt1, dt2, dt3, file = '~/Documents/GitHub/INFIMA/data/example-10-genes.rda')

data("example-10-genes")

raw_data <- raw_input_data(dt1, dt2, dt3)
model_data <- model_input_data(raw_data)
pseudocount <- compute_pseudocount(raw_data, model_data)
infima <- model_fitting(model_data, pseudocount, verbose = T)
effector_genes <- snp_link_gene(infima, model_data, fdr = 0.05, cum.pprob = 0.8, cred.set = 0.5)

context("Expected output")
test_that("snp_link_gene",{
  expect_equal(nrow(effector_genes), 136)
})

context("S3 Classes")
test_that("S3 Classes",{
  expect_equal(class(infima), "infima")
  expect_equal(class(pseudocount), "pseudocount")
  expect_equal(class(model_data), "model_data")
  expect_equal(class(raw_data), "raw_data")
})