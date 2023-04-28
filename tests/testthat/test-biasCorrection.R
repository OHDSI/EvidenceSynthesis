library(testthat)
library(EvidenceSynthesis)

# load the example likelihoods
data("ncLikelihoods")
data("ooiLikelihoods")

test_that('fit bias distribution', {
  # normal model
  biasDist5 = fitBiasDistribution(ncLikelihoods[[5]], numsamps = 1000)

  expect_type(biasDist5, "list")
  expect_named(biasDist5, c('mean', 'scale', 'bias'))
  expect_type(biasDist5$bias, 'double')

  # t model
  biasDist5 = fitBiasDistribution(ncLikelihoods[[5]], numsamps = 1000, robust = TRUE)

  expect_type(biasDist5, "list")
  expect_named(biasDist5, c('mean', 'scale', 'bias'))
  expect_type(biasDist5$bias, 'double')
})

test_that('sequentially fit bias distribution', {
  biasDist = sequentialFitBiasDistribution(ncLikelihoods, numsamps = 1000)

  expect_type(biasDist, "list")
  expect_named(biasDist, c('mean', 'scale', 'bias', 'Id'))
  expect_type(biasDist$bias, 'double')
})

test_that('Bayesian bias correction', {
  bbcSequential = biasCorrectionInference(ooiLikelihoods,
                                          ncLikelihoods,
                                          numsamps = 1000,
                                          doCorrection = TRUE)

  expect_type(bbcSequential, "list")
  expect_named(bbcSequential, c("median", "mean", "ci95Lb", "ci95Ub", "p1", "Id"))
  expect_gte(bbcSequential$ci95Ub[1], bbcSequential$median[1])
  expect_lte(bbcSequential$ci95Lb[1], bbcSequential$median[1])
})

test_that('plot bias distribution', {
  biasDist = sequentialFitBiasDistribution(ncLikelihoods, numsamps = 5000)
  biasPlot = plotBiasDistribution(biasDist)

  expect_s3_class(biasPlot, "ggplot")
})

test_that("plot bias correction", {
  bbcSequential = biasCorrectionInference(ooiLikelihoods,
                                          ncLikelihoods,
                                          numsamps = 5000,
                                          doCorrection = TRUE)
  bbcPlot = plotBiasCorrectionInference(bbcSequential)

  expect_s3_class(bbcPlot, "grob")
})
