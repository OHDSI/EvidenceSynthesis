library(testthat)

createApproximations <- function(populations, approximation) {
  fitModelInDatabase <- function(population, approximation) {
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
    approximation <-  approximateLikelihood(cyclopsFit, "x", approximation = approximation)
    return(approximation)
  }
  data <- lapply(populations, fitModelInDatabase, approximation = approximation)
  data <- do.call("rbind", data)
  return(data)
}

set.seed(1)
populations <- simulatePopulations(settings = createSimulationSettings(n = 2500, treatedFraction = 0.5, hazardRatio = 2))

test_that("Normal approximation", {
  data <- createApproximations(populations, "normal")
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_lt(estimate$lb, 2)
  expect_gt(estimate$ub, 2)
  
  estimate <- computeBayesianMetaAnalysis(data, chainLength = 1e5, burnIn = 1e4)
  expect_lt(estimate$mu95Lb, 2)
  expect_gt(estimate$mu95Ub, 2)
})

test_that("Skew normal approximation", {
  data <- createApproximations(populations, "skew normal")
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_lt(estimate$lb, 2)
  expect_gt(estimate$ub, 2)
  
  estimate <- computeBayesianMetaAnalysis(data, chainLength = 1e5, burnIn = 1e4)
  expect_lt(estimate$mu95Lb, 2)
  expect_gt(estimate$mu95Ub, 2)
})

test_that("Custom approximation", {
  data <- createApproximations(populations, "custom")
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_lt(estimate$lb, 2)
  expect_gt(estimate$ub, 2)
  
  estimate <- computeBayesianMetaAnalysis(data, chainLength = 1e5, burnIn = 1e4)
  expect_lt(estimate$mu95Lb, 2)
  expect_gt(estimate$mu95Ub, 2)
})

test_that("Grid approximation", {
  data <- createApproximations(populations, "grid")
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_lt(estimate$lb, 2)
  expect_gt(estimate$ub, 2)
  
  estimate <- computeBayesianMetaAnalysis(data, chainLength = 1e5, burnIn = 1e4)
  expect_lt(estimate$mu95Lb, 2)
  expect_gt(estimate$mu95Ub, 2)
})

test_that("Pooled", {
  estimate <- computeFixedEffectMetaAnalysis(populations)
  expect_lt(estimate$lb, 2)
  expect_gt(estimate$ub, 2)
  
  estimate <- computeBayesianMetaAnalysis(populations, chainLength = 1e4, burnIn = 1e3)
  expect_lt(estimate$mu95Lb, 2)
  expect_gt(estimate$mu95Ub, 2)
})

