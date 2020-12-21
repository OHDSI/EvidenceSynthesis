# Computing the gold standard takes a long time, so storing for quicker testing
recomputeGoldStandard <- FALSE
# setwd('tests/testthat')

library(testthat)
library(EvidenceSynthesis)
library(survival)

createApproximations <- function(populations, approximation) {
  fitModelInDatabase <- function(population, approximation) {
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                              data = population,
                                              modelType = "cox")
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData,
                                           fixedCoefficients = c(approximation != "normal"))
    approximation <- approximateLikelihood(cyclopsFit, "x", approximation = approximation)
    return(approximation)
  }
  data <- lapply(populations, fitModelInDatabase, approximation = approximation)
  data <- do.call("rbind", data)
  return(data)
}

if (recomputeGoldStandard) {
  set.seed(1)
  populations <- simulatePopulations(settings = createSimulationSettings(nSites = 10,
                                                                         n = 2500,
                                                                         treatedFraction = 0.25,
                                                                         hazardRatio = 2,
                                                                         randomEffectSd = 0.5))
  pooledFixedFxEstimate <- computeFixedEffectMetaAnalysis(populations)
  pooledRandomFxEstimate <- computeBayesianMetaAnalysis(populations)
  saveRDS(populations, "resources/populations.rds")
  saveRDS(pooledFixedFxEstimate, "resources/pooledFixedFxEstimate.rds")
  saveRDS(pooledRandomFxEstimate, "resources/pooledRandomFxEstimate.rds")
} else {
  populations <- readRDS("resources/populations.rds")
  pooledFixedFxEstimate <- readRDS("resources/pooledFixedFxEstimate.rds")
  pooledRandomFxEstimate <- readRDS("resources/pooledRandomFxEstimate.rds")
}

# seed <- round(runif(1, 0, 1e10))
seed <- 1

# Custom approximation
data <- createApproximations(populations, "custom")

test_that("Custom approximation: pooled matches fixed-effects meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_equal(estimate,
               pooledFixedFxEstimate,
               tolerance = 0.15,
               scale = 1,
               check.attributes = FALSE)
})

test_that("Custom approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  expect_equal(estimate[, c("mu", "tau", "logRr")],
               pooledRandomFxEstimate[, c("mu", "tau", "logRr")],
               tolerance = 0.10,
               scale = 1,
               check.attributes = FALSE)
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
               pooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
               tolerance = 0.30,
               scale = 1,
               check.attributes = FALSE)
})

# Grid approximation
data <- createApproximations(populations, "grid")

test_that("Grid approximation: pooled matches meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_equal(estimate,
               pooledFixedFxEstimate,
               tolerance = 0.15,
               scale = 1,
               check.attributes = FALSE)
})

test_that("Grid approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  expect_equal(estimate[, c("mu", "tau", "logRr")],
               pooledRandomFxEstimate[, c("mu", "tau", "logRr")],
               tolerance = 0.10,
               scale = 1,
               check.attributes = FALSE)
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
               pooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
               tolerance = 0.30,
               scale = 1,
               check.attributes = FALSE)
})

# Normal approximation
data <- createApproximations(populations, "normal")

test_that("Normal approximation: pooled matches meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  # Not really expecting normal approximation is close to gold standard:
  expect_equal(estimate, pooledFixedFxEstimate, tolerance = 10, check.attributes = FALSE)
})

test_that("Normal approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  # Not really expecting normal approximation is close to gold standard:
  expect_equal(estimate, pooledRandomFxEstimate, tolerance = 1, check.attributes = FALSE)
})

# Skew-normal approximation
data <- createApproximations(populations, "skew normal")

test_that("Skew-normal approximation: pooled matches meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  # Skew-normal is a poorer approximation, using higher tolerance:
  expect_equal(estimate,
               pooledFixedFxEstimate,
               tolerance = 0.30,
               scale = 1,
               check.attributes = FALSE)
})

test_that("Skew-normal approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  # Skew-normal is a poorer approximation, using higher tolerance:
  expect_equal(estimate[, c("mu", "tau", "logRr")],
               pooledRandomFxEstimate[, c("mu", "tau", "logRr")],
               tolerance = 0.50,
               scale = 1,
               check.attributes = FALSE)
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
               pooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
               tolerance = 1.00,
               scale = 1,
               check.attributes = FALSE)
})
