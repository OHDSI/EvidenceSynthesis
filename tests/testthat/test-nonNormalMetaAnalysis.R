# Computing the gold standard takes a long time, so storing for quicker testing
recomputeGoldStandard <- FALSE
# setwd('tests/testthat')

library(testthat)
library(EvidenceSynthesis)
library(survival)

if (recomputeGoldStandard) {
  set.seed(1)
  populations <- simulatePopulations(settings = createSimulationSettings(
    nSites = 10,
    n = 2500,
    treatedFraction = 0.25,
    hazardRatio = 2,
    randomEffectSd = 0.5
  ))
  pooledFixedFxEstimate <- computeFixedEffectMetaAnalysis(populations)
  pooledRandomFxEstimate <- computeBayesianMetaAnalysis(populations)
  saveRDS(populations, "resources/populations.rds")
  saveRDS(pooledFixedFxEstimate, "resources/pooledFixedFxEstimate.rds")
  saveRDS(pooledRandomFxEstimate, "resources/pooledRandomFxEstimate.rds")

  sccsPopulations <- simulatePopulations(settings = createSccsSimulationSettings(
    nSites = 10,
    n = 2500,
    atRiskTimeFraction = 0.25,
    timePartitions = 10,
    timeCovariates = 5,
    timeEffectSize = log(2),
    rateRatio = 2,
    randomEffectSd = 0.5
  ))
  sccsPooledFixedFxEstimate <- computeFixedEffectMetaAnalysis(sccsPopulations)
  sccsPooledRandomFxEstimate <- computeBayesianMetaAnalysis(sccsPopulations)
  saveRDS(sccsPopulations, "resources/sccsPopulations.rds")
  saveRDS(sccsPooledFixedFxEstimate, "resources/sccsPooledFixedFxEstimate.rds")
  saveRDS(sccsPooledRandomFxEstimate, "resources/sccsPooledRandomFxEstimate.rds")
} else {
  populations <- readRDS("resources/populations.rds")
  pooledFixedFxEstimate <- readRDS("resources/pooledFixedFxEstimate.rds")
  pooledRandomFxEstimate <- readRDS("resources/pooledRandomFxEstimate.rds")
  sccsPopulations <- readRDS("resources/sccsPopulations.rds")
  sccsPooledFixedFxEstimate <- readRDS("resources/sccsPooledFixedFxEstimate.rds")
  sccsPooledRandomFxEstimate <- readRDS("resources/sccsPooledRandomFxEstimate.rds")
}

# seed <- round(runif(1, 0, 1e10))
seed <- 1

# Custom approximation
data <- createApproximations(populations, "custom")

test_that("Custom approximation: pooled matches fixed-effects meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_equal(estimate[, c("rr", "logRr")],
    pooledFixedFxEstimate[, c("rr", "logRr")],
    tolerance = 0.15,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("lb", "ub", "seLogRr")],
    pooledFixedFxEstimate[, c("lb", "ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

test_that("Custom approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  expect_equal(estimate[, c("mu", "tau", "logRr")],
    pooledRandomFxEstimate[, c("mu", "tau", "logRr")],
    tolerance = 0.10,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    pooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

# Grid approximation
data <- createApproximations(populations, "grid")

test_that("Grid approximation: pooled matches fixed-effects meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_equal(estimate[, c("rr", "logRr")],
    pooledFixedFxEstimate[, c("rr", "logRr")],
    tolerance = 0.15,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("lb", "ub", "seLogRr")],
    pooledFixedFxEstimate[, c("lb", "ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

test_that("Grid approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  expect_equal(estimate[, c("mu", "tau", "logRr")],
    pooledRandomFxEstimate[, c("mu", "tau", "logRr")],
    tolerance = 0.10,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    pooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

# Adaptive grid approximation
data <- createApproximations(populations, "adaptive grid")

test_that("Adaptive grid approximation: pooled matches fixed-effects meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_equal(estimate[, c("rr", "logRr")],
    pooledFixedFxEstimate[, c("rr", "logRr")],
    tolerance = 0.15,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("lb", "ub", "seLogRr")],
    pooledFixedFxEstimate[, c("lb", "ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

test_that("Adaptive grid approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  expect_equal(estimate[, c("mu", "tau", "logRr")],
    pooledRandomFxEstimate[, c("mu", "tau", "logRr")],
    tolerance = 0.10,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    pooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

# Normal approximation
data <- createApproximations(populations, "normal")

test_that("Normal approximation: pooled matches fixed-effects meta-analysis", {
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

test_that("Normal approximation: pooled matches random-effects meta-analysis using tibble", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(dplyr::as_tibble(data), seed = seed)
  # Not really expecting normal approximation is close to gold standard:
  expect_equal(estimate, pooledRandomFxEstimate, tolerance = 1, check.attributes = FALSE)
})

# Skew-normal approximation
data <- createApproximations(populations, "skew normal")

test_that("Skew-normal approximation: pooled matches fixed-effects meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  # Skew-normal is a poorer approximation, using higher tolerance:
  expect_equal(estimate[, c("rr", "logRr")],
    pooledFixedFxEstimate[, c("rr", "logRr")],
    tolerance = 0.30,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("lb", "ub", "seLogRr")],
    pooledFixedFxEstimate[, c("lb", "ub", "seLogRr")],
    tolerance = 1.00,
    scale = 1,
    check.attributes = FALSE
  )
})

test_that("Skew-normal approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  # Skew-normal is a poorer approximation, using higher tolerance:
  expect_equal(estimate[, c("mu", "tau", "logRr")],
    pooledRandomFxEstimate[, c("mu", "tau", "logRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    pooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    tolerance = 1.00,
    scale = 1,
    check.attributes = FALSE
  )
})


# Grid with gradients approximation
data <- createApproximations(populations, "grid with gradients")

test_that("Grid with gradients approximation: pooled matches fixed-effects meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_equal(estimate[, c("rr", "logRr")],
    pooledFixedFxEstimate[, c("rr", "logRr")],
    tolerance = 0.15,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("lb", "ub", "seLogRr")],
    pooledFixedFxEstimate[, c("lb", "ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

test_that("Grid with gradients approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  expect_equal(estimate[, c("mu", "tau", "logRr")],
    pooledRandomFxEstimate[, c("mu", "tau", "logRr")],
    tolerance = 0.15,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    pooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})


# SCCS Adaptive grid approximation
data <- createApproximations(sccsPopulations, "adaptive grid")

test_that("SCCS adaptive grid approximation: pooled matches fixed-effects meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_equal(estimate[, c("rr", "logRr")],
    sccsPooledFixedFxEstimate[, c("rr", "logRr")],
    tolerance = 0.15,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("lb", "ub", "seLogRr")],
    sccsPooledFixedFxEstimate[, c("lb", "ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

test_that("SCCS adaptive grid approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  expect_equal(estimate[, c("mu", "tau", "logRr")],
    sccsPooledRandomFxEstimate[, c("mu", "tau", "logRr")],
    tolerance = 0.10,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    sccsPooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})


# SCCS grid with gradients approximation
data <- createApproximations(sccsPopulations, "grid with gradients")

test_that("SCCS adaptive grid approximation: pooled matches fixed-effects meta-analysis", {
  estimate <- computeFixedEffectMetaAnalysis(data)
  expect_equal(estimate[, c("rr", "logRr")],
    sccsPooledFixedFxEstimate[, c("rr", "logRr")],
    tolerance = 0.15,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("lb", "ub", "seLogRr")],
    sccsPooledFixedFxEstimate[, c("lb", "ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})

test_that("SCCS adaptive grid approximation: pooled matches random-effects meta-analysis", {
  skip_if_not(supportsJava8())
  estimate <- computeBayesianMetaAnalysis(data, seed = seed)
  expect_equal(estimate[, c("mu", "tau", "logRr")],
    sccsPooledRandomFxEstimate[, c("mu", "tau", "logRr")],
    tolerance = 0.10,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(estimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    sccsPooledRandomFxEstimate[, c("mu95Lb", "mu95Ub", "muSe", "tau95Lb", "tau95Ub", "seLogRr")],
    tolerance = 0.50,
    scale = 1,
    check.attributes = FALSE
  )
})
