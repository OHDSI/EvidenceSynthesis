# Computing the gold standard takes a long time, so storing for quicker testing
recomputeGoldStandard <- FALSE
# setwd('tests/testthat')

library(testthat)
library(EvidenceSynthesis)

if (recomputeGoldStandard) {
  set.seed(1)
  populations <- simulatePopulations(settings = createSimulationSettings(
    nSites = 10,
    n = 2500,
    treatedFraction = 0.25,
    hazardRatio = 2,
    randomEffectSd = 0.5
  ))
  data <- createApproximations(populations, "grid with gradients")
  estimate <- computeBayesianMetaAnalysis(data)
  traces <- attr(estimate, "traces")
  predictions <- do.call(c, lapply(seq_len(nrow(traces)), function(i) rnorm(100000, traces[i, 1], traces[i, 2])))
  gsPredictionInterval <- HDInterval::hdi(predictions, credMass = 0.95)
  saveRDS(traces, "resources/tracesForPi.rds")
  saveRDS(gsPredictionInterval, "resources/gsPredictionInterval.rds")
} else {
  traces <- readRDS("resources/tracesForPi.rds")
  gsPredictionInterval <- readRDS("resources/gsPredictionInterval.rds")
}

test_that("Prediction interval matches gold standard", {
  predictionInterval <- computePredictionInterval(traces)
  expect_equal(predictionInterval,
               gsPredictionInterval,
               tolerance = 0.15,
               scale = 1,
               check.attributes = FALSE
  )
})
