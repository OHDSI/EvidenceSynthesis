library(testthat)

createApproximations <- function(populations, approximation) {
  fitModelInDatabase <- function(population, approximation) {
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, fixedCoefficients = c(TRUE))
    approximation <-  approximateLikelihood(cyclopsFit, "x", approximation = approximation)
    return(approximation)
  }
  data <- lapply(populations, fitModelInDatabase, approximation = approximation)
  data <- do.call("rbind", data)
  return(data)
}

set.seed(1)
populations <- simulatePopulations(settings = createSimulationSettings(nSites = 10,
                                                                       n = 2000, 
                                                                       treatedFraction = 0.3, 
                                                                       hazardRatio = 2))
pooledEstimate <- computeFixedEffectMetaAnalysis(populations)

test_that("Custom approximation: pooled matches meta-analysis", {
  data <- createApproximations(populations, "custom")
  estimate <- computeFixedEffectMetaAnalysis(data)
  
  expect_equal(estimate, pooledEstimate, tolerance = 0.1, check.attributes = FALSE)  
})

test_that("Grid approximation: pooled matches meta-analysis", {
  data <- createApproximations(populations, "grid")
  estimate <- computeFixedEffectMetaAnalysis(data)
  
  expect_equal(estimate, pooledEstimate, tolerance = 0.1, check.attributes = FALSE)  
})

