library(survival)
library(testthat)

test_that("Simple Bayesian profile-normal model", {
  skip_if_not(supportsJava8())
  
  set.seed(666)
  population <- simulatePopulations(createSimulationSettings(nSites = 1))[[1]]
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                            data = population,
                                            modelType = "cox")
  cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
  likelihoodProfile <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "grid")
  
  traces <- approximateSimplePosterior(likelihoodProfile = likelihoodProfile, priorMean = 0, priorSd = 100, seed = 666)
  expect_equivalent(mean(traces$theta1), coef(cyclopsFit)[1], tolerance = 0.1)
})
