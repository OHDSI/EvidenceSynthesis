library(testthat)
library(survival)
library(EvidenceSynthesis)

# Simulate large balanced population. Should have normally distributed likelihood:
set.seed(1)
population <- simulatePopulations(settings = createSimulationSettings(nSites = 1,
                                                                      n = 10000, 
                                                                      treatedFraction = 0.5, 
                                                                      hazardRatio = 0.5))[[1]]
cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
x <- log(seq(0.1, 10, length.out = 100))
ll <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, "x", x)$value
ll <- ll - max(ll)

test_that("Normal approximation", {
  params <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "normal")
  # plotLikelihoodFit(data, cyclopsFit, parameter = "x")
  approxLl <- dnorm(x, mean = params$logRr, sd = params$seLogRr, log = TRUE)
  approxLl <- approxLl - max(approxLl)
  expect_equal(ll, approxLl, tolerance = 0.1)
})

test_that("Skew normal approximation", {
  params <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "skew normal")
  approxLl <- skewNormal(x, mu = params$mu, sigma = params$sigma, alpha = params$alpha)
  approxLl <- approxLl - max(approxLl)
  expect_equal(ll, approxLl, tolerance = 0.1)
})

test_that("Custom approximation", {
  params <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "custom")
  approxLl <- customFunction(x, mu = params$mu, sigma = params$sigma, gamma = params$gamma)
  approxLl <- approxLl - max(approxLl)
  expect_equal(ll, approxLl, tolerance = 0.1)
})

test_that("Grid approximation", {
  approxLl <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "grid")
  approxLl <- approxLl - max(approxLl)
  x <- as.numeric(names(approxLl))
  ll <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, "x", x)$value
  ll <- ll - max(ll)
  
  expect_equal(ll, approxLl, tolerance = 0.01, check.attributes = FALSE)
})
