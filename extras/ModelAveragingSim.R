# Simulation code to perform (unweighted) model averaging, when multiple methods have
# been used to estimate the same effect in the same database (but with different ways
# to construct the counterfactual).

library(survival)
library(dplyr)
# Also requires Cyclops, EvidenceSynthesis, and EmpiricalCalibration to be installed
settings <- list(
  nMethods = 4,
  nTarget = 1000,
  nComparator = c(500, 1000, 1000, 1500),
  nNegativeControls = 50,
  nullMean = c(0.1, 0, 0.2, 0),
  nullSd = c(0.1, 0.1, 0.2, 0.15),
  backgroundHazard = 0.001,
  censorHazard = 0.01,
  trueHazardRatio = 2
)

runSingleSimulation <- function(seed, settings) {
  set.seed(seed)

  profiles <- list()
  for (outcomeIdx in seq_len(settings$nNegativeControls + 1)) {
    profiles[[outcomeIdx]] <- list()
    # First outcome is outcome of interest (rest is negative controls):
    if (outcomeIdx == 1) {
      hazardRatio <- settings$trueHazardRatio
    } else {
      hazardRatio <- 1
    }
    hazardTarget <- settings$backgroundHazard * hazardRatio
    # Target is re-used by all methods:
    target <- tibble(
      exposure = 1,
      tOutcome = rexp(settings$nTarget, hazardTarget),
      tCensor = rexp(settings$nTarget, settings$censorHazard)
    )
    for (methodIdx in seq_len(settings$nMethods)) {
      sysError <- rnorm(1, settings$nullMean[methodIdx], settings$nullSd[methodIdx])
      hazardComparator <- settings$backgroundHazard / exp(sysError)
      # To avoid having power correlate with bias, adjust comparator size to standardize
      # overall hazard:
      nComparatorAdjusted <- round(settings$nComparator[methodIdx] * settings$backgroundHazard / hazardComparator)

      # Each method samples different counterfactual. Suboptimal counterfactual
      # selection leads to systematic error:
      comparator <- tibble(
        exposure = 0,
        tOutcome = rexp(nComparatorAdjusted, hazardComparator),
        tCensor = rexp(nComparatorAdjusted, settings$censorHazard)
      )
      data <- bind_rows(target, comparator) |>
        mutate(time = min(tOutcome, tCensor),
               outcome = tOutcome <= tCensor)
      cyclopsData <- Cyclops::createCyclopsData(Surv(time, outcome) ~ exposure, data = data, modelType = "cox")
      fit <- Cyclops::fitCyclopsModel(cyclopsData, startingCoefficients = 0.1)
      profile <- Cyclops::getCyclopsProfileLogLikelihood(fit, "exposure", bounds = c(log(0.1), log(10)))
      profiles[[outcomeIdx]][[methodIdx]] <- profile
    }
  }

  # Martijn's simple stupid frequentist two-stage approach. to be replaced by
  # fancy Bayesian one-stage approach:

  # Stage one: fit systematic error distributions:
  nulls <- list()
  for (methodIdx in seq_len(settings$nMethods)) {
    ncProfiles <- lapply(seq_len(settings$nNegativeControls) + 1,
                         function(outcomeIdx) profiles[[outcomeIdx]][[methodIdx]])
    nulls[[methodIdx]] <- EmpiricalCalibration::fitNullNonNormalLl(ncProfiles)
  }

  # Stage two: Construct likelihood function and find mode and 95% CI:
  # Normalize LLs:
  for (methodIdx in seq_len(settings$nMethods)) {
    profiles[[1]][[methodIdx]]$value <- profiles[[1]][[methodIdx]]$value -
      max(profiles[[1]][[methodIdx]]$value)
  }
  biasL <- function(bias, x, mu, sigma) {
    exp(EmpiricalCalibration:::gridLlApproximation(x - bias, profiles[[1]][[methodIdx]])) *
      dnorm(bias, mu, sigma)
  }
  methodL <- function(methodIdx, x) {
    null <- nulls[[methodIdx]]
    integrate(
      biasL,
      lower = null[1] - 5 * null[2],
      upper = null[1] + 5 * null[2],
      x = x,
      mu = null[1],
      sigma = null[2]
    )$value
  }
  logLikelihood <- function(x) {
    log(sum(sapply(seq_len(settings$nMethods), methodL, x = x)) / settings$nMethods)
  }
  estimate <- EvidenceSynthesis:::computeEstimateFromApproximation(logLikelihood)
  return(estimate)
}
cluster <- ParallelLogger::makeCluster(20)
ParallelLogger::clusterRequire(cluster, "dplyr")
ParallelLogger::clusterRequire(cluster, "survival")
estimates <- ParallelLogger::clusterApply(cluster, seq_len(100), runSingleSimulation, settings = settings)
estimates <- estimates |>
  bind_rows()
estimates |>
  summarise(
    coverage = mean(lb < settings$trueHazardRatio & ub > settings$trueHazardRatio),
    meanHr = exp(mean(logRr))
  )
# coverage   meanHr
# 1     0.99 1.940481
ParallelLogger::stopCluster(cluster)
