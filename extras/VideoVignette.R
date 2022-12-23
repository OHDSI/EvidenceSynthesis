# Code used in the video vignette:

# Simulate data ----------------------------------------------------------------
library(EvidenceSynthesis)
simulationSettings <- createSimulationSettings(nSites = 10,
                                               n = 10000,
                                               treatedFraction = 0.8,
                                               nStrata = 5,
                                               hazardRatio = 2,
                                               randomEffectSd = 0.5)
set.seed(1)
populations <- simulatePopulations(simulationSettings)

head(populations[[1]])

# Fit a model locally ----------------------------------------------------------
library(Cyclops)
# Assume we are at site 1:
population <- populations[[1]]


cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                 data = population,
                                 modelType = "cox")
cyclopsFit <- fitCyclopsModel(cyclopsData)

# Hazard ratio:
exp(coef(cyclopsFit))

# 95% confidence interval:
exp(confint(cyclopsFit, parm = "x")[2:3])

# Approximate the likelihood function at one site ------------------------------
approximation <-  approximateLikelihood(
  cyclopsFit = cyclopsFit,
  parameter = "x",
  approximation = "normal"
)
plotLikelihoodFit(approximation = approximation,
                  cyclopsFit = cyclopsFit,
                  parameter = "x")

approximation <-  approximateLikelihood(
  cyclopsFit = cyclopsFit,
  parameter = "x",
  approximation = "adaptive grid",
  bounds = c(log(0.1), log(10))
)
head(approximation)

plotLikelihoodFit(approximation = approximation,
                  cyclopsFit = cyclopsFit,
                  parameter = "x")

# Approximate at all sites -----------------------------------------------------
fitModelInDatabase <- function(population, approximation) {
  cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                   data = population,
                                   modelType = "cox")
  cyclopsFit <- fitCyclopsModel(cyclopsData)
  approximation <-  approximateLikelihood(cyclopsFit,
                                          parameter = "x",
                                          approximation = approximation)
  return(approximation)
}
adaptiveGridApproximations <- lapply(populations, fitModelInDatabase, approximation = "adaptive grid")
normalApproximations <- lapply(populations, fitModelInDatabase, approximation = "normal")
normalApproximations <- dplyr::bind_rows(normalApproximations)

# Synthesize evidence ----------------------------------------------------------

computeFixedEffectMetaAnalysis(adaptiveGridApproximations)
computeFixedEffectMetaAnalysis(normalApproximations)
computeFixedEffectMetaAnalysis(populations)

exp(computeBayesianMetaAnalysis(adaptiveGridApproximations)[, 1:3])
exp(computeBayesianMetaAnalysis(normalApproximations)[, 1:3])
exp(computeBayesianMetaAnalysis(populations)[, 1:3])

