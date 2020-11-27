library(testthat)
library(EvidenceSynthesis)

test_that("Plot forest using grid approximation", {
  populations <- simulatePopulations()
  labels <- paste("Data site", LETTERS[1:length(populations)])
  
  # Fit a Cox regression at each data site, and approximate likelihood function:
  fitModelInDatabase <- function(population) {
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                              data = population,
                                              modelType = "cox")
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
    approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "grid")
    return(approximation)
  }
  approximations <- lapply(populations, fitModelInDatabase)
  approximations <- do.call("rbind", approximations)
  
  # At study coordinating center, perform meta-analysis using per-site approximations:
  estimate <- computeBayesianMetaAnalysis(approximations)
  plot <- plotMetaAnalysisForest(approximations, labels, estimate)
  expect_s3_class(plot, "gtable")
})