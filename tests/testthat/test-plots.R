library(testthat)
library(EvidenceSynthesis)

set.seed(1234)
populations <- simulatePopulations()
labels <- paste("Data site", LETTERS[1:length(populations)])

# Fit a Cox regression at each data site, and approximate likelihood function:
fitModelInDatabase <- function(population) {
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
    data = population,
    modelType = "cox"
  )
  cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, control = Cyclops::createControl(convergenceType = "lange"))
  approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "grid")
  attr(approximation, "cyclopsFit") <- cyclopsFit
  return(approximation)
}
approximations <- lapply(populations, fitModelInDatabase)
approximation <- approximations[[1]]
approximations <- do.call("rbind", approximations)

tempFile <- tempfile(fileext = ".png")

# At study coordinating center, perform meta-analysis using per-site approximations:
estimate <- computeBayesianMetaAnalysis(approximations)

test_that("Plot forest using grid approximation", {
  plot <- plotMetaAnalysisForest(data = approximations, labels = labels, estimate = estimate, fileName = tempFile)
  expect_s3_class(plot, "gtable")
})

test_that("Plot MCMC traces using grid approximation", {
  plot <- plotMcmcTrace(estimate, fileName = tempFile)
  expect_s3_class(plot, "ggplot")
})

test_that("Plot MCMC traces per DB using grid approximation", {
  plot <- plotPerDbMcmcTrace(estimate, fileName = tempFile)
  expect_s3_class(plot, "ggplot")
})

test_that("Plot posterior density using grid approximation", {
  plot <- plotPosterior(estimate, fileName = tempFile)
  expect_s3_class(plot, "ggplot")
})

test_that("Plot posterior density per BD using grid approximation", {
  plot <- plotPerDbPosterior(estimate, fileName = tempFile)
  expect_s3_class(plot, "ggplot")
})

test_that("Plot likelihood fit", {
  suppressWarnings(
    plot <- plotLikelihoodFit(approximation, attr(approximation, "cyclopsFit"), parameter = "x", fileName = tempFile)
  )
  expect_s3_class(plot, "ggplot")
})

test_that("Plot empirical nulls", {
  site1 <- EmpiricalCalibration::simulateControls(n = 50, mean = 0, sd = 0.1, trueLogRr = 0)
  site1$label <- "Site 1"
  site2 <- EmpiricalCalibration::simulateControls(n = 50, mean = 0.1, sd = 0.2, trueLogRr = 0)
  site2$label <- "Site 2"
  site3 <- EmpiricalCalibration::simulateControls(n = 50, mean = 0.15, sd = 0.25, trueLogRr = 0)
  site3$label <- "Site 3"
  sites <- rbind(site1, site2, site3)

  plot <- plotEmpiricalNulls(logRr = sites$logRr, seLogRr = sites$seLogRr, labels = sites$label, fileName = tempFile)
  expect_s3_class(plot, "gtable")
})

unlink(tempFile, recursive = TRUE)
