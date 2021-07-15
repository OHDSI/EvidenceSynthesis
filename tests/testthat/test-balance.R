library(testthat)
library(EvidenceSynthesis)

set.seed(1234)
balance1 <- data.frame(beforeMatchingStdDiff = rnorm(1000, 0.1, 0.1),
                       afterMatchingStdDiff = rnorm(1000, 0, 0.01))
balance2 <- data.frame(beforeMatchingStdDiff = rnorm(1000, 0.2, 0.1),
                       afterMatchingStdDiff = rnorm(1000, 0, 0.05))
balance3 <- data.frame(beforeMatchingStdDiff = rnorm(1000, 0, 0.1),
                       afterMatchingStdDiff = rnorm(1000, 0, 0.03))

treatment <- rep(0:1, each = 100)
propensityScore <- c(rnorm(100, mean = 0.4, sd = 0.25), rnorm(100, mean = 0.6, sd = 0.25))
data <- data.frame(treatment = treatment, propensityScore = propensityScore)
data <- data[data$propensityScore > 0 & data$propensityScore < 1, ]

tempFile <- tempfile(fileext = ".png")

test_that("Plot covariate balance", {
  plot <- plotCovariateBalances(balances = list(balance1, balance2, balance3), labels = c("Site A", "Site B", "Site C"), 
                                fileName = tempFile)
  expect_s3_class(plot, "gtable")
})

test_that("Plot PS", {
  preparedPlot <- preparePsPlot(data)
  plot <- plotPreparedPs(preparedPsPlots = list(preparedPlot, preparedPlot, preparedPlot), 
                         labels = c("Site A", "Site B", "Site C"), 
                         fileName = tempFile)
  expect_s3_class(plot, "ggplot")
})

unlink(tempFile, recursive = TRUE)
