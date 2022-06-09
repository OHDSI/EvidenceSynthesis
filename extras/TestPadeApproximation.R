# Some code to test the Pade approximation
library(EvidenceSynthesis)
set.seed(1)
settings <- createSimulationSettings(nSites = 5,
                                     n = round(runif(5, 3000, 4000)),
                                     minBackgroundHazard = 0.00001,
                                     maxBackgroundHazard = 0.0001,
                                     treatedFraction = 0.1,
                                     hazardRatio = 1,
                                     randomEffectSd = 0,
                                     nStrata = 1)
populations <- simulatePopulations(settings)

# Explore individual approximations --------------------------------------------
population <- populations[[5]]

# Has zero count:
table(population[,c(3,5)])

approximation <- approximateLikelihoodUsingPade(population)
plotPadeApproximation(population, approximation) # See function def below

cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                          data = population,
                                          modelType = "cox")
cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
approximation <- approximateLikelihood(cyclopsFit = cyclopsFit,
                                       parameter = "x",
                                       approximation = "custom")
plotLikelihoodFit(approximation = approximation,
                  cyclopsFit = cyclopsFit,
                  parameter = "x")


# Test across databases using simulation ---------------------------------------
fitModelInDatabase <- function(population, method = "pade") {
  if (method == "pade") {
    approximation <- approximateLikelihoodUsingPade(population)
  } else {
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                              data = population,
                                              modelType = "cox")
    cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
    approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = method)
  }
  return(approximation)
}

# Pade
method <- "pade"
approximations <- lapply(populations, fitModelInDatabase, method = method)
approximations <- do.call("rbind", approximations)
estimate <- computeFixedEffectMetaAnalysis(approximations)
estimate
# rr        lb       ub      logRr   seLogRr
# 1 0.613201 0.1842669 1.519949 -0.4890625 0.5382872
estimate <- computeBayesianMetaAnalysis(approximations)
estimate
# mu   mu95Lb    mu95Ub      muSe       tau     tau95Lb   tau95Ub      logRr   seLogRr
# 1 -0.5563541 -1.61537 0.5136067 0.5455222 0.3438458 0.001198962 0.9381472 -0.5563541 0.5455222
plotPosterior(estimate)

# Custom function
method <- "custom"
approximations <- lapply(populations, fitModelInDatabase, method = method)
approximations <- do.call("rbind", approximations)
estimate <- computeFixedEffectMetaAnalysis(approximations)
estimate
# rr        lb       ub      logRr   seLogRr
# 1 0.624319 0.1909064 1.529607 -0.4710937 0.5308728
estimate <- computeBayesianMetaAnalysis(approximations)
estimate
# mu    mu95Lb    mu95Ub      muSe       tau     tau95Lb   tau95Ub      logRr   seLogRr
# 1 -0.5311875 -1.652897 0.5402607 0.5738145 0.3331693 0.001622732 0.9501842 -0.5311875 0.5738145
plotPosterior(estimate)

# Gold standard (data pooling)
estimate <- computeFixedEffectMetaAnalysis(populations)
estimate
# rr        lb       ub      logRr   seLogRr
# x 0.6292336 0.1906228 1.532602 -0.4632526 0.5317509
estimate <- computeBayesianMetaAnalysis(populations)
estimate
# mu    mu95Lb    mu95Ub      muSe       tau     tau95Lb   tau95Ub      logRr   seLogRr
# 1 -0.5613211 -1.605733 0.5343067 0.5551024 0.3563129 0.002501133 0.9717572 -0.5613211 0.5551024
plotPosterior(estimate)


# Function for plotting approximation ------------------------------------------
plotPadeApproximation <- function(survivalData,
                                  approximation,
                                  xLabel = "Hazard Ratio",
                                  limits = c(0.1, 10)) {
  x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                            data = survivalData,
                                            modelType = "cox")
  cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
  ll <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, 1, x)$value
  constraints <- EvidenceSynthesis:::getPadeConstraints(beta = approximation$beta,
                                                        a0 = approximation$a0,
                                                        a1 = approximation$a1,
                                                        a2 = approximation$a2,
                                                        b1 = approximation$b1,
                                                        b2 = approximation$b2)
  y <- EvidenceSynthesis:::padeConstrained(x, approximation$beta, approximation$a0, approximation$a1, approximation$a2, approximation$b1, approximation$b2, constraints$minBeta, constraints$maxBeta, constraints$minD1, constraints$maxD1)
  y <- y - max(y)
  ll <- ll - max(ll)

  plotData <- rbind(data.frame(x = x, ll = ll, type = "Likelihood"),
                    data.frame(x = x, ll = y, type = "Approximation"))
  yLabel <- "Log Likelihood"
  plotData$type <- factor(plotData$type, levels = c("Likelihood", "Approximation"))
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  plot <- ggplot2::ggplot(plotData, ggplot2::aes(x = .data$x, y = .data$ll)) +
    ggplot2::geom_vline(xintercept = log(breaks), color = "#AAAAAA", lty = 1, size = 0.5) +
    ggplot2::geom_line(ggplot2::aes(group = .data$type,
                                    color = .data$type,
                                    linetype = .data$type,
                                    size = .data$type), alpha = 0.7) +
    ggplot2::scale_size_manual(values = c(1, 2, 2, 2) * 0.6) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) +
    ggplot2::scale_color_manual(values = c("#000000", "#66c2a5", "#fc8d62", "#8da0cb")) +
    ggplot2::scale_y_continuous(yLabel) +
    ggplot2::scale_x_continuous(xLabel,
                                limits = log(limits),
                                breaks = log(breaks),
                                labels = breaks) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "top")
  plot
}
