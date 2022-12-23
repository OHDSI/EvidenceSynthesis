# Code for generating plot for video vignette:
library(EvidenceSynthesis)
simulationSettings <- createSimulationSettings(nSites = 10,
                                               n = 10000,
                                               treatedFraction = 0.8,
                                               nStrata = 5,
                                               hazardRatio = 2,
                                               randomEffectSd = 0.5)
set.seed(1)
populations <- simulatePopulations(simulationSettings)

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




limits = c(0.1, 10)
breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
plotData <- adaptiveGridApproximations[[10]]
plot <- ggplot2::ggplot(plotData, ggplot2::aes(x = point, y = value)) +
  ggplot2::geom_vline(xintercept = log(breaks), color = "#AAAAAA", lty = 1, linewidth = 0.5) +
  ggplot2::geom_line(alpha = 0.2) +
  ggplot2::geom_point(size = 0.5, alpha = 0.5) +
  ggplot2::scale_y_continuous("Log Likelihood") +
  ggplot2::scale_x_continuous("Hazard Ratio",
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
ggplot2::ggsave(filename = "d:/temp/plot.png", plot = plot, width = 6, height = 3, dpi = 200)
