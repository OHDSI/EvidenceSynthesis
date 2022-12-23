# Copyright 2022 Observational Health Data Sciences and Informatics
#
# This file is part of EvidenceSynthesis
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Plot the likelihood approximation
#'
#' @details
#' Plots the (log) likelihood and the approximation of the likelihood. Allows for reviewing the
#' approximation.
#'
#' @param approximation   An approximation of the likelihood function as fitted using the
#'                        [approximateLikelihood()] function.
#' @param cyclopsFit      A model fitted using the [Cyclops::fitCyclopsModel()] function.
#' @param parameter       The parameter in the `cyclopsFit` object to profile.
#' @param logScale        Show the y-axis on the log scale?
#' @param xLabel          The title of the x-axis.
#' @param limits          The limits on the x-axis.
#' @param fileName        Name of the file where the plot should be saved, for example 'plot.png'. See
#'                        the function [ggplot2::ggsave] in the ggplot2 package for supported file
#'                        formats.
#'
#' @return
#' A Ggplot object. Use the [ggplot2::ggsave] function to save to file.
#'
#' @examples
#' # Simulate a single database population:
#' population <- simulatePopulations(createSimulationSettings(nSites = 1))[[1]]
#'
#' # Approximate the likelihood:
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                           data = population,
#'                                           modelType = "cox")
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "custom")
#'
#' plotLikelihoodFit(approximation, cyclopsFit, parameter = "x")
#'
#' @export
plotLikelihoodFit <- function(approximation,
                              cyclopsFit,
                              parameter = "x",
                              logScale = TRUE,
                              xLabel = "Hazard Ratio",
                              limits = c(0.1, 10),
                              fileName = NULL) {
  if ("logRr" %in% colnames(approximation)) {
    inform("Detected normal approximation")
    x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
    y <- dnorm(x, mean = approximation$logRr, sd = approximation$seLogRr, log = TRUE)
  } else if ("gamma" %in% colnames(approximation)) {
    inform("Detected custom parameric approximation")
    x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
    y <- customFunction(x,
                        mu = approximation$mu,
                        sigma = approximation$sigma,
                        gamma = approximation$gamma)
  } else if ("alpha" %in% colnames(approximation)) {
    inform("Detected skew normal approximation")
    x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
    y <- skewNormal(x,
                    mu = approximation$mu,
                    sigma = approximation$sigma,
                    alpha = approximation$alpha)
  }  else if ("point" %in% names(approximation)) {
    inform("Detected adaptive grid approximation")
    x <- approximation$point
    y <- approximation$value
  } else {
    inform("Detected grid approximation")
    x <- as.numeric(names(approximation))
    if (any(is.na(x))) {
      abort("Expecting grid data, but not all column names are numeric")
    }
    y <- approximation

  }
  ll <- getLikelihoodProfile(cyclopsFit, parameter, x)
  x <- x[!is.nan(ll)]
  y <- y[!is.nan(ll)]
  ll <- ll[!is.nan(ll)]
  ll <- ll - max(ll)
  y <- y - max(y)
  plotData <- rbind(data.frame(x = x, ll = ll, type = "Likelihood"),
                    data.frame(x = x, ll = y, type = "Approximation"))
  if (logScale) {
    yLabel <- "Log Likelihood"
  } else {
    yLabel <- "Likelihood"
    plotData$ll <- exp(plotData$ll)
  }
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
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 5, height = 5, dpi = 400)
  return(plot)
}

#' Plot MCMC trace
#'
#' @details
#' Plot the samples of the posterior distribution of the mu and tau parameters. Samples are taken
#' using Markov-chain Monte Carlo (MCMC).
#'
#' @param estimate       An object as generated using the [computeBayesianMetaAnalysis()] function.
#' @param showEstimate   Show the parameter estimates (mode) and 95 percent confidence intervals?
#' @param dataCutoff     This fraction of the data at both tails will be removed.
#' @param fileName       Name of the file where the plot should be saved, for example 'plot.png'. See
#'                       the function [ggplot2::ggsave] in the ggplot2 package for supported file
#'                       formats.
#'
#' @seealso
#' [computeBayesianMetaAnalysis]
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' # Fit a Cox regression at each data site, and approximate likelihood function:
#' fitModelInDatabase <- function(population) {
#'   cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                             data = population,
#'                                             modelType = "cox")
#'   cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#'   approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "custom")
#'   return(approximation)
#' }
#' approximations <- lapply(populations, fitModelInDatabase)
#' approximations <- do.call("rbind", approximations)
#'
#' # At study coordinating center, perform meta-analysis using per-site approximations:
#' estimate <- computeBayesianMetaAnalysis(approximations)
#' plotMcmcTrace(estimate)
#'
#' @return
#' A Ggplot object. Use the [ggplot2::ggsave] function to save to file.
#'
#' @export
plotMcmcTrace <- function(estimate, showEstimate = TRUE, dataCutoff = 0.01, fileName = NULL) {
  traces <- attr(estimate, "traces")
  dataMu <- data.frame(x = 1:nrow(traces), trace = traces[, 1], var = "Mu")
  dataTau <- data.frame(x = 1:nrow(traces), trace = traces[, 2], var = "Tau")
  if (dataCutoff != 0) {
    limsMu <- quantile(traces[, 1], c(dataCutoff, 1 - dataCutoff))
    limsTau <- quantile(traces[, 2], c(dataCutoff, 1 - dataCutoff))
    dataMu <- dataMu[dataMu$trace > limsMu[1] & dataMu$trace < limsMu[2], ]
    dataTau <- dataTau[dataTau$trace > limsTau[1] & dataTau$trace < limsTau[2], ]
  }
  data <- rbind(dataMu, dataTau)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x,
                                             y = .data$trace)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::scale_x_continuous("Iterations") +
    ggplot2::facet_grid(var ~ ., scales = "free", switch = "both") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.background = ggplot2::element_blank())
  if (showEstimate) {
    mode <- data.frame(trace = c(estimate$mu, estimate$tau), var = c("Mu", "Tau"))
    ci <- data.frame(trace = c(estimate$mu95Lb,
                               estimate$mu95Ub,
                               estimate$tau95Lb,
                               estimate$tau95Ub), var = c("Mu", "Mu", "Tau", "Tau"))

    plot <- plot +
      ggplot2::geom_hline(ggplot2::aes(yintercept = trace), data = mode) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = trace), data = ci, linetype = "dashed")
  }
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  return(plot)
}

#' Plot MCMC trace for individual databases
#'
#' @details
#' Plot the samples of the posterior distribution of the theta parameter (the estimated log hazard
#' ratio) at each site. Samples are taken using Markov-chain Monte Carlo (MCMC).
#'
#' @param estimate       An object as generated using the [computeBayesianMetaAnalysis()] function.
#' @param showEstimate   Show the parameter estimates (mode) and 95 percent confidence intervals?
#' @param dataCutoff     This fraction of the data at both tails will be removed.
#' @param fileName       Name of the file where the plot should be saved, for example 'plot.png'. See
#'                       the function [ggplot2::ggsave] in the ggplot2 package for supported file
#'                       formats.
#'
#' @seealso
#' [computeBayesianMetaAnalysis]
#'
#' @return
#' A Ggplot object. Use the [ggplot2::ggsave] function to save to file.
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' # Fit a Cox regression at each data site, and approximate likelihood function:
#' fitModelInDatabase <- function(population) {
#'   cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                             data = population,
#'                                             modelType = "cox")
#'   cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#'   approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "custom")
#'   return(approximation)
#' }
#' approximations <- lapply(populations, fitModelInDatabase)
#' approximations <- do.call("rbind", approximations)
#'
#' # At study coordinating center, perform meta-analysis using per-site approximations:
#' estimate <- computeBayesianMetaAnalysis(approximations)
#' plotPerDbMcmcTrace(estimate)
#'
#' @export
plotPerDbMcmcTrace <- function(estimate, showEstimate = TRUE, dataCutoff = 0.01, fileName = NULL) {
  traces <- attr(estimate, "traces")
  getDbChain <- function(i) {
    data <- data.frame(x = 1:nrow(traces), trace = traces[, i], var = sprintf("Site %s", i - 2))
    if (dataCutoff != 0) {
      lims <- quantile(data$trace, c(dataCutoff, 1 - dataCutoff))
      data <- data[data$trace > lims[1] & data$trace < lims[2], ]
    }
    return(data)
  }
  data <- lapply(3:ncol(traces), getDbChain)
  data <- do.call(rbind, data)

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$trace)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::scale_x_continuous("Iterations") +
    ggplot2::scale_y_continuous(expression(theta)) +
    ggplot2::facet_grid(var ~ ., scales = "free")
  if (showEstimate) {
    getDbMedian <- function(i) {
      return(data.frame(trace = median(traces[, i]), var = sprintf("Site %s", i - 2)))
    }
    mode <- lapply(3:ncol(traces), getDbMedian)
    mode <- do.call(rbind, mode)

    getDbCi <- function(i) {
      return(data.frame(trace = HDInterval::hdi(traces[, i]), var = sprintf("Site %s", i - 2)))
    }
    ci <- lapply(3:ncol(traces), getDbCi)
    ci <- do.call(rbind, ci)

    plot <- plot +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$trace), data = mode) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$trace),
                          data = ci,
                          linetype = "dashed")
  }
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  return(plot)
}

#' Plot posterior density
#'
#' @details
#' Plot the density of the posterior distribution of the mu and tau parameters.
#'
#' @param estimate       An object as generated using the [computeBayesianMetaAnalysis()] function.
#' @param showEstimate   Show the parameter estimates (mode) and 95 percent confidence intervals?
#' @param dataCutoff     This fraction of the data at both tails will be removed.
#' @param fileName       Name of the file where the plot should be saved, for example 'plot.png'. See
#'                       the function [ggplot2::ggsave] in the ggplot2 package for supported file
#'                       formats.
#'
#' @seealso
#' [computeBayesianMetaAnalysis]
#'
#' @return
#' A Ggplot object. Use the [ggplot2::ggsave] function to save to file.
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' # Fit a Cox regression at each data site, and approximate likelihood function:
#' fitModelInDatabase <- function(population) {
#'   cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                             data = population,
#'                                             modelType = "cox")
#'   cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#'   approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "custom")
#'   return(approximation)
#' }
#' approximations <- lapply(populations, fitModelInDatabase)
#' approximations <- do.call("rbind", approximations)
#'
#' # At study coordinating center, perform meta-analysis using per-site approximations:
#' estimate <- computeBayesianMetaAnalysis(approximations)
#' plotPosterior(estimate)
#'
#' @export
plotPosterior <- function(estimate, showEstimate = TRUE, dataCutoff = 0.01, fileName = NULL) {
  traces <- attr(estimate, "traces")
  dataMu <- data.frame(x = 1:nrow(traces), trace = traces[, 1], var = "Mu")
  dataTau <- data.frame(x = 1:nrow(traces), trace = traces[, 2], var = "Tau")
  if (dataCutoff != 0) {
    limsMu <- quantile(traces[, 1], c(dataCutoff, 1 - dataCutoff))
    limsTau <- quantile(traces[, 2], c(dataCutoff, 1 - dataCutoff))
    dataMu <- dataMu[dataMu$trace > limsMu[1] & dataMu$trace < limsMu[2], ]
    dataTau <- dataTau[dataTau$trace > limsTau[1] & dataTau$trace < limsTau[2], ]
  }
  data <- rbind(dataMu, dataTau)
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$trace)) +
    ggplot2::geom_density(alpha = 0.4, fill = rgb(0, 0, 0), color = rgb(0, 0, 0)) +
    ggplot2::scale_y_continuous("Density") +
    ggplot2::facet_grid(~var, scales = "free", switch = "both") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   strip.placement = "outside",
                   strip.background = ggplot2::element_blank())
  if (showEstimate) {
    mode <- data.frame(trace = c(estimate$mu, estimate$tau), var = c("Mu", "Tau"))
    ci <- data.frame(trace = c(estimate$mu95Lb,
                               estimate$mu95Ub,
                               estimate$tau95Lb,
                               estimate$tau95Ub), var = c("Mu", "Mu", "Tau", "Tau"))

    plot <- plot +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$trace), data = mode) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$trace),
                          data = ci,
                          linetype = "dashed")
  }
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  return(plot)
}

#' Plot posterior density per database
#'
#' @details
#' Plot the density of the posterior distribution of the theta parameter (the estimated log hazard
#' ratio) at each site.
#'
#' @param estimate       An object as generated using the [computeBayesianMetaAnalysis()] function.
#' @param showEstimate   Show the parameter estimates (mode) and 95 percent confidence intervals?
#' @param dataCutoff     This fraction of the data at both tails will be removed.
#' @param fileName       Name of the file where the plot should be saved, for example 'plot.png'. See
#'                       the function [ggplot2::ggsave] in the ggplot2 package for supported file
#'                       formats.
#'
#' @return
#' A Ggplot object. Use the [ggplot2::ggsave] function to save to file.
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' # Fit a Cox regression at each data site, and approximate likelihood function:
#' fitModelInDatabase <- function(population) {
#'   cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'                                             data = population,
#'                                             modelType = "cox")
#'   cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#'   approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "custom")
#'   return(approximation)
#' }
#' approximations <- lapply(populations, fitModelInDatabase)
#' approximations <- do.call("rbind", approximations)
#'
#' # At study coordinating center, perform meta-analysis using per-site approximations:
#' estimate <- computeBayesianMetaAnalysis(approximations)
#' plotPerDbPosterior(estimate)
#'
#' @export
plotPerDbPosterior <- function(estimate, showEstimate = TRUE, dataCutoff = 0.01, fileName = NULL) {
  traces <- attr(estimate, "traces")
  getDbChain <- function(i) {
    data <- data.frame(x = 1:nrow(traces), trace = traces[, i], var = sprintf("Site %s", i - 2))
    if (dataCutoff != 0) {
      lims <- quantile(data$trace, c(dataCutoff, 1 - dataCutoff))
      data <- data[data$trace > lims[1] & data$trace < lims[2], ]
    }
    return(data)
  }
  data <- lapply(3:ncol(traces), getDbChain)
  data <- do.call(rbind, data)

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$trace)) +
    ggplot2::geom_density(alpha = 0.4, fill = rgb(0, 0, 0), color = rgb(0, 0, 0)) +
    ggplot2::scale_x_continuous(expression(theta)) +
    ggplot2::scale_y_continuous("Density") +
    ggplot2::facet_grid(var ~ ., scales = "free")
  if (showEstimate) {
    getDbMedian <- function(i) {
      return(data.frame(trace = median(traces[, i]), var = sprintf("Site %s", i - 2)))
    }
    mode <- lapply(3:ncol(traces), getDbMedian)
    mode <- do.call(rbind, mode)

    getDbCi <- function(i) {
      return(data.frame(trace = HDInterval::hdi(traces[, i]), var = sprintf("Site %s", i - 2)))
    }
    ci <- lapply(3:ncol(traces), getDbCi)
    ci <- do.call(rbind, ci)

    plot <- plot +
      ggplot2::geom_vline(ggplot2::aes(xintercept = trace), data = mode) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = trace), data = ci, linetype = "dashed")
  }
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  return(plot)
}
