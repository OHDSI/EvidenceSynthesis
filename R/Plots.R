# Copyright 2025 Observational Health Data Sciences and Informatics
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
#'   data = population,
#'   modelType = "cox"
#' )
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
  coords <- getLikelihoodCoordinates(approximation, limits)
  coords$ll <- getLikelihoodProfile(cyclopsFit, parameter, coords$x)
  coords <- coords[!is.nan(coords$ll), ]
  coords$ll <- coords$ll - max(coords$ll)
  plotData <- rbind(
    data.frame(x = coords$x, ll = coords$ll, type = "Likelihood"),
    data.frame(x = coords$x, ll = coords$y, type = "Approximation")
  )
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
    ggplot2::geom_line(ggplot2::aes(
      group = .data$type,
      color = .data$type,
      linetype = .data$type,
      size = .data$type
    ), alpha = 0.7) +
    ggplot2::scale_size_manual(values = c(1, 2, 2, 2) * 0.6) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) +
    ggplot2::scale_color_manual(values = c("#000000", "#66c2a5", "#fc8d62", "#8da0cb")) +
    ggplot2::scale_y_continuous(yLabel) +
    ggplot2::scale_x_continuous(xLabel,
                                limits = log(limits),
                                breaks = log(breaks),
                                labels = breaks
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "top"
    )
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 5, height = 5, dpi = 400)
  }
  return(plot)
}

getLikelihoodCoordinates <- function(approximation, limits, verbose = TRUE) {
  type <- detectApproximationType(approximation, verbose = verbose)
  if (type == "normal") {
    x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
    y <- dnorm(x, mean = approximation$logRr, sd = approximation$seLogRr, log = TRUE)
  } else if (type == "custom") {
    x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
    y <- customFunction(x,
                        mu = approximation$mu,
                        sigma = approximation$sigma,
                        gamma = approximation$gamma
    )
  } else if (type == "skew normal") {
    x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
    y <- skewNormal(x,
                    mu = approximation$mu,
                    sigma = approximation$sigma,
                    alpha = approximation$alpha
    )
  } else if (type == "adaptive grid") {
    x <- approximation$point
    y <- approximation$value
  } else if (type == "grid") {
    x <- as.numeric(names(approximation))
    y <- approximation
  } else if (type == "grid with gradients") {
    x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
    y <- hermiteInterpolation(x,
                              profile = approximation)
  } else {
    abort(sprintf("Approximation type '%s' not supported by this function", type))
  }
  y <- y - max(y)
  return(data.frame(x = x, y = y))
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
#'     data = population,
#'     modelType = "cox"
#'   )
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
  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data$x,
    y = .data$trace
  )) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::scale_x_continuous("Iterations") +
    ggplot2::facet_grid(var ~ ., scales = "free", switch = "both") +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank()
    )
  if (showEstimate) {
    mode <- data.frame(trace = c(estimate$mu, estimate$tau), var = c("Mu", "Tau"))
    ci <- data.frame(trace = c(
      estimate$mu95Lb,
      estimate$mu95Ub,
      estimate$tau95Lb,
      estimate$tau95Ub
    ), var = c("Mu", "Mu", "Tau", "Tau"))

    plot <- plot +
      ggplot2::geom_hline(ggplot2::aes(yintercept = trace), data = mode) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = trace), data = ci, linetype = "dashed")
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }
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
#'     data = population,
#'     modelType = "cox"
#'   )
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
                          linetype = "dashed"
      )
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }
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
#'     data = population,
#'     modelType = "cox"
#'   )
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
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.background = ggplot2::element_blank()
    )
  if (showEstimate) {
    mode <- data.frame(trace = c(estimate$mu, estimate$tau), var = c("Mu", "Tau"))
    ci <- data.frame(trace = c(
      estimate$mu95Lb,
      estimate$mu95Ub,
      estimate$tau95Lb,
      estimate$tau95Ub
    ), var = c("Mu", "Mu", "Tau", "Tau"))

    plot <- plot +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$trace), data = mode) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data$trace),
                          data = ci,
                          linetype = "dashed"
      )
  }
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }
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
#'     data = population,
#'     modelType = "cox"
#'   )
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
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }
  return(plot)
}


#' Plot bias distributions
#'
#' @details
#' Plot empirical bias distributions learned from analyzing negative controls.
#'
#' @param biasDist       A bias distribution object generated by the [fitBiasDistribution()] or [sequentialFitBiasDistribution()] function.
#' @param limits         The lower and upper limits in log-RR to plot.
#' @param logScale       Whether or not to show bias in log-RR; default FALSE (shown in RR).
#' @param numericId      (For sequential or group case only) whether or not to treat `Id` as a numeric variable; default: TRUE.
#' @param fileName       Name of the file where the plot should be saved, for example 'plot.png'. See
#'                       the function [ggplot2::ggsave] in the ggplot2 package for supported file
#'                       formats.
#'
#' @seealso
#' [fitBiasDistribution], [sequentialFitBiasDistribution]
#'
#' @return
#' A `ggplot` object. Use the [ggplot2::ggsave] function to save to file.
#'
#' @examples
#' # Fit a bias distribution for this example:
#' data("ncLikelihoods")
#' # NOT RUN
#' # singleBiasDist = fitBiasDistribution(ncLikelihoods[[5]], seed = 1)
#'
#' # Plot it
#' # NOT RUN
#' # plotBiasDistribution(singleBiasDist)
#'
#' @export
plotBiasDistribution <- function(biasDist,
                                 limits = c(-2, 2),
                                 logScale = FALSE,
                                 numericId = TRUE,
                                 fileName = NULL) {
  # basic settings
  ypadding <- 0.03

  xlims <- limits
  if (logScale) {
    xbreaks <- c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3)
    xlabels <- as.character(xbreaks)
    xtext <- "Log rate ratio"
  } else {
    RRbreaks <- c(0.05, 0.1, 0.5, 1, 2, 5, 20, 50)
    xbreaks <- log(RRbreaks)
    xlabels <- as.character(RRbreaks)
    xtext <- "Rate ratio"
  }

  if ("Id" %in% names(biasDist)) {
    # plot sequential/group bias distributions

    if (numericId) {
      if (any(gsub("[0-9.-]", "", biasDist$Id) != "")) {
        stop("Id column cannot be converted to numerics. Try with `numericId = FALSE` instead?")
      }

      biasDist$Id <- as.integer(biasDist$Id)
    }

    medianDat <- data.frame(
      x = sapply(
        unique(biasDist$Id),
        function(id) median(biasDist[biasDist$Id == id, "bias"])
      ),
      Id = unique(biasDist$Id)
    )

    medianDat$y <- ypadding

    plot <- ggplot2::ggplot(biasDist, ggplot2::aes(x = .data$bias)) +
      ggplot2::geom_vline(
        xintercept = 0, linetype = 2,
        linewidth = 1, color = "gray40"
      ) +
      ggplot2::geom_density(fill = "gray80") +
      ggplot2::geom_point(
        data = medianDat,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape = 4
      ) +
      ggplot2::scale_x_continuous(
        limits = xlims,
        breaks = xbreaks,
        labels = xlabels
      ) +
      ggplot2::labs(
        x = sprintf("%s estimates for negative controls", xtext),
        y = "",
        caption = "Analysis period"
      ) +
      ggplot2::scale_y_continuous(breaks = NULL) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::facet_grid(. ~ Id, switch = "both") +
      ggplot2::theme(
        strip.placement = "inside",
        strip.text.y.left = ggplot2::element_text(angle = 0),
        panel.border = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(size = 14, hjust = 0.5)
      )
  } else {
    # plot only one bias distribution
    plot <- ggplot2::ggplot(biasDist, ggplot2::aes(x = .data$bias)) +
      ggplot2::geom_vline(
        xintercept = 0, linetype = 2,
        linewidth = 1, color = "gray40"
      ) +
      ggplot2::geom_density(fill = "gray80") +
      ggplot2::scale_x_continuous(
        limits = xlims,
        breaks = xbreaks,
        labels = xlabels
      ) +
      ggplot2::labs(
        x = sprintf("%s estimates for negative controls", xtext),
        y = ""
      ) +
      ggplot2::scale_y_continuous(breaks = NULL) +
      ggplot2::theme_bw(base_size = 14)
  }

  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 6, height = 4.5, dpi = 400)
  }

  return(plot)
}

#' Plot bias correction inference
#'
#' @details
#' Plot empirical bias distributions learned from analyzing negative controls.
#'
#' @param bbcResult      A (sequential) analysis object generated by the [biasCorrectionInference()] function.
#' @param type           The type of plot. Must be one of `c("corrected", "raw", "compare")`.
#' @param ids            IDs of the periods/groups to plot result for; default is all IDs.
#' @param limits         The limits on log RR for plotting.
#' @param logScale       Whether or not to show bias in log-RR; default FALSE (shown in RR).
#' @param numericId      Whether or not to treat `Id` as a numeric variable; default: TRUE.
#' @param fileName       Name of the file where the plot should be saved, for example 'plot.png'. See
#'                       the function [ggplot2::ggsave] in the ggplot2 package for supported file
#'                       formats.
#'
#' @seealso
#' [biasCorrectionInference]
#'
#' @return
#' A `ggplot` object. Use the [ggplot2::ggsave] function to save to file.
#'
#' @examples
#' # Perform sequential analysis using Bayesian bias correction for this example:
#' data("ncLikelihoods")
#' data("ooiLikelihoods")
#' # NOT RUN
#' # bbcSequential = biasCorrectionInference(ooiLikelihoods, ncLikelihoodProfiles = ncLikelihoods)
#'
#' # Plot it
#' # NOT RUN
#' # plotBiasCorrectionInference(bbcSequential, type = "corrected")
#'
#' @export
plotBiasCorrectionInference <- function(bbcResult,
                                        type = "raw",
                                        ids = bbcResult$Id,
                                        limits = c(-3, 3),
                                        logScale = FALSE,
                                        numericId = TRUE,
                                        fileName = NULL) {
  if (!type %in% c("corrected", "raw", "compare")) {
    stop("Unsupported plotting type! Make sure `type` is one of `c(\"corrected\", \"raw\", \"compare\")`.")
  }

  if (!attr(bbcResult, "corrected") && type %in% c("corrected", "compare")) {
    stop("Results are not bias corrected! Cannot plot the corrected results or compare.")
  }

  # basic settings
  ylims <- limits
  if (logScale) {
    ybreaks <- c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3)
    ylabels <- as.character(ybreaks)
    ytext <- "Log rate ratio"
  } else {
    RRbreaks <- c(0.05, 0.1, 0.5, 1, 2, 5, 20)
    ybreaks <- log(RRbreaks)
    ylabels <- as.character(RRbreaks)
    ytext <- "Rate ratio"
  }

  if (type == "compare") {
    dat <- dplyr::bind_rows(bbcResult, attr(bbcResult, "summaryRaw"))
    dat$biasCorrected <- rep(c("Yes", "No"), each = nrow(bbcResult))
    dat$biasCorrected <- factor(dat$biasCorrected, levels = c("Yes", "No"))

    dat <- dat[dat$Id %in% ids, ]

    if (numericId) {
      if (any(gsub("[0-9.-]", "", dat$Id) != "")) {
        stop("Id column cannot be converted to numerics. Try with `numericId = FALSE` instead?")
      }

      dat$Id <- as.integer(dat$Id)
      xbreaks <- seq(from = min(dat$Id), to = max(dat$Id), by = 2)
    }

    plot <- ggplot2::ggplot(dat, ggplot2::aes(
      x = .data$Id,
      y = .data$median,
      color = .data$biasCorrected
    )) +
      ggplot2::geom_hline(
        yintercept = 0, color = "gray50",
        linewidth = 1, linetype = 2
      ) +
      ggplot2::geom_pointrange(
        ggplot2::aes(
          ymin = .data$ci95Lb,
          ymax = .data$ci95Ub
        ),
        position = ggplot2::position_dodge(width = 0.3),
        size = 0.5, linewidth = 1.2
      ) +
      ggplot2::labs(
        x = "Sequential analysis period",
        y = sprintf("%s estimate (95%% CI)", ytext),
        color = "Bias correction?"
      ) +
      ggplot2::scale_y_continuous(
        limits = ylims, breaks = ybreaks,
        labels = ylabels
      ) +
      ggplot2::scale_color_manual(values = c("#046C9A", "#D69C4E")) +
      ggplot2::theme_bw(base_size = 14)

    if (numericId) {
      plot <- plot + ggplot2::scale_x_continuous(breaks = xbreaks)
    }
  } else {
    if (type == "corrected") {
      dat <- attr(bbcResult, "samples")
      capt <- "(with bias correction)"

      datProbs <- dplyr::bind_rows(
        data.frame(prob = bbcResult$p1, Id = bbcResult$Id, label = "P(H1)"),
        data.frame(prob = 1 - bbcResult$p1, Id = bbcResult$Id, label = "P(H0)")
      )
    } else {
      dat <- attr(bbcResult, "samplesRaw")
      capt <- "(without bias correction)"

      resultRaw <- attr(bbcResult, "summaryRaw")
      datProbs <- dplyr::bind_rows(
        data.frame(prob = resultRaw$p1, Id = resultRaw$Id, label = "P(H1)"),
        data.frame(prob = 1 - resultRaw$p1, Id = resultRaw$Id, label = "P(H0)")
      )
    }

    dat <- dat[dat$Id %in% ids, ]
    datProbs <- datProbs[datProbs$Id %in% ids, ]

    if (numericId) {
      if (any(gsub("[0-9.-]", "", dat$Id) != "")) {
        stop("Id column cannot be converted to numerics. Try with `numericId = FALSE` instead?")
      }

      dat$Id <- as.integer(dat$Id)
      datProbs$Id <- as.integer(datProbs$Id)
      xbreaks <- seq(from = min(datProbs$Id), to = max(datProbs$Id), by = 1)
    }

    medianDat <- data.frame(
      x = sapply(
        unique(dat$Id),
        function(id) median(dat[dat$Id == id, "beta"])
      ),
      Id = unique(dat$Id)
    )
    medianDat$y <- 0.05

    plot1 <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$beta)) +
      ggplot2::geom_vline(
        xintercept = 0, linetype = 2,
        linewidth = 1, color = "gray40"
      ) +
      ggplot2::geom_density(fill = "gray80") +
      ggplot2::geom_point(
        data = medianDat,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape = 4
      ) +
      ggplot2::scale_x_continuous(
        limits = ylims,
        breaks = ybreaks,
        labels = ylabels
      ) +
      ggplot2::labs(
        x = sprintf("%s", ytext),
        y = ""
      ) +
      ggplot2::scale_y_continuous(breaks = NULL) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::facet_grid(. ~ Id, switch = "both") +
      ggplot2::theme(
        strip.placement = "inside",
        strip.text.y.left = ggplot2::element_text(angle = 0),
        panel.border = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(size = 14, hjust = 0.5)
      )

    plot2 <- ggplot2::ggplot(
      datProbs,
      ggplot2::aes(
        x = .data$Id,
        y = .data$prob,
        color = .data$label
      )
    ) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        y = "Posterior probability",
        x = sprintf("Analysis period %s", capt),
        color = ""
      ) +
      ggplot2::scale_color_manual(values = c("#046C9A", "#D69C4E")) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(legend.position = "bottom")

    if (numericId) {
      plot2 <- plot2 + ggplot2::scale_x_continuous(
        breaks = xbreaks,
        expand = ggplot2::expansion(add = 0.3)
      )
    }

    plot <- gridExtra::grid.arrange(plot1, plot2, nrow = 2, heights = c(2.5, 1.5))
  }

  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 7.5, height = 7, dpi = 400)
  }

  return(plot)
}
