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

#' Create a forest plot
#'
#' @description
#' Creates a forest plot of effect size estimates, including the summary estimate.
#'
#' @details
#' Creates a forest plot of effect size estimates, including a meta-analysis estimate.
#'
#' @param data      A data frame containing either normal, skew-normal, custom parametric, or grid
#'                  likelihood data. One row per database.
#' @param labels    A vector of labels for the data sources.
#' @param estimate  The meta-analytic estimate as created using either ['computeFixedEffectMetaAnalysis()`]
#'                  or [`computeBayesianMetaAnalysis()`] function.
#' @param xLabel    The label on the x-axis: the name of the effect estimate.
#' @param summaryLabel The label for the meta-analytic estimate.
#' @param limits    The limits of the effect size axis.
#' @param alpha     The alpha (expected type I error).
#' @param showLikelihood Show the likelihood curve for each estimate?
#' @param fileName  Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                  function [ggplot2::ggsave] ifor supported file formats.
#'
#' @return
#' A Ggplot object. Use the [ggplot2::ggsave] function to save to file.
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#' labels <- paste("Data site", LETTERS[1:length(populations)])
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
#' plotMetaAnalysisForest(approximations, labels, estimate)
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
plotMetaAnalysisForest <- function(data,
                                   labels,
                                   estimate,
                                   xLabel = "Relative risk",
                                   summaryLabel = "Summary",
                                   limits = c(0.1, 10),
                                   alpha = 0.05,
                                   showLikelihood = TRUE,
                                   fileName = NULL) {
  if (is.data.frame(data)) {
    data <- split(data, seq_len(nrow(data)))
  }
  if (length(data) != length(labels)) {
    abort("The number of labels should be equal to the number of approximatins in `data`")
  }

  d1 <- data.frame(logRr = -100,
                   logLb95Ci = -100,
                   logUb95Ci = -100,
                   type = "header",
                   label = "Source")
  getEstimate <- function(approximation) {
    ci <- suppressMessages(computeConfidenceInterval(approximation = approximation,
                                                     alpha = alpha))
    return(data.frame(logRr = ci$logRr,
                      logLb95Ci = log(ci$lb),
                      logUb95Ci = log(ci$ub),
                      type = "db"))
  }
  d2 <- lapply(data, getEstimate)
  d2 <- do.call(rbind, d2)
  d2$label <- labels
  if ("rr" %in% colnames(estimate)) {
    d3 <- data.frame(logRr = estimate$logRr,
                     logLb95Ci = log(estimate$lb),
                     logUb95Ci = log(estimate$ub),
                     type = "ma",
                     label = summaryLabel)
  } else {
    d3 <- data.frame(logRr = estimate$mu,
                     logLb95Ci = estimate$mu95Lb,
                     logUb95Ci = estimate$mu95Ub,
                     type = "ma",
                     label = sprintf("%s (tau = %.2f)", summaryLabel, estimate$tau))
  }
  d <- rbind(d1, d2, d3)
  d$y <- seq(nrow(d), 1)
  # d$name <- factor(d$name, levels = c(d3$name, rev(as.character(labels)), "Source"))

  # ggplot puts whisker for infinite values, but not large values:
  plotD <- d
  plotD$logLb95Ci[is.infinite(plotD$logLb95Ci)] <- -10
  plotD$logUb95Ci[is.infinite(plotD$logUb95Ci)] <- 10

  rowHeight <- 0.8
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  yLimits <- c(min(d$y) - rowHeight / 2, max(d$y) + rowHeight / 2)
  p <- ggplot2::ggplot(plotD, ggplot2::aes(x = exp(.data$logRr), y = .data$y))

  if (showLikelihood) {
    yLimit <- -qchisq(1 - alpha, df = 1)/2
    transformY <- function(coords, yOffset) {
      coords$y <- coords$y - max(coords$y)
      coords <- coords[coords$y >= yLimit, ]
      coords$y <- rowHeight/2 - (coords$y * rowHeight / yLimit) + yOffset
      coords$y2 <- yOffset - rowHeight/2
      return(coords)
    }
    getCoords <- function(i) {
      coords <- getLikelihoodCoordinates(data[[i]], limits, verbose = FALSE)
      coords <- transformY(coords, length(data) - i + 2)
      coords$group <- length(data) - i + 2
      return(coords)
    }
    llData1 <- lapply(seq_along(data), getCoords)
    llData1 <- do.call(rbind, llData1)

    if ("rr" %in% colnames(estimate)) {
      llData2 <- getLikelihoodCoordinates(estimate, limits, verbose = FALSE)
      llData2 <- transformY(llData2, 1)
      llData2$group <- 1
    } else {
      trace <- attr(estimate, "traces")[, 1]
      trace <- trace[trace > log(limits[1]) & trace < log(limits[2])]
      coords <- density(trace, from = log(limits[1]), to = log(limits[2]), n = 100)
      llData2 <- data.frame(x = coords$x, y = log(coords$y))
      llData2 <- transformY(llData2, 1)
      llData2$group <- 1
    }

    llData <- rbind(llData1,llData2)
    uniqueYs <- unique(llData$y2)
    nSteps <- 10
    ys <- rep(uniqueYs,  nSteps) + rep((rowHeight/2) * ((seq_len(nSteps)-1) / nSteps), each = length(uniqueYs))
    fadeData <- data.frame(xmin = limits[1],
                           xmax = limits[2],
                           ymin = ys,
                           ymax = ys + rowHeight / 2 / nSteps,
                           alpha = 1 - rep(seq_len(nSteps) / nSteps, each = length(uniqueYs)),
                           logRr = 0,
                           y = 1)
    p <- p +  ggplot2::geom_ribbon(
      ggplot2::aes(x = exp(.data$x), ymax = .data$y, ymin = .data$y2, group = .data$group),
      fill = "#0000AA",
      alpha = 0.25,
      data = llData) +
      ggplot2::geom_rect(ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax, alpha = .data$alpha), fill = "#FFFFFF", data = fadeData)
  }
  p <- p +
    ggplot2::geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
    ggplot2::geom_vline(xintercept = 1, size = 0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = exp(.data$logLb95Ci),
                                         xmax = exp(.data$logUb95Ci)), height = 0.15) +
    ggplot2::geom_point(size = 3, shape = 23, ggplot2::aes(fill = .data$type)) +
    ggplot2::scale_fill_manual(values = c("#000000", "#000000", "#FFFFFF")) +
    ggplot2::scale_x_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) +
    ggplot2::coord_cartesian(xlim = limits, ylim = yLimits) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.border = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines"))

  p
  d$logLb95Ci[is.infinite(d$logLb95Ci)] <- NA
  d$logUb95Ci[is.infinite(d$logUb95Ci)] <- NA
  labels <- sprintf("%0.2f (%0.2f - %0.2f)", exp(d$logRr), exp(d$logLb95Ci), exp(d$logUb95Ci))
  labels <- gsub("NA", "", labels)
  labels <- data.frame(y = rep(d$y, 2),
                       x = rep(1:2, each = nrow(d)),
                       label = c(as.character(d$label), labels),
                       stringsAsFactors = FALSE)
  labels$label[nrow(d) + 1] <- paste(xLabel, "(95% CI)")
  data_table <- ggplot2::ggplot(labels, ggplot2::aes(x = .data$x,
                                                     y = .data$y,
                                                     label = .data$label)) +
    ggplot2::geom_text(size = 4, hjust = 0, vjust = 0.5) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = nrow(d) - 0.5)) +
    ggplot2::scale_y_continuous(limits = yLimits) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(colour = "white"),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_line(colour = "white"),
                   plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines")) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::coord_cartesian(xlim = c(1, 3))

  plot <- gridExtra::grid.arrange(data_table, p, ncol = 2)

  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 7, height = 1 + nrow(data) * 0.3, dpi = 400)
  invisible(plot)
}
