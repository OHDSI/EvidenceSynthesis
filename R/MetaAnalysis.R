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
#' @param estimate  The meta-analytic estimate as created using either [computeFixedEffectMetaAnalysis]
#'                  or [computeBayesianMetaAnalysis] function.
#' @param xLabel    The label on the x-axis: the name of the effect estimate.
#' @param summaryLabel The label for the meta-analytic estimate.
#' @param limits    The limits of the effect size axis.
#' @param alpha     The alpha (expected type I error).
#' @param showPredictionInterval Show the prediction interval (for random effects models).
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
#'     data = population,
#'     modelType = "cox"
#'   )
#'   cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#'   approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = "grid with gradients")
#'   return(approximation)
#' }
#' approximations <- lapply(populations, fitModelInDatabase)
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
                                   showPredictionInterval = TRUE,
                                   showLikelihood = TRUE,
                                   fileName = NULL) {
  if (is.numeric(data)) {
    data <- apply(data, 1, function(row) row, simplify = FALSE)
  }

  if (is.data.frame(data)) {
    data <- split(data, seq_len(nrow(data)))
  }

  if (length(data) != length(labels)) {
    abort("The number of labels should be equal to the number of approximations in `data`")
  }

  d1 <- tibble(
    logRr = -100,
    logLb95Ci = -100,
    logUb95Ci = -100,
    type = "header",
    label = "Source"
  )
  getEstimate <- function(approximation) {
    ci <- suppressMessages(computeConfidenceInterval(
      approximation = approximation,
      alpha = alpha
    ))
    return(tibble(
      logRr = ci$logRr,
      logLb95Ci = log(ci$lb),
      logUb95Ci = log(ci$ub),
      type = "db"
    ))
  }
  d2 <- lapply(data, getEstimate)
  d2 <- bind_rows(d2) |>
    mutate(label = labels)

  if ("rr" %in% colnames(estimate)) {
    # Estimate produced by computeFixedEffectMetaAnalysis
    d3 <- tibble(
      logRr = estimate$logRr,
      logLb95Ci = log(estimate$lb),
      logUb95Ci = log(estimate$ub),
      type = "ma",
      label = summaryLabel
    )
  } else if ("mu" %in% colnames(estimate)) {
    # Estimate produced by computeBayesianMetaAnalysis
    d3 <- tibble(
      logRr = c(estimate$logRr, NA),
      logLb95Ci = c(estimate$mu95Lb, NA),
      logUb95Ci = c(estimate$mu95Ub, NA),
      type = c("ma", "maSub"),
      label = c("Bayesian random effects", sprintf("\u03C4 = %.2f (%.2f - %.2f)", estimate$tau, estimate$tau95Lb, estimate$tau95Ub))
    )
    if (showPredictionInterval) {
      d3 <- bind_rows(
        d3,
        tibble(
          logRr = NA,
          logLb95Ci = estimate$predictionInterval95Lb,
          logUb95Ci = estimate$predictionInterval95Ub,
          type = "pi",
          label = "Prediction interval"
        )
      )
    }
  } else {
    stop("Unknown summary estimate type")
  }
  d <- bind_rows(d1, d2, d3)
  d <- d |>
    mutate(y = if_else(.data$type == "maSub", 0.5, 1)) |>
    mutate(y = sum(.data$y) - cumsum(.data$y) + 1)
  maEstimates <- d |>
    filter(.data$type == "ma", !is.na(.data$logRr))
  diamondData <- tibble(
    x = exp(c(maEstimates$logLb95Ci, maEstimates$logRr, maEstimates$logUb95Ci, maEstimates$logRr)),
    y = c(maEstimates$y, maEstimates$y + 0.2, maEstimates$y, maEstimates$y - 0.2),
    group = rep(maEstimates$y, 4)
  )

  maBoundaryY <- d |>
    filter(.data$type == "ma") |>
    summarise(max(.data$y)) |>
    pull() + 0.5

  # ggplot puts whisker for infinite values, but not large values:
  plotD <- d
  plotD$logLb95Ci[is.infinite(plotD$logLb95Ci)] <- -10
  plotD$logUb95Ci[is.infinite(plotD$logUb95Ci)] <- 10
  plotD <- plotD |>
    filter(.data$type != "ma")

  rowHeight <- 0.5
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  yLimits <- c(min(d$y) - rowHeight / 2, max(d$y) + rowHeight / 2)
  rightPlot <- ggplot2::ggplot(plotD, ggplot2::aes(x = exp(.data$logRr), y = .data$y)) +
    ggplot2::geom_rect(xmin = -10, xmax = 10, ymin = 0, ymax = maBoundaryY, size = 0, fill = "#69AED5", alpha = 0.25, data = tibble(logRr = 1, y = 1)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, y = .data$y, xend = .data$x, yend = .data$yend), color = "#AAAAAA",  size = 0.2, data = data.frame(x = breaks, y = 0, yend = max(d$y) - 0.5)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, y = .data$y, xend = .data$x, yend = .data$yend), size = 0.5, data = data.frame(x = 1, y = 0, yend = max(d$y) - 0.5)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = max(d$y) - 0.5)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = exp(.data$logLb95Ci), xmax = exp(.data$logUb95Ci)), height = 0.15) +
    ggplot2::geom_point(size = 3, shape = 16) +
    ggplot2::geom_polygon(ggplot2::aes(x = .data$x, y = .data$y, group = .data$group), data = diamondData) +
    ggplot2::scale_x_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) +
    ggplot2::coord_cartesian(xlim = limits, ylim = yLimits) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(),
      plot.margin = ggplot2::unit(c(0, 0, 0, -0.19), "lines")
    )
  rightPlot
  d$logLb95Ci[is.infinite(d$logLb95Ci)] <- NA
  d$logUb95Ci[is.infinite(d$logUb95Ci)] <- NA
  d$logRr[exp(d$logRr) < limits[1] | exp(d$logRr) > limits[2]] <- NA
  estimateLabels <- sprintf("%0.2f", exp(d$logRr))
  estimateLabels <- gsub("NA", "-", estimateLabels)
  estimateLabels[d$type %in% c("maSub", "pi", "header")] <- ""
  intervalLabels <- sprintf("(%0.2f - %0.2f)", exp(d$logLb95Ci), exp(d$logUb95Ci))
  intervalLabels <- gsub("NA", "", intervalLabels)
  intervalLabels <- gsub("\\( - \\)", "", intervalLabels)
  intervalLabels[d$type %in% c("maSub", "header")] <- ""

  textTable <- data.frame(
    y = rep(d$y, 3),
    x = rep(c(1, 2, 2.2), each = nrow(d)),
    label = c(as.character(d$label), estimateLabels, intervalLabels),
    fontface = rep(if_else(d$type %in% c("ma", "header", "pi"), "bold", "plain"), 3)
  )
  textTable$label[nrow(d) + 1] <- paste(xLabel, "(95% CI)")
  xLimits <- c(1, 2.75)
  widths <- c(1.5, 1)
  width <- 7
  leftPlot <- ggplot2::ggplot(textTable, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label)) +
    ggplot2::geom_rect(xmin = -10, xmax = 10, ymin = 0, ymax = maBoundaryY, size = 0, fill = "#69AED5", alpha = 0.25, data = tibble(x = 1, y = 1, label = "NA")) +
    ggplot2::geom_text(ggplot2::aes(fontface = .data$fontface), size = 4, hjust = 0, vjust = 0.5) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = max(d$y) - 0.5)) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::coord_cartesian(xlim = xLimits, ylim = yLimits) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "white"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(),
      plot.margin = ggplot2::unit(c(0, 0, 0, 0), "lines")
    )
  plot <- gridExtra::grid.arrange(leftPlot, rightPlot, ncol = 2, widths = widths, padding = ggplot2::unit(0, "line"))
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = width, height = 1 + max(d$y) * 0.3, dpi = 300)
  }
  invisible(plot)
}
