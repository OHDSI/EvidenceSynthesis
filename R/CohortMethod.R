# Copyright 2023 Observational Health Data Sciences and Informatics
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

#' Plot covariate balances
#'
#' @description
#' Plots the covariate balance before and after matching for multiple data sources.
#'
#' @details
#' Creates a plot showing the covariate balance before and after matching. Balance distributions are
#' displayed as box plots combined with scatterplots.
#'
#' @param balances      A list of covariate balance objects as created using the
#'                      `computeCovariateBalance()` function in the `CohortMethod` package. Each
#'                      balance object is expected to be a data frame with at least these two columns:
#'                      `beforeMatchingStdDiff` and `afterMatchingStdDiff`.
#' @param labels        A vector containing the labels for the various sources.
#' @param threshold     Show a threshold value for the standardized difference.
#' @param beforeLabel   Label for before matching / stratification / trimming.
#' @param afterLabel    Label for after matching / stratification / trimming.
#' @param fileName      Name of the file where the plot should be saved, for example 'plot.png'. See
#'                      the function [ggplot2::ggsave] for supported file formats.
#'
#' @return
#' A Ggplot object. Use the [ggplot2::ggsave].
#'
#' @examples
#' # Some example data:
#' balance1 <- data.frame(
#'   beforeMatchingStdDiff = rnorm(1000, 0.1, 0.1),
#'   afterMatchingStdDiff = rnorm(1000, 0, 0.01)
#' )
#' balance2 <- data.frame(
#'   beforeMatchingStdDiff = rnorm(1000, 0.2, 0.1),
#'   afterMatchingStdDiff = rnorm(1000, 0, 0.05)
#' )
#' balance3 <- data.frame(
#'   beforeMatchingStdDiff = rnorm(1000, 0, 0.1),
#'   afterMatchingStdDiff = rnorm(1000, 0, 0.03)
#' )
#' plotCovariateBalances(
#'   balances = list(balance1, balance2, balance3),
#'   labels = c("Site A", "Site B", "Site C")
#' )
#'
#' @export
plotCovariateBalances <- function(balances,
                                  labels,
                                  threshold = 0,
                                  beforeLabel = "Before matching",
                                  afterLabel = "After matching",
                                  fileName = NULL) {
  if (length(balances) != length(labels)) {
    abort("The balances and labels arguments must have the same length")
  }
  prepare <- function(i) {
    before <- data.frame(
      stdDiff = balances[[i]]$beforeMatchingStdDiff,
      type = beforeLabel,
      label = labels[i],
      stringsAsFactors = FALSE
    )
    after <- data.frame(
      stdDiff = balances[[i]]$afterMatchingStdDiff,
      type = afterLabel,
      label = labels[i],
      stringsAsFactors = FALSE
    )
    return(rbind(before, after))
  }
  balances <- lapply(1:length(balances), prepare)
  balance <- do.call("rbind", balances)
  labelToY <- data.frame(label = labels, y = 1 + length(labels) - (1:length(labels)))
  balance <- merge(balance, labelToY)
  balance$yDither <- balance$y + runif(nrow(balance), -0.3, 0.3)
  combis <- unique(balance[, c("type", "label", "y")])
  prepareAgg <- function(i) {
    subset <- balance[balance$type == combis$type[i] & balance$label == combis$label[i], ]
    result <- quantile(subset$stdDiff, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    result <- data.frame(
      type = combis$type[i],
      label = combis$label[i],
      y = combis$y[i],
      ymin = result[1],
      lower = result[2],
      middle = result[3],
      upper = result[4],
      ymax = result[5],
      count = nrow(subset)
    )
    return(result)
  }
  agg <- lapply(1:nrow(combis), prepareAgg)
  agg <- do.call("rbind", agg)
  agg$type <- factor(agg$type, levels = c(beforeLabel, afterLabel))
  agg$label <- factor(agg$label, levels = labels)
  balance$type <- factor(balance$type, levels = c(beforeLabel, afterLabel))
  balance$label <- factor(balance$label, levels = labels)
  limits <- c(min(balance$stdDiff, na.rm = TRUE), max(balance$stdDiff, na.rm = TRUE))
  plot <- ggplot2::ggplot(balance, ggplot2::aes(x = .data$y, group = .data$label)) +
    ggplot2::geom_point(ggplot2::aes(x = .data$yDither, y = .data$stdDiff),
      color = rgb(0, 0, 0.8, alpha = 0.3),
      shape = 16,
      size = 0.5
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ymin, ymax = .data$ymin), data = agg) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ymax, ymax = .data$ymax), data = agg) +
    ggplot2::geom_boxplot(ggplot2::aes(
      ymin = .data$ymin,
      lower = .data$lower,
      middle = .data$middle,
      upper = .data$upper,
      ymax = .data$ymax
    ), outlier.shape = NA, data = agg, stat = "identity", alpha = 0.65) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_x_continuous(limits = c(0.5, length(unique(agg$label)) + 1.75)) +
    ggplot2::scale_y_continuous("Standardized difference of mean", limits = limits) +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(~type) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#AAAAAA"),
      panel.background = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 11),
      axis.title.x = ggplot2::element_text(size = 11),
      axis.ticks.x = ggplot2::element_line(color = "#AAAAAA"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 11),
      plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines")
    )

  if (threshold != 0) {
    plot <- plot + ggplot2::geom_hline(yintercept = c(threshold, -threshold), linetype = "dotted")
  }
  after <- agg[agg$type == afterLabel, ]
  after$max <- pmax(abs(after$ymin), abs(after$ymax))
  text <- data.frame(
    y = rep(c(after$y, nrow(after) + 1.25), 3),
    x = rep(c(1, 2, 3), each = nrow(after) + 1),
    label = c(c(
      as.character(after$label),
      "Source",
      formatC(after$count,
        big.mark = ",",
        format = "d"
      ),
      "Covariate\ncount",
      formatC(after$max,
        digits = 2,
        format = "f"
      ),
      paste(afterLabel,
        "max(absolute)",
        sep = "\n"
      )
    )), dummy = ""
  )

  data_table <- ggplot2::ggplot(text, ggplot2::aes(
    x = .data$x,
    y = .data$y,
    label = .data$label
  )) +
    ggplot2::geom_text(size = 4, hjust = 0, vjust = 0.5) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = nrow(after) + 0.5)) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = "white"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = "white"),
      strip.background = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines")
    ) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::facet_grid(~dummy) +
    ggplot2::coord_cartesian(xlim = c(1, 4), ylim = c(0.5, length(unique(agg$label)) + 1.75))

  plot <- gridExtra::grid.arrange(data_table, plot, ncol = 2)

  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 8, height = 1 + length(balances) * 0.4, dpi = 400)
  }
  invisible(plot)
}


computePreferenceScore <- function(data, unfilteredData = NULL) {
  if (is.null(unfilteredData)) {
    proportion <- sum(data$treatment) / nrow(data)
  } else {
    proportion <- sum(unfilteredData$treatment) / nrow(unfilteredData)
  }
  propensityScore <- data$propensityScore
  propensityScore[propensityScore > 0.9999999] <- 0.9999999
  x <- exp(log(propensityScore / (1 - propensityScore)) - log(proportion / (1 - proportion)))
  data$preferenceScore <- x / (x + 1)
  return(data)
}

#' Prepare to plot the propensity score distribution
#'
#' @description
#' Prepare to plot the propensity (or preference) score distribution. It computes the distribution, so
#' the output does not contain person-level data.
#'
#' @param data             A data frame with at least the two columns described below
#' @param unfilteredData   To be used when computing preference scores on data from which subjects have
#'                         already been removed, e.g. through trimming and/or matching. This data frame
#'                         should have the same structure as `data`.
#' @param scale            The scale of the graph. Two scales are supported: `scale = 'propensity'` or
#'                         `scale = 'preference'`. The preference score scale is defined by Walker et
#'                         al. (2013).
#'
#' @details
#' The data frame should have a least the following two columns:
#' - **treatment** (integer): Column indicating whether the person is in the treated (1) or comparator
#' (0) group. - **propensityScore** (numeric): Propensity score.
#'
#' @seealso
#' [plotPreparedPs]
#'
#' @return
#' A data frame describing the propensity score (or preference score) distribution at 100
#' equally-spaced points.
#'
#' @examples
#' # Simulate some data for this example:
#' treatment <- rep(0:1, each = 100)
#' propensityScore <- c(rnorm(100, mean = 0.4, sd = 0.25), rnorm(100, mean = 0.6, sd = 0.25))
#' data <- data.frame(treatment = treatment, propensityScore = propensityScore)
#' data <- data[data$propensityScore > 0 & data$propensityScore < 1, ]
#'
#' preparedPlot <- preparePsPlot(data)
#'
#' @references
#' Walker AM, Patrick AR, Lauer MS, Hornbrook MC, Marin MG, Platt R, Roger VL, Stang P, and
#' Schneeweiss S. (2013) A tool for assessing the feasibility of comparative effectiveness research,
#' Comparative Effective Research, 3, 11-20
#'
#' @export
preparePsPlot <- function(data, unfilteredData = NULL, scale = "preference") {
  if (!("treatment" %in% colnames(data))) {
    abort("Missing column treatment in data")
  }
  if (!("propensityScore" %in% colnames(data))) {
    abort("Missing column propensityScore in data")
  }
  if (!is.null(unfilteredData)) {
    if (!("treatment" %in% colnames(unfilteredData))) {
      abort("Missing column treatment in unfilteredData")
    }
    if (!("propensityScore" %in% colnames(unfilteredData))) {
      abort("Missing column propensityScore in unfilteredData")
    }
  }
  if (scale != "propensity" && scale != "preference") {
    abort(paste0("Unknown scale '", scale, "', please choose either 'propensity' or 'preference'"))
  }

  if (scale == "preference") {
    data <- computePreferenceScore(data, unfilteredData)
    d1 <- density(data$preferenceScore[data$treatment == 1], from = 0, to = 1, n = 100)
    d0 <- density(data$preferenceScore[data$treatment == 0], from = 0, to = 1, n = 100)
    d <- data.frame(
      preferenceScore = c(d1$x, d0$x),
      y = c(d1$y, d0$y),
      treatment = c(rep(1, length(d1$x)), rep(0, length(d0$x)))
    )
    d$y <- d$y / max(d$y)
  } else {
    d1 <- density(data$propensityScore[data$treatment == 1], from = 0, to = 1, n = 100)
    d0 <- density(data$propensityScore[data$treatment == 0], from = 0, to = 1, n = 100)
    d <- data.frame(
      propensityScore = c(d1$x, d0$x),
      y = c(d1$y, d0$y),
      treatment = c(rep(1, length(d1$x)), rep(0, length(d0$x)))
    )
    d$y <- d$y / max(d$y)
  }
  return(d)
}

#' Plot the propensity score distribution
#'
#' @param preparedPsPlots   list of prepared propensity score data as created by the [preparePsPlot()]
#'                          function.
#' @param labels            A vector containing the labels for the various sources.
#' @param treatmentLabel    A label to us for the treated cohort.
#' @param comparatorLabel   A label to us for the comparator cohort.
#' @param fileName          Name of the file where the plot should be saved, for example 'plot.png'.
#'                          See the function [ggplot2::ggsave] for supported file formats.
#'
#' @examples
#' # Simulate some data for this example:
#' treatment <- rep(0:1, each = 100)
#' propensityScore <- c(rnorm(100, mean = 0.4, sd = 0.25), rnorm(100, mean = 0.6, sd = 0.25))
#' data <- data.frame(treatment = treatment, propensityScore = propensityScore)
#' data <- data[data$propensityScore > 0 & data$propensityScore < 1, ]
#' preparedPlot <- preparePsPlot(data)
#'
#' # Just reusing the same data three times for demonstration purposes:
#' preparedPsPlots <- list(preparedPlot, preparedPlot, preparedPlot)
#' labels <- c("Data site A", "Data site B", "Data site C")
#'
#' plotPreparedPs(preparedPsPlots, labels)
#'
#' @seealso
#' [preparePsPlot]
#'
#' @return
#' A ggplot object. Use the [ggplot2::ggsave] function to save to file in a different format.
#'
#' @export
plotPreparedPs <- function(preparedPsPlots,
                           labels,
                           treatmentLabel = "Target",
                           comparatorLabel = "Comparator",
                           fileName = NULL) {
  if (length(preparedPsPlots) != length(labels)) {
    abort("The preparedPsPlots and labels arguments must have the same length")
  }

  for (i in 1:length(preparedPsPlots)) {
    preparedPsPlots[[i]]$label <- labels[i]
  }
  data <- do.call("rbind", preparedPsPlots)
  data$GROUP <- treatmentLabel
  data$GROUP[data$treatment == 0] <- comparatorLabel
  data$GROUP <- factor(data$GROUP, levels = c(treatmentLabel, comparatorLabel))
  if (!is.null(data$propensityScore)) {
    xLabel <- "Propensity score"
    data$x <- data$propensityScore
  } else {
    xLabel <- "Preference score"
    data$x <- data$preferenceScore
  }
  data$label <- factor(data$label, levels = labels)

  plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data$x,
    y = .data$y,
    color = .data$GROUP,
    group = .data$GROUP,
    fill = .data$GROUP
  )) +
    ggplot2::geom_density(stat = "identity") +
    ggplot2::scale_fill_manual(values = c(
      rgb(0.8, 0, 0, alpha = 0.5),
      rgb(0, 0, 0.8, alpha = 0.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      rgb(0.8, 0, 0, alpha = 0.5),
      rgb(0, 0, 0.8, alpha = 0.5)
    )) +
    ggplot2::scale_x_continuous(xLabel, limits = c(0, 1)) +
    ggplot2::scale_y_continuous("Density") +
    ggplot2::facet_grid(label ~ ., switch = "both") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = "top",
      strip.text.y = ggplot2::element_text(angle = 180, size = 11),
      strip.background = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 11),
      axis.title.x = ggplot2::element_text(size = 11),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = rgb(0.7, 0.7, 0.7))
    )
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName,
      plot,
      width = 4,
      height = 1 + length(preparedPsPlots) * 0.75,
      dpi = 400
    )
  }
  invisible(plot)
}
