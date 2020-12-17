# Copyright 2020 Observational Health Data Sciences and Informatics
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

#' Plot empirical null distributions
#'
#' @description
#' Plot the empirical null distribution for multiple data sources.
#'
#' @details
#' Creates a plot showing the empirical null distributions. Distributions are shown as mean plus minus
#' one standard deviation, as well as a distribution plot.
#'
#' @param logRr      A numeric vector of effect estimates for the negative controls on the log scale.
#' @param seLogRr    The standard error of the log of the effect estimates. Hint: often the standard
#'                   error = (log(lower bound 95 percent confidence interval) - l og(effect
#'                   estimate))/qnorm(0.025).
#' @param labels     A vector containing the labels for the various sources. Should be of equal length
#'                   as `logRr` and `seLogRr`.
#' @param xLabel     The label on the x-axis: the name of the effect estimate.
#' @param limits     The limits of the effect size axis.
#' @param showCis    Show the 95 percent confidence intervals on the null distribution and distribution
#'                   parameter estimates?
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See the
#'                   function [ggplot2::ggsave()] for supported file formats.
#'
#' @seealso
#' [EmpiricalCalibration::fitNull], [EmpiricalCalibration::fitMcmcNull]
#'
#' @return
#' A Ggplot object. Use the [ggplot2::ggsave()] function to save to file.
#'
#' @examples
#' # Some example data:
#' site1 <- EmpiricalCalibration::simulateControls(n = 50, mean = 0, sd = 0.1, trueLogRr = 0)
#' site1$label <- "Site 1"
#' site2 <- EmpiricalCalibration::simulateControls(n = 50, mean = 0.1, sd = 0.2, trueLogRr = 0)
#' site2$label <- "Site 2"
#' site3 <- EmpiricalCalibration::simulateControls(n = 50, mean = 0.15, sd = 0.25, trueLogRr = 0)
#' site3$label <- "Site 3"
#' sites <- rbind(site1, site2, site3)
#'
#' plotEmpiricalNulls(logRr = sites$logRr, seLogRr = sites$seLogRr, labels = sites$label)
#'
#' @export
plotEmpiricalNulls <- function(logRr,
                               seLogRr,
                               labels,
                               xLabel = "Relative risk",
                               limits = c(0.1, 10),
                               showCis = TRUE,
                               fileName = NULL) {
  if (length(logRr) != length(seLogRr))
    abort("The logRr and seLogRr arguments must have the same length")
  if (length(logRr) != length(labels))
    abort("The logRr and labels arguments must have the same length")
  
  d <- data.frame(label = unique(labels),
                  mean = NA,
                  sd = NA,
                  xMin = NA,
                  xMax = NA,
                  meanLabel = "",
                  sdLabel = "",
                  y = length(unique(labels)) - (1:length(unique(labels))) + 1,
                  stringsAsFactors = FALSE)
  dist <- data.frame(label = rep(unique(labels), each = 100),
                     x = seq(log(limits[1]), log(limits[2]), length.out = 100),
                     yMax = NA,
                     yMaxUb = NA,
                     yMaxLb = NA,
                     yMin = NA)
  if (isRmdCheck()) {
    iter <- 1000
  } else {
    iter <- 10000
  }
  for (i in 1:nrow(d)) {
    idx <- labels == d$label[i]
    if (showCis) {
      null <- EmpiricalCalibration::fitMcmcNull(logRr = logRr[idx], seLogRr = seLogRr[idx], iter = iter)
      mcmc <- attr(null, "mcmc")
      lb95Mean <- quantile(mcmc$chain[, 1], 0.025)
      ub95Mean <- quantile(mcmc$chain[, 1], 0.975)
      ub95Sd <- 1/sqrt(quantile(mcmc$chain[, 2], 0.025))
      lb95Sd <- 1/sqrt(quantile(mcmc$chain[, 2], 0.975))
      d$mean[i] <- null[1]
      d$sd[i] <- 1/sqrt(null[2])
      d$xMin[i] <- d$mean[i] - d$sd[i]
      d$xMax[i] <- d$mean[i] + d$sd[i]
      d$meanLabel[i] <- sprintf("% 1.2f (% 1.2f - % 1.2f)", d$mean[i], lb95Mean, ub95Mean)
      d$sdLabel[i] <- sprintf("%1.2f (%1.2f - %1.2f)", d$sd[i], lb95Sd, ub95Sd)
      idx <- dist$label == d$label[i]
      compute <- function(x) {
        yMcmc <- dnorm(rep(x,
                           nrow(mcmc$chain)), mean = mcmc$chain[, 1], sd = 1/sqrt(mcmc$chain[, 2]))
        return(quantile(yMcmc, c(0.025, 0.5, 0.975)))
      }
      ys <- sapply(dist$x[idx], compute)
      y <- ys[2, ]
      yMaxLb <- ys[1, ]
      yMaxUb <- ys[3, ]
      normFactor <- max(ys[2, ])
      yMaxUb[yMaxUb > normFactor] <- normFactor
      y <- 0.7 * y/normFactor
      yMaxLb <- 0.7 * yMaxLb/normFactor
      yMaxUb <- 0.7 * yMaxUb/normFactor
      dist$yMax[idx] <- d$y[i] - 0.35 + y
      dist$yMaxLb[idx] <- d$y[i] - 0.35 + yMaxLb
      dist$yMaxUb[idx] <- d$y[i] - 0.35 + yMaxUb
      dist$yMin[idx] <- d$y[i] - 0.35
    } else {
      null <- EmpiricalCalibration::fitNull(logRr = logRr[idx], seLogRr = seLogRr[idx])
      d$mean[i] <- null[1]
      d$sd[i] <- null[2]
      d$xMin[i] <- null[1] - null[2]
      d$xMax[i] <- null[1] + null[2]
      d$meanLabel[i] <- sprintf("% 1.2f", d$mean[i])
      d$sdLabel[i] <- sprintf("%1.2f", d$sd[i])
      idx <- dist$label == d$label[i]
      y <- dnorm(dist$x[idx], mean = null[1], sd = null[2])
      y <- y/max(y)
      y <- y * 0.7
      dist$yMax[idx] <- d$y[i] - 0.35 + y
      dist$yMin[idx] <- d$y[i] - 0.35
    }
  }
  
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  plot <- ggplot2::ggplot(d, ggplot2::aes(group = .data$label)) +
    ggplot2::geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.2) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_ribbon(ggplot2::aes(x = .data$x, ymax = .data$yMax, ymin = .data$yMin),
                         fill = rgb(1, 0, 0),
                         alpha = 0.6,
                         data = dist)
  
  if (showCis) {
    plot <- plot + ggplot2::geom_ribbon(ggplot2::aes(x = .data$x,
                                                     ymax = .data$yMaxUb,
                                                     ymin = .data$yMaxLb),
                                        fill = rgb(0.8, 0.2, 0.2),
                                        alpha = 0.3,
                                        data = dist)
  }
  plot <- plot + ggplot2::geom_errorbarh(ggplot2::aes(xmax = .data$xMax, xmin = .data$xMin, y = .data$y),
                                         height = 0.5,
                                         color = rgb(0, 0, 0),
                                         size = 0.5) + 
    ggplot2::geom_point(ggplot2::aes(x = mean, y = y), shape = 16, size = 2) + 
    ggplot2::coord_cartesian(xlim = log(limits), ylim = c(0.5, (nrow(d) + 1))) + 
    ggplot2::scale_x_continuous(xLabel, breaks = log(breaks), labels = breaks) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
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
                   plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines"))
  
  text <- data.frame(y = rep(c(d$y, nrow(d) + 1), 3),
                     x = rep(c(1, 2, 3.2), each = nrow(d) + 1),
                     label = c(c(as.character(d$label),
                                 "Source",
                                 d$meanLabel,
                                 " Mean",
                                 d$sdLabel,
                                 "SD")),
                     dummy = "")
  
  data_table <- ggplot2::ggplot(text, ggplot2::aes(x = .data$x,
                                                   y = .data$y,
                                                   label = .data$label)) + 
    ggplot2::geom_text(size = 4, hjust = 0, vjust = 0.5) + 
    ggplot2::geom_hline(ggplot2::aes(yintercept = nrow(d) + 0.5)) + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(colour = "white"),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_line(colour = "white"),
                   strip.background = ggplot2::element_blank(),
                   plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines")) + 
    ggplot2::labs(x = "", y = "") + 
    ggplot2::coord_cartesian(xlim = c(1, 4), ylim = c(0.5, (nrow(d) + 1)))
  
  plot <- gridExtra::grid.arrange(data_table, plot, ncol = 2)
  
  if (!is.null(fileName))
    ggplot2::ggsave(fileName,
                    plot,
                    width = 8 + showCis * 2,
                    height = 1 + length(d) * 0.1,
                    dpi = 400)
  invisible(plot)
}
