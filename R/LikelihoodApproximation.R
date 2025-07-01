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

#' Approximate a likelihood function
#'
#' @description
#' Approximate the likelihood function using a parametric (normal, skew-normal, or custom parametric),
#' or grid approximation. The approximation does not reveal person-level information, and can
#' therefore be shared among data sites. When counts are low, a normal approximation might not be
#' appropriate.
#'
#' @param cyclopsFit      A model fitted using the [Cyclops::fitCyclopsModel()] function.
#' @param parameter       The parameter in the `cyclopsFit` object to profile.
#' @param approximation   The type of approximation. Valid options are `'normal'`, `'skew normal'`,
#'                        `'custom'`, `'grid'`, `'adaptive grid'`, or `'grid with gradients'`.
#' @param bounds          The bounds on the effect size used to fit the approximation.
#'
#' @seealso
#' [computeConfidenceInterval], [computeFixedEffectMetaAnalysis], [computeBayesianMetaAnalysis]
#'
#' @return
#' A vector of parameters of the likelihood approximation.
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'   data = populations[[1]],
#'   modelType = "cox"
#' )
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' approximation <- approximateLikelihood(cyclopsFit, "x")
#' approximation
#'
#' # (Estimates in this example will vary due to the random simulation)
#'
#' @export
approximateLikelihood <- function(cyclopsFit,
                                  parameter = 1,
                                  approximation = "custom",
                                  bounds = c(log(0.1), log(10))) {
  if (!approximation %in% c("normal", "skew normal", "custom", "grid", "adaptive grid", "grid with gradients")) {
    stop("'approximation' argument should be 'normal', 'skew normal', 'custom', 'grid', 'adaptive grid', or 'grid with gradients'.")
  }
  if (!is(cyclopsFit, "cyclopsFit")) {
    stop("'cyclopsFit' argument should be of type 'cyclopsFit'")
  }

  if (approximation == "grid with gradients") {
    profile <- getGridWithGradients(cyclopsFit, parameter, bounds)
    return(profile)
  } else if (approximation == "grid") {
    x <- seq(bounds[1], bounds[2], length.out = 1000)
    result <- getLikelihoodProfile(cyclopsFit, parameter, x)
    names(result) <- x
    return(result)
  } else if (approximation == "adaptive grid") {
    result <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, parameter, bounds = bounds, tolerance = 0.001)
    return(result)
  } else if (approximation == "normal") {
    if (cyclopsFit$return_flag != "SUCCESS") {
      return(data.frame(rr = NA, ci95Lb = NA, ci95Ub = NA, logRr = NA, seLogRr = NA))
    } else {
      mode <- coef(cyclopsFit)[parameter]
      ci95 <- confint(cyclopsFit, parameter, level = 0.95)
      return(data.frame(
        rr = exp(mode),
        ci95Lb = exp(as.numeric(ci95[2])),
        ci95Ub = exp(as.numeric(ci95[3])),
        logRr = mode,
        seLogRr = (as.numeric(ci95[3]) - as.numeric(ci95[2])) / (2 * qnorm(0.975))
      ))
    }
  } else {
    x <- seq(bounds[1], bounds[2], length.out = 100)
    ll <- getLikelihoodProfile(cyclopsFit, parameter, x)
    if (approximation == "skew normal") {
      fun <- skewNormal
    } else {
      fun <- customFunction
    }
    result <- fitLogLikelihoodFunction(x, ll, fun = fun)
    if (approximation == "skew normal") {
      names(result)[names(result) == "gamma"] <- "alpha"
    }
    return(result)
  }
}

getGridWithGradients <- function(cyclopsFit, parameter, bounds) {
  # Should probably move this code to Cyclops at some point:
  mleProfile <- data.frame()
  if (cyclopsFit$return_flag == "SUCCESS") {
    mle <- coef(cyclopsFit)[parameter]
    if (mle > bounds[1] && mle < bounds[2] && mle != 0) {
      # There appears to be a MLE, so save that. Note: if the model was fitted
      # with the parameter fixed the coefficient may not actually be the MLE.
      mleProfile <- Cyclops::getCyclopsProfileLogLikelihood(object = cyclopsFit,
                                                         parm = parameter,
                                                         x = mle,
                                                         returnDerivatives = TRUE)
    }
  }
  x <- seq(bounds[1], bounds[2], length.out = 8)
  profile <- Cyclops::getCyclopsProfileLogLikelihood(object = cyclopsFit,
                                                     parm = parameter,
                                                     x = x,
                                                     returnDerivatives = TRUE)
  profile <- bind_rows(mleProfile, profile) |>
    arrange(.data$point)
  invalid <- is.nan(profile$value) |
    is.infinite(profile$value) |
    is.nan(profile$derivative) |
    is.infinite(profile$derivative)
  if (any(invalid)) {
    if (all(invalid)) {
      warning("Failing to compute likelihood at entire grid.")
      return(NULL)
    }

    start <- min(which(!invalid))
    end <- max(which(!invalid))
    if (start == end) {
      warning("Failing to compute likelihood at entire grid except one. Giving up")
      return(NULL)
    }
    profile <- profile[start:end, ]
    invalid <- invalid[start:end]
    if (any(invalid)) {
      warning("Failing to compute likelihood in non-extreme regions. Giving up.")
      return(NULL)
    }
    warning("Failing to compute likelihood at extremes. Truncating bounds.")
  }
  return(profile)
}

#' A custom function to approximate a log likelihood function
#'
#' @details
#' A custom parametric function designed to approximate the shape of the Cox log likelihood function.
#' When `gamma = 0` this function is the normal distribution.
#'
#' @param x       The log(hazard ratio) for which to approximate the log likelihood.
#' @param mu      The position parameter.
#' @param sigma   The scale parameter.
#' @param gamma   The skew parameter.
#'
#' @examples
#' customFunction(x = 0:3, mu = 0, sigma = 1, gamma = 0)
#'
#' @return
#' The approximate log likelihood for the given x.
#'
#' @export
customFunction <- function(x, mu, sigma, gamma) {
  return(((exp(gamma * (x - mu)))) * ((-(x - mu)^2) / (2 * sigma^2)))

  # Derivative: -(exp(gamma * (x - mu)) * (gamma * (x - mu) + 2) * (x - mu))/(2 * sigma^2)
}

#' The skew normal function to approximate a log likelihood function
#'
#' @details
#' The skew normal function. When `alpha = 0` this function is the normal distribution.
#'
#' @param x       The log(hazard ratio) for which to approximate the log likelihood.
#' @param mu      The position parameter.
#' @param sigma   The scale parameter.
#' @param alpha   The skew parameter.
#'
#' @examples
#' skewNormal(x = 0:3, mu = 0, sigma = 1, alpha = 0)
#'
#' @return
#' The approximate log likelihood for the given x.
#'
#' @references
#' Azzalini, A. (2013). The Skew-Normal and Related Families. Institute of Mathematical Statistics
#' Monographs. Cambridge University Press.
#'
#' @export
skewNormal <- function(x, mu, sigma, alpha) {
  if (is.infinite(sigma)) {
    return(rep(0, length(x)))
  }
  return(log(2) + dnorm(x, mu, sigma, log = TRUE) + pnorm(alpha * (x - mu), 0, sigma, log.p = TRUE))
}


#' Cubic Hermite interpolation using both values and gradients to approximate a log likelihood function
#'
#' @param x       The log(hazard ratio) for which to approximate the log likelihood.
#' @param profile A profile as created with `approximateLikelihood()` with
#'                `approximation = "grid with gradients"`. This is a data frame
#'                with 3 columns: `point`, `value`, and `derivative`, sorted by
#'                `point`.
#'
#' @details
#' Performs spline interpolation using cubic Hermite polynomials (Catmull et al.  1974)
#' between the points specified in the profile. We use linear extrapolation
#' outside the points.
#'
#' @examples
#' profile <- data.frame(point = c(1.1, 2.1), value = c(1, 1), derivative = c(0.1, -0.1))
#' hermiteInterpolation(x = 0:3, profile = profile)
#'
#' @return
#' The approximate log likelihood for the given x.
#'
#' @references
#' Catmull, Edwin; Rom, Raphael (1974), "A class of local interpolating splines",
#' in Barnhill, R. E.; Riesenfeld, R. F. (eds.), Computer Aided Geometric Design,
#' New York: Academic Press, pp. 317–326
#'
#' @export
hermiteInterpolation <- function(x, profile) {
  n <- length(profile$point)
  result <- numeric(length(x))

  # Hermite interpolation for points within the range of profile$point:
  for (i in seq_len(n - 1)) {
    t <- (x - profile$point[i]) / (profile$point[i+1] - profile$point[i])
    h00 <- 2 * t^3 - 3 * t^2 + 1
    h10 <- t^3 - 2 * t^2 + t
    h01 <- -2 * t^3 + 3 * t^2
    h11 <- t^3 - t^2

    idx <- which(x >= profile$point[i] & x <= profile$point[i+1])
    result[idx] <- h00[idx] * profile$value[i] +
      h10[idx] * (profile$point[i+1] - profile$point[i]) * profile$derivative[i] +
      h01[idx] * profile$value[i+1] +
      h11[idx] * (profile$point[i+1] - profile$point[i]) * profile$derivative[i+1]
  }

  # Extrapolate to the left (linear):
  leftIdx <- which(x < profile$point[1])
  result[leftIdx] <- profile$value[1] + (x[leftIdx] - profile$point[1]) * profile$derivative[1]

  # Extrapolate to the right (linear):
  rightIdx <- which(x > profile$point[n])
  result[rightIdx] <- profile$value[n] + (x[rightIdx] - profile$point[n]) * profile$derivative[n]

  return(result)
}

fitLogLikelihoodFunction <- function(beta, ll, weighByLikelihood = TRUE, fun = customFunction) {
  sumSquares <- function(p, maxAnchor = TRUE, idx = NULL) {
    approxLl <- fun(beta, p[1], p[2], p[3])
    if (maxAnchor) {
      approxLl <- approxLl - max(approxLl)
    } else {
      approxLl <- approxLl - approxLl[idx]
    }
    if (weighByLikelihood) {
      result <- sum(weights * (approxLl - ll)^2)
    } else {
      result <- sum((approxLl - ll)^2)
    }
    return(result)
  }
  idx <- !is.nan(ll)
  if (!any(idx)) {
    return(data.frame(mu = as.numeric(NA), sigma = as.numeric(NA), gamma = as.numeric(NA)))
  }
  beta <- beta[idx]
  ll <- ll[idx]

  # Scale to standard (translate in log space so max at 0):
  ll <- ll - max(ll)
  weights <- exp(ll)
  # Weights shouldn't be too small:
  weights[weights < 0.001] <- 0.001

  mode <- beta[ll == 0][1]
  if (mode == min(beta) || mode == max(beta)) {
    mode <- 0
  }

  if (min(ll) > -1e-06) {
    return(data.frame(mu = 0, sigma = Inf, gamma = 0))
  }

  fit <- tryCatch(
    {
      suppressWarnings(nlm(sumSquares, c(mode, 1, 0), maxAnchor = TRUE))
    },
    error = function(e) {
      list(estimate = c(0, 0, 0), minimum = Inf)
    }
  )
  result <- data.frame(mu = fit$estimate[1], sigma = fit$estimate[2], gamma = fit$estimate[3])
  minimum <- fit$minimum

  fit <- tryCatch(
    {
      suppressWarnings(optim(c(mode, 1, 0), sumSquares, maxAnchor = TRUE))
    },
    error = function(e) {
      list(par = c(0, 0, 0), value = Inf)
    }
  )
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1], sigma = fit$par[2], gamma = fit$par[3])
    minimum <- fit$value
  }

  # Scale to standard (translate in log space so intersects at (0,0)):
  idx <- which(abs(beta) == min(abs(beta)))
  ll <- ll - ll[idx]
  fit <- tryCatch(
    {
      suppressWarnings(nlm(sumSquares, c(mode, 1, 0), maxAnchor = FALSE, idx = idx))
    },
    error = function(e) {
      list(estimate = c(0, 0, 0), minimum = Inf)
    }
  )
  if (fit$minimum < minimum) {
    result <- data.frame(mu = fit$estimate[1], sigma = fit$estimate[2], gamma = fit$estimate[3])
    minimum <- fit$minimum
  }
  fit <- tryCatch(
    {
      suppressWarnings(optim(c(mode, 1, 0), sumSquares, maxAnchor = FALSE, idx = idx))
    },
    error = function(e) {
      list(par = c(0, 0, 0), value = Inf)
    }
  )
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1], sigma = fit$par[2], gamma = fit$par[3])
    minimum <- fit$value
  }
  threshold <- 0.05
  if (minimum / length(beta) > threshold) {
    warn(paste("Mean squared error greater than ", threshold, ". Probably bad fit."))
  }
  result$sigma <- abs(result$sigma)
  return(result)
}

getLikelihoodProfile <- function(cyclopsFit, parameter, x) {
  ll <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, parameter, x)$value
  return(ll)
}

#' Compute the point estimate and confidence interval given a likelihood function approximation
#'
#' @details
#' Compute the point estimate and confidence interval given a likelihood function approximation.
#'
#' @param approximation   An approximation of the likelihood function as fitted using the
#'                        [approximateLikelihood()] function.
#' @param alpha           The alpha (expected type I error).
#'
#' @return
#' A data frame containing the point estimate, and upper and lower bound of the confidence interval.
#'
#' @examples
#' # Simulate some data for this example:
#' populations <- simulatePopulations()
#'
#' cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#'   data = populations[[1]],
#'   modelType = "cox"
#' )
#' cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
#' approximation <- approximateLikelihood(cyclopsFit, "x")
#' computeConfidenceInterval(approximation)
#'
#' @export
computeConfidenceInterval <- function(approximation, alpha = 0.05) {
  type <- detectApproximationType(approximation)
  if (type == "normal") {
    estimate <- data.frame(
      rr = exp(approximation$logRr),
      lb = exp(approximation$logRr + qnorm(alpha / 2) * approximation$seLogRr),
      ub = exp(approximation$logRr + qnorm(1 - alpha / 2) * approximation$seLogRr),
      logRr = approximation$logRr,
      seLogRr = approximation$seLogRr
    )
    return(estimate)
  } else if (type == "custom") {
    estimate <- computeEstimateFromApproximation(
      approximationFuntion = customFunction,
      a = alpha,
      mu = approximation$mu,
      sigma = approximation$sigma,
      gamma = approximation$gamma
    )
    return(estimate)
  } else if (type == "skew normal") {
    estimate <- computeEstimateFromApproximation(
      approximationFuntion = skewNormal,
      a = alpha,
      mu = approximation$mu,
      sigma = approximation$sigma,
      alpha = approximation$alpha
    )
    return(estimate)
  } else if (type == "adaptive grid") {
    temp <- approximation$value
    names(temp) <- approximation$point
    estimate <- computeEstimateFromGrid(temp, alpha = alpha)
    return(estimate)
  } else if (type == "grid") {
    estimate <- computeEstimateFromGrid(approximation, alpha = alpha)
    return(estimate)
  } else if (type == "grid with gradients") {
    estimate <- computeEstimateFromApproximation(
      approximationFuntion = hermiteInterpolation,
      a = alpha,
      profile = approximation
    )
    return(estimate)
  } else {
    abort(sprintf("Approximation type '%s' not supported by this function", type))
  }
}

#' Detect the type of likelihood approximation based on the data format
#'
#' @param data    The approximation data. Can be a single approximation, or approximations
#'                from multiple sites.
#' @param verbose Should the detected type be communicated to the user?
#'
#' @return
#' A character vector with one of the following values: "normal", "custom", "skew normal",
#' "pooled", "grid", or "adaptive grid".
#'
#' @examples
#' detectApproximationType(data.frame(logRr = 1, seLogRr = 0.1))
#'
#' @export
detectApproximationType <- function(data, verbose = TRUE) {
  if (is.list(data) && !is.data.frame(data)) {
    columnNames <- names(data[[1]])
  } else {
    columnNames <- colnames(data)
  }

  if ("logRr" %in% columnNames) {
    if (verbose) {
      inform("Detected data following normal distribution")
    }
    return("normal")
  } else if ("gamma" %in% columnNames) {
    if (verbose) {
      inform("Detected data following custom parameric distribution")
    }
    return("custom")
  } else if ("alpha" %in% columnNames) {
    if (verbose) {
      inform("Detected data following skew normal distribution")
    }
    return("skew normal")
  } else if ("stratumId" %in% columnNames) {
    if (any(grepl("^x[0-9]+$", columnNames))) {
      if (verbose) {
        inform("Detected (pooled) patient-level SCCS data")
      }
      return("pooled sccs")
    } else {
      if (verbose) {
        inform("Detected (pooled) patient-level data")
      }
      return("pooled")
    }
  } else if ("derivative" %in% columnNames) {
    if (verbose) {
      inform("Detected data following grid with gradients distribution")
    }
    return("grid with gradients")
  } else if ("point" %in% columnNames) {
    if (verbose) {
      inform("Detected data following adaptive grid distribution")
    }
    return("adaptive grid")
  } else {
    if (verbose) {
      inform("Detected data following grid distribution")
    }
    x <- as.numeric(columnNames)
    if (any(is.na(x))) {
      abort("Expecting grid data, but not all column names are numeric")
    }
    return("grid")
  }
}

cleanApproximations <- function(data) {
  type <- detectApproximationType(data, verbose = FALSE)
  if (type == "normal") {
    data <- cleanData(data, c("logRr", "seLogRr"), minValues = c(-100, 1e-05))
  } else if (type == "custom") {
    data <- cleanData(data, c("mu", "sigma", "gamma"), minValues = c(-100, 1e-05, -100))
  } else if (type == "skew normal") {
    data <- cleanData(data,
                      c("mu", "sigma", "alpha"),
                      minValues = c(-100, 1e-05, -10000),
                      maxValues = c(100, 10000, 10000)
    )
  } else if (type == "adaptive grid") {
    for (i in 1:length(data)) {
      cleanedData <- as.data.frame(data[[i]])
      cleanedData$value <- cleanedData$value - max(cleanedData$value)
      cleanedData <- cleanData(cleanedData,
                               c("point", "value"),
                               minValues = c(-100, -1e6),
                               maxValues = c(100, 0)
      )
      data[[i]] <- cleanedData
    }
  } else if (type == "grid with gradients") {
    for (i in 1:length(data)) {
      cleanedData <- as.data.frame(data[[i]])
      cleanedData$value <- cleanedData$value - max(cleanedData$value)
      cleanedData <- cleanData(cleanedData,
                               c("point", "value", "derivative"),
                               minValues = c(-100, -1e6, -1e6),
                               maxValues = c(100, 0, 1e6)
      )
      data[[i]] <- cleanedData
    }
  } else if (type == "grid") {
    if (is.list(data) && !is.data.frame(data)) {
      for (i in 1:length(data)) {
        data <- cleanData(data,
                          colnames(data),
                          minValues = rep(-1e6, ncol(data)),
                          maxValues = rep(0, ncol(data)),
                          grid = TRUE
        )
      }
    } else {
      data <- cleanData(as.data.frame(data),
                        colnames(data),
                        minValues = rep(-1e6, ncol(data)),
                        maxValues = rep(0, ncol(data)),
                        grid = TRUE
      )
    }
  }
  return(data)
}

cleanData <- function(data,
                      columns,
                      minValues = rep(-100, length(columns)),
                      maxValues = rep(100, length(columns)),
                      grid = FALSE) {
  for (i in 1:length(columns)) {
    column <- columns[i]
    if (any(is.infinite(data[[column]]))) {
      if (grid) {
        warn(paste("Estimate(s) with infinite log-likelihood detected. Removing before computing meta-analysis."))
      } else {
        warn(paste(
          "Estimate(s) with infinite",
          column,
          "detected. Removing before computing meta-analysis."
        ))
      }
      data <- data[!is.infinite(data[[column]]), ]
    }
    if (any(is.na(data[[column]]))) {
      if (grid) {
        warn(paste("Estimate(s) with NA log-likelihood detected. Removing before computing meta-analysis."))
      } else {
        warn(paste(
          "Estimate(s) with NA",
          column,
          "detected. Removing before computing meta-analysis."
        ))
      }
      data <- data[!is.na(data[[column]]), ]
    }
    if (any(data[[column]] > maxValues[i])) {
      if (grid) {
        warn(paste("Estimate(s) with positive log-likelihood detected. Removing before computing meta-analysis."))
      } else {
        warn(sprintf(
          "Estimate(s) with extremely high %s (>%s) detected. Removing before computing meta-analysis.",
          column,
          maxValues[i]
        ))
      }
      data <- data[data[[column]] <= maxValues[i], ]
    }
    if (any(data[[column]] < minValues[i])) {
      if (grid) {
        warn(paste("Estimate(s) with extremely low log-likelihood detected. Removing before computing meta-analysis."))
      } else {
        warn(sprintf(
          "Estimate(s) with extremely low %s (<%s) detected. Removing before computing meta-analysis.",
          column,
          minValues[i]
        ))
      }
      data <- data[data[[column]] >= minValues[i], ]
    }
  }
  if (nrow(data) == 0) {
    warn("No estimates left after removing estimates with NA, infinite or extreme values")
  }
  return(data)
}


#' Construct `DataModel` objects from approximate likelihood or profile likelihood data
#'
#' @param data             The likelihood data. Can be a single approximation, approximations
#'                         from multiple sites, or (adaptive) gride profile likelihoods.
#' @param labelReferences  Optional parameter that provides a reference list that
#'                         maps string names to integer indices; only applies to
#'                         "grid" or "adaptive grip" type of data.
#'
#' @examples
#' data("likelihoodProfileLists")
#' dataModel = constructDataModel(likelihoodProfileLists[[1]])
#'
#' @export
constructDataModel <- function(data, labelReferences = NULL){
  type <- detectApproximationType(data)
  data <- cleanApproximations(data)
  if (type == "normal") {
    if (nrow(data) == 0) {
      return(NULL)
    }
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.NormalDataModel")
    for (i in 1:nrow(data)) {
      dataModel$addLikelihoodParameters(
        as.numeric(c(data$logRr[i], data$seLogRr[i])),
        as.numeric(c(NA, NA))
      )
    }
    dataModel$finish()
  } else if (type == "custom") {
    if (nrow(data) == 0) {
      return(NULL)
    }
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.ParametricDataModel")
    for (i in 1:nrow(data)) {
      dataModel$addLikelihoodParameters(
        as.numeric(c(data$mu[i], data$sigma[i], data$gamma[i])),
        as.numeric(c(NA, NA))
      )
    }
    dataModel$finish()
  } else if (type == "skew normal") {
    if (nrow(data) == 0) {
      return(NULL)
    }
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.SkewNormalDataModel")
    for (i in 1:nrow(data)) {
      dataModel$addLikelihoodParameters(
        as.numeric(c(data$mu[i], data$sigma[i], data$alpha[i])),
        as.numeric(c(NA, NA))
      )
    }
    dataModel$finish()
  } else if (type == "pooled") {
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.CoxDataModel")
    for (i in 1:length(data)) {
      dataModel$addLikelihoodData(
        as.integer(data[[i]]$stratumId),
        as.integer(data[[i]]$y),
        as.numeric(data[[i]]$time),
        as.numeric(data[[i]]$x)
      )
    }
    dataModel$finish()
  } else if (type == "pooled sccs") {
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.SccsDataModel")
    for (i in 1:length(data)) {
      x <- data[[i]] %>%
        select(matches("^x[0-9]+$")) %>%
        apply(2, as.numeric)
      dataModel$addLikelihoodData(
        as.integer(data[[i]]$y),
        as.double(data[[i]]$a),
        rJava::.jarray(x, dispatch = TRUE),
        as.integer(data[[i]]$stratumId),
        as.double(data[[i]]$time)
      )
    }
    dataModel$finish()
  } else if (type == "adaptive grid") {
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.ExtendingEmpiricalDataModel")
    for (i in 1:length(data)) {
      if(!is.null(labelReferences) && !is.null(names(data))){
        id = labelReferences[[names(data)[i]]]
        cat(sprintf("Data source %s, with label %s ...\n",
                    names(data)[i], id))
      }else{
        id = i
      }
      dataModel$addLikelihoodParameters(data[[i]]$point, data[[i]]$value, id) # specify identifier
    }
    dataModel$finish()
  } else if (type == "grid") {
    if (nrow(data) == 0) {
      return(NULL)
    }
    x <- as.numeric(colnames(data))
    data <- as.matrix(data)
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.ExtendingEmpiricalDataModel")
    for (i in 1:nrow(data)) {
      if(!is.null(labelReferences) && !is.null(row.names(data))){
        id = labelReferences[[row.names(data)[i]]]
      }else{
        id = i
      }
      dataModel$addLikelihoodParameters(x, data[i, ], id) # specify identifier
    }
    dataModel$finish()
  } else if (type == "grid with gradients") {
    dataModel <- rJava::.jnew("org.ohdsi.metaAnalysis.GridWithGradientsDataModel")
    for (i in 1:length(data)) {
      dataModel$addLikelihoodData(data[[i]]$point, data[[i]]$value, data[[i]]$derivative)
    }
    dataModel$finish()
  }else {
    abort(sprintf("Approximation type '%s' not supported by this function", type))
  }
  return(dataModel)
}


#' Build a list of references that map likelihood names to integer labels for later use
#'
#' @param data             The likelihood data. Can be a single approximation, approximations
#'                         from multiple sites, or (adaptive) grid profile likelihoods.
#'
#' @examples
#' data("likelihoodProfileLists")
#' refLabs = buildLabelReferences(likelihoodProfileLists)
#'
#' @export
buildLabelReferences <- function(data){
  type <- detectApproximationType(data[[1]])
  if(type == "grid" || type == "adaptive grid"){

    labelRefs = list()
    counter = 1
    for(i in 1:length(data)){
      if(type == "grid"){
        thisLabels = row.names(data[[i]])
      }else{
        thisLabels = names(data[[i]])
      }
      # if(i == 1 && is.null(thisLabels)){
      #   return(NULL)
      # }
      for(l in thisLabels){
        if(!l %in% names(labelRefs)){
          labelRefs[[l]] = counter
          counter = counter + 1
        }
      }
    }

    if(length(labelRefs) == 0) {labelRefs = NULL}

  }else{
    labelRefs = NULL
  }

  return(labelRefs)
}

#' Prepare SCCS interval data for pooled analysis
#'
#' @param sccsIntervalData An object of type `SccsIntervalData` as created using the `createSccsIntervalData`
#'                         function in the OHDSI `SelfControlledCaseSeries` package.
#' @param covariateId      The ID of the covariate of interest, for which the estimate will be synthesized.
#'                         All other coveriates will be considered nuisance variables.
#'
#' @return
#' A tibble that can be used in the `computeBayesianMetaAnalysis()` function.
#'
#' @export
prepareSccsIntervalData <- function(sccsIntervalData, covariateId) {
  ensureInstalled("tidyr")
  covariates <- tidyr::pivot_wider(collect(sccsIntervalData$covariates),
                                   names_from = "covariateId",
                                   names_prefix = "x",
                                   values_from = "covariateValue",
                                   values_fill = 0)
  data <- covariates %>%
    rename(a = paste0("x", covariateId)) %>%
    inner_join(collect(sccsIntervalData$outcomes), by = join_by("rowId", "stratumId")) %>%
    select(-"rowId")
  return(data)
}
