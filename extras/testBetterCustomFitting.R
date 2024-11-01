library(EvidenceSynthesis)
library(Cyclops)
library(survival)
# Cox -------------------------
settings <- createSimulationSettings(nSites = 1, n = 3000)
population <- simulatePopulations(settings)[[1]]
cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                 data = population,
                                 modelType = "cox"
)
param <- "x"

# SCCS -----------------------
settings <- createSccsSimulationSettings(nSites = 1, n = 1000)
population <- simulatePopulations(settings)[[1]]
cyclopsData <- createCyclopsData(y ~ a + x1 + x2 + x3 + x4 + x5 + strata(stratumId),
                                 data = population,
                                 modelType = "cpr"
)
param <- "a"



fit <- fitCyclopsModel(cyclopsData)

# Plot likelihood with gradients --------------------------------
x <- seq(from = log(0.1), to = log(10), length.out = 20)
profile <- getCyclopsProfileLogLikelihood(object = fit,
                                          parm = param,
                                          x = x,
                                          returnDerivatives = TRUE)
profile
library(ggplot2)
library(dplyr)
delta <- 0.1
gradientLineData <- profile |>
  transmute(xMin = point - delta,
            xMax = point + delta,
            yMin = value + derivative * delta,
            yMax = value - derivative * delta)

ggplot(profile, aes(x = point, y = value)) +
  geom_point(shape = 16) +
  geom_segment(aes(x = xMin, y = yMin, xend = xMax, yend = yMax), data = gradientLineData, color = "red") +
  scale_x_continuous("Log effect size") +
  scale_y_continuous("Log likelihood")



  if (fit$return_flag == "SUCCESS" && coef(fit)[param] > log(0.1) && coef(fit)[param] < log(10)) {
    # If mode found and in range, set mu to range to ensure gradient = 0

    logRr <- coef(fit)[param]
    x <- seq(log(0.1), log(10), length.out = 100)
    y <- getCyclopsProfileLogLikelihood(fit, param, x)$value
    y <- y - max(y)
    sse <- function(p) {
      error <- customFunction(x, mu, p[1], p[2]) - y
      weights <- exp(y)
      weights[weights < 0.001] <- 0.001
      return(sum(weights * error ^ 2))
    }

    mu <- logRr
    p <- optim(c(1,0), sse)
    sigma <- p$par[1]
    gamma <- p$par[2]
    approximation <- data.frame(mu = mu, sigma = sigma, gamma = gamma)
  } else {
    # Otherwise, use data to fit all 3 parameters
    approximation <- approximateLikelihood(fit, approximation = "custom")
  }
plotLikelihoodFit(approximation = approximation,
                  cyclopsFit = fit,
                  parameter = param)
print(approximation$mu)
approximation <- approximateLikelihood(fit)
plotLikelihoodFit(approximation = approximation,
                  cyclopsFit = fit,
                  parameter = param)
print(approximation$mu)

approximation <- approximateLikelihood(fit, parameter = param, approximation = "adaptive grid")
plotLikelihoodFit(approximation = approximation,
                  cyclopsFit = fit,
                  parameter = param)
print(nrow(approximation))



# Exp taylor series
fitLogLikelihoodFunction <- function(beta, ll, weighByLikelihood = TRUE, pDim = 3, fun = customFunction) {
  sumSquares <- function(p, maxAnchor = TRUE, idx = NULL) {
    approxLl <- fun(beta, p)
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

  beta <- beta[!is.nan(ll)]
  ll <- ll[!is.nan(ll)]

  # Scale to standard (translate in log space so max at 0):
  ll <- ll - max(ll)
  weights <- exp(ll)
  # Weights shouldn't be too small:
  weights[weights < 0.001] <- 0.001
  pInit <- rep(0, pDim)

  fit <- tryCatch(
    {
      suppressWarnings(nlm(sumSquares, pInit, maxAnchor = TRUE))
    },
    error = function(e) {
      list(estimate = pInit, minimum = Inf)
    }
  )
  result <- fit$estimate
  minimum <- fit$minimum

  fit <- tryCatch(
    {
      suppressWarnings(optim(pInit, sumSquares, maxAnchor = TRUE))
    },
    error = function(e) {
      list(par = pInit, value = Inf)
    }
  )
  if (fit$value < minimum) {
    result <- par
    minimum <- fit$value
  }

  # Scale to standard (translate in log space so intersects at (0,0)):
  idx <- which(abs(beta) == min(abs(beta)))
  ll <- ll - ll[idx]
  fit <- tryCatch(
    {
      suppressWarnings(nlm(sumSquares, pInit, maxAnchor = FALSE, idx = idx))
    },
    error = function(e) {
      list(estimate = pInit, minimum = Inf)
    }
  )
  if (fit$minimum < minimum) {
    result <- fit$estimate
    minimum <- fit$minimum
  }
  fit <- tryCatch(
    {
      suppressWarnings(optim(pInit, sumSquares, maxAnchor = FALSE, idx = idx))
    },
    error = function(e) {
      list(par = pInit, value = Inf)
    }
  )
  if (fit$value < minimum) {
    result <- fit$par
    minimum <- fit$value
  }
  threshold <- 0.05
  if (minimum / length(beta) > threshold) {
    warn(paste("Mean squared error greater than ", threshold, ". Probably bad fit."))
  }
  return(result)
}

population <- simulatePopulations(settings)[[1]]
cyclopsData <- createCyclopsData(y ~ a + x1 + x2 + x3 + x4 + x5 + strata(stratumId),
                                 data = population,
                                 modelType = "cpr"
)
param <- "a"

fit <- fitCyclopsModel(cyclopsData)

x <- seq(log(0.01), log(100), length.out = 100)
#idx <- x > log(0.1) & x < log(10)
idx <- seq_along(x)
ll <- EvidenceSynthesis:::getLikelihoodProfile(fit, param, x)
expTaylor <- function(x, p) {
  return(exp(p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3 + p[5]*x^4))
  #return(exp(exp(p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3 + p[5]*x^4)))
}
pet <- fitLogLikelihoodFunction(x[idx], ll[idx], fun = expTaylor, pDim = 5, weighByLikelihood = T)
yet <- expTaylor(x, pet)
p <- approximateLikelihood(fit, parameter = param)
y <- customFunction(x, p$mu, p$sigma, p$gamma)
plot(x,ll - max(ll))
lines(x, yet - max(yet))
lines(x, y-max(y), col = "red")

# newFun <- function(x, mu, sigma, gamma) {
#   return(-(-customFunction(x, mu, sigma, gamma))^0.84)
# }
# newFun <- function(x, mu, sigma, gamma) {
#   return(-(-customFunction(x, mu, sigma, gamma)))
# }
# newFun <- function(x, mu, sigma, gamma) {
#   return(-(abs(x-mu)^gamma)*sigma)
# }
#
# sse <- function(p) {
#   error <- newFun(x, mu, p[1], p[2]) - y
#   return(sum(error ^ 2))
# }
# sse <- function(p) {
#   error1 <- newFun(ci[2], mu, p[1], p[2]) - ll
#   error2 <- newFun(ci[3], mu, p[1], p[2]) - ll
#   return(error1 ^ 2 + error2 ^ 2)
# }
#
# mu <- logRr
# p <- optim(c(1,1), sse)
# sigma <- p$par[1]
# gamma <- p$par[2]
# y2 <- newFun(x, mu, sigma, gamma)
# y2 <- y2 - max(y2)
# library(ggplot2)
# ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
#   geom_line() +
#   geom_line(color = "green", data = data.frame(x =x, y = y2))




