library(EvidenceSynthesis)
library(Cyclops)
library(survival)
# Cox -------------------------
# settings <- createSimulationSettings(nSites = 1, n = 3000)
# population <- simulatePopulations(settings)[[1]]
# cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
#                                  data = population,
#                                  modelType = "cox"
# )
# param <- "x"

# SCCS -----------------------
settings <- createSccsSimulationSettings(nSites = 1, n = 1000)
population <- simulatePopulations(settings)[[1]]
cyclopsData <- createCyclopsData(y ~ a + x1 + x2 + x3 + x4 + x5 + strata(stratumId),
                                 data = population,
                                 modelType = "cpr"
)
param <- "a"



fit <- fitCyclopsModel(cyclopsData)

# Hermite interpolation -----------------------------------------
hermiteInterpolation <- function(theta, ll, llGradient, x) {
  n <- length(theta)
  result <- numeric(length(x))

  for (i in seq_len(n - 1)) {
    t <- (x - theta[i]) / (theta[i+1] - theta[i])
    h00 <- 2 * t^3 - 3 * t^2 + 1
    h10 <- t^3 - 2 * t^2 + t
    h01 <- -2 * t^3 + 3 * t^2
    h11 <- t^3 - t^2

    idx <- which(x >= theta[i] & x <= theta[i+1])
    result[idx] <- h00[idx] * ll[i] +
      h10[idx] * (theta[i+1] - theta[i]) * llGradient[i] +
      h01[idx] * ll[i+1] +
      h11[idx] * (theta[i+1] - theta[i]) * llGradient[i+1]
  }

  return(result)
}

# Compute profile as grid of points, including MLE if exists:
sampleX <- seq(from = log(0.1), to = log(10), length.out = 5)
profile <- getCyclopsProfileLogLikelihood(object = fit,
                                          parm = param,
                                          x = sampleX,
                                          returnDerivatives = TRUE)
profile$derivative <- -profile$derivative # Bug in current Cyclops version
if (fit$return_flag == "SUCCESS" && coef(fit)[param] > log(0.1) && coef(fit)[param] < log(10)) {
  profile <- profile |>
    bind_rows(data.frame(
      point = coef(fit)[param],
      value = fit$log_likelihood,
      derivative = 0
    ))
}

library(dplyr)

isConvex <- function(df) {
  df <- df %>% arrange(point)  # Ensure points are ordered

  # Compute second derivative approximation
  secondDerivatives <- diff(df$derivative) / diff(df$point)

  # Check if second derivative is non-negative
  convexCheck <- all(secondDerivatives < 0 | abs(secondDerivatives) < .Machine$double.eps)

  # Check if value changes when derivative is nonzero
  valueChangeCheck <- all(df$derivative == 0 | c(diff(df$value) != 0, TRUE))

  return(convexCheck && valueChangeCheck)
}
if (!isConvex(profile)) {
  stop("Issue with convexity")
}

# Plot profile, true LL, and approximation:
x <- seq(from = log(0.1), to = log(10), length.out = 200)
vizData <- data.frame(
  x = x,
  approxY = hermiteInterpolation(profile$point, profile$value, profile$derivative, x),
  trueY = getCyclopsProfileLogLikelihood(object = fit, parm = param, x = x, returnDerivatives = FALSE)$value
)
delta <- 0.1
gradientLineData <- profile |>
  transmute(xMin = point - delta,
            xMax = point + delta,
            yMin = value - derivative * delta,
            yMax = value + derivative * delta)

ggplot(vizData, aes(x = x, y = trueY)) +
  geom_line() +
  geom_point(aes(x = point, y = value), shape = 16, size = 3, color = "blue", data = profile)+
  geom_segment(aes(x = xMin, y = yMin, xend = xMax, yend = yMax), size = 1, color = "blue", alpha = 0.7, data = gradientLineData) +
  geom_line(aes(y = approxY), color = "red", linetype = "dashed", size = 2, alpha = 0.6) +
  scale_x_continuous("Log effect size") +
  scale_y_continuous("Log likelihood")

print(sprintf("MSE: %0.4f", mean(sum((vizData$approxY - vizData$trueY)^2))))





