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
settings <- createSccsSimulationSettings(nSites = 1, n = 10000)
population <- simulatePopulations(settings)[[1]]
cyclopsData <- createCyclopsData(y ~ a + x1 + x2 + x3 + x4 + x5 + strata(stratumId),
                                 data = population,
                                 modelType = "cpr"
)
param <- "a"

fit <- fitCyclopsModel(cyclopsData)
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




