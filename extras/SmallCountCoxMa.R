# Exploring how bad violations of normality assumption are when synthesizing across small-count stratified Cox regressions
library(Cyclops)
library(ggplot2)
library(survival)
library(meta)

simulateData <- function(nstrata = 10,
                         nrows = 10000,
                         rr = 2,
                         treatedFraction = 0.2,
                         minBackgroundP = 0.0000002,
                         maxBackgroundP = 0.00002) {
  population <- data.frame(rowId = 1:nrows,
                           stratumId = round(runif(nrows,min = 1,max = nstrata)),
                           y = 0,
                           x = as.numeric(runif(nrows) < treatedFraction))
  strataBackgroundProb <- runif(nstrata, min = minBackgroundP, max = maxBackgroundP )
  population$rate <-  strataBackgroundProb[population$stratumId]
  population$rate[population$x == 1] <- population$rate[population$x == 1] * rr
  population$timeToOutcome <- 1 + round(rexp(n = nrows, population$rate))
  population$timeToCensor <- 1 + round(runif(n = nrows, min = 0, max = 499))
  population$time <- population$timeToOutcome
  population$time[population$timeToCensor < population$timeToOutcome] <- population$timeToCensor[population$timeToCensor < population$timeToOutcome]
  population$y <- as.integer(population$timeToCensor > population$timeToOutcome)
  writeLines(paste("People with the outcome:", sum(population$y), ", in target:", sum(population$y[population$x == 1]), ", in comparator:", sum(population$y[population$x == 0])))
  return(population)
}

fitIndividualModel <- function(population, useCyclops = FALSE) {
  if (useCyclops) {
    # Using Cyclops (doesn't really assume normality)
    cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
    fit <- fitCyclopsModel(cyclopsData)
    mode <- coef(fit)
    ci95 <- confint(fit, 1, level = .95)
    return(data.frame(rr = exp(mode),
                      ci95Lb = exp(ci95[2]),
                      ci95Ub = exp(ci95[3]),
                      logRr = mode,
                      seLogRr = (ci95[3] - ci95[2])/(2 * qnorm(0.975))))
  } else {
    # Using vanilla Cox regression:
    fit <- coxph(Surv(time, y) ~ x + strata(stratumId), population)
    ci95 <- confint(fit)
    return(data.frame(rr = exp(fit$coefficients[1]),
                      ci95Lb = ci95[1],
                      ci95Ub = ci95[2],
                      logRr = fit$coefficients[1],
                      seLogRr = sqrt(fit$var[1,1])))
  }
}

plotLikelihood <- function(population) {
  cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
  fit <- fitCyclopsModel(cyclopsData)
  mode <- coef(fit)
  ci95 <- confint(fit, 1, level = .95)
  
  
  # Plot likelihood distribution ----------------------------------------------------
  ci99 <- confint(fit, 1, level = .99)
  if (is.na(ci99[2])) {
    ci99[2] <- mode - (ci99[3] - mode) * 2
  }
  if (is.na(ci99[3])) {
    ci99[3] <- mode + (mode - ci99[2]) * 2
  }
  x <- seq(from = ci99[2], to = ci99[3], length.out = 100)
  y <- rep(0, 100)
  for (i in 1:100) {
    # Set starting coefficient, then tell Cyclops that it is fixed:
    temp <- fitCyclopsModel(cyclopsData, startingCoefficients = x[i], fixedCoefficients = 1)
    y[i] <- temp$log_likelihood
  }
  plot <- ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_line() +
    geom_vline(aes(xintercept = mode), data = data.frame(mode = mode)) +
    geom_vline(linetype = "dotted", aes(xintercept = ci), data = data.frame(ci = ci95[2:3])) +
    xlab("Beta") +
    ylab("Log likelihood")
  return(plot)
}

combinePopulations <- function(populations) {
  highestId <- 0
  for (i in 1:length(populations)) {
    # Making sure stratum IDs are unique:
    populations[[i]]$stratumId <- populations[[i]]$stratumId + highestId
    highestId <- max(populations[[i]]$stratumId)
  }
  pooledPop <- do.call("rbind", populations) 
  return(pooledPop)
}

# Simulation -------------------------------------------------------------------


population <- simulateData()

nDatabases <- 5
rr <- 2

for (i in 1:10) {
  writeLines(paste("Iteration", i))
  set.seed(i)
  populations <- lapply(X = 1:nDatabases, 
                        FUN = simulateData, 
                        nrows = max(1000, round(rnorm(1, 6000, 4000))), 
                        rr = rr)
  
  # Fit individual models, then perform meta-analysis:
  estimates <- lapply(populations, fitIndividualModel)
  estimates <- do.call("rbind", estimates)
  estimates <- estimates[!is.na(estimates$seLogRr) & estimates$seLogRr < 100, ]
  meta <- metagen(TE = estimates$logRr, seTE = estimates$seLogRr)
  # ma <- summary(meta)$random
  ma <- summary(meta)$fixed
  maEstimate <- data.frame(rr = exp(ma$TE),
                           ci95Lb = ma$lower,
                           ci95Ub = ma$upper,
                           logRr = ma$TE,
                           seLogRr = ma$seTE)

    # Pool data, then fit single model:
  pooledPop <- combinePopulations(populations)
  fit <- fitIndividualModel(pooledPop)
  writeLines("Meta-analysis:")
  print(maEstimate, row.names = FALSE)
  writeLines("Pooled:")
  print(fit, row.names = FALSE)
  writeLines("")
}


# Real data -----------------------------------------------------------
folder <- "C:/Home/Research/SmallCountMetaAnalysis/ExampleEnaHtzd"
files <- c("EnaHtzd_CCAE.rds", "EnaHtzd_MDCD.rds", "EnaHtzd_MDCR.rds", "EnaHtzd_Optum.rds", "EnaHtzd_Panther.rds")

# Cleanup
# for (file in files) {
#   x <- readRDS(file.path(folder, file)) 
#   x <- data.frame(rowId = x$rowId,
#                   stratumId = x$stratumId,
#                   x = x$treatment,
#                   y = as.integer(x$outcomeCount != 0),
#                   time = x$survivalTime)
#   saveRDS(x, file.path(folder, file)) 
# }


populations <- lapply(files, function(x) readRDS(file.path(folder, x)))

# Fit individual models, then perform meta-analysis:
estimates <- lapply(populations, fitIndividualModel)
estimates <- do.call("rbind", estimates)
estimates <- estimates[!is.na(estimates$seLogRr) & estimates$seLogRr < 100, ]
meta <- metagen(TE = estimates$logRr, seTE = estimates$seLogRr)
# ma <- summary(meta)$random
ma <- summary(meta)$fixed
maEstimate <- data.frame(rr = exp(ma$TE),
                         ci95Lb = ma$lower,
                         ci95Ub = ma$upper,
                         logRr = ma$TE,
                         seLogRr = ma$seTE)

# Pool data, then fit single model:
pooledPop <- combinePopulations(populations)
fit <- fitIndividualModel(pooledPop)
print(maEstimate)
print(fit)



