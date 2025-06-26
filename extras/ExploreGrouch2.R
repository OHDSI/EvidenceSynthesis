library(dplyr)
library(Cyclops)

# Simulation parameters ----------------------------------------------------------------------------

nDatabases <- 40
nAgeBins <- 3
nPersonsPerDatabase <- round(runif(nDatabases, 100, 10000))
pFemale <- runif(nDatabases)
pAge <- rep(1/100, 100) # Equal probability for every age
baselineIr <- 0.001
betaFemale <- rnorm(1, mean = 0, sd = 1)
betaAge <- splinefun(c(1, 25, 50, 75, 100), c(2, 0, 1, 1.5, 3))
betaDatabase <- rep(0, nDatabases)

# Create two outlier databases:
betaDatabase[3] <- -1
betaDatabase[40] <- 1

# Simulate data ------------------------------------------------------------------------------------
simulateDatabase <- function(i) {
  nPersons <- nPersonsPerDatabase[i]
  population <- tibble(
    female = rbinom(nPersons, 1, pFemale[i]),
    age = sample.int(100, nPersons, replace = TRUE, prob = pAge),
    daysObserved = round(rbeta(nPersons, 2, 3) * 365)
  ) |>
    mutate(
      ir = exp(log(baselineIr) + betaFemale * female + betaAge(age) + betaDatabase[i])
    ) |>
    mutate (
      events = rpois(nPersons, ir * daysObserved)
    )
  summaryStats <- list(
    nPersons = nrow(population),
    nFemale = sum(population$female),
    nAge = population |> group_by(age) |> summarise(count = n()) |> arrange(age),
    nDaysObserved = population |> mutate(daysBin = round(daysObserved / 30)) |> group_by(daysBin) |> summarise(count = n()) |> arrange(daysBin),
    nEvents = population |> summarise(sum(events)) |> pull(),
    nEventPersons = population |> filter(events > 0) |> count() |> pull(),
    databaseId = i
  )
  return(summaryStats)
}

summaryStatsList <- lapply(seq_len(nDatabases), simulateDatabase)

# Prepare data for model fitting -------------------------------------------------------------------

# Define spline
nInternalKnots <- 4
nAge <- bind_rows(lapply(summaryStatsList, function(x) x$nAge)) |>
  group_by(age) |>
  summarise(count = sum(count)) |>
  arrange(age) |>
  mutate(cumCount = cumsum(count))
nTotal <- sum(nAge$count)
knotsAtCumCount <- (1:nInternalKnots / (nInternalKnots+1)) * nTotal
internalKnots <- sapply(knotsAtCumCount, function(x) nAge$age[max(which(nAge$cumCount <= x))])
boundaryKnots <- c(min(nAge$age), max(nAge$age))
designMatrix <- splines::bs(x = boundaryKnots[1]:boundaryKnots[2],
                            knots = internalKnots,
                            Boundary.knots = boundaryKnots,
                            degree = 3)
# Drop one variable to avoid reduncancy:
# designMatrix <- designMatrix[, -1]

summaryStats <- summaryStatsList[[4]]
prepareDatabaseData <- function(summaryStats) {
  # Length of observation is not in the model, and since we assume it is independent from age and sex
  # we just need the total days observed:
  daysObserved <- summaryStats$nDaysObserved |>
    mutate(daysObserved = (15 + daysBin * 30) * count) |>
    summarise(daysObserved = sum(daysObserved)) |>
    pull()

  ageVars <- apply(summaryStats$nAge, 1, function(x) designMatrix[x[1], ] * x[2] / summaryStats$nPersons, simplify = FALSE) |>
    bind_rows() |>
    apply(2, sum) |>
    t() |>
    as_tibble()
  ageVars <- ageVars[, -1]
  colnames(ageVars) <- paste0("age", colnames(ageVars) )
  row <- bind_cols(
    ageVars,
    tibble(
      female = summaryStats$nFemale / summaryStats$nPersons,
      daysObserved = daysObserved,
      databaseId = summaryStats$databaseId,
      nEvents = summaryStats$nEvents
    )
  )
  return(row)
}
data <- lapply(summaryStatsList, prepareDatabaseData)
data <- bind_rows(data) |>
  mutate(across(starts_with("age"), ~ifelse(is.na(.), 0, .)))
data$databaseId <- as.factor(data$databaseId)
ageVarNames <- colnames(data)[grep("^age", colnames(data))]
formula <- as.formula(sprintf("nEvents ~ %s + female + offset(log(daysObserved))",
                              paste(ageVarNames, collapse = " + ")))

# cyclopsData <- createCyclopsData(formula, modelType = "pr", data = data)
# fit <- fitCyclopsModel(cyclopsData, prior = createPrior("laplace", 0.1))
# coefs <- coef(fit)
# # confint(fit, parm = "female")
# betaFemale
#
# y <- apply(designMatrix[, -1] %*% coefs[ageVarNames], 1, sum)
# plot(1:100, y)
# plot(1:100, betaAge(1:100))

# Create prediction intervals using bootstrap ------------------------------------------------------
doPrediction <- function(dummy, data) {
  sampledData <- data[sample.int(nrow(data), nrow(data), replace = TRUE), ]
  cyclopsData <- createCyclopsData(formula, modelType = "pr", data = sampledData)
  fit <- fitCyclopsModel(cyclopsData, prior = createPrior("laplace", 0.1))
  if (fit$return_flag != "SUCCESS") {
    return()
  }
  # Cyclops predict() does not support predicting for new data, except when using sparse representation,
  # so using own implementation instead:
  coefs <- coef(fit)
  prediction <- exp(coefs[1] + apply(t((t(data[ageVarNames]) * coefs[ageVarNames])), 1, sum) + coefs["female"]*data$female) * data$daysObserved
  return(prediction)
}

# t((t(data[ageVarNames]) * coefs[ageVarNames]))[2, 1]
# data[ageVarNames][2, 1] * coefs[ageVarNames][1]


bootstrap <- lapply(1:10000, doPrediction, data = data)
bootstrap <- do.call(rbind, bootstrap)
# predictions <- doPrediction(1, data)

alpha <- 0.01 / nrow(data)
cis <- apply(bootstrap, 2, function(x) quantile(x, c(0.5, alpha/2, 1 - alpha/2)))
rownames(cis) <- c("predictedEvents", "lb", "ub")
dataWithPredictions <- data |>
  bind_cols(as_tibble(t(cis)))

dataWithPredictions <- dataWithPredictions |>
  mutate(outlier = nEvents < lb | nEvents > ub)

which(dataWithPredictions$outlier)

# hist(bootstrap[, 10])
