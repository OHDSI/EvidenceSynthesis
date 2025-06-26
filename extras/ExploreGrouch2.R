library(dplyr)
library(Cyclops)

nDatabases <- 40
nAgeBins <- 3
nPersonsPerDatabase <- round(runif(nDatabases, 100, 10000))
pFemale <- runif(nDatabases)
pAgeBin <- t(apply(matrix(runif(nDatabases * nAgeBins), ncol = nAgeBins), 1, function(x) x / sum(x)))
baselineIr <- 0.001
betaFemale <- rnorm(1, mean = 0, sd = 1)
betaAgeBin <- rnorm(nAgeBins, mean = 0, sd = 1)
betaDatabase <- c(rep(0, nDatabases - 1), 1)

# Simulation
simulateDatabase <- function(i) {
  nPersons <- nPersonsPerDatabase[i]
  population <- tibble(
    female = rbinom(nPersons, 1, pFemale[i]),
    ageBin = sample(seq_len(nAgeBins), nPersons, replace = TRUE, prob = pAgeBin[i, ]),
    daysObserved = round(rbeta(nPersons, 2, 3) * 365)
  ) |>
    mutate(
      ir = exp(log(baselineIr) + betaFemale * female + betaAgeBin[ageBin] + betaDatabase[i])
    ) |>
    mutate (
      events = rpois(nPersons, ir * daysObserved)
    )
  #mean(population$female)
  summaryStats <- list(
    nPersons = nrow(population),
    nFemale = sum(population$female),
    nAgeBin = population |> group_by(ageBin) |> summarise(count = n()) |> arrange(ageBin),
    nDaysObserved = population |> mutate(daysBin = round(daysObserved / 30)) |> group_by(daysBin) |> summarise(count = n()) |> arrange(daysBin),
    nEvents = population |> summarise(sum(events)) |> pull(),
    nEventPersons = population |> filter(events > 0) |> count() |> pull(),
    databaseId = i
  )
  return(summaryStats)
}

summaryStatsList <- lapply(seq_len(nDatabases), simulateDatabase)

summaryStats <- summaryStatsList[[4]]
prepareDatabaseData <- function(summaryStats) {
  # Length of observation is not in the model, and since we assume it is independent from age and sex
  # we just need the total days observed:
  daysObserved <- summaryStats$nDaysObserved |>
    mutate(daysObserved = (15 + daysBin * 30) * count) |>
    summarise(daysObserved = sum(daysObserved)) |>
    pull()

  # Remove first age bin to avoid redundancy:
  nAgeBins <- summaryStats$nAgeBin |>
    filter(ageBin != 1)
  ageVars <- t(nAgeBins$count) / summaryStats$nPersons
  colnames(ageVars) <- paste0("age", nAgeBins$ageBin)
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

# cyclopsData <- createCyclopsData(nEvents ~ age2 + age3 + female + databaseId + offset(log(daysObserved)), modelType = "pr", data = data)
# fit <- fitCyclopsModel(cyclopsData)
# coef(fit)
# estimates <- coef(fit)[paste0("databaseId", 2:length(summaryStatsList))]
# cis <- confint(fit, parm = paste0("databaseId", 2:length(summaryStatsList)))
# ses <- (cis[, 3] - cis[, 2]) /  2*qnorm(0.975)
# zs <- estimates / ses
# ps <- 2 * pmin(pnorm(zs), 1 - pnorm(zs))
# ps
#
cyclopsData <- createCyclopsData(nEvents ~ age2 + age3 + female + offset(log(daysObserved)), modelType = "pr", data = data)
fit <- fitCyclopsModel(cyclopsData)
coef(fit)
confint(fit, parm = 2:4)
betaAgeBin[2:length(betaAgeBin)] - betaAgeBin[1]
betaFemale
#
# plot(data$female, data$nEvents / data$daysObserved)
# betaFemale
# plot(data$age2, data$nEvents / data$daysObserved)
# plot(data$age3, data$nEvents / data$daysObserved)
# betaAgeBin[2:length(betaAgeBin)] - betaAgeBin[1]


# Create prediction intervals using bootstrap ------------------------------------------------------
doPrediction <- function(dummy, data, doSample = FALSE) {
  sampledData <- data[sample.int(nrow(data), nrow(data), replace = TRUE), ]
  cyclopsData <- createCyclopsData(nEvents ~ age2 + age3 + female + offset(log(daysObserved)), modelType = "pr", data = sampledData)
  fit <- fitCyclopsModel(cyclopsData)
  # Cyclops predict() does not support predicting for new data, except when using sparse representation,
  # so using own implementation instead:
  coefs <- coef(fit)
  prediction <- exp(coefs[1] + coefs["age2"]*data$age2 + coefs["age3"]*data$age3 + coefs["female"]*data$female) * data$daysObserved
  return(prediction)
}

bootstrap <- lapply(1:10000, doPrediction, data = data, doSample = TRUE)
bootstrap <- do.call(rbind, bootstrap)
# predictions <- doPrediction(1, data)

alpha <- 0.05 / nrow(data)
cis <- apply(bootstrap, 2, function(x) quantile(x, c(0.5, alpha/2, 1 - alpha/2)))
rownames(cis) <- c("predictedEvents", "lb", "ub")
dataWithPredictions <- data |>
  bind_cols(as_tibble(t(cis)))

dataWithPredictions <- dataWithPredictions |>
  mutate(outlier = nEvents < lb | nEvents > ub)

which(dataWithPredictions$outlier)

hist(bootstrap[, 10])
