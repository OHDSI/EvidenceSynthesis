library(Cyclops)
library(survival)
library(meta)
n <- c(1000, 100, 2000, 500)
dbs <- 4
hr <- 1.5
baselineHazard <- 0.0001

perDbCoverage <- rep(0, dbs*100)
maCoverage <- rep(0, 100)
smaller <- rep(0,100)
for (i in 1:100) {
  logRr <- rep(0, dbs)
  seLogRr <- rep(0, dbs)
  for (d in 1:dbs) {
    # set.seed(i)
    population <- data.frame(rowId = 1:n[d],
                             treatment = rep(c(FALSE, TRUE), n[d]/2),
                             stratumId = rep(1:(n[d]/2), each=2))
    hazard <- baselineHazard * (1+(hr-1)*population$treatment)
    timeToEvent <- round(rexp(n[d], hazard))
    timeToCensor <- round(rexp(n[d], 0.01))
    population$survivalTime <- timeToEvent
    population$survivalTime[timeToEvent > timeToCensor] <- timeToCensor[timeToEvent > timeToCensor]
    population$y <- 1
    population$y[timeToEvent > timeToCensor] <- 0
    # print(sum(population$y[population$treatment]))
    # print(sum(population$y[!population$treatment]))
    print(sum(population$y))
    cyclopsData <- createCyclopsData(Surv(survivalTime, y) ~ treatment + strata(stratumId), data = population,  modelType = "cox")
    fit <- fitCyclopsModel(cyclopsData)
    if (fit$return_flag == "SUCCESS") {
      ci <- exp(confint(fit, "treatmentTRUE")[2:3])
      logRr[d] <- coef(fit)
      seLogRr[d] <- (log(ci[2]) - log(ci[1])) / (2*qnorm(0.975))
      if (is.na(seLogRr[d])) 
        perDbCoverage[i*d] <- NA
      else
        perDbCoverage[i*d] <- hr >= ci[1] & hr <= ci[2]
    } else {
      logRr[d] <- NA
      seLogRr[d] <- NA
      perDbCoverage[i*d] <- NA
    }
  }
  logRr <- logRr[!is.na(seLogRr)]
  seLogRr <- seLogRr[!is.na(seLogRr)]
  if (length(logRr) != 0) {
    meta <- metagen(TE = logRr, seTE = seLogRr, hakn = TRUE)
    rnd <- summary(meta)$random
    maCoverage[i] <- hr >= exp(rnd$lower) & hr <= exp(rnd$upper)
    fix <- summary(meta)$fixed
    smaller[i] <- rnd$seTE < fix$seTE
  } else {
    maCoverage[i] <- NA
    smaller[i] <- NA
  }
}
mean(perDbCoverage, na.rm = TRUE)
mean(maCoverage, na.rm = TRUE)
mean(smaller, na.rm = TRUE)
