library(Cyclops)
library(survival)
library(CohortMethod)
library(ggplot2)

createPopulation <- function(n) {
  population <- data.frame(rowId = 1:n,
                           treatment = rep(c(0, 1), n/2))
  return(population)
}

computeTimeToEvent <- function(population, baselineHazard, hrs) {
  population$timeToEvent <- 99999
  start <- 0
  for (i in 1:nrow(hrs)) {
    hazard <- baselineHazard * (1+(hrs$hr[i]-1)*population$treatment)
    timeToEvent <- round(rexp(nrow(population), hazard)) + start
    timeToEvent[timeToEvent > hrs$endDay[i]] <- 99999
    population$timeToEvent[population$timeToEvent > start] <- timeToEvent[population$timeToEvent > start]
    start <- hrs$endDay[i]
  }
  return(population)
}

censorFollowup <- function(population, shape = 3, scale = 180, max = 730) {
  # population$timeToCensor <- round(rnorm(nrow(population), mean, sd))
  population$timeToCensor <- round(rweibull(nrow(population), shape, scale))
  population$timeToCensor[population$timeToCensor < 0] <- 0
  population$timeToCensor[population$timeToCensor > max] <- max
  population$survivalTime <- population$timeToEvent
  population$survivalTime[population$timeToEvent > population$timeToCensor] <- population$timeToCensor[population$timeToEvent > population$timeToCensor]
  population$y <- 1
  population$y[population$timeToEvent > population$timeToCensor] <- 0
  population$outcomeCount <- population$y
  population <- population[population$survivalTime > 0, ]
  return(population)  
}

estimateHr <- function(population, weighted = FALSE) {
  if (weighted) {
    cyclopsData <- createCyclopsData(Surv(survivalTime, y) ~ treatment, weights = population$weight, data = population,  modelType = "cox")
  } else {
    cyclopsData <- createCyclopsData(Surv(survivalTime, y) ~ treatment, data = population,  modelType = "cox")
  }
  fit <- fitCyclopsModel(cyclopsData)
  if (fit$return_flag == "SUCCESS") {
    ci <- exp(confint(fit, "treatment")[2:3])
    hr <- data.frame(hr = exp(coef(fit)), lb = ci[1], ub = ci[2])
  } else {
    hr <- data.frame(hr = NA, lb = NA, ub = NA)
  }
  return(hr)
  # if (weighted) {
  # fit <- coxph(Surv(survivalTime, y) ~ treatment, weights = population$weight, data = population)
  # } else {
  # fit <- coxph(Surv(survivalTime, y) ~ treatment, data = population)
  # } 
  # ci <- exp(confint(fit, "treatment")[1:2])
  # hr <- data.frame(hr = exp(coef(fit)), lb = ci[1], ub = ci[2])
  # return(hr)
}

compareKm <- function(population1, population2, fileName = NULL) {
  sv1 <- survival::survfit(survival::Surv(survivalTime, y) ~ treatment, population1, conf.int = TRUE)
  data1 <- data.frame(time = sv1$time,
                      n.censor = sv1$n.censor,
                      s = sv1$surv,
                      strata = summary(sv1, censored = T)$strata,
                      upper = sv1$upper,
                      lower = sv1$lower,
                      study = "Trial")
  cutoff1 <- quantile(population1$survivalTime, 0.9)
  data1 <- data1[data1$time <= cutoff1, ]
  sv2 <- survival::survfit(survival::Surv(survivalTime, y) ~ treatment, population2, conf.int = TRUE)
  data2 <- data.frame(time = sv2$time,
                      n.censor = sv2$n.censor,
                      s = sv2$surv,
                      strata = summary(sv2, censored = T)$strata,
                      upper = sv2$upper,
                      lower = sv2$lower,
                      study = "Observational")
  cutoff2 <- quantile(population2$survivalTime, 0.9)
  data2 <- data2[data2$time <= cutoff2, ]
  data <- rbind(data1, data2)
  levels(data$strata)[levels(data$strata) == "treatment=0"] <- "Comparator"
  levels(data$strata)[levels(data$strata) == "treatment=1"] <- "Target"
  data$strata <- factor(data$strata, levels = c("Target", "Comparator"))
  
  xLabel <- "Time in days"
  yLabel <- "Survival probability"
  cutoff <- max(cutoff1, cutoff2)
  xlims <- c(0, cutoff)
  
  if (cutoff <= 300) {
    xBreaks <- seq(0, cutoff, by = 50)
  } else if (cutoff <= 600) {
    xBreaks <- seq(0, cutoff, by = 100)
  } else {
    xBreaks <- seq(0, cutoff, by = 250)
  }
  
  ylims <- c(quantile(data$lower, 0.1), 1)
  
  hr1 <- estimateHr(population1)
  hr1$study <- "Trial"
  hr2 <- estimateHr(population2)
  hr2$study <- "Observational"
  hrData <- rbind(hr1, hr2)
  hrData$time <- 1
  hrData$s <- ylims[2] +  (ylims[1] - ylims[2] )* 0.9
  hrData$label <- sprintf("HR = %.2f (%.2f-%.2f)", hrData$hr, hrData$lb, hrData$ub)
  
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = time, y = s)) +
    ggplot2::geom_ribbon(ggplot2::aes(fill = strata,
                                      ymin = lower,
                                      ymax = upper), color = rgb(0, 0, 0, alpha = 0)) +
    ggplot2::geom_step(ggplot2::aes(color = strata), size = 1) +
    ggplot2::geom_label(ggplot2::aes(label = label), hjust = 0, data = hrData) +
    ggplot2::scale_color_manual(values = c(rgb(0.8, 0, 0, alpha = 0.8),
                                           rgb(0, 0, 0.8, alpha = 0.8))) +
    ggplot2::scale_fill_manual(values = c(rgb(0.8, 0, 0, alpha = 0.3),
                                          rgb(0, 0, 0.8, alpha = 0.3))) +
    ggplot2::scale_x_continuous(xLabel, breaks = xBreaks) +
    ggplot2::scale_y_continuous(yLabel) +
    ggplot2::coord_cartesian(ylim=ylims, xlim=xlims) +
    ggplot2::facet_grid(study~.) +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "top")
  if (!is.null(fileName))
    ggsave(plot = plot, filename = fileName, width = 4, height = 4, dpi = 200)
  return(plot)
}

plotAtRiskDistribution <- function(population1, population2, fileName = NULL) {
  computeWeightedAtRisk <- function(time, population) {
    return(sum(population$weight[population$survivalTime >= time]))
  }
  weightedAtRisk2 <- sapply(timesOfInterest, computeWeightedAtRisk, population = population2)
  weightedAtRisk2 <- weightedAtRisk2/max(weightedAtRisk2)
  data <- data.frame(time = rep(timesOfInterest, 3),  atRisk = c(atRisk1, atRisk2, weightedAtRisk2), population = as.factor(rep(c("Trial", "Observational", "Obs. reweighted"), each = length(timesOfInterest))))
  plot <- ggplot(data, aes(x = time, y = atRisk, group = population, color = population)) +
    geom_line(size = 1, alpha = 0.8) +
    scale_x_continuous("Time in days") +
    scale_y_continuous("Fraction at risk") +
    theme(legend.position = "top",
          legend.title = element_blank())
  if (!is.null(fileName))
    ggsave(plot = plot, filename = fileName, width = 4, height = 4, dpi = 200)
  return(plot)
}

computeWeights <- function(population1, population2, nQUantiles = 10) {
  # # Weighting option 1:
  # timesOfInterest <- unique(population2$survivalTime)
  # timesOfInterest <- timesOfInterest[order(timesOfInterest)]
  # computeAtRisk <- function(time, population) {
  #   return(sum(population$survivalTime >= time))
  # }
  # atRisk1 <- sapply(timesOfInterest, computeAtRisk, population = population1)
  # atRisk1 <- atRisk1/max(atRisk1)
  # atRisk2 <- sapply(timesOfInterest, computeAtRisk, population = population2)
  # atRisk2 <- atRisk2/max(atRisk2)
  # weight <- atRisk1 / atRisk2
  # population2$weight <- weight[match(population2$survivalTime, timesOfInterest)]
  # sum(population2$weight > 10)/nrow(population2)
  # population2$weight[population2$weight > 10] <- 10
  # population2 <- population2[population2$weight > 0, ]
  
  # Weighting option 2: 
  q <- c(0, quantile(population2$survivalTime, 1:nQUantiles/nQUantiles))
  x <- cut(population1$survivalTime, q)
  counts1 <- aggregate(count ~ x, data = data.frame(count = 1, x = x), sum)
  counts1$count <- counts1$count/sum(counts1$count)
  x <- cut(population2$survivalTime, q)
  return(counts1$count[match(x, counts1$x)])
}

hrs <- data.frame(hr = c(1.01, 3, 1.5), endDay = c(90, 365, 99999))

baselineHazard <- 0.0001

population1 <- createPopulation(n = 100000)
population1 <- computeTimeToEvent(population = population1, baselineHazard = baselineHazard, hrs = hrs)
population1 <- censorFollowup(population = population1, shape = 5, scale = 365, max = 730)

population2 <- createPopulation(n = 100000)
population2 <- computeTimeToEvent(population = population2, baselineHazard = baselineHazard, hrs = hrs)
population2 <- censorFollowup(population = population2, shape = 1, scale = 180, max = 730)

compareKm(population1, population2, "c:/temp/Km100k.png")

population2$weight <- computeWeights(population1, population2, nQUantiles = 10)

hr <- estimateHr(population2, weighted = TRUE)
sprintf("HR = %.2f (%.2f-%.2f)", hr$hr, hr$lb, hr$ub)

plotAtRiskDistribution(population1, population2, "c:/temp/AtRisk100k.png")

population1 <- createPopulation(n = 4000)
population1 <- computeTimeToEvent(population = population1, baselineHazard = baselineHazard, hrs = hrs)
population1 <- censorFollowup(population = population1, shape = 5, scale = 365, max = 730)

population2 <- createPopulation(n = 40000)
population2 <- computeTimeToEvent(population = population2, baselineHazard = baselineHazard, hrs = hrs)
population2 <- censorFollowup(population = population2, shape = 1, scale = 180, max = 730)

compareKm(population1, population2, "c:/temp/Km4k40k.png")

population2$weight <- computeWeights(population1, population2, nQUantiles = 10)

hr <- estimateHr(population2, weighted = TRUE)
sprintf("HR = %.2f (%.2f-%.2f)", hr$hr, hr$lb, hr$ub)

plotAtRiskDistribution(population1, population2, "c:/temp/AtRisk4k40k.png")
