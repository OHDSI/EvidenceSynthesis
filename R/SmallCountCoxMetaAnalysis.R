profileCoxLikelihood <- function(population) {
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), 
                                            data = population, 
                                            modelType = "cox")
  fit <- Cyclops::fitCyclopsModel(cyclopsData)
  mode <- coef(fit)
  ci95 <- confint(fit, 1, level = .95)
  ci99 <- confint(fit, 1, level = .99)
  if (is.na(ci99[2])) {
    ci99[2] <- mode - (ci99[3] - mode) * 2
  }
  if (is.na(ci99[3])) {
    ci99[3] <- mode + (mode - ci99[2]) * 2
  }
  beta <- seq(from = ci99[2], to = ci99[3], length.out = 100)
  ll <- rep(0, 100)
  for (i in 1:100) {
    # Set starting coefficient, then tell Cyclops that it is fixed:
    temp <- Cyclops::fitCyclopsModel(cyclopsData, startingCoefficients = beta[i], fixedCoefficients = 1)
    ll[i] <- temp$log_likelihood
  }
  result <- data.frame(beta = beta, ll = ll)
  attr(result, "mode") <- mode
  attr(result, "se") <- (ci95[3] - ci95[2]) / (2*qnorm(0.975))
  return(result)
}

pseudoCoxLogLikelihood <- function(x, position, scale, skew, shape) {
  return((1/scale)*(x - position) - log(exp(skew*(1/scale)*(x - position)) + shape))
}

fitPseudoCox <- function(coxProfile) {
  
  sumSquares <- function(p, coxProfile) {
    pseudoCoxLl <- pseudoCoxLogLikelihood(coxProfile$beta, p[1], p[2], p[3], p[4])
    coxLl <- coxProfile$ll
    
    # Scale to standard:
    coxLl <- coxLl - min(coxLl)
    coxLl <- coxLl / max(coxLl)
    pseudoCoxLl <- pseudoCoxLl - min(pseudoCoxLl)
    pseudoCoxLl <- pseudoCoxLl / max(pseudoCoxLl)
    
    result <- sum((pseudoCoxLl - coxLl)^2)
    print(paste(paste(p, collapse = ","), result))
    return(result)
  }
  
  fit <- optim(c(3,12,2.2,1), sumSquares, coxProfile = coxProfile)
  
  result <- data.frame(position = fit$par[1],
                       scale  = fit$par[2],
                       skew = fit$par[3],
                       shape = fit$par[4])
  return(result)
}

plotLikelihoodFit <- function(coxProfile, pseudoCoxFit) {
  pseudoCoxLl <- pseudoCoxLogLikelihood(x = coxProfile$beta, 
                                        position = pseudoCoxFit$position,
                                        scale = pseudoCoxFit$scale,
                                        skew = pseudoCoxFit$skew,
                                        shape = pseudoCoxFit$shape)
  mean <- attr(coxProfile, "mode")
  se <- attr(coxProfile, "se")
  normalLl <- dnorm(coxProfile$beta, mean, se, log = TRUE) 
  coxLl <- coxProfile$ll
  
  # Scale to standard:
  coxLl <- coxLl - min(coxLl)
  coxLl <- coxLl / max(coxLl)
  normalLl <- normalLl - min(normalLl)
  normalLl <- normalLl / max(normalLl)
  pseudoCoxLl <- pseudoCoxLl - min(pseudoCoxLl)
  pseudoCoxLl <- pseudoCoxLl / max(pseudoCoxLl)
  
  
  plotData <- rbind(data.frame(beta = coxProfile$beta,
                               ll = coxLl,
                               label = "Cox"),
                    data.frame(beta = coxProfile$beta,
                               ll = pseudoCoxLl,
                               label = "Pseudo Cox"),
                    data.frame(beta = coxProfile$beta,
                               ll = normalLl,
                               label = "Normal"))
  ggplot2::ggplot(plotData, ggplot2::aes(x = beta, y = -ll, group = label, color = label)) +
    ggplot2::geom_line(size = 1)
}