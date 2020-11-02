test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

library(survival)

test_that("AML example", {
  gold <- coxph(Surv(time, status) ~ x, data = aml,
                ties = "breslow")
  
  data <- rJava::.jnew("org.ohdsi.data.CoxData", 
                       as.integer(aml$status),
                       as.double(aml$time),
                       as.double(aml$x) - 1)
  
  parameter <- rJava::.jnew("dr.inference.model.Parameter$Default", coef(gold))
  
  likelihood <- rJava::.jnew("org.ohdsi.likelihood.CoxPartialLikelihood",
                             rJava::.jcast(parameter, "dr.inference.model.Parameter"), 
                             data$getSortedData())
  
  tolerance <- 1E-10
  
  expect_equal(likelihood$getLogLikelihood(), 
               as.numeric(logLik(gold)), tolerance = tolerance)
})

test_that("Cox likelihood with failure ties and strata", {

  test <- read.table(header=T, sep = ",", text = "
start, length, event, x1, x2
0, 4,  1,0,0
0, 3,  1,2,0
0, 3,  0,0,1
0, 2,  1,0,1
0, 2,  1,1,1
0, 1,  0,1,0
0, 1,  1,1,0
")
  
  gold <- coxph(Surv(length, event) ~ x1 + strata(x2), test, ties = "breslow")
  
  data <- rJava::.jnew("org.ohdsi.data.CoxData", 
                       as.integer(test$x2),
                       as.integer(test$event),
                       as.double(test$length),
                       as.double(test$x1))
  
  parameter <- rJava::.jnew("dr.inference.model.Parameter$Default", coef(gold))
  
  likelihood <- rJava::.jnew("org.ohdsi.likelihood.CoxPartialLikelihood",
                             rJava::.jcast(parameter, "dr.inference.model.Parameter"), 
                             data$getSortedData())
  
  tolerance <- 1E-10
  
  expect_equal(likelihood$getLogLikelihood(), 
               as.numeric(logLik(gold)), tolerance = tolerance)
})

test_that("normal logPdf", {
  normal <- rJava::.jnew("org.ohdsi.metaAnalysis.NormalDataModel")
  normal$addLikelihoodParameters(
    as.numeric(c(2.0, 3.0)), 
    as.numeric(c(NA, NA)))
  normal$finish()
  likelihood <- normal$getLikelihood()
  parameter <- normal$getCompoundParameter()
  parameter$setParameterValue(as.integer(0), 0.5)
  logPdf <- likelihood$getLogLikelihood()
  
  expect_equal(logPdf, dnorm(0.5, 2, 3, log = TRUE))
})

test_that("skew-normal logPdf and gradLogPdf", {
  skewNormal <- rJava::.jnew("org.ohdsi.metaAnalysis.SkewNormalDataModel")
  skewNormal$addLikelihoodParameters(
    as.numeric(c(2.0, 3.0, 4.0)), 
    as.numeric(c(NA, NA, NA)))
  skewNormal$finish()
  likelihood <- skewNormal$getLikelihood()
  parameter <- skewNormal$getCompoundParameter()
  parameter$setParameterValue(as.integer(0), 0.5)
  logPdf <- likelihood$getLogLikelihood()
  
  expect_equal(logPdf, sn::dsn(0.5, 2, 3, 4, log = TRUE))
  
  analytic <- likelihood$getGradientLogDensity(rJava::.jarray(as.numeric(0.5)))
  
  dx <- function(x, f, eps = 1E-7, ...) {
    (f(x + eps, ...) - f(x - eps, ...)) / (2 * eps)
  }  
  numeric <- dx(x = 0.5, sn::dsn, 
                xi = 2, omega = 3, alpha = 4, log = TRUE)
  
  expect_equal(analytic, numeric, tolerance = 1E-6)
})
