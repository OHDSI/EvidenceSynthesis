# @author Fan Bu
# batch simulation experiments to test hierarchical meta analysis functions
# set up to run in parallel

rJava::.jinit(parameters="-Xmx32g", force.init = TRUE)
options(java.parameters = c("-Xms200g", "-Xmx200g"))

library(EvidenceSynthesis)
library(dplyr)

#### settings of the simulation experiments -----
#setwd("./extras")
cachepath = "cache-6"
if(!dir.exists(cachepath)){dir.create(cachepath)}

# number of repetitions of the experiments
numrep = 50

# list of true RRs
#trueRRs = c(1, 2, 4, 6)
trueRRs = c(1, 1.5, 2, 4)

# simulation setup
nSites = 10
sitePop = 10000
mNegativeControls = 50
meanBias = 0.5
biasStd = 0.1
meanSiteEffect = 0 # has to be fixed at 0 !
siteEffectStd = 0.15

# set random seed
seed0 = 666
set.seed(seed0)

## simulation function
simulateAndApproximate <- function(trueRR, seed = 42, cachePath = NULL, repId = 1){
  metaPopulations = simulateMetaAnalysisWithNegativeControls(meanExposureEffect = log(trueRR),
                                                             nSites = nSites,
                                                             mNegativeControls = mNegativeControls,
                                                             sitePop = sitePop,
                                                             treatedFraction = 0.3,
                                                             minBackgroundHazard = 0.05,
                                                             maxBackgroundHazard = 0.05,
                                                             nStrata = 1,
                                                             meanBias = meanBias,
                                                             biasStd = biasStd,
                                                             meanSiteEffect = meanSiteEffect,
                                                             siteEffectStd = siteEffectStd,
                                                             seed = seed)
  cat("Simulation data generated.\n")

  metaLPs = lapply(metaPopulations, createApproximations, "adaptive grid")
  cat("Profile likelihoods extracted using adaptive grids.\n\n")

  if(!is.null(cachePath)){
    saveRDS(metaLPs, file.path(cachePath, sprintf("metaLPs-%s-%s.rds", trueRR, repId)))
  }

  return(metaLPs)
}

## function to shift likelihood profiles (only do it on the exposure profiles (last entry of meta list))
shiftLP <- function(LPlist, oldRR = 1, newRR = 2, cachePath = NULL, repId = 1){
  shift = log(newRR) - log(oldRR)
  exposureLP = LPlist[[length(LPlist)]]
  for(l in 1:length(exposureLP)){
    exposureLP[[l]]$point = exposureLP[[l]]$point + shift
  }
  LPlist[[length(LPlist)]] = exposureLP

  if(!is.null(cachePath)){
    saveRDS(LPlist, file.path(cachePath, sprintf("metaLPs-%s-%s.rds", newRR, repId)))
  }

  return(LPlist)
}


# summarize estimates in string
estimateAsString <- function(v, digits = 3, exponentiate = TRUE){
  if(exponentiate){v = exp(v)}
  v = round(v, digits = digits)
  s = sprintf("%s (%s, %s)", v[1], v[2], v[3])
  return(s)
}

# summarize estimates in (postMedian, lb95, ub95) tuples, with a method label column
estimatesAsVector <- function(v, label = "standard", exponentiate = TRUE){
  if(exponentiate){v = exp(v)}

  res = as.list(v)
  names(res) = c("Estimate", "LB", "UB")
  res$label = label
  return(res)
}


## fit four different models
## to really update this, need to update `fitModels` function in `hmaSimulationUtils.R`!!!
# fitModels <- function(trueRR, metaLPs, seed = 666, exportAsString = TRUE,
#                       cachePath = "cache2", repId = 1,
#                       ...){
#
#   if(exportAsString){
#     res = list(RR = trueRR)
#   }else{
#     res = NULL
#   }
#
#   print("fit models one by one! ")
#
#
#   # 1 ## fit a joint model with main exposure effect included (uninformative mean):
#   resFile = file.path(cachePath, sprintf("fittedModel-uninformed-%s-%s.rds", trueRR, repId))
#   if(file.exists(resFile)){
#     maWithExposure = readRDS(resFile)
#     cat("loading from fitted model object...\n")
#   }else{
#     settings = generateBayesianHMAsettings(globalExposureEffectPriorMean = c(0),
#                                            globalExposureEffectPriorStd = c(10.0),
#                                            exposureEffectCount = 1,
#                                            ...)
#     maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
#                                                      settings = settings,
#                                                      seed = seed)
#     saveRDS(maWithExposure, file = resFile)
#   }
#
#   # final row is main exposure...
#   estRow = nrow(maWithExposure)
#   cat(estRow, "\n")
#
#   if(exportAsString){
#     res$uninformed = estimateAsString(maWithExposure[estRow,2:4])
#   }else{
#     res = bind_rows(res, as.data.frame(estimatesAsVector(maWithExposure[estRow,2:4], label = "uninformed")))
#   }
#   rm(maWithExposure)
#
#   # 2. try using a "friendlier" prior for the meta-analytic effect
#   resFile = file.path(cachePath, sprintf("fittedModel-informed-%s-%s.rds", trueRR, repId))
#   if(file.exists(resFile)){
#     maWithExposure = readRDS(resFile)
#   }else{
#     settings = generateBayesianHMAsettings(globalExposureEffectPriorMean = c(log(trueRR)),
#                                            globalExposureEffectPriorStd = c(2.0),
#                                            exposureEffectCount = 1,
#                                            ...)
#     maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
#                                                      settings = settings,
#                                                      seed = seed)
#     saveRDS(maWithExposure, file = resFile)
#   }
#
#   estRow = nrow(maWithExposure)
#   if(exportAsString){
#     res$informed = estimateAsString(maWithExposure[estRow,2:4])
#   }else{
#     res = bind_rows(res, as.data.frame(estimatesAsVector(maWithExposure[estRow,2:4], label = "informed")))
#   }
#   rm(maWithExposure)
#
#   #3.(try with the option that does the separate prior on biased effect estimand
#   resFile = file.path(cachePath, sprintf("fittedModel-separate-%s-%s.rds", trueRR, repId))
#   if(file.exists(resFile)){
#     maWithExposure = readRDS(resFile)
#   }else{
#     settings = generateBayesianHMAsettings(globalExposureEffectPriorMean = c(log(trueRR)),
#                                            globalExposureEffectPriorStd = c(10.0),
#                                            exposureEffectCount = 1,
#                                            separateExposurePrior = TRUE, ...)
#     maWithExposure = computeHierarchicalMetaAnalysis(data = metaLPs,
#                                                      settings = settings,
#                                                      seed = seed)
#     saveRDS(maWithExposure, file = resFile)
#   }
#
#   estRow = nrow(maWithExposure)
#
#   cat(estRow, "\n")
#   print(maWithExposure[estRow,2:4])
#   if(exportAsString){
#     res$separate = estimateAsString(maWithExposure[estRow,2:4])
#   }else{
#     res = bind_rows(res, as.data.frame(estimatesAsVector(maWithExposure[estRow,2:4], label = "separate")))
#   }
#   rm(maWithExposure)
#
#   #4. two-stage approach
#   resFile = file.path(cachePath, sprintf("fittedModel-twoStage-%s-%s.rds", trueRR, repId))
#   if(file.exists(resFile)){
#     adjustedMainEffectSamps= readRDS(resFile)
#   }else{
#     settings = generateBayesianHMAsettings(includeExposureEffect = FALSE, ...)
#     maNCsOnly = computeHierarchicalMetaAnalysis(data = metaLPs[1:(length(metaLPs)-1)],
#                                                 seed = seed,
#                                                 settings = settings)
#     traces = attr(maNCsOnly, "traces")
#
#     # newBiasSamps = rnorm(nrow(traces),
#     #                      mean = traces[,"outcome.mean"] + traces[,"source.mean"],
#     #                      sd = 1/sqrt(traces[,"outcome.scale"]))
#     newBiasSamps = traces[,"outcome.mean"] + traces[,"source.mean"] # try the mean bias samples only...
#
#     maExposure = computeBayesianMetaAnalysis(data = metaLPs[[length(metaLPs)]],
#                                              seed = seed, ...)
#     tracesExposure = attr(maExposure, "traces")
#     mainEffectSamps = tracesExposure[,1]
#     adjustedMainEffectSamps = mainEffectSamps - newBiasSamps
#
#     saveRDS(adjustedMainEffectSamps, file = resFile)
#   }
#   if(exportAsString){
#     res$twoStage = estimateAsString(summarizeChain(adjustedMainEffectSamps)[2:4])
#   }else{
#     res = bind_rows(res, as.data.frame(estimatesAsVector(summarizeChain(adjustedMainEffectSamps)[2:4],
#                                                          label = "twoStage")))
#   }
#   rm(adjustedMainEffectSamps)
#
#   res$id = repId
#   if(exportAsString){
#     res = data.frame(res)
#   }else{
#     res$RR = trueRR
#   }
#   return(res)
#
# }



#### start running simulation experiments ----

exportStrings = FALSE # export numeric results for later processing


## parallel function to do simulations ----
paraSimulate <- function(id,
                         trueRRs = c(1, 1.5, 2, 4),
                         cachepath = "cache-6",
                         ...){

  #source("extras/hmaSimulationUtils.R")
  source("./hmaSimulationUtils.R") # working directory should be under "extras/"

  res = NULL

  cat("ROUND", id, "....\n\n")
  this.seed = sample(1000, 1)
  for(trueRR in trueRRs){
    cat("\nSimulation experiments with true RR =", trueRR, "...\n")
    # simulate data
    LPfile = file.path(cachepath, sprintf("metaLPs-%s-%s.rds", trueRR, id))
    if(file.exists(LPfile)){
      metaLPs = readRDS(LPfile)
    }else{
      # make it easier & faster: shifting the exposure effect LPs directly intead of simulating everything
      if(trueRR == 1){
        metaLPs = simulateAndApproximate(trueRR, seed = this.seed, cachePath = cachepath, repId = id)
      }else{
        nullEffectLPsFile = file.path(cachepath, sprintf("metaLPs-1-%s.rds", id))
        if(file.exists(nullEffectLPsFile)){
          nullEffectLPs = readRDS(nullEffectLPsFile)
          metaLPs = shiftLP(nullEffectLPs, oldRR = 1, newRR = trueRR, cachePath = cachepath, repId = id)
          rm(nullEffectLPs)
        }else{
          metaLPs = simulateAndApproximate(trueRR, seed = this.seed, cachePath = cachepath, repId = id)
        }
      }
    }
    # fit models
    this.res = fitModels(trueRR, metaLPs,
                         seed = this.seed,
                         cachePath = cachepath,
                         repId = id,
                         exportAsString = FALSE,
                         chainLength = 300000,
                         burnIn = 5e+04, ...)

    res = dplyr::bind_rows(res, this.res)
    cat("Done!\n\n")

    rm(metaLPs)
  }

  return(res)
}

### run simulations in parallel ----
numrep = 50

## experiments using HMC...
# cluster <- ParallelLogger::makeCluster(14)
# ParallelLogger::clusterRequire(cluster, c("dplyr", "EvidenceSynthesis"))
# simRes <- ParallelLogger::clusterApply(cluster, 1:numrep, paraSimulate)
# ParallelLogger::stopCluster(cluster)
# simRes <- bind_rows(simRes)
#
# saveRDS(simRes, file.path(cachepath, "hmaSimulationsRes-3.rds"))

## experiments using homoscedastic/heteroscedastic model or regular MH/HMC
cluster <- ParallelLogger::makeCluster(14)
ParallelLogger::clusterRequire(cluster, c("dplyr", "EvidenceSynthesis"))
simRes <- ParallelLogger::clusterApply(cluster, 1:numrep, paraSimulate,
                                       cachepath="cache-6",
                                       useHeteroscedasticModel = FALSE, # = TRUE for heteroscedastic model
                                       useHMC = FALSE) # = TRUE is using HMC; HMC is better but slower!
ParallelLogger::stopCluster(cluster)
simRes <- bind_rows(simRes)

saveRDS(simRes, file.path(cachepath, "hmaSimulationsRes-6.rds"))



#### check out results and make plots ----

## load pre-saved results ...
cachepath = "cache-6"
# cachepath = "cache-hetero"
# res = readRDS(file.path(cachepath, "hmaSimulationsRes-3.rds"))
# res = readRDS(file.path(cachepath, "hmaSimulationsRes-4.rds"))
# res = readRDS(file.path(cachepath, "hmaSimulationsRes-5.rds"))
res = readRDS(file.path(cachepath, "hmaSimulationsRes-6.rds"))

## process numeric results in the big dataframe

res = res %>%
  mutate(label = case_when(
    label == "uninformed" ~ "diffuse",
    label == "separate" ~ "separable",
    label == "twoStage" ~ "two-stage",
    TRUE ~ label
  ))

## check out range of estimates (posterior median)
summ = res %>% group_by(label, RR) %>%
  mutate(covered = (LB <= RR) & (RR <= UB)) %>%
  summarise(avgEst = mean(Estimate),
            estIqr1 = quantile(Estimate, 0.1),
            estIqr2 = quantile(Estimate, 0.9),
            coverage = mean(covered)) %>%
  ungroup()


## make plots
library(ggplot2)
library(scales)
library(wesanderson)

## try OHDSI color palette
ohdsiPaletteOrig <- c(
  "#336B91", "#69AED5", "#11A08A", "#FBC511", "#EB6622"
)

ohdsiPalette = scales::alpha(ohdsiPaletteOrig, alpha = 0.6)

## (1) estimates range
ggplot(res, aes(x = label, y = Estimate, fill = as.factor(RR))) +
  geom_hline(yintercept = trueRRs, linetype = 2, size = 0.7, color = "gray60") +
  geom_boxplot(outlier.alpha = 0) +
  scale_y_continuous(trans = log_trans(),
                     breaks = trueRRs,
                     labels = c("1", "1.5", "2", "4"))+
  # scale_x_discrete(limits = c("uninformed", "informed", "separate", "twoStage"),
  #                  labels = c("diffuse", "informed", "separate", "two-stage")) +
  scale_fill_manual(values = ohdsiPalette[2:5], #wes_palette("Royal2")[c(2:5)]
                    breaks = c("4", "2", "1.5", "1")) +
  labs(y = "Posterior median", x = "", fill = "HR") +
  theme_bw(base_size = 16)

## (2) coverage of 95% credible intervals
ggplot(summ, aes(x = label, y = coverage, fill = as.factor(RR))) +
  geom_hline(yintercept = .95, linetype = 2, size = 0.7, color = "gray60") +
  geom_bar(stat = "identity", position = position_dodge()) +
  # scale_x_discrete(limits = c("uninformed", "informed", "separate", "twoStage"),
  #                  labels = c("diffuse", "informed", "separate", "two-stage")) +
  scale_fill_manual(values = wes_palette("Royal2")[c(2:5)]) +
  labs(y = "coverage", x = "", fill = "HR") +
  theme_bw(base_size = 16)

## move the "two-stage" coverage results around a little so can show
summ1 = summ %>%
  mutate(coverage = if_else(label == "two-stage", coverage - 0.01, coverage))

ggplot(summ1, aes(x = RR, y = coverage, color = label)) +
  geom_hline(yintercept = .95, linetype = 2, size = 1, color = "gray60") +
  geom_line(size = 1.5) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = trueRRs) +
  scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values = ohdsiPaletteOrig[2:5]) +
  labs(y = "coverage", x = "HR", color = "method") +
  theme_bw(base_size = 16)

# (3) show 95% CIs by simulations...

## use one RR value as example
pickedRR = 1.5

ggplot(res %>% filter(RR == pickedRR), aes(color = label, x= id, y = Estimate)) +
  geom_hline(yintercept =pickedRR, color = "gray60", size = 1)+
  geom_point(size = 0.8) +
  geom_errorbar(aes(ymax = UB, ymin = LB), width = 0.2) +
  facet_grid(label~.) +
  labs(x = "simulations", y = "HR (95% CI)", color = "variant") +
  theme_bw(base_size = 16)











