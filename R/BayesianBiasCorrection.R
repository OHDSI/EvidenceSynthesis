# Copyright 2023 Observational Health Data Sciences and Informatics
#
# This file is part of EvidenceSynthesis
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Fit Bias Distribution
#'
#' @description
#' Learn an empirical distribution on estimation bias by simultaneously analyzing
#' a large set of negative control outcomes by a Bayesian hierarchical model through MCMC.
#' Analysis is based on a list of extracted likelihood profiles.
#'
#' @param likelihoodProfiles   A list of grid profile likelihoods regarding negative controls.
#' @param priorSds             A two-dimensional vector with the standard deviation of the prior for
#'                             the average bias and the sd/scale parameter, respectively.
#' @param numsamps             Total number of MCMC samples needed.
#' @param thin                 Thinning frequency: how many iterations before another sample is obtained?
#' @param minNCs               Minimum number of negative controls needed to fit a bias distribution;
#'                             default (also recommended): 5.
#' @param robust               Whether or not to use a t-distribution model; default: FALSE.
#' @param seed                 Seed for the random number generator.
#'
#' @seealso
#' [computeBayesianMetaAnalysis]
#'
#' @return
#' A dataframe with three columns and `numsamps` number of rows.
#' Column `mean` includes MCMC samples for the average bias,
#' `scale` for the sd/scale parameter,
#' and `bias` for predictive samples of the bias.
#'
#' @export
fitBiasDistribution <- function(likelihoodProfiles,
                                priorSds = c(2,0.5),
                                numsamps = 10000,
                                thin = 10,
                                minNCs = 5,
                                robust = FALSE,
                                seed = 1){

  # check if there are at least `minNC` likelihoods available
  if(length(likelihoodProfiles) < minNCs){
    message(sprintf('Less than %s NC profile likelihoods available! Not fitting a bias distribution...\n',
                    minNCs))
    return(NULL)
  }

  # fit bias distribution
  null = EvidenceSynthesis::computeBayesianMetaAnalysis(data = likelihoodProfiles,
                                                        chainLength = numsamps + numsamps * thin,
                                                        burnIn = numsamps,
                                                        subSampleFrequency = thin,
                                                        priorSd = priorSds,
                                                        robust = robust,
                                                        seed = seed)
  traces = attr(null, "traces")
  means = traces[,1]
  scales = traces[,2]

  # generate predictive samples of the new bias
  # robust --> t distribution
  # not robust ---> normal
  if(robust){
    biases = ggdist::rstudent_t(numsamps, df = 4, mu = means, sigma = scales)
  }else{
    biases = rnorm(numsamps, means, scales)
  }

  data.frame(mean = means, scale = scales, bias = biases)
}


#' Fit Bias Distribution Sequentially or in Groups
#'
#' @description
#' Learn empirical bias distributions sequentially or in groups; for each sequential
#' step or analysis group, bias distributions is learned by by simultaneously analyzing
#' a large set of negative control outcomes by a Bayesian hierarchical model through MCMC.
#'
#' @param LikelihoodProfileList   A list of lists, each of which is a set of grid profile likelihoods regarding negative controls,
#'                                indexed by analysis period ID for sequential analyses or group ID for group analyses.
#'
#' @seealso
#' [fitBiasDistribution], [computeBayesianMetaAnalysis]
#'
#' @return
#' A (long) dataframe with four columns.
#' Column `mean` includes MCMC samples for the average bias,
#' `scale` for the sd/scale parameter,
#' `bias` for predictive samples of the bias, and
#' `Id` for the period ID or group ID.
#'
#' @export
sequentialFitBiasDistribution <- function(LikelihoodProfileList,
                                          ...){
  if(typeof(LikelihoodProfileList) != 'list'){
    stop('Data format incompatible. Should be a list!')
  }else if(typeof(ncLikelihoods[[1]][[1]]) != 'list'){
    stop('Data format incompatible. Make sure `LikelihoodProfileList` is a list of lists of profile likelihoods.')
  }else if(length(ncLikelihoods) == 1 &&
           typeof(ncLikelihoods[[1]][[1]]) == 'list' &&
           typeof(ncLikelihoods[[1]][[1]][[1]]) == 'double'){
    warning('Looks like a single period/group of profile likelihoods! Deferring to `fitBiasDistribution`...')
    return(fitBiasDistribution(LikelihoodProfileList[[1]], ...))
  }

  # if `LikelihoodProfileList` is not a named list, then give it names...
  if(is.null(names(LikelihoodProfileList))){
    names(LikelihoodProfileList) = as.character(1:length(LikelihoodProfileList))
  }

  groups = names(LikelihoodProfileList)

  res = list()
  for(g in groups){
    res.g = fitBiasDistribution(LikelihoodProfileList[[g]], ...)
    res.g$Id = g
    res[[g]] = res.g
  }

  dplyr::bind_rows(res)

}


#' Bias Correction with Inference
#'
#' @description
#' Perform Bayesian posterior inference regarding an outcome of interest with bias correction
#' using negative control analysis. There is an option to not perform bias correction so that
#' un-corrected results can be obtained.
#'
#' @param likelihoodProfiles      A list of grid profile likelihoods for the outcome of interest.
#' @param ncLikelihoodProfiles    Likelihood profiles for the negative control outcomes.
#'                                Must be a list of lists of profile likelihoods;
#'                                if there is only one analysis period, then this must be a length-1 list,
#'                                with the first item as a list all outcome-wise profile likelihoods.
#' @param biasDistributions       Pre-saved bias distribution(s), formatted as the output
#'                                from [fitBiasDistribution()] or [sequentialFitBiasDistribution()].
#'                                If NULL, then `ncLikelihoodProfiles` must be provided.
#' @param priorMean               Prior mean for the effect size (log rate ratio).
#' @param priorSd                 Prior standard deviation for the effect size (log rate ratio).
#' @param numsamps                Total number of MCMC samples needed.
#' @param thin                    Thinning frequency: how many iterations before another sample is obtained?
#' @param doCorrection            Whether or not to perform bias correction; default: TRUE.
#' @param seed                    Seed for the random number generator.
#'
#' @seealso
#' [approximateSimplePosterior], [fitBiasDistribution]
#'
#' @return
#' A dataframe with five columns, including posterior `median` and `mean` of log RR
#' effect size estimates, 95% credible intervals (`ci95Lb` and `ci95Ub`),
#' posterior probability that log RR > 0 (`p1`), and the period or group ID (`Id`).
#'
#' It is accompanied by the following attributes:
#'
#'   - `samplesCorrected`: all MCMC samples for the bias corrected log RR effect size estimate.
#'   - `samplesRaw`: all MCMC samples for log RR effect size estimate, without bias correction.
#'   - `biasDistributions`: the learned empirical bias distribution from negative control analysis.
#'   - `summaryRaw`: a summary dataframe (same format as in the main result) without bias correction.
#'   - `corrected`: a logical flag indicating if bias correction has been performed; = TRUE if `doCorrection = TRUE`.
#'
#' @export
biasCorrectionInference <- function(likelihoodProfiles,
                                    ncLikelihoodProfiles = NULL,
                                    biasDistributions = NULL,
                                    priorMean = 0,
                                    priorSd = 1,
                                    numsamps = 10000,
                                    thin = 10,
                                    doCorrection = TRUE,
                                    seed = 1, ...){

  if(is.null(ncLikelihoodProfiles) && is.null(biasDistributions)){
    stop("At least one between `ncLikelihoodProfiles` and `biasDistributions` should be provided!")
  }

  # if likelihoodProfiles is not a named list, give it names
  if(is.null(names(likelihoodProfiles))){
    names(likelihoodProfiles) = as.character(1:length(likelihoodProfiles))
  }

  if(!is.null(biasDistributions)){
    # check if bias distributions are usable
    if(length(likelihoodProfiles) == 1 && 'Id' %in% names(biasDistributions)){
      if(names(likelihoodProfiles) %in% unique(biasDistributions$Id)){
        biasDistributions = biasDistributions[biasDistributions$Id == names(likelihoodProfiles),]
      }else{
        stop('Cannot find matching results in the provided `biasDistributions`!')
      }
    }
    if(length(likelihoodProfiles) > 1 && any(!names(likelihoodProfiles) %in% unique(biasDistributions$Id))){
      stop("Likelihood profile IDs and bias distributions IDs do not match!")
    }
  }else{
    # check if all needed periods/groups do have NCs
    if(length(likelihoodProfiles) < length(ncLikelihoodProfiles)){
      if(any(!names(likelihoodProfiles) %in% names(ncLikelihoodProfiles))){
        stop("Outcome of interest likelihoods IDs and negative control likelihoods IDs cannot be matched!")
      }
    }else if(length(likelihoodProfiles) > length(ncLikelihoodProfiles)){
      stop("Insufficient number of negative control profile likelihood groups provided!")
    }


    # fit the bias distribution
    biasDistributions = sequentialFitBiasDistribution(ncLikelihoodProfiles,
                                                      seed = seed,
                                                      numsamps = numsamps,
                                                      thin = thin, ...)
  }

  groups = names(likelihoodProfiles)

  res = list()
  summaryData = list()

  resRaw = list()
  summaryDataRaw = list()

  for(g in groups){
    # inference
    lik = as.numeric(likelihoodProfiles[[g]][[1]]$value)
    names(lik) = as.character(likelihoodProfiles[[g]][[1]]$point)

    mcmc = approximateSimplePosterior(lik,
                                      chainLength = numsamps * thin + 1e5,
                                      burnIn = 1e5,
                                      subSampleFrequency = thin,
                                      priorMean = priorMean,
                                      priorSd = priorSd,
                                      seed = seed)
    samps = mcmc$theta1

    if(!is.null(samps) && length(samps) > 1){
      resRaw[[g]] = data.frame(beta = samps, Id = g)

      betaHDI = HDInterval::hdi(samps, credMass = .95)
      summaryDataRaw[[g]] = data.frame(median = median(samps),
                                       mean = mean(samps),
                                       ci95Lb = betaHDI[1],
                                       ci95Ub = betaHDI[2],
                                       p1 = mean(samps > 0),
                                       Id = g)
    }else{
      warning(sprintf("Cannot obtain MCMC samples for Id=%s! Skipped.", g))
      next
    }

    if(doCorrection){
      if(is.null(biasDistributions) || length(biasDistributions) == 0 || nrow(biasDistributions) == 0){
        warning("Data insufficient for fitting bias distributions. Cannot perform bias correction! Using un-corrected results...")
      }else{
        if("Id" %in% names(biasDistributions)){
          biases = biasDistributions[biasDistributions$Id == g, "bias"]
        }else if(length(likelihoodProfiles) == 1){
          biases = biasDistributions$bias
        }else{
          stop("More than one periods/groups of bias distributions or negative control likelihoods are needed!")
        }

        if(!is.null(biases) && length(biases) > 1){
          # check numsamps match-up
          if(length(samps) <= length(biases)){
            samps = samps - biases[1:length(samps)]
          }else{
            warning('Insufficent bias samples! Trimming effect estimate sample size to fit. Use with caution.')
            samps = samps[1:length(biases)] - biases
          }
        }else{
          warning("Data insufficient for fitting bias distributions. Cannot perform bias correction! Using un-corrected results...")
        }

      }
    }

    if(!is.null(samps) && length(samps) > 1){
      res[[g]] = data.frame(beta = samps, Id = g)

      betaHDI = HDInterval::hdi(samps, credMass = .95)
      summaryData[[g]] = data.frame(median = median(samps),
                                    mean = mean(samps),
                                    ci95Lb = betaHDI[1],
                                    ci95Ub = betaHDI[2],
                                    p1 = mean(samps > 0),
                                    Id = g)
    }else{
      warning(sprintf("Cannot obtain MCMC samples for Id=%s! Skipped.", g))
    }

  }

  samples = dplyr::bind_rows(res)
  summaryData = dplyr::bind_rows(summaryData)
  row.names(summaryData) = NULL

  samplesRaw = dplyr::bind_rows(resRaw)
  summaryDataRaw = dplyr::bind_rows(summaryDataRaw)
  row.names(summaryDataRaw) = NULL

  attr(summaryData, 'samples') = samples
  attr(summaryData, 'samplesRaw') = samplesRaw
  attr(summaryData, 'summaryRaw') = summaryDataRaw
  attr(summaryData, 'biasDistributions') = biasDistributions
  attr(summaryData, 'corrected') = doCorrection

  return(summaryData)
}
