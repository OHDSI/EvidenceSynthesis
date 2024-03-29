% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesianBiasCorrection.R
\name{sequentialFitBiasDistribution}
\alias{sequentialFitBiasDistribution}
\title{Fit Bias Distribution Sequentially or in Groups}
\usage{
sequentialFitBiasDistribution(LikelihoodProfileList, ...)
}
\arguments{
\item{LikelihoodProfileList}{A list of lists, each of which is a set of grid profile likelihoods
regarding negative controls, indexed by analysis period ID
for sequential analyses or group ID for group analyses.}

\item{...}{Arguments passed to the \code{\link[=fitBiasDistribution]{fitBiasDistribution()}} function.}
}
\value{
A (long) dataframe with four columns.
Column \code{mean} includes MCMC samples for the average bias,
\code{scale} for the sd/scale parameter,
\code{bias} for predictive samples of the bias, and
\code{Id} for the period ID or group ID.
}
\description{
Learn empirical bias distributions sequentially or in groups; for each sequential
step or analysis group, bias distributions is learned by by simultaneously analyzing
a large set of negative control outcomes by a Bayesian hierarchical model through MCMC.
}
\examples{
# load example data
data("ncLikelihoods")

# fit bias distributions over analysis periods
# NOT RUN
# biasDistributions = sequentialFitBiasDistribution(ncLikelihoods, seed = 42)

}
\seealso{
\link{fitBiasDistribution}, \link{computeBayesianMetaAnalysis}
}
