\name{PLMIX}
\alias{PLMIX}
\docType{package}
\title{
Bayesian Analysis of Finite Mixtures of Plackett-Luce Models for Partial Rankings/Orderings
}
\description{
The PLMIX package for R provides functions to fit and analyze finite mixtures of Plackett-Luce models for partial top rankings/orderings within the Bayesian framework. It provides MAP point estimates via EM algorithm and posterior MCMC simulations via Gibbs Sampling. It also fits MLE as a special case of the noninformative Bayesian analysis with vague priors.
}
\details{
The Plackett-Luce model is one of the most popular and frequently applied parametric distributions to analyze partial top rankings/orderings of a finite set of items. The present package allows to account for unobserved sample heterogeneity of partially ranked data with a model-based analysis relying on Bayesian finite mixtures of Plackett-Luce models. The package provides a suite of functions that covers the fundamental phases of a model-based analysis: 
\describe{
    \item{Ranking data manipulation}{\code{\link{rank_ord_switch}} switches data format from rankings to orderings and vice-versa;  \code{\link{make_complete}} performs random completion of partial orderings/rankings data; \code{\link{make_partial}} performs censoring of complete data; \code{\link{binary_group_ind}} constructs binary group membership matrix from mixture component labels; \code{\link{freq_to_unit}} transforms from frequency distribution representation to individual orderings/rankings; \code{\link{unit_to_freq}} does the vice-versa}.
    \item{Ranking data simulation}{\code{\link{rPLMIX}} draws a random sample of complete rankings/orderings according to a finite mixture of plackett-Luce model specification}.
    \item{Ranking data description}{\code{\link{rank_summaries}} provides basic summary statistics of ranking/ordering data, \code{\link{paired_comparisons}} computes the paired comparison frequencies}.
    \item{Model estimation}{\code{\link{mapPLMIX}} performs MAP estimation via EM algorithm; \code{\link{mapPLMIX_multistart}} does the same with multiple stating values; \code{\link{gibbsPLMIX}} provides Bayesian analysis with MCMC posterior simulation via Gibbs sampling; \code{\link{likPLMIX}} and \code{\link{loglikPLMIX}} provide likelihood and loglikelihood evaluation for a mixture of Placket-Luce model}.
    \item{Model selection}{\code{\link{selectPLMIX}} computes alternative Bayesian model selection criteria;  \code{\link{bicPLMIX}} computes the BIC value}.
    \item{Model assessment}{\code{\link{ppcheckPLMIX}} and \code{\link{ppcheckPLMIX_cond}} provide posterior predictive diagnostics}.
}
Data have to be supplied as an object of class \code{matrix} with missing positions denotes as zero entries. For a more efficient implementation of the methods, partial sequences with a single missing entry should be preliminarily filled in, as they correspond to complete rankings/orderings. In the present setting, ties are not allowed. Some quantities frequently recalled in the manual are the following: 
\describe{
    \item{\eqn{N}}{Sample size}.
    \item{\eqn{K}}{Number of possible items}.
    \item{\eqn{G}}{Number of mixture components}.
    \item{\eqn{L}}{Size of the final posterior MCMC sample (after burn-in phase)}.
}
}
\author{
Cristina Mollica and Luca Tardella

Maintainer: Cristina Mollica <\url{cristina.mollica@uniroma1.it}>
}
\references{
Mollica, C., Tardella, L. (2016). Bayesian Plackett-Luce mixture models for partially ranked data. \emph{Psychometrika} (published online): \url{http://dx.doi.org/10.1007/s11336-016-9530-0}
}
\keyword{
Ranking data, Plackett-Luce model, Mixture models, Bayesian inference, MAP estimation, MCMC approximation
}
