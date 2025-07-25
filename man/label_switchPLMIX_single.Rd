% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLMIXfunctions.R
\name{label_switchPLMIX_single}
\alias{label_switchPLMIX_single}
\title{Label switching adjustment for mixtures of Plackett-Luce models}
\usage{
label_switchPLMIX_single(pi_inv, G, MCMCsampleP, MCMCsampleW, MAPestP, MAPestW)
}
\arguments{
\item{pi_inv}{An object of class \code{top_ordering}, collecting the numeric \eqn{N}\eqn{\times}{x}\eqn{K} data matrix of partial orderings, or an object that can be coerced with \code{\link{as.top_ordering}}.}

\item{G}{Number of mixture components.}

\item{MCMCsampleP}{Numeric \eqn{L}\eqn{\times}{x}\eqn{G*K} matrix with the MCMC samples of the component-specific support parameters to be processed.}

\item{MCMCsampleW}{Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix with the MCMC samples of the mixture weights to be processed.}

\item{MAPestP}{Numeric \eqn{G}\eqn{\times}{x}\eqn{K} matrix of MAP component-specific support parameter estimates to be used as pivot in the PRA method.}

\item{MAPestW}{Numeric vector of the \eqn{G} MAP estimates of the mixture weights as pivot in the PRA method.}
}
\value{
A list of named objects:

 \item{\code{final_sampleP}}{ Numeric \eqn{G}\eqn{\times}{x}\eqn{K}\eqn{\times}{x}\eqn{L} array MCMC samples of the component-specific support parameters adjusted for label switching.}
 \item{\code{final_sampleW}}{ Numeric \eqn{L}\eqn{\times}{x}\eqn{G} matrix of MCMC samples of the mixture weights adjusted for label switching.}
}
\description{
Remove the label switching phenomenon from the MCMC samples of Bayesian mixtures of Plackett-Luce models with a different number of components.
}
\details{
The \code{label_switchPLMIX} function performs the label switching adjustment of the MCMC samples via the Pivotal Reordering Algorithm (PRA) described in Marin et al (2005), by recalling the \code{\link[label.switching]{pra}} function from the \code{\link[label.switching]{label.switching}} package.
}
\author{
Cristina Mollica and Luca Tardella
}
