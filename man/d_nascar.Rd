\name{d_nascar}
\alias{d_nascar}
\docType{data}
\title{NASCAR Data (partial orderings)}
\description{The NASCAR dataset (\code{d_nascar}) collects the results of the 2002 season of stock car racing held in the United States. The 2002 championship consisted of \eqn{N=36} races, with 43 car drivers competing in each race. A total of \eqn{K=87} drivers participated in the 2002 season, taking part to a different number of races: some of them competed in all the races, some others in only one. The results of the entire 2002 season were collected in the form of top-43 orderings, where the position of the not-competing drivers in each race is assumed lower than 43th, but undetermined. Missing positions are denoted with zero entries.}
\usage{data(d_nascar)}
\format{
 Data matrix of partial orderings with \eqn{N=36} rows and \eqn{K=87} columns. Each row lists the car drivers from the top position (\code{rank1}) to the bottom one (\code{rank87}) in a given race. Columns from 44th to 87th are filled with zeros, because only 43 drivers competed in each race.
}
\references{
Caron, F., Doucet, A. (2012). Efficient Bayesian inference for Generalized Bradley-Terry models. \emph{J. Comput. Graph. Statist.}, \bold{21}(1), pages 174--196.

Guiver, J., Snelson, E. (2009). Bayesian inference for Plackett-Luce ranking models. In Bottou, L. and Littman, M., editors, \emph{Proceedings of the 26th International Conference on Machine Learning - ICML 2009}, pages 377--384. Omnipress.

Hunter, D. R. (2004). MM algorithms for Generalized Bradley-Terry models. \emph{Ann. Statist.}, \bold{32}(1), pages 384--406.
}
\source{
The NASCAR dataset in the MATLAB format used by Hunter, D. R. (2004) can be downloaded from \url{http://sites.stat.psu.edu/~dhunter/code/btmatlab/}. At the same link, a .xls file with drivers' names is also available.}
\examples{
data(d_nascar)
head(d_nascar)
}
\keyword{datasets}
