\name{d_german}
\alias{d_german}
\docType{data}
\title{German sample Data (complete rankings)}
\description{The German sample dataset (\code{d_german}) is part of a comparative cross-sectional study on political actions and mass participation involving five Western countries. The dataset regards a sample of \eqn{N=2262} German respondents who were asked to rank \eqn{K=4} political goals in order of desirability, namely: 1 = maintaining order in the nation, 2 = giving people more say in the decisions of government, 3 = fighting rising prices and 4 = protecting freedom of speech. The dataset is composed of complete orderings. 
}
\usage{data(d_german)}
\format{
Data matrix of complete orderings with \eqn{N=2262} rows and \eqn{K=4} columns, named as:
  \describe{
    \item{\code{rank1}}{Political goal occupying position 1 (most desirable)}.
    \item{\code{rank2}}{Political goal occupying position 2}.
    \item{\code{rank3}}{Political goal occupying position 3}.
    \item{\code{rank4}}{Political goal occupying position 4 (least desirable)}.
  }
}
\references{
Barnes, S. H. et al. (1979). Political action. Mass participation in five Western democracies. London: Sage.

Croon, M. A. (1989). \emph{Latent class models for the analysis of rankings}. Advances in psychology, 60, pages 99--121.
}
\examples{
data(d_german)
head(d_german)
}
\keyword{datasets}