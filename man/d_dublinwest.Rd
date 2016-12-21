\name{d_dublinwest}
\alias{d_dublinwest}
\docType{data}
\title{Dublin West Data (partial orderings)}
\description{The Dublin West dataset (\code{d_dublinwest}) contains some results of the voting ballots of the 2002 Irish general election from the Dublin West constituency. The Irish voting system allows voters to rank the candidates in order of preferences, rather than only specify the favorite one. In the Dublin West constituency, \eqn{N=29988} voters ranked a maximum of \eqn{K=9} candidates, labeled as: 1 = Bonnie R., 2 = Burton J., 3 = Doherty-Ryan D., 4 = Higgins J., 5 = Lenihan B., 6 = McDonald M., 7 = Morrissey T., 8 = Smyth J. and 9 = Terry S..
The dataset is composed of partial top orderings of varying lengths. Missing positions are denoted with zero entries.}
\usage{data(d_dublinwest)}
\format{
Data matrix of partial orderings with \eqn{N=29988} rows and \eqn{K=9} columns. Each row lists the candidates from the most-liked (\code{rank1}) to the least-liked (\code{rank9}) in a given voting ballot. 
}
\references{
Gormley, I. C., Murphy, T. B. (2008). Exploring Voting Blocs Within the Irish Electorate: A Mixture Modeling Approach. \emph{Journal of the America Statistical Association}, \bold{103}(483), pages 1014--1027.

Gormley, I. C., Murphy, T. B. (2009). A grade of membership model for rank data. \emph{Bayesian
Analysis}, \bold{4}(2),  pages 65--295.

Mattei, N., Walsh, T. (2013) PrefLib: A Library of Preference Data. \emph{Proceedings of Third International Conference on Algorithmic Decision Theory} (ADT 2013). Springer, Lecture Notes in Artificial Intelligence, November 13-15, 2013.

}
\source{
The 2002 Dublin West data have been downloaded from \url{http://www.preflib.org/} Preflib: A Library for Preferences. In that repository, preferences with ties are also included. The original source was publicly available from the Dublin County Returning Officer
at the following URL: \url{https://dublincountyreturningofficer.com/}.}
\examples{
data(d_dublinwest)
head(d_dublinwest)
}
\keyword{datasets}
