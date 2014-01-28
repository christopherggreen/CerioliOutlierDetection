\name{hr05CriticalValue}
\alias{hr05CriticalValue}
\title{ Hardin and Rocke (2005) Critical Value for Testing MCD-based Mahalanobis Distances
}
\description{
  Hardin and Rocke (2005) provide an approximate \eqn{F} distribution for
  testing whether Mahalanobis distances calculated using the MCD dispersion
  estimate are unusually large, and hence, indicative of outliers in the data.
}
\usage{
hr05CriticalValue(em, p.dim, signif.alpha)
}
\arguments{
  \item{em}{(Numeric) Degrees of freedom for Wishart distribution 
    approximation to the MCD scatter matrix.} 
  \item{p.dim}{(Integer) Dimension of the data, i.e.,
    number of variables. }
  \item{signif.alpha}{ (Numeric) Significance level for testing the 
    null hypothesis }
}
\details{
  Hardin and Rocke (2005) derived an \eqn{F} distributional approximation
  for the Mahalanobis distances of the observations that were excluded from
  the MCD calculation; see equation 3.2 on page 938 of the paper.

  It is assumed here that the MCD covariance estimate used in the Mahalanobis
  distance calculation was adjusted by the consistency factor, so it is not 
  included in the calculation here. (If one needs the consistency factor it
  is returned by the function \code{\link{ch99AsymptoticDF}} in this package
  or by the private function \code{MCDcons} in the \code{robustbase} 
  package.
}
\value{ The appropriate cutoff value (from the \eqn{F} 
  distributional approximation) for testing whether a
  Mahalanobis distance is unusually large at the specified
  significance level.
}
\note{ It can happen that one of the \eqn{F} distribution paramaters, 
  \eqn{m - p + 1}, is non-positive, in which case \code{qf} will return
  NaN. \code{hr05CriticalValue} will issue a warning in this case, and
  return NA.
}
\references{
J. Hardin and D. M. Rocke. The distribution of robust distances. Journal of Computational
and Graphical Statistics, 14:928-946, 2005.
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\seealso{
\code{\link{hr05AdjustedDF}}, \code{\link{hr05CutoffMvnormal}}
}
\examples{
hr05CriticalValue( hr05AdjustedDF( 1000, 20 ), 20, 0.05 ) 
}
\keyword{ robust }
\keyword{ multivariate }
