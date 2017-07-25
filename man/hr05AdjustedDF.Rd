\name{hr05AdjustedDF}
\alias{hr05AdjustedDF}
\concept{Mahalanobis Distances}
\concept{Outliers}
\title{Adjusted Degrees of Freedom for Testing Robust Mahalanobis Distances for Outlyingness}
\description{
  Computes the degrees of freedom for the adjusted F distribution for testing
  Mahalanobis distances calculated with the minimum covariance 
  determinant (MCD) robust dispersion estimate (for data with 
  a model normal distribution) as described in Hardin and 
  Rocke (2005) or in Green and Martin (2014).
}
\usage{
hr05AdjustedDF( n.obs, p.dim, mcd.alpha, m.asy, method = c("HR05", "GM14"))
}
\arguments{
  \item{n.obs}{(Integer) Number of observations}
  \item{p.dim}{(Integer) Dimension of the data, i.e.,
    number of variables. }
  \item{mcd.alpha}{(Numeric) Value that determines the
    fraction of the sample used to compute the MCD estimate. 
	Default value corresponds to the maximum breakdown point
	case of the MCD.}
  \item{m.asy}{(Numeric) Asymptotic Wishart degrees of freedom. 
  	The default value uses \code{\link{ch99AsymptoticDF}} to 
	obtain the the finite-sample asymptotic value, but the user can also 
	provide a pre-computed value.}
  \item{method}{Either "HR05" to use the method of Hardin and Rocke (2005),
    or "GM14" to use the method of Green and Martin (2014).}
}
\details{
  Hardin and Rocke (2005) derived an approximate \eqn{F} distribution
  for testing robust Mahalanobis distances, computed using the MCD
  estimate of dispersion, for outlyingness. This distribution improves
  upon the standard \eqn{\chi^2} distribution for identifying outlying
  points in data set. The method of Hardin and Rocke was designed to work
  for the maximum breakdown point case of the MCD, where 
  \deqn{\alpha = \lfloor (n.obs + p.dim + 1)/2 \rfloor/n.obs.}{alpha = \lfloor (n.obs + p.dim + 1)/2 \rfloor/n.obs.}
  Green and Martin (2014) extended 
  this result to MCD(\eqn{\alpha}), where \eqn{\alpha} controls the
  size of the sample used to compute the MCD estimate, as well as the
  breakdown point of the estimator.

  With argument \code{method = "HR05"} the function returns 
  \eqn{m_{pred}}{m_pred} as given in Equation 3.4 of Hardin and Rocke (2005).
  The Hardin and Rocke method is only supported for the maximum breakdown
  point case; an error will be generated for other values of \code{mcd.alpha}.

  The argument \code{method = "GM14"} uses the extended methodology 
  described in Green and Martin (2014) and is available for all values
  of \code{mcd.alpha}.
}
\value{ Returns the adjusted F degrees of freedom based on the
  asymptotic value, the dimension of the data, and the sample size.
}
\note{ This function is typically not called directly by users; rather
  it is used in the construction of other functions.
}
\references{
C. G. Green and R. Douglas Martin. An extension of a method of Hardin and Rocke, with
an application to multivariate outlier detection via the IRMCD method of Cerioli.
Working Paper, 2017. Available from 
\url{http://christopherggreen.github.io/papers/hr05_extension.pdf}

J. Hardin and D. M. Rocke. The distribution of robust distances. Journal of Computational
and Graphical Statistics, 14:928-946, 2005. \doi{10.1198/106186005X77685}
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\seealso{
\code{\link{ch99AsymptoticDF}}
}
\examples{
hr05tester <- function(n,p) {
	a <- floor( (n+p+1)/2 )/n
	hr05AdjustedDF( n, p, a, ch99AsymptoticDF(n,p,a)$m.hat.asy, method="HR05" )
}
# compare to m_pred in table on page 941 of Hardin and Rocke (2005)
hr05tester(  50, 5)
hr05tester( 100,10)
hr05tester( 500,10)
hr05tester(1000,20)

# using default arguments
hr05tester <- function(n,p) {
	hr05AdjustedDF( n, p, method="HR05" )
}
# compare to m_pred in table on page 941 of Hardin and Rocke (2005)
hr05tester(  50, 5)
hr05tester( 100,10)
hr05tester( 500,10)
hr05tester(1000,20)

# Green and Martin (2014) improved method
hr05tester <- function(n,p) {
	hr05AdjustedDF( n, p, method="GM14" )
}
# compare to m_sim in table on page 941 of Hardin and Rocke (2005)
hr05tester(  50, 5)
hr05tester( 100,10)
hr05tester( 500,10)
hr05tester(1000,20)
}
\keyword{ robust }
\keyword{ multivariate }
