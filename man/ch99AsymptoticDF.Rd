\name{ch99AsymptoticDF}
\alias{ch99AsymptoticDF}
\title{
  Croux and Haesbroeck (1999) finite-sample asymptotic approximation
  parameters for the MCD estimate
}
\description{
  Computes the asymptotic Wishart degrees of freedom and
  consistency constant for the MCD robust dispersion estimate 
  (for data with a model normal distribution) as described in 
  Hardin and Rocke (2005) and using the formulas described in 
  Croux and Haesbroeck (1999).
}
\usage{
ch99AsymptoticDF(n.obs, p.dim, mcd.alpha)
}
\arguments{
  \item{n.obs}{(Integer) Number of observations}
  \item{p.dim}{(Integer) Dimension of the data, i.e.,
    number of variables. }
  \item{mcd.alpha}{ (Numeric) Value that  
    determines the fraction of the sample used to 
    compute the MCD estimate. \eqn{1 - mcd.alpha}{1 - mcd.alpha}
    will be the fraction of observations that are
    omitted in computing the MCD estimate. The
    default value is 
    \deqn{ \lfloor (n.obs + p.dim + 1)/2 \rfloor/n.obs,}{\lfloor (n.obs + p.dim + 1)/2 \rfloor/n.obs,}  
    which yields the MCD estimate with the maximum possible
    breakdown point.
}
}
\value{
  \item{c.alpha}{the asymptotic consistency coefficient for the MCD
    estimate of the dispersion matrix}
  \item{m.hat.asy}{the asymptotic degrees of freedom for the Wishart
    distribution approximation to the distribution of the MCD dispersion
    estimate}
}
\details{
  The consistency factor \code{c.alpha} is already available in the
  \code{robustbase} library as the function 
  \code{.MCDcons}. (See the code for \code{covMcd}.) \code{ch99AsymptoticDF}
  uses the result of \code{.MCDcons} for consistency.

  The computation of the asymptotic Wishart degrees of freedom parameter \code{m}
  follows the Appendix of Hardin and Rocke (2005).
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\references{
Christopher Croux and Gentiane Haesbroeck. Influence function and efficiency of the minimum
covariance determinant scatter matrix estimator. Journal of Multivariate Analysis,
71:161-190, 1999.

J. Hardin and D. M. Rocke. The distribution of robust distances. Journal of Computational
and Graphical Statistics, 14:928-946, 2005.
}
\examples{
# compare to table from p941 of Hardin and Rocke (2005)
ch99AsymptoticDF(  50, 5)
ch99AsymptoticDF( 100,10)
ch99AsymptoticDF( 500,10)
ch99AsymptoticDF(1000,20)
}
\keyword{ robust }
\keyword{ multivariate }
