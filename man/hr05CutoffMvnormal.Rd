\name{hr05CutoffMvnormal}
\alias{hr05CutoffMvnormal}
\title{ Corrected Critical Values for Testing MCD-based Mahalanobis Distances }
\description{
  Provides critical values for testing for outlyingness using MCD-based 
  Mahalanobis distances and the \eqn{F} distributional approximation 
  developed by Hardin and Rocke (2005) or the enhancement by Green and Martin (2014).
}
\usage{ hr05CutoffMvnormal(n.obs, p.dim, mcd.alpha, signif.alpha, 
  method = c("GM14", "HR05"), use.consistency.correction = FALSE)
}
\arguments{
  \item{n.obs}{(Integer) Number of observations}
  \item{p.dim}{(Integer) Dimension of the data, i.e.,
    number of variables. }
  \item{mcd.alpha}{ (Numeric) Value that  
    determines the fraction of the sample used to 
    compute the MCD estimate. Defaults to the value
    used in the maximum breakdown point case of the
    MCD.
  }
  \item{signif.alpha}{ (Numeric) Significance level for testing the 
    null hypothesis. Default value is 0.05. 
  }
  \item{method}{Either "HR05" to use the method of Hardin and Rocke (2005),
    or "GM14" to use the method of Green and Martin (2014).}
  \item{use.consistency.correction}{ (Logical) By default, the method does
    not multiply the cutoff values by the consistency correction for the MCD,
    under the assumption that the correction was applied during the calculation
    of the MCD-based Mahalanobis distances. Specify TRUE to add the correction
    factor if you need it for your application.}
}
\details{
  \code{hr05CutoffMvnormal} is the typical way in which a user will calculate
  critical values for testing outlyingness via MCD-based Mahalanobis distances.
  The critical values come from the \eqn{F} distributional approximation 
  derived by Hardin and Rocke (2005). One can use either the corrected degrees
  of freedom parameter derived in that paper (which was only shown to work for 
  the maximum breakdown point case of MCD), or the correction derived in 
  Green and Martin (2014) for arbitrary values of \code{mcd.alpha}.
}
\value{
  \item{ cutoff.pred }{ Critical value based on the predicted
    Wishart degrees of freedom \code{m.pred} }
  \item{ cutoff.asy  }{ Critical value based on the asymptotic
    Wishart degrees of freedom \code{m.asy} }
  \item{ c.alpha     }{ The value of the consistency correction
    factor, \eqn{c_{\alpha} } }
  \item{ m.asy       }{Asymptotic Wishart degrees of freedom
    parameter }
  \item{ m.pred      }{Predicted Wishart degrees of freedom 
    (using the method specified in \code{method}) }
  \item{ n.obs       }{Number of observations}
  \item{ p.dim       }{Number of variables}
}
\references{
C. G. Green and R. Douglas Martin. An extension of a method of Hardin and Rocke, with
an application to multivariate outlier detection via the IRMCD method of Cerioli.
Working Paper, 2014. Available from 
\url{http://students.washington.edu/cggreen/uwstat/papers/cerioli_extension.pdf}

J. Hardin and D. M. Rocke. The distribution of robust distances. Journal of Computational
and Graphical Statistics, 14:928-946, 2005.
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\seealso{
\code{\link{hr05CriticalValue}}, \code{\link{hr05AdjustedDF}}
}
\examples{
# examples from page 941 of Hardin and Rocke
hr05CutoffMvnormal(n.obs=50  , p.dim=5 , signif.alpha=0.05)
hr05CutoffMvnormal(n.obs=100 , p.dim=10, signif.alpha=0.05)
hr05CutoffMvnormal(n.obs=500 , p.dim=10, signif.alpha=0.05)
hr05CutoffMvnormal(n.obs=1000, p.dim=20, signif.alpha=0.05)
}
\keyword{ robust }
\keyword{ multivariate }
