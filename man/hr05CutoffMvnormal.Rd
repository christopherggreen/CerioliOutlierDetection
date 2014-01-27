\name{hr05CutoffMvnormal}
\alias{hr05CutoffMvnormal}
\title{ Corrected Critical Values for Testing MCD-based Mahalanobis Distances }
\description{
  Provides critical values for testing for outlyingness using MCD-based 
  Mahalanobis distances and the \eqn{F} distributional approximation 
  developed by Hardin and Rocke (2005) or the enhancement by Green (2014).
}
\usage{ hr05CutoffMvnormal(n.obs, p.dim, mcd.alpha, signif.alpha, 
  method = c("CG", "HR05"), use.consistency.correction = FALSE)
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
    null hypothesis 
  }
  \item{method}{Either "HR05" to use the method of Hardin and Rocke (2005),
    or "CG" to use the method of Green (2014).}
  \item{use.consistency.correction}{ (Logical) By default, the method does
    not multiply the cutoff values by the consistency correction for the MCD,
    under the assumption that the correction was applied during the calculation
    of the MCD-based Mahalanobis distances. Specify TRUE to add the correction
    factor if you need it for your application.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
}
\references{
J. Hardin and D. M. Rocke. The distribution of robust distances. Journal of Computational
and Graphical Statistics, 14:928–946, 2005.

Christopher G. Green
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
"hr05CutoffMvnormal"
}
\keyword{ robust }
\keyword{ multivariate }
