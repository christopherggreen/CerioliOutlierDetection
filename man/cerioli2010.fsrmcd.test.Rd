\name{cerioli2010.fsrmcd.test}
\alias{cerioli2010.fsrmcd.test}
\title{ Finite-Sample Reweighted MCD Outlier Detection Test of Cerioli (2010) }
\description{
  Given a set of observations, this function tests whether 
  there are outliers in the data set and identifies outlying
  points. Outlier testing/identification is done using the
  Mahalanobis-distances based on the MCD dispersion estimate. 
  The finite-sample reweighted MCD method of Cerioli (2010) is used 
  to test for unusually large distances.
}
\usage{
cerioli2010.fsrmcd.test(datamat, 
  mcd.alpha = max.bdp.mcd.alpha(n,v), 
  signif.alpha = 0.05, nsamp = 500, 
  nmini = 300, trace = FALSE)
}
\arguments{
  \item{datamat}{ (Data Frame or Matrix)
    Data set to test for outliers (rows = observations, 
      columns = variables). \code{datamat} cannot have
    missing values; please deal with them prior to
    calling this function. \code{datamat} will be converted
    to a matrix.
  }
  \item{mcd.alpha}{ (Numeric)
    Value to control the fraction
    of observations used to compute the determinant.
    Default value is corresponds to the maximum breakpoint
    case of the MCD; valid values are between
    0.5 and 1. See the \code{\link{covMcd}} documentation 
    in the \code{robustbase} library.
  }
  \item{signif.alpha}{ (Numeric)
    Desired nominal size of the individual outlier test 
	(e.g., 0.05). Equivalently, significance level at which 
	to test observations for outlyingness. To test the
	intersection hypothesis of no outliers in the data,
	specify 
	\eqn{\alpha = 1 - (1 - \gamma)^{(1/n)}}{\alpha = 1 - (1 - \gamma)^(1/n)}, 
	where \eqn{\gamma} is the nominal size of the test and 
	\eqn{n} is the number of observations.
  }
  \item{nsamp}{ (Integer)
    Number of subsamples to use in computing the MCD. see 
    \code{\link{covMcd}} documentation in the
    \code{robustbase} library.
  }
  \item{nmini}{ (Integer)
    see \code{\link{covMcd}} documentation in the
    \code{robustbase} library.
  }
  \item{trace}{ (Logical)
    see \code{\link{covMcd}} documentation in the
    \code{robustbase} library.
   }
}
%%  ~~ \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  \item{mu.hat}{Location estimate from the MCD calculation}
  \item{sigma.hat}{Dispersion estimate from the MCD calculation}
  \item{mahdist}{Mahalanobis distances calculated using the MCD estimate} 
  \item{DD}{Hardin-Rocke critical values for testing MCD distances. Used
    to produce weights for reweighted MCD. See Equation (16) in Cerioli (2010).}
  \item{weights}{Weights used in the reweighted MCD. See Equation (16) in Cerioli
    (2010).} 
  \item{mu.hat.rw}{Location estimate from the reweighted MCD calculation} 
  \item{sigma.hat.rw}{Dispersion estimate from the reweighted MCD calculation}
  \item{mahdist.rw}{a matrix of dimension \code{nrow(datamat)} by
    \code{length(signif.alpha)} of Mahalanobis distances computed
    using the finite-sample reweighted MCD methodology in
    Cerioli (2010). Even though the distances do not depend
    on \code{signif.alpha}, there is one column per entry in 
    \code{signif.alpha} for user convenience.
  }
  \item{critvalfcn}{Function to compute critical values for Mahalanobis distances
    based on the reweighted MCD; see Equations (18) and (19) in Cerioli (2010). The
    function takes a signifance level as its only argument.} 
  \item{signif.alpha}{Significance levels used in testing.} 
  \item{mcd.alpha}{Fraction of the observations used to compute the MCD estimate} 
  \item{outliers}{A matrix of dimension \code{nrow(datamat)} by 
    \code{length(signif.alpha)} indicating 
    whether each row of \code{datamat} is an outlier. The \code{i-th} 
    column corresponds to the result of testing observations
    for outlyingness at significance level \code{signif.alpha[i]}.}
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>
}
\references{
Andrea Cerioli. Multivariate outlier detection with high-breakdown estimators. 
Journal of the American Statistical Association, 105(489):147-156, 2010.

Andrea Cerioli, Marco Riani, and Anthony C. Atkinson. Controlling the size of multivariate
outlier tests with the mcd estimator of scatter. Statistical Computing, 19:341-353, 2009.
}
\seealso{
  \code{\link{cerioli2010.irmcd.test}}
}
\examples{

require(mvtnorm, quiet=TRUE)

############################################
# dimension v, number of observations n
v <- 5
n <- 200
simdata <- array( rmvnorm(n*v, mean=rep(0,v), 
    sigma = diag(rep(1,v))), c(n,v) )
#
# detect outliers with nominal sizes 
# c(0.05,0.01,0.001)
#
sa <- 1. - ((1. - c(0.05,0.01,0.001))^(1./n))
results <- cerioli2010.fsrmcd.test( simdata, 
    signif.alpha=sa )
# count number of outliers detected for each 
# significance level
colSums( results$outliers )


#############################################
# add some contamination to illustrate how to 
# detect outliers using the fsrmcd test
# 10/200 = 5% contamination
simdata[ sample(n,10), ] <- array( 
  rmvnorm( 10*v, mean=rep(2,v), sigma = diag(rep(1,v))),
  c(10,v)
)
results <- cerioli2010.fsrmcd.test( simdata, 
  signif.alpha=sa )
colMeans( results$outliers )


\dontrun{
#############################################
# example of how to ensure the size of the test is correct

  n.sim <- 5000
  simdata <- array( 
    rmvnorm(n*v*n.sim, mean=rep(0,v), sigma=diag(rep(1,v))),
    c(n,v,n.sim)
  )
  # in practice we'd do this using one of the parallel processing
  # methods out there
  sa <- 1. - ((1. - 0.01)^(1./n))
  results <- apply( simdata, 3, function(dm) {
    z <- cerioli2010.fsrmcd.test( dm, 
      signif.alpha=sa )
    # true if outliers were detected in the data, false otherwise
    any(z$outliers[,1,drop=TRUE])
  })
  # count the percentage of samples where outliers were detected;
  # should be close to the signif.alpha value used (0.01) in these
  # samples
  mean(results)

}
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ htest }
