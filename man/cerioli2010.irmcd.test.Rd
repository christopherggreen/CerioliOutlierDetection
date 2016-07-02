\name{cerioli2010.irmcd.test}
\alias{cerioli2010.irmcd.test}
\title{
Iterated RMCD test of Cerioli (2010)
}
\description{
  Given a set of observations, this function tests whether 
  there are outliers in the data set and identifies outlying
  points. Outlier testing/identification is done using the
  Mahalanobis-distances based on the MCD dispersion estimate. 
  The iterated reweighted MCD method of Cerioli (2010) is used 
  to ensure the test has the specified size.
}
\usage{
cerioli2010.irmcd.test(datamat, 
  mcd.alpha = max.bdp.mcd.alpha(n,v), 
  signif.gamma = 0.05, nsamp = 500, 
  nmini = 300, trace = FALSE, hrdf.method=c("GM14","HR05"))
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
  \item{signif.gamma}{ (Numeric)
    Desired nominal size of the \emph{intersection} outlier 
    test (e.g., 0.05), i.e., a test that there are no outliers 
    in the data. This is the \eqn{\gamma} parameter in 
    Cerioli (2010). The corresponding \eqn{\alpha} parameter
    for testing individual observations for outlyingness will
    be calculated from \eqn{\gamma} as 
  \deqn{ \alpha = 1 - ( 1 - \gamma )^{(1/n)}. }
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
  \item{hrdf.method}{ (String)
    Method to use for computing degrees of freedom
    and cutoff values for the non-MCD subset. The
    original method of Hardin and Rocke (2005) and
    the expanded method of Green and Martin (2014)
    are available as the options ``HR05'' and ``GM14'',
    respectively. ``GM14'' is the default, as it 
    is more accurate across a wider range of \code{mcd.alpha}
    values.
  }
}
\details{
  Calls the finite-sample reweighted MCD (FSRMCD) outlier 
  detection function \code{\link{cerioli2010.fsrmcd.test}}
  first to test for the existence of any outliers in the
  data. If the FSRMCD method rejects the null hypothesis of
  no outliers in the data, individual observations are then
  tested for outlyingness using the critical value function
  returned by \code{\link{cerioli2010.fsrmcd.test}} with 
  significance \eqn{\gamma}.
}
\value{
  \item{outliers}{A matrix of dimension \code{nrow(datamat)} by 
    \code{length(signif.gamma)} indicating 
    whether each row of \code{datamat} is an outlier. The \code{i-th} 
    column corresponds to the result of testing observations
    for outlyingness at significance level \code{signif.gamma[i]}.}
  \item{mahdist.rw}{a matrix of dimension \code{nrow(datamat)} by
    \code{length(signif.gamma)} of Mahalanobis distances computed
    using the finite-sample reweighted MCD methodology in
    Cerioli (2010). Even though the distances do not depend
    on \code{signif.gamma}, there is one column per entry in 
    \code{signif.gamma} for user convenience.
  }
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
  \code{\link{cerioli2010.fsrmcd.test}}
}
\examples{
require(mvtnorm, quiet=TRUE)

############################################
# dimension v, number of observations n
v <- 5
n <- 200
simdata <- array( rmvnorm(n*v, mean=rep(0,v), 
    sigma = diag(rep(1,v))), c(n,v) )
# detect outliers
results <- cerioli2010.irmcd.test( simdata, 
    signif.gamma=c(0.05,0.01,0.001) )
# count number of outliers detected for each 
# significance level
colSums( results$outliers )


#############################################
# add some contamination to illustrate how to 
# detect outliers using the irmcd test
# 10/200 = 5% contamination
simdata[ sample(n,10), ] <- array( 
  rmvnorm( 10*v, mean=rep(2,v), sigma = diag(rep(1,v))),
  c(10,v)
)
results <- cerioli2010.irmcd.test( simdata, 
  signif.gamma=0.01 )
mean( results$outliers[,1,drop=TRUE] )


#############################################
# example of how to ensure the size of the test is correct
\dontrun{
  n.sim <- 5000
  simdata <- array( 
    rmvnorm(n*v*n.sim, mean=rep(0,v), sigma=diag(rep(1,v))),
    c(n,v,n.sim)
  )
  # in practice we'd do this using one of the parallel processing
  # methods out there
  results <- apply( simdata, 3, function(dm) {
    z <- cerioli2010.irmcd.test( dm, 
      signif.gamma=0.01 )
    # true if outliers were detected in the data, false otherwise
    any(z$outliers[,1,drop=TRUE])
  })
  # count the percentage of samples where outliers were detected;
  # should be close to the signif.gamma value used (0.01) in these
  # samples
  mean(results)
}
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ htest }
