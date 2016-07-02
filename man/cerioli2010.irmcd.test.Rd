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
# banknote example from Cerioli (2010)
\dontrun{

  require(rrcov) # for CovMcd
  require(alr3)  # banknote data set lives here
  data(banknote, package="alr3")
  # length, width of left edge, width of right edge,
  # width of bottom edge, width of top edge, length
  # of image diagonal, counterfeit (1=counterfeit)

  bnk.gamma <- 0.01
  # genuine banknotes
  # classical mean and covariance
  banknote.real <- banknote[ banknote[,"Y"]==0, 1:6 ]
  cov.cls <- CovClassic( banknote.real  )
  # 1 - (1 - 0.01)^(1/100) quantile of scaled-Beta distribution
  # with m=100 and v=6
  bnk.m <- nrow( banknote.real )
  bnk.v <- ncol( banknote.real ) 
  bnk.alpha <- 1. - ((1. - bnk.gamma)^(1./bnk.m))
  cutoff.cls <- (bnk.m-1.)*(bnk.m-1.)*qbeta( 1. - bnk.alpha, bnk.v/2., 
	(bnk.m - bnk.v - 1.)/2.)/bnk.m
  # Figure 4 (left) in Cerioli (2010)
  plot( getDistance( cov.cls ), xlab="Index number", 
	ylab="Squared Mahalanobis Distance", type="p", 
	ylim=c(0,45)
  )
  abline( h=cutoff.cls )

  # reweighted MCD, maximum breakdown point case
  cov.rob <- CovMcd( banknote.real, 
	alpha=floor((bnk.m + bnk.v + 1.)/2.)/bnk.m, nsamp="best" )
  # cutoff using chi-squared individually
  cutoff.rmcdind <- qchisq(1. - bnk.gamma, df=bnk.v)
  # cufoff using simultaneous chi-square
  cutoff.rmcdsim <- qchisq(1. - bnk.alpha, df=bnk.v)
  # scaled-F cutoff using FSRMCD
  # cutoff value is returned by critvalfcn for observations
  # with weight=0
  tmp.fsrmcd <- cerioli2010.fsrmcd.test( banknote.real, 
	signif.alpha=bnk.alpha )
  cutoff.fsrmcd <- unique(tmp.fsrmcd$critvalfcn( bnk.alpha )[tmp.fsrmcd$weights==0])
  # Figure 4 (right)
  plot( getDistance( cov.rob ), xlab="Index number",
	ylab="Squared Robust Reweighted Distance", type="p",
	ylim=c(0,45)
  )
  abline( h=cutoff.rmcdind, lty="dotted" )
  abline( h=cutoff.rmcdsim, lty="dashed" )
  abline( h=cutoff.fsrmcd, lty="solid" )
  legend( "topright", c("RMCD_ind","RMCD","FSRMCD"), 
    lty=c("dotted","dashed","solid") )

  # forged banknotes
  # classical mean and covariance
  banknote.fake <- banknote[ banknote[,"Y"]==1, 1:6 ]
  cov.cls <- CovClassic( banknote.fake  )
  # 1 - (1 - 0.01)^(1/100) quantile of scaled-Beta distribution
  # with m=100 and v=6
  bnk.m <- nrow( banknote.fake )
  bnk.v <- ncol( banknote.fake ) 
  bnk.alpha <- 1. - ((1. - bnk.gamma)^(1./bnk.m))
  cutoff.cls <- (bnk.m-1.)*(bnk.m-1.)*qbeta( 1. - bnk.alpha, bnk.v/2., 
  	(bnk.m - bnk.v - 1.)/2.)/bnk.m
  # Figure 5 (left) in Cerioli (2010)
  plot( getDistance( cov.cls ), xlab="Index number", 
  	ylab="Squared Mahalanobis Distance", type="p", 
  	ylim=c(0,45)
  )
  abline( h=cutoff.cls )
  
  
  # reweighted MCD, maximum breakdown point case
  cov.rob <- CovMcd( banknote.fake, 
  	alpha=floor((bnk.m + bnk.v + 1.)/2.)/bnk.m, nsamp="best" )
  # cutoff using chi-squared individually
  cutoff.rmcdind <- qchisq(1. - bnk.gamma, df=bnk.v)
  # scaled-F cutoff using FSRMCD
  # cutoff value is returned by critvalfcn for observations
  # with weight=0
  tmp.fsrmcd <- cerioli2010.fsrmcd.test( banknote.fake, 
  	signif.alpha=bnk.alpha )
  cutoff.fsrmcd <- unique(tmp.fsrmcd$critvalfcn( bnk.alpha )[tmp.fsrmcd$weights==0])
  cutoff.irmcd <- unique(tmp.fsrmcd$critvalfcn( bnk.gamma )[tmp.fsrmcd$weights==0])
  # Figure 5 (right) in Cerioli (2010)
  plot( getDistance( cov.rob ), xlab="Index number",
  	ylab="Squared robust reweighted Distance", type="p",
  	ylim=c(0,150)
  )
  abline( h=cutoff.rmcdind, lty="dotted" )
  abline( h=cutoff.fsrmcd, lty="dashed" )
  abline( h=cutoff.irmcd, lty="solid" )
  legend( "topright", c("RMCD_ind","FSRMCD","IRMCD"), 
  	lty=c("dotted","dashed","solid") )

}

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
