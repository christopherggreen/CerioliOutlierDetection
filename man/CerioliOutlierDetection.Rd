\name{CerioliOutlierDetection}
\alias{CerioliOutlierDetection}
\alias{CerioliOutlierDetection-package}
\docType{package}
\title{
CerioliOutlierDetection: package for implementing the Iterated Reweighted 
MCD outlier detection method of Cerioli (2010) 
}
\description{
Implements the outlier detection methodology of Cerioli (2010) based on 
Mahalanobis distances and the minimum covariance determinant (MCD) estimate 
of dispersion. Also provides critical values for testing outlyingness of
MCD-based Mahalanobis distances using the distribution approximations 
developed by Hardin and Rocke (2005) and Green and Martin (2014).
}
\details{
  \tabular{ll}{
    Package: \tab CerioliOutlierDetection\cr
    Type: \tab Package\cr
    Version: \tab 1.1.5\cr
  }

  The function \code{cerioli2010.irmcd.test()} provides the outlier detection 
  methodology of Cerioli (2010), and is probably the best place for a new user
  of this package to start. See the documentation for that function for examples.

  This package was also used to produce the results presented in Green and 
  Martin (2014). There is a companion \code{R} package, \code{HardinRockeExtension},
  that provides code that can be used to replicate the results of that paper. The
  package \code{HardinRockeExtension} is available from Christopher G. Green's 
  GitHub: \url{https://github.com/christopherggreen/HardinRockeExtension} .
}
\author{
Written and maintained by Christopher G. Green <christopher.g.green@gmail.com>,
with advice and support from Doug Martin.
}
\references{
Andrea Cerioli. Multivariate outlier detection with high-breakdown estimators. 
Journal of the American Statistical Association, 105(489):147-156, 2010.

C. G. Green and R. Douglas Martin. An extension of a method of Hardin and Rocke, with
an application to multivariate outlier detection via the IRMCD method of Cerioli.
Working Paper, 2017. Available from 
\url{http://christopherggreen.github.io/papers/hr05_extension.pdf}

J. Hardin and D. M. Rocke. The distribution of robust distances. Journal of Computational
and Graphical Statistics, 14:928-946, 2005.
}
\keyword{ package }
