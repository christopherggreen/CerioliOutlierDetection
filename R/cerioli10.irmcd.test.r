cerioli2010.irmcd.test <- 
# 
# implements the iterated reweighted MCD
# outlier detection method presented in Cerioli (2010)
#
# Author: Christopher G. Green
# Date: 2011-06-23
#
# Update 2016-05-27
# changing argument from signif.alpha to signif.gamma
# to make it more clear what the parameterization is
#
function( datamat, mcd.alpha=max.bdp.mcd.alpha(n,v), 
  signif.gamma=0.05, nsamp = 500, nmini = 300, trace=FALSE,
  delta = 0.025, hrdf.method=c("GM14","HR05")) 
{

  datamat <- as.matrix(datamat)
  if ( any(is.na(datamat)) ) 
    stop("datamat cannot have missing values.")
  n <- nrow(datamat) # number of observations
  v <- ncol(datamat) # dimension

  hrdf.method <- match.arg(hrdf.method)

  # steps 1-4: compute the FSRMCD of Cerioli (2010) 
  # compute alpha needed to ensure Type I error rate
  # of signif.gamma for the intersection test
  alpha.ind    <- 1. - ((1. - signif.gamma)^(1./n))
  fsout    <- cerioli2010.fsrmcd.test( datamat, mcd.alpha=mcd.alpha, 
    signif.alpha=alpha.ind, nsamp=nsamp, nmini=nmini, trace=trace,
	delta=delta, hrdf.method=hrdf.method)
  # test whether any point is an outlier
  n.siggam <- length(signif.gamma)
  outliers <- fsout$outliers
  for ( i in 1:n.siggam ) {
    if ( any(outliers[,i]) ) {
      # test each mahalanobis distance at the gamma[i] level
      outliers[,i] <- fsout$mahdist.rw[,i] > fsout$critvalfcn(signif.gamma[i])
    } else {
      # accept null hypothesis: no outliers
      # don't need to do anything to outliers[,i]
    }
  }

  list(outliers=outliers, mahdist.rw=fsout$mahdist.rw)
}
