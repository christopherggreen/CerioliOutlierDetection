hr05AdjustedDF <- 
#
# Calculate the adjusted degrees of freedom parameter
# for the F distribution given in Hardin and Rocke (2005)
# for testing Mahalanobis distances calculated with the 
# MCD
# 
# Christopher G. Green
# 2011
#
function( mcd.alpha, m.asy, p, n, method=c("HR05","CG") ) {

  method <- match.arg(method)
  retval <- numeric(0)

  retval <- if ( method == "HR05" ) {
    if ( mcd.alpha == floor((n+p+1)/2)/n ) {
      # original equation from Hardin and Rocke 2005
      hr05.predict.050.hr05(m.asy,p,n)
    } else {
      stop("HR05 unsupported for alpha other than maximum breakdown case.")
    }
  } else if ( method == "CG" ) {
    # use fitted models to adjust asymptotic degrees of freedom to 
    # simulated values for small samples
    predictfunc( m.asy, p, n, mcd.alpha )
  }

  retval

}

# internal functions

hr05.predict.050.hr05 <- function(m.asy, p, n)  { 
  # original equation from Hardin and Rocke 2005
  m.asy * exp( 0.725  - 0.00663*p - 0.0780*log(n)) 
}

predictfunc <- function(m.asy, p, n, alpha) {
  # use fitted models to adjust asymptotic degrees of freedom to 
  # simulated values for small samples
  z1 <- 13.1072594 - 15.0079984 * alpha + 0.1347435 * p
  z2 <- n^(0.5119686 + 0.2235726 * alpha)
  m.asy * exp( z1/z2 )
}

