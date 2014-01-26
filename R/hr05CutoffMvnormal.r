hr05CutoffMvnormal <- 
# 
# corrected cutoff values, degrees of freedom, etc.
# from Hardin Rocke 2005 paper
# assumes MCD subset comes from MV normal
#
# Christopher G. Green, 2/5/2010
#
function( n.obs, p.dim, mcd.alpha=1/2, signif.alpha, 
	method=c("CG","HR05"), use.consistency.correction=FALSE ) {

	method          <- match.arg(method)

  # constants from Croux and Haesbroeck 1999 paper
#  one.minus.alpha  <- robustbase::h.alpha.n(mcd.alpha,n.obs,p.dim)/n.obs
#  q.alpha          <- qchisq(one.minus.alpha, df=p.dim)
#  p.alpha          <- pchisq(q.alpha, df=p.dim+2)
#  c.alpha          <- one.minus.alpha / p.alpha
#  c2               <- -0.5*p.alpha
#  c3               <- -0.5*pchisq(q.alpha, df=p.dim+4)
#  #c4              <- 3*c3
#  #b1              <- c.alpha * (c3 - c4 )/one.minus.alpha
#  b1               <- c.alpha * (-2. * c3) /one.minus.alpha
#  b2               <- 0.5 + (c.alpha/one.minus.alpha)*(c3 - 
#                     (q.alpha/p.dim)*(c2 + 0.5*one.minus.alpha))
#  z                <- b1 - p.dim*b2
#  v1               <- one.minus.alpha * b1 * b1 * ( 
#                     mcd.alpha * ( (c.alpha*q.alpha)/p.dim - 1. )^2 - 1.) - 
#                     2.*c3*c.alpha*c.alpha*(3.*z*z +
#                     (p.dim + 2)*b2*(b1 + z))
#
#  v2               <- n.obs*c.alpha*c.alpha*(b1*z*one.minus.alpha)^2
#  v                <- v1/v2
#  m.hat.asy        <- 2. / (c.alpha*c.alpha*v)
  ch99             <- ch99AsymptoticDF( n.obs, p.dim, mcd.alpha )
  # ch99$c.alpha, ch99$m.hat.asy

  # use hardin and rocke 2005 results or my simulated versions
  # to get predicted degrees of freedom value

  #m.pred           <- 0.725 - 0.00663*p.dim - 0.0780*log(n.obs)
  m.pred           <- hr05AdjustedDF(n.obs=n.obs, p.dim=p.dim, mcd.alpha=mcd.alpha, 
     m.asy=ch99$m.hat.asy, method=method)
  # cgg this is now handled by hr05.predict
  #m.pred           <- ch99$m.hat.asy * exp(m.pred)

  cutoff.pred      <- hr05CriticalValue(m.pred        ,p.dim,signif.alpha)
  cutoff.asy       <- hr05CriticalValue(ch99$m.hat.asy,p.dim,signif.alpha)

  if ( use.consistency.correction ) {
    cutoff.pred <- cutoff.pred * ch99$c.alpha
	 cutoff.asy  <- cutoff.asy  * ch99$c.alpha
  }

  list( cutoff.pred = cutoff.pred, cutoff.asy = cutoff.asy, c.alpha = ch99$c.alpha, 
	m.asy = ch99$m.hat.asy, m.pred = m.pred, n = n.obs, p = p.dim )
}
