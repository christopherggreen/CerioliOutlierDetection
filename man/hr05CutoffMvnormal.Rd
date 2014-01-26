\name{hr05CutoffMvnormal}
\alias{hr05CutoffMvnormal}
\title{ Hardin Rocke Cutoff }
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hr05CutoffMvnormal(n.obs, p.dim, mcd.alpha = 1/2, signif.alpha, method = c("CG", "HR05"), use.consistency.correction = FALSE)
}
\arguments{
  \item{n.obs}{(Integer) Number of observations}
  \item{p.dim}{(Integer) Dimension of the data, i.e.,
    number of variables. }
  \item{mcd.alpha}{ (Numeric) Alpha value that  
    determines the fraction of the sample used to 
  compute the MCD estimate. }
  \item{signif.alpha}{
%%     ~~Describe \code{signif.alpha} here~~
}
  \item{method}{Either "HR05" to use the method of Hardin and Rocke (2005),
    or "CG" to use the method of Green (2014).}
  \item{use.consistency.correction}{
%%     ~~Describe \code{use.consistency.correction} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
J. Hardin and D. M. Rocke. The distribution of robust distances. Journal of Computational
and Graphical Statistics, 14:928–946, 2005.
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
