% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_outer.R
\name{lps.bma}
\alias{lps.bma}
\title{Log Predictive Score}
\usage{
lps.bma(object, realized.y, newdata = NULL)
}
\arguments{
\item{object}{an object of class \code{\link{pred.density}}, or class
\code{bma} (cf. \code{\link{bms}}), or class \code{\link{zlm}}}

\item{realized.y}{a vector with realized values of the dependent variables
to be plotted in addition to the predictive density, must have its length
conforming to \code{newdata}}

\item{newdata}{Needs to be provided if \code{object} is not of class
\code{\link{pred.density}}: a data.frame, matrix or vector containing
variables with which to predict.}
}
\value{
A scalar denoting the log predictive score
}
\description{
Computes the Log Predictive Score to evaluate a forecast based on a bma
object
}
\details{
The log predictive score is an indicator for the likelihood of several
forecasts.\cr It is defined as minus the arithmethic mean of the logarithms
of the point densities for \code{realized.y} given \code{newdata}.\cr Note
that in most cases is more efficient to first compute the predictive density
object via a call to \code{\link{pred.density}} and only then pass the
result on to \code{lps.bma}.
}
\examples{

 data(datafls)
 mm=bms(datafls,user.int=FALSE,nmodel=100)
 
 #LPS for actual values under the used data (static forecast)
 lps.bma(mm, realized.y=datafls[,1] , newdata=datafls[,-1])
 
 #the same result via predicitve.density
 pd=pred.density(mm, newdata=datafls[,-1])
 lps.bma(pd,realized.y=datafls[,1])
 
 # similarly for a linear model (not BMA)
 zz = zlm(datafls)
 lps.bma(zz, realized.y=datafls[,1] , newdata=datafls[,-1])
 
}
\seealso{
\code{\link{pred.density}} for constructing predictive densities,
\code{\link{bms}} for creating \code{bma} objects, \code{\link{density.bma}}
for plotting coefficient densities

Check \url{http://bms.zeugner.eu} for additional help.
}
\keyword{utilities}
