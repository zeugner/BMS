% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_outer.R
\name{is.bma}
\alias{is.bma}
\title{Tests for a 'bma' Object}
\usage{
is.bma(bmao)
}
\arguments{
\item{bmao}{a 'bma' object: see 'value'}
}
\value{
Returns \code{TRUE} if bmao is of class 'bma', \code{FALSE}
otherwise.
}
\description{
tests for objects of class "bma"
}
\examples{

 data(datafls)
 mm=bms(datafls,burn=1000, iter=4000)
 is.bma(mm)
}
\seealso{
'Output' in \code{\link{bms}} for the structure of a 'bma' object

Check \url{http://bms.zeugner.eu} for additional help.
}
\keyword{classes}
