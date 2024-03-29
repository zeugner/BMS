% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BMS-package.R
\docType{package}
\name{BMS-package}
\alias{BMS}
\alias{BMS-package}
\title{BMS: Bayesian Model Averaging Library}
\description{
Bayesian Model Averaging for linear models with a wide choice of (customizable) priors. Built-in priors include coefficient priors (fixed, hyper-g and empirical priors), 5 kinds of model priors, moreover model sampling by enumeration or various MCMC approaches. Post-processing functions allow for inferring posterior inclusion and model probabilities, various moments, coefficient and predictive densities. Plotting functions available for posterior model size, MCMC convergence, predictive and coefficient densities, best models representation, BMA comparison. Also includes Bayesian normal-conjugate linear model with Zellner's g prior, and assorted methods.
}
\details{
The key function you need is  \code{\link{bms}}.
}
\references{
\url{http://bms.zeugner.eu}: BMS package homepage with help and tutorials

Feldkircher, M. and S. Zeugner (2015): Bayesian Model Averaging Employing 
Fixed and Flexible Priors: The BMS Package for R, Journal of Statistical Software 68(4).

Feldkircher, M. and S. Zeugner (2009): Benchmark Priors
Revisited: On Adaptive Shrinkage and the Supermodel Effect in Bayesian Model
Averaging, IMF Working Paper 09/202.
}
\seealso{
\code{\link{coef.bma}}, \code{\link{plotModelsize}} and
\code{\link{density.bma}} for some operations on the resulting 'bma' object, 
as well as 
\code{\link{predict.bma}} or \code{\link{gdensity}}, or
\code{\link{zlm}} for individual Zellner regression models.

Check \url{http://bms.zeugner.eu} for additional help.
}
\author{
Martin Feldkircher, Paul Hofmarcher, and Stefan Zeugner
}
\keyword{internal}
