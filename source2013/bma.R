######################################
# Bayesian Model Averaging Program   #
######################################
# This version: adjusted on 2011-05-05
# Martin Feldkircher
# martin.feldkircher@gzpace.net, http://feldkircher.gzpace.net
# Stefan Zeugner
# stefan.zeugner@ulb.ac.be, http://www.zeugner.eu
#####################
# The main code starts at line 169 with the function "bms=function(....)" and is written
# by Martin Feldkircher and Stefan Zeugner as part of their work at the Institute for Advanced Studies (IHS),
# Vienna in 2006/2007. Descriptions of the algorithms and priors used can be found in
# Gary Koop ("Bayesian Econometrics", Wiley & Sons), Fernandez, C., E. Ley and M.F.J. Steel (2001b) 
# ?Model Uncertainty in CrossCountry Growth Regressions,? Journal of Applied Econometrics,
# Fernandez, C., E. Ley and M.F.J. Steel (2001a) ?Benchmark Priors for Bayesian Model Averaging,?
# Journal of Econometrics, 100: 381?427. and 
# Liang, F., R. Paulo, G. Molina, M.A. Clyde, and J.O. Berger (2008): "Mixtures of g Priors for Bayesian Variable Selection",
# Journal of the American Statistical Association, 103: 410-423

####################
# USAGE            # 
###################################################################################################
# bms <-function(X.data,burn=1000,iter=NA,nmodel=100,mcmc="bd",g="UIP",mprior="random",mprior.size=NA,user.int=TRUE,
#                start.value=NA,g.stats=TRUE,logfile=FALSE,logstep=10000,force.full.ols=FALSE)
#
######################
# FUNCTION ARGUMENTS #
###################################################################################################
#X.data       data.frame or matrix: submit a data frame or a matrix, where the first column corresponds to the dependent variable
#             followed by the covariates, including a constant term is not necessary since y and X is 
#             demeaned automatically.
#burn         integer >=0: is the number of burn-in draws, default 1000    (not taken into account if mcmc="enumerate")
#iter         integer >=0: is the number of posterior draws, default 3000  (not taken into account if mcmc="enumerate"); 
#             if mcmc="enumerate", then iter is the number of models to be sampled, starting from 0 (default 2^K-1) - cf. start.value
#nmodel       integer >=0: is the number of best models for which information is stored. Convergence analysis of
#             the sampler by means of the correlation between analytic posterior  model probabilities
#             (pmp's) and hat of the MCMC sampler is based on the number you  have set in nmodels. Also
#             if you want to save the regression coefficients (beta.save=T), they are taken from the
#             nmodel best models. Setting nmodel500 slows down the MCMC sampler. Note that posterior
#             inclusion probabilites, and mean calculations are based on the MCMC frequencies as opposed
#             to exact analytical calculations (as in Fernandez, Ley and Steel).
#             Set nmodel=0 to speed up sampling (if topmodel information is not needed)
#mcmc         character: * MC3 samplers: default is "bd" which corresponds to a birth / death MCMC alogrithm. You can choose
#             also "rev.jump" where we have implemented a true reversible jump algorithm adding a "move" step to the birth / death steps from "bd".
#             * enumeraton sampler: moreover there is mcmc="enumerate"/"enum" which will enumerate all possible regressor combinations (not advisable for K>23) 
#             * interaction sampler: adding an ".int" to an MC3 sampler (e.g. "mcmc="bd.int") does an interaction sampler
#             interaction terms will only be sampled along with their component variables, in the colnumn names of X.data interaction terms need to be 
#             denominated by names consisting of the base terms separated by "#" (e.g. an interaction term of base variables "A", "B" and "C" needs column name "A#B#C")
#g            character or integer: is the hyperparameter on Zellner's g-prior for the regression coefficients. You can specify g="UIP", corresponding to g=N (default)
#             g="bric" corresponding to the benchmark prior suggestion from FLS (2001), i.e   g=max(N, K^2)
#             with K denoting the total number of regressors and N the number of observations, g="EBL" estimates a 
#             local empirical Bayes g-parameter (as in Liang et a. (2008), JASA). g="hyper", takes the 'hyper-g' 
#             prior distribution (As in Liang et al.) with the default hyper-parameter a=3; This hyperparameter can 
#             be adjusted (between 2<a<=4) by setting g="hyper=2.9", for instance.
#mprior       regards the prior on model size. It can be either "random" (Default), fixed", "uniform", "customk", or "pip" and is based on the 
#             working paper "On the Effect of Prior Assumptions in Bayesian Model Averaging with Applications
#             to Growth Regression", by  Ley and Steel (2008). Mprior denotes the a priori inclusion probability
#             of a regressor. Their suggestion is to use a binomial-beta hyperprior on mprior (i.e. mprior="random")
#             in order to be noninformative on model size. You can use mprior="fixed" if you have strong prior
#             information on model size. mprior="uniform" employs the uniform model prior; 
#             mprior="customk" allows for custom model size priors (cf. mprior.size)
#             mprior="pip" allows for custom prior inclusion probabilities (cf mprior.size)
#             Note that the prior on models with more than N-3 regressors is automatically zero; these models will not be sampled.
#mprior.size  corresponds to the expected value of the model size prior (default K/2). For mprior=random there is little
#             impact on results by varying mprior.size. For fixed mprior (i.e. informative) prior this is a
#             sensible choice. 
#             if mprior="customk" then a custom model size prior can be provided as K+1 vector detailing the priors from model size 0 to K 
#             (e.g. rep(1,K+1) for the uniform model prior)
#             if mprior="pip" then custom prior inclusion probabilities can be provided as a vector of size K, with elements in the interval (0,1] 
#user.int     'interactive mode' print out results to command line after ending the routine and do two charts
#start.value  specifies the starting model. You can choose either a specific model by the corresponding
#             column indices (e.g. starting.model=numeric(K) starts from the null model including 
#             solely a constant term) or you set a number (e.g. starting.model=20). In the latter case
#             randomly 20 covariates are chosen and the starting model is identified by those regressors
#             with a t-statistics>0.2. 
#             The default value start.value=NA corresponds to start.value=min(nrow(X)-2,ncol(X),nrow(X)-3)
#             start.value=0 or start.value=NULL starts from the null model
#             * if mcmc="enumerate" then start.value is the index to start the iteration (default=0) . Any number between 0 and K^2-1 is admissible
#g.stats      TRUE or FALSE whether statistics on the shrinkage factor should be collected. default=TRUE
#             set g.stats=FALSE for faster iteration
#logfile      setting logfile=TRUE produces a logfile named "test.log" in your current working directory,
#             in order to keep track of the sampling procedure. 
#             setting logfile equal to some filepath (like logfile="subfolder/bla.txt") puts the logfile 
#             into that specified position. (default: logfile=FALSE)
#             Note that logfile="" implies log printouts on the console
#logstep      specifies at which number of posterior draws information is written to the log file; default: 100 000 iterations
#force.full.ols  default FALSE. If force.full.ols=TRUE, the ols estimation part of the sampling procedure relies on slower matrix inversion, 
#             instead of streamlined routines. force.full.ols=TRUE can slow down sampling but may deal better with highly collinear data
#beta.save    obsolete
#int          obsolete, see "mcmc"
#exact        obsolete: see estimates.bma()
#ask.set      obsolete, with no replacement
#printRes     obsolete; see user.int
#return.g.stats renamed intop 'g.stats'
#theta        renamed into 'mprior'
#prior.msize  renamed into 'mprior.size'
#
################
# OUTPUT       #
###################################################################################################
# A call to bms returns a list/"bma" object with the following elements:
#info  a list containing some posterior statistics. it has the following elements
#   iter      (numeric(1)) The number of iteration runs in the sampler (without burn-ins)
#   burn      (numeric(1)) The number of burn-in   iterations
#   inccount  (numeric(K)) The (weighted) cumulative sum of the binary vector with regressor inclusions
#   models.visited (numeric(1)) the number of model candidates that have been accepted (including burn-ins)
#   b1mo      ((numeric(K)) the (weighted) cumulative sum of first posterior moment of beta: Sum p(M)E(beta|X,y,M)
#   b2mo      ((numeric(K)) the (weighted) cumulative sum of second posterior moment of beta: Sum p(M)E(beta^2|X,y,M)
#   add.otherstats (numeric(0 or some integer)) cumulative sum of some additional statistics (such as posterior moments of g)
#   cumsumweights (numeric(1)) the denominator to turn cumulative sums into (weighted) averages
#   K          (numeric(1)) number of covariates in X.data
#   N          (numeric(1)) number of observations in X.data
#   corr.pmp   (numeric(1)) correlation between analytical and MCMC frequencies for the best nmodel models; if mcmc="enum", then this is NA
#   msize      (numeric(1)) the cumulative sum of model sizes
#   timed      (difftime) time taken for the main sampling routine in seconds
#   k.vec      (numeric(K+1)) cumulative sum of the different model siyes, from zero to K
#   cons       (numeric(1)) scalar value of E(constant|X,y), same as BMAOBJECT$cons
#   pos.sign   (numeric(K)) cumulative sum of the coefficients>0, from 1 to K
# arguments    named list with the all the function arguments to the function bms, after being adjusted for inconsistencies
# topmod       object/list containing the best nmodel models; for mor information see help on function .top10
# start.pos    (numeric(<=K)) the indexes of covariates in the actual starting model, that started the MCMC sampling chain
# gprior.info  (list) A list containing information on the gprior: its type ("BRIC", "EBL", "hyper" or "numeric"), whether it is constant, and its 
#                     first moment (in case of the hyper-prior, also the second moment)
# X.data       (data.frame) the data in the bms function argument X.data; the same as BMAOBJECT$arguments$X.data, retained for backward compatibility
# reg.names    (character(K)) the covariates' column names (generic ones, if no names were provided with X.data)
# bms.call     (call) the function call to bms() as it was entered into the command line (regularized by match.call())
#
###############################
# Deprecated output elements  #
###################################################################################################
# in order to convert your bma object into a bma object of versions before this version: 20090612,
# use the function as.oldbma(BMAOBJECT)
#     ----- The following are deprecated elements of BMAOBJECT =bms(...) -------
# estimates   use estimates.bma(BMAOBJECT)
# estimates when exact=TRUE: use estimates.bma(BMAOBJECT, exact=TRUE)
# info        use info.bma(BMAOBJECT)
# topmodels   use topmodels.bma(BMAOBJECT)
# beta.draws  use beta.draws.bma(BMAOBJECT)
# pmp.10      use pmp.bma(BMAOBJECT)
#
#
#

###### BMS MAIN FUNCTION #########################









#' Bayesian Model Sampling and Averaging
#' 
#' Given data and prior information, this function samples all possible model
#' combinations via MC3 or enumeration and returns aggregate results.
#' 
#' Ad \code{mcmc}: \cr Interaction sampler: adding an ".int" to an MC3 sampler
#' (e.g. "mcmc="bd.int") provides for special treatment of interaction terms.
#' Interaction terms will only be sampled along with their component variables:
#' In the colnumn names of X.data, interaction terms need to be denominated by
#' names consisting of the base terms separated by \code{#} (e.g. an
#' interaction term of base variables \code{"A"}, \code{"B"} and \code{"C"}
#' needs column name \code{"A#B#C"}). Then variable \code{"A#B#C"} will only be
#' included in a model if all of the component variables ("A", "B", and "C")
#' are included.
#' 
#' The MC3 samplers "\code{bd}", "\code{rev.jump}", "\code{bd.int}" and
#' "\code{rev.jump.int}", iterate away from a starting model by adding,
#' dropping or swapping (only in the case of rev.jump) covariates.
#' 
#' In an MCMC fashion, they thus randomly draw a candidate model and then move
#' to it in case its marginal likelihood (marg.lik.) is superior to the
#' marg.lik. of the current model.
#' 
#' In case the candidate's marg.lik is inferior, it is randomly accepted or
#' rejected according to a probability formed by the ratio of candidate
#' marg.lik over current marg.lik.  Over time, the sampler should thus converge
#' to a sensible distribution. For aggregate results based on these MC3
#' frequencies, the first few iterations are typically disregarded (the
#' 'burn-ins').
#' 
#' Ad \code{g} and the hyper-g prior: The hyper-g prior introduced by Liang et
#' al. (2008) puts a prior distribution on the shrinkage factor \eqn{g/(1+g)},
#' namely a Beta distribution \eqn{ Beta(1, 1/2-1)} that is governed by the
#' parameter \eqn{a}. \eqn{a=4} means a uniform prior distribution of the
#' shrinkage factor, while \eqn{a>2} close to 2 concentrates the prior
#' shrinkage factor close to one. \cr The prior expected value is
#' \eqn{E(g/1+g)) = 2/a}. In this sense \code{g="hyper=UIP"} and
#' \code{g="hyper=BRIC"} set the prior expected shrinkage such that it conforms
#' to a fixed UIP-g (eqng=N) or BRIC-g (\eqn{g=max(K^2,N)} ).
#' 
#' 
#' @param X.data a data frame or a matrix, with the dependent variable in the
#' first column, followed by the covariates (alternatively, \code{X.data} can
#' also be provided as a \code{\link{formula}}).  Note that \code{bms}
#' automatically estimates a constant, therefore including constant terms is
#' not necessary.
#' @param burn The (positive integer) number of burn-in draws for the MC3
#' sampler, defaults to 1000. (Not taken into account if mcmc="enumerate")
#' @param iter If mcmc is set to an MC3 sampler, then this is the number of
#' iteration draws to be sampled (ex burn-ins), default 3000 draws. \cr If
#' \code{mcmc="enumerate"}, then iter is the number of models to be sampled,
#' starting from 0 (defaults to \eqn{2^K-1}) - cf. \code{start.value}.
#' @param nmodel the number of best models for which information is stored
#' (default 500). Best models are used for convergence analysis between
#' likelihoods and MCMC frequencies, as well as likelihood-based inference.\cr
#' Note that a very high value for \code{nmodel} slows down the sampler
#' significantly. Set nmodel=0 to speed up sampling (if best model information
#' is not needed).
#' @param mcmc a character denoting the model sampler to be used.\cr The MC3
#' sampler \code{mcmc="bd"} corresponds to a birth/death MCMC algogrithm.
#' \code{mcmc="rev.jump"} enacts a reversible jump algorithm adding a "swap"
#' step to the birth / death steps from "bd".\cr Alternatively, the entire
#' model space may be fully enumerated by setting \code{mcmc="enumerate"} which
#' will iterate all possible regressor combinations (Note: consider that this
#' means \eqn{2^K} iterations, where K is the number of covariates.)\cr Default
#' is full enumeration (\code{mcmc="enumerate"}) with less then 15 covariates,
#' and the birth-death MC3 sampler (\code{mcmc="bd"}) with 15 covariates or
#' more. Cf. section 'Details' for more options.
#' @param g the hyperparameter on Zellner's g-prior for the regression
#' coefficients.\cr \code{g="UIP"} corresponds to \eqn{g=N}, the number of
#' observations (default);\cr \code{g="BRIC"} corresponds to the benchmark
#' prior suggested by Fernandez, Ley and Steel (2001), i.e \eqn{g=max(N, K^2)},
#' where K is the total number of covariates;\cr \code{g="RIC"} sets
#' \eqn{g=K^2} and conforms to the risk inflation criterion by George and
#' Foster (1994)\cr \code{g="HQ"} sets \eqn{g=log(N)^3} and asymptotically
#' mimics the Hannan-Quinn criterion with \eqn{C_{HQ}=3} (cf. Fernandez, Ley
#' and Steel, 2001, p.395)\cr \code{g="EBL"} estimates a local empirical Bayes
#' g-parameter (as in Liang et al. (2008));\cr \code{g="hyper"} takes the
#' 'hyper-g' prior distribution (as in Liang et al., 2008) with the default
#' hyper-parameter \eqn{a} set such that the prior expected shrinkage factor
#' conforms to 'UIP';\cr This hyperparameter \eqn{a} can be adjusted (between
#' \eqn{2<a<=4}) by setting \code{g="hyper=2.9"}, for instance.\cr
#' Alternatively, \code{g="hyper=UIP"} sets the prior expected value of the
#' shrinkage factor equal to that of UIP (default), \code{g="hyper=BRIC"} sets
#' it according to BRIC \cr cf section 'Details' fro more on the hyper-g prior
#' @param mprior a character denoting the model prior choice, defaulting to
#' "random":\cr \code{mprior="fixed"} denotes fixed common prior inclusion
#' probabilities for each regressor as e.g. in Sala-i-Martin, Doppelhofer, and
#' Miller(2004) - for their fine-tuning, cf. \code{mprior.size}. Preferable to
#' \code{mcmc="random"} if strong prior information on model size exists;\cr
#' \code{mprior="random"} (default) triggers the 'random theta' prior by Ley
#' and Steel (2008), who suggest a binomial-beta hyperprior on the a priori
#' inclusion probability;\cr \code{mprior="uniform"} employs the uniform model
#' prior;\cr \code{mprior="customk"} allows for custom model size priors (cf.
#' \code{mprior.size});\cr \code{mprior="pip"} allows for custom prior
#' inclusion probabilities (cf. \code{mprior.size});\cr Note that the prior on
#' models with more than N-3 regressors is automatically zero: these models
#' will not be sampled.
#' @param mprior.size if \code{mprior} is "fixed" or "random",
#' \code{mprior.size} is a scalar that denotes the prior expected value of the
#' model size prior (default K/2).\cr If \code{mprior="customk"} then a custom
#' model size prior can be provided as a K+1 vector detailing the priors from
#' model size 0 to K (e.g. rep(1,K+1) for the uniform model prior);\cr if
#' \code{mprior="pip"}, then custom prior inclusion probabilities can be
#' provided as a vector of size K, with elements in the interval (0,1)
#' @param user.int 'interactive mode': print out results to console after
#' ending the routine and plots a chart (default TRUE).
#' @param start.value specifies the starting model of the iteration chain. For
#' instance a specific model by the corresponding column indices (e.g.
#' starting.model=numeric(K) starts from the null model including solely a
#' constant term) or \code{start.value=c(3,6)} for a starting model only
#' including covariates 3 and 6.\cr If \code{start.model} is set to an integer
#' (e.g. \code{start.model=15}) then that number of covariates (here: 15
#' covariates) is randomly chosen and the starting model is identified by those
#' regressors with an OLS t-statistic>0.2.\cr The default value
#' \code{start.value=NA} corresponds to
#' \code{start.value=min(ncol(X.data),nrow(X.data)-3)}. Note that
#' \code{start.value=0} or \code{start.value=NULL} starts from the null
#' model.\cr If \code{mcmc="enumerate"} then \code{start.value} is the index to
#' start the iteration (default: 0, the null model) . Any number between 0 and
#' \eqn{K^2-1} is admissible.
#' @param g.stats \code{TRUE} if statistics on the shrinkage factor g/(1+g)
#' should be collected, defaulting to TRUE (Note: set \code{g.stats=FALSE} for
#' faster iteration.)
#' @param logfile setting \code{logfile=TRUE} produces a logfile named
#' \code{"test.log"} in your current working directory, in order to keep track
#' of the sampling procedure. \code{logfile} equal to some filepath (like
#' \code{logfile="subfolder/log.txt"}) puts the logfile into that specified
#' position. (default: \code{logfile=FALSE}). Note that \code{logfile=""}
#' implies log printouts on the console.
#' @param logstep specifies at which number of posterior draws information is
#' written to the log file; default: 10 000 iterations
#' @param force.full.ols default FALSE. If \code{force.full.ols=TRUE}, the OLS
#' estimation part of the sampling procedure relies on slower matrix inversion,
#' instead of streamlined routines. \code{force.full.ols=TRUE} can slow down
#' sampling but may deal better with highly collinear data
#' @param fixed.reg indices or variable names of \code{X.data} that are fixed
#' regressors to be always included in every sampled model. Note: the parameter
#' \code{mprior.size} refers to prior model size including these fixed
#' regressors.
#' @return A list of class \code{bma}, that may be displayed using e.g.
#' \code{\link{summary.bma}} or \code{\link{coef.bma}}. The list contains the
#' following elements: \item{info}{a list of aggregate statistics: \code{iter}
#' is the number of iterations, \code{burn} the number of burn-ins.\cr The
#' following have to be divided by \code{cumsumweights} to get posterior
#' expected values: \code{inccount} are the posterior inclusion probabilities,
#' \code{b1mo} and \code{b2mo} the first and second moment of coefficients,
#' \code{add.otherstats} other statistics of interest (typically the moments of
#' the shrinkage factor), \code{msize} is the post. expected model size,
#' \code{k.vec} the posterior model size distribution, \code{pos.sign} the
#' unconditional post. probability of positive coefficients, \code{corr.pmp} is
#' the correlation between the best models' MCMC frequencies and their marg.
#' likelihoods.\cr \code{timed} is the time that was needed for MCMC sampling,
#' \code{cons} is the posterior expected value of the constant. \code{K} and
#' \code{N} are the maximum number of covariates and the sample size,
#' respectively.} \item{arguments}{a list of the evaluated function arguments
#' provided to \code{bms} (see above)} \item{topmod}{a 'topmod' object
#' containing the best drawn models. see \code{\link{topmod}} for more details}
#' \item{start.pos}{the positions of the starting model. If bmao is a'bma'
#' object this corresponds to covariates bmao$reg.names[bmao$start.pos]. If
#' bmao is a chain that resulted from several starting models (cf.
#' \code{\link{c.bma}}, then \code{start.pos} is a list detailing all of them.}
#' \item{gprior.info}{a list of class \code{\link{gprior-class}}, detailing
#' information on the g-prior: \code{gtype} corresponds to argument \code{g}
#' above, \code{is.constant} is FALSE if \code{gtype} is either "hyper" or
#' "EBL", \code{return.g.stats} corresponds to argument \code{g.stats} above,
#' \code{shrinkage.moments} contains the first and second moments of the
#' shrinkage factor (only if \code{return.g.stats==TRUE}), \code{g} details the
#' fixed g (if \code{is.constant==TRUE}), \code{hyper.parameter} corresponds to
#' the hyper-g parameter \eqn{a} as in Liang et al. (2008) }
#' \item{mprior.info}{a list of class \code{\link{mprior-class}}, detailing
#' information on the model prior: \code{origargs} lists the original arguments
#' to \code{mprior} and \code{mprior.size} above; \code{mp.msize} denotes the
#' prior mode size; \code{mp.Kdist} is a (K+1) vector with the prior model size
#' distribution from 0 to K} \item{X.data}{data.frame or matrix: corresponds to
#' argument \code{X.data} above, possibly cleaned for NAs}
#' \item{reg.names}{character vector: the covariate names to be used for X.data
#' (corresponds to \code{\link{variable.names.bma}} } \item{bms.call}{the
#' original call to the \code{bms} function}
#' @note There are several ways to speed-up sampling: \code{nmodel=10} saves
#' only the ten best models, at most a marginal improvement. \code{nmodels=0}
#' does not save the best (500) models, however then posterior convergence and
#' likelihood-based inference are not possible.  %\code{beta.save=FALSE} saves
#' the best models, but not their coefficients, which renders the use of
#' \code{image.bma} and the paramer \code{exact=TRUE} in functions such as
#' \code{coef.bma} infeasible.  \code{g.stats=FALSE} saves some time by not
#' retaining the shrinkage factors for the MC3 chain (and the best models).
#' \code{force.fullobject=TRUE} in contrast, slows sampling down significantly
#' if \code{mcmc="enumerate"}.
#' @section Theoretical background: The models analyzed are Bayesian
#' normal-gamma conjugate models with improper constant and variance priors
#' akin to Fernandez, Ley and Steel (2001): A model \eqn{M} can be described as
#' follows, with \eqn{\epsilon} ~ \eqn{N(0,\sigma^2 I)}: \deqn{latex}{ y=
#' \alpha + X \beta + \epsilon} \deqn{f(\beta | \sigma, M, g) ~ N(0, g \sigma^2
#' (X'X)^-1) }
#' 
#' Moreover, the (improper) prior on the constant \eqn{f(\alpha)} is put
#' proportional to 1. Similarly, the variance prior \eqn{f(\sigma)} is
#' proportional to \eqn{1/\sigma}.
#' @author Martin Feldkircher and Stefan Zeugner
#' @seealso \code{\link{coef.bma}}, \code{\link{plotModelsize}} and
#' \code{\link{density.bma}} for some operations on the resulting 'bma' object,
#' \code{\link{c.bma}} for integrating separate MC3 chains and splitting of
#' sampling over several runs.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @references Feldkircher, M. and S. Zeugner (2009): Benchmark Priors
#' Revisited: On Adaptive Shrinkage and the Supermodel Effect in Bayesian Model
#' Averaging, IMF Working Paper 09/202.
#' 
#' Fernandez, C. E. Ley and M. Steel (2001): Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics 100(2), 381--427
#' 
#' Ley, E. and M. Steel (2008): On the Effect of Prior Assumptions in Bayesian
#' Model Averaging with Applications to Growth Regressions. working paper
#' 
#' Liang, F., Paulo, R., Molina, G., Clyde, M. A., and Berger, J. O. (2008).
#' Mixtures of g Priors for Bayesian Variable Selection. Journal of the
#' American Statistical Association 103, 410-423.
#' 
#' Sala-i-Martin, X. and G. Doppelhofer and R.I. Miller (2004): Determinants of
#' long-term growth: a Bayesian averaging of classical estimates (BACE)
#' approach. American Economic Review 94(4), 813--835
#' @keywords models
#' @examples
#' 
#'   data(datafls)
#'   #estimating a standard MC3 chain with 1000 burn-ins and 2000 iterations and uniform model priors
#'   bma1 = bms(datafls,burn=1000, iter=2000, mprior="uniform")
#' 
#'   ##standard coefficients based on exact likelihoods of the 100 best models:
#'   coef(bma1,exact=TRUE, std.coefs=TRUE) 
#'   
#'   #suppressing user-interactive output, using a customized starting value, and not saving the best 
#'   #  ...models for only 19 observations (but 41 covariates)
#'   bma2 = bms(datafls[20:39,],burn=1000, iter=2000, nmodel=0, start.value=c(1,4,7,30),
#'      user.int=FALSE)
#'   coef(bma2)
#'   
#'   #MC3 chain with a hyper-g prior (custom coefficient a=2.1), saving only the 20 best models, 
#'   # ...and an alternative sampling procedure; putting a log entry to console every 1000th step
#'   bma3 = bms(datafls,burn=1000, iter=5000, nmodel=20, g="hyper=2.1", mcmc="rev.jump",
#'       logfile="",logstep=1000)
#'   image(bma3) #showing the coefficient signs of the 20 best models
#'   
#'   #enumerating with 10 covariates (= 1024 models), keeping the shrinkage factors 
#'   #  ...of the best 200 models
#'   bma4 = bms(datafls[,1:11],mcmc="enumerate",nmodel=200,g.stats=TRUE)
#' 
#'   #using an interaction sampler for two interaction terms
#'   dataint=datafls
#'   dataint=cbind(datafls,datafls$LifeExp*datafls$Abslat/1000,
#'         datafls$Protestants*datafls$Brit-datafls$Muslim)
#'   names(dataint)[ncol(dataint)-1]="LifeExp#Abslat"
#'   names(dataint)[ncol(dataint)]="Protestants#Brit#Muslim"
#'   bma5 = bms(X.data=dataint,burn=1000,iter=9000,start.value=0,mcmc="bd.int") 
#'   
#'   density(bma5,reg="English") # plot posterior density for covariate "English"
#'   
#'   # a matrix as X.data argument
#'   bms(matrix(rnorm(1000),100,10))
#'   
#'   # keeping a set of fixed regressors:
#'   bms(datafls, mprior.size=7, fixed.reg = c("PrScEnroll", "LifeExp", "GDP60"))
#'   # Note that mprior.size=7 means prior model size of 3 fixed to 4 'uncertain' regressors
#'   
#'   
bms <-function(X.data,burn=1000,iter=NA,nmodel=500,mcmc="bd",g="UIP",mprior="random",mprior.size=NA,user.int=TRUE,
                start.value=NA,g.stats=TRUE,logfile=FALSE,logstep=10000,force.full.ols=FALSE, fixed.reg=numeric(0)) {
#                beta.save=TRUE,exact=NA,int=NA,printRes=NA,ask.set=NA,return.g.stats=NA,theta=NULL,prior.msize=NULL #deprecated function arguments, retained for compatibility with older versions
                 


### getting data dimensions ####################
  if (class(X.data)[[1]]=="formula") { X.data=stats::model.frame(X.data); if (!is.null(ncol(X.data[[2]]))) X.data=cbind(X.data[[1]],X.data[[2]][,-1]) }

  if (any(is.na(X.data))) {
     X.data=na.omit(X.data)
     if (nrow(X.data)<3) {stop("Too few data observations. Please provide at least three data rows without NA entries.") }
     warning("Argument 'X.data' contains NAs. The corresponding rows have not been taken into account.")
  }                                                  

  
  N<-nrow(X.data)
  K=ncol(X.data)-1
  maxk=N-3  #maximum number of admissible k per model

  
  
############################################################################################################################################
#### User Checks: ##########################################################################################################################
############################################################################################################################################

                
  # check for deprecated arguments
#  if (!is.na(exact)) { warning("Function argument 'exact' has been deprecated, please refer to function 'estimates.bma' instead.") }
#  if (!is.na(int)) { mcmc=paste(mcmc,".int",sep=""); warning("Function argument 'int' has been deprecated, please add an 'int' to the argument 'mcmc' instead.") }
#  if (!is.na(printRes)) { user.int=printRes; warning("Function argument 'printRes' has been deprecated, please refer to the argument 'user.int' instead.") }  
#  if (!is.na(ask.set)) { warning("Function argument 'ask.set' has been deprecated, with no replacement.") }    
#  if (!is.na(return.g.stats)) { g.stats=return.g.stats; warning("Function argument 'return.g.stats' has been renamed into 'g.stats'.") }      
#  if (!is.null(theta)) { mprior=theta; warning("Function argument 'theta' has been renamed into 'mprior'.") }  
#  if (!is.null(prior.msize)) { mprior.size=prior.msize; warning("Function argument 'prior.msize' has been renamed into 'mprior.size'.") }      
#  return.g.stats=g.stats; #theta=mprior; prior.msize=mprior.size
  
  if (is.null(nmodel[1])||is.na(nmodel[1])||nmodel[1]<=0) {dotop=FALSE;nmodel=0} else {dotop=TRUE}
###########################################################################

    nameix=1:K; names(nameix)=colnames(X.data[,-1,drop=FALSE]); fixed.pos=nameix[fixed.reg]; rm(nameix)

######################################################################################################################################      
    #assign the sampling procedure                                                                                              
   if (missing(mcmc)&&((K-length(fixed.pos))<15)) {mcmc="enum"} 

   int=FALSE; is.enum=FALSE #int: whether interaction sampler is wanted; is.enum: whether the sampler is enumeration
   if (is.function(mcmc)) { 
      samplingfun=mcmc
      mcmc="custom"
   } else {
    if (length(grep("int",mcmc,ignore.case=TRUE))) {int=TRUE}   
   
    if (length(grep("enum",mcmc,ignore.case=TRUE)))  {
      is.enum=TRUE; samplingfun=.iterenum
      if (K>maxk) samplingfun=.iterenum.KgtN
    } else if(length(grep("bd",mcmc,ignore.case=TRUE))){
      samplingfun=switch(int+1,.fls.samp,.fls.samp.int)
    } else {
      samplingfun=switch(int+1,.rev.jump,.rev.jump.int)
    }
   }
   if (int&&(length(fixed.pos)>0L)) { warning("interaction sampler does not allow for non-zero argument fixed.pos; consequently it was set fixed.pos=0"); fixed.pos=numeric(0); }
    sampling=.fixedset.sampler(samplingfun,fullK=K,fixed.pos=fixed.pos, X.data=X.data); 
######################################################################################################################################    
  # specific enumeration user checks & init
  if (is.enum) {      
     #check for a start.value index to start enumeration from seomewhere in between (and not do all possible models)  
     burn=0;  int=FALSE; mcmc="enum"; is.enum=TRUE
     tmp=.enum_startend(iter=iter, start.value=start.value, K=K, maxk=maxk, fixed.pos=fixed.pos); iter=tmp$iter; start.value=tmp$start.value
  } else {
     if (is.na(iter)) {iter=3000}; #if no enumeration and iter not given, set to default value 3000
  }

######################################################################################################################################
   # generate logfile if desired
    if(logfile!=FALSE){
        if (is.character(logfile)) {
          sfilename=logfile}
        else {
          sfilename="test.log"
        }
        if (nchar(sfilename)>0) file.create(sfilename)
        logfile=TRUE
      cat(as.character(Sys.time()),": starting loop ... \n",append=TRUE, file=sfilename)  #write one line
      if (logstep!=10000) fact=logstep else fact=max(floor((burn+iter)/100),logstep)
    }  

   


######################################################################################################################################  
######################################################################################################################################  
 # subtract mean from all regressors as in FLS
  y<-as.matrix(X.data[,1])
  X<-as.matrix(X.data[,2:ncol(X.data)])
  y<-y-mean(y)
  X<-X-matrix(colMeans(X),N,K,byrow=TRUE)


 # multiply the whole matrix stuff out before going into the simulation loops
  XtX.big=crossprod(X)
  Xty.big=crossprod(X,y)
  yty = as.vector(crossprod(y))

 # user check: whether we have to use force.full.ols
 coreig=eigen(cor(X),symmetric=TRUE,only.values=TRUE)$values
 if (!force.full.ols) { #this line was added due to feature requests
    if (sum(coreig>1e-7)<min(K,(N-1))) { force.full.ols=TRUE }
 }
 if (any(coreig[1:min(K,(N-1))] <1e-16))  { warning(paste("data seems to be rank-deficient: its rank seems to be only ", sum(coreig>1e-13))) }

######################################################################################################################################
 # for the case that X contains interaction terms
  if(int){
      if(length(grep("#",colnames(X.data),fixed=TRUE))==0) stop("Please separate column names of interaction terms by # (e.g. A#B)")
      mPlus=.constr.intmat(X,K)
  }
  else{ mPlus<-NA }                                
  
######################################################################################################################################
    #Prior Initialization
    
    #model prior: outsourced to stupid function for cleanliness
    pmplist=.choose.mprior(mprior,mprior.size,K=K,X.data=X.data,fixed.pos=fixed.pos)
    mprior=pmplist$mp.mode;
    
    #g-prior
    gprior.info = .choose.gprior(g,N,K,return.g.stats=g.stats,yty=yty,X.data=X.data) # gprior.info is a list that summarizes info about the choice of the g-prior
    lprobcalc = gprior.info$lprobcalc

    
######################################################################################################################################
#The function Starter selects randomly a start matrix and runs a
#regression. From this regression, the
#start Design matrix is that for which the t-stats are >0.2. So we
#can circumvent starting from a very bad start point.
  start.list=.starter(K,start.value,y,N=N,XtX.big=XtX.big,Xty.big=Xty.big,X=X,fixed.pos=fixed.pos)
  molddraw=start.list$molddraw; start.position=start.list$start.position
  kold=sum(molddraw)
  position=(1:K)[molddraw==1]

########################################################################################################################################
 


################################################################################################
#    Initializing                                                                              #
################################################################################################


  # initialize sampler-specific variables    ########################################
  # these are to select additional statistics (such as g)
  collect.otherstats=FALSE
  otherstats=numeric(0)
  add.otherstats=numeric(0)
  # initialize the vector for collecting the empirical shrinkage factor moments
  if (gprior.info$return.g.stats & !(gprior.info$is.constant)) { add.otherstats=gprior.info$shrinkage.moments; collect.otherstats=TRUE } 
  cumsumweights=iter
  null.lik=lprobcalc$just.loglik(ymy=yty,k=0) # calculate Likelihood for NullModel
  
  #adding up posterior stats has been outsourced to sub-functions for speed reasons
  if (collect.otherstats) {
    addup<-  function() {

      inccount <<- inccount + molddraw #PIPs
      msize<<-msize + kold   # average size of models

      #for speed reasons, iterative adding with indexing should be done in one stacked vector      
      if (kold!=0) {
        bm[c(position,K+position,2*K+kold,3*K+position)]=c(b1,b2,1,b1>0); bmo <<- bmo+bm
        #bmo is partitioned: first K entries have cum. b1 ("b1mo"), second K entries have cum. b2 ("b2mo"), third K entries have model size dist ("k.vec"), and fourth K entries are like inccount for positive betas (add up pos. sign covariate selections)
      } else {
        null.count<<-null.count+1
      }
      
      # collect e.g. estimated g-priors, etc   
      otherstats<<-lik.list[["otherstats"]]; add.otherstats<<-add.otherstats + otherstats        
    }


  } else {

    addup <- function() {

      inccount <<- inccount + molddraw #PIPs
      msize<<-msize + kold   # average size of models

      #for speed reasons, iterative adding with indexing should be done in one stacked vector      
      if (kold!=0) {
        bm[c(position,K+position,2*K+kold,3*K+position)]=c(b1,b2,1,b1>0); bmo <<- bmo+bm
        #bmo is partitioned: first K entries have cum. b1 ("b1mo"), second K entries have cum. b2 ("b2mo"), third K entries have model size dist ("k.vec"), and fourth K entries are like inccount for positive betas (add up pos. sign covariate selections)
      } else {
        null.count<<-null.count+1
      }
    }
  }
  if (is.enum) {
    cumsumweights=0
    if (collect.otherstats) {
      addup<- function() {

      weight=  exp(pmpold+lprobold-null.lik)
      inccount <<- inccount + weight*molddraw #PIPs
      msize<<-msize + weight*kold   # average size of models
      cumsumweights<<-cumsumweights+weight #denominator to get at sum of PMPs=1     

      #for speed reasons, iterative adding with indexing should be done in one stacked vector      
      if (kold!=0) {
        bm[c(position,K+position,2*K+kold,3*K+position)]=weight*c(b1,b2,1,b1>0); bmo <<- bmo+bm
        #bmo is partitioned: first K entries have cum. b1 ("b1mo"), second K entries have cum. b2 ("b2mo"), third K entries have model size dist ("k.vec"), and fourth K entries are like inccount for positive betas (add up pos. sign covariate selections)
      } else {
        null.count<<-null.count+weight
      }
      otherstats<<-lik.list[["otherstats"]]; add.otherstats<<-add.otherstats + weight*otherstats  
      }
    } else {
      addup <- function() {
      weight=  exp(pmpold+lprobold-null.lik)
      
      inccount <<- inccount + weight*molddraw #PIPs
      msize<<-msize + weight*kold   # average size of models
      cumsumweights<<-cumsumweights+weight #denominator to get at sum of PMPs=1     

      #for speed reasons, iterative adding with indexing should be done in one stacked vector      
      if (kold!=0) {
        bm[c(position,K+position,2*K+kold,3*K+position)]=weight*c(b1,b2,1,b1>0); bmo <<- bmo+bm
        #bmo is partitioned: first K entries have cum. b1 ("b1mo"), second K entries have cum. b2 ("b2mo"), third K entries have model size dist ("k.vec"), and fourth K entries are like inccount for positive betas (add up pos. sign covariate selections)
      } else {
        null.count<<-null.count+weight
      }
      }

    }
  }
  environment(addup) <- environment() 
  ##################################################################################

 
 
  ##initialize model varibales with starter model ###################################
  ols.object=.ols.terms2(positions=(1:K)[molddraw==1],yty=yty,k=kold,N,K=K,XtX.big=XtX.big,Xty.big=Xty.big) #OLS results from starter model

  lik.list=lprobcalc$lprob.all(ymy=ols.object$ymy, k=kold, bhat=ols.object$bhat, diag.inverse=ols.object$diag.inverse) #likelihood and expected values for starter model
  lprobold=lik.list$lprob
  b1=lik.list$b1new
  b2=lik.list$b2new
  
  ## calculate the posterior model probability for the first model
  pmpold=pmplist$pmp(ki=kold,mdraw=molddraw)
  ##################################################################################
  
  
  
  
  ## initialize top 10 function ####################################################
  #topmods=.top10(nmaxregressors=K,nbmodel=nmodel,bbeta=FALSE,bbeta2=FALSE,lengthfixedvec=length(add.otherstats))
  topmods=topmod(nbmodels=nmodel,nmaxregressors=K,bbeta=FALSE,lengthfixedvec=length(add.otherstats))
  if (mcmc=="enum") { try(topmods$duplicates_possible(FALSE), silent=TRUE) }
  if (dotop && (burn==0L)) topmods$addmodel(mylik=pmpold+lprobold,vec01=molddraw,fixedvec=lik.list$otherstats)
  ##################################################################################




  ## Initialize the rest  ###########################################################
  null.count=0             #number the null model has been drawn
  models.visited=0         #how often a model has been accepted (in burn-ins and iterations)
  inccount=numeric(K)      #how often the respective covariate has been included
  msize=0                  #average model size
  k.vec=numeric(K)         #how often the respective model size has been accepted
  b1mo=numeric(K)          #holds aggregate first moment of all coefficients
  ab=numeric(K)            #Initialize them here
  b2mo=numeric(K)          #holds aggregate second moment of all coefficients
  bb=numeric(K)              
  possign=inccount         # the number of times the respective coefficent has been positive
  mnewdraw=numeric(K)      #holds the binary vector denoting the proposed model
  if (force.full.ols) {candi.is.full.object=TRUE} else {candi.is.full.object=FALSE} #candi.is.full: if TRUE, standard OLS, else OLS via Frisch-Waugh tricks
  bmo=numeric(4*K); bm=bmo #common placeholder for b1mo, b2mo, k.vec and possign
  if (is.enum) { addup() } # in case the sampler is enumeration then count the starting value as well (no burn-ins)
  if (!is.finite(pmpold)) pmpold = -1e90 # this is if in case of an MCMC sampler the starting model got 0 PMP

###############################################################################################################
###############################################################################################################





#############################################################################################
    set.seed(as.numeric(Sys.time()))              #Set Seed randomly for number generator
                                                 
    t1<-Sys.time()                                #Save time before going into the loop 
###########################################################################################
#START MAIN LOOP
###########################################################################################
nrep=burn+iter; i=0;
while(i<nrep) {
      i=i+1;
      if(logfile){ if (i %% fact==0) { cat(as.character(Sys.time()),":",i,"current draw \n",append=TRUE, file=sfilename)} } #write one line  
    
##########################################################################################
#Start sampling program
###########################################################################################

      #sample a model                                                                                           
      a=sampling(molddraw=molddraw,K=K,mPlus=mPlus,maxk=maxk,oldk=kold)
      mnewdraw=a[["mnewdraw"]]; positionnew=a[["positionnew"]]; knew=length(positionnew)

      #calculate prior model prob
      pmpnew=pmplist[["pmp"]](ki=knew,mdraw=mnewdraw) # get the (log) model prior prob
      

      if (!is.enum) {
        if (int) {if (length(c(a$dropi,a$addi))>2|i<3|force.full.ols) {candi.is.full.object=TRUE} else {candi.is.full.object=FALSE}}
        #candi.is.full.object = TRUE if there were multiple regs dropped or added due to interaction terms

        if (candi.is.full.object) {
            ols.candidate = .ols.terms2(positions=positionnew,yty=yty,k=knew,N,K=K,XtX.big=XtX.big,Xty.big=Xty.big) #in case of changing interaction terms, draw the big OLS stuff
            ymy.candi =ols.candidate[["ymy"]]
        } else {
            ymy.candi=ols.object[["child.ymy"]](a$addi,a$dropi,k=knew) #if standard sampling, use Frisch-Waugh to get the new ResidSS (faster)
        }
    
        if ( (ymy.candi<0) | is.na(ymy.candi) ) stop(paste("stumbled on rank-deficient model" ))
        lprobnew = lprobcalc[["just.loglik"]](ymy=ymy.candi,k=knew) # get the log-likelihood out of the ResidSS
        
        #Now decide whether to accept candidate draw
        
        accept.candi = as.logical(log(runif(1,0,1))< lprobnew-lprobold + pmpnew-pmpold)

      } else {
        accept.candi=TRUE
        candi.is.full.object=FALSE
      }



      if(accept.candi){
          if (!candi.is.full.object) {
               #in case one has used Frisch-Waugh and the new model got accepted,
               #calculate the 'real' inverse in order not to make copying mistakes
              ols.res = ols.object[["mutate"]](addix=a$addi, dropix=a$dropi, newpos=positionnew, newk=knew)
          } else {
              ols.object = ols.candidate
              ols.res = ols.candidate[["full.results"]]()
          }

          lik.list = lprobcalc[["lprob.all"]](max(0,ols.res$ymy), knew, ols.res$bhat, ols.res$diag.inverse)


          lprobold=lik.list[["lprob"]]
          position = positionnew
          pmpold=pmpnew  #get posterior odds for new model  if accepted
          molddraw=mnewdraw
          kold=knew
          models.visited=models.visited+1 #does not account for revisiting models
      }


# Collect Posterior Draws
########################################################################
    if (i>burn){   
          b1=lik.list[["b1new"]]; b2=lik.list[["b2new"]];
          addup() #addup does iterative, cumulative sums of quantities of interest (betas, model size, etc.)
          
          # add log(lik)*p(M) to topmodels
          if (dotop) topmods[["addmodel"]](mylik=pmpold+lprobold,vec01=molddraw,fixedvec=otherstats)
    }
}
###########################################################################################
#END MAIN LOOP
###########################################################################################


###########################################################################################
   #adjust the topmod object and calculate all the betas after sampling
   #similar to havving set bbeta=T, and bbeta2=T in the call to .top10 above
   if (dotop) topmods=.topmod.as.bbetaT(topmods,gprior.info,X.data)

###########################################################################################

###########################################################################################
  
  timed<-difftime(Sys.time(),t1)

  # do aggregating calculations
  if (is.enum) {iter=iter+1; models.visited=models.visited+1}
  bmo=matrix(bmo,4,byrow=TRUE); b1mo=bmo[1,]; b2mo=bmo[2,]; k.vec=bmo[3,]; possign=bmo[4,]; rm(bmo)
  post.inf=.post.calc(gprior.info,add.otherstats,k.vec,null.count,X.data,topmods,b1mo,b2mo,iter,burn,inccount,models.visited,K,N,msize,timed,cumsumweights,mcmc,possign)
  
  result=list(info=post.inf$info,arguments=.construct.arglist(bms, environment()),topmod=topmods,start.pos=sort(start.position),gprior.info=post.inf$gprior.info,mprior.info=pmplist,X.data=X.data,reg.names=post.inf$reg.names,bms.call=try(match.call(bms,sys.call(0)),silent=TRUE))
  class(result)=c("bma")
  
  
###########################################################################################  

  # print results to console
  if(user.int){
    print(result)
    print(timed)
    plot.bma(result) # do modelsize plot
  }

  return(invisible(result))
}



