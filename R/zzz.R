




#' Class "topmod"
#' 
#' An updateable list keeping the best x models it encounters in any kind of
#' model iteration
#' 
#' 
#' @name topmod-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls to
#' \code{\link{topmod}}, or indirectly by calls to \code{\link{bms}}.\cr
#' 
#' A 'topmod' object (as created by \code{topmod}) holds three basic vectors:
#' \code{lik} (for the (log) likelihood of models or similar), \code{bool()}
#' for a hexcode presentation of the model binaries (cf. \code{\link{bin2hex}})
#' and ncount() for the times the models have been drawn.\cr All these vectors
#' are sorted descendantly by \code{lik}, and are of the same length. The
#' maximum length is limited by the argument \code{nbmodels}.
#' 
#' If \code{tmo} is a topmod object, then a call to \code{tmo$addmodel} (e.g.
#' \code{tmo$addmodel(mylik=4,vec01=c(T,F,F,T))} updates the object \code{tmo}
#' by a model represented by \code{vec01} (here the one including the first and
#' fourth regressor) and the marginal (log) likelihood \code{lik} (here: 4).\cr
#' If this model is already part of \code{tmo}, then its respective
#' \code{ncount} entry is incremented by one; else it is inserted into a
#' position according to the ranking of \code{lik}.\cr In addition, there is
#' the possibility to save (the first moments of) coefficients of a model
#' (\code{betas}) and their second moments (\code{betas2}), as well as an
#' arbitrary vector of statistics per model (\code{fixed_vector}).\cr
#' @author Martin Feldkircher and Stefan Zeugner
#' @seealso \code{\link{topmod}} to create \code{topmod} objects and a more
#' detailed description,
#' \code{\link{is.topmod}} to test for this class
#' @references \url{http://bms.zeugner.eu}
#' @keywords classes
#' @examples
#' 
#'   tm= topmod(2,4,TRUE,0) #should keep a  maximum two models
#'   tm$addmodel(-2.3,c(1,1,1,1),1:4,5:8) #update with some model
#'   tm$addmodel(-2.2,c(0,1,1,1),1:3,5:7) #add another model
#'   tm$addmodel(-2.2,c(0,1,1,1),1:3,5:7) #add it again -> adjust ncount
#'   tm$addmodel(-2.5,c(1,0,0,1),1:2,5:6) #add another model
#'   
#'   #read out
#'   tm$lik()
#'   tm$ncount()
#'   tm$bool_binary()
#'   tm$betas()
#' 
setClass("topmod",representation(addmodel="function", lik="function", bool="function", ncount="function", nbmodels="function", nregs="function", betas_raw="function", betas2_raw="function", kvec_raw="function", bool_binary="function", betas="function", betas2="function", fixed_vector="function"))


#' Class "bma"
#' 
#' A list holding results from a BMA iteration chain
#' 
#' 
#' @name bma-class
#' @docType class
#' @section Objects from the Class: Objects can be created via calls to
#' \code{\link{bms}}, but indirectly also via \code{\link{c.bma}}\cr A
#' \code{bma} object is a list whose elements hold information on input and
#' output for a Bayesian Model Averaging iteration chain, such as from a call
#' to \code{\link{bms}}:
#' @author Martin Feldkircher and Stefan Zeugner
#' @seealso \code{\link{bms}} for creating \code{bma} objects,\cr or
#' \code{\linkS4class{topmod}} for the topmod object
#' @references \url{http://bms.zeugner.eu}
#' @keywords classes
#' @examples
#' 
#'  data(datafls)
#'  mm=bms(datafls)
#'  #show posterior model size
#'  print(mm$info$msize/mm$info$cumsumweights)
#'  #is the same number as in
#'  summary(mm)
#'  
#' 
setClass("bma",representation(info="list",arguments="list",topmod="topmod",start.pos="integer",gprior.info="list",mprior.info="list",X.data="data.frame",reg.names="character",bms.call="call"))



#' Class "zlm"
#' 
#' A list holding output from the Bayesian Linar Model under Zellner's g prior
#' akin to class 'lm'
#' 
#' 
#' @name zlm-class
#' @docType class
#' @section Objects from the Class: Objects can be created via calls to
#' \code{\link{zlm}}, but indirectly also via \code{\link{as.zlm}}.\cr
#' \code{\link{zlm}} estimates a Bayesian Linear Model under Zellner's g prior
#' - its output is very similar to objects of class \code{\link{lm}} (cf.
#' section 'Value')
#' @author Martin Feldkircher and Stefan Zeugner
#' @seealso \code{\link{zlm}} and \code{\link{as.zlm}} for creating \code{zlm}
#' objects,\cr \code{\link{density.zlm}}, \code{\link{predict.zlm}} and
#' \code{\link{summary.zlm}} for other posterior results
#' @references \url{http://bms.zeugner.eu}
#' @keywords classes
setClass("zlm",representation(coefficients="numeric",residuals="numeric",rank="numeric",fitted.values="numeric",df.residual="numeric",call="call",terms="formula",model="data.frame",coef2moments="numeric",marg.lik="numeric",gprior.info="list"))


#' Class "gprior"
#' 
#' An object pertaining to a coefficient prior
#' 
#' 
#' @name gprior-class
#' @docType class
#' @section Objects from the Class: A \code{gprior} object holds descriptions
#' and subfunctions pertaining to coefficient priors. Functions such as
#' \code{\link{bms}} or \code{\link{zlm}} rely on this class to 'convert' the
#' output of OLS results into posterior expressions for a Bayesian Linear
#' Model. Post-processing functions such as \code{\link{density.bma}} also
#' resort to gprior objects.\cr There are currently three coefficient prior
#' structures built into the BMS package, generated by the following functions
#' (cf. Feldkircher and Zeugner, 2009) : \cr \code{gprior.constg.init}: creates
#' a Zellner's g-prior object with constant \code{g}.\cr
#' \code{gprior.eblocal.init}: creates an Empricial Bayes Zellner's g-prior.\cr
#' \code{gprior.hyperg.init}: creates a hyper g-prior with a Beta-prior on the
#' shrinkage parameter.\cr The following describes the necessary slots
#' @author Martin Feldkircher and Stefan Zeugner
#' @seealso \code{\link{bms}} and \code{\link{zlm}} for creating \code{bma} or
#' \code{zlm} objects. \cr Check the appendix of \code{vignette(BMS)} for a
#' more detailed description of built-in priors.\cr Check
#' \url{http://bms.zeugner.eu/custompriors.php} for examples.
#' @references Feldkircher, M. and S. Zeugner (2009): Benchmark Priors
#' Revisited: On Adaptive Shrinkage and the Supermodel Effect in Bayesian Model
#' Averaging, IMF Working Paper 09/202.
#' @keywords classes
#' @examples
#' 
#' 
#' data(datafls)
#' mm1=bms(datafls[,1:10], g="EBL")
#' gg=mm1$gprior.info # is the g-prior object, augmented with some posterior statistics
#' 
#' mm2=bms(datafls[,1:10], g=gg) #produces the same result
#' 
#' mm3=bms(datafls[,1:10], g=BMS:::.gprior.eblocal.init) 
#' #this passes BMS's internal Empirical Bayes g-prior object as the coefficient prior 
#' # - any other obejct might be used as well
#' 
#' 
#' 
setClass("gprior")


#' Class "mprior"
#' 
#' An object pertaining to a BMA model prior
#' 
#' 
#' @name mprior-class
#' @docType class
#' @section Objects from the Class: An \code{mprior} object holds descriptions
#' and subfunctions pertaining to model priors. The BMA functions
#' \code{\link{bms}} and post-processing functions rely on this class. \cr
#' There are currently five model prior structures built into the BMS package,
#' generated by the following functions (cf. the appendix of
#' \code{vignette(BMS)}): \cr \code{mprior.uniform.init}: creates a uniform
#' model prior object.\cr \code{mprior.fixedt.init}: creates the popular
#' binomial model prior object with common inclusion probabilities.\cr
#' \code{mprior.randomt.init}: creates a beta-binomial model prior object.\cr
#' \code{mprior.pip.init}: creates a binomial model prior object that allows
#' for defining individual prior inclusion probabilities.\cr
#' \code{mprior.customk.init}: creates a model prior object that allows for
#' defining a custom prior for each model parameter size.\cr The following
#' describes the necessary slots:
#' @author Martin Feldkircher and Stefan Zeugner
#' @seealso \code{\link{bms}} for creating \code{bma} objects. \cr Check the
#' appendix of \code{vignette(BMS)} for a more detailed description of built-in
#' priors.\cr Check \url{http://bms.zeugner.eu/custompriors.php} for examples.
#' @keywords classes
setClass('mprior')






#' @importFrom stats density


.flsresultlist =function(item=NULL) {
  if (is.null(item)) return(.flsresults)
  return(.flsresults[[item]])
}