###########################################
# This version: adjusted on 2022-07-09    #
###########################################
# it includes all the auxiliary functions that should only be called inside as well as OUTSIDE the bms function


 #######################
 # FUNCTIONS FOR USERS #  
#########################################################################
 # SEE ALSO PLOTTING FUCNTIONS AT THE END 
 # beta.draws.bma, pmp.bma, estimates.bma, info.bma
 # cf. as well the functions combine_chains & [.bma farther below
 # moreover: bin2hex, hex2bin, f21simple, fullmodel.ssq, is.bma, print.bma, summary.bma, print.topmod, is.topmod
 # plotting: plotModelsize, plotComp, plotConv, plotDensity, plot.bma, image.bma




 



#' Coefficients of the Best Models
#' 
#' Returns a matrix whose columns are the (expected value or standard
#' deviations of) coefficients for the best models in a bma object.
#' 
#' 
#' @param bmao a 'bma' object (as e.g. resulting from \code{\link{bms}})
#' @param stdev if \code{stdev=FALSE} then \code{beta.draws.bma} returns the
#' (conditional) posterior expected values of the coefficients (i.e. 'Bayesian
#' coefficients'). If \code{stdev=TRUE} it returns their posterior standard
#' deviations.
#' @return Each column presents the coefficients for the model indicated by its
#' column name. The zero coefficients are the excluded covariates per model.
#' Note that the coefficients returned are only those of the best (100) models
#' encountered by the \code{bma} object (cf. argument \code{nmodels} of
#' \code{\link{bms}}).
#' 
#' For aggregate coefficients please refer to \code{\link{coef.bma}}.
#' @note Note that the elements of \code{beta.draws.bma(bmao)} correspond to
#' \code{bmao$topmod$betas()}
#' @seealso \code{\link{bms}} for creating bms objects, \code{\link{coef.bma}}
#' for aggregate coefficients
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#' 
#'   #sample a bma object:
#'   data(datafls)
#'   mm=bms(datafls,burn=500,iter=5000,nmodel=20)
#'   
#'   #coefficients for all
#'   beta.draws.bma(mm) 
#'   
#'   #standard deviations for the fourth- to eight best models
#'   beta.draws.bma(mm[4:8],TRUE); 
#' 
#' @export
beta.draws.bma <- function(bmao,stdev=FALSE) {
      # constructs a nice matrix of the betas of the best models stored in topmods
      # bmao: bma object
  if (!is.bma(bmao)) {stop("you need to provide a BMA object"); return()}

  resmat=.post.beta.draws(bmao$topmod, bmao$reg.names,FALSE)
  if (stdev) {
   mom2=.post.beta.draws(bmao$topmod, bmao$reg.names,TRUE)
   resmat=sqrt(mom2-resmat^2)
  }
  return(resmat)
  
}

#' Posterior Model Probabilities
#' 
#' Returns the posterior model probabilites for the best models encountered by
#' a 'bma' object
#' 
#' A call to bms with an MCMC sampler (e.g.
#' \code{bms(datafls,mcmc="bd",nmodel=100)} uses a Metropolis-Hastings
#' algorithm to sample through the model space - and the frequency of how often
#' models are drawn converges to the distribution of their posterior marginal
#' likelihoods.  While sampling, each 'bma' object stores the best models
#' encountered by its sampling chain with their marginal likelihood and their
#' MCMC frequencies.
#' 
#' \code{pmp.bma} then allows for comparing the posterior model probabilities
#' (PMPs) for the two different methods, similar to \code{\link{plotConv}}.  It
#' calculates the PMPs based on marginal likelihoods (first column) and the
#' PMPs based on MCMC frequencies (second column) for the best x models stored
#' in the bma object.
#' 
#' The correlation of the two columns is an indicator of how well the MCMC
#' sampler has converged to the actual PMP distribution - it is therefore also
#' given in the output of \code{\link{summary.bma}}.
#' 
#' The second column is slightly different in case the \code{\link{bms}}
#' argument \code{mcmc} was set to \code{mcmc="enumeration"}: In this case, the
#' second column is also based on marginal likelihoods. The correlation between
#' the two columns is therefore one.
#' 
#' @param bmao A bma object (see argument \code{nmodel} in \code{\link{bms}}),
#' alternatively an object of class \code{\link{topmod}}
#' @param oldstyle For normal use, leave this at \code{FALSE}. It is an
#' argument for compatibility with older BMS versions - see section 'Notes'
#' @return the result is a matrix, its row names describe the model binaries\cr
#' There are two columns in the matrix: \item{PMP (Exact)}{posterior model
#' probabilities based on the posterior likelihoods of the best models in
#' \code{bmao} } \item{PMP (MCMC)}{posterior model probabilities of the best
#' models in \code{bmao} based on their MCMC frequencies, relative to all
#' models encountered by \code{bmao} - see 'Details' }
#' @note The second column thus shows the PMPs of the best models relative to
#' all models the call to \code{\link{bms}} has sampled through (therefore
#' typically the second column adds up to less than one).  The first column
#' relates to the likelihoods of the best models, therefore it would add up to
#' 1.  In order estimate for their marginal likelihoods with respect to the
#' other models (the ones not retained in the best models), these PMP aadding
#' up to one are multiplied with the sum of PMP of the best models accroding to
#' MCMC frequencies.  Therefore, the two columns have the same column sum.
#' 
#' CAUTION: In package versions up to \code{BMS 0.2.5}, the first column was
#' indeed set always equal to one. This behaviour can still be mimicked by
#' setting \code{oldstyle=TRUE}.
#' 
#' @seealso \code{\link{plotConv}} for plotting \code{pmp.bma},
#' \code{\link{pmpmodel}} to obtain the PMP for any individual model,
#' \code{\link{bms}} for sampling bma objects
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#' 
#' ## sample BMA for growth dataset, MCMC sampler
#' data(datafls)
#' mm=bms(datafls[,1:10],nmodel=20, mcmc="bd")
#' 
#' ## mmodel likelihoods and MCMC frequencies of best 20 models
#' print(mm$topmod)
#' 
#' pmp.bma(mm)
#' #first column: posterior model prob based on model likelihoods,
#' #  relative to best models in 'mm'
#' #second column: posterior model prob based MCMC frequencies,
#' #  relative to all models encountered by 'mm'
#' 
#' #consequently, first column adds up to one
#' #second column shows how much of the sampled model space is
#' # contained in the best models
#' colSums(pmp.bma(mm))
#' 
#' 
#' #correlation betwwen the two shows how well the sampler converged
#' cor(pmp.bma(mm)[,1],pmp.bma(mm)[,2])
#' 
#' #is the same as given in summary.bma
#' summary(mm)["Corr PMP"]
#' 
#' #plot the two model probabilites
#' plotConv(mm)
#' 
#' #equivalent to the following chart
#' plot(pmp.bma(mm)[,2], type="s")
#' lines(pmp.bma(mm)[,1],col=2)
#' 
#' 
#' #moreover, note how the first column is constructed
#' liks=exp(mm$top$lik())
#' liks/sum(liks)
#' pmp.bma(mm)[,1] #these two are equivalent
#' 
#' 
#' 
#' #the example above does not converge well,
#' #too few iterations and best models
#' # this is already better, but also not good
#' mm=bms(datafls[,1:10],burn=2000,iter=5000,nmodel=200)
#' 
#' 
#' # in case the sampler has been 'enumeration' instead of MCMC,
#' # then both matrix columns are of course equivalent
#' mm=bms(datafls[,1:10],nmodel=512,mcmc="enumeration")
#' cor(pmp.bma(mm)[,1],pmp.bma(mm)[,2])
#' colSums(pmp.bma(mm))
#' 
#' 
#' @export
pmp.bma <- function(bmao, oldstyle=FALSE) {
     # constructs nice matrix with PMP analytical and PMP MCMC for best models stored in topmods
     # Prob of top "nmodel" models, analytical (after deleting rest of models)'
     # bmao: either the topmodel object or a "bma" object

     if (!(is.bma(bmao)||is.topmod(bmao))) stop("bmao needs to be a 'bma' object!")
     if (is.topmod(bmao)) {
       topmods=bmao
       was.enum=FALSE
       cumsumweights=sum(topmods$ncount())
       log.null.lik=0
     } else {
       topmods=bmao$topmod 
       log.mid.lik=bmao$info$log.mid.lik; if (is.null(log.mid.lik)) {log.mid.lik=0}
       log.null.lik=(1-bmao$info$N)/2*log(as.vector(crossprod(bmao$arguments$X.data[,1]-mean(bmao$arguments$X.data[,1])))) - log.mid.lik
       cumsumweights=bmao$info$cumsumweights
       was.enum=(bmao$arguments$mcmc=="enum")       
     }
    lt1=suppressWarnings(topmods$lik() - max(topmods$lik()))    # do this to get only positive probabilities
    lt1=exp(lt1)/sum(exp(lt1))
     
    
    if (was.enum) {
       lt2=exp(topmods$lik())/cumsumweights #Prob of top "nmodel" models, (loglik based)
    } else {
       lt2=topmods$ncount()/cumsumweights #MCMC Prob of top "nmodel" models, numerical 
    }
    
     cpoint=min(length(lt1),length(lt2))
     lt1=lt1[1:cpoint]; lt2=lt2[1:cpoint]
     if (!oldstyle) lt1 <- lt1*sum(lt2)
    #rbind the probs to the topmodmatrix
    topmodout=rbind(lt1, lt2)
    
    
    
    rownames(topmodout)=c("PMP (Exact)","PMP (MCMC)")
    colnames(topmodout)=topmods$bool()
    

     return(t(topmodout))
}





#' Posterior Model Probability for any Model
#' 
#' Returns the posterior model probability for any model based on bma results
#' 
#' If the model as provided in \code{model} is the null or the full model, or
#' is contained in \code{bmao}'s topmod object (cf. argument \code{nmodel} in
#' \code{\link{bms}}), \cr then the result is the same as in
#' \code{\link{pmp.bma}}.\cr If not and \code{exact=TRUE}, then \code{pmpmodel}
#' estimates the model based on comparing its marginal likelihood (times model
#' prior) to the likelihoods in the \code{topmod} object and multiplying by
#' their sum of PMP according to MCMC frequencies,
#' 
#' @param bmao A bma object as created by \code{\link{bms}}.
#' @param model A model index - either variable names, or a logical with model
#' binaries, or the model hexcode (cf. \code{\link{hex2bin}}, or a numeric with
#' positions of the variables to be included.
#' @param exact If \code{TRUE}, then the resulting PMP is based on analytical
#' model likelihoods (works for any model). \cr If \code{FALSE}, the the
#' resulting PMP is derived from MCMC frequencies (works only for the null and
#' fullmodel, as well as for models contained in \code{bmao}'s topmod
#' object.\cr If \code{bmao} is based on enumeration (cf. argument \code{mcmc}
#' in \code{\link{bms}}, then \code{exact} does not matter.
#' @return A scalar with (an estimate of) the posterior model probability for
#' \code{model}
#' @seealso \code{\link{pmp.bma}} for similar
#' functions
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#' 
#' ## sample BMA for growth dataset, enumeration sampler
#' data(datafls)
#' mm=bms(datafls[,1:10],nmodel=5)
#' 
#' #show the best 5 models:
#' pmp.bma(mm)
#' #first column: posterior model prob based on model likelihoods,
#' #second column: posterior model prob based MCMC frequencies,
#' 
#' ### Different ways to get the same result: #########
#' 
#' #PMP of 2nd-best model (hex-code representation)
#' pmpmodel(mm,"00c")
#' 
#' #PMP of 2nd-best model (binary representation)
#' incls=as.logical(beta.draws.bma(mm)[,2])
#' pmpmodel(mm,incls)
#' 
#' #PMP of 2nd-best model (via variable names)
#' #names of regressors in model "00c": 
#' names(datafls[,2:10])[incls]
#' pmpmodel(mm,c("SubSahara", "LatAmerica"))
#' 
#' #PMP of 2nd-best model (via positions)
#' pmpmodel(mm,c(6,7))
#' 
#' ####PMP of another model #########
#' pmpmodel(mm,1:5)
#' 
#' 
#' @export
pmpmodel= function(bmao, model=numeric(0), exact=TRUE) {
  #returns the PMP of any model 
  #if model is null, full or contained in topmod, returns the PMP as in pmp.bma 
  #if model is out of topmod but bmao is enumeration, calculates the marg-lik and compares
  #if model is out of topmod and bmao is MCMC-sampled, then estimates PMP by comparing marg-lik and ncounts
  
  if (!is.bma(bmao)) stop("bmao needs to be a bma object")
  if (!is.vector(model)) stop("model needs to be vector denoting a single model.")
  K=bmao$info$K
  was.enum=(bmao$arguments$mcmc=="enum")
  emptyindex=logical(K)
  modelhex=""
  
  #Conversion of user input
  if (length(model)==0L) model=numeric(0)
  if ((is.character(model))&&(all(model %in% bmao$reg.names))) {
      mix=match(model, bmao$reg.names)
      if (any(is.na(mix))) stop("Provided variable names do not conform to bma object")
      emptyindex[mix]=TRUE; model=emptyindex
  } else if ((length(model)==1L)&&all(is.character(model))) {
      modelhex=model
      model=as.logical(hex2bin(model)); if (length(model)>K) model=model[-(1:(length(model)-K))]
  } else if (is.logical(model)||((length(model)==K)&&(is.numeric(model)&& max(model)<2))) {
    if (length(model)>K) model=model[-(1:(length(model)-K))]
      model=as.logical(model)
  } else if (is.numeric(model)) {
      emptyindex[model]=TRUE; model=emptyindex
  } else stop("model needs to be an integer, logical or character model index representation (hexcode or variable names)")
  if (any(is.na(model))) stop("Provided model index does not seem to exist.")
  
  if (modelhex=="") modelhex=bin2hex(model)
  # now model is a logical vector, modelhex its hex representation
  
  #last user check
  fixed.pos = bmao$mprior.info$fixed.pos
  if (is.null(fixed.pos)) fixed.pos=numeric(0)
  if (any(model[fixed.pos]!=TRUE)) stop("Such a model was excluded by bmao's argument fixed.reg")
  
  # prior calculations
  bools=bmao$topmod$bool()
  liks=bmao$topmod$lik()
  ncounts=bmao$topmod$ncount()
  cumsumweights=bmao$info$cumsumweights
  yty=as.vector(crossprod(bmao$arguments$X.data[,1,drop=TRUE]-mean(bmao$arguments$X.data[,1,drop=TRUE])))
  log.mid.lik=bmao$info$log.mid.lik; if (is.null(log.mid.lik)) {log.mid.lik=0}
  log.null.lik= bmao$gprior.info$lprobcalc$just.loglik(ymy=yty,k=0) - log.mid.lik
  
  
  
  #look up whether model is in topmodels
  mix=match(modelhex,bools) 
  
  # first treat case when MCMC sampler and exact =FALSE
  if ((!exact) && (!was.enum)) {
    if (!is.na(mix)) {
      return(ncounts[[mix]]/cumsumweights)
    } else if (!any(model)||all(model)) {
      return(bmao$info$k.vec[sum(model)+1]/cumsumweights)   
    } else {
      stop("Model MCMC-based PMP cannot be found. Try exact=TRUE .")
    }
  } 
  
  # now treat the cases where exact = TRUE or was.enum=TRUE 
  if (!is.na(mix)) {
    loglik=liks[mix]  #if yes and exact return its PMP
  } else if (was.enum && (!any(model)||all(model))) {
    loglik = log(bmao$info$k.vec[sum(model)+1])+log.null.lik # stuff in info is saved as sum(Lik_model / Lik_null)
  } else {
    if (!was.enum && (length(liks)==0L)) stop("bmao needs to contain more than 0 top models to provide an estimate for your model index.")
    if (sum(model)==0L) {
      loglik= log.null.lik + bmao$mprior.info$pmp(ki=0, mdraw=rep(0,K), ymy=yty)
    } else {
      zz=zlm(bmao$arguments$X.data[,c(TRUE,model),drop=FALSE], g=bmao$gprior.info)
      loglik=zz$marg.lik+bmao$mprior.info$pmp(ki=sum(model), mdraw=as.numeric(model), ymy=zz$olsres$ymy)
    }
  }
  
  if (was.enum) {
    return(exp(loglik-log.null.lik)/cumsumweights)
  }
  
  #estimate pmp of model compared to top models
  # then multply with a factor by comparing ncount and liks so to estimate overall pmp of model
    pmp_withintop=exp(loglik-log.null.lik)/sum(exp(liks-log.null.lik))
    return(pmp_withintop*sum(ncounts)/cumsumweights) 
  
}


#' @rdname coef.bma
#' @export
estimates.bma <- function(bmao,exact=FALSE,order.by.pip=TRUE,include.constant=FALSE,incl.possign=TRUE,std.coefs=FALSE,condi.coef=FALSE) {
   # constructs a nice estimates matrix with 5 columns: 1) PIP, 2) E(beta|Y), 3) Var(beta|Y), 4) pos. coef. sign (cond. on inclusion, optional), 5) Index in X.data
   # bmao: bma object; 
   # exact: if True, then calcs posterior stats based on analytical PMPs (likelihoods) in topmods object; if False, then does it as a weighted average (typically with MCMC frequencies)
   # order.by.pip: if true then the matrix is sorted according to the PIPs (first column)
   # include.constant: if true, then add the constant in the last row of resultant matrix
   # incl.possign: if true, then the fouth column details how often the coefficent sign was positive (conditional on inclusion)
   # std.coefs: if true, then coefficents are standardized
   # condi.coefs: if true, then the coefficents and standard deviations are not given as overall expected values, but as expected values given inclusion
   
       
  if (!is.bma(bmao)) {stop("you need to provide a BMA object"); return()}
  if (exact) {
   #if (bmao$arguments$beta.save==FALSE) stop("exact=TRUE needs betas from the draws: Run estimation again and save betas via setting beta.save=TRUE")
   if (bmao$topmod$nbmodels==0) stop("exact=TRUE needs at least one 'top model': Run estimation again and set nmodel>0")
  }
  bmaest=.post.estimates(bmao$info$b1mo,bmao$info$b2mo,bmao$info$cumsumweights,bmao$info$inccount,bmao$topmod,bmao$arguments$X.data,bmao$reg.names,bmao$info$pos.sign,exact,order.by.pip,include.constant,incl.possign,std.coefs,condi.coef)
  return(bmaest)
}

#' @rdname summary.bma
#' @export
info.bma <- function(bmao) {
    # constructs an 'info' character matrix with 1 row and the following columns:
    # bmao: bma object
    # output: a list that contains information about the "Mean nr. of Regressors" (not counting the constant term),
    #             "Draws"=posterior draws, "Burnins"=nr. of burnins taken, "Time" denotes total time elapsed
    #             since calling the "bms" function, "Nr. of models visited" counts each time a model is accepted.
    #             Note that we do not take into account the case of revisiting models by the sampler.
    #             Modelspace is simply indicating the whole model space (2^K) and percentage visited is 
    #             the nr. of models visited as a percentage of 2^K. "Corr PMP" is the correlation between
    #             analytic and MCMC posterior model probabilites, where a correlation of 0.99 indicates 
    #             excellent convergence. For "nmodel=100" the best 100 models are considered for the correlation
    #             analysis. Finally, Nr. of Observations is given in the output as well.

    
    foo=bmao$info
    iter=foo$iter; burn=foo$burn; timed=foo$timed; models.visited=foo$models.visited; corr.pmp=foo$corr.pmp; K=foo$K; N=foo$N; msize=foo$msize; cumsumweights=foo$cumsumweights
    if (is.element("mprior.info",names(bmao))) {
      prior= paste(bmao$mprior.info$mp.mode, "/", bmao$mprior.info$mp.msize)
    } else {
      if (is.element("theta",names(bmao$arguments))&&is.element("prior.msize",names(bmao$arguments))) {
         if (!is.null(bmao$arguments$theta)&!is.null(bmao$arguments$prior.msize)) prior=paste(bmao$arguments$theta, "/", bmao$arguments$prior.msize) else prior=NA
      } else {
        prior=paste(bmao$arguments$mprior, "/", bmao$arguments$mprior.size)
      }
    }
    gprior.info=bmao$gprior.info
    gprior.choice=gprior.info$gtype
    
    model.space=2^K
    fraction.model=models.visited/model.space*100
    fraction.topmodel=sum(bmao$topmod$ncount())/iter*100 # this is the number of visits the topmodels accounted for relative to



    if (gprior.info$gtype=="hyper") {gprior.choice=paste(gprior.choice," (a=",2+signif(gprior.info$hyper.parameter-2,digits=4),")",sep="")}
    nr.reg=msize/cumsumweights

   info<-as.character(c(format(round(nr.reg,4),nsmall=4),format(iter,nsmall=0),format(burn,nsmall=0),
                      format(timed,nsmall=4),models.visited,format(model.space,digits=2),
                      format(fraction.model,digits=2),format(fraction.topmodel,digits=2),format(round(.cor.topmod(bmao$topmod),4),nsmall=4),
                      format(N,nsmall=4),prior, gprior.choice))



    names(info)<-c("Mean no. regressors", "Draws","Burnins", "Time", "No. models visited",
                 "Modelspace 2^K", "% visited","% Topmodels","Corr PMP","No. Obs.", "Model Prior", "g-Prior")

    if (gprior.info$return.g.stats) {
      gpriorav=gprior.info$shrinkage.moments[1]
      gstatsprint= paste("Av=", format(gpriorav,digits=4),sep="")
      if (length(gprior.info$shrinkage.moments)>1) { 
        gpriorstdev=sqrt(gprior.info$shrinkage.moments[2]-gprior.info$shrinkage.moments[1]^2)
        gstatsprint = paste(gstatsprint, ", Stdev=", format(gpriorstdev,digits=2),sep="") 
      }
      info <- c(info, gstatsprint)
      names(info)[13]<- "Shrinkage-Stats"
    }
  return(info)
}

#' Predict Method for bma Objects
#' 
#' Expected value of prediction based on 'bma' object
#' 
#' 
#' @param object a bma object - see \code{\link{bms}}
#' @param newdata An optional data.frame, matrix or vector containing variables
#' with which to predict. If omitted, then (the expected values of) the fitted
#' values are returned.
#' @param exact If \code{FALSE} (default), then prediction is based on all
#' models (i.e. on their MCMC frequencies in case the \code{\link{bms}}
#' parameter \code{mcmc} was set to an mcmc sampler.\cr If \code{TRUE}, then
#' prediction is based on analytical likelihoods of the best models retained in
#' \code{object} - cf. \code{\link{bms}} parameter \code{nmodel}.
#' @param topmodels index of the models with whom to predict: for instance,
#' \code{topmodels=1} predicts based solely on the best model, whereas
#' \code{topmodels=1:5} predicts based on a combination of the five best
#' models.\cr Note that setting \code{topmodels} triggers \code{exact=TRUE}.
#' @param \dots further arguments passed to or from other methods.
#' @return A vector with (expected values of) fitted values.
#' @seealso \code{\link{coef.bma}} for obtaining coefficients,
#' \code{\link{bms}} for creating bma objects, \code{\link{predict.lm}} for a
#' comparable function
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'  data(datafls)
#'  mm=bms(datafls,user.int=FALSE)
#'  
#'  predict(mm) #fitted values based on MCM frequencies
#'  predict(mm, exact=TRUE) #fitted values based on best models
#'  
#'  predict(mm, newdata=1:41) #prediction based on MCMC frequencies 
#'  
#'  predict(mm, newdata=datafls[1,], exact=TRUE) #prediction based on a data.frame
#'  
#'  # the following two are equivalent:
#'  predict(mm, topmodels=1:10)
#'  predict(mm[1:10], exact=TRUE)
#'  
#'  
#' @export
predict.bma <- function(object, newdata=NULL, exact=FALSE, topmodels=NULL, ...) {
   # does basic fitting in expected values
   # object: a bma object
   # newdata: newdata to be supplied (just eas in predict.lm)
   # topmodels: The index of the bestmodels to be included (e.g. 1 for best, 1:3 for three best); setting this parameter triggers exact=TRUE
   # exact: TRUE if fit should be based on exact likelihood of best models, FALSE if should be based on MCMC freqs
   # output: a vector with fitted values
   
    if (!is.bma(object)) {stop("you need to provide a BMA object"); return()}

    # check the topmodels argument
    if (!is.null(topmodels)) {
      if (!(is.numeric(topmodels)&&is.vector(topmodels))) {
        stop("topmodels must denote the models to take into account, e.g. 1:5 for the best five.")
      } else if (object$topmod$nbmodels < max(topmodels)) {
        stop(paste("Only",object$topmod$nbmodels,"best models are available, but you asked to take the", max(topmodels), "-best model into account."))
      }
      object=object[unique(topmodels)]
    }
    if ((!missing(topmodels))&&missing(exact)) exact=TRUE

    #get the betas as required
    betas=estimates.bma(object,exact=exact,order.by.pip=FALSE,include.constant=FALSE,std.coefs=FALSE,condi.coef=FALSE)[,2]

    #check the newdata argument
    if (is.null(newdata)) {
       newX<-as.matrix(object$arguments$X.data[,-1, drop=FALSE])
    } else {
       newX=as.matrix(newdata)
       if (!is.numeric(newX)) stop("newdata must be numeric!")
       if (is.vector(newdata)) newX=matrix(newdata,1)
       if (ncol(newX)!=length(betas)) {
         if (ncol(newX)==length(betas)+1) {
             newX=newX[,-1,drop=FALSE] # this is to achieve a bevavior similar to predict.lm in this case
         } else {
           stop("newdata must be a matrix or data.frame with", length(betas), "columns.")
         }
       }
       orinames=colnames(object$arguments$X.data[,-1, drop=FALSE])
       if (!is.null(colnames(newX))&& !is.null(orinames)) { #this is a user check whether columns had been submitted in the wrong  order
        if (all(orinames %in% colnames(newX) ) && !all(orinames == colnames(newX))  ) {
            warning("argument newdata had to be reordered according to its column names. Consider submitting the columns of newdata in the right order.")
            newX=newX[,orinames, drop=FALSE] 
        }
       }
       
    }
    cons=.post.constant(object$arguments$X.data,betas)
    return(as.vector(newX%*%betas)+cons)
}


#' OLS Statistics for the Full Model Including All Potential Covariates
#' 
#' A utility function for reference: Returns a list with R2 and sum of squares
#' for the OLS model encompassing all potential covariates that are included in
#' a bma object.
#' 
#' 
#' @param yX.data a bma object (cf. \code{\link{bms}}) - alternatively a
#' \link{data.frame} or \link{matrix} whose first column is the dependent
#' variable
#' @return Returns a list with some basic OLS statistics \item{R2}{The
#' R-squared of the full model} \item{ymy}{The sum of squares of residuals of
#' the full model} \item{ypy}{The explained sum of squares of the full model}
#' \item{yty}{The sum of squares of the (demeaned) dependent variable}
#' \item{Fstat}{The F-statistic of the full model}
#' @note This function is just for quick comparison; for proper OLS estimation
#' consider \code{\link{lm}}
#' @seealso \code{\link{bms}} for creating bma objects, \code{\link{lm}} for
#' OLS estimation
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#' data(datafls)
#' mm=bms(datafls)
#' 
#' fullmodel.ssq(mm)
#' 
#' #equivalent:
#' fullmodel.ssq(datafls)
#' 
#' 
#' @export
fullmodel.ssq <- function(yX.data) {
  # yX.data: a dataframe
  # returns the OLS sums of sqares for yX.data, where the first column is the dependent:
  # R2: r-squared; ymy: resid SS; ypy: explained SS; yty: (y-ymean)'(y-ymean)

  if (is.bma(yX.data)) {yX.data <- yX.data$arguments$X.data}
  y<-as.matrix(yX.data[,1])
  X<-as.matrix(yX.data[,2:ncol(yX.data)])
  N<-nrow(X)
  K=ncol(X)

  y.mean=mean(y)
  y<-y-matrix(y.mean,N,1,byrow=TRUE)
  X.mean=colMeans(X)
  X<-X-matrix(X.mean,N,K,byrow=TRUE)
  
  Xqr<-qr(X)
  yty=as.numeric(crossprod(y))
  ymy=as.numeric(crossprod(qr.resid(Xqr,y)))
  ypy=as.numeric(crossprod( qr.fitted(Xqr,y)))
  R2=ypy/yty
  return(list(R2=R2,ymy=ymy,ypy=ypy,yty=yty,Fstat=(R2/(1-R2))*(N-K-1)/K))
  
}

#' @export
print.bma <- function(x,...) {
  #defines how to print a bmao object (e.g. to the console)
  if (!is.bma(x)) {return(print(x))} 
  print(estimates.bma(x),include.constant=TRUE,...)
  cat("\n")
  print(info.bma(x),...)
  cat("\n")
}

#' Summary Statistics for a 'bma' Object
#' 
#' Returns a vector with summary statistics for a 'bma' object
#' 
#' \code{info.bma} is equivalent to \code{summary.bma}, its argument
#' \code{bmao} conforms to the argument \code{object}
#' 
#' @aliases summary.bma info.bma
#' @param object a list/object of class 'bma' that typically results from the
#' function \code{bms} (see \code{\link{bms}} for details)
#' @param \dots further arguments passed to or from other methods
#' @param bmao same as \code{object}
#' @return A character vector summarizing the results of a call to \code{bms}
#' \item{Mean no. of Regressors}{ the posterior mean of model size}
#' \item{Draws}{the number of iterations (ex burn-ins)} \item{Burnins}{the
#' number of burn-in iterations} \item{Time}{the time spent on iterating
#' through the model space} \item{No. of models visited}{the number of times a
#' model was accepted (including burn-ins)} \item{Modelspace }{the total model
#' space \eqn{2^K}}\item{list(list("2^K"))}{the total model space \eqn{2^K}}
#' \item{Percentage visited}{\code{No. of models visited/Modelspace*100}}
#' \item{Percentage Topmodels}{number of times the best models were drawn in
#' percent of \code{Draws}} \item{Corr. PMP}{the correlation between the MCMC
#' frequencies of the best models (the number of times they were drawn) and
#' their marginal likelihoods.} \item{No. Obs.}{Number of observations}
#' \item{Model Prior}{a character conforming to the argument \code{mprior} of
#' \code{bms}, and the expected prior model size} \item{g-prior}{a character
#' corresponding to argument \code{g} of function \code{bms}}
#' \item{Shrinkage-Stats}{Posterior expected value und standard deviation (if
#' applicable) of the shrinkage factor. Only included if argument
#' \code{g.stats} of function \code{bms} was set to TRUE}
#' @note All of the above statistics can also be directly extracted from the
#' bma object (\code{bmao}). Therefore \code{summary.bma} only returns a
#' character vector.
#' @seealso \code{\link{bms}} and \code{\link{c.bma}} for functions creating
#' bma objects, \code{print.bma} makes use of \code{summary.bma}.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'   data(datafls)
#' 
#'   m_fixed=bms(datafls,burn=1000,iter=3000,user.int=FALSE, )
#'   summary(m_fixed)
#'    
#'   m_ebl=bms(datafls,burn=1000,iter=3000,user.int=FALSE, g="EBL",g.stats=TRUE)
#'   info.bma(m_ebl)
#' 
#' @export
summary.bma <-function(object,...) {
  #just an alias for info.bma
  info.bma(object)
}

#' Posterior Inclusion Probabilities and Coefficients from a 'bma' Object
#' 
#' Returns a matrix with aggregate covariate-specific Bayesian model Averaging:
#' posterior inclusion probabilites (PIP), post. expected values and standard
#' deviations of coefficients, as well as sign probabilites
#' 
#' More on the argument \code{exact}: \cr In case the argument
#' \code{exact=TRUE}, the PIPs, coefficient statistics and conditional sign
#' probabilities are computed on the basis of the (500) best models the
#' sampling chain encountered (cf. argument \code{nmodel} in
#' \code{\link{bms}}). Here, the weights for Bayesian model averaging (BMA) are
#' the posterior marginal likelihoods of these best models. \cr In case
#' \code{exact=FALSE}, then these statistics are based on all accepted models
#' (except burn-ins): If \code{mcmc="enumerate"} then this are simply all
#' models of the traversed model space, with their marginal likelihoods
#' providing the weights for BMA.\cr If, however, the bma object \code{bmao}
#' was based on an MCMC sampler (e.g. when \code{\link{bms}} argument
#' \code{mcmc="bd"}), then BMA statistics are computed differently: In contrast
#' to above, the weights for BMA are MCMC frequencies, i.e. how often the
#' respective models were encountered by the MCMC sampler. (cf. a comparison of
#' MCMC frequencies and marginal likelihoods for the best models via the
#' function \code{\link{pmp.bma}}).
#' 
#' @aliases estimates.bma coef.bma
#' @param object,bmao a 'bma' object (cf. \code{\link{bms}})
#' @param exact if \code{exact=FALSE}, then PIPs, coefficients, etc. will be
#' based on aggregate information from the sampling chain with posterior model
#' distributions based on MCMC frequencies (except in case of enumeration - cf.
#' 'Details');\cr if \code{exact=TRUE}, estimates will be based on the
#' \code{\link[=bms]{nmodel}} best models encountered by the sampling chain,
#' with the posterior model distribution based on their \emph{exact} marginal
#' likelihoods - cf. 'Details' below.
#' @param order.by.pip \code{order.by.pip=TRUE} orders the resulting matrix
#' according to posterior inclusion probabilites, \code{order.by.pip=FALSE}
#' ranks them according to the original data (order of the covariates as in
#' provided in \code{X.data} to \code{\link{bms}}), default \code{TRUE}
#' @param include.constant If \code{include.constant=TRUE} then the resulting
#' matrix includes the expected value of the constant in its last row. Default
#' \code{FALSE}
#' @param incl.possign If \code{incl.possign=FALSE}, then the sign probabilites
#' column (cf. 'Values' below) is omitted from the result. Default \code{TRUE}
#' @param std.coefs If \code{std.coefs=TRUE} then the expected values and
#' standard deviations are returned in standardized form, i.e. as if the
#' original data all had mean zero and variance 1. If \code{std.coefs=FALSE}
#' (default) then both expected values and standard deviations are returned 'as
#' is'.
#' @param condi.coef If \code{condi.coef=FALSE} (default) then coefficients
#' \eqn{\beta_i} and standard deviations are unconditional posterior expected
#' values, as in standard model averaging; if \code{condi.coef=FALSE} then they
#' are given as conditional on inclusion (equivalent to \eqn{\beta_i / PIP_i}).
#' @param ...  further arguments for other \code{\link{coef}} methods
#' @return A matrix with five columns (or four if \code{incl.possign=FALSE})
#' \item{Column 'PIP'}{Posterior inclusion probabilities \eqn{\sum p(\gamma|i
#' \in \gamma, Y) / sum p(\gamma|Y) }} \item{Column 'Post Mean'}{posterior
#' expected value of coefficients, unconditional \eqn{E(\beta|Y)=\sum
#' p(\gamma|Y) E(\beta|\gamma,Y)}, where \eqn{E(\beta_i|\gamma,i \notin \gamma,
#' Y)=0} if \code{condi.coef=FALSE}, or conditional on inclusion
#' (\eqn{E(\beta|Y) / \sum p(\gamma|Y, i \in \gamma) } ) if
#' \code{condi.coef=TRUE}} \item{Column 'Post SD'}{posterior standard deviation
#' of coefficients, unconditional or conditional on inclusion, depending on
#' \code{condi.coef}} \item{Column 'Cond.Pos.Sign'}{The ratio of how often the
#' coefficients' expected values were positive conditional on inclusion. (over
#' all visited models in case \code{exact=FALSE}, over the best models in case
#' \code{exact=TRUE})} \item{Column 'Idx'}{the original order of covariates as
#' the were used for sampling. (if included, the constant has index 0)}
#' 
#' @seealso \code{\link{bms}} for creating bma objects, \code{\link{pmp.bma}}
#' for comparing MCMC frequencies and marginal likelihoods.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords models
#' @examples
#' 
#' #sample, with keeping the best 200 models:
#' data(datafls)
#' mm=bms(datafls,burn=1000,iter=5000,nmodel=200)
#' 
#' #standard BMA PIPs and coefficients from the MCMC sampling chain, based on 
#' #  ...how frequently the models were drawn
#' coef(mm)
#' 
#' #standardized coefficients, ordered by index
#' coef(mm,std.coefs=TRUE,order.by.pip=FALSE)
#' 
#' #coefficients conditional on inclusion:
#' coef(mm,condi.coef=TRUE)
#' 
#' #same as
#' ests=coef(mm,condi.coef=FALSE)
#' ests[,2]/ests[,1]
#' 
#' #PIPs, coefficients, and signs based on the best 200 models
#' estimates.bma(mm,exact=TRUE)
#' 
#' #... and based on the 50 best models
#' coef(mm[1:50],exact=TRUE)
#' 
#' 
#' @method coef bma
#' @export
coef.bma <-function(object,exact = FALSE, order.by.pip = TRUE, include.constant = FALSE,
    incl.possign = TRUE, std.coefs = FALSE, condi.coef = FALSE, ...) {
  #just an alias for estimates.bma
  estimates.bma(object, exact=exact, order.by.pip = order.by.pip, include.constant = include.constant,
    incl.possign = incl.possign, std.coefs = std.coefs, condi.coef = condi.coef)
}

#' Tests for a 'bma' Object
#' 
#' tests for objects of class "bma"
#' 
#' 
#' @param bmao a 'bma' object: see 'value'
#' @return Returns \code{TRUE} if bmao is of class 'bma', \code{FALSE}
#' otherwise.
#' 
#' @seealso 'Output' in \code{\link{bms}} for the structure of a 'bma' object
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords classes
#' @examples
#' 
#'  data(datafls)
#'  mm=bms(datafls,burn=1000, iter=4000)
#'  is.bma(mm)
#' @export
is.bma <-function(bmao) {
  #returns true if the class of the object is either bma or a related class
   if (any(is.element(class(bmao),c("bma","bma.fcast","bma.sar","oldbma","bmav0")))) return(TRUE) else return(FALSE)
}




#' @export
is.topmod <-function(tmo) {
  #returns true if the class of the object is a "topmod" list
   if (is.element("topmod",class(tmo))) return(TRUE) else return(FALSE)
}



#' Gaussian Hypergeometric Function F(a,b,c,z)
#' 
#' Computes the value of a Gaussian hypergeometric function \eqn{ F(a,b,c,z) }
#' for \eqn{-1 \leq z \leq 1} and \eqn{a,b,c \geq 0}
#' 
#' The function \code{f21hyper} complements the analysis of the 'hyper-g prior'
#' introduced by Liang et al. (2008).\cr For parameter values, compare cf.
#' \url{https://en.wikipedia.org/wiki/Hypergeometric_function}.
#' 
#' @param a The parameter \code{a} of the Gaussian hypergeometric function,
#' must be a positive scalar here
#' @param b The parameter \code{b} of the Gaussian hypergeometric function,
#' must be a positive scalar here
#' @param c The parameter \code{c} of the Gaussian hypergeometric function,
#' must be a positive scalar here
#' @param z The parameter \code{z} of the Gaussian hypergeometric function,
#' must be between -1 and 1 here
#' @return The value of the Gaussian hypergeometric function \eqn{ F(a,b,c,z) }
#' @note This function is a simple wrapper function of sped-up code that is
#' intended for sporadic application by the user; it is neither efficient nor
#' general; for a more general version cf. the package '\code{hypergeo}'
#' 
#' @seealso package \code{hypergeo} for a more proficient implementation.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @references Liang F., Paulo R., Molina G., Clyde M., Berger J.(2008):
#' Mixtures of g-priors for Bayesian variable selection. J. Am. Statist. Assoc.
#' 103, p. 410-423
#' 
#' \url{https://en.wikipedia.org/wiki/Hypergeometric_function}
#' @keywords utilities
#' @examples
#' 
#'   
#'   f21hyper(30,1,20,.8) #returns about 165.8197
#'   
#'   f21hyper(30,10,20,0) #returns one
#'   
#'   f21hyper(10,15,20,-0.1) # returns about 0.4872972
#' @export
f21hyper = function(a,b,c,z) {
  if ((length(a)!=1)|(length(b)!=1)|(length(c)!=1)|(length(z)!=1)) stop("All function arguments need to be scalars")
  if ((a<0)|(b<0)|(c<0)) stop("Arguments a, b, and c need to be non-negative")
  if ((z>1)|(z<=(-1))) stop("Argument z needs to be between -1 and 1")
  nmx=max(100,3*floor(((a+b)*z-c-1)/(1-z)))
  if (nmx>10000) warning("Power series probably does not converge")
  serie=0:nmx
  return(1+sum(cumprod((a+serie)/(c+serie)*(b+serie)/(1+serie)*z)))
}


########################################################################
# Auxilary functions on the final bma object or its components #########
########################################################################

.post.constant <- function(X.data,Ebeta) {
    # calculates E(constant|Y): X.data is dataframe, Ebeta a vector
      Xmeans= colMeans(X.data)
     cons=Xmeans[1]-crossprod(Ebeta,Xmeans[-1])
     return(as.vector(cons))
}



.post.beta.draws <- function(topmods,reg.names,moment2=FALSE) {
      # constructs a nice matrix of the betas of the best models stored in topmods
      # topmods: topmod-object; reg.names: character vector (like colnames(X))
      # moment2: TRUE: return betas2(), FALSE: return betas()

      if(moment2) beta.draws=as.matrix(topmods$betas2()) else beta.draws=as.matrix(topmods$betas())
      if(sum(beta.draws)==0){
          stop("The tompod object provided does not have saved betas. cf. bbeta argument in function topmod")
        }

      
      
      
      if(nrow(beta.draws)!=length(reg.names)){
        rownames(beta.draws)=c(reg.names,"W-Index")
      }
      else{
        rownames(beta.draws)=c(reg.names)
      }
      beta.names=topmods$bool()

      if(length(which(beta.names=="0"))>0){
        colnames(beta.draws)=beta.names[-c(which(beta.names=="0"))]
      } else {
        colnames(beta.draws)=beta.names
      }
      return(beta.draws)
}



.post.topmod.includes <- function(topmods,reg.names) {
      # constructs nice 0-1-matrix with ionclusion vectors for the best models stored in topmods
      # topmods: topmodel object, reg.names: character-vector (like colnames(X))
      topmod = topmods$bool_binary()
      colnames(topmod)<-topmods$bool()
      rownames(topmod)=reg.names

    return(topmod)

}


.post.topmod.bma <- function(topmods,reg.names=numeric(0)) {
    # constructs nice 0-1 matrix for regressors in model, with the last two rows 
    #     being the analytical and MCMC PMPs
    # topmods: either the topmodel object or a "bma" object
    
    pmps = pmp.bma(topmods)
  if (is.bma(topmods)) {
    reg.names=topmods$reg.names; topmods=topmods$topmod
  }
  # constructs nice matrix combining the 0-1- includes matrix with info on PMP

  rbind(.post.topmod.includes(topmods,reg.names),t(pmps))
}




#' Model Binaries and their Posterior model Probabilities
#' 
#' Returns a matrix whose columns show which covariates were included in the best models in a 'bma' object. The last two columns detail posterior model probabilities. 
#' 
#' @param bmao an object of class 'bma' - see \code{\link{bma-class}}
#' @return Each column in the resulting matrix corresponds to one of the 'best' models in bmao: the first column for the best model, the second for the second-best model, etc.
#' The model binaries have elements 1 if the regressor given by the row name was included in the respective models, and 0 otherwise.
#' The second-last row shows the model's posterior model probability based on marginal likelihoods (i.e. its marginal likelihood over the sum of likelihoods of all best models)
#' The last row shows the model's posterior model probability based on MCMC frequencies (i.e. how often the model was accepted vs sum of acceptance of all models) 
#' Note that the column names are hexcode representations of the model binaries (e.g. "03" for c(0,0,0,1,0,0)) 
#' 
#' 
#' @details
#' Each bma class (the result of bms) contains 'top models', the x models with tthe best 
#' analytical likelihood  that bms had encountered while sampling
#' 
#' See \code{\link{pmp.bma}} for an explanation of likelihood vs. MCMC
#' frequency concepts
#' 
#' 
#' 
#' @seealso \code{\link{topmod}} for creating topmod objects, \code{\link{bms}}
#' for their typical use, \code{\link{pmp.bma}} for comparing posterior model
#' probabilities
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @examples
#' 
#' data(datafls)
#' #sample with a limited data set for demonstration
#' mm=bms(datafls[,1:12],nmodel=20)
#' 
#' #show binaries for all
#' topmodels.bma(mm)
#' 
#' #show binaries for 2nd and 3rd best model, without the model probs
#' topmodels.bma(mm[2:3])[1:11,]
#' 
#' #access model binaries directly
#' mm$topmod$bool_binary()
#'  
#' @export
topmodels.bma <-function(bmao)  {# function alias
  if (!is.bma(bmao)) {stop("you need to provide a bma object")}
  return(.post.topmod.bma(bmao))
}




.post.estimates <- function(b1mo=NULL,b2mo=NULL,cumsumweights=NULL,inccount=NULL,topmods=NULL,X.data=NULL,reg.names=NULL,pos.sign=NULL,exact=FALSE,order.by.pip=TRUE,include.constant=FALSE,incl.possign=TRUE,std.coefs=FALSE,condi.coef=FALSE) {
   # b1mo: (weighted) cumsum of first moment of beta
   # b2mo: (weighted) cumsum of second moment of beta
   # cumsumweights: sum of the weights/model probs - typically just the number of draws "iter"
   # inccount: (weighted) number of inclusions per coefficient
   # reg.names: character vector like colnames(X.data)
   # topmods: topmodl object (best n models)
   # exact: if True, then calcs posterior stats based on analytical PMPs (likelihoods) in topmods object; if False, then does it as a weighted average (typically with MCMC frequencies)
   # order.by.pip: if true then the matrix is sorted according to the PIPs (first column)
   # include.constant: if true, then add the constant in the last row of resultant matrix
   
   idx=1:(length(b1mo))

   if(exact){
      #calculate with analytical PMPs from topmods object (best n models)
      lt1=topmods$lik() - max(topmods$lik())    # do this to get only positive probabilities
      exact.pmp=as.vector(exp(lt1)/sum(exp(lt1)))
      pip=as.vector(topmods$bool_binary()%*%exact.pmp)
      idx=1:(length(pip))
      betas=topmods$betas()
      betas2=topmods$betas2()
      K=nrow(betas)
      Eb1=tcrossprod(betas,t(exact.pmp))[,1]   #multiply betas with corr. pmp's and sum up
      Eb2=tcrossprod(betas2,t(exact.pmp))[,1]         
      Ebsd=sqrt(Eb2-Eb1^2)           
      possign=round(tcrossprod((betas>0),t(exact.pmp))[,1]/pip,8)
      possign[is.nan(possign)]=NA
    }

    else{
      # Posterior mean and stand dev of each coefficient'
      # calculate by taking b1mo etc. as a weighted sum (typically MCMC draws), divided by the sum of the weights (typically number of runs - "iter")
      pip=inccount/cumsumweights
      Eb1 = b1mo/cumsumweights
      Eb2 = b2mo/cumsumweights
      Ebsd = sqrt(Eb2-Eb1^2)
      possign=round(pos.sign/inccount,8)
      possign[is.nan(possign)]=NA
    }

    if (include.constant) constterm=.post.constant(X.data,Eb1)

    if (condi.coef) {
      Eb1=Eb1/pip; Eb2=Eb2/pip; Ebsd=sqrt(Eb2-Eb1^2);
      Eb1[is.nan(Eb1)]=0; Ebsd[is.nan(Ebsd)]=0;
    }


    if (std.coefs) {
      #if standardized coefficients, then adjust them by SDs
      #important that this is done after all other Eb1 and Ebsd computations
      sddata=apply(as.matrix(X.data),2,stats::sd)
      Eb1=Eb1/sddata[1]*sddata[-1]
      Ebsd=Ebsd/sddata[1]*sddata[-1]   
      if (include.constant) constterm=constterm/sddata[1]
    }



    if (incl.possign) {
      post.mean<-cbind(pip,Eb1,Ebsd,possign,idx)
      rownames(post.mean)<-reg.names
      colnames(post.mean)<-c("PIP","Post Mean", "Post SD","Cond.Pos.Sign","Idx")
    } else {
      post.mean<-cbind(pip,Eb1,Ebsd,idx)
      #post.mean<-cbind(pip,Eb1,Ebsd)
      rownames(post.mean)<-reg.names
      colnames(post.mean)<-c("PIP","Post Mean", "Post SD","Idx")
      #colnames(post.mean)<-c("PIP","Post Mean", "Post SD")
    }





    if (order.by.pip) {
      post.mean<-post.mean[order(-post.mean[,1]),]   #order the results table according to PIP
    }

    if (include.constant) {
        constrow=matrix(c(1,constterm,NA,rep(NA,incl.possign),0),1)
        rownames(constrow)="(Intercept)"
        post.mean=rbind(post.mean,constrow)
    }

    return(post.mean)
}





#####################################################################################################

########################################################################
# functions on retrieving a FULL sys.call(0) ###########################
########################################################################




.construct.arglist <- function(funobj,envir=NULL) {
    # evaluates all function arguments at envir (default: calling environment) 
    # (they might have changed during running the function, or may be in variables)
    # construct.arglist gets rid of the variables and so on, and returns the argument list as a list
  
   namedlist=formals(funobj)
   argnames=names(namedlist)
   if (!is.environment(envir)) envir=sys.frame(-1)
  for (argn in 1:length(namedlist)) {
    # the following test is to cater for non.existent arguments
    testval=as.logical(try(exists(argnames[argn],envir=envir),silent=TRUE)); if (is.na(testval)) testval=FALSE
    if  (testval) {
          
          namedlist[[argn]]=try(get(argnames[argn],envir=envir))
    }
  }
  return(namedlist)
} 


#####################################################################################################

########################################################################
# functions for combining bma objects ##################################
########################################################################


.top10=function(nmaxregressors=10,nbmodels=10,bbeta=FALSE,lengthfixedvec=0,bbeta2=FALSE,...,
              inivec_lik=numeric(0),inivec_bool=character(0),inivec_count=numeric(0),inivec_vbeta=numeric(0), inivec_vbeta2=numeric(0),inivec_veck=0,inivec_fixvec=numeric(0)){
    #object used by bms to save the best models
    #set up the variables to be augmented in the process
    #use .top10(...) to initalise this object
    # top-level INPUT arguments:
    # nmaxregressors: maximum number of regreswsors possible in the models (ncol(X.data)-1) for bms()
    # nmodels: integer >=0, maximum number of models to store
    # bbeta: logical whether additionally, the models` beta coefficents should be stored as well
    # lengthfixedvec: integer >=0; a vector of fixed length =lengthfixedvec will be stored for each model (e.g. posterior moments of g, etc.);
    #                 if lengthfixedvec=0, nothing will be stored
    # inivec_...: for the advanced user: vectors to initialize the top10 object in case one wants to 'add' several models right at the beginning:
    #      inivec_lik: vector of likelihoods (length: nb of models) corresponds to: lik()
    #      inivec_bool: vector of model binaries in hexcode (length: nb of models); althernatively, a numeric vector of 1s and 0s (length: nmaxregressors)
    #                  correspoinds to: bool()
    #      inivec_count: vector of ncounts  (length: nb of models), corresponds to: ncount()
    #      inivec_vbeta: vector of betas (stripped of zeroes), corresponds to betas_raw()
    #      inivec_veck: vector of ks (length: nb of models), i.e. the number of coefs per model, corresponds to: kvec_raw()
    #      inivec_fixedvec: vector of fixed vectors (length: nb of models times lengthfixedvec)
    # OUTPUT: list with follwoing elements
    #  addmodel(mylik,vec01,vbeta=numeric(0),fixedvec=numeric(0)): Function that adds a model to be stored in the object.
    #       if the model is already there, only the model counter will be incremented: Reffeerence: see below
    #  lik(), bool(), ncount(), nbmodels, nregs, betas_raw(), kvec_raw(), bool_binary(), betas(), fixed_vector():
    #    these are final output functions, for a reference, refer to these directly


     findex=function() {
       seq_incl=seq_len(nbmodel)
       if (nbmodel==nbmodels) {
         seq_incl[indices]=seq_incl
       } else {
         truncindex=indices; truncindex[(nbmodel+1):nbmodels]=0L
         seq_incl[truncindex]=seq_incl
       }
       return(seq_incl)
     }


    betamat = function(top10_betavec) { # return a matrix: rows: betas per model (including zeros); columns: model
        bins=(sapply(as.list(top10_bool[findex()]),hexobject$as.binvec))

        betamatx=matrix(0,nmaxregressors,nbmodel)
        if (length(top10_betavec)>0) {betamatx[which(bins==1)]=top10_betavec} else betamatx=betamatx[,1]
        return(betamatx)
    }

#    sortall = function() {
#       if (!is_sorted) {
#          if (nbmodel<nbmodels) {
#             inclindex=seq_len(nbmodel)
#             indices <<- indices[inclindex];
#             top10_lik <<- top10_lik[inclindex];
#             top10_count <<- top10_count[inclindex];
#             top10_bool <<- top10_bool[inclindex];
#             top10_fixedvec <<- top10_count[seq_len(lengthfixedvec*nbmodel)];
#          }
#          top10_lik[indices]<<-top10_lik;
#          top10_count[indices]<<-top10_count;
#          top10_count[indices]<<-top10_count;
#
#       if (lengthfixedvec>0) { top10_fixvec <<- as.vector(matrix(lastm_fixvec,lengthfixedvec)[,order_index]) }
#       last_visited <<- c(lastnewmodel:1,length(lastm_lik):(lastnewmodel+1))[orderindex]
#       is_sorted <<- TRUE
#       }
#    }




    hexobject<-.hexcode.binvec.convert(nmaxregressors) #initialize object for hexcode to logical vector conversion
    if (nbmodels<0) {nbmodels=0}


    #declare needed objects in full length now to optimize memory allocation later
    indices=integer(nbmodels); top10_lik=rep(-Inf,nbmodels)
    top10_bool=character(nbmodels); top10_count=integer(nbmodels)
    top10_fixvec=numeric(lengthfixedvec*nbmodels)
    if (bbeta) lbetas=vector("list",nbmodels)
    if (bbeta2) lbetas2=vector("list",nbmodels)
    seq_nbmodel=seq_len(nbmodels); ix_of_mybool=logical(nbmodels)
    #is_sorted = FALSE

    # read in initial data
    nbmodel = length(inivec_lik)
    top10_lik[seq_len(nbmodel)]=inivec_lik; top10_count[seq_len(nbmodel)]=inivec_count;
    #read in initial binaries: character, or a list of vectors, a single vectors or a matrix, whose columns are the binaries
    if (is.character(inivec_bool)) {top10_bool[seq_len(nbmodel)]=inivec_bool} else {
      if (is.vector(inivec_bool)&(length(inivec_bool)==nmaxregressors)) {
        top10_bool[seq_len(nbmodel)]=hexobject$as.hexcode(inivec_bool)
      } else if (is.list(inivec_bool)) {
        top10_bool[seq_len(nbmodel)]=sapply(inivec_bool, hexobject$as.hexcode)
      } else if (is.matrix(inivec_bool)) {
        top10_bool[seq_len(nbmodel)]=sapply(as.list(as.data.frame(inivec_bool)), hexobject$as.hexcode)
      } else stop("inivec_bool is wrong format!")
    }
    top10_fixvec=inivec_fixvec;

    if (is.na(inivec_veck[1])) {inivec_veck=0}

    #read in initial beta information
    if (bbeta|bbeta2) {
      veck_ix=c(0,cumsum(inivec_veck))
      veckix_aux=as.list(seq_len(nbmodel)); veckix_aux=lapply(veckix_aux,function(x) { if (veck_ix[[x]]==veck_ix[[x+1]]) c(0,0) else c(veck_ix[[x]]+1,veck_ix[[x+1]]) } )
    }
    if (bbeta) { lbetas[seq_len(nbmodel)]=lapply(veckix_aux,function(x) inivec_vbeta[x[[1]]:x[[2]]]) } else lbetas=list(numeric(0))
    if (bbeta2) { lbetas2[seq_len(nbmodel)]=lapply(veckix_aux,function(x) inivec_vbeta2[x[[1]]:x[[2]]]) } else lbetas2=list(numeric(0))






    lastvec01=integer(nmaxregressors);
    modidx=length(top10_lik);

    indices[seq_len(nbmodel)]=order(inivec_lik,decreasing=TRUE)
    min.index = which.max(indices)

    if (length(min.index)>0) {
       min.top10_lik=top10_lik[[min.index]]
    } else {
      if (nbmodels>0) min.top10_lik=-Inf else min.top10_lik=Inf
    }


    index.of.mybool = function(mybool) {
      ix_of_mybool <<- (mybool==top10_bool)
    }
    check4dupl = index.of.mybool; dupl.possible=TRUE;

    retlist=list(
    addmodel=function(mylik,vec01,vbeta=numeric(0),vbeta2=numeric(0),fixedvec=numeric(0)) {
        #mylik: scalar likelihood, vec01: numeric vector of 0s and 1s, vbeta: small vector of betas (if bbeta was set to TRUE) that does NOT contain restricted zeros, vbeta2: small vector of betas^2 (is assumed to have same length as bvbeta)
        #use this function to add a model:
        #if its already among the best "nbmodels" models, its counter will be incremented by one
        #if it is not already in the best "nbmodels" models though its likelihood justifies that, the model will be added to the list


        if (mylik>=min.top10_lik|nbmodel<nbmodels) {

          #look whether the model is 'better' than the least model in the list
          if (identical(lastvec01,vec01)) {
              #if the model is the same as the model before, just adjust the counter
              top10_count[[modidx]]<<-top10_count[[modidx]]+1

          } else  {
                #look whether the model is already in the bestof list
                lastvec01<<-vec01
                mybool=hexobject$as.hexcode(vec01)

                check4dupl(mybool)
            if (!any(ix_of_mybool)) {

                #the model is not yet contained in the bestof list, but should be in there -> add model to list
                if (nbmodel<nbmodels) {
                    nbmodel <<- nbmodel+1
                    modidx <<- nbmodel
                } else {
                    modidx<<-min.index;
                }

                ltmylik=(top10_lik<=mylik) #adjusted on 2011-04-19, in order to cope with mylik=-Inf



                indices <<- indices + ltmylik
                indices[[modidx]] <<- nbmodels-sum(ltmylik)+1

                top10_lik[[modidx]] <<- mylik
                top10_bool[[modidx]] <<- mybool
                top10_count[[modidx]] <<- 1

                min.index <<- which.max(indices)
                min.top10_lik <<- top10_lik[[min.index]]

                if (lengthfixedvec>0) {
                    top10_fixvec[(modidx-1)*lengthfixedvec+seq_len(lengthfixedvec)] <<- fixedvec
                }


                if (bbeta) { lbetas[[modidx]] <<- vbeta }
                if (bbeta2) { lbetas2[[modidx]] <<- vbeta2 }





            } else {
                #the model is already contained in the bestof list -> just adjust counter
                modidx <<- seq_nbmodel[ix_of_mybool]
                top10_count[[modidx]]<<-top10_count[[modidx]]+1

            }


          }
        }
    },
    lik = function(){return(top10_lik[findex()])}, #return a vector of the best "nbmodels" likelihoods
    bool = function(){return(top10_bool[findex()])}, #return a vector of the best "nbmodels" codes as hexadecimal (e.g. model c(0,1,0,0) as "4")
    ncount = function(){return(top10_count[findex()])}, #return a vector of how each of the best models was chosen
    #counters = function(){return(c(tcalls,inlik,added,maintained))}, # for programming checks
    nbmodels = nbmodels, # return the maximum number of best mdoels saved in this object
    nregs = nmaxregressors, # return K, the maximum number of regressors overall
    betas_raw = function(){return(unlist(lbetas[findex()]))}, # return the vector of beta coefs. of the best models in one line without the zeros
    betas2_raw = function(){return(unlist(lbetas2[findex()]))}, # return the vector of beta^2 coefs. of the best models in one line without the zeros
    kvec_raw = function(){return(sapply(lbetas,length)[findex()])}, #return a vector that details how many coefs. each of the best models has
    bool_binary = function(){return(sapply(as.list(top10_bool[findex()]),hexobject$as.binvec))}, #return a matrix: each column: one of the best models; rows: logical of coeficcient inclusion
    betas = function() {
       betamat(unlist(lbetas[findex()]))
    },
    betas2 = function() {
       betamat(unlist(lbetas2[findex()]))
    },
    fixed_vector = function(){
        if (lengthfixedvec<=0) {return(matrix(0,0,0))} else
        {
        findex_base=(findex()-1)*lengthfixedvec
        findex_fixvec=numeric(0)
        for (xx in 1:lengthfixedvec) findex_fixvec=rbind(findex_fixvec,findex_base+xx)
        return(matrix(top10_fixvec[c(findex_fixvec)],lengthfixedvec))
        }
    }, # return the fixed vector that may contain additional statistics
    duplicates_possible = function(possible=NULL) {
      if (!is.logical(possible)) return(dupl.possible)
      if (possible) {
         check4dupl <<- index.of.mybool; dupl.possible <<- TRUE; ix_of_mybool <<- logical(nbmodels)
      } else {
         check4dupl <<- function(mybool){}; dupl.possible <<- FALSE; ix_of_mybool <<- FALSE
      }
    }
    )

    class(retlist)="topmod"
    return(retlist)
}



############################################################
#auxiliary functions for topmod object

#' @export
"[.topmod" <- function(tm,idx) {
# this function (applied as topmod[idx] ) provides a topmodel object with only the models indicated by idx
# e.g. topmod[1] contains only the best model, topmod[-(90:100)] eliminates the models ranked 90 to 100
if (any(is.na(suppressWarnings(as.integer(idx))))) idx=1:length(tm$lik())


  if (length(tm$betas_raw())>1) {
    bbeta=TRUE 
    bet=as.vector(tm$betas()[,idx])
    bet=bet[bet!=0]
  } else {
    bbeta=FALSE
    bet=numeric(0)
  }
  
  if (length(tm$betas2_raw())>1) {
    bbeta2=TRUE 
    bet2=as.vector(tm$betas2()[,idx])
    bet2=bet2[bet2!=0]                
  } else { 
    bbeta2=FALSE  
    bet2=numeric(0)
  }
  
  fixvec=tm$fixed_vector()
  if (!length(as.vector(fixvec))) fixvec=numeric(0) else fixvec=as.vector(t(fixvec[,idx]))
  
   .top10(nmaxregressors=tm$nregs,nbmodels=tm$nbmodels,bbeta=bbeta,lengthfixedvec=nrow(tm$fixed_vector()),bbeta2=bbeta2,
        inivec_lik=tm$lik()[idx],inivec_bool=tm$bool()[idx],inivec_count=tm$ncount()[idx],inivec_vbeta=bet, 
        inivec_vbeta2=bet2,inivec_veck=tm$kvec_raw()[idx],inivec_fixvec=fixvec)
}

#' @export
"[.bma" <- function(bmao,idx) {
# bma[idx] should have the same effect as applying the index to the topmod, for convenience
  bmao$topmod <- bmao$topmod[idx]
  return(bmao)
}

#' Printing topmod Objects
#' 
#' Print method for objects of class 'topmod', typically the best models stored
#' in a 'bma' object
#' 
#' See \code{\link{pmp.bma}} for an explanation of likelihood vs. MCMC
#' frequency concepts
#' 
#' @param x an object of class 'topmod' - see \code{\link{topmod}}
#' @param \dots additional arguments passed to \code{link{print}}
#' @return if \code{x} contains more than one model, then the function returns
#' a 2-column matrix: \item{Row Names}{show the model binaries in hexcode 
#' } \item{Column 'Marg.Log.Lik'}{shows the
#' marginal log-likelihoods of the models in \code{x}} \item{Column 'MCMC
#' Freq'}{shows the MCMC frequencies of the models in \code{x}}
#' 
#' if \code{x} contains only one model, then more detailed information is shown
#' for this model: \item{first line}{'Model Index' provides the model binary in
#' hexcode, 'Marg.Log.Lik' its marginal log likelhood, 'Sampled Freq.' how
#' often it was accepted (function \code{ncount()} in \code{\link{topmod}})}
#' \item{Estimates}{first column: covariate indices included in the model,
#' second column: posterior expected value of the coefficients, third column:
#' their posterior standard deviations (excluded if no coefficients were stored
#' in the topmod object - cf. argument \code{bbeta} in \code{\link{topmod}}) }
#' \item{Included Covariates}{the model binary} \item{Additional
#' Statistics}{any custom additional statistics saved with the model}
#' 
#' @seealso \code{\link{topmod}} for creating topmod objects, \code{\link{bms}}
#' for their typical use, \code{\link{pmp.bma}} for comparing posterior model
#' probabilities
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords print
#' @examples
#' 
#' # do some small-scale BMA for demonstration
#' data(datafls)
#' mm=bms(datafls[,1:10],nmodel=20)
#' 
#' #print info on the best 20 models
#' print(mm$topmod)
#' print(mm$topmod,digits=10)
#' 
#' #equivalent:
#' cbind(mm$topmod$lik(),mm$topmod$ncount())
#' 
#' 
#' 
#' #now print info only for the second-best model:
#' print(mm$topmod[2])
#' 
#' #compare 'Included Covariates' to:
#' topmodels.bma(mm[2])
#' 
#' #and to
#' as.vector(mm$topmod[2]$bool_binary())
#' 
#' 
#' @export
print.topmod <- function(x,...) {
#this function prints the matrix of logliks and MCMC frequencies
#if topmod contains only one model, more detiled infomration is given
# try e.g. topmod, or topmod[1:3] or topmod[1]
  tm=x
  if (length(tm$lik())==1) {
    infomat=c(tm$bool(), tm$lik(),tm$ncount())
    names(infomat)=c("Model Index","Marg.Log.Lik.","Sampled Freq.")
    print(infomat)    
    betamat=cbind(as.vector(tm$betas_raw()),sqrt(as.vector(tm$betas2_raw())-as.vector(tm$betas_raw())^2))
    if (nrow(betamat)!=0) {
     if (ncol(betamat)==1) {colnames(betamat)="Coef."} else {colnames(betamat)=c("Coef.","Std.Dev.")}
     rownames(betamat)=which(as.logical(as.vector(tm$bool_binary())))
     cat("\nEstimates:\n")
     print(betamat)
    }
    bin=as.vector(tm$bool_binary())
    names(bin)=1:length(bin)
    cat("\nIncluded Covariates:\n")
    print(bin)
    cat("\nAdditional Statistics:\n")    
    print(as.vector(tm$fixed_vector()))                  
  } else {
    mout=cbind(tm$lik(),tm$ncount())   
    colnames(mout)=c("Marg.Log.Lik", "MCMC Freq")
    rownames(mout)=tm$bool()           
    print(mout,...)                        
  }
}

#' Topmodel Object
#' 
#' Create or use an updateable list keeping the best x models it encounters
#' (for advanced users)
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
#' fourth regressor) and the marginal (log) likelihood \code{lik} (here: 4).
#' 
#' If this model is already part of \code{tmo}, then its respective
#' \code{ncount} entry is incremented by one; else it is inserted into a
#' position according to the ranking of \code{lik}.
#' 
#' In addition, there is the possibility to save (the first moments of)
#' coefficients of a model (\code{betas}) and their second moments
#' (\code{betas2}), as well as an arbitrary vector of statistics per model
#' (\code{fixed_vector}).
#' 
#' \code{is.topmod} returns \code{TRUE} if the argument is of class 'topmod'
#' 
#' @aliases topmod is.topmod
#' @param nbmodels The maximum number of models to be retained by the topmod
#' object
#' @param nmaxregressors The maximum number of covariates the models in the
#' topmod object are allowed to have
#' @param bbeta if \code{bbeta=TRUE}, then first and second moments of model
#' coefficients are stored in addition to basic model statistics (Note: if
#' \code{bbeta<0} then only the first moments are saved)
#' @param lengthfixedvec The length of an optional fixed vector adhering to
#' each model (for instance R-squared, etc). If \code{lengthfixedvec=0} then no
#' additonal fixed vector will be stored.
#' @param liks optional vector of log-likelihoods to initialize topmod object
#' with (length must be \code{<=nbmodels}) - see example below
#' @param ncounts optional vector of MCMC frequencies to initialize topmod
#' object with (same length as \code{liks}) - see example below
#' @param modelbinaries optional matrix whose columns detail model binaries to
#' initialize topmod object with (same nb columns as \code{liks}, nb rows as
#' \code{nmaxregressors}) - see example below
#' @param betas optional matrix whose columns are coefficients to initialize
#' topmod object with (same dimensions as \code{modelbinaries}) - see example
#' below
#' @param betas2 optional matrix whose columns are coefficients' second moments
#' to initialize topmod object with (same dimensions as \code{modelbinaries}) -
#' see example below
#' @param fixed_vector optional matrix whose columns are a fixed vector
#' initialize topmod object with (same \code{ncol} as \code{modelbinaries}) -
#' see example below
#' @return a call to \code{topmod} returns a list of class "topmod" with the
#' following elements:
#' \item{addmodel(mylik,vec01,vbeta=numeric(0),vbeta2=numeric(0),fixedvec=numeric(0))}{function
#' that adjusts the list of models in the 'topmod' object (see Details).
#' \code{mylik} is the basic selection criterion (usually log likelihood),
#' \code{vec01} is the model binary (logical or numeric) indicating which
#' regressors are included.\cr \code{vbeta} is a vector of length equal to
#' \code{sum(vec01)}, contianing only the non-zero coefficients (only accounted
#' for if \code{bbeta!=FALSE}). \code{vbeta2} is a similar vector of second
#' moments etc. (only accounted for if \code{bbeta=TRUE}); \code{fixedvec} is
#' an arbitrary vector of length \code{lengthfixedvec} (see above)}
#' \item{lik()}{A numeric vector of the best models (log) likelihoods, in
#' decreasing order} \item{bool()}{A character vector of hexmode expressions
#' for the model binaries (cf. \code{\link{bin2hex}}), sorted by \code{lik()} }
#' \item{ncount()}{A numeric vector of MCMC frequencies for the best models
#' (i.e. how often the respective model was introduced by \code{addmodel})}
#' \item{nbmodels}{Returns the argument \code{nbmodel}} \item{nregs}{Returns
#' the argument \code{nmaxregressors}} \item{bool_binary()}{Returns a matrix
#' whose columns present the models conforming to \code{lik()} in binary form}
#' \item{betas()}{a matrix whose columns are the coefficients conforming to
#' \code{bool_binary()} (Note that these include zero coefficients due to
#' non-inclusion of covariates); Note: if \code{bbeta=FALSE} this returns an
#' empty matrix} \item{betas2()}{similar to \code{betas} , for the second
#' moments of coefficients Note: if \code{bbeta<=0}, this returns an empty
#' matrix} \item{fixed_vector()}{The columns of this matrix return the
#' \code{fixed_vector} statistics conforming to \code{lik()} (see Details);
#' Note: if \code{lengthfixedvec=0} this returns an empty matrix}
#' @note \code{topmod} is rather intended as a building block for programming;
#' it has no direct application for a user of the BMS package.
#' 
#' @seealso the object resulting from \code{\link{bms}} includes an element of
#' class 'topmod'
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords models
#' @examples
#' 
#' #standard use
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
#'   is.topmod(tm)
#'   
#'   #extract a topmod oobject only containing the second best model
#'   tm2=tm[2]
#'   
#'   
#'   
#'   #advanced: should return the same result as
#'   #initialize
#'   tm2= topmod(2,4,TRUE,0, liks = c(-2.2,-2.3), ncounts = c(2,1), 
#'           modelbinaries = cbind(c(0,1,1,1),c(1,1,1,1)), betas = cbind(0:3,1:4), 
#'           betas2 = cbind(c(0,5:7),5:8)) 
#' 
#'   #update 
#'   tm$addmodel(-2.5,c(1,0,0,1),1:2,5:6) #add another model
#'   
#'   #read out
#'   tm$lik()
#'   tm$ncount()
#'   tm$bool_binary()
#'   tm$betas()
#' 
#' @export
topmod <- function(nbmodels,nmaxregressors=NA,bbeta=FALSE,lengthfixedvec=0,liks=numeric(0),
                   ncounts=numeric(0),modelbinaries=matrix(0,0,0),betas=matrix(0,0,0),
                   betas2=matrix(0,0,0),fixed_vector=matrix(0,0,0)) {
    #user-friendly function to create a 'topmod' object
    #nbmodels: integer, maxmium number of models to be retained
    #nmaxregressors: integer, maximum possible number of covariates (optional if arguments modelbinaries or betas are provided)
    #bbeta: If TRUE, model coefficients (first & second moments) are to saved by the resulting object; If FALSE not; (If bbeta<0 then only the first mometns are to be saved)
    #lengthfixedvec: If lengthfixedvec>0 then the resulting topmod object saves a numeric of length 'lengthfixedvec' for each added model; if lengthfixedvec=0, it does not; can be omitted if fixed_vector is supplied
    #
    #optional arguments to initialize the topmod object directly
    #liks: a vector of (log) likelihoods (does not need to be sorted)
    #ncounts: an (optional) vector of MCMC frrequencies for each model (if provided, needs to have the same length as liks)
    #modelbinaries: a logical or numeric matrix denoting the models to initialize the topmod object with; With nrow(modelbinaries) equal to nmaxregressors and ncol(modelbinaries) equal to length(liks) (can be omitted if argument betas is provided)
    #betas: a numeric matrix with the (expecteed values of) coefficents of the models; With nrow(betas) equal to nmaxregressors and ncol(betas) equal to length(liks); Note that providing betas automatically sets bbeta=-1 
    #betas2: a numeric matrix with the second moments coefficents of the models; With nrow(betas2) equal to nmaxregressors and ncol(betas2) equal to length(liks); Note that providing betas automatically sets bbeta=-1 
    #fixed_vector: an optional fixed vector to be saved with each model
    
    if (!is.numeric(nbmodels)) stop("argument 'nbmodels' needs to be an integer>0")
    nbmodels=floor(nbmodels[[1]])
    if (nbmodels[[1]]<0) stop("argument 'nbmodels' needs to be an integer>0")

    if (bbeta>0) bbeta2=TRUE else bbeta2=FALSE
    bbeta=as.logical(bbeta)

    if (!bbeta&(length(betas)>0)) bbeta=TRUE
    if (!bbeta2&(length(betas2)>0)) bbeta2=TRUE

    if (is.na(lengthfixedvec[1])) lengthfixedvec=0
    if ((lengthfixedvec==0)&length(fixed_vector)>0) {lengthfixedvec=nrow(fixed_vector)}
    
    if (length(liks)>nbmodels) stop("liks longer than nbmodels allows")
    if (length(ncounts)>nbmodels) stop("ncounts longer than nbmodels allows")

    if ((length(modelbinaries)==0)&(length(betas)>0)) {modelbinaries=as.logical(betas);dim(modelbinaries)=dim(betas) }
    if (ncol(modelbinaries)>nbmodels) stop("modelbinaries has more columns than than nbmodels allows")
    bindim=dim(modelbinaries); modelbinaries=as.logical(modelbinaries); dim(modelbinaries)=bindim;
    
    if ((is.na(nmaxregressors[1]))&(length(modelbinaries)>0)) {nmaxregressors=nrow(modelbinaries) }
    if (is.na(nmaxregressors)) stop("argument 'nmaxregressors' is missing")
    nmaxregressors=floor(nmaxregressors[[1]])
    if (nmaxregressors<=0) stop("argument 'nmaxregressors' needs to be a positive integer")

    if ((length(ncounts)==0)&(length(liks)>0)) {ncounts=rep(1,length(liks)) }

    #check conformance with K
    if (length(modelbinaries)>0) if (nmaxregressors!=nrow(modelbinaries)) stop("nrow() of modelbinaries does not match nmaxregressors")
    if (bbeta&(length(betas)>0)) if (nmaxregressors!=nrow(betas)) stop("nrow() of betas does not match nmaxregressors")
    if (bbeta2&(length(betas2)>0)) if (nmaxregressors!=nrow(betas2)) stop("nrow() of betas2 does not match nmaxregressors")
    
    #check that all have the same model length
    N=length(liks)
    if (length(ncounts)!=length(liks)) stop("lengths of arguments 'liks' and 'ncounts' do not conform")
    if (ncol(modelbinaries)!=length(liks)) stop("nrow of argument 'modelbinaries' does not conform to length of argument 'liks'")
    if (bbeta) if (ncol(betas)!=length(liks)) stop("nrow of argument 'betas' does not conform to length of argument 'liks'")
    if (bbeta2) if (ncol(betas2)!=length(liks)) stop("nrow of argument 'betas2' does not conform to length of argument 'liks'")
    if (lengthfixedvec) if (ncol(fixed_vector)!=length(liks)) stop("nrow of argument 'fixed_vector' does not conform to length of argument 'liks'")


    morder=order(liks,decreasing=TRUE)
    liks=liks[morder]
    modelbinaries=modelbinaries[,morder,drop=FALSE]
    ncounts=ncounts[morder]
    if (bbeta) { betas=betas[,morder,drop=FALSE] }
    if (bbeta2) { betas2=betas2[,morder,drop=FALSE] }
    if (lengthfixedvec) { fixed_vector=fixed_vector[,morder,drop=FALSE] }

    hexobj=.hexcode.binvec.convert(nmaxregressors)
    bools=as.vector(sapply(as.list(as.data.frame(modelbinaries)),hexobj$as.hexcode))
    if (length(bools)==0) {bools=character(0)}

    veck=numeric(0); betas_raw=numeric(0); betas2_raw=numeric(0);
    if (bbeta&(length(bbeta)>0)) {
      veck=colSums(modelbinaries)
      betas_raw=as.vector(betas)[as.vector(modelbinaries)]
    }

    if (bbeta2&(length(bbeta2)>0)) {
      veck=colSums(modelbinaries)
      betas2_raw=as.vector(betas2)[as.vector(modelbinaries)]
    }

      fixedvec=as.vector(fixed_vector)

      .top10(nmaxregressors=nmaxregressors,nbmodels=nbmodels,bbeta=bbeta,lengthfixedvec=lengthfixedvec,bbeta2=bbeta2,
              inivec_lik=liks,inivec_bool=bools,inivec_count=ncounts,inivec_vbeta=betas_raw, inivec_vbeta2=betas2_raw,inivec_veck=veck,inivec_fixvec=fixedvec)
    
}

.cor.topmod <- function(tmo) {
    if (is.bma(tmo)) tmo=tmo$topmod
    pmp.10=pmp.bma(tmo,oldstyle=TRUE)
      if (nrow(pmp.10)==1|suppressWarnings(length(grep("error",class(try(cor(pmp.10[,1],pmp.10[,2]),silent=TRUE)))))) {
         corr.pmp=NA
      } else {
        if (var(pmp.10[,2])==0) corr.pmp=NA else corr.pmp=cor(pmp.10[,1],pmp.10[,2])
      }
    return(corr.pmp)
}

.topmod.as.bbetaT <- function (tm,gprior.info=NULL,yXdata=NULL, addr2=FALSE) {
       #this is a small function to convert a topmod object with bbeta=FALSE into on ewith bbeta=TRUE
       #CAUTION: this is not a generic function, but tailored to very specific topmod objects in conjunction with bms()
       # tm: a topmod object as the one resulting from a call to bms()
       # gprior.info: an object as the one resulting from a call to bms(); 
       #   optionally you can easily create one with the function .choose.gprior()
       # yXdata: the original data, just like the object X.data resulting from a call to bms()
       # if addr2=TREU, then the fixed_vector is appended to a first row containing the R-squareds
       # returns amn adjusted topmod object
       #
       # Optionally, tm can be a bma object, then this function returns an adjusted bma object
       
       is.bmao=FALSE
       if (is.bma(tm)) {
         #in case tm is a bma object...
         is.bmao=TRUE
         bmao=tm;
         yXdata=bmao$arguments$X.data; 
         gprior.info=bmao$gprior.info; 
         tm =bmao$topmod 
       }

      #retrieve necessary info
       yXdata=as.matrix(yXdata); 
       N=nrow(yXdata); K=ncol(yXdata)-1
       yXdata=yXdata-matrix(colMeans(yXdata),N,K+1,byrow=TRUE)
       if (length(tm$lik())<1) { if (is.bmao) return(bmao) else return(tm)}
       if (!addr2) if ((length(tm$betas_raw())>0)&(ncol(as.matrix(tm$betas()))==length(tm$lik()))) { if (is.bmao) return(bmao) else return(tm)}
       
       bools=(tm$bool_binary())
       yty=c(crossprod(yXdata[,1]))
       positions=lapply(lapply(as.list(as.data.frame(bools)),as.logical),which)
       
      #calculate all the OLS results for all the models
       olsmodels=lapply(lapply(positions,.ols.terms2,yty=yty,N=N,K=K,XtX.big=crossprod(yXdata[,-1]),Xty.big=c(crossprod(yXdata[,-1],yXdata[,1]))),function (x) x$full.results())
      
      #initialize the right gprior-function
#        if (gprior.info$is.constant) {
#          lprobo=.lprob.constg.init(g=gprior.info$g,N=N,K=K,yty=yty)
#        } else if (gprior.info$gtype=="EBL") {
#          lprobo=.lprob.eblocal.init(g=gprior.info$g,N=N,K=K,yty=yty)
#        } else if (gprior.info$gtype=="hyper") {
#          lprobo=.lprob.hyperg.init(g=gprior.info$g,N=N,K=K,yty=yty,f21a=gprior.info$hyper.parameter)
#        } else {
#           stop("gprior not recognizeable")
#        }
      lprobo=gprior.info$lprobcalc
       
      #calculate the posterior statistics
       lpl=lapply(olsmodels,function(x) lprobo$lprob(x$ymy,length(x$bhat),x$bhat,x$diag.inverse)) # caution: lprobs are not the same due to missing prior model probs
       veck=as.vector(unlist(lapply(lapply(lpl,"[[","b1new"),length)))
       b1raw=as.vector(unlist(lapply(lpl,"[[","b1new")))
       b2raw=as.vector(unlist(lapply(lpl,"[[","b2new")))
       fixedvecmat=tm$fixed_vector()
       if (addr2) { # add R-squared as the first elements among the fixed vectors
        r2=1-sapply(olsmodels,function(x) x$ymy)/yty
        if (nrow(fixedvecmat)==0) {
          fixedvecmat=matrix(0,0,length(veck))
        } else if ( mean(abs(r2-fixedvecmat[1,]))< 1e-17 ) {
          fixedvecmat=fixedvecmat[-1,,drop=FALSE]
        }
        fixedvecmat = rbind(r2,fixedvecmat)
       }
       lengthfixedvec=nrow(fixedvecmat)
       
      #now create new topmod with bbeta=TRUE
       tm<-.top10(nmaxregressors=tm$nregs,nbmodels=tm$nbmodels,bbeta=TRUE,lengthfixedvec=lengthfixedvec,bbeta2=TRUE,
              inivec_lik=tm$lik(),inivec_bool=tm$bool(),inivec_count=tm$ncount(),inivec_vbeta=b1raw, inivec_vbeta2=b2raw,inivec_veck=veck,inivec_fixvec=c(fixedvecmat))
              
              
              
       if (is.bmao) {
         bmao$topmod <- tm
#         bmao$arguments$beta.save=TRUE
         return(bmao)
       }
       return(tm)
}




combine_chains <- function(...) {
    #to combine outputs of the the function bms
    #EXAMPLE:
    #bma1<-bms(X.data=t5.within[,1:20],burn=100,iter=1000,g=TRUE,nmodel=10,logfile=TRUE,beta.save=FALSE,start.value=41,step=1000)
    #bma2<-bms(X.data=X.data,burn=1000,iter=10000,g=TRUE,nmodel=10,logfile=TRUE,beta.save=FALSE,start.value=41,step=1000)
    #out=combine_chains(bma1,bma2)
    ### or: combine_chains(bma1,bma2,bma3,bma4)...
    # output is a standard bma object


  combine_topmods <- function(topmodobj1,topmodobj2) {
    #to combine top10 models objects of function bms()
    #e.g. ppp=combine_topmods(test1$topmod,test2$topmod)
    # output: a topmodel object

    #retrieve the necessary properties
    nregs1=topmodobj1$nregs
    nregs2=topmodobj2$nregs
    if (nregs1!=nregs2) {stop("The number of regressors in both BMA chains has to be the same!")}
    k1=length(topmodobj1$ncount())
    k2=length(topmodobj2$ncount())

    #read out the importnant stuff
    nbmodels1=topmodobj1$nbmodels
    nbmodels2=topmodobj2$nbmodels
    ncount1=topmodobj1$ncount()
    ncount2=topmodobj2$ncount()
    lik1=topmodobj1$lik()
    lik2=topmodobj2$lik()
    bool1=topmodobj1$bool()
    bool2=topmodobj2$bool()
    betas1=topmodobj1$betas()
    betas2=topmodobj2$betas()
    betas2_1=topmodobj1$betas2()    
    betas2_2=topmodobj2$betas2()    
    fv1=topmodobj1$fixed_vector()
    fv2=topmodobj2$fixed_vector()
    if (all(betas1==0)|all(betas2==0)) {dobetas=FALSE} else {dobetas=TRUE}
    if (all(betas2_1==0)|all(betas2_2==0)) {dobetas2=FALSE} else {dobetas2=TRUE}
    
    #first look which models of 1 are already in 2 and
    #for these just update the ncounts (note: this is quite easy, since this subset in 1 has the same order as in 2)
    idxin2_boolof1in2=match(bool1,bool2)
    idxin1_boolof1in2=which(!is.na(idxin2_boolof1in2))
    idxin2_boolof1in2=idxin2_boolof1in2[!is.na(idxin2_boolof1in2)]
    ncount2[idxin2_boolof1in2]=ncount2[idxin2_boolof1in2]+ncount1[idxin1_boolof1in2]

    if (any(idxin1_boolof1in2)) { # in case there are models in 1 that also show up in 2
        #strip 1 of all the models that were already in 2
        ncount1=ncount1[-idxin1_boolof1in2]
        lik1=lik1[-idxin1_boolof1in2]
        bool1=bool1[-idxin1_boolof1in2]
    }
    #now do A u (B\(AnB))
    lika=c(lik2,lik1)
    orderlika=order(lika,decreasing=TRUE)
    lika=lika[orderlika]
    ncounta=c(ncount2,ncount1)[orderlika]
    boola=c(bool2,bool1)[orderlika]

    if (dobetas) {
       # if there are betas, do the same for betas
        if (any(idxin1_boolof1in2)) betas1=betas1[,-idxin1_boolof1in2]
        betasa=cbind(betas2,betas1)[,orderlika]
        betasa_not0=betasa!=0
        vecka=colSums(betasa_not0)
        vbetaa=as.vector(betasa[as.vector(betasa_not0)])
    } else {
        vecka=0;vbetaa=numeric(0)
    }


    if (dobetas2) {                                  
       # if there are betas^2, do the same for betas^2  
        if (any(idxin1_boolof1in2)) betas2_1=betas2_1[,-idxin1_boolof1in2] 
        betasa2=cbind(betas2_2,betas2_1)[,orderlika]                        
        vbetaa2=as.vector(betasa2[as.vector(betasa_not0)]);            
    } else {                                        
        vbetaa2=numeric(0)                            
    }                                               
    
    
    fva=numeric(0); lfv=0;
    if ( (nrow(fv1)==nrow(fv2)) & ((nrow(fv1)>0) & (nrow(fv2)>0)) ) {
        # if there is a fixed vector then combine it the same way    
        if (any(idxin1_boolof1in2)) fv1=fv1[,-idxin1_boolof1in2]
        if (!is.matrix(fv1)) fv1=t(fv1)
        fva=as.vector(cbind(fv2,fv1)[,orderlika])
        lfv=nrow(fv1)
    }
    return(.top10(nmaxregressors=nregs1,nbmodels=length(lika),bbeta=dobetas,lengthfixedvec=lfv,bbeta2=dobetas2,inivec_lik=lika,inivec_bool=boola,inivec_count=ncounta,inivec_vbeta=vbetaa,inivec_vbeta2=vbetaa2,inivec_veck=vecka,inivec_fixvec=fva))
}


  combine_2chains <- function(flso1,flso2) {
#      if (!exists("sPath")) sPath=""
#      source(paste(sPath,"aux_inner.r",sep=""),local=TRUE)

 #   if (flso1$arguments$beta.save & flso2$arguments$beta.save) {
 #       beta.save = TRUE
 #   } else {
 #       beta.save = FALSE
 #   }

    # combine the topmodel opbjects
    topmod.combi=combine_topmods(flso1$topmod,flso2$topmod)
    
    #prepare the gprior.info obejct for post.calc
    gpi <- flso1$gprior.info
    gpi$shrinkage.moments=numeric(length(gpi$shrinkage.moments))    

    # use post.calc to compute info, gprior.info, and reg.names
    io1=flso1$info; io2 = flso2$info
    obj.combi=.post.calc(gprior.info=gpi,add.otherstats=io1$add.otherstats + io2$add.otherstats,k.vec=(io1$k.vec[-1]+io2$k.vec[-1]),null.count=(io1$k.vec[1]+io2$k.vec[1]),
       flso1$arguments$X.data,topmods=topmod.combi,b1mo=io1$b1mo + io2$b1mo,b2mo=io1$b2mo + io2$b2mo,iter=io1$iter + io2$iter,burn=io1$burn + io2$burn,
       inccount=io1$inccount + io2$inccount,models.visited=io1$models.visited + io2$models.visited,K=io1$K,N=io1$N,msize=io1$msize + io2$msize,
       timed=io1$timed + io2$timed,cumsumweights=io1$cumsumweights + io2$cumsumweights,mcmc=flso1$arguments$mcmc,possign=io1$pos.sign+io2$pos.sign)
    
    # concatenate start.positions by cbinding them   
    stpos1=as.matrix(flso1$start.pos);stpos2=as.matrix(flso2$start.pos)
    startpos.combi =cbind(rbind(stpos1,matrix(0,max(0,nrow(stpos2)-nrow(stpos1)),ncol(stpos1))),rbind(stpos2,matrix(0,max(0,nrow(stpos1)-nrow(stpos2)),ncol(stpos2))))

    # concatenate bms.calls in list
    call.combi= c(flso1$bms.call,flso2$bms.call)

    # combine arguments
    args.combi = flso1$arguments; args2=flso2$arguments
    args.combi$burn = args.combi$burn + args2$burn
    args.combi$iter = args.combi$iter + args2$iter    
    if ((length(args.combi$mprior.size)==1)|(length(args.combi$mprior.size)==1)) {args.combi$mprior.size = mean(c(args.combi$mprior.size,args2$mprior.size))}
    args.combi$nmodel = topmod.combi$nbmodels
#    args.combi$beta.save = (args.combi$beta.save & args2$beta.save)
    args.combi$user.int = (args.combi$user.int & args2$user.int)
    args.combi$g.stats = (args.combi$g.stats & args2$g.stats)    

    #model prior object
    mp1=flso1$mprior.info; mp2=flso2$mprior.info
    if (mp1$mp.mode!=mp2$mp.mode) {mpall=list()} else {
       mpall=mp1
       mpall$mp.msize = .5*mp1$mp.msize+.5*mp2$mp.msize
       mpall$origargs$mpparam = .5*mp1$origargs$mpparam+.5*mp2$origargs$mpparam
       mpall$mp.Kdist= .5*mp1$mp.Kdist + .5*mp2$mp.Kdist    
    }

      result=list(info=obj.combi$info,arguments=args.combi, topmod=topmod.combi,start.pos=startpos.combi,
      gprior.info=obj.combi$gprior.info, mprior.info=mpall, X.data=flso1$arguments$X.data,reg.names=obj.combi$reg.names,bms.call=call.combi)
    class(result)="bma"
    
    return(result)
  }
  #############################################################################################################
  #this is the rest of the combine function; the combine function is iteratively used to combine as many chains
  #as are specified by (...)
  ############################################################################################################
    arglist=list(...)
                                                 
    if ( !all(unlist(lapply(arglist,is.bma)))) stop("All of the input arguments must be BMA objects!")
    if ( !all(lapply(arglist,function (xx) xx$info$K)==arglist[[1]]$info$K) ) stop("All of the input BMA objects must have an equal number of max regressors (i.e. equal (X.data))!")
    if ( !all(lapply(arglist,function (xx) xx$info$N)==arglist[[1]]$info$N) ) stop("All of the input BMA objects must have equal X.data!")
    if ( !all(lapply(arglist,function (xx) xx$gprior.info$gtype)==arglist[[1]]$gprior.info$gtype) ) stop("All of the input BMA objects must have the same type of g-prior (bms-argument g)")
    if ( length(arglist)==1) return(arglist[[1]])

    combined_output <- combine_2chains(arglist[[1]],arglist[[2]])
  if (nargs()>2) {
    for (inarg in 3:nargs()) {
        combined_output <- combine_2chains(arglist[[inarg]],combined_output)
    }
  }
  ############################################################################################################
  return(combined_output)
}

#' Concatenate bma objects
#' 
#' Combines bma objects (resulting from \code{\link{bms}}). Can be used to
#' split estimation over several machines, or combine the MCMC results obtained
#' from different starting points.
#' 
#' Aggregates the information obtained from several chains. The result is a
#' 'bma' object (cf. 'Values' in \code{\link{bms}}) that can be used just as a
#' standard 'bma' object.\cr Note that \code{combine_chains} helps in
#' particular to paralllelize the enumeration of the total model space: A model
#' with \eqn{K} regressors has \eqn{2^K} potential covariate combinations: With
#' \eqn{K} large (more than 25), this can be pretty time intensive.  With the
#' \code{\link{bms}} arguments \code{start.value} and \code{iter}, sampling can
#' be done in steps: cf. example 'enumeration' below.
#' 
#' @aliases combine_chains c.bma
#' @param \dots At least two 'bma' objects (cf. \code{\link{bms}})
#' @param recursive retained for compatibility with \code{\link{c}} method
#' 
#' @seealso \code{\link{bms}} for creating bma objects
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords models
#' @examples
#' 
#'  data(datafls)
#'   
#'  #MCMC case ############################
#'  model1=bms(datafls,burn=1000,iter=4000,mcmc="bd",start.value=c(20,30,35))
#'  model2=bms(datafls,burn=1500,iter=7000,mcmc="bd",start.value=c(1,10,15))
#'  
#'  model_all=c(model1,model2)
#'  coef(model_all)
#'  plot(model_all)
#'  
#'  
#'  
#'  #splitting enumeration ########################
#'  
#'  #standard case with 12 covariates (4096 differnt combinations):
#'  enum0=bms(datafls[,1:13],mcmc="enumerate")
#'  
#'  # now split the task:
#'  # enum1 does everything from model zero (the first model) to model 1999
#'  enum1=bms(datafls[,1:13],mcmc="enumerate",start.value=0,iter=1999)
#'  
#'  # enum2 does models from index 2000 to the index 3000 (in total 1001 models)
#'  enum2=bms(datafls[,1:13],mcmc="enumerate",start.value=2000,iter=1000)
#'  
#'  # enum3 does models from index 3001 to the end
#'  enum3=bms(datafls[,1:13],mcmc="enumerate",start.value=3001)
#'  
#'  enum_combi=c(enum1,enum2,enum3)
#'  coef(enum_combi)
#'  coef(enum0)
#'  #both enum_combi and enum0 have exactly the same results 
#'  #(one difference: enum_combi has more 'top models' (1500 instead of 500))
#' 
#' @export
c.bma <- function(...,recursive=FALSE) {
  #simple wrapper, recursive has no meaning and ist retained for compatibility
  combine_chains(...) 
}


########################################################################
# auxiliary functions for the topmodel object ##########################
########################################################################


.hexcode.binvec.convert <- function(length.of.binvec) {
    #function to initialise conversion betwwen logical vector (such as c(1,0,0,0)) and character hexcode vector (such as "f")
    #length.of.binvec is the desired length of the inserted and resulting logical vectors; 
    #the initialisation will fit some leading zeros to make it convertible into hexcode (length of bin. vector as a multiple of 4)
    if (length(length.of.binvec)>1) length.of.binvec=length(length.of.binvec)
    addpositions=4-length.of.binvec%%4; positionsby4=(length.of.binvec+addpositions)/4;
    hexvec=c(0:9,"a","b","c","d","e","f"); #lookup list for converting from binary to hexadecimal
    hexcodelist=list("0"=numeric(4),"1"=c(0,0,0,1),"2"=c(0,0,1,0),"3"=c(0,0,1,1),"4"=c(0,1,0,0),"5"=c(0,1,0,1),"6"=c(0,1,1,0), 
      "7"=c(0,1,1,1),"8"=c(1,0,0,0),"9"=c(1,0,0,1),"a"=c(1,0,1,0),"b"=c(1,0,1,1),"c"=c(1,1,0,0),"d"=c(1,1,0,1),"e"=c(1,1,1,0),"f"=c(1,1,1,1));
      #lookup list for converting from hexadecimal to binary
    
    return(list(
    as.hexcode = function(binvec) {
        #convert logical vector to hexcode character
        incl=c(numeric(addpositions),binvec);dim(incl)=c(4,positionsby4); #split into elements of four positions
        return(paste(hexvec[crossprod(incl,2L^(3:0))+1],collapse=""));
    },
    as.binvec = function(hexcode) {
        #convert hexcode character to numeric vector (e.g. "a" to c(1,0,1,0))
        return(unlist(hexcodelist[unlist(strsplit(hexcode,"",fixed=TRUE),recursive=FALSE,use.names=FALSE)],recursive=FALSE,use.names=FALSE)[-(1:addpositions)])    
    }))
}


#' @rdname bin2hex
#' @export
hex2bin<-function(hexcode) {
    #user-friendly function to convert some hexcode character to numeric vector (e.g. "a" to c(1,0,1,0))
    if (!is.character(hexcode)) stop("please input a character like '0af34c'");
    hexcode <- paste("0",tolower(hexcode),sep="")
    hexobj<-.hexcode.binvec.convert(length(hexcode)*16L)
    return(hexobj$as.binvec(hexcode))
}


#' Converting Binary Code to and from Hexadecimal Code
#' 
#' A simple-to-use function for converting a logical ('binary') vector into hex
#' code and reverse.
#' 
#' The argument is an integer in binary form (such as "101"), provided as a
#' logical (\code{c(T,F,T)}) or numeric vector (\code{c(1,0,1)}).\cr
#' \code{bin2hex} then returns a character denoting this number in hexcode (in
#' this case "5").
#' 
#' The function \code{hex2bin} does the reverse operation, e.g.
#' \code{hex2bin("5")} gives (\code{c(1,0,1)}).
#' 
#' @aliases bin2hex hex2bin
#' @param binvec a logical vector (alternatively a vector coercible into
#' logical)
#' @param hexcode a single-element character denoting an integer in hexcode
#' (admissible character: 0 to 9, ato f)
#' @return \code{bin2hex} returns a single element character; \code{hex2bin}
#' returns a numeric vector equivalent to a logical vector
#' 
#' @seealso \code{\link{hex2bin}} for converting hexcode into binary vectors,
#' \code{\link{format.hexmode}} for a related R function.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords arith
#' @examples
#' 
#'   bin2hex(c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE))
#'   bin2hex(c(1,0,1,0,1,1))
#'   hex2bin("b8a")
#'   bin2hex(hex2bin("b8a"))
#' 
#' @export
bin2hex<-function(binvec) {
    #user-friendly function to convert some logical vector to hexcode character(e.g. c(1,0,1,0) or c(T,F,T,F) to "a")
    if (!is.logical(binvec)) {if (is.numeric(binvec)) {binvec=as.logical(binvec)} else {stop("need to supply a logical vector like c(T,F) or c(1,0)")} }
    hexobj<-.hexcode.binvec.convert(length(binvec))
    hexcode=hexobj$as.hexcode(binvec)    
    if (nchar(hexcode)>(floor((length(binvec)-1)/4)+1)) {hexcode=substring(hexcode,2)}
    return(hexcode)    
}


 #######################
 # FUNCTIONS FOR USERS #  
#########################################################################

#' Plot Model Size Distribution
#' 
#' Plots posterior and prior model size distribution
#' 
#' 
#' @param bmao a 'bma' object (cf. \code{\link{bms}})
#' @param exact if \code{TRUE}, then the posterior model distribution is based
#' on the best models of \code{bmao} and their marginal likelihoods;\cr if
#' \code{FALSE} (default) then the distribution is based on all encountered
#' models and their MCMC frequencies (cf. 'Details' in \code{\link{coef.bma}})
#' @param ksubset integer vector detailing for which model sizes the plot
#' should be done
#' @param include.legend if \code{TRUE}, a small legend is included via the
#' low-level command \code{\link{legend}}
#' @param do.grid if \code{TRUE}, a \code{\link{grid}} is added to the plot
#' (with a simple \code{grid()}).
#' @param \dots parameters passed on to \code{\link{matplot}} with sensible
#' defaults
#' @return As a default, \code{plotModelsize} plots the posterior model size
#' distribution as a blue line, and the prior model distribution as a dashed
#' red line.\cr In addition, it returns a list with the following elements:
#' \item{mean}{The posterior expected value of model size} \item{var}{The
#' variance of the posterior model size distribution} \item{dens}{A vector
#' detailing the posterior model size distribution from model size \eqn{0} (the
#' first element) to \eqn{K} (the last element)}
#' 
#' @seealso See also \code{\link{bms}}, \code{\link{image.bma}},
#' \code{\link{density.bma}}, \code{\link{plotConv}}
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords hplot
#' @examples
#' 
#' data(datafls)
#' mm=bms(datafls,burn=1500, iter=5000, nmodel=200,mprior="fixed",mprior.size=6)
#' 
#' #plot Nb.1 based on aggregate results
#' postdist= plotModelsize(mm)
#' 
#' #plot based only on 30 best models
#' plotModelsize(mm[1:30],exact=TRUE,include.legend=FALSE)
#' 
#' #plot based on all best models, but showing distribution only for model sizes 1 to 20
#' plotModelsize(mm,exact=TRUE,ksubset=1:20)
#' 
#' # create a plot similar to plot Nb. 1
#' plot(postdist$dens,type="l") 
#' lines(mm$mprior.info$mp.Kdist)
#' 
#' 
#' @export
plotModelsize<-function(bmao,exact=FALSE,ksubset=NULL,include.legend=TRUE, do.grid=TRUE, ...) { #,lwd=1.5,xaxt="n",col=c("steelblue3","tomato"),main=NULL,cex.main=0.8,xlab="Model Size",ylab=""){
    #plots posterior vs. prior model size distribution
    # bma: bma object
    # exact: whether posterior distribution is to be based on MC3 results or on best 'top' models (as in bmao$topmod)
    # ksubset: an integer vector detailing for which model sizes the plot should be done
    # include.legend: if TRUE, a small legend is included via legend()
    # ... parameters passed on to ?matplot (Note: defaults are as above, default 'main' is 'Posterior Model Size Distribution Mean: x' where x is the actual mean)   
    # besides plotting, the function returns a list with the mean and variance of the posterior model size distribution (not altered by ksubset)
        
    dotargs = match.call(expand.dots=FALSE)$...
    
    if (length(exact)>1) { topmodidx=exact; exact=TRUE} else {topmodidx=NA}
    # get required information
    K=bmao$info$K
    if (is.element("mprior.info",names(bmao))) m=bmao$mprior.info$mp.msize else m=bmao$arguments$prior.msize        
    pmp.10=pmp.bma(bmao$topmod[topmodidx],oldstyle=TRUE)

    if(exact){
      modelSmean=sum(apply(.post.topmod.bma(bmao$topmod[topmodidx]),2,function(x) length(which(x==1)))*pmp.10[,1])
      modelS.var=sum(apply(.post.topmod.bma(bmao$topmod[topmodidx]),2,function(x) length(which(x==1)))^2*pmp.10[,1])-modelSmean^2
      x = apply(.post.topmod.bma(bmao$topmod[topmodidx]),2,function(x) length(which(x==1)))
      y = pmp.10[,1]
      result = c()
      for( i in sort(unique(x)) )
         result = c(result, sum(y[which(x==i)]))
      names(result) = sort(unique(x))
      kvec=rep(0,(K+1))
      kvec[(as.numeric(names(result))+1)]=result
    }
    else{
      k.vec=bmao$info$k.vec
     #calculate expected value
      summi=sum(k.vec)
     #modelSmean -1 because we have to deduct the constant from the mean model size
      modelSmean=sum((1:length(k.vec))*(k.vec/summi))-1
      kvec=k.vec/sum(k.vec)
      modelSmean.sq=sum(((1:length(k.vec))^2)*(k.vec/summi))
      modelS.var=modelSmean.sq-modelSmean^2   #var(x)=E(X^2)-(E(X))^2
    }

    upper=min(ceiling(modelSmean+5*modelS.var), K)
    lower=max(floor(modelSmean-5*modelS.var),0)

    
    if (is.element("mp.Kdist",names(bmao$mprior.info))) { 
      prior = bmao$mprior.info$mp.Kdist
    } else if (is.element("theta", names(bmao$arguments))) { #kept for historical reasons
       theta=bmao$arguments$theta
       if(theta=="random"){
            beta.bin=function(a=1,b=(K-m)/m,K=K,w=0:K){
              return(lgamma(a+b)-(lgamma(a)+lgamma(b)+lgamma(a+b+K))+log(choose(K,w))+lgamma(a+w)+lgamma(b+K-w))
            }
          prior=exp(beta.bin(a=1,b=(K-m)/m,K=K,w=0:K))
       }
       if(theta!="random"){ prior=stats::dbinom(x=0:K,size=K,prob=m/K,log=FALSE)}
    } else {
      prior = rep(NA,length(kvec)) 
    }
     
    mat=cbind(kvec,prior)
    
    #do the plot
        upper.ylim=max(kvec,prior, na.rm=TRUE)
        if(is.null(ksubset)){ksubset=(lower:upper)}
         dotargs=.adjustdots(dotargs,type="l",ylim=c(0,1.1*upper.ylim),lwd=1.5,xaxt="n",col=c("steelblue3","tomato"),main=paste("Posterior Model Size Distribution","\n","Mean:",round(modelSmean,4)),cex.main=0.8,xlab="Model Size",ylab="", lty=1:2, pch=4, cex.axis=.9)


        matsubset=mat[ksubset+1,]
        eval(as.call(c(list(as.name("matplot"),as.name("matsubset")),as.list(dotargs))))
        #matplot(mat[ksubset,],type="l",xaxt=xaxt,col=col,
        #        main=ifelse(is.null(main),paste("Posterior Model Size Distribution","\n","Mean:",
        #        round(modelSmean,4)),main),cex.main=cex.main,xlab=xlab,
        #        ylim=c(0,1.1*upper.ylim),lwd=lwd,ylab=ylab,...)
        if (as.logical(do.grid)) graphics::grid()
        graphics::points(kvec[ksubset+1],cex=0.8,pch=eval(dotargs$pch))
        
        #if(lower==0){
        #    graphics::axis(1, las=1, at =1:length(lower:upper), label = c(0:K)[lower:(upper+1)],cex.axis=0.7)
        #}
        #else{
        #    axis(1, las=1, at =1:length(lower:upper), label = c(0:K)[lower:upper],cex.axis=0.7)
        #}
        graphics::axis(1, las=1, at =1:length(ksubset), labels = ksubset, cex.axis=eval(dotargs$cex.axis))
        if (include.legend) {
          if (is.null(prior)||all(is.na(prior))) {
            graphics::legend(x="topright",lty=eval(dotargs$lty),legend=c("Posterior"),col=eval(dotargs$col),ncol=1,bty="n",lwd=eval(dotargs$lwd))
          } else {
            graphics::legend(x="topright",lty=eval(dotargs$lty),legend=c("Posterior","Prior"),col=eval(dotargs$col),ncol=2,bty="n",lwd=eval(dotargs$lwd))
          }
        }
   return(invisible(list(mean=modelSmean,var=modelS.var,dens=kvec)))
}


#' Coefficient Marginal Posterior Densities
#' 
#' Calculates the mixture marginal posterior densities for the coefficients
#' from a BMA object and plots them
#' 
#' The argument \code{addons} specifies what additional information should be
#' added to the plot(s) via the low-level commands \code{\link{lines}} and
#' \code{\link{legend}}:\cr \code{"e"} for the posterior expected value (EV) of
#' coefficients conditional on inclusion (see argument \code{exact=TRUE} in
#' \code{\link{coef.bma}}),\cr \code{"s"} for 2 times posterior standard
#' deviation (SD) bounds,\cr \code{"m"} for the posterior median,\cr \code{"b"}
#' for posterior expected values of the individual models whom the density is
#' averaged over,\cr \code{"E"} for posterior EV under MCMC frequencies (see
#' argument \code{exact=FALSE} in \code{\link{coef.bma}}),\cr \code{"S"} for
#' the corresponding SD bounds (MCMC),\cr \code{"p"} for plotting the Posterior
#' Inclusion Probability above the density plot,\cr \code{"l"} for including a
#' \code{\link{legend}}, \code{"z"} for a zero line, \code{"g"} for adding a
#' \code{\link{grid}}
#' 
#' Any combination of these letters will give the desired result. Use
#' \code{addons=""} for not using any of these.\cr In case of
#' \code{density.zlm}, only the letters \code{e}, \code{s}, \code{l}, \code{z},
#' and \code{g} will have an effect.
#' 
#' @aliases density.bma density.zlm
#' @param x A bma object (see \code{\link{bms}}) or a \code{\link{zlm}} object.
#' @param reg A scalar integer or character detailing which covariate's
#' coefficient should be plotted. If \code{reg=NULL} (default), then all
#' regressors are plotted one after the other, waiting for user interaction.
#' @param addons character. Specifies which additional information should be
#' added to the plot via low-level commands (see 'Details' below).
#' @param std.coefs logical. If \code{TRUE} then the posterior density is
#' estimated for standardized coefficients (representing the case where all
#' variables have mean zero and standard deviation 1) - default is
#' \code{FALSE}.
#' @param n numeric. the number of equally spaced points at which the density
#' is to be estimated.
#' @param plot logical.  If \code{TRUE} (default), the density is plotted; if
#' \code{FALSE} then \code{density.bma} only returns the estimated posterior
#' densities without plotting.
#' @param hnbsteps even integer, default 30. The number of numerical
#' integration steps to be used in case of a hyper-g prior (cf. argument
#' \code{g} in \code{\link{bms}}). Increase this number to increase accuracy.
#' @param addons.lwd scalar, default 1.5. Line width to be used for the
#' low-level plotting commands specified by \code{addons}. Cf. argument
#' \code{lwd} in \code{\link{par}}
#' @param \dots Additional arguments for \code{\link{plot.default}} with
#' sensible defaults
#' @return The function returns a list containing objects of the class
#' \code{\link{density}} detailing the marginal posterior densities for each
#' coefficient provided in \code{reg}.\cr In case of \code{density.zlm}, simple
#' marginal posterior coefficient densities are computed, while
#' \code{density.bma} calculates there mixtures over models according to
#' posterior model probabilities.\cr These densities contain only the density
#' points apart from the origin. (see 'Note' below)
#' 
#' As long as \code{plot=TRUE}, the densities are plotted too.  Note that (for
#' \code{density.bma}) if the posterior inclusion probability of a covariate is
#' zero, then it will not be plotted, and the returned density will be
#' \code{list(x=numeric(n),y=numeric(n))}.
#' @note The computed marginal posterior densities from \code{density.bma} are
#' a Bayesian Model Averaging mixture of the marginal posterior densities of
#' the individual models.  The accuracy of the result therefore depends on the
#' number of 'best' models contained in \code{x} (cf. argument \code{nmodel} in
#' \code{\link{bms}}).
#' 
#' The marginal posterior density can be interpreted as 'conditional on
#' inclusion': If the posterior inclusion probability of a variable is smaller
#' than one, then some of its posterior density is Dirac at zero.  Therefore
#' the integral of the returned density vector adds up to the posterior
#' inclusion probability, i.e. the probability that the coefficient is not
#' zero.
#' 
#' Correspondingly, the posterior EV and SD specified by \code{addons="es"} are
#' based on 'best' model likelihoods ('exact') and are conditional on
#' inclusion.  They correspond to the results from command
#' \code{coef.bma(x,exact=TRUE,condi.coef=TRUE,order.by.pip=FALSE)} (cf. the
#' example below).
#' 
#' The low-level commands enacted by the argument \code{addons} rely on colors
#' of the \code{\link{palette}}: color 2 for \code{"e"} and \code{"s"}, color 3
#' for \code{"m"}, color 8 for \code{"b"}, color 4 for \code{"E"} and
#' \code{"S"}. The default colors may be changed by a call to
#' \code{\link{palette}}.
#' 
#' Up to BMS version 0.3.0, \code{density.bma} may only cope with built-in
#' \code{gprior}s, not with any user-defined priors.
#' 
#' @seealso \code{\link{quantile.coef.density}} for extracting quantiles,
#' \code{\link{coef.bma}} for similar concepts, \code{\link{bms}} for creating
#' bma objects
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords aplot utilities
#' @examples
#' 
#' 
#'  data(datafls)
#'  mm=bms(datafls)
#' 
#'  density(mm,reg="SubSahara")
#'  density(mm,reg=7,addons="lbz") 
#'  density(mm,1:9)
#'  density(mm,reg=2,addons="zgSE",addons.lwd=2,std.coefs=TRUE)
#' 
#' # plot the posterior density only for the very best model
#'  density(mm[1],reg=1,addons="esz")
#' 
#' 
#' #using the calculated density for other purposes...
#'  dd=density(mm,reg="SubSahara")
#'  plot(dd) 
#' 
#'  dd_list=density(mm,reg=1:3,plot=FALSE,n=400)
#'  plot(dd_list[[1]])
#' 
#' 
#' #Note that the shown density is only the part that is not zero
#'  dd=density(mm,reg="Abslat",addons="esl")
#'  pip_Abslat=sum(dd$y)*diff(dd$x)[1]
#' 
#'  #this pip and the EV conform to what is done by the follwing command
#'  coef(mm,exact=TRUE,condi.coef=TRUE)["Abslat",]
#' 
#' @export
density.bma <- function(x,reg=NULL,addons="lemsz",std.coefs=FALSE,n=300,plot=TRUE,hnbsteps=30,addons.lwd=1.5,...) {
# x: a bma object
# reg: the covariate to be calculated; can be character or integer
# addons: low-level additions to the plot: "e" for post. expected value, "s" for 2 times Std.Dev. bounds,
#  "l" for including a legend, "m" for the median, "z" for a zero line
#  "b" for posterior exp. values of the individual models   
#  "E" for post. exp. value under MCMC frequencies, "S" for the corresponding SD bounds (MCMC)
#  ?"g" for adding a grid()
#   "p" for an additional box drawing the PIP
# n: the number of points for which the density should be estimated
# std.coefs: If TRUE then plot post. dist. of standardized coefficients
# plot: if FALSE, then only the density is returned; 
# hnbsteps: steps for numerical integration in case of a hyper-g prior
# addons.lwd: lty for addons lines
#...: parameters passed on to plot
#
#returns a list of class "density" (cf. ?density)

  dtcm=function(x,df,ncp,varp) {
     #a wrapper for univariate t-dist with non-centrality parameter and variance parameter 
     # (variance is df/(df-2)*varp)
     sqvarp=sqrt(varp)
     stats::dt((x-ncp)/sqvarp,df=df)/sqvarp
  }

  dsgivenykernel <- function(sf,kpa,N,z) {
    #the post. density of the shrinkge factor f(s|Y)*F((N-1)/2,1,(k+a)/2,R2)
    (kpa-2)/2*(1-sf)^((kpa-4)/2)*(1-sf*z)^(-(N-1)/2)
  }


  
    #check user input and get basic info
    dotargs = match.call(expand.dots=FALSE)$...
    bmao=x
    if (!is.bma(bmao)) stop("Argument bmao needs to be a bma object")
    if (hnbsteps%%2) stop("Argument nbsteps needs to be an even integer")
    nbsteps=max(hnbsteps,2)
    n=max(ceiling(n),1)
  
    N=bmao$info$N; K=bmao$info$K
    if(is.null(reg)) reg=1:K
    nameix=1:K; names(nameix)=bmao$reg.names; reg=nameix[reg]
    #if (is.na(reg)) stop("Argument reg is out of bounds")
    ishyper=(bmao$gprior$gtype=="hyper")
    tm=bmao$topmod
    bools=(tm$bool_binary())
    betas=tm$betas()
    betas2=tm$betas2()    
    
    
           
    if (std.coefs) {
      #if standardized coefficients are wanted, then adjust moments accordingly
      sddata=apply(as.matrix(bmao$arguments$X.data),2,stats::sd)
      betas=diag(sddata[-1])%*%betas/sddata[1]
      betas2=diag(sddata[-1]^2)%*%betas2/sddata[1]^2
    }
        
    sigmadiag=(betas2-betas^2)*(N-3)/(N-1) #the variance parameters for the t-dist
    
    pmps=pmp.bma(bmao$topmod,oldstyle=TRUE)[,1] #calc post model probs (exact)

    #the stuff below is similar to estimates.bma(bmao, condi.coef=TRUE,exact=TRUE)
    pips=c(tcrossprod(bools,t(pmps)))   
    Eb1=c(tcrossprod(betas,t(pmps)))/pips
    Ebsd=sqrt(c(tcrossprod(betas2,t(pmps)))/pips-Eb1^2); Ebsd[is.nan(Ebsd)]=0; Eb1[is.nan(Eb1)]=0
    Eball=cbind(Eb1,Ebsd) #conditional coefficients
    
    
    
    if ((any(grep("E",addons,ignore.case=FALSE)))|(any(grep("S",addons,ignore.case=FALSE)))) {
      #in case the user wants to include MCMC results (see estimates.bma(,exact=FALSE),
      #then compute them here
      Eb1.mcmc = bmao$info$b1mo/bmao$info$inccount
      Ebsd.mcmc = sqrt(bmao$info$b2mo/bmao$info$inccount-Eb1.mcmc^2)
      if (std.coefs) {
         sddata=apply(as.matrix(bmao$arguments$X.data),2,stats::sd)
         Eb1.mcmc=Eb1.mcmc*sddata[-1]/sddata[1];
         Ebsd.mcmc=Ebsd.mcmc*sddata[-1]/sddata[1];
      }      
    }
    

    if (ishyper) {
       #in case of hyper-g, we cannot rely on sigmadiag, but have to numerically integrate over
       #differnt shrinkages, for that we need this
       yXdata=as.matrix(bmao$arguments$X.data); yXdata=yXdata-matrix(colMeans(yXdata),N,K+1,byrow=TRUE)
       if (std.coefs) yXdata=yXdata%*%diag(1/sddata)
       yty=c(crossprod(yXdata[,1]))
       positions=lapply(lapply(as.list(as.data.frame(bools)),as.logical),which)
       olsmodels=lapply(lapply(positions,.ols.terms2,yty=yty,N=N,K=K,XtX.big=crossprod(yXdata[,-1]),Xty.big=c(crossprod(yXdata[,-1],yXdata[,1]))),function (x) x$full.results())       
       f21a=bmao$gprior.info$hyper.parameter
    }


 plotndens <- function(ix,doplot=FALSE) {
      #this function does marginal density and plot for a specific covariate
      #depends heavily on parent scope!
 
      sss=function(lbound,uboundp1,nbsteps) {
       #simple simpson integration over shrinkage factor
       #this is just a convenience function and depends on variables in parent scope, only used in case of hyper-g
       s.seq=seq(lbound,uboundp1,(uboundp1-lbound)/nbsteps)[-nbsteps]
       tmat=sapply(as.list(s.seq),function(ss) { dtcm(seqs,N-1,ss*bhati,invdiagi*ss*(1-ss*z)/(N-1)*yty)}) #matrix of t-densities for different s
       smat=sapply(as.list(s.seq),dsgivenykernel, kpa=k+f21a,N=N,z=z)  #vector of posterior densities for the differnet s
       if (any(is.infinite(smat))) smat[is.infinite(smat)]=0
       intconst=(4*sum(smat[c(FALSE,TRUE)])+2*sum(smat[c(TRUE,FALSE)])-3*smat[nbsteps]-smat[1])*(s.seq[nbsteps]-s.seq[1])/nbsteps/3 #calc the value of F((N-1)/2,1,(k+a)/2,R2)
       return(list(dv=c(4*tmat[,c(FALSE,TRUE)]%*%smat[c(FALSE,TRUE)]+2*tmat[,c(TRUE,FALSE)]%*%smat[c(TRUE,FALSE)]-3*tmat[,nbsteps]*smat[nbsteps]-tmat[,1]*smat[1])*(s.seq[nbsteps]-s.seq[1])/nbsteps/3, ic=intconst))
       #return the estimated density and a normalization constant
       #is done in parts because it may be recomposed later
     }

 

    if (pips[ix]==0) { 
        reslist=list(x=numeric(n),y=numeric(n),n=n,call=sys.call(),data.name=names(nameix)[ix],has.na=FALSE); class(reslist)=c("density", "coef.density")  
        return(reslist)
    }
    #the x vector
    lbound=min(betas[ix,as.logical(bools[ix,])])-3*Eball[ix,2]; ubound=max(betas[ix,as.logical(bools[ix,])])+3*Eball[ix,2]; 
    seqs=seq(lbound,ubound,(ubound-lbound)/(n-1))    


    densvec=numeric(length(seqs))
    
    #loop through models and calc post dens
    for (m in 1:length(pmps)) {
       if (bools[ix,m]) {
        if (ishyper) {
            ixadj=sum(bools[1:ix,m])
            bhati=olsmodels[[m]]$bhat[[ixadj]]; invdiagi=olsmodels[[m]]$diag.inverse[[ixadj]]; k=sum(bools[,m]);
            Esf=betas[ix,m]/bhati; z=1-olsmodels[[m]]$ymy/yty;  midpoint=1-(1-Esf)*4
          if (midpoint<.5) {
            dvl=sss(.0001,.9999999,nbsteps*2); addvec=dvl$dv/dvl$ic
          } else {
            dvl1=sss(.0001,midpoint,nbsteps); dvl2=sss(midpoint,1,nbsteps)
            addvec=(dvl1$dv+dvl2$dv)/(dvl1$ic+dvl2$ic)
          }
        } else {
          addvec=dtcm(seqs,N-1,betas[ix,m],sigmadiag[ix,m])
        }
         densvec=densvec+pmps[m]*addvec
     
       }
    }


    reslist=list(x=seqs,y=densvec,bw=NULL,n=n,call=sys.call(),data.name=names(nameix)[ix],has.na=FALSE); class(reslist)="density"


  if (!doplot) {  return(reslist) }
     #plot stuff  
     main_default=paste("Marginal Density:",names(nameix)[ix],"(PIP",round(c(crossprod(pmps,bools[ix,]))*100,2),"%)")
     if (any(grep("p",addons,ignore.case=TRUE))) { 
        decr=.12; 
        parplt=graphics::par()$plt; 
        parplt_temp=parplt; parplt_temp[4]=(1-decr)*parplt[4] +decr*parplt[3]; graphics::par(plt=parplt_temp)
        main_temp=main_default; main_default=NULL
#         graphics::layout(1:2,heights=c(.1,1))
#         opm=par()$mar
#         par(mar=c(0,opm[2],1,opm[4]))
#         plot(0,type="n",xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
#         mtext("PIP",side=2,las=2)
#         rect(0,0,1,1,col=8)
#         rect(0,0,pips[ix],1,col=9)
#         par(mar=opm)
     }
     

    
    dotargs=.adjustdots(dotargs,type="l", col="steelblue4", main=main_default,
        xlab=if (std.coefs) "Standardized Coefficient" else "Coefficient", ylab = "Density");
#    if (!is.element('type',names(dotargs))) { dotargs$type <- "l" }
#    if (!is.element('col',names(dotargs))) { dotargs$col <- "steelblue4" }
#    if (!is.element('main',names(dotargs))) { dotargs$main=paste("Marginal Density:",names(nameix)[ix],"(PIP",round(c(crossprod(pmps,bools[ix,]))*100,2),"%)") }
#    if (!is.element('xlab',names(dotargs))) { dotargs$xlab <- if (std.coefs) "Standardized Coefficient" else "Coefficient" }
#    if (!is.element('ylab',names(dotargs))) { dotargs$ylab <- "Density" }
#

    eval(as.call(c(list(as.name("plot"),x=as.name("seqs"),y=as.name("densvec")),as.list(dotargs))))
    leg.col=numeric(0);leg.lty=numeric(0); leg.legend=character(0)
    
    if (any(grep("g",addons,ignore.case=TRUE))) { # grid
       graphics::grid()
    }

    if (any(grep("b",addons,ignore.case=TRUE))) { # post exp values of the individual models
      for (m in 1:length(pmps)) {
         Ebm=betas[ix,m] ;
         if (as.logical(Ebm)) {
           Ebheight=min(densvec[max(sum(seqs<Ebm),1)],densvec[sum(seqs<Ebm)+1])
           graphics::lines(x=rep(Ebm,2),y=c(0,Ebheight),col=8)
         }
      }
      leg.col=c(leg.col,8);leg.lty=c(leg.lty,1);leg.legend=c(leg.legend,"EV Models")
    }

    if (any(grep("e",addons,ignore.case=FALSE))) {  # posterior mean
       graphics::abline(v=Eball[ix,1],col=2,lwd=addons.lwd)
       leg.col=c(leg.col,2);leg.lty=c(leg.lty,1); leg.legend=c(leg.legend,"Cond. EV")
    }
    if (any(grep("s",addons,ignore.case=FALSE))) { # posterior SD bounds
       graphics::abline(v=Eball[ix,1]-2*Eball[ix,2],col=2,lty=2,lwd=addons.lwd)
       graphics::abline(v=Eball[ix,1]+2*Eball[ix,2],col=2,lty=2,lwd=addons.lwd)    
       leg.col=c(leg.col,2);leg.lty=c(leg.lty,2);leg.legend=c(leg.legend,"2x Cond. SD")       
    }
    if (any(grep("m",addons,ignore.case=TRUE))) { # posterior median
       median_index=sum(cumsum(densvec)<sum(densvec)/2)
       graphics::abline(v=(seqs[median_index]+seqs[median_index+1])/2,col=3,lwd=addons.lwd)
       leg.col=c(leg.col,3);leg.lty=c(leg.lty,1);leg.legend=c(leg.legend,"Median")       
    }
    
    if (any(grep("z",addons,ignore.case=TRUE))) { #zero line
      graphics::abline(h=0,col="gray",lwd=addons.lwd)
    }
    
    if (any(grep("E",addons,ignore.case=FALSE))) { #post exp value of MCMC results (see estimates.bma(,exact=F)
      graphics::abline(v=Eb1.mcmc[ix],col=4,lwd=addons.lwd)
       leg.col=c(leg.col,4);leg.lty=c(leg.lty,1); leg.legend=c(leg.legend,"Cond. EV (MCMC)")
    }
    if (any(grep("S",addons,ignore.case=FALSE))) { #2 times post SD of MCMC results (see estimates.bma(,exact=F)
      graphics::abline(v=Eb1.mcmc[ix]-2*Ebsd.mcmc[ix],col=4,lty=2,lwd=addons.lwd)
      graphics::abline(v=Eb1.mcmc[ix]+2*Ebsd.mcmc[ix],col=4,lty=2,lwd=addons.lwd)       
       leg.col=c(leg.col,4);leg.lty=c(leg.lty,2); leg.legend=c(leg.legend,"2x SD (MCMC)")
    }
    
    
    if (any(grep("l",addons,ignore.case=TRUE))&(length(leg.col)>0)) { #legend
      leg.pos="topright"; if (Eball[ix,1]>seqs[floor(n/2)]) leg.pos="topleft"; 
      graphics::legend(x=leg.pos,lty=leg.lty,col=leg.col,legend=leg.legend,box.lwd=0,bty="n",lwd=addons.lwd)
    }
    
    if (any(grep("p",addons,ignore.case=TRUE))) { 
      pusr=graphics::par()$usr
      graphics::rect(pusr[1],pusr[4]*(1+decr*.2), pusr[2], pusr[4]*(1+decr),xpd=TRUE,col=8)
      graphics::rect(pusr[1],pusr[4]*(1+decr*.2), pips[ix]*pusr[2]+(1-pips[ix])*pusr[1], pusr[4]*(1+decr),xpd=TRUE,col=9)
      graphics::mtext("PIP:",side=2, las=2,line=1, at=pusr[4]*(1+decr*.6))
      graphics::par(plt=parplt)
      graphics::title(main_temp)
    }
    return(reslist)
  }
  
  
  densres=list()
  oldask=graphics::par()$ask
  plots=0
  for (vbl in 1:length(reg)) {
    doplot=(if (as.logical(pips[reg[vbl]])) plot else FALSE)
    plots=plots+doplot
    if (plots==2) {graphics::par(ask=TRUE)}
    densres[[nameix[vbl]]]=plotndens(reg[vbl],doplot)
    densres[[nameix[vbl]]]$call=sys.call()   #call("density.bma",bmao=bmao,reg=reg,n=300,hnbsteps=30)
  }
  
  graphics::par(ask=oldask)
  if (length(densres)==1) densres=densres[[1]] else class(densres) = c("coef.density",class(densres))
  if (!plot) return(densres)
  if (plot&(plots==0)) {warning("No plot produced as PIPs of provided variables are zero under 'exact' estimation.")}
  return(invisible(densres))
}


#' Posterior Density of the Shrinkage Factor
#' 
#' Calculates the mixture marginal posterior density for the shrinkage factor
#' (g/(1+g)) from a BMA object under the hyper-g prior and plots it
#' 
#' The function \code{gdensity} estimates and plots the posterior density for
#' the shrinkage factor \eqn{g/(1+g)}\cr This is evidently only possible if the
#' shrinkage factor if not fixed, i.e. if the bma object \code{x} was estimated
#' with a hyper-g prior - cf. argument \code{g} in \code{\link{bms}}\cr The
#' density is based only on the best models retained in the bma object
#' \code{x}, cf. argument \code{nmodel} in \code{\link{bms}}\cr A note on
#' argument \code{n}: The points at which the density is estimated start at
#' \eqn{max(0,E-5*SD)}, where \eqn{E} and \eqn{SD} are the expected value and
#' standard deviation of the shrinkage factor, respectively. For plotting the
#' entire domain \eqn{(0,1)} use \code{xlim=c(0,1)} as an argument for
#' \code{gdensity}.
#' 
#' The argument \code{addons} specifies what additional information should be
#' added to the plot(s) via the low-level commands \code{\link{lines}} and
#' \code{\link{legend}}:\cr \code{"e"} for the posterior expected value (EV) of
#' the shrinkage factor,\cr \code{"s"} for 2 times posterior standard deviation
#' (SD) bounds,\cr \code{"m"} for the posterior median,\cr \code{"f"} for
#' posterior expected values of the individual models whom the density is
#' averaged over,\cr \code{"z"} for a zero line, \code{"l"} for including a
#' \code{\link{legend}}\cr The following two are only possible if the bma
#' object collected statistics on shrinkage, cf. argument \code{g.stats} in
#' \code{\link{bms}} \code{"E"} for posterior expected value under MCMC
#' frequencies (see argument \code{exact} in \code{\link{coef.bma}}),\cr
#' \code{"S"} for the corresponding 2 times standard deviation bounds
#' (MCMC),\cr
#' 
#' Any combination of these letters will give the desired result. Use
#' \code{addons=""} for not using any of these.
#' 
#' @param x A bma object (see \code{\link{bms}}).
#' @param n The integer number of equally spaced points at which the density is
#' to be estimated. see 'Details' below
#' @param addons character, defaulting to \code{"zles"}. Specifies which
#' additional information should be added to the plot via low-level commands
#' (see 'Details' below).
#' @param plot logical.  If \code{TRUE} (default), the density is plotted; if
#' \code{FALSE} then \code{gdensity} only returns the estimated posterior
#' density without plotting.
#' @param addons.lwd scalar, default 1.5. Line width to be used for the
#' low-level plotting commands specified by \code{addons}. Cf. argument
#' \code{lwd} in \code{\link{par}}
#' @param \dots Additional arguments for \code{\link{plot.default}} with
#' sensible defaults
#' @return \code{gdensity} returns an object of the class \code{\link{density}}
#' detailing the posterior mixture density of the shrinkage factor.
#' @note The computed marginal posterior density is a Bayesian Model Averaging
#' mixture of the marginal posterior densities of the shrinkage factor under
#' individual models.  The accuracy of the result therefore depends on the
#' number of 'best' models contained in \code{x} (cf. argument \code{nmodel} in
#' \code{\link{bms}}).
#' 
#' Correspondingly, the posterior EV and SD specified by \code{addons="es"} are
#' based on 'best' model likelihoods ('exact') and are conditional on
#' inclusion.
#' 
#' The low-level commands enacted by the argument \code{addons} rely on colors
#' of the \code{\link{palette}}: color 2 for \code{"e"} and \code{"s"}, color 3
#' for \code{"m"}, color 8 for \code{"f"}, color 4 for \code{"E"} and
#' \code{"S"}. The default colors may be changed by a call to
#' \code{\link{palette}}.
#' 
#' @seealso \code{\link{density.bma}} for computing coefficient densities,
#' \code{\link{bms}} for creating bma objects, \code{\link{density}} for the
#' general method
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords aplot utilities
#' @examples
#' 
#' 
#'  data(datafls)
#'  mm=bms(datafls,g="hyper=UIP")
#' 
#'  gdensity(mm) # default plotting
#'  
#'  # the grey bars represent expected shrinkage factors of the individual models
#'  gdensity(mm,addons="lzfes") 
#'  
#'  # #plotting the median 'm' and the posterior mean and bounds based on MCMC results:
#'  gdensity(mm,addons="zSEm",addons.lwd=2)
#' 
#' # plot the posterior shrinkage density only for the very best model
#'  gdensity(mm[1],addons="esz")
#' 
#' 
#' #using the calculated density for other purposes...
#'  dd=gdensity(mm,plot=FALSE)
#'  plot(dd) 
#' 
#' @export
gdensity <- function(x,n=512,plot=TRUE,addons="zles",addons.lwd=1.5,...) { #main="Posterior Density of the Shrinkage Factor",type="l",xlab="Shrinkage factor",ylab="Density",col="steelblue4"
# plots posterior density of shrinkage factor for hyper-g bma objects
# INPUTS:
#  x: bma object
#  n: number of equally spaced points where density is computed
#  plot: whether to do a plot in addition or not
#  addons: low-level plot commands: draw "z"=zero line, "l"=legend, "e"=exact exp. val, "s"=2x St dev bounds (exact), "E"=MCMC exp. val., "S"=MCMC SD bounds, "g"=shrinkage Exp vals for individual models
#  addons.lwd: linwe width for addons stuff
#  ... commands pased on to plot.default

  dsgivenykernel <- function(kpazvec,sf,N) {
    #the post. density of the shrinkge factor f(s|Y)*F((N-1)/2,1,(k+a)/2,R2)
    #kpazvec is a vector with two elements: first element is k+a, second is z (the R-squared)
    #sf: a vector of shrinkage values
    #N: sample size
    (kpazvec[[1]]-2)/2*(1-sf)^((kpazvec[[1]]-4)/2)*(1-sf*kpazvec[[2]])^(-(N-1)/2)
    #mode is at sf=(N-k-a+3)/(N-1-z*(k+a-4))
  }

    #user checks
     if (!is.bma(x)) stop("argument needs to an object of class 'bma'")
     if (!(x$gprior$gtype=="hyper")) stop("g prior density makes only sense for hyper-g prior.")
     if (n<2) stop("n needs to be at least 2")
     n=floor(n); #if (!(n%%2)) {n=n+1}
     dotargs = match.call(expand.dots=FALSE)$...
     
    #extract info
      N=x$info$N; K=x$info$K
      tm=x$topmod
      bools=tm$bool_binary()
      betas=tm$betas()
      betas2=tm$betas2()
      smoments = tm$fixed_vector()


    #re-compute the R-squareds and other stuff for the topmodels
       yXdata=as.matrix(x$arguments$X.data); yXdata=yXdata-matrix(colMeans(yXdata),N,K+1,byrow=TRUE)
       yty=c(crossprod(yXdata[,1]))
       positions=lapply(lapply(as.list(as.data.frame(bools)),as.logical),which) #vector of who is in where
       ymyvec=unlist(lapply(lapply(positions,.ols.terms2,yty=yty,N=N,K=K,XtX.big=crossprod(yXdata[,-1]),Xty.big=c(crossprod(yXdata[,-1],yXdata[,1]))),function (x) x$full.results()$ymy)) # vector of SSResid
       kvec=tm$kvec_raw() #vector of parameter number
       zvec=1-ymyvec/yty #vector of r.-squared
       pmpexact=pmp.bma(x,oldstyle=TRUE)[,1] #vector of 'exact' posterior model probs


       f21a=x$gprior.info$hyper.parameter #hyper parameter

       if (length(smoments)==0) { #if not there, re-compute individual moments for models
         #lprob=.lprob.hyperg.init(N=N,K=K,yty=yty,f21a=f21a, return.gmoments=TRUE)
         lprob = x$gprior.info$lprobcalc
         smoments = sapply(lapply(as.list(as.data.frame(rbind(kvec,ymyvec))),function(x) lprob$lprob.all(ymy=x[2],k=x[1],bhat=numeric(x[1]),diag.inverse=rep(1,x[1]))),"[[","otherstats")
       }

      Es=c(crossprod(smoments[1,],pmpexact)) #exp val
      Es2=c(crossprod(smoments[2,],pmpexact))
      Esd = sqrt(Es2-Es^2) #st dev



      nbsteps=n
      cutoff=max(0,Es-5*Esd) #this is to concentrate on where the mass is
      sdiff=(1-cutoff)/(nbsteps+1)
      s.seq=seq(sdiff+cutoff,cutoff+nbsteps*sdiff,sdiff)
      sdensl=lapply(as.list(as.data.frame(rbind(kvec+f21a,zvec))),dsgivenykernel,sf=s.seq,N=N)
      intconsts=lapply(lapply(sdensl,sum),"*",sdiff)       # a crude numerical integration to save time
      sdensvecs=mapply("/",sdensl,intconsts) #normalize by integration constants
      sdens=sdensvecs%*%pmpexact #mixture density

      reslist = list(x=s.seq,y=sdens,bw=NULL,n=n,call=sys.call(),data.name="Shrinkage",has.na=FALSE);
      class(reslist)="density"

      if (!plot) {return(reslist)}


   ##### PLOTTING #########################

    #Main plot
    dotargs=.adjustdots(dotargs,ylab="Density", xlab = "Shrinkage factor", main = "Posterior Density of the Shrinkage Factor", type="l", col="steelblue4")

    eval(as.call(c(list(as.name("plot"),as.name("s.seq"),as.name("sdens")),as.list(dotargs))))


    #preparing stuff for addons
      leg.col=numeric(0);leg.lty=numeric(0); leg.legend=character(0)


    if (any(grep("f",addons,ignore.case=TRUE))) { # post exp values of the individual models
      for (m in 1:length(pmpexact)) {
         Esm=smoments[1,m] ;
         if (as.logical(Esm)) {
           ixlower=max(sum(s.seq<Esm),1)
           Esheight=(sdens[ixlower+1]-sdens[ixlower])*(Esm-s.seq[ixlower])+sdens[ixlower]
           graphics::lines(x=rep(Esm,2),y=c(0,Esheight),col=8,lwd=addons.lwd)
         }
      }
      leg.col=c(leg.col,8);leg.lty=c(leg.lty,1);leg.legend=c(leg.legend,"EV Models")
    }



    if (any(grep("e",addons,ignore.case=FALSE))) {  # posterior mean
      graphics::abline(v=Es,col=2,lwd=addons.lwd)
       leg.col=c(leg.col,2);leg.lty=c(leg.lty,1); leg.legend=c(leg.legend,"EV")
    }
    if (any(grep("s",addons,ignore.case=FALSE))) { # posterior SD bounds
       if (!(Es-2*Esd)<0) graphics::abline(v=Es-2*Esd,col=2,lty=2,lwd=addons.lwd)
       if (!(Es+2*Esd)>1) graphics::abline(v=Es+2*Esd,col=2,lty=2,lwd=addons.lwd)
       leg.col=c(leg.col,2);leg.lty=c(leg.lty,2);leg.legend=c(leg.legend,"2x SD")
    }
    if (any(grep("m",addons,ignore.case=TRUE))) { # posterior median
       median_index=sum(cumsum(sdens)<sum(sdens)/2)
       graphics::abline(v=(s.seq[median_index]+s.seq[median_index+1])/2,col=3,lwd=addons.lwd)
       leg.col=c(leg.col,3);leg.lty=c(leg.lty,1);leg.legend=c(leg.legend,"Median")
    }

    if (any(grep("z",addons,ignore.case=TRUE))) { #zero line
      graphics::abline(h=0,col="gray",lwd=addons.lwd)
    }



    if (any(grep("E",addons,ignore.case=FALSE))) { #post exp value of MCMC results (see estimates.bma(,exact=F)
       if (all(x$gprior.info$shrinkage.moments==0)) warning("bma object needs to contain posterior g statistics - cf. argument 'g.stats' in 'help(bms)'") else {
         graphics::abline(v=x$gprior.info$shrinkage.moments[1],col=4,lwd=addons.lwd)
          leg.col=c(leg.col,4);leg.lty=c(leg.lty,1); leg.legend=c(leg.legend,"EV (MCMC)")
       }
    }
    if (any(grep("S",addons,ignore.case=FALSE))) { #2 times post SD of MCMC results (see estimates.bma(,exact=F)
      if (!all(x$gprior.info$shrinkage.moments==0)) {
         ES=x$gprior.info$shrinkage.moments[1]; SDs=sqrt(x$gprior.info$shrinkage.moments[2]-x$gprior.info$shrinkage.moments[1]^2)
         if (ES-2*SDs>0) graphics::abline(v=ES-2*SDs,col=4,lty=2,lwd=addons.lwd)
         if (ES+2*SDs<1) graphics::abline(v=ES+2*SDs,col=4,lty=2,lwd=addons.lwd)
         leg.col=c(leg.col,4);leg.lty=c(leg.lty,2); leg.legend=c(leg.legend,"2x SD (MCMC)")
       }
    }


    if (any(grep("l",addons,ignore.case=TRUE))&(length(leg.col)>0)) { #legend
      leg.pos="topleft";
      graphics::legend(x=leg.pos,lty=leg.lty,col=leg.col,legend=leg.legend,box.lwd=0,bty="n")
    }




    return(invisible(reslist))
}


#' @rdname quantile.pred.density
#' @export
quantile.density = function(x, probs=seq(.25,.75,.25), names=TRUE, normalize=TRUE, ...) {
# a generic function for objects of class density or lists whose elements are densities

  # the actual subfunction for object of class  "density"
  my.quantile.density = function(x, probs, names, normalize, ...) {
    ycs=(cumsum(x$y)-(x$y-x$y[[1]])/2)*diff(x$x[1:2])
    if (normalize) ycs=ycs/(ycs[[length(ycs)]])
    xin=x$x; maxi=length(ycs)
    qqs=sapply(as.list(probs), function(qu) {iii=sum(ycs<=qu); if (iii==maxi) return(Inf) else if (iii==0L) return(-Inf) else {  return(
      xin[[iii+1]] + ( (ycs[[iii+1]]-qu)/(ycs[[iii+1]]-ycs[[iii]]) ) *(xin[[iii]]-xin[[iii+1]]) )
      }})
    if (as.logical(names)) names(qqs)= paste(format(100 * probs, trim = TRUE, digits = max(2L, getOption("digits"))), "%", sep = "")#paste(signif(probs,5),"%",sep="")
    return(qqs)
  }  
  
  # user checks
  probs=as.vector(probs)
  if (is.element("density",class(x))) return(my.quantile.density(x=x, probs=probs, names=names, normalize=normalize))
  if (!all(sapply(x,function(dd) is.element("density",class(dd))))) stop("x needs to be a density or list of densities")
  if (length(x)==1L) return(my.quantile.density(x=x[[1]], probs=probs, names=names, normalize=normalize))
  #combining a list of denisties
  qout=sapply(x, my.quantile.density, probs=probs, names=FALSE, normalize=normalize) 
  
  #formatting the output into a matrix
  if (!is.matrix(qout)) { #some formatting
    if (length(probs)>1) return(qout)
    qout=as.matrix(qout)
  }  else qout=t(qout)
  if (as.logical(names)) colnames(qout)= paste(format(100 * probs, trim = TRUE, digits = max(2L, getOption("digits"))), "%", sep = "")
  return(qout)
}


.quantile.density=quantile.density

#' @rdname quantile.pred.density
#' @export
quantile.coef.density = function(x, probs=seq(.25,.75,.25), names=TRUE, ...) {
  #customizing quantile.density to stuff resulting from density.bma
 quout= .quantile.density(x, probs=probs, names=names, normalize=TRUE) 
 if (is.matrix(quout)&&as.logical(names)) rownames(quout) <- sapply(x, function(lx) lx[["data.name"]])
 return(quout)
}

#' Extract Quantiles from 'density' Objects
#' 
#' Quantiles for objects of class "density", "pred.density" or "coef.density"
#' 
#' The methods \code{quantile.coef.density} and \code{quantile.pred.density}
#' both apply \code{quantile.density} to densities nested with object of class
#' \code{coef.density} or \code{pred.density}.\cr The function
#' \code{quantile.density} applies generically to the built-in class
#' \code{\link{density}} (as least for versions where there is no such method
#' in the pre-configured packages).\cr Note that \code{quantile.density} relies
#' on trapezoidal integration in order to compute the cumulative densities
#' necessary for the calculation of quantiles.
#' 
#' @aliases quantile.pred.density quantile.coef.density quantile.density
#' @param x a object of class \code{\link{pred.density}}, \code{coef.density},
#' \code{\link{density}}, or a list of densities.
#' @param probs numeric vector of probabilities with values in [0,1] - elements
#' very close to the boundaries return \code{Inf} or \code{-Inf}
#' @param names logical; if \code{TRUE}, the result has a \code{names}
#' attribute, resp. a \code{rownames} and \code{colnames} attributes. Set to
#' \code{FALSE} for speedup with many probs.
#' @param normalize logical; if \code{TRUE} then the values in \code{x$y} are
#' multiplied with a factor such that their integral is equal to one.
#' @param \dots further arguments passed to or from other methods.
#' @return If \code{x} is of class \code{density} (or a list with exactly one
#' element), a vector with quantiles.\cr If \code{x} is a \code{\link{list}} of
#' densities with more than one element (e.g. as resulting from
#' \code{pred.density} or \code{coef.density}), then the output is a matrix of
#' quantiles, with each matrix row corresponding to the respective density.
#' @author Stefan Zeugner
#' @seealso \code{\link{quantile.default}} for a comparable function,
#' \code{\link{pred.density}} and \code{\link{density.bma}} for the
#' BMA-specific objects.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'  data(datafls)
#'  mm = bms(datafls[1:70,], user.int=FALSE)
#'  
#'  #predict last two observations with preceding 70 obs:
#'  pmm = pred.density(mm, newdata=datafls[71:72,], plot=FALSE) 
#'  #'standard error' quantiles
#'  quantile(pmm, c(.05, .95))
#'  
#'  #Posterior density for Coefficient of "GDP60"
#'  cmm = density(mm, reg="GDP60", plot=FALSE) 
#'  quantile(cmm, probs=c(.05, .95))
#'  
#'  
#'  #application to generic density:
#'  dd1 = density(rnorm(1000))
#'  quantile(dd1)
#'  
#'\dontrun{
#'  #application to list of densities:
#'  quantile.density( list(density(rnorm(1000)), density(rnorm(1000))) )
#' }
#' 
#' @export
quantile.pred.density = function(x, probs=seq(.25,.75,.25), names=TRUE, ...) {
  #customizing quantile.density to stuff resulting from pred.density
  quout= .quantile.density(x$densities(), probs=probs, names=names, normalize=FALSE)
   if (is.matrix(quout)&&as.logical(names)) rownames(quout) <- names(x$fit)
  return(quout)
}


#' Plot Convergence of BMA Sampler
#' 
#' Plots the posterior model probabilites based on 1) marginal likelihoods and
#' 2) MCMC frequencies for the best models in a 'bma' object and details the
#' sampler's convergence by their correlation
#' 
#' A call to bms with a MCMC sampler (e.g.
#' \code{bms(datafls,mcmc="bd",nmodel=100)} uses a Metropolis-Hastings
#' algorithm to sample through the model space: the frequency of how often
#' models are drawn converges to the distribution of their posterior marginal
#' likelihoods.\cr While sampling, each 'bma' object stores the best models
#' encountered by its sampling chain with their marginal likelihood and their
#' MCMC frequencies.\cr \code{plotConv} compares the MCMC frequencies to
#' marginal likelihoods, and thus visualizes how well the sampler has
#' converged.
#' 
#' @param bmao an object of class 'bma' - see \code{\link{bms}}
#' @param include.legend whether to include a \code{\link{legend}} in the plot
#' @param add.grid whether to include a \code{\link{grid}} in the plot
#' @param \dots other parameters for \code{\link{matplot}}
#' @note \code{plotConv} is also used by \code{\link{plot.bma}}
#' 
#' @seealso \code{\link{pmp.bma}} for posterior model probabilites based on the
#' two concepts, \code{\link{bms}} for creating objects of class 'bma'
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords aplot
#' @examples
#' 
#' 
#' data(datafls)
#' mm=bms(datafls[,1:12],user.int=FALSE)
#' 
#' plotConv(mm)
#' 
#' #is similar to
#' matplot(pmp.bma(mm),type="l")
#' 
#' @export
plotConv<-function(bmao,include.legend=TRUE,add.grid=TRUE,...){
 #function that compares MCMC and exact PMP's for the models in bmao$topmod
 # bmao: bma object
 # include.legend: TRUE or FALSE
 # main: if NULL then default is "Posterior Model Probabilites (Corr: x)", where x is the correlation betw. exact PMPs and MCMC freqs
 # other parameters: cf. help(plot.default)
 
  if (!is.bma(bmao)) stop("submit an object of class bma")
  # now get PMP analytical and MCMC
  mat=pmp.bma(bmao,oldstyle=TRUE)
  norm_const = sum(mat[,1])/sum(mat[,2])
  mat=cbind(mat[,2]*norm_const,mat[,1])
  if (length(bmao$topmod$lik())==0L) {
    stop("plotConv needs at least one model stored in topmod in order to produce a plot")
  }
  cor.pmp=format(round(.cor.topmod(bmao$topmod),4),nsmall=4)
  dotargs = match.call(graphics::plot,expand.dots=FALSE)$...
  dotargs=.adjustdots(dotargs, lwd=2,main=paste("Posterior Model Probabilities\n(Corr: ",cor.pmp,")",sep=""),lty=1,col=c("steelblue3","tomato"),cex.main=0.8,xlab="Index of Models",ylab="",type="l")
  eval(as.call(c(list(as.name("matplot"),as.name("mat")),as.list(dotargs))))

    
#   if (is.null(main)) main=paste("Posterior Model Probabilities\n(Corr: ",cor.pmp,")",sep="")
#  matplot(mat,type="l",lty=lty,col=col,lwd=lwd,main=main,xlab=xlab,cex.main=cex.main,ylab=ylab,...)
  if (as.logical(add.grid)) graphics::grid()
  if (as.logical(include.legend)) graphics::legend("topright",lty=eval(dotargs$lty),legend=c("PMP (MCMC)", "PMP (Exact)"),col=eval(dotargs$col),ncol=2,bty="n",cex=1,lwd=eval(dotargs$lwd));
}











#' Compare Two or More bma Objects
#' 
#' Plots a comparison of posterior inclusion probabilites, coefficients or
#' their standard deviation between various bma objects
#' 
#' 
#' @param \dots one or more objects of class 'bma' to be compared.
#' \code{plotComp} passes on any other parameters in \code{\dots{}} to
#' \code{\link{matplot}}.
#' @param varNr optionally, covariate indices to be included in the plot, can
#' be either integer vector or character vector - see examples
#' @param comp a character denoting what should be compared: \code{comp="PIP"}
#' (default) for posterior inclusion probabilities, \code{comp="Post Mean"} for
#' coefficients, \code{comp="Post SD"} for their standard deviations,
#' \code{comp="Std Mean"} or standardized coefficients, or \code{comp="Std SD"}
#' for standardized standard deviations
#' @param exact if \code{FALSE}, the statistics to be compared are based on
#' aggregate bma statistics, if \code{TRUE}, they are based solely on the best
#' models retained in the bma objects
#' @param include.legend whether to include a default legend in the plot
#' (custom legends can be added with the command \code{\link{legend}})
#' @param add.grid whether to add a \code{\link{grid}} to the plot
#' @param do.par whether to adjust \code{par("mar")} in order to fit in the
#' tick labels on the x-axis
#' @param cex.xaxis font size scaling parameter for the x-axis - cf. argument
#' \code{cex.axis} in \code{\link{par}}
#' 
#' @seealso \code{\link{coef.bma}} for the underlying function
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords hplot
#' @examples
#' 
#' ## sample two simple bma objects
#' data(datafls)
#' mm1=bms(datafls[,1:15])
#' mm2=bms(datafls[,1:15])
#' 
#' #compare PIPs
#' plotComp(mm1,mm2)
#' 
#' #compare standardized coefficeitns
#' plotComp(mm1,mm2,comp="Std Mean")
#' 
#' #...based on the lieklihoods of best models 
#' plotComp(mm1,mm2,comp="Std Mean",exact=TRUE)
#' 
#' #plot only PIPs for first four covariates
#' plotComp(mm1,mm2,varNr=1:4, col=c("black","red"))
#' 
#' #plot only coefficients for covariates 'GDP60 ' and 'LifeExp'
#' plotComp(mm1,mm2,varNr=c("GDP60", "LifeExp"),comp="Post Mean")
#' 
#' 
#' 
#' @export
plotComp <-function(...,varNr=NULL,comp="PIP",exact=FALSE,include.legend=TRUE,add.grid=TRUE,do.par=TRUE,cex.xaxis=0.8) { #,main=NULL,type="p",lty=1:5, lwd=1.5, pch=NULL,col=NULL,cex=NULL,bg=NA,xlab="",ylab=NULL){
# this plot compares results from different bma specifications, for example different W matrices
# for the SAR bma or if a new computer routine was coded we want to compare results of it with that
# of former computer routines
# in case you plug in results for different sets of regressors (e.g. once you use interaction sampling with 
# some interacted variables and compare that to results without interaction sampling) the plot looks only 
# at the coefficients / SD's for the set of variables that is contained in both specifications; so you can
# answer questions like, if I include this set of variables, how do they impact on remaining coefficients if
# the set is not included. 
# you can specify bmaList (list of bma results), varNr is the number of variables you want to see in the plot,
# comp is one of "PIP", "Post Mean", "Post SD" or "Std Mean" or "Std SD" (standardized coefficients),
# bmaNames are the names the routines get assigned in the legend (e.g. "new Model", "old Model" or "inv Dist.", "inv Dist Squared"),
# main is the title of the graph and cexM is the size of the labeling of the regressor names

  col_default =c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E" ,"#E6AB02", "#A6761D", "#666666")

    
    #bmaList = list(...)
  bmaList=list(...); bmaix=sapply(bmaList,is.bma); bmaList=bmaList[bmaix]
  
  #user check:
  bmaNr=length(bmaList)
  
  if(!all(sapply(bmaList,is.bma))){
    stop("Submit only bma objects to compare results (no other objects)")
  }
  dotargs = match.call(expand.dots=FALSE)$...;
  dotargs=.adjustdots(dotargs,ylab=paste(comp),pch=1:bmaNr,col=col_default,type="p",lty=1:5,lwd=1.5,xlab="",xaxt="n")
  dotargs=dotargs[!c(bmaix,logical(length(dotargs)-length(bmaix)))]

  
  #  # take care of name order and look whether we compare results for same set of variables!
  xMat=lapply(bmaList,function(x) rownames(estimates.bma(x,exact=exact)))
  xNames=xMat[[1]]
  ind=as.numeric(unlist(lapply(xMat,function(x) length(x))))

  # in case we do not have the same set of vars in each submitted specification
  if(length(unique(ind)>1)){
     smallestSet=which.min(ind)
     indMat=array(0:0,dim=c(length(xMat[[smallestSet]]),length(ind)))
     for(i in 1:length(ind)){
      indMat[,i]=as.numeric(xMat[[smallestSet]] %in% xMat[[i]])
     }
     xNamesInd=which(rowSums(indMat)==bmaNr)
     xNames=xMat[[smallestSet]][xNamesInd]
  }

  
  compNames=c(colnames(estimates.bma(bmaList[[1]])),"Std Mean", "Std Coef")

  if(is.null(xNames)){
    stop("the bma objects have to have (the same) rownames attached to them")
  }
  if(!(comp %in% compNames)){
    stop("Please specify comp as one of PIP, Post Mean, Post SD, Std Mean, or Std Coef")
  }

  if(comp=="Std Mean"){
    compMatrix=sapply(bmaList,function(x) estimates.bma(x,std.coefs=TRUE,exact=exact)[xNames,"Post Mean"])
    comp="Standardized Coefficients"
  } else if (comp=="Std SD") {
    compMatrix=sapply(bmaList,function(x) estimates.bma(x,std.coefs=TRUE,exact=exact)[xNames,"Post SD"])
    comp="Standardized SD"    
  } else{
    compMatrix=sapply(bmaList,function(x) estimates.bma(x,exact=exact)[xNames,comp])
  }
  
  bmaNames=names(list(...))[bmaix]
  colnames(compMatrix)=paste("Model", 1:bmaNr)
  if(!is.null(bmaNames) && (length(bmaNames)==ncol(compMatrix))) {
    for (bix in 1:bmaNr) { 
      colnames(compMatrix)[[bix]] <- ifelse(bmaNames[[bix]]=="",paste("Model", bix), bmaNames[[bix]])
    }
  }
  
  # in case you do not want to plot the whole stuff but only the first varNr regressors
  if(!is.null(varNr)){
    compMatrix=compMatrix[varNr,,drop=FALSE]
  }
  
  
  # do the plot ####################
  if (as.logical(do.par)) {
    oldmar=graphics::par()$mar
    spaceforxaxis=graphics::strwidth(rownames(compMatrix)[which.max(nchar(rownames(compMatrix)))],units="inches", cex=cex.xaxis)*(graphics::par("mar")/graphics::par("mai"))[[2]]
    tempmar=oldmar; tempmar[1]=min(max(oldmar[1],spaceforxaxis+oldmar[1]/3), .5*graphics::par("fin")[[2]]*(graphics::par("mar")/graphics::par("mai"))[[1]])
    graphics::par(mar=tempmar)
  }
  

  eval(as.call(c(list(as.name("matplot"),as.name("compMatrix")),as.list(dotargs))))
  
  #matplot(compMatrix,main=main,type=type,col=col,cex=cex,bg=bg,ylab=ylab,xlab=xlab,xaxt="n",pch=pch,lty=lty,lwd=lwd)
  if (as.logical(include.legend)) {
    extractfromdotargs=function(...) { dal=list(...); return(list(col=dal$col,pch=dal$pch)) }
    myargs=eval(as.call(c(list(as.name("extractfromdotargs")),as.list(dotargs))))
    graphics::legend("topright", colnames(compMatrix),pch=myargs$pch,col=myargs$col,bty="n")
  }
  if (as.logical(add.grid)) graphics::grid()
  graphics::axis(1, las=2, at = 1:nrow(compMatrix), labels = rownames(compMatrix),cex.axis=cex.xaxis)
  #graphics::layout(matrix(1))
  if (as.logical(do.par)) graphics::par(mar=oldmar)
}



    
#' Plot Posterior Model Size and Model Probabilities
#' 
#' Produces a combined plot: upper row shows prior and posterior model size
#' distribution, lower row shows posterior model probabilities for the best
#' models
#' 
#' 
#' @param x an object of class 'bma'
#' @param \dots additional arguments for \code{\link{matplot}}
#' @return combines the plotting functions \code{\link{plotModelsize}} and
#' \code{\link{plotConv}}
#' @note The upper plot shows the prior and posterior distribution of model
#' sizes (\code{\link{plotModelsize}}).\cr The lower plot is an indicator of
#' how well the bma object has converged (\code{\link{plotConv}}).
#'  and Paul
#' @seealso \code{\link{plotModelsize}} and \code{\link{plotConv}}
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords hplot
#' @examples
#' 
#' data(datafls)
#' mm=bms(datafls,user.int=FALSE)
#' 
#' plot(mm)
#' @export
plot.bma <-function(x,...) {
   # does a combined plot of plotConv and plotModelsize for bma object bmao
    # bmao: bma object
    # ... passes on parameters to these two functions
    
    if (!is.bma(x)) stop("Need to provide object of class 'bma'!")
    if (x$arguments$nmodel<3) {
      try(plotModelsize(x,...),silent=TRUE)
    } else {
     graphics::layout(matrix(1:2,2,1))
     try(plotModelsize(x,...),silent=TRUE)
     try(plotConv(x,...),silent=TRUE)
     graphics::layout(1)
    }
}


#' Plot Signs of Best Models
#' 
#' Plots a grid with signs and inclusion of coefficients vs. posterior model
#' probabilities for the best models in a 'bma' object:
#' 
#' Under default settings, blue corresponds to positive sign, red to a negative
#' sign, white to non-inclusion.
#' 
#' @param x a list of class bma (cf. \code{\link{bms}} for further details)
#' @param yprop2pip if \code{yprop2pip=TRUE} then the grid lines on the
#' vertical axis are scaled according to the coefficients' inclusion
#' probabilites.\cr If \code{yprop2pip=FALSE} (default) then the grid lines on
#' the vertical axis are equidistant.
#' @param order.by.pip with \code{order.by.pip=TRUE} (default), coefficients
#' are sorted according to their posterior inclusion probabilites along the
#' vertical axis. If \code{order.by.pip=FALSE} they are ordered as they were
#' provided to \code{\link{bms}}.
#' @param do.par Defaults to \code{do.par=TRUE}, which adjusts
#' \code{\link{par}()$mar} for optimal positioning. Set \code{do.par=FALSE} for
#' customizing \code{par} yourself.
#' @param do.grid \code{do.grid=TRUE} (default) plots grid lines among the
#' chart's boxes, akin to the low level command \code{\link{grid}}.
#' \code{do.grid=FALSE} omits the grid lines.
#' @param do.axis \code{do.axis=TRUE} (default) plots axis tick marks and
#' labels (cf. \code{\link{axis}}). \code{do.axis=FALSE} omits them.
#' @param cex.axis font size for the axes (cf. \code{\link{axis}}), defaults to
#' 1
#' @param \dots Parameters to be passed on to \code{\link{image.default}}.
#' 
#' @seealso \link{coef.bma} for the coefficients in matrix form, \link{bms} for
#' creating 'bma' objects.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords hplot
#' @examples
#' 
#'  data(datafls)
#'  
#'  model=bms(datafls,nmodel=200)
#'  
#'  #plot all models
#'  image(model,order.by.pip=FALSE)
#'  image(model,order.by.pip=TRUE,cex.axis=.8)
#'  
#'  #plot best 7 models, with other colors
#'  image(model[1:7],yprop2pip=TRUE,col=c("black","lightgrey"))
#'  
#' @export
image.bma <- function(x,yprop2pip=FALSE,order.by.pip=TRUE,do.par=TRUE,do.grid=TRUE,do.axis=TRUE,cex.axis=1,...) { #,main=NULL,col=c("tomato","blue"),xlab="Cumulative Model Probabilities",ylab=""
 #does a 'grid plot' of bmao's best models (contained in bmao$topmod) 
 #putting the variables' (on the y-axis) coefficent signs per model vs. the post. model probs (x-axis)
  # bmao: bma object
  # yprop2pip: if TRUE, then horizontal grid lines (on the y axis) are proportional to the variables PIPs; if FALSE, they are equidistanrt
  # do.par: if TRUE adjust par temporarily to fit labels properly; if FALSE, please take care of the par()$mar parameter yourself
  # do.grid: if TRUE, highlights the boundaries between coefficient areas via the grid() function; if FALSE, it does not
  # do.axis: if TRUE, plots default axis tickmarks and labels; if FALSE, it does not (consider using axis() afterwards)
  # cex.axis: denotes label font size (cf. ?axis )  
  # ... parameters passed on to image.default (Note: default colors are col=c(4,2), default axis labels are empty strings)  

  dotargs = match.call(expand.dots=FALSE)$...
  ests=estimates.bma(x,exact=TRUE,order.by.pip=order.by.pip,include.constant=FALSE)
  ests=ests[nrow(ests):1,]
  pips=ests[,"PIP"]
  idx=ests[,"Idx"]
  pmp.res=pmp.bma(x,oldstyle=TRUE)
  pmps=pmp.res[,1]
  normali_factor=pmin(sum(pmp.res[,2]),1) #this mnormali_factor business only matters for MCMCl, not for iteration
  betasigns=beta.draws.bma(x)[idx,,drop=FALSE]
  betasigns=betasigns[as.logical(pips),]
  betasigns=sign(betasigns)/2+.5
  betasigns[betasigns==.5]=NA
  pips=pips[as.logical(pips)]
  if (yprop2pip) {
    pipbounds=(c(0,cumsum(pips)))  
  } else {
    pipbounds=0:length(pips)
    names(pipbounds)=c("",names(pips))
  }
  pmpbounds=(c(0,cumsum(pmps)))  

  if (do.par) {
    oldmar=graphics::par()$mar    
    spaceforyaxis=graphics::strwidth(names(pipbounds)[which.max(nchar(names(pipbounds)))],units="inches")*(graphics::par("mar")/graphics::par("mai"))[[2]]
    tempmar=oldmar; tempmar[2]=min(spaceforyaxis+oldmar[2]/2, .5*graphics::par("fin")[[1]]*(graphics::par("mar")/graphics::par("mai"))[[2]])
    
    graphics::par(mar=tempmar)
  }

  dotargs=.adjustdots(dotargs, ylab="", xlab="Cumulative Model Probabilities", col=c("tomato", "blue"), main = paste("Model Inclusion Based on Best ",length(pmps), " Models"))
  dotargs$axes <- FALSE
#  if (!is.element("ylab",names(dotargs))) { dotargs$ylab <- "" }
#  if (!is.element("xlab",names(dotargs))) { dotargs$xlab <- "Cumulative Model Probabilities" }
#  if (!is.element("col",names(dotargs))) { dotargs$col <- c("tomato", "blue") }
#  if (!is.element("main",names(dotargs))) { dotargs$main <- paste("Model Inclusion Based on Best ",length(pmps), " Models") }


  tbetasigns=t(betasigns)
  eval(as.call(c(list(as.name("image.default"),as.name("pmpbounds"),as.name("pipbounds"),as.name("tbetasigns")),as.list(dotargs))))

  #if (is.null(main)) { main=paste("Model Inclusion Based on Best ",length(pmps), " Models") }
  #image.default(pmpbounds,pipbounds,t(betasigns),col=col,axes=FALSE,xlab=xlab,ylab=ylab,main=main,...)

  if (do.axis) {
    graphics::axis(1,at=pmpbounds, labels=round(normali_factor*pmpbounds,2),cex.axis=cex.axis)
    graphics::axis(2,at=pipbounds,labels=FALSE,line=FALSE)    
    graphics::axis(2,at=pipbounds[-1]-diff(pipbounds)/2,labels=names(pipbounds[-1]),tick=FALSE,las=1,cex.axis=cex.axis)
  }

  if (do.grid) {
    graphics::abline(v=round(pmpbounds,2),lty="dotted",col="grey")
    graphics::abline(h=round(pipbounds,2),lty="dotted",col="grey")
  }
  if (do.par) {graphics::par(mar=oldmar)}
}






   ##########################
   # LINEAR MODELS (NO BMA) #
###########################################################################


#' Bayesian Linear Model with Zellner's g
#' 
#' Used to fit the Bayesian normal-conjugate linear model with Zellner's g
#' prior and mean zero coefficient priors. Provides an object similar to the
#' \code{\link{lm}} class.
#' 
#' \code{zlm} estimates the coefficients of the following model \eqn{y = \alpha
#' + X \beta + \epsilon} where \eqn{\epsilon} ~ \eqn{N(0,\sigma^2)} and \eqn{X}
#' is the design matrix\cr The priors on the intercept \eqn{\alpha} and the
#' variance \eqn{\sigma} are improper: \eqn{alpha \propto 1}, \eqn{sigma
#' \propto \sigma^{-1}} \cr Zellner's g affects the prior on coefficients:
#' \eqn{beta} ~ \eqn{N(0, \sigma^2 g (X'X)^{-1})}. \cr Note that the prior mean
#' of coefficients is set to zero by default and cannot be adjusted. Note
#' moreover that \code{zlm} always includes an intercept.
#' 
#' @param formula an object of class "formula" (or one that can be coerced to
#' that class), such as a data.frame - cf. \code{\link{lm}}
#' @param data an optional \code{\link{data.frame}} (or one that can be coerced
#' to that class): cf. \code{\link{lm}}
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param g specifies the hyperparameter on Zellner's g-prior for the
#' regression coefficients.\cr \code{g="UIP"} corresponds to \eqn{g=N}, the
#' number of observations (default); \code{g="BRIC"} corresponds to the
#' benchmark prior suggested by Fernandez, Ley and Steel (2001), i.e
#' \eqn{g=max(N, K^2)}, where K is the total number of covariates;\cr
#' \code{g="EBL"} estimates a local empirical Bayes g-parameter (as in Liang et
#' al. (2008));\cr \code{g="hyper"} takes the 'hyper-g' prior distribution (as
#' in Liang et al., 2008) with the default hyper-parameter \eqn{a=3}; This
#' hyperparameter can be adjusted (between \eqn{2<a<=4}) by setting
#' \code{g="hyper=2.9"}, for instance.\cr Alternatively, \code{g="hyper=UIP"}
#' sets the prior expected value of the shrinkage factor equal to that of UIP
#' (above), \code{g="hyper=BRIC"} sets it according to BRIC
#' @return Returns a list of class \code{zlm} that contains at least the
#' following elements (cf. \code{\link{lm}}):
#' 
#' \item{coefficients}{a named vector of posterior coefficient expected values}
#' \item{residuals}{the residuals, that is response minus fitted values}
#' \item{fitted.values}{the fitted mean values} \item{rank}{the numeric rank of
#' the fitted linear model} \item{df.residual}{the residual degrees of freedom}
#' \item{call}{the matched call} \item{terms}{the \code{\link{terms}} object
#' used} \item{model}{the model frame used} \item{coef2moments}{a named vector
#' of coefficient posterior second moments} \item{marg.lik}{the log marginal
#' likelihood of the model} \item{gprior.info}{a list detailing information on
#' the g-prior, cf. output value \code{gprior.info} in \code{\link{bms}}}
#' @author Stefan Zeugner
#' @seealso The methods \code{\link{summary.zlm}} and \code{\link{predict.lm}}
#' provide additional insights into \code{zlm} output.\cr The function
#' \code{\link{as.zlm}} extracts a single out model of a \code{bma} object (as
#' e.g. created through\code{\link{bms}}).\cr Moreover, \code{\link{lm}} for
#' the standard OLS object, \code{\link{bms}} for the application of \code{zlm}
#' in Bayesian model averaging.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @references The representation follows Fernandez, C. E. Ley and M. Steel
#' (2001): Benchmark priors for Bayesian model averaging. Journal of
#' Econometrics 100(2), 381--427
#' 
#' See also \url{http://bms.zeugner.eu} for additional help.
#' @keywords models
#' @examples
#' 
#' 
#' data(datafls)
#' 
#' #simple example
#' foo = zlm(datafls)
#' summary(foo)
#' 
#' #example with formula and subset
#' foo2 = zlm(y~GDP60+LifeExp, data=datafls, subset=2:70) #basic model, omitting three countries
#' summary(foo2)
#' 
#' 
#' @export
zlm <- function(formula, data=NULL, subset=NULL, g="UIP") {
    #does a normal-gamma linear Bayesian model with Zellner's g prior
    #INPUTS:
    # formula, data, subset: cf help(lm)
    # g: a g prior character (cf help(bms)), or a gprior object as from .choose.gprior
    #OUTPUT:
    # an object class c("zlm","lm") that is comparable to the lm object

    thiscall=match.call()

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if (is.matrix(formula)) {
      mf <- model.frame(as.data.frame(formula,drop.unused.levels=TRUE))
    } else {
      mf <- eval(mf, parent.frame())
    }


    yXdata=as.matrix(mf);
    N=nrow(yXdata); K=ncol(yXdata)-1
    dmdata=yXdata-matrix(colMeans(yXdata),N,K+1,byrow=TRUE)

    yty=c(crossprod(dmdata[,1]));


    olsres=.ols.terms2(positions=rep(TRUE,K),yty=yty,N=N,K=K,XtX.big=crossprod(dmdata[,-1,drop=FALSE]),Xty.big=c(crossprod(dmdata[,-1,drop=FALSE],dmdata[,1])))$full.results()

    #get gprior
    if (is.list(g)) {
      if (any(is.element(names(g),"gtype"))) gprior.info=g else stop("Please provide a proper g-prior. see help(zlm)")
    }
    gprior.info=.choose.gprior(g=g,N=N,K=K,return.g.stats=TRUE,yty=yty)
    lprobcalc=gprior.info$lprobcalc
    
#     if (gprior.info$gtype=="EBL")  {
#        lprobcalc=.lprob.eblocal.init(N=N,K=K,yty=yty,return.g=gprior.info$return.g.stats)
#     } else if (gprior.info$gtype=="hyper")  {
#        lprobcalc=.lprob.hyperg.init(N=N,K=K,yty=yty,f21a=gprior.info$hyper.parameter,return.gmoments=gprior.info$return.g.stats)
#     } else {
#        lprobcalc=.lprob.constg.init(g=gprior.info$g,N=N,K=K,yty=yty)
#     }

    zres = lprobcalc$lprob.all(ymy=olsres$ymy, k=K, bhat=olsres$bhat, diag.inverse=olsres$diag.inverse)
    betas = c(zres$b1)
    betas2 = c(zres$b2)
    alpha = mean(yXdata[,1]) - c(crossprod(betas, colMeans(yXdata)[-1]))
    fitval = c(yXdata[, -1,drop=FALSE]%*%betas)+alpha
    resids = yXdata[,1]-fitval

    if (gprior.info$is.constant) {
       gprior.info$shrinkage.moments = 1 - 1/(1+gprior.info$g)
    } else {
       gprior.info$shrinkage.moments = zres$otherstats
    }

    #lm-like stuff
    mt <- attr(mf, "terms")
    alphabeta=c(alpha,betas)
    names(alphabeta) <- c("(Intercept)", attr(mt,"term.labels"))

    res=list()
    res$coefficients <- alphabeta
    res$residuals <- resids
    res$rank <- K+1
    res$fitted.values <- fitval
    res$df.residual <- N-K-1
    res$xlevels <- stats::.getXlevels(mt, mf)
    res$call <- thiscall
    res$terms <- mt
    res$model <- mf
    res$na.action <- attr(mf,"na.action")
    res$coef2moments <- c(NA,betas2)
    res$marg.lik <- zres$lprob
    res$gprior.info <- gprior.info
    res$olsres <- olsres
    res$zres <- zres

    class(res)=c("zlm","lm")
    return(res)

}

#' Summarizing Linear Models under Zellner's g
#' 
#' summary method for class "\code{zlm}"
#' 
#' \code{summary.zlm} prints out coefficients expected values and their
#' standard deviations, as well as information on the gprior and the log
#' marginal likelihood. However, it invisibly returns a list with elements as
#' described below:
#' 
#' @param object an object of class \code{zlm}: see "Examples" below
#' @param printout If \code{TRUE} (default, then information is printed to
#' console in a neat form
#' @param \dots further arguments passed to or from other methods
#' @return A \code{\link{list}} with the following elements \item{residuals}{
#' The expected value of residuals from the model} \item{coefficients}{The
#' posterior expected values of coefficients (including the intercept) }
#' \item{coef.sd}{Posterior standard deviations of the coefficients (the
#' intercept SD is \code{NA}, since an improper prior was used)}
#' \item{gprior}{The g prior as it has been submitted to \code{object}}
#' \item{E.shrinkage}{the shrinkage factor \eqn{g/(1+g)}, respectively its
#' posterior expected value in case of a hyper-g prior}
#' \item{SD.shrinkage}{(Optionally) the shrinkage factor's posterior standard
#' deviation (in case of a hyper-g prior)} \item{log.lik}{The log marginal
#' likelihood of the model}
#' @author Stefan Zeugner
#' @seealso \code{\link{zlm}} for creating \code{zlm} objects,
#' \code{link{summary.lm}} for a similar function on OLS models
#' 
#' See also \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#' data(datafls)
#' 
#' #simple example
#' foo = zlm(datafls)
#' summary(foo)
#' 
#' sfoo = summary(foo,printout=FALSE)
#' print(sfoo$E.shrinkage)
#' 
#' @export
summary.zlm <- function(object, printout=TRUE, ...) {
    #prints a summary for a zlm object: coefficients, std devs, shinkage stats, marg.lik.

    betas=object$coefficients
    betas2= object$coef2moments
    sds= sqrt(betas2-betas^2)
    ests=cbind(betas,sds)
    gi=object$gprior.info
    gi.choice = gi$gtype
    if (gi$gtype=="hyper") {gi.choice=paste(gi$gtype," (a=",2+signif(gi$hyper.parameter-2,digits=4),")",sep="")}
    gi.sd=-1; gi.sdtext=""
    if (length(gi$shrinkage.moments)>1) {
       gi.sd=sqrt(gi$shrinkage.moments[[2]]-gi$shrinkage.moments[[1]]^2)
       gi.sdtext=paste(" St.Dev.:", round(gi.sd,3))
    }


    rownames(ests) = c("(Intercept)",attr(object$terms,"term.labels"))
    colnames(ests)=c("Exp.Val.","St.Dev.")
    cat("Coefficients\n")
    print(ests)
    cat("\n Log Marginal Likelihood:\n")
    cat(object$marg.lik)
    cat(paste("\n g-Prior:", gi.choice,"\n"))
    cat(paste("Shrinkage Factor",ifelse(gi$is.constant,": ", " Exp.Val: "), round(gi$shrinkage.moments[[1]],3), gi.sdtext, "\n",sep=""))
    
    res=list()
    res$residuals <- object$residuals
    res$coefficients <- object$coefficients
    res$coef.sd <- sds
    res$gprior <- gi.choice
    res$E.shrinkage <- gi$shrinkage.moments[[1]]
    if (gi.sd>-1) {res$SD.shrinkage <- gi.sd}
    res$log.lik <- object$marg.lik
    return(invisible(res))
}



#' Extract a Model from a bma Object
#' 
#' Extracts a model out of a \code{bma} object's saved models and converts it
#' to a \code{\link{zlm}} linear model
#' 
#' A bma object stores several 'best' models it encounters (cf. argument
#' \code{nmodel} in \code{\link{bms}}). \code{as.zlm} extracts a single model
#' and converts it to an object of class \code{\link{zlm}}, which represents a
#' linear model estimated under Zellner's g prior.\cr The utility
#' \code{\link{model.frame}} allows to transfrom a \code{zlm} model into an OLS
#' model of class \code{\link{lm}}.
#' 
#' @param bmao A \code{bma} object, e.g. resulting from a call to
#' \code{\link{bms}}
#' @param model The model index, in one of the following forms:\cr An integer,
#' denoting the rank of the model (1 for best, 2 for second-best, ...)\cr A
#' numeric or logical vector of length K describing which covariates are
#' contained in the model\cr A hexcode character describing which covariates
#' are contained in the model
#' @return a list of class \code{\link{zlm}}
#' @author Stefan Zeugner
#' @seealso \code{\link{bms}} for creating \code{bma} objects,
#' \code{\link{zlm}} for creating \code{zlm} objects,
#' \code{\link{pmp.bma}} for displaying the
#' topmodels in a \code{bma} object
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords models
#' @examples
#' 
#' data(datafls)
#' 
#' mm=bms(datafls[,1:6],mcmc="enumeration") # do a small BMA chain
#' topmodels.bma(mm)[,1:5] #display the best 5 models
#' 
#' m2a=as.zlm(mm,4) #extract the fourth best model
#' summary(m2a)
#' 
#' # Bayesian Model Selection:
#' # transform the best model into an OLS model:
#' lm(model.frame(as.zlm(mm)))
#' 
#' # extract the model only containing the 5th regressor
#' m2b=as.zlm(mm,c(0,0,0,0,1)) 
#' 
#' # extract the model only containing the 5th regressor in hexcode
#' print(bin2hex(c(0,0,0,0,1)))
#' m2c=as.zlm(mm,"01")
#' 
#' 
#' 
#' 
#' @export
as.zlm <- function(bmao, model=1) {
   #this function extracts a single topmodel from a bma object and converts it to zlm format
   #Inputs:
   # bmao: bma object
   # model: index of the model
   #Output: zlm object

   thiscall=match.call()
   if (!is.bma(bmao)) stop("bmao needs to be a bma object")

   bools=bmao$topmod$bool()
   if (all(is.character(model))&&length(model)==1) {
       model=(1:length(bools))[bools==model[[1]]]
       if (length(model)==0) stop("Provided model hex-index was not found in bmao object topmodels")
   } else if (all(is.character(model))&&(length(model)>1) ) {
       mix=match(model, bmao$reg.names)
       if (any(is.na(mix))) stop("Provided variable names do not conform to bma object")
       ll=logical(bmao$info$K); ll[mix]=TRUE; model=(1:length(bools))[bools==bin2hex(ll)]; rm(ll,mix)
       if (length(model)==0) stop("Model conforming to provided variable names was not found in bmao object topmodels")
   } else if ((length(model)==bmao$info$K)&&(is.numeric(model)||is.logical(model))) {
       model=(1:length(bools))[bools==bin2hex(model)]
       if (length(model)==0) stop("Provided binary model index was not found in bmao object topmodels")
   } else if ((length(model)==1)&&(is.numeric(model)||is.logical(model))) {
       if (model<1|model>length(bools)) stop("Provided numeric model index was not found in bmao object topmodels")
   } else stop("model needs to be an integer, logical or character model index representation (hexcode or variable names)")


   inclvbls= as.logical(bmao$topmod$bool_binary()[,model, drop=TRUE])

   yXdf =as.data.frame(bmao$arguments$X.data)

   zlmres=zlm(as.formula(yXdf[,c(TRUE,inclvbls)]),data=yXdf,g=bmao$gprior.info)

   zlmres$call <- thiscall
   return(zlmres)
}


#' Predict Method for zlm Linear Model
#' 
#' Expected value (And standard errors) of predictions based on 'zlm' linear
#' Bayesian model under Zellner's g prior
#' 
#' 
#' @param object a zlm linear model object - see \code{\link{zlm}}
#' @param newdata An optional data.frame, matrix or vector containing variables
#' with which to predict. If omitted, then (the expected values of) the fitted
#' values are returned.
#' @param se.fit A switch indicating if the standard deviations for the
#' predicted varaibles are required.
#' @param \dots further arguments passed to or from other methods.
#' @return A vector with (expected values of) fitted values.\cr If
#' \code{se.fit} is \code{TRUE}, then the output is a list with the following
#' elements: \item{fit}{ a vector with the expected values of fitted values}
#' \item{std.err}{ a vector with the standard deviations of fitted values}
#' \item{se.fit}{ a vector with the standard errors without the residual scale
#' akin to \code{se.fit} in \code{\link{predict.lm}} } \item{residual.scale}{
#' The part from the standard deviations that involves the identity matrix.
#' Note that \code{sqrt(se.fit^2+residual.scale^2)} yields \code{std.err}. }
#' 
#' @seealso \code{\link{bms}} for creating zlm objects,
#' \code{\link{predict.lm}} for a comparable function,
#' \code{\link{predict.bma}} for predicting with bma objects
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'  data(datafls)
#'  mm=zlm(datafls,g="EBL")
#'  
#'  predict(mm) #fitted values 
#'  predict(mm, newdata=1:41) #prediction based on a 'new data point' 
#'  
#'  #prediction based on a 'new data point', with 'standard errors'
#'  predict(mm, newdata=datafls[1,], se.fit=TRUE) 
#'  
#' @export
predict.zlm <- function(object, newdata=NULL, se.fit=FALSE, ...) {
   # does basic fitting in expected values, cf. predict.lm
   # object: a zlm object
   # newdata: newdata to be supplied (just eas in predict.lm)
   # se.fit:  whether standard erros should be calculated
   # output: a vector with fitted values, if se.fit=T, then a list with fit, standard errors and residual scale

    if (!is(object,"zlm")) {stop("you need to provide a zlm object"); return()}


    #get the betas as required
    betas=object$coefficients[-1,drop=FALSE]
    alpha=object$coefficients[[1]]

    #check the newdata argument
    if (is.null(newdata)) {
       newX<-as.matrix(object$model[,-1,drop=FALSE])
    } else {
       newX=as.matrix(newdata)
       if (!is.numeric(newX)) stop("newdata must be numeric!")
       if (is.vector(newdata)) newX=matrix(newdata,1)
       if (ncol(newX)!=length(betas)) {
         if (ncol(newX)==length(betas)+1) {
             newX=newX[,-1,drop=FALSE] # this is to achieve a behavior similar to predict.lm in this case
         } else {
           stop("newdata must be a matrix or data.frame with ", length(betas), " columns.")
         }
       }
    }

    #if se.fit==FALSE, we are done now:
    if (!se.fit) return(as.vector(newX%*%betas)+alpha)

    #################################################
    #if se.fit==TRUE, get additional info to compute the standard errors
    yXdata <- as.matrix(object$model)
    oldXraw <- yXdata[,-1,drop=FALSE]
    if (!is.null(colnames(newX))&& !is.null(colnames(oldXraw))) { if ( all(colnames(oldXraw) %in% colnames(newX)) && !all(colnames(oldXraw) == colnames(newX))  ) { #this is another user check, when the stuff came in malformed
      warning("argument newdata had to be reordered according to its column names. Consider submitting the columns of newdata in the right order.")
      newX=newX[,colnames(oldXraw), drop=FALSE] 
    } }
    yraw <- yXdata[,1,drop=TRUE]
    N <- length(yraw)
    k <- ncol(oldXraw)
    oldXmeans <- colMeans(oldXraw)
    oldXdm <- oldXraw-matrix(oldXmeans,N,k,byrow=TRUE)
    newXdm <- newX-matrix(oldXmeans,nrow(newX),k,byrow=TRUE)

    #compute x_f' (X'X)^-1 x_f
    xtxinv=chol2inv(chol(crossprod(oldXdm)))
    xtxinv_xf= tcrossprod(xtxinv,newXdm)
    xf_xx_xf=unlist(lapply(1:nrow(newXdm),function(x) {crossprod(newXdm[x,],xtxinv_xf[,x])[[1L]]} ) )

    #these are some factors multipliying the above
    bvar=object$coef2moments[-1]-object$coefficients[-1]^2
    bvar_factor=bvar[[1L]]/xtxinv[[1L]]
    yty=as.vector(crossprod(yraw)-N*mean(yraw)^2)
    r2 = 1-object$olsres$ymy/yty

    if (object$gprior.info$gtype=="hyper") {
      f21a=object$gprior.info$hyper.parameter
      f21_recover = exp( (object$marg.lik) + (N-1)/2*log(yty) + log((k+f21a-2)/(f21a-2)) ) #this might be not so accurate due to numerical errors, but in this case the part involving F21 is usually very very small...
      res_scale = yty/(N-3)/(N-1-k-f21a)*( (N-3)*(1-r2) - (k+f21a-2)/f21_recover)
      svar_woscale = res_scale/N +  bvar_factor*xf_xx_xf
      #svar_woscale is the standard error without the residual scale (the identiy matrix part)
    } else {
      sf=object$gprior.info$shrinkage.moments[[1]]
      res_scale = (1-sf*r2)*yty/(N-3)
      svar_woscale = res_scale/N + bvar_factor*xf_xx_xf
    }

    #produce output
    reslist=list()
    reslist$fit <- as.vector(newX%*%betas)+alpha
    reslist$std.err <- sqrt(svar_woscale+res_scale)
    reslist$se.fit <- sqrt(svar_woscale)
    reslist$residual.scale <- sqrt(res_scale)

    return(reslist)
}


#' @export
density.zlm <- function(x,reg=NULL,addons="lesz",std.coefs=FALSE,n=300,plot=TRUE,hnbsteps=30,addons.lwd=1.5,...) {
    #this function does just the same as density.bma, but for an object of class zlm
    #permitted values for addons: e, s, l, z, g

    addons=gsub("E","",addons,ignore.case=FALSE)
    addons=gsub("S","",addons,ignore.case=FALSE)
    addons=gsub("b","",addons,ignore.case=FALSE)
    addons=gsub("m","",addons,ignore.case=FALSE)
    addons=gsub("p","",addons,ignore.case=FALSE)

    N=length(x$residuals); K=length(x$coefficients)-1

    tmo=topmod(1,nmaxregressors=K, bbeta=TRUE,liks=x$marg.lik, ncounts=1, modelbinaries=matrix(rep(1,K),K,1), betas=matrix(as.vector(x$coefficients[-1]),K), betas2=matrix(as.vector(x$coef2moments[-1]),K))
    tokenbma=list(info=list(K=K,N=N),arguments=list(),topmod=tmo,start.pos=integer(0), gprior.info=x$gprior.info, X.data=x$model, reg.names=names(x$coefficients)[-1], bms.call=new("call"))
    class(tokenbma)="bma"
    return(density.bma(tokenbma,reg=reg, addons=addons, std.coefs=std.coefs, n=n, plot=plot, hnbsteps=hnbsteps, addons.lwd=addons.lwd,...))
}



 .adjustdots = function(dotargs,...) {
             # helper function to adjust default arguments for passing on to other functions
             # in particular to adjust type, col, main etc. when passing on to plotting functions
            defargs=list(...); defargnames=names(defargs);
            dotargs=as.list(dotargs);
            
            if (is.null(dotargs)) { dotargs= list() }
            for (di in seq_len(length(defargs))) {
                if (!is.element(defargnames[[di]],names(dotargs))) {
                   dotargs[[defargnames[[di]]]] <- defargs[[di]];
                }
            }
            return(dotargs);
 }


#####################################
# PREDICTIVE DENSITY           ######
#####################################

#' Predictive Densities for bma Objects
#' 
#' Predictive densities for conditional forecasts
#' 
#' The predictive density is a mixture density based on the \code{nmodels} best
#' models in a \code{bma} object (cf. \code{nmodel} in \code{\link{bms}}).\cr
#' The number of 'best models' to retain is therefore vital and should be set
#' quite high for accuracy.
#' 
#' @aliases pred.density pred.density-class print.pred.density
#' @param object a bma object - see \code{\link{bms}}, alternativel a
#' \code{\link{zlm}} object
#' @param newdata A data.frame, matrix or vector containing variables with
#' which to predict.
#' @param n The integer number of equally spaced points at which the density is
#' to be estimated.
#' @param hnbsteps The number of numerical integration steps to be used in case
#' of a hyper-g prior (cf. argument \code{g} in \code{\link{bms}}). Increase
#' this number to increase accuracy. Must be an even integer.
#' @param \dots arguments to be passed on to \code{\link{plot.density}}.
#' @return \code{pred.density} returns a list of class \code{pred.density} with
#' the following elements \item{densities()}{a list whose elements each contain
#' the estimated density for each forecasted observation} \item{fit}{a vector
#' with the expected values of the predictions (the 'point forecasts')}
#' \item{std.err}{a vector with the standard deviations of the predictions (the
#' 'standard errors')} \item{dyf(realized.y, predict_index=NULL)}{Returns the
#' densities of realized response variables provided in \code{realized.y}. \cr
#' If \code{realized.y} is a matrix, then each row corresponds to a forecast
#' observation in \code{newdata}\cr if not left empty, \code{predict.index}
#' specifies to which observations in newdata the realized.y should apply}
#' \item{lps(realized.y, predict_index=NULL)}{Computes the log predictive score
#' for the response varaible provided in \code{realized.y} (cf.
#' \code{\link{lps.bma}}) -\cr Note that the LPS equals minus the mean of the
#' logarithmized results from \code{dyf}) } \item{plot((x, predict_index =
#' NULL, addons = "eslz", realized.y = NULL, addons.lwd = 1.5, ...)}{the same
#' as \code{plot.pred.density}} \item{n}{The number of equally spaced
#' points for which the density (under \code{densities()} was computed.}
#' \item{nmodel}{The number of best models predictive densities are based
#' upon.} \item{call}{the call that created this \code{pred.density} object}
#' @note In BMS version 0.3.0, \code{pred.density} may only cope with built-in
#' \code{gprior}s, not with any user-defined priors.
#' 
#' @seealso \code{\link{predict.bma}} for simple point forecasts,
#' \code{plot.pred.density} for plotting predictive densities,
#' \code{\link{lps.bma}} for calculating the log predictive score
#' independently, \code{\link{quantile.pred.density}} for extracting quantiles
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'  data(datafls)
#'  mm=bms(datafls,user.int=FALSE)
#'  
#'  #predictive densityfor two 'new' data points
#'  pd=pred.density(mm,newdata=datafls[1:2,]) 
#'  
#'  
#'  #fitted values based on best models, same as predict(mm, exact=TRUE)
#'  pd$fit
#'  
#'  #plot the density for the first forecast observation
#'  plot(pd,1)  
#'  
#'  # the same plot ' naked'
#'  plot(pd$densities()[[1]])
#'  
#'  
#'  #predict density for the first forecast observation if the dep. variable is 0
#'  pd$dyf(0,1) 
#'  
#'  #predict densities for both forecasts for the realizations 0 and 0.5
#'  pd$dyf(rbind(c(0,.5),c(0,.5)))
#'  
#'  # calc. Log Predictive Score if both forecasts are realized at 0:
#'  lps.bma(pd,c(0,0))
#'  
#'  
#' @export
pred.density <- function(object, newdata=NULL, n=300, hnbsteps=30, ...) {
    
    dtcm=function(x,df,ncp,varp) {
      #a wrapper for univariate t-dist with non-centrality parameter and variance parameter
      # (variance is df/(df-2)*varp)
      sqvarp=sqrt(varp)
      stats::dt((x-ncp)/sqvarp,df=df)/sqvarp
    }

    dsgivenykernel <- function(sf,kpa,N,z) {
      #the post. density of the shrinkge factor f(s|Y)*F((N-1)/2,1,(k+a)/2,R2)
      (kpa-2)/2*(1-sf)^((kpa-4)/2)*(1-sf*z)^(-(N-1)/2)
    }

    #user checks
    nbsteps=max(hnbsteps,2)
    n=max(ceiling(n),1)

    #read out data from object
    is.hyper=(object$gprior.info$gtype=="hyper")
    if (is.hyper) f21a=object$gprior.info$hyper.parameter
    if (is.bma(object)) {
      K=object$info$K; N=object$info$N
      yXdata=as.matrix(object$arguments$X.data)
      tmo <- object$topmod
    } else if (is(object,"zlm")) {
      yXdata=as.matrix(object$model)
      K=ncol(yXdata)-1; N=nrow(yXdata)
      tmo <- topmod(1,nmaxregressors=K, bbeta=TRUE,liks=object$marg.lik, ncounts=1, modelbinaries=matrix(rep(1,K),K,1), betas=matrix(as.vector(object$coefficients[-1]),K), betas2=matrix(as.vector(object$coef2moments[-1]),K))
    } else stop("argument 'object' requires class 'bma' or 'zlm'")

    rm(object)


    #check the newdata argument, checks are small as newdata has already survived predict.zlm
    if (missing(newdata)) {
       stop("You must provide the argument newdata")
    } else {
       newX=as.matrix(newdata)
       if (!is.numeric(newX)) stop("newdata must be numeric!")
       if (is.vector(newdata)) newX=matrix(newdata,1)
       if (ncol(newX)!=K) {
         if (ncol(newX)==K+1) {
             newX=newX[,-1,drop=FALSE] # this is to achieve a bevavior similar to predict.lm in this case
         } else {
           stop("newdata must be a matrix or data.frame with ", K, " columns.")
         }
       }
       orinames=colnames(yXdata[,-1,drop=FALSE]) #this is a user check whether columns had been submitted in the wrong  order
       if (!is.null(colnames(newX)) && !is.null(orinames)) { 
        if (all(orinames %in% colnames(newX) ) && !all(orinames == colnames(newX))  ) {
            warning("argument newdata had to be reordered according to its column names. Consider submitting the columns of newdata in the right order.")
            newX=newX[,orinames, drop=FALSE] 
        }
       }
    }

    if(!is.null(rownames(newX))) {
        newXnames=rownames(newX)
    } else {
        newXnames=as.character(1:nrow(newX))
    }
    rnew = nrow(newX)




    y.mean=mean(yXdata[,1])
    y<-yXdata[,1]-matrix(y.mean,N,1,byrow=TRUE)
    X<-yXdata[,-1,drop=FALSE]-matrix(colMeans(yXdata[,-1,drop=FALSE]),N,K,byrow=TRUE)
    XtX.big=crossprod(X)
    Xty.big=as.vector(crossprod(X,y))
    yty = crossprod(y)[[1]]
    newXdm= newX-matrix(colMeans(yXdata[,-1,drop=FALSE]),rnew,K,byrow=TRUE)


    hexobject<-.hexcode.binvec.convert(K)

    make_xfxxxf = function(hex) {

       syminv <- function(symmat, ndim=ncol(symmat)) {
              #this does the same as chol2inv(chol.default(x)), but is stripped-down for speed purposes
              # Caution: symmat must always have length(symmat)>0!!!
              if (!is.matrix(symmat)) {symmat=as.matrix(symmat)}
              if (dim(symmat)[[1]]==0) return(matrix(numeric(0),0,0))
              return( chol2inv(chol(symmat), size=ndim) )
      }


      boolvec = as.logical(hexobject$as.binvec(hex))
      if (!any(boolvec)) return(c(numeric(rnew),numeric(rnew),Inf,Inf,0))
      newXsub= newXdm[,boolvec,drop=FALSE]
      xtxinv = syminv(XtX.big[boolvec,boolvec,drop=FALSE])
      xty=Xty.big[boolvec]
      betas = as.vector(crossprod(xtxinv,xty),mode="numeric")
      r2=crossprod(xty,betas)[[1]]/yty
      xtxinv_xf = tcrossprod(xtxinv,newXsub)
      xf_xx_xf = unlist(lapply(1:nrow(newXsub),function(x) {crossprod(newXsub[x,],xtxinv_xf[,x])[[1L]]} ) )
      xf_bhat = as.vector(newXsub%*%betas)
      return(c(xf_xx_xf, xf_bhat, xtxinv[[1L]], betas[[1L]],r2))


    }

    pmps=pmp.bma(tmo,oldstyle=TRUE)[,1,drop=TRUE]
    bools=tmo$bool()
    nmodel=length(bools)
    linvres=lapply(bools,make_xfxxxf)
    mat_xfxxxf = array(unlist(lapply(linvres,"[",1:rnew)),dim=c(rnew,nmodel))
    mat_xfbhat = array(unlist(lapply(linvres,"[",rnew+(1:rnew))),dim=c(rnew,nmodel))
    xtxinv_elem1 = unlist(lapply(linvres,"[[",rnew*2+1))
    betahat_elem1 = unlist(lapply(linvres,"[[",rnew*2+2))
    r2 = unlist(lapply(linvres,"[[",rnew*2+3))

    kvec=tmo$kvec_raw()
    kvec_cs=c(1,cumsum(kvec)+1); kvec_cs=kvec_cs[-length(kvec_cs)]
    firstbetas=tmo$betas_raw()[kvec_cs]
    firstbetas2=tmo$betas2_raw()[kvec_cs]
    Es = firstbetas / betahat_elem1
    varmult = (firstbetas2-firstbetas^2) / xtxinv_elem1
    if (is.hyper) {
      first_factor= yty/(N-3)*(N+1)/N* (1+2/(N-kvec-f21a-1)-r2*Es )
    } else {
      first_factor= yty/(N-3)*(1-Es*r2)*(N+1)/N
    }



    Sigmas = (matrix((N-3)/(N-1)*first_factor,rnew,nmodel,byrow=TRUE)+t(t(mat_xfxxxf)*((N-3)/(N-1)*varmult))) #by recyling this multplies each column with varmult
    Evals_minusy = t(t(mat_xfbhat)*Es)
    Eyf = as.vector(Evals_minusy %*%pmps+y.mean)
    Varyf = as.vector(Sigmas%*%pmps)*(N-1)/(N-3)

    premultfactor=yty/(N-1) #needed for hyper
    interceptfactor=(N+1)/N #needed for hyper



    ####################
    calcdensvec = function(xf_index, seqy, m_index) {


        sss=function(lbound,uboundp1,nbsteps,seqs,xf.index) {
         #simple simpson integration over shrinkage factor
         #this is just a convenience function and depends on variables in parent scope, only used in case of hyper-g
         #caution: seqs needs to have at least two elements!
         s.seq=seq(lbound,uboundp1,(uboundp1-lbound)/nbsteps)[-nbsteps]
         tmat=array(unlist(lapply(as.list(s.seq),function(ss) { dtcm(seqs,N-1,y.mean+ss*myev,premultfactor*(1-ss*myr2)*(interceptfactor+ss*myxfxxxf))})),dim=c(length(seqs),nbsteps)) #matrix of t-densities for different s
         smat=sapply(as.list(s.seq),dsgivenykernel, kpa=myk+f21a,N=N,z=myr2)  #vector of posterior densities for the different s
         if (any(is.infinite(smat))) smat[is.infinite(smat)]=0
         intconst=(4*sum(smat[c(FALSE,TRUE)])+2*sum(smat[c(TRUE,FALSE)])-3*smat[nbsteps]-smat[1])*(s.seq[nbsteps]-s.seq[1])/nbsteps/3 #calc the integration constant
         return(list(dv=c(4*tmat[,c(FALSE,TRUE)]%*%smat[c(FALSE,TRUE)]+2*tmat[,c(TRUE,FALSE)]%*%smat[c(TRUE,FALSE)]-3*tmat[,nbsteps]*smat[nbsteps]-tmat[,1]*smat[1])*(s.seq[nbsteps]-s.seq[1])/nbsteps/3, ic=intconst))
         #return the estimated density and a normalization constant
         #is done in parts because it may be recomposed later
       }



      if (any(is.na(newX[xf_index,]))) {
         densvec=numeric(0)
      }


      if (is.hyper) {
            myev=mat_xfbhat[xf_index,m_index]; myxfxxxf=mat_xfxxxf[xf_index,m_index];
            myk=kvec[[m_index]]; myr2=r2[[m_index]]



            midpoint=1-(1-Es[[m_index]])*4
           if (midpoint<.5) {
              dvl=sss(.0001,.9999999,nbsteps*2,seqy,xf_index); densvec=dvl$dv/dvl$ic
           } else {
             dvl1=sss(.0001,midpoint,nbsteps,seqy,xf_index); dvl2=sss(midpoint,1,nbsteps,seqy,xf_index)
             densvec=(dvl1$dv+dvl2$dv)/(dvl1$ic+dvl2$ic)
           }
      } else {

        densvec=dtcm(seqy,N-1,Evals_minusy[xf_index,m_index]+y.mean,Sigmas[xf_index,m_index])
      }



      return(densvec)
    }
    ##########################



    ##########################
    dens_yf = function(yfr,xf_indices=NULL) {
       if (is.null(xf_indices)) xf_indices=seq_len(rnew)
       yfdens=array(NA,dim=dim(yfr))
       for (myxf in 1:length(xf_indices)) {
        allm_dens=sapply(seq_len(nmodel), function(x) calcdensvec(xf_indices[[myxf]],yfr[myxf,],x))
        yfdens[myxf,]=as.vector(allm_dens%*%pmps)
       }
       yfdens[!is.finite(yfdens)]=NA
       if (ncol(yfdens)==1) dim(yfdens) <- NULL
       return(yfdens)
    }


    ################################################

    emptydens=list(x=numeric(0), y=numeric(0), bw=NULL, n=0, has.na=TRUE)
    class(emptydens)="density"
    dlist = lapply(vector("list",nrow(newX)),function (x)  emptydens)
    densities_calculated <- FALSE

    #if (!exists("xf_index")) xf_index=NULL

    calc_alldens = function() {
       if (densities_calculated) return(NULL)
       for (xf.index in 1:rnew) {
          if (!any(is.na(newX[xf.index,]))) {
 

 #assign("xf_index",xf.index,envir=env) # IS THIS NECESSARY?
          #xf_index <<- xf.index
          lbound=Eyf[[xf.index]]-sqrt(Varyf[[xf.index]])*4; ubound=Eyf[[xf.index]]+sqrt(Varyf[[xf.index]])*4
          seqs=seq(lbound,ubound,(ubound-lbound)/(n-1))
          allm_dens=sapply(seq_len(nmodel), function(x) calcdensvec(xf.index,seqs,x))

          myy=as.vector(tcrossprod(t(as.matrix(pmps)),allm_dens))
          mydens=list(x=seqs,y=myy,bw=NULL, n=n, has.na=FALSE)
          class(mydens)="density"
          dlist[[xf.index]] <<-  mydens
          }
       }
       densities_calculated <<- TRUE
    }
    ########################################

    consistent.yf = function(yf, xf.indices=NULL) {
       #a user check function whether a realized y conforms to requirements
       xf_series=seq_len(rnew)
       wasnull=FALSE
       if (is.null(xf.indices)) {
          wasnull=TRUE
          xf.indices=xf_series
       } else {
          if (!all(xf.indices %in% xf_series)) stop(paste("predict_index needs to be an integer between 1 and ",rnew,"!",sep=""))
       }

       if (!is.numeric(yf)) stop("realized.y must be a numeric matrix or vector!")
       if (!is.matrix(yf)) yf <- as.matrix(yf)
       if ((length(xf.indices)==1)&(nrow(yf)>1)&(ncol(yf)==1)) yf <- t(yf)
       if (nrow(newX[xf.indices,,drop=FALSE])!=nrow(yf)) {
           if (wasnull) stop(paste("realized.y must have", rnew , "elements/rows corresponding to newdata"))
           else stop("The number of rows/elements in realized.y must have the same length as predict_index!")
       }

       return(yf)

    }
    
    consistent_predict_index =function(pix)  {
      # a user check function to convert characgter predict_index into numeric
        if (is.character(pix)) {
          if (all(pix %in%  newXnames)) {
            return( match( pix, newXnames)  )  
          } else {
            stop("Forecast IDs provided in predict_index do not conform to rownames of predicted data");
          }
        } else return(pix)
    }

    ########################################



    plot.preddens = function(xf.index=1, addons="eslz", yf.addons=NULL, predict_index=NULL, addons.lwd=1.5, ...) { #, main=NULL,col="steelblue4", xlab="Response variable") {
      dotargs = match.call(expand.dots=FALSE)$...
      if (rnew>1) {
        main_default <- paste("Predictive Density Obs ", newXnames[[xf.index]], " (", nmodel, " Models)", sep="")
      } else {
        main_default <- paste("Predictive Density", " (", nmodel, " Models)",sep="")
      }
      
      dotargs=.adjustdots(dotargs,xlab="Response variable",main=main_default,col=4, zero.line=FALSE)
      
#       if (is.null(main)) {
#       if (rnew>1) {
#         main<- paste("Predictive Density Obs ", newXnames[[xf.index]], " (", nmodel, " Models)", sep="")
#       } else {
#         main<- paste("Predictive Density", " (", nmodel, " Models)",sep="")
#       }}
#       
#       plot.density(dlist[[xf.index]], main=main, xlab=xlab, col=col, zero.line=any(grep("z",addons,ignore.case=TRUE)), ...)

      thingy=dlist[[xf.index]]
      eval(as.call(c(list(as.name("plot"),as.name("thingy")),as.list(dotargs))))
      leg.col=numeric(0);leg.lty=numeric(0); leg.legend=character(0)
      
      if (any(grep("g",addons,ignore.case=TRUE))) { # grid
        graphics::grid()
      }

      if (any(grep("e",addons,ignore.case=FALSE))) {  # posterior mean
        graphics::abline(v=fit[[xf.index]],col=2,lwd=addons.lwd)
         leg.col=c(leg.col,2);leg.lty=c(leg.lty,1); leg.legend=c(leg.legend,"Exp. Value")
      }
      if (any(grep("s",addons,ignore.case=FALSE))) { # standard error bounds
        graphics::abline(v=fit[[xf.index]]-2*stderrs[[xf.index]],col=2,lty=2,lwd=addons.lwd)
        graphics::abline(v=fit[[xf.index]]+2*stderrs[[xf.index]],col=2,lty=2,lwd=addons.lwd)
         leg.col=c(leg.col,2);leg.lty=c(leg.lty,2);leg.legend=c(leg.legend,"2x Std.Errs")
      }
      
      if (any(grep("z",addons,ignore.case=TRUE))) { #zero line
        graphics::abline(h=0,col="gray",lwd=addons.lwd)
      }

      if (!is.null(yf.addons)&&is.numeric(yf.addons)) { #yf actual y realization line

         yfs=as.vector(yf.addons)
         #if (length(yfs)==rnew) {
         if (!is.na(yfs[[xf.index]])) {
           graphics::abline(v=yfs[[xf.index]],col=1,lwd=addons.lwd, lty=2)
           leg.col=c(leg.col,1);leg.lty=c(leg.lty,2);leg.legend=c(leg.legend,"Realized y")
         } else warning("yf.addons must be a vector with the same number of elements as rows in newdata!")
         #}

      }

      if (any(grep("l",addons,ignore.case=TRUE))&(length(leg.col)>0)) { #legend
        graphics::legend(x="topright",lty=leg.lty,col=leg.col,legend=leg.legend,box.lwd=0,bty="n",lwd=addons.lwd)
      }

    }

    #######################################




    fit=Eyf; names(fit)=newXnames
    stderrs=sqrt(Varyf); names(stderrs)=newXnames
    
    

    reslist= list()
    reslist$densities = function() { calc_alldens(); return(dlist) }
    reslist$fit = fit
    reslist$std.err = stderrs
    reslist$dyf = function(realized.y, predict_index=NULL) {
        predict_index=consistent_predict_index(predict_index)
        if (missing(realized.y)) { stop("You must provide a realization of the dependent variable in realized.y")}
        return(dens_yf(consistent.yf(realized.y, predict_index), predict_index))
    }

    reslist$lps = function(realized.y, predict_index=NULL) {
        predict_index=consistent_predict_index(predict_index)
        if (missing(realized.y)) { stop("You must provide a realization of the dependent variable in realized.y")}
        yf=consistent.yf(realized.y, predict_index); if (ncol(yf)!=1) stop("realized.y must have only one column!")
        yf.dens=dens_yf(yf, predict_index)
        return(-sum(log(yf.dens[!is.na(yf.dens)]))/length(yf))
    }

    reslist$plot = function(predict_index=NULL, addons="eslz", realized.y=NULL, addons.lwd=1.5, ...) {

            dotargs = match.call(expand.dots=FALSE)$...
            xf_series=seq_len(rnew)
            
            #user checks            
            predict_index=consistent_predict_index(predict_index)              
    
            if (is.null(predict_index)) { predict_index=xf_series
            } else if (!all(predict_index%in%xf_series)) stop(paste("predict_index needs to be an integer between 1 and ",rnew,"!",sep=""))
            if (!(is.null(realized.y))) {
                if (length(realized.y)!=length(predict_index)) {
                   stop("realized.y must be a vector with the same number of elements as rows in newdata (or predict_index)!")
                }
            }
            if (!is.null(realized.y)) realized.y <- consistent.yf(realized.y,predict_index);
            
      calc_alldens()
	    #temp=reslist$densities()
            oldask = graphics::par()$ask
            plotnb=0

            for (xf_index in predict_index) {
                  doplot=!dlist[[xf_index]]$has.na; plotnb=plotnb+doplot
                  if (plotnb==2) graphics::par(ask=TRUE)
                  dotargs=.adjustdots(dotargs,main=NULL,col="steelblue4", xlab="Response variable")
                  if (doplot) {
                    eval(as.call(c(list(as.name("plot.preddens"),as.name("xf_index"),addons=as.name("addons"),yf.addons=as.name("realized.y"),addons.lwd=as.name("addons.lwd")),as.list(dotargs))))
                    #plot.preddens(xf_index, addons=addons, yf.addons=realized.y, addons.lwd = addons.lwd, ..., main=main, col = col, xlab=xlab)
                  }
            }
            graphics::par(ask=oldask)
    }
    reslist$n=n
    reslist$nmodel=nmodel
    reslist$call=sys.call(0)
    class(reslist) ="pred.density"

    rm(betahat_elem1, bools, emptydens, firstbetas, firstbetas2, linvres)

    return(reslist)
}


#' Log Predictive Score
#' 
#' Computes the Log Predictive Score to evaluate a forecast based on a bma
#' object
#' 
#' The log predictive score is an indicator for the likelihood of several
#' forecasts.\cr It is defined as minus the arithmethic mean of the logarithms
#' of the point densities for \code{realized.y} given \code{newdata}.\cr Note
#' that in most cases is more efficient to first compute the predictive density
#' object via a call to \code{\link{pred.density}} and only then pass the
#' result on to \code{lps.bma}.
#' 
#' @param object an object of class \code{\link{pred.density}}, or class
#' \code{bma} (cf. \code{\link{bms}}), or class \code{\link{zlm}}
#' @param realized.y a vector with realized values of the dependent variables
#' to be plotted in addition to the predictive density, must have its length
#' conforming to \code{newdata}
#' @param newdata Needs to be provided if \code{object} is not of class
#' \code{\link{pred.density}}: a data.frame, matrix or vector containing
#' variables with which to predict.
#' @return A scalar denoting the log predictive score
#' 
#' @seealso \code{\link{pred.density}} for constructing predictive densities,
#' \code{\link{bms}} for creating \code{bma} objects, \code{\link{density.bma}}
#' for plotting coefficient densities
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'  data(datafls)
#'  mm=bms(datafls,user.int=FALSE,nmodel=100)
#'  
#'  #LPS for actual values under the used data (static forecast)
#'  lps.bma(mm, realized.y=datafls[,1] , newdata=datafls[,-1])
#'  
#'  #the same result via predicitve.density
#'  pd=pred.density(mm, newdata=datafls[,-1])
#'  lps.bma(pd,realized.y=datafls[,1])
#'  
#'  # similarly for a linear model (not BMA)
#'  zz = zlm(datafls)
#'  lps.bma(zz, realized.y=datafls[,1] , newdata=datafls[,-1])
#'  
#' @export
lps.bma <- function(object, realized.y, newdata=NULL) {
  if (!any(class(object) %in% c("pred.density","bma","zlm"))) stop("object must be of class 'pred.density', 'bma' or 'zlm'!")
  if  (any(class(object) %in% c("bma","zlm"))) {
      if (is.null(newdata)) stop("newdata must be provided if object is of class 'bma' or 'zlm'.")
      object=pred.density(object,newdata=newdata)
  }
  return(object$lps(realized.y))
}

#' @export
plot.pred.density <- function(x,  predict_index=NULL, addons="eslz", realized.y=NULL, addons.lwd=1.5, ...) {
    if (!is(x,"pred.density")) stop("x must be of class 'pred.density'!")
    x$plot(predict_index, realized.y=realized.y, addons=addons, addons.lwd=addons.lwd, ...)
}

#' @export
print.pred.density <-function(x, digits=NULL, ...) {
     outmat=matrix(numeric(0),length(x$fit),2)
     colnames(outmat)=c("Exp.Val.","Std.Err.")
     rownames(outmat)=names(x$fit)
     outmat[,1]=x$fit; outmat[,2]=x$std.err
     cat("Call:\n")
     print(x$call)
     cat(paste("\nDensities for conditional forecast(s)\n",x$n, " data points, based on ", x$nmodel, " models;\n",sep=""))
     print(outmat, digits=digits, ...)
}



#####################################
# NEW UTILITIES                ######
#####################################

#' @export
deviance.bma = function(object, exact=FALSE, ...) {
 #calculates (N-1)*posterior variance of a bma object = effective residual sum of squares
 # also works for objects of class zlm (and in principle for lm)
 #akin to method 'deviance'
 if (is.bma(object)) {
  xx=as.matrix(object$arguments$X.data);  
  ebeta = estimates.bma(object,order.by.pip=FALSE,exact=exact)[,2,drop=TRUE] 
 } else if (is(object,"lm")) {
   xx=as.matrix(object$model)
   ebeta = coef(object); if (length(ebeta)==ncol(xx) ) ebeta=ebeta[-1]
 } else stop("Required input is an object of class 'bma' or 'lm'/'zlm'.")
  
 xx =xx - matrix(colMeans(xx), nrow(xx), ncol(xx),byrow=TRUE)
 ess=as.vector(crossprod(ebeta,as.vector(crossprod(xx[,-1,drop=FALSE],xx[,1]))))
 return((as.vector(crossprod(xx[,1,drop=TRUE]))-ess))
}


#' @export
deviance.zlm=function(object, ...) deviance.bma(object)


#' @export
model.frame.bma = function(formula, ...) { 
  #akin to method 'model.frame'
  if (!is.bma(formula)) stop("argument 'formula' needs to be a bma object")
  return(as.data.frame(formula$arguments$X.data))
}

#' Variable names and design matrix
#' 
#' Simple utilities retrieving variable names and design matrix from a bma
#' object
#' 
#' All functions are \code{bma}-functions for the generic methods
#' \code{\link{variable.names}}, \code{\link{deviance}}, and
#' \code{\link{model.frame}}.
#' 
#' @aliases variable.names.bma model.frame.bma
#' @param object A \code{bma} object (as produced by \code{\link{bms}})
#' @param ... further arguments passed to or from other methods
#' 
#' @seealso \code{\link{bms}} for creating bma objects
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'  data(datafls)
#'  bma_enum=bms(datafls[1:20,1:10])
#'  
#'  model.frame(bma_enum) # similar to 
#'  bma_enum$arguments$X.data
#'  
#'  variable.names(bma_enum)[-1] # is equivalent to
#'  bma_enum$reg.names
#'  
#' @export
variable.names.bma = function(object, ...) {
  #akin to method 'variable.names'
  if (!is.bma(object)) stop("argument 'object' needs to be a bma object")
  return(c("(Intercept)", object$reg.names))
}

#' Variable names and design matrix
#' 
#' Simple utilities retrieving variable names and design matrix from a bma
#' object
#' 
#' \code{variable.names.zlm}: method \code{\link{variable.names}} for a
#' \code{\link{zlm}} model. \cr \code{vcov.zlm}: the posterior
#' variance-covariance matrix of the coefficients of a \code{\link{zlm}} model
#' - cf. \code{\link{vcov}} \cr \code{logLik.zlm}: a \code{\link{zlm}} model's
#' log-likelihood \code{p(y|M)} according to the implementation of the
#' respective coefficent prior \cr
#' 
#' @aliases variable.names.zlm vcov.zlm logLik.zlm
#' @param object A \code{bma} object (as produced by \code{\link{bms}})
#' @param ... further arguments passed to or from other methods
#' 
#' @seealso \code{\link{zlm}} for creating \code{zlm} objects
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'  data(datafls)
#'   
#'  zz=zlm(datafls)
#'  variable.names(zz)
#'  vcov(zz)
#'  logLik(zz)
#'  
#' @export
variable.names.zlm = function(object, ...) {
  #akin to method 'variable.names'
  if (!is(object,"zlm")) stop("argument 'object' needs to be zlm object")
  return(names(object$coefficients))
}

#' @export
logLik.zlm = function(object,...) {
  #marginal likelihood of a 'zlm' model, akin to method 'logLik'
  if (!is(object,"zlm")) stop("argument 'formula' needs to be zlm object")
  ret = object$marg.lik
  attr(ret, "df") = object$rank+1
  attr(ret, "nbobs") = object$rank+object$df.residual
  class(ret)="logLik"
  return(ret)
}


#' @export
vcov.zlm = function(object, include.const = FALSE, ...) {
  #akin to vcov.lm
  
  #get initial stuff
  Xmat=as.matrix(model.frame(object)[,-1,drop=FALSE]);
  if (ncol(Xmat)<1) stop("Needs at least one non-constant regressor")
  regnames= colnames(Xmat)
  Xmat=Xmat-matrix(colMeans(Xmat),nrow(Xmat),ncol(Xmat),byrow=TRUE)
  
  #complement the diagonal VCOV with off-diagonal elements, which are proportional to OLS-VCOV
  xxinv=chol2inv(chol(crossprod(Xmat)))
  outmat=((object$coef2moments[[2]]-object$coefficients[[2]]^2)/xxinv[[1]]) * xxinv
  if (include.const) {
    outmat=rbind(rep(NA,nrow(outmat)+1), cbind(rep(NA,ncol(outmat)),outmat))
    regnames=c("(Intercept)", regnames)
  }
  colnames(outmat)<-rownames(outmat) <- regnames
  return(outmat)
}


#' Posterior Variance and Deviance
#' 
#' Returns posterior residual variance, deviance, or pseudo R-squared,
#' according to the chosen prior structure
#' 
#' \code{post.var}: Posterior residual variance as according to the prior
#' definitions contained in \code{object} \cr \code{post.pr2}: A
#' pseudo-R-squared corresponding to unity minus posterior variance over
#' dependent variance. \cr \code{deviance.bma}: returns the
#' \code{\link{deviance}} of a \code{bma} model as returned from
#' \code{\link{bms}}. \cr \code{deviance.zlm}: returns the
#' \code{\link{deviance}} of a \code{\link{zlm}} model.
#' 
#' @aliases post.var post.pr2 deviance.bma deviance.zlm
#' @param object A \code{bma} object (as produced by \code{\link{bms}}) or a
#' \code{\link{zlm}} object.
#' @param exact When \code{exact=FALSE}, then \code{deviance} will be based on
#' MCMC frequencies, if \code{exact=TRUE} then it will be based on\cr
#' analytical posterior model probabilities - cf. argument \code{exact} in
#' \code{\link{coef.bma}}.
#' @param ... further arguments passed to or from other methods
#' 
#' @seealso \code{\link{bms}} for creating \code{bma} objects and priors,
#' \code{\link{zlm}} object.
#' 
#' Check \url{http://bms.zeugner.eu} for additional help.
#' @keywords utilities
#' @examples
#' 
#'  data(datafls)
#'   
#'  mm=bms(datafls[,1:10])
#'  deviance(mm)/nrow(datafls) # is equivalent to
#'  post.var(mm)
#'  
#'  post.pr2(mm) # is equivalent to
#'  1 - post.var(mm) / ( var(datafls[,1])*(1-1/nrow(datafls)) )
#'  
#' @export
post.var= function(object,exact=FALSE,...) {
  #calculates the expected posterior standard error based on effective Residual sum of squares, for objects of class 'bma', 'zlm', 'lm', ...
  if (!(is.bma(object) | is(object,"lm"))) stop("Required input is an object of class 'bma' or 'lm'/'zlm'.")
  od=deviance(object, exact=exact)  
  oy=model.frame(object)[,1,drop=TRUE]
  ret=od/length(oy)
  attr(ret,"nobs") = length(oy)
  return(ret)
}

#' @export
post.pr2= function(object,exact=FALSE) {
  #calculates a pseudo-R-squared based on effective Residual sum of squares, for objects of class 'bma', 'zlm', 'lm', ...
  if (!(is.bma(object) | is(object,"lm"))) stop("Required input is an object of class 'bma' or 'lm'/'zlm'.")
  od=deviance(object, exact=exact)  
  oy=model.frame(object)[,1,drop=TRUE]
  return(1-(od/crossprod(oy-mean(oy))[[1]]))
}
