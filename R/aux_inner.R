###########################################
# This version: adjusted on 2011-05-05    #
###########################################
# it includes all the auxiliary functions that should only be called INSIDE the bms function



  #START: SUBFUNCTIONS ################



.ols.terms2<-function(positions,yty,k=NULL,N=N,K=K,XtX.big=XtX.big,Xty.big=Xty.big,...){
# function calculates ols terms uhat'uhat ("ymy"), beta, (X'X)^-1,...
# it's child function child.ymy and mutate are most used


   syminv <- function(symmat, ndim=ncol(symmat)) {
     #this does the same as chol2inv(chol.default(x)), but is stripped-down for speed purposes
     # Caution: symmat must always have length(symmat)>0!!!
     if (!is.matrix(symmat)) {symmat=as.matrix(symmat)}
     return( chol2inv(chol(symmat), size=ndim) )
   }

      if(is.null(k)) k=length(positions)
      XtXinv.return=numeric(0)
      # if nullmodel
      if(sum(k)==0){
        Xty=numeric(0);XtXinv=matrix(0,0,0);bhat=numeric(0);ymy=yty;positions=0;
      } else {
        XtX<-XtX.big[positions,positions,drop=FALSE]
        Xty<-Xty.big[positions]
        #do cholesky split: A=XtX=LL'
        #get lower triangular matrix from cholesky split
        XtXinv<-syminv(XtX,ndim=k)
        bhat<-crossprod(XtXinv,Xty); ymy<-yty-crossprod(Xty,bhat)[[1]]
      }

      return(list(
      full.results = function() {

        return(list(ymy=ymy, bhat=bhat, diag.inverse=XtXinv[1:k+0:(k-1)*k]))
        },

      child.ymy = function(addix=0,dropix=0,...) {
          if (!any(as.logical(c(addix,dropix)))) {return(ymy)}
          if (all(as.logical(c(addix,dropix)))) {  #swap
              jhere={1:k}[positions==dropix]; poshere=positions[-jhere];Xj=XtXinv[,jhere];Xtxi=XtX.big[poshere,addix]
              bxlessj=crossprod(XtXinv,XtX.big[positions,addix])-Xj*XtX.big[addix,dropix]; bhatx=bxlessj[-jhere]-Xj[-jhere]*bxlessj[jhere]/Xj[jhere]
              child.ymy = ymy+bhat[jhere]^2/Xj[jhere]-{Xty.big[addix]-crossprod(Xty.big[poshere],bhatx)[[1]]}^2/{XtX.big[addix,addix]-crossprod(bhatx,Xtxi)[[1]]}
              return(child.ymy)
          } else {
            if (addix==0) { #drop
              jhere={1:k}[positions==dropix]
              child.ymy=ymy+bhat[jhere]^2/XtXinv[jhere,jhere]
              return(child.ymy)
            } else { #add
              Xtxi=XtX.big[positions,addix]
              bhatx=crossprod(XtXinv,Xtxi)[,1]
              child.ymy = ymy - {Xty.big[addix]-crossprod(bhatx,Xty)[[1]]}^2 /{XtX.big[addix,addix]-crossprod(bhatx,Xtxi)[[1]]}
              return(child.ymy)
            }

          }

      },
      mutate= function(addix=0,dropix=0,newpos=numeric(0),newk=0,...) {
            #return(ols.terms2(newpos,yty,length(k),N,K=K,XtX.big=XtX.big,Xty.big=Xty.big,...,return.inverse=F))
         if (newk==0) {
             XtXinv<<-matrix(0,0,0); Xty<<-numeric(0)
         } else {
#         if (newk<7|dropix[[1]]*addix[[1]]!=0|length(c(dropix,addix))>2) { #
         if (newk<7|addix[[1]]!=0|length(c(dropix,addix))>2) { #
             Xty<<-Xty.big[newpos]
             XtXinv<<- syminv(XtX.big[newpos,newpos,drop=FALSE],ndim=newk)
         } else {

           if (dropix[1]>0) {
               jhere=sum(positions<=dropix)
               Xty<<-Xty[-jhere]
               Xj=XtXinv[,jhere];
               XtXinv<<- {XtXinv-tcrossprod(Xj/Xj[jhere],Xj)}[-jhere,-jhere]

           } else {
               jhere=sum(positions<addix)+1
               Xtxx=XtX.big[addix,newpos];Xtx=Xtxx[-jhere];Xty<<-Xty.big[newpos]
               bhatx=crossprod(XtXinv,Xtx)[,1]; bhatxadj=c(bhatx[0:(jhere-1)],-1,bhatx[jhere:k]); if (jhere==newk) bhatxadj=bhatxadj[-(jhere+1:2)]
               newinv=tcrossprod(bhatxadj,bhatxadj/(Xtxx[jhere]-crossprod(Xtx,bhatx)[[1]]))
               newinv[-jhere,-jhere]=newinv[-jhere,-jhere]+XtXinv
               XtXinv<<- newinv;
           }
         }}
#         if (any(diag(XtXinv)<0)) browser()
         positions<<-newpos; k<<-newk; bhat<<-crossprod(XtXinv,Xty)[,1]; ymy<<-yty-crossprod(Xty,bhat)[[1]]
         return(list(ymy=ymy, bhat=bhat, diag.inverse=XtXinv[1:k+0:{k-1}*k]))
      },
      return.inverse= function() XtXinv,
      ymy=ymy, bhat=bhat, diag.inverse=XtXinv[1:k+0:{k-1}*k]
      ))
 }





  




## MODEL PRIOR DEFINITIONS #####################
#the following functions define an mprior.info object
# out of the provided parameters
# mpparam: corresponds to argument 'mprior.size' in function bms   
# K: the number of covariates
# arguments passed from bms to .mprior.*.init:
# mpmode=mprior,mpparam=mprior.size,K=K,X.data=X.data,fixed.pos=fixed.pos
# Example:
# bms(attitude,mprior=.mprior.uniform.init)
# bms(attitude,mprior=.mprior.pip.init(K=ncol(attitude)-1,mpparam=seq(.1,.6,.1)))

.mprior.uniform.init = function(K,...) {
  #defines uniform model prior
    return(list(
      mp.mode="uniform", #name of the model prior
      mp.msize=K/2, #prior expected model size
      pmp=function(...) return(0), #the actual model prior function
      mp.Kdist=exp(lchoose(K,0:K)-K*log(2)) #a K+1 vector of prior model probailtiy for each model size - optional, but necessary for some processing functions
    ))
}


.mprior.fixedt.init = function(K, mpparam, ...) {
  #Defines binomial model prior ("fixed theta")
  #user checks:
  if (is.na(mpparam[1])) mpparam<-K/2
  if((mpparam[[1]]>=K)&(length(mpparam)==1)){
    warning("Submitted prior model size is >= than the nr. of   regressors\n, used K/2 instead\n\n")
    mpparam<-K/2
  }

  # actual model prior
  m=mpparam[[1]]
  return(list(
    mp.mode="fixed",
    mp.msize=m,
    pmp=function(ki,...) {
         post.odds1=ki*log(m/K)+{K-ki}*log(1-m/K)
         return(post.odds1)
    },
    mp.Kdist=stats::dbinom(x=0:K,size=K,prob=m/K,log=FALSE)
  ))
}



.mprior.randomt.init = function(K, mpparam, ...) {
  #model prior definition for beta-binmomial (random theta) model prior
  #user checks:
  if (is.na(mpparam[1])) mpparam<-K/2
  if((mpparam[[1]]>=K)&(length(mpparam)==1)){
    warning("Submitted prior model size is >= than the nr. of   regressors\n, used K/2 instead\n\n")
    mpparam<-K/2
  }
  
  # necessary initialization variables for speeding up sampling
  m=mpparam[[1]]
  vecofpriors=lgamma(1+0:K)+lgamma({K-m}/m+K-0:K)
  beta.bin=function(a=1,b=(K-m)/m,K=K,w=0:K){lgamma(a+b)-{lgamma(a)+lgamma(b)+lgamma(a+b+K)}+log(choose(K,w))+lgamma(a+w)+lgamma(b+K-w)}
  
  #actual mpinfo object
  return(list(
    mp.mode="random",
    mp.msize=m,
    pmp=function(ki,...) {
       return(vecofpriors[[ki+1]])
    },     
    mp.Kdist=exp(beta.bin(a=1,b={K-m}/m,K=K,w=0:K))
  ))

}

.getpolycoefs <- function(polyroots) {
   # helper function for mprior.pip.init
   #given the roots of a homogenous polynomial, 
   # this function finds the polynomial coefficients (ordered from highest exponent term to constant)
  if (length(polyroots)==1) return(c(1,polyroots))
  restterms=.getpolycoefs(polyroots[-1])
  c(restterms,0)+c(0,polyroots[1]*restterms)
}


.mprior.pip.init = function(K,mpparam,...) {
  #defines model prior for custom prior inclusion probabilities
  if (any(is.na(mpparam))) mpparam=rep(.5,K);
  if (!is.numeric(mpparam)) stop("For prior inclusion probabilites, you need to provide a K vector with elements between 0 and 1 for argument 'mprior.size'.")
  mpparam=as.vector(mpparam)
  if (!((length(mpparam)==K)&all(mpparam>0)&all(mpparam<=1))) stop("For prior inclusion probabilites, you need to provide a K vector with elements between 0 and 1 for argument 'mprior.size'.")
  if (any(mpparam==1L)) warning("Prior Inclsuion Prob. = 1 are impractical. Try using the argument fixed.reg")
  inclfacts=log(mpparam/(1-mpparam))

  return(list(
    mp.mode="pip",
    mp.msize=sum(mpparam),
    pmp=function(mdraw,...) {
      return(sum(inclfacts[as.logical(mdraw)]))
    },
    mp.Kdist=.getpolycoefs(mpparam/{1-mpparam})*prod(1-mpparam)
    ))
}

.mprior.customk.init = function(K, mpparam, ...) {
  # defines model prior for custom size-based model priors
  if (any(is.na(mpparam))) mpparam=rep(.5,K);
  if (!is.numeric(mpparam)) stop("For custom model size priors, you need to provide a K+1 vector with positive elements for argument 'mprior.size'.")
  mpparam=as.vector(mpparam)
  if (!((length(mpparam)==(K+1))&all(mpparam>0))) {     stop("For custom model size priors, you need to provide a K+1 vector with positive elements for argument 'mprior.size'.")     }
  mpkvec=log(mpparam)

  return(list(
    mp.mode="custom",
    mp.msize=sum(choose(K,0:K)*mpparam*{0:K})/sum(choose(K,0:K)*mpparam),      
    pmp=function(ki,...) {
      return(mpkvec[[ki+1]])
    },
    mp.Kdist=choose(K,0:K)*mpparam/sum(choose(K,0:K)*mpparam)
  ))
}





.fixedset.mprior = function(mprior.function,fullK,fixed.pos=numeric(0),K=NA,...) {
  #CAUTION: NO adjustment for mpparam yet!
  if (length(fixed.pos)==0) return(mprior.function(K=fullK,...))
  fixed.pos={1:fullK}[fixed.pos] #to convert from possible binary vector
  flexpos={1:fullK}[-fixed.pos]; flexk=length(flexpos)
  mprior=mprior.function(K=flexk,...)
  fixk=length(fixed.pos)
  mpl=list(
    mp.mode=mprior$mp.mode,
    mp.msize=mprior$mp.msize+fixk,
    pmp=function(ki,mdraw,...) {
      return(mprior$pmp(ki=ki-fixk,mdraw=mdraw[flexpos]))
    },
    mp.Kdist = c(numeric(fixk),mprior$mp.Kdist)
   )
  return(mpl)
}


.choose.mprior <- function(mpmode,mpparam,K,...,fixed.pos=numeric(0)) {
   #this function fetches an mprior.info object
  # mpmode: corresponds to argument 'mprior' in function bms
  # mpmode determines one of 5 generic model prior classes:
  # uniform, fixed equal inclusion probs, random 'theta' equal inclusion probs, 
  # custom model size priors and custom inclusion probs
  # but can also be a custom function or list such as the indivdual model prior definition functions
  # mpparam: corresponds to argument 'mprior.size' in function bms   
  # K: the number of covariates
  # fixed.pos: the positions of the variables in X the need to be kept in any model
   
   origargs=list(mpmode=mpmode,mpparam=mpparam) #read input arguments
   fixed.pos={1:K}[fixed.pos]; fixed.exist=as.logical(length(fixed.pos)); fixk=length(fixed.pos)
 
 
   
  if (!( is.character(mpmode)|| is.function(mpmode) || is.list(mpmode))) stop("'mprior' parameter must be character! (or function/list)")
  
  #switch list: choose model prior - Caution: mpparam is adjusted for fixed.set variables 
   if (is.function(mpmode) || is.list(mpmode)) { # custom model prior provided as function
     mpinfo=mpmode; 
   } else if (any(grep("fix",mpmode,ignore.case=TRUE))) { #fixed theta prior
     mpinfo=.mprior.fixedt.init     
     if (is.numeric(mpparam)) mpparam=mpparam[[1]]-fixk
  } else if (any(grep("unif",mpmode,ignore.case=TRUE))) { #uniform
     mpinfo=.mprior.uniform.init
  } else if (any(grep("custom",mpmode,ignore.case=TRUE))) { #custom model size prior 
     mpinfo=.mprior.customk.init
     if (fixed.exist && is.numeric(mpparam)) if (length(mpparam)==K+1) mpparam=mpparam[(fixk+1):length(mpparam)]
  } else if (any(grep("pip",mpmode,ignore.case=TRUE))) { #prior inclusion probabilities
     mpinfo=.mprior.pip.init  
     if (fixed.exist && is.numeric(mpparam)) if (length(mpparam)==K) mpparam=mpparam[-fixed.pos]
  } else { #random theta prior
     mpinfo=.mprior.randomt.init 
     if (is.numeric(mpparam)) mpparam=mpparam[[1]]-fixk
  } 
  
  if (is.function(mpinfo)) mpinfo=.fixedset.mprior(mpinfo,fullK=K,fixed.pos=fixed.pos,K=NA,mpparam=mpparam, mpmode=mpmode, ...)
  if (!all( c("mp.mode","mp.msize", "pmp") %in% names(mpinfo))) stop("The provided custom-built model prior is deficient.")
  if (!("origargs" %in% names(mpinfo))) mpinfo$origargs=origargs; 
  if (length(fixed.pos)>0) mpinfo$fixed.pos=fixed.pos
  class(mpinfo) <- c("mprior", class(mpinfo))
  return(mpinfo) 
}


## end: model priors




.starter=function(K,start.value,y,N=N,XtX.big=XtX.big,Xty.big=Xty.big,X=X,fixed.pos=numeric(0)){
#"starter" draws random start vector of size start.value, and keeps those regressors with a t-stat>0.2
    # in case user submitted a single number as start.value
    # we randomly draw a model with start.value regressors
    # finally only regressors with t-stats>0.2 are kept and build   the starting model
    
    #some input checks
    if (is.na(start.value[1])) {start.value=min((N-3),K)}
    if (any(start.value<0)|!(is.numeric(start.value)|is.logical(start.value))) { 
      start.value=min((N-3),K)
      warning("Argument 'start.value' did not conform to required format. start.value has been changed to default - min(N-3,K)")
    }
    if (length(start.value)==0) {start.value=numeric(K)}
    if (length(start.value)==1) {if (start.value==0) {start.value=numeric(K)}}      
    if (length(start.value)>1 && any(start.value>1)) {sv=numeric(K); sv[start.value]=1; start.value=sv; rm(sv) } 
    
    if(length(start.value)==1){
        if (start.value>min((N-3),K)){
          cat("Submitted Start value is too large, used\n min(N-3,K) as starting model size instead\n\n")
          start.value=min((N-3),K)
        }
        # draw randomly
        sorter=stats::runif(K)
        start.position=order(sorter,seq(1:K))[1:start.value]


        # calculate bhats and t-stats
        XtX.start<-XtX.big[start.position,start.position]
        XtXinv.start<-chol2inv(chol(XtX.start))
        bhat=XtXinv.start%*%Xty.big[start.position]
        e=y-X[,start.position]%*%bhat
        sse=crossprod(e)
        s2=as.numeric(sse/(N-length(start.position)))
        bcov=s2*XtXinv.start
        bt=bhat/sqrt(diag(bcov))
        # choose only regressors with t-stat>0.2
        molddraw=rep(0,K)
        goodguy=as.numeric(abs(bt)>.2)
        molddraw[start.position]=goodguy
        start.position = (1:K)[as.logical(molddraw)]
    outstart=list(molddraw=molddraw,bhat=bhat,start.position=start.position)
  }
    # else we start with the user specified starting model
    # in this case we do not test whether the t-stat's are greater than 0.2
    # but start with exactly the selected model no matter whether it is a
    # good starting model or not

    # in case you want to start with the null model
    if(length(start.value)>1 && sum(start.value)==0){
      outstart=list(molddraw=rep(0,K),bhat=rep(0,K),start.position=integer(0))
    }

    if(length(start.value)>1 && sum(start.value)>0){
          if(length(start.value)!=K){
            # in case user specified a too big model
            stop("Starting Model contains unequal to K regressors,please respecify")
          }

        start.position=which(as.logical(start.value))
        XtX.start<-XtX.big[start.position,start.position]
        XtXinv.start<-chol2inv(chol(XtX.start))
        bhat=XtXinv.start%*%Xty.big[start.position]
        molddraw=rep(0,K); molddraw[start.position]=1

      outstart=list(molddraw=molddraw,bhat=bhat,start.position=start.position)
    }
    fixed.pos=(1:K)[fixed.pos]
    if (length(fixed.pos)>0) { 
      outstart$molddraw[fixed.pos]=1 
      outstart$start.position = (1:K)[as.logical(outstart$molddraw)]
    }
    
    return(outstart)
}


###########################################################################################################################
#Sample Functions
############################################################################################################################
#First, we have implemented the original FLS Sample Function; here a variable is drawn from the set of
# K regressors, then
#conditional on whether it is included in the current model it can be discarded or added.


 .fls.samp=function(molddraw=molddraw,K=K,...,maxk=Inf,oldk=0){
 # the original FLS Sample Function; here a variable is drawn from the set of K regressors, then
 #conditional on whether it is included in the current model it can be discarded or added.


    indch<-ceiling(stats::runif(1,0,K)) #rounding to the smallest integer part by floor, uniform distr. [0,1]  
    bdropit<-as.logical(molddraw[[indch]]); 
    if (oldk==maxk) if (!bdropit) {indch=(1:K)[molddraw==1][[ceiling(stats::runif(1,0,sum(molddraw)))]]; bdropit=molddraw[[indch]]}
     if (bdropit){  #dropping
         addvar<-0;dropvar<-indch;
         molddraw[[indch]]<-0
     }
     else {                    #adding
          addvar<-indch;dropvar<-0;
          molddraw[[indch]]<-1
     }
#    addvar<-(!bdropit)*indch; dropvar<-bdropit*indch;
#    molddraw[[indch]]<- !bdropit
  
    positionnew <- {1:K}[molddraw==1]
  return(list(mnewdraw=molddraw,positionnew=positionnew,addi=addvar,dropi=dropvar))
 }
#############################################################################################################################
#Second, a reversible jump algorithm, where we
#have added a move step. See below
##############################################################################################################################
 #Reversible Jump Algorithm
 .rev.jump=function(molddraw=molddraw,K=K,...,maxk=Inf,oldk=0){
  #Reversible Jump Algorithm: with equal prob decides between swap step or add/drop step (.fls.samp)
    rev.idx=ceiling(stats::runif(1,0,2))  #rev.idx is a flag that indicates  the three possible steps of
                                   #the reversible jump algorithm, 1=birth or death and 2=move.
   # Perform Death, Birth or Move Step
   # if rev.idx is 1, do the same as in fls sampler (i.e. increase or  decrease depending
   # on variables already included in the model
   if(rev.idx==1){
     birth.death=.fls.samp(molddraw=molddraw,K=K,maxk=maxk,oldk=oldk)
     mnewdraw=birth.death[["mnewdraw"]]
     positionnew=birth.death[["positionnew"]]
     addvar=birth.death[["addi"]]; dropvar=birth.death[["dropi"]]
   }
   #move step
   if(rev.idx==2){
      var.in=(1:K)[as.logical(molddraw)]  #positions of the variables that are currently in the model
      var.out=(1:K)[!as.logical(molddraw)] #positions of the variables that are currently out of the model
      var.in.rand=ceiling(length(var.in)*stats::runif(1,0,1)); 
      addvar=var.out[ceiling(length(var.out)*stats::runif(1,0,1))]
      dropvar=var.in[var.in.rand]
      mnewdraw=molddraw; mnewdraw[addvar]=1; mnewdraw[dropvar]=0;
      positionnew=(1:K)[as.logical(mnewdraw)]
      dropvar=max(dropvar,0); addvar=max(addvar,0) # in case one of the indexes is integer(0) (borderline case)
   }
   return(list(mnewdraw=mnewdraw,positionnew=positionnew,addi=addvar,dropi=dropvar))
 }
#############################################################################################################################
#Third, there is a 'contiguity enumeration' sampler, that enumerates all possible combinations of models. 
#In particular it moves between continguous models (without repeating itself). i.e. there is always only an "add" or "drop" opration in the move to the next model.
##############################################################################################################################

.iterenum <- function(molddraw=numeric(0),K=length(molddraw),...) { 
 #Contiguity Enumeration Sampler
 # takes a binary vector (like molddraw) and iterates it such that it performs 'contiguity enumeration', 
 # i.e. enumerating all possible combinations by always changing only one entry in molddraw
         even.lead1={1:K}[!{cumsum(molddraw)%%2}]; i=even.lead1[length(even.lead1)]
        # i ist the last entry where the number of 1's up to entry i is even:
        molddraw[i]=!molddraw[i]; # then change entry i (either T->F or F->T)
        addi=molddraw[i]*i;dropi={!molddraw[i]}*i; #indch=i    
          return(list(mnewdraw=molddraw,positionnew={1:K}[as.logical(molddraw)],addi=addi,dropi=dropi))
} 

##### Enumeration in case K>N-2 ######################

.iterenum.bone = function(molddraw=numeric(0),maxk=Inf) {
# an auxiliary function for iterenum.KgtN
        even.lead1=((1:length(molddraw))[!(cumsum(molddraw)%%2)]); i=even.lead1[length(even.lead1)]
        # i ist the last entry where the number of 1's up to entry i is even:
        molddraw[i]=!molddraw[i]; # then change entry i (either T->F or F->T)
        if (sum(molddraw)>maxk) return(.iterenum.bone(molddraw,maxk)) else return(molddraw)
}

.iterenum.KgtN <- function(molddraw=numeric(0),maxk=Inf,oldk=0,...) { 
#iterenum.KgtN is slightly slower than iterenum and itererates only through models with k<maxk; else similar to iterenum
         mnewdraw=.iterenum.bone(molddraw=molddraw,maxk)
         addi=(1:length(mnewdraw))[molddraw<mnewdraw]; if (length(addi)==0) addi=0
         dropi=(1:length(mnewdraw))[molddraw>mnewdraw]; if (length(dropi)==0) dropi=0
           return(list(mnewdraw=mnewdraw,positionnew=(1:length(mnewdraw))[as.logical(mnewdraw)],addi=addi,dropi=dropi))
} 

#####  ############
#get the enumeration sampler vector for a certain iteration index
.enum_fromindex <- function(lindex){
   #lindex: an integer index; this function returns a logical vector that corresponds to the enumeration sampler draw for this index (without leading zeros)
   lindex=lindex[[1]]
   if (lindex==0) return(FALSE)
   log2=ceiling(log(lindex+1,2))
   return( as.logical((lindex+2^((log2-1):0))%/%(2^(log2:1))%% 2))
}


.enum_startend = function(iter=NA,start.value=0,K=1,maxk=K, fixed.pos=numeric(0)) {
  #does user checks and returns the starting model for enumeration and the number of iterations
  # in its most basic version, it returns list(numeric(K), 2^K-1)
  # however, it also adjusts for N-3<K and fixed variables
  # moreover, it converts to a subspace of the enumeration space if iter and start.value say so
  
  fixed.pos={1:K}[fixed.pos]
  effk=K-length(fixed.pos)
  flexpos={1:K}; if (length(fixed.pos)>0) flexpos={1:K}[-fixed.pos]
   
  start.value2=0
   
   
  if (length(start.value)==1) { 
    start.value2=suppressWarnings(as.integer(start.value)) 
    if (any(is.na(start.value2))|start.value2[[1]]<0|start.value2[[1]]>=(2^effk-1)) {
       start.value=0;start.value2=0
    } 
  } else { 
    start.value=0 
  }  
  if (length(start.value)==1) { 
     #if startvalue is an integer index satysfying above conditions then convert it into a 'draw' to start from
     start.value_cut=.enum_fromindex(start.value2)
     start.value=rep(1,K)
     start.value[flexpos]=c(numeric(effk-length(start.value_cut)),start.value_cut)   
    
  }
 
  #   do all the other enumeration stuff
  if (K>maxk) {lastindex=2^effk-1-sum(choose(effk,(maxk+1):effk))} else lastindex=2^effk-1
  if (is.na(iter)) {iter=lastindex-start.value2} #default iter is the rest to complete the enumeration from start.value
  iter=min(iter,2^effk-1-start.value2); # don't do too much!
  
  return(list(start.value=start.value, iter=iter))

} 

.fixedset.sampler = function(sampler.function, fullK, fixed.pos=numeric(0),...) {
  #converts any model sampler function into one that retains fixed variables defined in fixed.pos
  if (length(fixed.pos)==0) return(sampler.function)
  
  #initialize
  fixed.pos={1:fullK}[fixed.pos] #to convert from possible binary vector
  flexpos={1:fullK}[-fixed.pos]; flexk=length(flexpos)
  outdraw=rep(1,fullK)
  
  #now create a 'stoepsler' function that plugs in the draw of flexible varaibles into the entire set
  outfun = function(molddraw=molddraw,K=flexk,...) {
    flexdraw=sampler.function(molddraw=molddraw[flexpos],K=flexk,...)
    outdraw[flexpos]=flexdraw[["mnewdraw"]]
    addi=flexdraw[["addi"]]; dropi=flexdraw[["dropi"]];
    #indch= flexpos[flexdraw[["indch"]]]; #is this really necessary?
    if (is.numeric(addi)||is.numeric(dropi)) {
    if (addi>0) addi= flexpos[addi] else addi=0
    if (dropi>0) dropi= flexpos[dropi] else dropi=0
    }
    return(list(mnewdraw=outdraw,positionnew={1:fullK}[as.logical(outdraw)],addi=addi,dropi=dropi))
  }
  
  return(outfun)
  
}


###########################################################################################################################
#SAMPLERS WITH INTERACTION TERMS
############################################################################################################################
#initialization
.constr.intmat = function(X,K) {
    #this function identifies the columns of X named with multiple terms (sep'd by "?") as interactions
    #and constructs a matrix whose rows identify the base terms corresponding to each interaction term
   intix=grep("#",colnames(X),fixed=TRUE) # indices of columsn with interaction terms
   mPlus=diag(K)
   colnames(mPlus)<-colnames(X)
   for (jj in 1:length(intix)) {
     cix= intix[jj]
     mPlus[cix,unlist(strsplit(colnames(mPlus)[cix],"#",fixed=TRUE))]=1 # put a one in row (interaction term) and columns (all base terms)
   }
   return(mPlus)      
}


#First, we have implemented the original FLS Sample Function; here a variable is drawn from the set of
# K regressors, then
#conditional on whether it is included in the current model it can be discarded or added.


 .fls.samp.int=function(molddraw=molddraw,K=K,mPlus=mPlus,maxk=Inf,oldk=0){
      #interactions sampler for .fls.samp
      indch=ceiling(stats::runif(1,0,1)*K) #rounding to the smallest integer  part by floor, uniform distr. [0,1]
       # have to make sure that we delete all the regressors                                            
      
      if (molddraw[indch]==1){  #dropping
          mnewdraw=as.numeric(molddraw>mPlus[,indch])
          dropvar = (1:K)[xor(molddraw,mnewdraw)]; addvar=0;
      }
       else{                    #adding
          mnewdraw=as.numeric(molddraw|mPlus[indch,])
          addvar = (1:K)[xor(molddraw,mnewdraw)]; dropvar=0;
      }
      positionnew=which(mnewdraw==1)
      if (length(positionnew)>maxk) {
        return(.fls.samp.int(molddraw=molddraw,K=K,mPlus=mPlus,maxk,oldk))
      } else {
        return(list(mnewdraw=mnewdraw,positionnew=positionnew,addi=addvar,dropi=dropvar))
      }
 }
#############################################################################################################################
#Second, we have implemented a reversible jump algorithm, where we
#have added a move step. See below
##############################################################################################################################
 #Reversible Jump Algorithm
 .rev.jump.int=function(molddraw=molddraw,K=K,mPlus=mPlus,maxk=Inf,oldk=0){
       #interactions sampler for .rev.jump
    rev.idx=floor(stats::runif(1,0,1)*2)  #rev.idx is a flag that indicates  the three possible steps of
                                   #the reversible jump algorithm, 1=birth or death and 2=move.
   # Perform Death, Birth or Move Step
   # if rev.idx is 1, do the same as in fls sampler (i.e. increase or  decrease depending
   # on variables already included in the model
   if((rev.idx)|oldk==0){
     birth.death=.fls.samp.int(molddraw=molddraw,K=K,mPlus=mPlus,maxk,oldk)
     mnewdraw=birth.death$mnewdraw
     positionnew=birth.death$positionnew
     addvar=birth.death$addi; dropvar=birth.death$dropi
   } else {
      var.in=(1:K)[as.logical(molddraw)]  #positions of the variables that are currently in the model
      var.out=(1:K)[!as.logical(molddraw)] #positions of the variables that are currently out of the model
      mnewdraw=(molddraw>mPlus[,var.in[ceiling(length(var.in)*stats::runif(1,0,1))]])
      mnewdraw=mnewdraw|mPlus[var.out[ceiling(length(var.out)*stats::runif(1,0,1))],]
      positionnew=(1:K)[mnewdraw]
      addvar = (1:K)[molddraw<mnewdraw]; dropvar = (1:K)[molddraw>mnewdraw];
      if (length(dropvar)==0) dropvar=0;  if (length(addvar)==0) addvar=0
      
   }
#return(list(mnewdraw=as.numeric(mnewdraw),positionnew=positionnew,addi=addvar,dropi=dropvar,indch=rev.idx))   
   if (length(positionnew)>maxk) {
        return(.rev.jump.int(molddraw=molddraw,K=K,mPlus=mPlus,maxk,oldk))
   } else {
        return(list(mnewdraw=as.numeric(mnewdraw),positionnew=positionnew,addi=addvar,dropi=dropvar))
   }


   
 }





.post.calc <- function(gprior.info,add.otherstats,k.vec,null.count,X.data,topmods,b1mo,b2mo,iter,burn,inccount,models.visited,K,N,msize,timed,cumsumweights=NA,mcmc="bd",possign=NA) {
    #customized function for posterior results from bms
    postad.k.vec <- function(k.vec, null.count) c(null.count,k.vec) #concatenates the vector of freqencies for 1:K model sizes with freq for null model
    
    postad.gprior.info <- function(gprior.info,add.otherstats=numeric(0),cumsumweights=1) {
      # adjusts the gprior.info object resulting from choose.gprior():
      # add.otherstats: vector, cumsumweights: scalar by which otherstats are to be divided (typically nb of draws 'iter')
      if (gprior.info$return.g.stats) {
        if (length(add.otherstats)>0) {
          gprior.info$shrinkage.moments=add.otherstats/cumsumweights
        } else {
          gprior.info$shrinkage.moments=1/(1+1/gprior.info$g)
        }
      }
      return(gprior.info)
    }
    
    
    
    postad.reg.names <- function(X.data) {
      # extracts the column names of covariates or constructs ones: X.data is data.frame
      xdcn=colnames(X.data)
      if(is.null(xdcn)) {xdcn<- rep("",NCOL(X.data))}
      if(anyNA(xdcn)) {xdcn[is.na(xdcn)]<- rep("",NCOL(X.data))[is.na(xdcn)]}
      if(any(trimws(xdcn)=='')) {xdcn[trimws(xdcn)=='']<-paste0("Vbl",0:K)[trimws(xdcn)==''] }
      reg.names=xdcn[-1]
      return(reg.names)
    }


      gprior.info=postad.gprior.info(gprior.info,add.otherstats,cumsumweights)
      k.vec = postad.k.vec(k.vec,null.count)
      #if (is.na(cumsumweights)) cumsumweights=iter
      cons=.post.constant(X.data,b1mo/cumsumweights)
      pmp.10=pmp.bma(topmods,oldstyle=TRUE)
      
      if (nrow(pmp.10)==1|suppressWarnings(length(grep("error",class(try(cor(pmp.10[,1],pmp.10[,2]),silent=TRUE)))))) {
         corr.pmp=NA
      } else {
        if (var(pmp.10[,2])==0) corr.pmp=NA else corr.pmp=cor(pmp.10[,1],pmp.10[,2])
      }
      
      if (is.na(possign[[1]])) possign=numeric(K)
      
      
      info.object=list(iter=iter,burn=burn,inccount=inccount,models.visited=models.visited,b1mo=b1mo,b2mo=b2mo,
        add.otherstats=add.otherstats,cumsumweights=cumsumweights,K=K,N=N,corr.pmp=corr.pmp,msize=msize,timed=timed,k.vec=k.vec,cons=cons,pos.sign=possign)

      reg.names=postad.reg.names(X.data)
      return(list(info=info.object,k.vec=k.vec,cons=cons,gprior.info=gprior.info,pmp.10=pmp.10,reg.names=reg.names))
}
  




 
         





 ################ g prior functions #########################################

# a gprior object is a list, or a function that returns that list
# such a gprior object can be passed to bms (cf example below)
#
# EXAMPLE: a simple constant g prior list is returned by the following function:
# 
# simpleg <- function(g=NA,N=100,K=2,yty=1,...,myg=NULL) {
#  # bms/zlm calls this function with the following arguments: g, N (nb obs), K (nb of total variables), yty (TSS), null.lik (usually NA) and return.g.stats (same as bms argument g.stats), y and X
#  # additional arguments may be used by providing bms with the output of this function, rather than the function itself
#  if (is.null(myg)) myg=10 #stuff to be valid for all evaluations of the sub-functions below
#  
#  gprior.info=list( #general information for ppost-processing
#     gtype = "wappler-g", # a name for the g prior sub-type
#     is.constant = TRUE, # whether g is a predifined scalar, or changes with models (important for post-processing)
#     g=myg, #a scalar - only necessary if is.constant==TRUE
#     return.g.stats = FALSE # whether g-statistics should be collected; TRUE possible in principle, but mor bug-prone
#     )
#     
#   #  the actual gprior-computation routines are in a sub-element of the gprior object - note that myg was defined upfront   
#  gprior.info$lprobcalc = list(
#     just.loglik = function(ymy, k, ...) { 
#       #for speed reasons, this function is used when evaluationg the liklihood of a candiade model in bms
#       return(.5*{ (N-1-k)*log(1+myg)-(N-1)*log(myg*ymy + yty) } ) #the (scalar) log-likelihood of the model
#     },
#     lprob.all=function(ymy,k,bhat,diag.inverse,...) {
#       return(list(
#         lprob=.5*{ (N-1-k)*log(1+myg) - (N-1)*log(myg*ymy + yty)}, # log-likelihood
#         b1new = myg/(1+myg)*bhat, # the posterior expected coefficents
#         b2new ={(yty/myg+ymy) * (myg/(1+myg))^2 /(N-3)}*diag.inverse + (myg/(1+myg)*bhat)^2, #second posteriro moment of the coefficients
#         otherstats=numeric(0) #addtional stats to be counted when return.g.stats==TRUE
#         ))
#     }    
#   )
#   
#   return(gprior.info)
# }
# bms(attitude,g=simpleg)    
# zlm(attitude,g=simpleg)
# bms(attitude,g=simpleg(myg=20,N=nrow(attitude),K=ncol(attitude)-1,yty=var(attitude[,1])*(nrow(attitude)-1)))

.gprior.constg.init <- function(g=NA,return.g.stats=TRUE,N=N,K=K,yty=1,null.lik=NA,...) {
  gg=NULL
  if (!(is.character(g)||is.numeric(g))) g="UIP"
  
  if (any(grep("BRIC",g,ignore.case=TRUE))|any(grep("FLS",g,ignore.case=TRUE))) {
    if (N<=(K^2)){gg=(K^2)} else { gg=N }
    gtype="BRIC"
  }

  if (any(grep("RIC",g,ignore.case=TRUE))&& (!any(grep("BRIC",g,ignore.case=TRUE)))) {
    gg=(K^2)
    gtype="RIC"
  }

  if (any(grep("HQ",g,ignore.case=TRUE))|any(grep("Hannan",g,ignore.case=TRUE))) {
    gg=(log(N))^3
    gtype="Hannan-Quinn"
  }

  if(is.numeric(g)){
    gg=g
    gtype="numeric"
  }
  
  if (is.null(gg)) {
    if (!(any(grep("UIP",g,ignore.case=TRUE))|any(grep("BIC",g,ignore.case=TRUE)))) warning("The provided g prior could not be identified. Therefore the default g prior (UIP) has been selected.")
    gg=N
    gtype="UIP"
  }
 
  gprior.info=list(gtype=gtype, is.constant=TRUE, return.g.stats=return.g.stats, shrinkage.moments=gg/(gg+1), g=gg)
  g=gg
  if (!is.numeric(null.lik))  { null.lik={1-N}/2*log(yty) }
  g2=g/{g+1}; l1g=log(1+g) ; g2sq=g2^2; n1=N-1
  
  gprior.info$lprobcalc <-  list(
     #estimates the standard posterior stats under a cfixed g-prior
     #lprob.constg.init is called only once: it calculates constant terms to reduce redundancy
     #each lprob..init function has the two subfunctions:
     #   * just.loglik: this just calculates the log likelihood from given terms - as the sampler mostly needs just that
     #   * lprob.all: this calculates the log lik, as well as b1new (E(beta|Y) the expected value (normal-gamma) of coefficients), and b2new (the Expected value of coefficents squared E(beta^2|Y)= Var(beta|Y)+E(beta|Y)^2)
   
    #function(ymy, k) { gg=N; l1g=log(1+gg); n1=N-1; g2=g0/(1+g0); .5*{-k*l1g - n1*log(1-g2*(1-ymy/yty))} } -n1
    
     just.loglik=function(ymy,k,...) {
         return(.5*{{n1-k}*l1g-n1*log(g*ymy + yty)})
     },
     lprob.all=function(ymy,k,bhat,diag.inverse,...) {
         b1new = g2*bhat
        return(list(lprob=.5*{{n1-k}*l1g-n1*log(g*ymy + yty)},b1new = b1new,b2new ={{yty/g+ymy}*g2sq/{N-3}}*diag.inverse+b1new^2,otherstats=numeric(0)))
     }
     
     )
  class(gprior.info) <- c("gprior",class(gprior.info))
  return(gprior.info)
 }

 
 
.gprior.eblocal.init <- function(g=NA,return.g.stats=TRUE,N=N,K=K,yty=1,null.lik=NA,...) {
  gprior.info=list(gtype="EBL", is.constant=FALSE, return.g.stats=return.g.stats, shrinkage.moments=numeric(1), g=NA)
  
  if (!is.numeric(null.lik))  { null.lik=(1-N)/2*log(yty) }
  ymy.current=-1; k.current=-1; loglik=null.lik; Fstat=numeric(0); g2=0
  if (return.g.stats) { otherstats=numeric(1) } else { otherstats=numeric(0) }
  
  gprior.info$lprobcalc <- list(
    # estimates the local Empirical Bayes g prior as given in Liang et al (2008): "Mixtures of g Priors for Bayesian Variable Selection"; JASA
    #lprob.eblocal.init is called only once: it initializes an object that links just.loglik and lprob.all
    #each lprob..init function has the two subfunctions:
    #   * just.loglik: this just calculates the log likelihood from given terms - as the sampler mostly needs just that
    #   * lprob.all: this calculates the log lik, as well as b1new (E(beta|Y) the expected value (normal-gamma) of coefficients), and b2new (the Expected value of coefficents squared E(beta^2|Y)= Var(beta|Y)+E(beta|Y)^2)
    #           if at initalization return.g.stats=TRUE, then it additionally returns the estimated value of the shrinkage factor g2
  just.loglik=function(ymy,k,...) {
      ymy.current<<-ymy; k.current<<- k; 
      if (k==0) { return(null.lik) }
      Fstat<<-(N-k-1)/k*(yty-ymy)/ymy
      if (Fstat>1) { #g0=1/max(Fstat-1,0)
        g0<-1/(Fstat-1); g2<<-1/(g0+1)  # g0 corresponds to 1/g
      } else {
        g2<<-0; return(null.lik)
      }
      lFstat=log(Fstat); lg02=-lFstat; #Note that: g2=1-1/F, lg02=log(g0*g2)=log(1/Fstat), loglik=.5* (-k*log(Fstat)-(N-1)*log(ymy+(1/Fstat-1)*yty)-(N-1)*(log(F-1)-log(Fstat))        
      loglik<<-.5*{k*lg02-{N-1}*{log(ymy + g0*yty) + {log(Fstat-1)-lFstat }}}
      return(loglik)
  },
  
  lprob.all=function(ymy,k,bhat,diag.inverse,...) {
      if (k==0) { return(list(lprob=null.lik, b1new=numeric(0), b2new=numeric(0),otherstats=otherstats)) }    
      if ((ymy!=ymy.current)|(k!=k.current)) { #if just.loglik was already called just before with the same parameters, no need to calculate the F-stat stuff again
        Fstat<<- {N-k-1}/k*(yty-ymy)/ymy
        if (Fstat>1) { 
          g0=1/{Fstat-1}; g2<<-1/(g0+1)
          lFstat=log(Fstat); lg02=-lFstat; 
          loglik<<-.5*{k*lg02-{N-1}*{log(ymy + g0*yty) + {log(Fstat-1)-lFstat }}}          
        } else { 
          g0=0; g2<<-0;
          loglik<<-null.lik
        }
      }       
      if (return.g.stats) { otherstats=g2 }
      b1new = g2*bhat
      if (g2>0) { b2new ={{(1/g2-1)*yty+ymy}*{g2^2}/{N-3}}*diag.inverse+b1new^2 } else {b2new = numeric(k)}
      return(list(lprob=loglik,b1new = b1new,b2new=b2new,otherstats=otherstats)) 
  }
  
  )
  class(gprior.info) <- c("gprior",class(gprior.info))
  return(gprior.info)
}
    
    
.gprior.hyperg.init <- function(g=NA,return.g.stats=TRUE,N=N,K=K,yty=1,null.lik=NA,...) {
  #user checks
  if (!is.character(g)) g="hyper"
  
  if (any(grep("=",g))) { 
    f21a=suppressWarnings(as.numeric(unlist(strsplit(g,"="))[2]))
    if (!is.numeric(f21a)|is.na(f21a)) {
      f21a.char=suppressWarnings(as.character(unlist(strsplit(g,"="))[2]))
      if (any(grep("bric",f21a.char,ignore.case=TRUE))) {
          f21a=2+2/max(N,K^2)
      } else if (any(grep("uip",f21a.char,ignore.case=TRUE))) {
          f21a=2+2/N
      } else {
          warning("You did not supply a proper 'a' parameter for the hyper g prior (like e.g. the format g='hyperg=3.1' or g='hyper=UIP') - thus set to default value 'hyper=UIP' instead.")
          f21a=2+2/N
      }
    } else {
      if (f21a<=2|f21a>4) {
        f21a=2+2/N
        warning("You provided an 'a' parameter for the hyper g prior that is not element of (2,4]. I chose the default value 'hyper=UIP' instead.")
      }
    }
  } else {
    f21a=2+2/N
  }
   gprior.info=list(gtype="hyper", is.constant=FALSE, return.g.stats=return.g.stats, shrinkage.moments=numeric(2),g=NA, hyper.parameter=f21a)
  
  
  
  #estimates the hyper g prior as given in Liang et al (2008): "Mixtures of g Priors for Bayesian Variable Selection"; JASA
  #initializing an object that links just.loglik and lprob.all
  #   argument f21a corresponds to the hyper-parameter "a" in Liang et al. (2008): any value in ]2,4], default 3 
  #each lprob..init function has the two subfunctions:
  #   * just.loglik: this just calculates the log likelihood from given terms - as the sampler mostly needs just that
  #   * lprob.all: this calculates the log lik, as well as b1new (E(beta|Y) the expected value (normal-gamma) of coefficients), and b2new (the Expected value of coefficents squared E(beta^2|Y)= Var(beta|Y)+E(beta|Y)^2)
  
    
  #initalizing global values
  if (!is.numeric(null.lik))  { null.lik={1-N}/2*log(yty) }
  gmoments=numeric(2)
  N12={N-1}/2; la2=log(f21a-2)
  log.lik = null.lik; ymy.current=-1; k.current=-1; intconstinv=f21a-2; #intconstinv is the inverse of 1/(a-2) times the integration constant, i.e. (k+a-2)/ 2F1( (N-1)/2,1, (k+a)/2, R^2 )
  
  #initalizing the 2F1 hypergeometric function object
  f21o=.f21_4hyperg(N,K,f21a)
  
  gprior.info$lprobcalc <- list(
  
  just.loglik=function(ymy,k,...) {
      if (k==0) {return(null.lik)}
      ymy.current <<- ymy; k.current <<- k
      intconstinv <<- {k+f21a-2}/f21o[["calcit"]](1-ymy/yty,k)
      if (intconstinv<0) {intconstinv <<- k+f21a-2} #this may happen because of numerical inaccuracies: e.g. ymy marginally greater than yty
      log.lik <<- null.lik + la2 - log(intconstinv)
      return(log.lik)
  },
  
  lprob.all=function(ymy,k,bhat,diag.inverse,...) {
      if (k==0) {return(list(lprob=null.lik,b1new=numeric(0), b2new=numeric(0),otherstats=c(2/f21a,8/f21a/(f21a+2))))}
      N3=N-3; ka2=k+f21a-2; R2=1-ymy/yty; #collect terms
      
      if ((ymy!=ymy.current)|(k!=k.current)) { #if just.loglik was already called just before with the same parameters, no need to calculate the F-stat stuff again
        intconstinv <<- ka2/f21o[["calcit"]](R2,k)
        log.lik <<- null.lik + la2 - log(intconstinv)
      }
      g2hyper= {intconstinv-ka2+N3*R2}/{R2*{N3-ka2}} # E(g/(1+g)|Y)
      gbetavar = {{1+2/N3*R2/{1-R2}}*intconstinv+{N3-2}*R2-ka2} *N3*{1-R2}/{N3-ka2}/{N3-ka2-2}/R2 * yty/N3*diag.inverse #Cov(beta|Y)
      
      if (return.g.stats) {
          ka=ka2+2;
          Eg22= { {{N3-2}*R2-ka}*intconstinv + {N3*R2-ka2}^2-2*{N3*R2^2-ka2} }/R2^2/{N3-ka2}/{N3-ka}
          gmoments=c(g2hyper,Eg22)
      }
                
      return(list(lprob=log.lik, b1new = g2hyper*bhat ,b2new =gbetavar+g2hyper^2*bhat^2, otherstats=gmoments)) 
  }
  
  )
  
  class(gprior.info) <- c("gprior",class(gprior.info))
  return(gprior.info)
  
}



### used for gprior.hyperg.init #####
.f21_4hyperg=function(N,K,f21a,ltermbounds=c(200,600,1400,3000)) {
  # this function calculates the value of a Gaussian hypergeometric function 2F1((N-1)/2,1,(f21a+k)/2,z)
  # as given in Liang et al (2008): "Mixtures of g Priors for Bayesian Variable Selection"; JASA, formula (18)
  # first initialize the object by calling myobject = .f21_4hyoperg() : this calcs recycable terms to reduce redundancy
  # 
  # then evaluate by myobject$calcit(z,k)

  #Initialization:
  create.lterms=function(cc) lapply(as.list(ltermbounds),function(x) (a+0:x)/(cc+0:x)) # a subfunction
  a=(N-1)/2; cv=(f21a+0:K)/2 #create the a and c terms for the (K+1) different model sizes
  lterms=(lapply(cv,create.lterms)) #create the Pochhammer terms in vectors of different lengths, lterms is a list: each element has a list of which each element contains Pochhammer terms up to the index specified in ltermbounds (see function argument)
  ltermbounds=c(0,ltermbounds[-length(ltermbounds)]) #this is just used for deciding which vector length to take
  
 return(list(
  calcit=function(z,k) {
    # this function actually calculates 2F1((N-1)/2,1,(f21a+k)/2,z) in its power series formulation
    nbterms=sum(ceiling(abs({a-cv[k+1]}/{1-z})*1.5)>=ltermbounds) # l=(a-c)/(1-z) is the index from which on the summands decline 
                                                                 # this times 1.5 decides which lvector length (lbounds) to take
    return(sum(cumprod(z*lterms[[k+1]][[nbterms]]))+1) #calculate it
  }                   
 ))
}

.f21simple=function(a,c,z) {
# this function calculates the value of the Gaussian hypergeometric function 2F1(a,1,c,z), where 0=<z=<1
# this function is a simple wrapper destined for user application
  f21o=.f21_4hyperg(2*a+1,0,c*2)
  f21o$calcit(z,0)
}






.choose.gprior <- function(g,N,K,return.g.stats=FALSE,yty=N,...) {
  #chooses the g-prior subject to the two parameters given from bms() user input
  #returns a list with gtype, is.constant, g, as well as prior-specific elements
  # use benchmark uip criterion as standard setting
  # g corresponds to argument "g" in function bms, return.g.stats to "g.stats"
  # N is number of obs, K max number of regressors
  #
  # returns a list conforming to the gprior.info in function 'bms'

  # first check for user-defined g-priors
  if (is.list(g)) { 
    if (!all(c("gtype","is.constant","return.g.stats", "lprobcalc") %in% names(g))) stop("The provided g-prior list (in argument 'g') does not conform to the standards of a g-prior object.")  
    if (!("g" %in% names(g))) {g$is.constant = FALSE; g$g=NA}
    if (!("shrinkage.moments" %in% names(g))) {g$shrinkage.moments=ifelse((g$is.constant&&is.numeric(g)), g/(1+g), 0)}
    if (!all(sapply(g$lprobcalc,is.function))) stop("The slot 'lprobcalc' in the provided g-prior list (in argument 'g') does not conform to the standards of a g-prior object.")
    return(g)
  }
  
  
  if (is.function(g)) {
    return(g(g=g,return.g.stats=return.g.stats,N=N,K=K,yty=yty,...))
  }


  # now guess a standard g-prior from user's arguments
  
  if(is.numeric(g)){
      return(.gprior.constg.init(g=g,return.g.stats=return.g.stats,N=N,K=K,yty=yty))
  }

  # if user specifies empirical bayes estimation (local)
  if (any(grep("EBL",g,ignore.case=TRUE))) {
     return(.gprior.eblocal.init(g=g,return.g.stats=return.g.stats,N=N,K=K,yty=yty))
  }

  # if user specifies hyper-prior (local)
  if (any(grep("hyper",g,ignore.case=TRUE))) {   
       return(.gprior.hyperg.init(g=g,return.g.stats=return.g.stats,N=N,K=K,yty=yty))
  }

  return(.gprior.constg.init(g=g,return.g.stats=return.g.stats,N=N,K=K,yty=yty))
  
}

# end g-Stuff ##########################################




# FUNCTIONS DESIGNED FOR BEING CALLED IN bms() #######################################
# these functions are only subfunction to be called inside bms()
# [therefore they necessitate the statement environemnt(SUBFUNCTION) <- environment() inside bms() ]
# the are defined outside of bms() for readability and modularity purposes


# compatibility with package version 0.2.5
.lprob.constg.init = function(...) {
 gpo=.gprior.constg.init(...)
 return(gpo$lprobcalc)
}

.lprob.eblocal.init = function(...) {
 gpo=.gprior.eblocal.init(...)
 return(gpo$lprobcalc)
}


.lprob.hyperg.init = function(...) {
 gpo=.gprior.hyperg.init(...)
 return(gpo$lprobcalc)
}
