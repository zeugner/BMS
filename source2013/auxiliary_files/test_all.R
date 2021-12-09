# This file is meant to test all instances of the 'Fls' bma code

rm(list=ls(all.names=TRUE)) #delete all

library(BMS)
data(datafls)

##### BASIC TESTS ##########################

#run the sampling code
 m_default =bms(X.data=datafls,burn=1000,iter=9000,user.int=T)

# some results:
 estimates.bma(m_default,exact=T,std.coefs=T) #standardized coefficeint estimates based on top 100 models
 pmp.bma(m_default) #post. model probs. for top 100 models based on MCMC freqs and likelihoods
 plotConv(m_default[1:15]) #plot post model probs for best 15 models
 m_default[1]$topmod #best model
 density(m_default,reg=1) #plot coef. density for regressor 1
 oo=density(m_default,plot=FALSE,std.coefs=TRUE)
 info.bma(m_default) #some stats on m_default
 image(m_default)
 plot(m_default)

#### MORE DETAILED TESTS ##################

options(error=recover)

test.post.cmds = function(model4test) {
 #post level
 test= pmp.bma(model4test) #
 test= topmodels.bma(model4test)
 test= estimates.bma(model4test,exact=F,order.by.pip=T,include.constant=T)
 test= estimates.bma(model4test,exact=T,order.by.pip=F,include.constant=F)
 test= estimates.bma(model4test,exact=F,order.by.pip=T,include.constant=T,std.coefs=T,condi.coef=T,incl.possign=F)
 test= info.bma(model4test)
 test= print(model4test)


 #plots
 cat("plotModelsize()\n")
 plotModelsize(model4test)
 cat("plotModelsize(,exact=T)\n")
 plotModelsize(model4test,exact=T)
 cat(" density.bma(,reg=best)\n")
 oo=density.bma(model4test,reg=rownames(estimates.bma(model4test))[1]);
 plot(oo)
 cat("plotConv()\n")
 plotConv(model4test[1:10])
 cat("image()\n")
 image(model4test,yprop2pip=F)
 }
 

 #testing default routine
 m_default =bms(X.data=datafls,burn=1000,iter=9000,user.int=T)
 test.post.cmds(m_default)
 
 m_default2 =bms(X.data=datafls,burn=1000,iter=9000,start.value=6,user.int=F)
 test.post.cmds(m_default2)
 
 is.bma(combine_chains(m_default,m_default2)) #check whether combining chains works
 plotComp(m_default,m_default2) #compare results
 
 m_fls = bms(X.data=datafls,burn=10000,iter=90000,nmodel=2000,g="BRIC",mcmc="bd",mprior="uniform",user.int=F)
 estimates.bma(m_fls,exact=T) #similar to what FLS (2001) did

# check different settings for g 
 m_g1 = bms(X.data=datafls,burn=1000,iter=9000,g=nrow(datafls))
 m_g2 = bms(X.data=datafls,burn=1000,iter=9000,g="EBL",g.stats=T)
 m_g3 = bms(X.data=datafls,burn=1000,iter=9000,g="hyper=2.01",g.stats=T)
 m_g4 = bms(X.data=datafls,burn=1000,iter=9000,g="hyper=BRIC",g.stats=T) 
 test.post.cmds(m_g3)

 #check differnt starting values
 m_start0 = bms(X.data=datafls,burn=1000,iter=9000,g="hyper",start.value=0)
 m_start1 = bms(X.data=datafls,burn=1000,iter=9000,g="hyper",start.value=NA)
 m_start2 = bms(X.data=datafls,burn=1000,iter=9000,g="hyper",start.value=300)
 test.post.cmds(m_start0)
 
 #check different model priors
 m_mprior0 = bms(X.data=datafls,burn=1000,iter=9000,g="EBL",mprior="fixed",mprior.size=4) #similar to Doppelhofer, Miller Sala-i-Martin
 m_mprior1 = bms(X.data=datafls,burn=1000,iter=9000,g="EBL",mprior="random") #similar to Ley and Steel (2007)
 m_mprior2 = bms(X.data=datafls,burn=1000,iter=9000,g="EBL",mprior="random",mprior.size=ncol(datafls)/4)
 m_mprior3 = bms(X.data=datafls,burn=1000,iter=9000,g="hyper",mprior="uniform") #uniform model prior
 m_mprior4 = bms(X.data=datafls,burn=1000,iter=9000,g="UIP",mprior="customk",mprior.size=(rep(1,ncol(datafls)))) #custom model size prior (actually uniform)
 m_mprior5 = bms(X.data=datafls,burn=1000,iter=9000,g="UIP",mprior="pip",mprior.size=c(0.99,rep(0.5,ncol(datafls)-2))) #custom prior inclusion probabilities (actually uniform) 
 
 test.post.cmds(m_mprior2)
 
 
 m_mcmc0 = bms(X.data=datafls,burn=1000,iter=9000,mcmc="rev.jump") #reversible jump MC3 sampler
 m_mcmc1 = bms(X.data=datafls[,1:12],burn=1000,iter=9000,mcmc="enumerate") #total enumeration instead of MC3 sampling
 
 m_nmodel0 = bms(X.data=datafls[,1:12],burn=1000,iter=9000,g="EBL",nmodel=0) #saving no top models



 
  m_logfile = bms(X.data=datafls,burn=1000,iter=9000,logfile="",logstep=500) #print out logfile to console
 
  m_all0 = bms(X.data=attitude,burn=1000,iter=9000,g="EBL",mprior="uniform", mprior.size=3,mcmc="bd",nmodel=300,user.int=T,g.stats=T) #try a lot
  m_all1 = bms(X.data=datafls[,1:9],burn=1000,iter=9000,g="hyper=2.5",mprior="uniform", mprior.size=7,mcmc="enumerate",nmodel=300,user.int=T,g.stats=T)

 m_default3 =bms(X.data=datafls[1:35,],burn=1000,iter=9000,start.value=28,user.int=F) #try out with more K than N
 test.post.cmds(m_default3)


 ##INTERACTION SAMPLER 
 #create a dataset with two interaction terms
 dataint=datafls
 dataint=cbind(datafls,datafls$LifeExp*datafls$Abslat/1000,datafls$Protestants*datafls$Brit-datafls$Muslim) 
 names(dataint)[ncol(dataint)-1]="LifeExp#Abslat"
 names(dataint)[ncol(dataint)]="Protestants#Brit#Muslim"

 #check interaction sampler
 m_int0 = bms(X.data=dataint,burn=1000,iter=9000,start.value=0) 
 m_int1 = bms(X.data=dataint,burn=1000,iter=9000,start.value=(ncol(dataint)-1),g="hyper",mprior="random",mcmc="bd.int", nmodel=400)
 m_int2 = bms(X.data=dataint[20:39,],burn=1000,iter=9000,start.value=0,mprior="fixed",mcmc="bd.int",mprior.size=4)  #case K>N
 rm(dataint) 

 
 
##Tests on ENUM & combine
  menum=bms(X.data=datafls[,1:12],mcmc="enum",nmodel=2^11)  #standard enum

  #menum a should be exactly the same as menum
  menum1=bms(X.data=datafls[,1:12],mcmc="enum",iter=1000,nmodel=1020)
  menum2=bms(X.data=datafls[,1:12],mcmc="enum",start.value=1001,nmodel=1020)  
  menuma=combine_chains(menum1,menum2)

  if (any(abs(estimates.bma(menum,order.by.pip=F)-estimates.bma(menuma,order.by.pip=F))>2E-8)) stop("Problem A")
  if (any(abs(estimates.bma(menum,order.by.pip=F,exact=T)-estimates.bma(menuma,order.by.pip=F,exact=T))>2E-8)) stop("Problem B")  
 

## additional functions
  if (bin2hex(hex2bin("fa04"))!="fa04") stop("Problem with bin2hex")
  fullmodel.ssq(m_default)
  f21hyper(40,1,5,.9) #should be 1.853011e+31
  print(m_default)
  print(m_default[1])
  print(m_default[1]$topmod)
  print(m_default[2:10]$topmod)
 




 
 options(error=NULL)
 
 