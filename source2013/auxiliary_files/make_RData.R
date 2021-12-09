
# for making an RData out of the BMS package
#!!! THIS ONLY WORKS ON VERSIONS 2.5 to 2.9 !!!

#### LOAD NECESSARY DATA #########################
rm(list=ls(all.names=T))
setwd("C:/Dokumente und Einstellungen/zeugner/Eigene Dateien/ULB/bmastuff/current_toolbox")

setClass("topmod",representation(addmodel="function", lik="function", bool="function", ncount="function", nbmodels="function", nregs="function", betas_raw="function", betas2_raw="function", kvec_raw="function", bool_binary="function", betas="function", betas2="function", fixed_vector="function"))
setClass("bma",representation(info="list",arguments="list",topmod="topmod",start.pos="integer",gprior.info="list",mprior.info="list",X.data="data.frame",reg.names="character",bms.call="call"))
setClass("zlm",representation(coefficients="numeric",residuals="numeric",rank="numeric",fitted.values="numeric",df.residual="numeric",call="call",terms="formula",model="data.frame",coef2moments="numeric",marg.lik="numeric",gprior.info="list"))
source("aux_outer.R")
source("aux_inner.R")
source("bma.R")
dataFrame=read.table("fls_data.txt",header=T)
remofls=c("Country","pop50","Latitude","Longitude","Port", "OECD")
datafls=dataFrame[,-charmatch(remofls,colnames(dataFrame))]
rm(dataFrame); rm(remofls);

save(list = ls(all=TRUE), file = "../package/BMS.RData")


