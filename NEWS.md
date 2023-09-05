# BMS 0.3.5

Last build 2022-08-05. Rebuilding package with modern means in order to prepare expansion of package.
Coding started in Jan 2022, fixing small issues...
## Main changes
* License changed from artistic to BSD3
* Co-author added: Paul Hofmarcher, maintainer email changed.
* Source coude released at github: github.com/zeugner/BMS, no with roxygen documentation
* For issues, check bms.zeugner.eu
* 
## Bug fixes
* `bms()`: Fixed numerical tolerance bugs relating to enumeration over datasets with 1000+ observations
* Fixed coercion to logical(1) issue
* Reduced file size of bma objects somewhat


# BMS 0.3.3

Build 2013-11-22. Maintenance release to cope with changes in R version 3 compared to version 2. Alse see release notes at https://modelaveraging.wordpress.com/category/news/


# BMS 0.3.2

Build 2013-11-17. Maintenance release to cope with changes in R version 3 compared to version 2

# BMS 0.3.1 

Build 2012-08-28. Maintenance release to fix inconsistencies in package version 0.3.0 (Warnings in package checks) under new R versions >=2.14


# BMS 0.3.0 
Build 2011-05-04 

##Changes
* Changed behavior of `pmp.bma()` 
    
## New features: 
* `bms()` now can hold a set of fixed regressors (argument fixed.reg)
* New post-processing functions functions `pred.density()`, `plot.pred.density()`, `print.pred.density()`, `predict.bma()`, `predict.zlm()`
* New functions `quantile.density()`, `quantile.pred.density()`, `quantile.coef.density()`
* New function `pmpmodel()`
* New post-processing functions `post.var()` and `post.pr2()` 
* New S3 additionms `model.frame.bma()`, `model.frame.zlm()`, `deviance.bma()`, `deviance.zlm()`, `variable.names.bma()`, `variable.names.zlm()`, `logLik.zlm()`, `vcov.zlm()`
* bms can now accommodate user-defined coefficient (g) and model-priors, as well as samplers
    
## Bug fixes: 
* `bms()`: 
        * bug with accept.candi = NA fixed  
        * problems with environment in non-standard settings fixed (in internal function .construct.arglist)
        * mprior =="custompip" is now working properly (fixed in internal function .mprior.pip.init)
        * fixed problem with sys.call in non-standard calling environments
           
* `plotComp()`:
        * also works for single bma object now
        * fixed display when variable names are too long, fixed many problems with graphical parameters - completely redone
* `topmod()`: 
        * bug when likelihood==-Inf, fixed (in internal function .top10) 
        * initializing with bbeta==FALSE now works properly (fixed in .top10)
* `density.bma()`: bug with par and plot.new fixed when argument addons contains "p"
* `image.bma()`: fixed display problems when variable names too long
* `plotConv()` (resp. plot.bma): fixed warning 'In min(x) : no non-missing arguments to min; returning Inf'
* `plotModelsize()`: 
        * bugfixes for parameter ksubset and graphical parameters
        * can now cope with case when model prior is not available

## Changes to internal objects: 
* New classes: gprior, mprior, pred.density, implicit class coef.density
* Streamlining and customizability of gprior required changes to following functions:
        `bms()`, `zlm()`, .choose.gprior, gdensity, .topmod.as.bbetaT
        New functions gprior.constantg.init, gprior.eblocal.init and .gprior.hyperg.init now provide for generic handling of g-priors
        Old functions lprob.*.init are still available for compatibility
* Streamlining and customizability of mprior required changes to following functions:
        .choose.mprior, and all .mprior.*.init functions
* Streamlining and customizability of sampling functions led to streamlining of:
        .iterenum, .fls.samp, .rev.jump, .enum_startend
* Introduction of bms argument fixed.reg is handled generically by functions    
        .fixedset.mprior, .fixedset.sampler, .starter, .choose.mprior, bms
* Speeding-up of syntax (replacing $ with [[, directly calling internal functions) yields 5% speed gain in bms vs version 0.2.5
* .construct.arglist has been vastly changed



# BMS 0.2.5 

Build 2010-07-30, fixed a bug affecting with coefficient appearance 

# BMS 0.2.4

Build 2010-07-25, first release on CRAN