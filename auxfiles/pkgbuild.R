rmpkg=function(pname,silent=TRUE) { if (any(grepl(paste0(pname,"$"),searchpaths()))) try(detach(paste0("package:",pname),character.only = TRUE),silent=silent);try(unloadNamespace(pname),silent=silent);try(remove.packages(pname),silent=silent); }
rmpkg('BMS')

s3m=list(stats=c("density", "quantile", "variable.names", "vcov", "predict", "deviance", "coef", "model.frame"))

library(devtools)
if (!grepl('err',class(try(utils:::win.version(),silent=TRUE)))) {
   setwd("C:/Users/zeugnst/rpackages/BMS/BMS/BMS")
} else {
  setwd("/home/admin-b1/test/BMS")
}

.flsresults= readRDS("auxfiles/flsresultsvignette.rds")
usethis::use_data(.flsresults,overwrite = TRUE, internal = TRUE)
dataFrame=read.table(paste("source2013/auxiliary_files/fls_data.txt",sep=""),header=T)
remofls=c("Country","pop50","Latitude","Longitude","Port", "OECD")
datafls=dataFrame[,-charmatch(remofls,colnames(dataFrame))];
ctryiso=c("DZ", "AR", "AU", "AT", "BE", "BO", "BW", "BR", "CM", "CA", "CL", "CO", "CG", "CR", "CY", "DK", "DO", "EC", "SV", "ET", "FI", "FR", "DE", "GH", "GR", "GT", "HT", "HN", "HK", "IN", "IE", "IL", "IT", "JM", "JP", "JO", "KE", "KR", "MG", "MW", "MY", "MX", "MA", "NL", "NI", "NG", "NO", "PK", "PA", "PY", "PE", "PH", "PT", "SN", "SG", "ES", "LK", "SE", "CH", "TW", "TZ", "TH", "TN", "TR", "UG", "UK", "US", "UY", "VE", "ZR", "ZM", "ZW")
rownames(datafls)=ctryiso
rm(remofls); rm(dataFrame); rm(ctryiso)
usethis::use_data(datafls,overwrite = TRUE)


devtools::document()
devtools::load_all(export_all = FALSE)
devtools::document()
PNS=readLines("NAMESPACE")
cat(c(PNS[1],"",paste0("importFrom(stats,",s3m$stats,")"),
      'importFrom("methods", "is", "new")',
      'importFrom("stats", "as.formula", "cor", "na.omit", "var")',
      PNS[-1]),sep = '\n', file="NAMESPACE")
rm(s3m,PNS)

devtools::install()


hh=Sys.getenv('PATH')
if (!grepl("TinyTeX",hh)) {
Sys.setenv(PATH=paste0(hh,";C:\\ProgramData\\Microsoft\\AppV\\Client\\Integration\\E5372EBE-A735-49DA-9EB2-A1B54EE9FBEB\\Root\\RNonStandardPkgs\\TEX\\TinyTeX\\bin\\win32"))
}
rm(hh)


if (!grepl('err',class(try(utils:::win.version(),silent=TRUE)))) {
  devtools::build(vignettes = F)
 
} else {
  devtools::build()
 }




#%ALLUSERSPROFILE%\Microsoft\AppV\Client\Integration\602C1D12-9652-4057-A07E-CB0228E7F453\Root\R\R-4.1.2\bin\x64\R.exe CMD check --as-cran C:\Users\zeugnst\rpackages\BMS\BMS\BMS_0.3.5.tar.gz
#shell("%ALLUSERSPROFILE%\Microsoft\AppV\Client\Integration\602C1D12-9652-4057-A07E-CB0228E7F453\Root\R\R-4.1.2\bin\x64\R.exe CMD check --as-cran C:\Users\zeugnst\rpackages\BMS\BMS\BMS_0.3.5.tar.gz")


# importFrom(stats,density)
# importFrom(stats,quantile)
# importFrom(stats,variable.names)
# importFrom(stats,vcov)
# importFrom(stats,predict)
# importFrom(stats,deviance)
# importFrom(stats,coef)
# importFrom(stats,model.frame)
# 

