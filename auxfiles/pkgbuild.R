rmpkg=function(pname) { try(detach(paste0("package:",pname)));try(unloadNamespace(pname));try(remove.packages(pname)); }
rmpkg('BMS')


library(devtools)
setwd("C:/Users/zeugnst/rpackages/BMS/BMS/BMS")

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
devtools::install()
devtools::build()
