#--------------------------------------------#
RLIBS<-"/home/genouest/irset/archives/softs/R/3.2.3/libs"
setwd("/home/genouest/irset/archives/projects/ChemPSy")
.libPaths(RLIBS)
#--------------------------------------------#

#--------------------------------------------#
# Custom CDF file linked to Entrez Gene IDs must be downloaded from the BRAINARRAY website: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
CDF.File  <- commandArgs(TRUE)[1]
#CDF.File <- "/home/genouest/irset/archives/projects/ChemPSy/scripts/CDF/Rat2302_Rn_ENTREZG.cdf"
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
if("makecdfenv" %in% rownames(installed.packages()) == FALSE) { biocLite("makecdfenv")}
library(makecdfenv)
#--------------------------------------------#

#--------------------------------------------#
CDF.Env <- make.cdf.env(basename(CDF.File),cdf.path=dirname(CDF.File))
#--------------------------------------------#


#--------------------------------------------#
save(CDF.Env,file=sprintf("%s.RData",CDF.File))
#--------------------------------------------#

#--------------------------------------------#
quit()