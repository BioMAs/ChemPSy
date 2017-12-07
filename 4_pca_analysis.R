#--------------------------------------------#
# DrugMatrix downloaded at ftp://anonftp.niehs.nih.gov/drugmatrix
# folder organization is similar
#--------------------------------------------#

#--------------------------------------------#
RLIBS<-"/home/genouest/irset/archives/softs/R/3.2.3/libs"
setwd("/home/genouest/irset/archives/projects/ChemPSy")
.libPaths(RLIBS)
#--------------------------------------------#


#--------------------------------------------#
RData.File <- commandArgs(TRUE)[1]
Output.Dir <- commandArgs(TRUE)[2]

#RData.File       <- "processed_data/Rattus_norvegicus/KIDNEY/3_reduce_data/data.reduced.RData"
#RData.File       <- "test/data.reduced.RData"
#Output.Dir       <- "test"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
if("FactoMineR" %in% rownames(installed.packages()) == FALSE) {install.packages("FactoMineR", lib=RLIBS)}
library(FactoMineR)
load(RData.File)
#--------------------------------------------#

#--------------------------------------------#
#Data.bin.ncp.res = estim_ncp(t(Data.bin),scale=FALSE)
#Data.bin.pca.res = PCA(t(Data.bin), scale.unit=FALSE, ncp=Data.bin.ncp.res$ncp, graph=F)

#Data.log2fc.ncp.res = estim_ncp(t(Data.log2fc),scale=FALSE)
#Data.log2fc.pca.res = PCA(t(Data.log2fc), scale.unit=FALSE, ncp=Data.log2fc.ncp.res$ncp, graph=F)

Reduced.Data.bin.ncp.res = estim_ncp(t(Reduced.Data.bin),scale=FALSE)
Reduced.Data.bin.pca.res = PCA(t(Reduced.Data.bin), scale.unit=FALSE, ncp=Reduced.Data.bin.ncp.res$ncp, graph=F)

Reduced.Data.log2fc.ncp.res = estim_ncp(t(Reduced.Data.log2fc),scale=FALSE)
Reduced.Data.log2fc.pca.res = PCA(t(Reduced.Data.log2fc), scale.unit=FALSE, ncp=Reduced.Data.log2fc.ncp.res$ncp, graph=F)
#--------------------------------------------#

#--------------------------------------------#
InformationPlot <- function(ncp,PDF.File) {
  pdf(file=PDF.File, width=7, height=7, title="Information plot")
  plot(1:length(ncp$criterion),ncp$criterion,xlab="Components",ylab="Mean error")
  lines(1:length(ncp$criterion),ncp$criterion)
  lines(c(ncp$ncp,ncp$ncp),c(0,ncp$criterion[ncp$ncp]),col="red")
  lines(c(1,ncp$ncp),c(ncp$criterion[ncp$ncp],ncp$criterion[ncp$ncp]),col="red")
  text(ncp$ncp,ncp$criterion[ncp$ncp],labels=ncp$ncp,pos=3)
  dev.off()
  return(PDF.File)
}
#--------------------------------------------#

#--------------------------------------------#
#InformationPlot(Data.bin.ncp.res    , sprintf("%s/PCA.Data.bin.InformationPlot.pdf"   ,Output.Dir))
#InformationPlot(Data.log2fc.ncp.res , sprintf("%s/PCA.Data.log2fc.InformationPlot.pdf",Output.Dir))
InformationPlot(Reduced.Data.bin.ncp.res    , sprintf("%s/PCA.reduced.Data.bin.InformationPlot.pdf"   ,Output.Dir))
InformationPlot(Reduced.Data.log2fc.ncp.res , sprintf("%s/PCA.reduced.Data.log2fc.InformationPlot.pdf",Output.Dir))
#--------------------------------------------#

#--------------------------------------------#
save(list=ls()[grep("^Reduced.|^Data.",ls())], file = sprintf("%s/PCA.RData",Output.Dir))
#--------------------------------------------#

quit()
