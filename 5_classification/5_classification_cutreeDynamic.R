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
PCA.RData.File <- commandArgs(TRUE)[1]
Output.Dir <- commandArgs(TRUE)[2]

#PCA.RData.File       <- "test/PCA.RData"
#Output.Dir       <- "test"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
if("dynamicTreeCut" %in% rownames(installed.packages()) == FALSE) {install.packages("dynamicTreeCut", lib=RLIBS)}
if("flashClust" %in% rownames(installed.packages()) == FALSE) {install.packages("flashClust", lib=RLIBS)}
library(dynamicTreeCut)
library(flashClust)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
FUNCTION <- function(Data,DistMethod,ClustMethod,dSValue,Prefix) {
	 if (ClustMethod == "tree" && dSValue > 1) {return ()}

	 cat(Prefix,"is running...\n")
	 dist.res <- dist(Data,method=DistMethod)
	 hclu.res <- hclust(dist.res, method='ward')
	 
   	 cutree.res <- cutreeDynamic(hclu.res,method=ClustMethod,distM=as.matrix(dist.res),deepSplit=dSValue)
  	 K <- length(unique(cutree.res))
	 k <- 0	 
	 for (ki in unique(cutree.res)) {
	     k <- k+1
   	     Individuals <- rownames(Data)[which(cutree.res == ki)]
	     write.table(matrix(Individuals,ncol=1),file = sprintf("%s_C%s.%s.txt",Prefix,K,k), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
	 }
	 return()
}
#--------------------------------------------#

#--------------------------------------------#
#for (DistMethod in c("euclidean")) {
#    for (ClustMethod in c("hybrid","tree")) {
#    	for (dSValue in 0:4) {
#	    Data        <-  t(Reduced.Data.bin)
#	    DataType    <- "bin_reduced"
#            Prefix      <- sprintf("%s/%s_cutreeDynamic_%s_%s_deepSplit%s",Output.Dir,DataType,DistMethod,ClustMethod,dSValue)
#	    FUNCTION(Data,DistMethod,ClustMethod,dSValue,Prefix)
#
#	    Data        <-  t(Reduced.Data.log2fc)
#	    DataType    <- "log2fc_reduced"
#            Prefix      <- sprintf("%s/%s_cutreeDynamic_%s_%s_deepSplit%s",Output.Dir,DataType,DistMethod,ClustMethod,dSValue)
#	    FUNCTION(Data,DistMethod,ClustMethod,dSValue,Prefix)
#	}
#    }
#}
#--------------------------------------------#

#--------------------------------------------#
for (DistMethod in c("euclidean","canberra")) {
    for (ClustMethod in c("hybrid","tree")) {
    	for (dSValue in 0:4) {
	    Data        <-  Reduced.Data.bin.pca.res$ind$coord
	    DataType    <- "PCA_bin_reduced"
            Prefix      <- sprintf("%s/%s_cutreeDynamic_%s_%s_deepSplit%s",Output.Dir,DataType,DistMethod,ClustMethod,dSValue)
	    FUNCTION(Data,DistMethod,ClustMethod,dSValue,Prefix)

	    Data        <-  Reduced.Data.log2fc.pca.res$ind$coord
	    DataType    <- "PCA_log2fc_reduced"
            Prefix      <- sprintf("%s/%s_cutreeDynamic_%s_%s_deepSplit%s",Output.Dir,DataType,DistMethod,ClustMethod,dSValue)
	    FUNCTION(Data,DistMethod,ClustMethod,dSValue,Prefix)
	}
    }
}
#--------------------------------------------#

quit()
