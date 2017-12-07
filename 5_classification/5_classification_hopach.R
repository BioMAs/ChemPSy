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

#PCA.RData.File <- "test/PCA.reduced.RData"
#Output.Dir     <- "test"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
if("hopach" %in% rownames(installed.packages()) == FALSE) {install.packages("hopach", lib=RLIBS)}
if("flashClust" %in% rownames(installed.packages()) == FALSE) {install.packages("flashClust", lib=RLIBS)}
library(hopach)
library(flashClust)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
FUNCTION <- function(Data,DistMethod,ClustMethod,Prefix) {
	 cat(Prefix,"is running...\n")
	 dist.res <- dist(Data,method=DistMethod)
	 
   	 hopach.res <- hopach(Data,dmat=as.matrix(dist.res),clusters=ClustMethod)
  	 K <- length(unique(hopach.res$clustering$labels))
	 k <- 0	 
	 for (ki in unique(hopach.res$clustering$labels)) {
	     k <- k+1
   	     Individuals <- rownames(Data)[which(hopach.res$clustering$labels == ki)]
	     write.table(matrix(Individuals,ncol=1),file = sprintf("%s_C%s.%s.txt",Prefix,K,k), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
	 }
	 return()
}
#--------------------------------------------#

#--------------------------------------------#
for (DistMethod in c("euclidean")) {
    for (ClustMethod in c("best","greedy")) {
        Data        <-  t(Reduced.Data.bin)
	DataType    <- "bin_reduced"
        Prefix      <- sprintf("%s/%s_hopach_%s_%s",Output.Dir,DataType,DistMethod,ClustMethod)
	try(FUNCTION(Data,DistMethod,ClustMethod,Prefix),silent=TRUE)

	Data        <-  t(Reduced.Data.log2fc)
	DataType    <- "log2fc_reduced"
        Prefix      <- sprintf("%s/%s_hopach_%s_%s",Output.Dir,DataType,DistMethod,ClustMethod)
	try(FUNCTION(Data,DistMethod,ClustMethod,Prefix),silent=TRUE)

        Data        <-  Reduced.Data.bin.pca.res$ind$coord
	DataType    <- "PCA_bin_reduced"
        Prefix      <- sprintf("%s/%s_hopach_%s_%s",Output.Dir,DataType,DistMethod,ClustMethod)
	try(FUNCTION(Data,DistMethod,ClustMethod,Prefix),silent=TRUE)

	Data        <-  Reduced.Data.log2fc.pca.res$ind$coord
	DataType    <- "PCA_log2fc_reduced"
        Prefix      <- sprintf("%s/%s_hopach_%s_%s",Output.Dir,DataType,DistMethod,ClustMethod)
	try(FUNCTION(Data,DistMethod,ClustMethod,Prefix),silent=TRUE)
    }
}
#--------------------------------------------#

quit()
