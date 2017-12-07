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

#PCA.RData.File   <- "processed_data/Rattus_norvegicus/KIDNEY/4_pca_analysis/PCA.reduced.RData"
#Output.Dir       <- "test"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
if("dbscan" %in% rownames(installed.packages()) == FALSE) {install.packages("dbscan", lib=RLIBS)}
library(dbscan)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
FUNCTION <- function(Data,Prefix) {
	 cat(Prefix,"is running...\n")

	 # Default parameters dbscan(x, eps, minPts = 5, weights = NULL,borderPoints = TRUE, search = "kdtree", bucketSize = 10,splitRule = "suggest", approx = 0, ...)
	 dbscan.res <- dbscan(Data,eps=.4)

  	 K <- length(unique(dbscan.res$cluster))
	 k <- 0	 
	 for (ki in unique(dbscan.res$cluster)) {
	     k <- k+1
   	     Individuals <- rownames(Data)[which(dbscan.res$cluster == ki)]
	     write.table(matrix(Individuals,ncol=1),file = sprintf("%s_C%s.%s.txt",Prefix,K,k), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
	 }
	 return()
}
#--------------------------------------------#

#--------------------------------------------#
Data        <-  t(Reduced.Data.bin)
DataType    <- "bin_reduced"
Prefix      <- sprintf("%s/%s_dbscan",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  t(Reduced.Data.log2fc)
DataType    <- "log2fc_reduced"
Prefix      <- sprintf("%s/%s_dbscan",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  Reduced.Data.bin.pca.res$ind$coord
DataType    <- "PCA_bin_reduced"
Prefix      <- sprintf("%s/%s_dbscan",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  Reduced.Data.log2fc.pca.res$ind$coord
DataType    <- "PCA_log2fc_reduced"
Prefix      <- sprintf("%s/%s_dbscan",Output.Dir,DataType)
FUNCTION(Data,Prefix)
#--------------------------------------------#

quit()
