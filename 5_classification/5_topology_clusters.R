#--------------------------------------------#
# DrugMatrix downloaded at ftp://anonftp.niehs.nih.gov/drugmatrix
# folder organization is similar
#--------------------------------------------#

#--------------------------------------------#
setwd("/home/genouest/irset/archives/projects/ChemPSy")
.libPaths("scripts/RLIBRARIES")
source("scripts/heatmap.R")
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
if("mclust" %in% rownames(installed.packages()) == FALSE) {install.packages("mclust", lib="scripts/RLIBRARIES")}
if("flashClust" %in% rownames(installed.packages()) == FALSE) {install.packages("flashClust", lib="scripts/RLIBRARIES")}
library(mclust)
library(flashClust)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
FUNCTION <- function(Data,Prefix) {
	 cat(Prefix,"is running...\n")

	 mclust.res <- Mclust(Data)

  	 K <- length(unique(mclust.res$classification))
	 k <- 0	 
	 for (ki in unique(mclust.res$classification)) {
	     k <- k+1
   	     Individuals <- rownames(Data)[which(mclust.res$classification == ki)]
	     write.table(matrix(Individuals,ncol=1),file = sprintf("%s_C%s.%s.txt",Prefix,K,k), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
	 }
	 return()
}
#--------------------------------------------#

#--------------------------------------------#
Data        <-  t(Reduced.Data.bin)
DataType    <- "bin_reduced"
Prefix      <- sprintf("%s/%s_mclust",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  t(Reduced.Data.log2fc)
DataType    <- "log2fc_reduced"
Prefix      <- sprintf("%s/%s_mclust",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  Reduced.Data.bin.pca.res$ind$coord
DataType    <- "PCA_bin_reduced"
Prefix      <- sprintf("%s/%s_mclust",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  Reduced.Data.log2fc.pca.res$ind$coord
DataType    <- "PCA_log2fc_reduced"
Prefix      <- sprintf("%s/%s_mclust",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  t(Data.bin)
DataType    <- "bin"
Prefix      <- sprintf("%s/%s_mclust",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  t(Data.log2fc)
DataType    <- "log2fc"
Prefix      <- sprintf("%s/%s_mclust",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  Data.bin.pca.res$ind$coord
DataType    <- "PCA_bin"
Prefix      <- sprintf("%s/%s_mclust",Output.Dir,DataType)
FUNCTION(Data,Prefix)

Data        <-  Data.log2fc.pca.res$ind$coord
DataType    <- "PCA_log2fc"
Prefix      <- sprintf("%s/%s_mclust",Output.Dir,DataType)
FUNCTION(Data,Prefix)
#--------------------------------------------#

quit()
