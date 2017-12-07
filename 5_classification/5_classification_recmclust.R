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

#PCA.RData.File   <- "processed_data/Rattus_norvegicus/KIDNEY/4_pca_analysis/PCA.RData"
#Output.Dir       <- "test"
#--------------------------------------------#

#--------------------------------------------#
nMax <- 50
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
if("mclust" %in% rownames(installed.packages()) == FALSE) {install.packages("mclust", lib=RLIBS)}
library(mclust)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
Rec.Mclust <- function(Data,nMax=50,Classes=c()) {
   if ( length(Classes) == 0 ) { Classes <- rep(1,nrow(Data)) }

   mclust.res <- c()
   try(mclust.res <- Mclust(Data),silent=TRUE)
   if (is.null(mclust.res) == TRUE) {return(Classes)}
   if (length(unique(mclust.res$classification)) <= 1) {return(Classes)}

   ki <- 0
   for (class in unique(mclust.res$classification)) {
       ki <- ki+1
       Indexes <- which(mclust.res$classification == class)
       Classes[Indexes] <- sprintf("%s.%s",Classes[Indexes],ki)

        if (length(Indexes)<=nMax) {next()}
       Classes[Indexes] <- Rec.Mclust(Data[Indexes,],nMax,Classes[Indexes])
   }
   return(Classes)
}
#--------------------------------------------#

#--------------------------------------------#
FUNCTION <- function(Data,nMax,Prefix) {
	 cat(Prefix,"is running...\n")

	 recmclust.class <- Rec.Mclust(Data,nMax)

  	 K <- length(unique(recmclust.class))
	 k <- 0	 
	 for (class in unique(recmclust.class)) {
	     k <- k+1
   	     Individuals <- rownames(Data)[which(recmclust.class == class)]
	     write.table(matrix(Individuals,ncol=1),file = sprintf("%s_C%s.%s.txt",Prefix,K,k), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
	 }
	 return()
}
#--------------------------------------------#


#--------------------------------------------#
Data        <-  t(Reduced.Data.bin)
DataType    <- "bin_reduced"
Prefix      <- sprintf("%s/%s_recmclust",Output.Dir,DataType)
FUNCTION(Data,nMax,Prefix)

Data        <-  t(Reduced.Data.log2fc)
DataType    <- "log2fc_reduced"
Prefix      <- sprintf("%s/%s_recmclust",Output.Dir,DataType)
FUNCTION(Data,nMax,Prefix)

Data        <-  Reduced.Data.bin.pca.res$ind$coord
DataType    <- "PCA_bin_reduced"
Prefix      <- sprintf("%s/%s_recmclust",Output.Dir,DataType)
FUNCTION(Data,nMax,Prefix)

Data        <-  Reduced.Data.log2fc.pca.res$ind$coord
DataType    <- "PCA_log2fc_reduced"
Prefix      <- sprintf("%s/%s_recmclust",Output.Dir,DataType)
FUNCTION(Data,nMax,Prefix)
#--------------------------------------------#

quit()
