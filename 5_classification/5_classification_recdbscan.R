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
if("dbscan" %in% rownames(installed.packages()) == FALSE) {install.packages("dbscan", lib=RLIBS)}
library(dbscan)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
Rec.dbscan <- function(Data,nMax=50,Classes=c()) {
   if ( length(Classes) == 0 ) { Classes <- rep(1,nrow(Data)) }

   dbscan.res <- c()
   try(dbscan.res <- dbscan(Data,eps=.4),silent=TRUE)
   if (is.null(dbscan.res) == TRUE) {return(Classes)}
   if (length(unique(dbscan.res$cluster)) <= 1) {return(Classes)}

   ki <- 0
   for (class in unique(dbscan.res$cluster)) {
       ki <- ki+1
       Indexes <- which(dbscan.res$cluster == class)
       Classes[Indexes] <- sprintf("%s.%s",Classes[Indexes],ki)

        if (length(Indexes)<=nMax) {next()}
       Classes[Indexes] <- Rec.dbscan(Data[Indexes,],nMax,Classes[Indexes])
   }
   return(Classes)
}
#--------------------------------------------#

#--------------------------------------------#
FUNCTION <- function(Data,nMax,Prefix) {
	 cat(Prefix,"is running...\n")

	 recdbscan.class <- Rec.dbscan(Data,nMax)

  	 K <- length(unique(recdbscan.class))
	 k <- 0	 
	 for (class in unique(recdbscan.class)) {
	     k <- k+1
   	     Individuals <- rownames(Data)[which(recdbscan.class == class)]
	     write.table(matrix(Individuals,ncol=1),file = sprintf("%s_C%s.%s.txt",Prefix,K,k), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
	 }
	 return()
}
#--------------------------------------------#


#--------------------------------------------#
Data        <-  t(Reduced.Data.bin)
DataType    <- "bin_reduced"
Prefix      <- sprintf("%s/%s_recdbscan",Output.Dir,DataType)
FUNCTION(Data,nMax,Prefix)

Data        <-  t(Reduced.Data.log2fc)
DataType    <- "log2fc_reduced"
Prefix      <- sprintf("%s/%s_recdbscan",Output.Dir,DataType)
FUNCTION(Data,nMax,Prefix)

Data        <-  Reduced.Data.bin.pca.res$ind$coord
DataType    <- "PCA_bin_reduced"
Prefix      <- sprintf("%s/%s_recdbscan",Output.Dir,DataType)
FUNCTION(Data,nMax,Prefix)

Data        <-  Reduced.Data.log2fc.pca.res$ind$coord
DataType    <- "PCA_log2fc_reduced"
Prefix      <- sprintf("%s/%s_recdbscan",Output.Dir,DataType)
FUNCTION(Data,nMax,Prefix)
#--------------------------------------------#

quit()
