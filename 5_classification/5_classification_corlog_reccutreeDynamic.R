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

#PCA.RData.File   <- "test/PCA.RData"
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
if("dynamicTreeCut" %in% rownames(installed.packages()) == FALSE) {install.packages("dynamicTreeCut", lib=RLIBS)}
if("flashClust" %in% rownames(installed.packages()) == FALSE) {install.packages("flashClust", lib=RLIBS)}
library(dynamicTreeCut)
library(flashClust)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
Rec.CutTree <- function(Data,DistMethod,ClustMethod,dSValue,nMax=50,Classes=c()) {
   print("passage")
   if ( length(Classes) == 0 ) { Classes <- rep(1,nrow(Data)) }
   
   if (DistMethod == "correlation"){
       dist.res   <- cor(t(Data))
       dist.res   <- 1-dist.res
       dist.res   <- as.dist(dist.res)
       hclu.res   <- hclust(dist.res, method='ward')
       
   }
   
   if (DistMethod == "logEuclidean"){
       dist.res   <- dist(Data,method="euclidean")
       dist.res   <-log(dist.res)
       hclu.res   <- hclust(dist.res, method='ward')
   }
   

   cutree.res <- c()
   try(cutree.res <- as.vector(cutreeDynamic(hclu.res,method=ClustMethod,distM=as.matrix(dist.res),deepSplit=dSValue)),silent=TRUE)
   if (is.null(cutree.res) == TRUE) {return(Classes)}
   if (length(unique(cutree.res)) <= 1) {return(Classes)}
   

   ki <- 0
   for (class in unique(cutree.res)) {
       ki <- ki+1
       Indexes <- which(cutree.res == class)
       Classes[Indexes] <- sprintf("%s.%s",Classes[Indexes],ki)
       if (length(Indexes)<=nMax) {next()}
       	Classes[Indexes] <- Rec.CutTree(Data[Indexes,],DistMethod,ClustMethod,dSValue,nMax,Classes[Indexes])
   	   }
   return(Classes)
}
#--------------------------------------------#

#--------------------------------------------#
FUNCTION <- function(Data,DistMethod,ClustMethod,dSValue,nMax,Prefix) {
	 if (ClustMethod == "tree" && dSValue > 1) {return ()}
	 cat(Prefix,"is running...\n")

	 recct.class <- Rec.CutTree(Data,DistMethod,ClustMethod,dSValue,nMax)

  	 K <- length(unique(recct.class))
	 k <- 0	 
	 for (class in unique(recct.class)) {
	     k <- k+1
   	     Individuals <- rownames(Data)[which(recct.class == class)]
 	     write.table(matrix(Individuals,ncol=1),file = sprintf("%s_C%s.%s.txt",Prefix,K,k), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
	 }
	 return()
}
#--------------------------------------------#
ClustMethod <- "hybrid"
dSValue     <- 1
#--------------------------------------------#

#DistMethod  <- "logEuclidean"
#Data        <-  Reduced.Data.bin.pca.res$ind$coord
#DataType    <- "PCA_bin_reduced"
#Prefix      <- sprintf("%s/%s_reccutreeDynamic_%s_%s_deepSplit%s",Output.Dir,DataType,DistMethod,ClustMethod,dSValue)
#FUNCTION(Data,DistMethod,ClustMethod,dSValue,nMax,Prefix)


DistMethod  <- "correlation"
Data        <-  Reduced.Data.bin.pca.res$ind$coord
DataType    <- "PCA_bin_reduced"
Prefix      <- sprintf("%s/%s_reccutreeDynamic_%s_%s_deepSplit%s",Output.Dir,DataType,DistMethod,ClustMethod,dSValue)
FUNCTION(Data,DistMethod,ClustMethod,dSValue,nMax,Prefix)

#--------------------------------------------#

quit()
