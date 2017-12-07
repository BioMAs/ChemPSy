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
if("hopach" %in% rownames(installed.packages()) == FALSE) {install.packages("hopach", lib=RLIBS)}
if("flashClust" %in% rownames(installed.packages()) == FALSE) {install.packages("flashClust", lib=RLIBS)}
library(hopach)
library(flashClust)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
Rec.hopach <- function(Data,DistMethod="euclidean",ClustMethod="greedy",nMax=50,Classes=c()) {
   if ( length(Classes) == 0 ) { Classes <- rep(1,nrow(Data)) }
   dist.res <- dist(Data,method=DistMethod)
   hopach.res <- c()
   try(hopach.res <- hopach(Data,dmat=as.matrix(dist.res),clusters=ClustMethod),silent=TRUE)
   if (is.null(hopach.res) == TRUE) {return(Classes)}
   if (length(unique(hopach.res$clustering$labels)) <= 1) {return(Classes)}

   ki <- 0
   for (class in unique(hopach.res$clustering$labels)) {
       ki <- ki+1
       Indexes <- which(hopach.res$clustering$labels == class)
       Classes[Indexes] <- sprintf("%s.%s",Classes[Indexes],ki)
       if (length(Indexes)<=nMax) {next()}
       Classes[Indexes] <- Rec.hopach(Data[Indexes,],DistMethod,ClustMethod,nMax,Classes[Indexes])
   }
   return(Classes)
}
#--------------------------------------------#

#--------------------------------------------#
FUNCTION <- function(Data,DistMethod="euclidean",ClustMethod="best",nMax=50,Prefix) {
	 cat(Prefix,"is running...\n")

	 rechopach.class <- Rec.hopach(Data,DistMethod,ClustMethod,nMax)

  	 K <- length(unique(rechopach.class))
	 k <- 0	 
	 for (class in unique(rechopach.class)) {
	     k <- k+1
   	     Individuals <- rownames(Data)[which(rechopach.class == class)]
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
        Prefix      <- sprintf("%s/%s_rechopach_%s_%s",Output.Dir,DataType,DistMethod,ClustMethod)
	FUNCTION(Data,DistMethod,ClustMethod,nMax,Prefix)

	Data        <-  t(Reduced.Data.log2fc)
	DataType    <- "log2fc_reduced"
        Prefix      <- sprintf("%s/%s_rechopach_%s_%s",Output.Dir,DataType,DistMethod,ClustMethod)
	FUNCTION(Data,DistMethod,ClustMethod,nMax,Prefix)

        Data        <-  Reduced.Data.bin.pca.res$ind$coord
	DataType    <- "PCA_bin_reduced"
        Prefix      <- sprintf("%s/%s_rechopach_%s_%s",Output.Dir,DataType,DistMethod,ClustMethod)
	FUNCTION(Data,DistMethod,ClustMethod,nMax,Prefix)

	Data        <-  Reduced.Data.log2fc.pca.res$ind$coord
	DataType    <- "PCA_log2fc_reduced"
        Prefix      <- sprintf("%s/%s_rechopach_%s_%s",Output.Dir,DataType,DistMethod,ClustMethod)
	FUNCTION(Data,DistMethod,ClustMethod,nMax,Prefix)
    }
}
#--------------------------------------------#

quit()
