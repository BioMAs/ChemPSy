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
RData.File     <- commandArgs(TRUE)[1]
Output.Dir     <- commandArgs(TRUE)[2]

#RData.File       <- "processed_data/Rattus_norvegicus/KIDNEY/2_merge_data/data.RData"
#Output.Dir       <- "test"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
load(RData.File)
#--------------------------------------------#

#--------------------------------------------#
# Chemicals altering at least the expression of one gene are selected. Genes altered in at least one chemical are selected.
Reduced.Data.bin    <- Data.bin[,which( apply(abs(Data.bin),2,sum,na.rm=TRUE) >= 10)]
Reduced.Data.bin    <- Reduced.Data.bin[which( apply(abs(Reduced.Data.bin),1,sum,na.rm=TRUE) >= 1),]
Reduced.Data.log2fc <- Data.log2fc[rownames(Reduced.Data.bin),colnames(Reduced.Data.bin)]
Reduced.Data.pBH    <- Data.pBH[rownames(Reduced.Data.bin),colnames(Reduced.Data.bin)]
#--------------------------------------------#

#--------------------------------------------#
write.table(Reduced.Data.log2fc, file=sprintf("%s/data_log2fc.reduced.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
write.table(Reduced.Data.bin   , file=sprintf("%s/data_binary.reduced.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
write.table(Reduced.Data.pBH   , file=sprintf("%s/data_pBH.reduced.txt"   ,Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
save(list=ls()[grep("^Reduced.|^Data.",ls())], file = sprintf("%s/data.reduced.RData",Output.Dir))
#--------------------------------------------#

quit()
