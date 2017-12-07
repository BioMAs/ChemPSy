#--------------------------------------------#
# DrugMatrix downloaded at ftp://anonftp.niehs.nih.gov/drugmatrix
# folder organization is similar
#--------------------------------------------#

#--------------------------------------------#
RLIBS<-"/groups/irset/archives/softs/R/3.2.3/libs"
setwd("/groups/irset/archives/projects/ChemPSy")
.libPaths(RLIBS)
#--------------------------------------------#

#--------------------------------------------#
FiltrationFiles.File <- commandArgs(TRUE)[1]
Output.Dir       <- commandArgs(TRUE)[2]

#FiltrationFiles.File <- "/home/genouest/irset/archives/projects/ChemPSy/processed_data/Rattus_norvegicus/HEART/HEART.txt"
#Output.Dir           <- "/home/genouest/irset/archives/projects/ChemPSy/test"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
Filtration.Files <- as.vector(read.table(FiltrationFiles.File)[,1])
Filtration.Names <- basename(dirname(Filtration.Files))
#--------------------------------------------#

#--------------------------------------------#
PS.IDs                <- rownames(read.table(Filtration.Files[1],header=TRUE,row.names=1,sep='\t'))

Data.log2fc           <- matrix(0,ncol=length(Filtration.Names),nrow=length(PS.IDs))
rownames(Data.log2fc) <- PS.IDs
colnames(Data.log2fc) <- Filtration.Names

Data.bin              <- matrix(0,ncol=length(Filtration.Names),nrow=length(PS.IDs))
rownames(Data.bin)    <- PS.IDs
colnames(Data.bin)    <- Filtration.Names

Data.pBH              <- matrix(NA,ncol=length(Filtration.Names),nrow=length(PS.IDs))
rownames(Data.pBH)    <- PS.IDs
colnames(Data.pBH)    <- Filtration.Names

for (i in 1:length(Filtration.Files)) {
	Filtration.File <- Filtration.Files[i]
	Filtration.Name <- Filtration.Names[i]

	Filtration.Data               <- read.table(Filtration.File,header=TRUE,row.names=1,sep='\t')
	Filtration.IDs                <- intersect(PS.IDs,rownames(Filtration.Data)[which(Filtration.Data[,1] == 1 & abs(Filtration.Data[,3]) == 1 & abs(Filtration.Data[,4]) == 1)])

	cat(i,"->",Filtration.Name,"->",length(Filtration.IDs),"\n")

	Data.bin[Filtration.IDs,i]    <- Filtration.Data[Filtration.IDs,4]  
	Data.log2fc[Filtration.IDs,i] <- Filtration.Data[Filtration.IDs,5] 


	Filtration.IDs                <- intersect(PS.IDs,rownames(Filtration.Data)[which(is.na(Filtration.Data[,7]) == FALSE)])
	Data.pBH[Filtration.IDs,i]    <- Filtration.Data[Filtration.IDs,7]  
}

write.table(Data.log2fc, file=sprintf("%s/data_log2fc.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
write.table(Data.bin   , file=sprintf("%s/data_binary.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
write.table(Data.pBH   , file=sprintf("%s/data_pBH.txt"   ,Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
save(Data.log2fc,Data.bin,Data.pBH, file = sprintf("%s/data.RData",Output.Dir))
#--------------------------------------------#

quit()
