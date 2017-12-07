#--------------------------------------------#
# DrugMatrix downloaded at ftp://anonftp.niehs.nih.gov/drugmatrix
# folder organization is similar
#--------------------------------------------#
 
#--------------------------------------------#
RLIBS<-"/home/genouest/irset/archives/softs/R/3.2.3/libs"
setwd("/home/genouest/irset/archives/projects/ChemPSy")
.libPaths(RLIBS)
#--------------------------------------------#
source("scripts/heatmap.R")
#--------------------------------------------#

#--------------------------------------------#
STEP5.Dir                    <- commandArgs(TRUE)[1]
STEP6.Dir                    <- commandArgs(TRUE)[2]
CTD_diseases_data.RData.File <- commandArgs(TRUE)[3]
PCA.RData.File               <- commandArgs(TRUE)[4]
Output.Dir                   <- commandArgs(TRUE)[5]

print(STEP5.Dir)
print(STEP6.Dir)
print(CTD_diseases_data.RData.File)
print(PCA.RData.File)
print(Output.Dir)
#--------------------------------------------#

#--------------------------------------------#
print("Create DIR")
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
print("LOAD")
load(CTD_diseases_data.RData.File)
print("LOAD CTD OK")
load(PCA.RData.File)
print("LOAD PCA OK")
#--------------------------------------------#

#--------------------------------------------#
print("GROUP FILE")
All.Group.Files   <- list.files(path=STEP5.Dir,pattern=".txt")
All.Group.Files   <- sprintf("%s/%s",STEP5.Dir,All.Group.Files[grep("_C[0-9]+\\.[0-9]+\\.txt$",All.Group.Files)])
#--------------------------------------------#

#--------------------------------------------#
All.CTD.IDs       <- c()
for (Group.File in All.Group.Files) {
    Group.IDs     <- as.matrix(read.table(Group.File,header=FALSE,sep="\t",quote=""))[,1]
    Group.CTD.IDs <- unique(gsub("^MESH:","",CTD.Chemicals.Data[intersect(rownames(CTD.Chemicals.Data),Group.IDs),3]))
    All.CTD.IDs   <- c(All.CTD.IDs,Group.CTD.IDs)
}
All.CTD.IDs       <- unique(All.CTD.IDs)
#--------------------------------------------#

#--------------------------------------------#
N      <- length(All.CTD.IDs)
All.ns <- apply(CTD.Chemical2Disease.Data[All.CTD.IDs,],2,sum,na.rm=TRUE)
All.Ns <- rep(N,length(All.ns))
for (Group.File in All.Group.Files) {
    Group.IDs     <- as.matrix(read.table(Group.File,header=FALSE,sep="\t",quote=""))[,1]
    Group.CTD.IDs <- unique(gsub("^MESH:","",CTD.Chemicals.Data[intersect(rownames(CTD.Chemicals.Data),Group.IDs),3]))
    Group.Disease.Data <- CTD.Chemical2Disease.Data[Group.CTD.IDs,]
    R             <- length(Group.CTD.IDs)
    Rs            <- rep(R,length(All.ns))
    Ns            <- All.Ns
    ns            <- All.ns

    rs            <- apply(Group.Disease.Data,2,sum,na.rm=TRUE)
    rRs           <- rs/Rs
    MESH.IDs      <- colnames(CTD.Chemical2Disease.Data)

    indexes1       <- which(ns > 0)
    if (length(indexes1) <= 0) {next()}

    MESH.IDs      <- MESH.IDs[indexes1]	
    rs            <- rs[indexes1]
    Rs            <- Rs[indexes1]
    Ns            <- Ns[indexes1]
    ns            <- ns[indexes1]
    rRs           <- rRs[indexes1]
   
    pvalues       <- phyper(rs-1,ns,Ns-ns,Rs,lower.tail=FALSE)
    pvalues.BH    <- p.adjust(pvalues,method="BH")

    indexes2       <- which(rs >= 3 & rRs >= 0.25)
    if (length(indexes2) <= 0) {next()}
    MESH.IDs       <- MESH.IDs[indexes2]

    Included.CTD.IDs <- c()
    Excluded.CTD.IDs <- c()
    for (MESH.ID in MESH.IDs) {
    	incl.CTD.IDs <- rownames(Group.Disease.Data)[which(Group.Disease.Data[,MESH.ID] >0 )]
	excl.CTD.IDs <- setdiff(Group.CTD.IDs,incl.CTD.IDs)

	Included.CTD.IDs <- c(Included.CTD.IDs ,paste(incl.CTD.IDs,sep="|",collapse="|"))
	Excluded.CTD.IDs <- c(Excluded.CTD.IDs ,paste(excl.CTD.IDs,sep="|",collapse="|"))
    }
    
    Data           <- matrix(nrow=length(indexes2),ncol=10)
    colnames(Data) <- c("MESH","r","R","n","N","r/R","p","pBH","Included.CTD.IDs","Excluded.CTD.IDs")
    Data[,1]       <- MESH.IDs
    Data[,2]       <- rs[indexes2]
    Data[,3]       <- Rs[indexes2]
    Data[,4]       <- ns[indexes2]
    Data[,5]       <- Ns[indexes2]
    Data[,6]       <- rRs[indexes2]
    Data[,7]       <- pvalues[indexes2]
    Data[,8]       <- pvalues.BH[indexes2]
    Data[,9]       <- Included.CTD.IDs 
    Data[,10]      <- Excluded.CTD.IDs
    
    if (length(indexes2) > 1) {
	Data           <- Data[sort(as.double(Data[,7]),decreasing=FALSE,index.return=TRUE)$ix,]
    }
    cat(basename(Group.File),"\n")
    write.table(Data,file = gsub(".txt$",".chem2enr.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)
}
#--------------------------------------------#


#--------------------------------------------#
N.HGIDs <- unique(Gene2HomoloGene.Data[intersect(rownames(Data.bin),rownames(Gene2HomoloGene.Data)),3])
N       <- length(N.HGIDs)
All.ns  <- apply(CTD.HomoloGene2Disease.Data,2,sum,na.rm=TRUE)
All.Ns <- rep(N,length(All.ns))
for (Group.File in All.Group.Files) {
    Group.Gene.File <- sprintf("%s/%s",STEP6.Dir,gsub(".txt$",".GeneIDs.txt",basename(Group.File)))
    if (file.exists(Group.Gene.File) == FALSE || file.info(Group.Gene.File) <= 0) {next()}
    Group.GeneIDs       <- as.character(as.matrix(read.table(Group.Gene.File,header=FALSE,sep="\t",quote=""))[,1])
    Group.HGIDs         <- as.character(unique(Gene2HomoloGene.Data[intersect(Group.GeneIDs,rownames(Gene2HomoloGene.Data)),3]))
    Group.Disease.Data  <- CTD.HomoloGene2Disease.Data[intersect(Group.HGIDs,rownames(CTD.HomoloGene2Disease.Data)),]
    R             <- length(Group.HGIDs)
    Rs            <- rep(R,length(All.ns))
    Ns            <- All.Ns
    ns            <- All.ns

    rs            <- apply(Group.Disease.Data,2,sum,na.rm=TRUE)
    rRs           <- rs/Rs
    MESH.IDs      <- colnames(CTD.HomoloGene2Disease.Data)

    indexes1       <- which(ns > 0)
    if (length(indexes1) <= 0) {next()}

    MESH.IDs      <- MESH.IDs[indexes1]	
    rs            <- rs[indexes1]
    Rs            <- Rs[indexes1]
    Ns            <- Ns[indexes1]
    ns            <- ns[indexes1]
    rRs           <- rRs[indexes1]
   
    pvalues       <- phyper(rs-1,ns,Ns-ns,Rs,lower.tail=FALSE)
    pvalues.BH    <- p.adjust(pvalues,method="BH")

    indexes2       <- which(rs >= 3 & rRs >= 0.25)
    if (length(indexes2) <= 0) {next()}
    MESH.IDs       <- MESH.IDs[indexes2]

    Included.Gene.IDs <- c()
    Excluded.Gene.IDs <- c()
    for (MESH.ID in MESH.IDs) {
    	incl.Gene.IDs <- rownames(Group.Disease.Data)[which(Group.Disease.Data[,MESH.ID] >0 )]
	excl.Gene.IDs <- setdiff(Group.HGIDs,incl.Gene.IDs)

	Included.Gene.IDs <- c(Included.Gene.IDs ,paste(incl.Gene.IDs,sep="|",collapse="|"))
	Excluded.Gene.IDs <- c(Excluded.Gene.IDs ,paste(excl.Gene.IDs,sep="|",collapse="|"))
    }

    Data           <- matrix(nrow=length(indexes2),ncol=10)
    colnames(Data) <- c("MESH","r","R","n","N","r/R","p","pBH","Included.Gene.IDs","Excluded.Gene.IDs")
    Data[,1]       <- MESH.IDs
    Data[,2]       <- rs[indexes2]
    Data[,3]       <- Rs[indexes2]
    Data[,4]       <- ns[indexes2]
    Data[,5]       <- Ns[indexes2]
    Data[,6]       <- rRs[indexes2]
    Data[,7]       <- pvalues[indexes2]
    Data[,8]       <- pvalues.BH[indexes2]
    Data[,9]       <- Included.Gene.IDs 
    Data[,10]      <- Excluded.Gene.IDs
    
    if (length(indexes2) > 1) {
	Data           <- Data[sort(as.double(Data[,7]),decreasing=FALSE,index.return=TRUE)$ix,]
    }

    cat(basename(Group.File),"\n")
    write.table(Data,file = gsub(".txt$",".gene2enr.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = NA , row.names = TRUE , quote=FALSE)
}
#--------------------------------------------#

quit()





#--------------------------------------------#
N.CTD.IDs       <- c()
for (Group.File in All.Group.Files) {
    Group.IDs     <- as.matrix(read.table(Group.File,header=FALSE,sep="\t",quote=""))[,1]
    Group.CTD.IDs <- unique(gsub("^MESH:","",CTD.Chemicals.Data[intersect(rownames(CTD.Chemicals.Data),Group.IDs),3]))
    N.CTD.IDs     <- c(N.CTD.IDs,Group.CTD.IDs)
}
N.CTD.IDs         <- unique(N.CTD.IDs)
N.CTD <- length(N.CTD.IDs)
#--------------------------------------------#

#--------------------------------------------#
N.HG.IDs <- unique(Gene2HomoloGene.Data[intersect(rownames(Data.bin),rownames(Gene2HomoloGene.Data)),3])
N.HG    <- length(N.HG.IDs)
#--------------------------------------------#
for (Group.File in All.Group.Files) {
    if (file.exists(Group.File) == FALSE || file.info(Group.File) <= 0) {next()}
    Group.Gene.File        <- gsub(".txt$",".GeneIDs.txt",Group.File)
    if (file.exists(Group.Gene.File) == FALSE || file.info(Group.Gene.File) <= 0) {next()}
    cat(Group.File,"\n")

    #--------------------------------------------#
    Group.Chemical.IDs     <- as.matrix(read.table(Group.File,header=FALSE,sep="\t",quote=""))[,1]
    Group.CTD.IDs          <- unique(gsub("^MESH:","",CTD.Chemicals.Data[intersect(rownames(CTD.Chemicals.Data),Group.IDs),3]))
    Group.Disease.CTD.Data <- CTD.Chemical2Disease.Data[Group.CTD.IDs,]
    R.CTD                  <- length(Group.CTD.IDs)

    rs.CTD                 <- apply(Group.Disease.CTD.Data,2,sum,na.rm=TRUE)
    Rs.CTD                 <- rep(R.CTD,length(rs.CTD))
    rRs.CTD                <- rs.CTD/Rs.CTD

    indexes.CTD            <- which(rs.CTD >= 3 & rRs.CTD >= 0.25)
    #if (length(indexes.CTD) <= 0) {next()}
    rs.CTD                 <- rs.CTD[indexes.CTD]
    if (length(indexes.CTD) > 1) {
        ns.CTD             <- apply(CTD.Chemical2Disease.Data[N.CTD.IDs,indexes.CTD],2,sum,na.rm=TRUE)
    } else {
        ns.CTD             <- sum(CTD.Chemical2Disease.Data[N.CTD.IDs,indexes.CTD])
    }
    Rs.CTD                 <- Rs.CTD[indexes.CTD]
    Ns.CTD                 <- rep(N.CTD,length(rs.CTD))
    rRs.CTD                <- rRs.CTD[indexes.CTD]

    pvalues.CTD            <- phyper(rs.CTD-1,ns.CTD,Ns.CTD-ns.CTD,Rs.CTD,lower.tail=FALSE)
    pvalues.BH.CTD         <- p.adjust(pvalues.CTD,method="BH")

    Selected.CTD           <- which(pvalues.BH.CTD<=0.05)
    cat(basename(Group.File),"-> CTD -> ",colnames(CTD.Chemical2Disease.Data)[indexes.CTD[Selected.CTD]],"\n")
    #--------------------------------------------#

    #--------------------------------------------#
    Group.Gene.IDs         <- as.character(as.matrix(read.table(Group.Gene.File,header=FALSE,sep="\t",quote=""))[,1])
    Group.HG.IDs           <- as.character(unique(Gene2HomoloGene.Data[intersect(Group.Gene.IDs,rownames(Gene2HomoloGene.Data)),3]))
    Group.Disease.HG.Data  <- CTD.HomoloGene2Disease.Data[intersect(Group.HG.IDs,rownames(CTD.HomoloGene2Disease.Data)),]
    R.HG                   <- length(Group.HG.IDs)

    rs.HG                  <- apply(Group.Disease.HG.Data,2,sum,na.rm=TRUE)
    indexes.HG             <- which(rs.HG >= 3)
    #if (length(indexes.HG) <= 0) {next()}
    rs.HG                  <- rs.HG   [indexes.HG]
    if (length(indexes.HG) > 1) {
        ns.HG              <- apply(CTD.HomoloGene2Disease.Data[,indexes.HG],2,sum,na.rm=TRUE)
    } else {
        ns.HG              <- sum(CTD.HomoloGene2Disease.Data[,indexes.HG])
    }
    Rs.HG                  <- rep(R.HG,length(rs.HG))
    Ns.HG                  <- rep(N.HG,length(rs.HG))

    pvalues.HG             <- phyper(rs.HG-1,ns.HG,Ns.HG-ns.HG,Rs.HG,lower.tail=FALSE)
    pvalues.BH.HG          <- p.adjust(pvalues.HG,method="BH")

    Selected.HG            <- which(pvalues.BH.HG<=0.05)
    cat(basename(Group.File),"-> HG -> ",colnames(CTD.HomoloGene2Disease.Data)[indexes.HG[Selected.HG]],"\n")
    #--------------------------------------------#

    #--------------------------------------------#
    Selected <- intersect(Selected.HG,Selected.CTD)
    if (length(Selected) > 0) {next()}
    cat(Group.File,"->",Selected,"\n")
    #--------------------------------------------#

}

