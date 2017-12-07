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
Input.Dir                    <- commandArgs(TRUE)[1]
PCA.RData.File               <- commandArgs(TRUE)[2]
Output.Dir                   <- commandArgs(TRUE)[3]

if (0==1) {
  Input.Dir                     <- "test_all_recmclust"
  PCA.RData.File                <- "processed_data/Rattus_norvegicus/ALL/3_reduce_data/4_pca_analysis/PCA.RData"
  Output.Dir                    <- "test"
}

rR.cutoff  <- 0.8
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
All.Group.Files   <- list.files(path=Input.Dir,pattern=".txt")
All.Group.Files   <- sprintf("%s/%s",Input.Dir,All.Group.Files[grep("_C[0-9]+\\.[0-9]+\\.txt$",All.Group.Files)])
#--------------------------------------------#

#--------------------------------------------#
All.Chemical.IDs       <- c()
for (Group.File in All.Group.Files) {
    Group.Chemical.IDs <- as.matrix(read.table(Group.File,header=FALSE,sep="\t",quote=""))[,1]
    All.Chemical.IDs   <- c(All.Chemical.IDs,Group.Chemical.IDs)
}
All.Chemical.IDs       <- unique(All.Chemical.IDs)
#--------------------------------------------#

#--------------------------------------------#
for (Group.File in All.Group.Files) {
    Group.Chemical.IDs <- unique(as.matrix(read.table(Group.File,header=FALSE,sep="\t",quote=""))[,1])
    if (length(Group.Chemical.IDs) < 3) {next()}

    R <- length(Group.Chemical.IDs)

    Sum <- apply(Reduced.Data.bin[,Group.Chemical.IDs],1,sum)

    Up.GeneIDs  <- c()
    Up.Indexes  <- which(Sum > 0)
    if (length(Up.Indexes)>0) {
       if (length(Up.Indexes)>1) {
       	  Up.rs       <- apply(abs(Reduced.Data.bin[Up.Indexes,Group.Chemical.IDs]),1,sum)
       } else {
       	  Up.rs       <- sum(abs(Reduced.Data.bin[Up.Indexes,Group.Chemical.IDs]))
       }
       Up.Rs       <- rep(R,length(Up.rs))
       Up.rRs      <- Up.rs/Up.Rs
 
       Up.Selected <- which(Up.rRs >= rR.cutoff)
       Up.GeneIDs  <- rownames(Reduced.Data.bin)[Up.Indexes[Up.Selected]]
    }

    Down.GeneIDs  <- c()
    Down.Indexes  <- which(Sum < 0)
    if (length(Down.Indexes)>0) {
       if (length(Down.Indexes)>1) {
          Down.rs       <- apply(abs(Reduced.Data.bin[Down.Indexes,Group.Chemical.IDs]),1,sum)
       } else {
          Down.rs       <- sum(abs(Reduced.Data.bin[Down.Indexes,Group.Chemical.IDs]))
       }
       Down.Rs       <- rep(R,length(Down.rs))
       Down.rRs      <- Down.rs/Down.Rs
 
       Down.Selected <- which(Down.rRs >= rR.cutoff)
       Down.GeneIDs  <- rownames(Reduced.Data.bin)[Down.Indexes[Down.Selected]]
    }
    Selected.GeneIDs <- unique(c(Up.GeneIDs,Down.GeneIDs))
    if (length(Selected.GeneIDs) > 0) {
       cat(basename(Group.File),"-> all ->",length(Selected.GeneIDs),"\n")
       write.table(matrix(Selected.GeneIDs,ncol=1),file = gsub(".txt$",".GeneIDs.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
    }
}

quit()



#--------------------------------------------#
rR.cutoff  <- 0.5
pBH.cutoff <- 0.05
N <- length(All.Chemical.IDs)
for (Group.File in All.Group.Files) {
    Group.Chemical.IDs <- unique(as.matrix(read.table(Group.File,header=FALSE,sep="\t",quote=""))[,1])
    if (length(Group.Chemical.IDs) < 3) {next()}

    R <- length(Group.Chemical.IDs)

    Sum <- apply(Reduced.Data.bin[,Group.Chemical.IDs],1,sum)

    Up.GeneIDs  <- c()
    Up.Indexes  <- which(Sum > 0)
    if (length(Up.Indexes)>0) {
       if (length(Up.Indexes)>1) {
       	  Up.rs       <- apply(abs(Reduced.Data.bin[Up.Indexes,Group.Chemical.IDs]),1,sum)
       	  Up.ns       <- apply(abs(Reduced.Data.bin[Up.Indexes,]),1,sum)
       } else {
       	  Up.rs       <- sum(abs(Reduced.Data.bin[Up.Indexes,Group.Chemical.IDs]))
       	  Up.ns       <- sum(abs(Reduced.Data.bin[Up.Indexes,]))
       }
       Up.Rs       <- rep(R,length(Up.rs))
       Up.rRs      <- Up.rs/Up.Rs
       Up.Ns       <- rep(N,length(Up.rs))
       Up.pvalues  <- p.adjust(phyper(Up.rs-1,Up.ns,Up.Ns-Up.ns,Up.Rs,lower.tail=FALSE),method="BH")

       Up.Selected <- which(Up.rRs >= rR.cutoff)
       Up.GeneIDs  <- rownames(Reduced.Data.bin)[Up.Indexes[Up.Selected]]

       SpecUp.Selected <- which(Up.rRs >= rR.cutoff & Up.pvalues <= pBH.cutoff)
       SpecUp.GeneIDs  <- rownames(Reduced.Data.bin)[Up.Indexes[SpecUp.Selected]]
    }

    Down.GeneIDs  <- c()
    Down.Indexes  <- which(Sum < 0)
    if (length(Down.Indexes)>0) {
       if (length(Down.Indexes)>1) {
          Down.rs       <- apply(abs(Reduced.Data.bin[Down.Indexes,Group.Chemical.IDs]),1,sum)
          Down.ns       <- apply(abs(Reduced.Data.bin[Down.Indexes,]),1,sum)
       } else {
          Down.rs       <- sum(abs(Reduced.Data.bin[Down.Indexes,Group.Chemical.IDs]))
          Down.ns       <- sum(abs(Reduced.Data.bin[Down.Indexes,]))
       }
       Down.Rs       <- rep(R,length(Down.rs))
       Down.rRs      <- Down.rs/Down.Rs
       Down.Ns       <- rep(N,length(Down.rs))
       Down.pvalues  <- p.adjust(phyper(Down.rs-1,Down.ns,Down.Ns-Down.ns,Down.Rs,lower.tail=FALSE),method="BH")

       Down.Selected <- which(Down.rRs >= rR.cutoff)
       Down.GeneIDs  <- rownames(Reduced.Data.bin)[Down.Indexes[Down.Selected]]

       SpecDown.Selected <- which(Down.rRs >= rR.cutoff & Down.pvalues <= pBH.cutoff)
       SpecDown.GeneIDs  <- rownames(Reduced.Data.bin)[Down.Indexes[SpecDown.Selected]]
    }
    Selected.GeneIDs <- unique(c(Up.GeneIDs,Down.GeneIDs))
    if (length(Selected.GeneIDs) > 0) {
       cat(basename(Group.File),"-> all ->",length(Selected.GeneIDs),"\n")
       write.table(matrix(Selected.GeneIDs,ncol=1),file = gsub(".txt$",".GeneIDs.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
    }

    Spec.Selected.GeneIDs <- unique(c(SpecUp.GeneIDs,SpecDown.GeneIDs))
    if (length(Spec.Selected.GeneIDs) > 0) {
       cat(basename(Group.File),"-> specific ->",length(Spec.Selected.GeneIDs),"\n")
       write.table(matrix(Spec.Selected.GeneIDs,ncol=1),file = gsub(".txt$",".Spec.GeneIDs.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
    }

    next()
    if (length(Up.GeneIDs) > 0) {
      write.table(matrix(Up.GeneIDs,ncol=1),file = gsub(".txt$",".Up.GeneIDs.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
    }
    if (length(Down.GeneIDs) > 0) {
      write.table(matrix(Down.GeneIDs,ncol=1),file = gsub(".txt$",".Down.GeneIDs.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
    }
    if (length(SpecUp.GeneIDs) > 0) {
      write.table(matrix(SpecUp.GeneIDs,ncol=1),file = gsub(".txt$",".SpecUp.GeneIDs.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
    }
    if (length(SpecDown.GeneIDs) > 0) {
      write.table(matrix(SpecDown.GeneIDs,ncol=1),file = gsub(".txt$",".SpecDown.GeneIDs.txt",sprintf("%s/%s",Output.Dir,basename(Group.File))), sep = '\t' , col.names = FALSE , row.names = FALSE , quote=FALSE)
    }
}
#--------------------------------------------#


quit()


