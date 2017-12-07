RLIBS<-"/home/genouest/irset/archives/softs/R/3.1.0/libs"
.libPaths(RLIBS)


source("http://bioconductor.org/biocLite.R")
if("scales" %in% rownames(installed.packages()) == FALSE) {install.packages("scales", lib=RLIBS)}
if("flashClust" %in% rownames(installed.packages()) == FALSE) {install.packages("flashClust", lib=RLIBS)}
library(scales)
library(flashClust)

#--------------------------------------------#
HeatMap <- function(Data,Group.Files,OUT.File) { 
  #--------------------------------------------#
  Q.Val  <- quantile(Data,p=c(0,0.05,0.95,1),na.rm=TRUE)
  q5.Val <- Q.Val[2]
  q95.Val <- Q.Val[3]
  min.Val <- floor(Q.Val[1])
	max.Val <- ceiling(Q.Val[4])
  q.Val <- max(abs(c(q5.Val,q95.Val)))
  m.Val <- max(abs(c(min.Val,max.Val)))
    
 	Scale <- c(-1*m.Val,extended_breaks(n=10)(floor(-1*q.Val):ceiling(q.Val)),m.Val)
	n <- (length(Scale)-1)/2
	nLeft  <- n
	nRight <- n
	Colors <- c(rgb((0:(nLeft-1))/nLeft,(0:(nLeft-1))/nLeft,(nLeft-1)/nLeft),rgb((nRight-1)/nRight,((nRight-1):0)/nRight,((nRight-1):0)/nRight))
  #--------------------------------------------#

  #--------------------------------------------#
  Old.Group.Files <- Group.Files
  Medoid.Data     <- matrix(ncol=ncol(Data),nrow=length(Group.Files))
  rownames(Medoid.Data) <- Old.Group.Files
  colnames(Medoid.Data) <- colnames(Data)
  Group.Files <- c()
  for (Group.File in Old.Group.Files) {
    Group.IDs <- as.vector(read.table(Group.File,header=FALSE,sep="\t",quote="")[,1])
    if (length(Group.IDs) < 1) {
      next()
    } else if (length(Group.IDs) == 1) {
      Medoid.Data[Group.File,] <- Data[Group.IDs,]
    } else {
      Medoid.Data[Group.File,] <- apply(Data[Group.IDs,],2,median)
    }
  }
  z      <- (1+apply(Medoid.Data, 1, which.max)) %/% 2
  ddr    <- reorder( as.dendrogram( hclust(dist(Medoid.Data), method='ward') ), z )
  rowInd <- order.dendrogram(ddr)
  Group.Files <- Old.Group.Files[rowInd]
  #--------------------------------------------#

  #--------------------------------------------#
  Tmp.Data <- matrix(Scale[2:length(Scale)],nrow=1)
  Labels <- sprintf("] %s, %s ]",Scale[1:(length(Scale)-1)],Scale[2:length(Scale)])
  pdf(file=gsub(".pdf",".Colorscale.pdf",OUT.File), width=5, height=2, title="HeatMAP")
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  layout(matrix(c(1,0), 2, 1, byrow = TRUE), widths  = c(1), heights = c(1,3), respect = FALSE)
  image(1:ncol(Tmp.Data),1,t(Tmp.Data),xlab='',ylab='',axes=FALSE,cex=0.8,breaks=Scale,col=Colors)
  axis(1,at=1:ncol(Tmp.Data),labels=Labels,las=2,tick=FALSE)
  dev.off()
  #--------------------------------------------#

  #--------------------------------------------#
  nGroups <- length(Group.Files)
  N <- 0
  for (Group.File in Group.Files) {N<-N+length(as.vector(read.table(Group.File,header=FALSE,sep="\t",quote="")[,1])) }
  Blank.Data <- matrix(NA,ncol=ncol(Data),nrow=round(N/80))
    
  New.Data <- matrix(nrow=0,ncol=ncol(Data))
  colnames(New.Data) <- colnames(Data)
  for (i in 1:nGroups) {
    Group.File        <- Group.Files[i]
    Group.Individuals <- as.vector(read.table(Group.File,header=FALSE,sep="\t",quote="")[,1])

    if (length(Group.Individuals)>1) {
        Group.Data        <- Data[Group.Individuals,]
   	z                 <- (1+apply(Group.Data, 1, which.max)) %/% 2
   	ddr               <- reorder( as.dendrogram( hclust(dist(Group.Data), method='ward') ), z )
   	rowInd            <- order.dendrogram(ddr)
   	Group.Data        <- Group.Data[rowInd,]
    } else {
        Group.Data        <- matrix(Data[Group.Individuals,],nrow=1)
	rownames(Group.Data) <- Group.Individuals
	colnames(Group.Data) <- colnames(New.Data)
    }	
    New.Data          <- rbind(New.Data , Group.Data)
    if (i < nGroups) {
        New.Data          <- rbind(New.Data , Blank.Data)
    }
  }
  
  #ddc    <- reorder( as.dendrogram( hclust(dist(t(Data)))), 1:dim(t(Data))[2])
  #colInd <- order.dendrogram(ddc)
  colInd  <- khorder(t(Data),200)
	
  #New.Data <- New.Data[,colInd]
  New.Data <- New.Data[nrow(New.Data):1,]
  #--------------------------------------------#
    		
  #--------------------------------------------#
  pdf(file=OUT.File, width=7, height=7, title="HeatMAP")
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  layout(matrix(c(1,0,0,0), 2, 2, byrow = TRUE), widths  = c(7,1), heights = c(7,1), respect = FALSE)
  image(1:ncol(New.Data),1:nrow(New.Data),t(New.Data),xlab='',ylab='',axes=FALSE,cex=0.8,col=Colors,breaks=Scale)
  axis(1, at=1:ncol(New.Data), labels=colnames(New.Data), las=3,tick=FALSE,cex.axis=0.5)
  axis(4, at=1:nrow(New.Data), labels=rownames(New.Data), las=2,tick=FALSE,cex.axis=0.5)
  dev.off()
  #--------------------------------------------#

  #--------------------------------------------#
  png(filename = sprintf("%s.png",OUT.File), width = 5000, height = 3000)
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  layout(matrix(c(1,0,0,0), 2, 2, byrow = TRUE), widths  = c(7,1), heights = c(7,1), respect = FALSE)
  image(1:ncol(New.Data),1:nrow(New.Data),t(New.Data),xlab='',ylab='',axes=FALSE,cex=0.8,col=Colors,breaks=Scale)
  axis(1, at=1:ncol(New.Data), labels=colnames(New.Data), las=3,tick=FALSE,cex.axis=0.5)
  axis(4, at=1:nrow(New.Data), labels=rownames(New.Data), las=2,tick=FALSE,cex.axis=0.5)
  dev.off()
  #--------------------------------------------#


  return(OUT.File)
}
#--------------------------------------------#

#--------------------------------------------#
HeatMapFC <- function(Data,Group.Files,OUT.File) {
  #--------------------------------------------#
  Data <- Data[,which(apply(Data,2,sd)!=0)]
  #--------------------------------------------#
 
  #--------------------------------------------#
  Q.Val  <- quantile(Data,p=c(0,0.05,0.95,1),na.rm=TRUE)
  min.Val <- floor(Q.Val[1])
  max.Val <- ceiling(Q.Val[4])
  m.Val <- max(abs(c(min.Val,max.Val)))
    
  Scale <- extended_breaks(n=20)(floor(-1*m.Val):ceiling(m.Val))
  n <- (length(Scale)-1)/2
  nLeft  <- n
  nRight <- n
  Colors <- c(rgb((0:(nLeft-1))/nLeft,(0:(nLeft-1))/nLeft,(nLeft-1)/nLeft),rgb((nRight-1)/nRight,((nRight-1):0)/nRight,((nRight-1):0)/nRight))
  #--------------------------------------------#

  #--------------------------------------------#
  Old.Group.Files <- Group.Files
  Medoid.Data     <- matrix(ncol=ncol(Data),nrow=length(Group.Files))
  rownames(Medoid.Data) <- Old.Group.Files
  colnames(Medoid.Data) <- colnames(Data)
  Group.Files <- c()
  for (Group.File in Old.Group.Files) {
    Group.IDs <- as.vector(read.table(Group.File,header=FALSE,sep="\t",quote="")[,1])
    if (length(Group.IDs) < 1) {
      next()
    } else if (length(Group.IDs) == 1) {
      Medoid.Data[Group.File,] <- Data[Group.IDs,]
    } else {
      Medoid.Data[Group.File,] <- apply(Data[Group.IDs,],2,median)
    }
  }
  z      <- (1+apply(Medoid.Data, 1, which.max)) %/% 2
  ddr    <- reorder( as.dendrogram( hclust(dist(Medoid.Data), method='ward') ), z )
  rowInd <- order.dendrogram(ddr)
  Group.Files <- Old.Group.Files[rowInd]
  #--------------------------------------------#

  #--------------------------------------------#
  Tmp.Data <- matrix(Scale[2:length(Scale)],nrow=1)
  Labels <- sprintf("] %s, %s ]",Scale[1:(length(Scale)-1)],Scale[2:length(Scale)])
  pdf(file=gsub(".pdf",".Colorscale.pdf",OUT.File), width=5, height=2, title="HeatMAP")
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  layout(matrix(c(1,0), 2, 1, byrow = TRUE), widths  = c(1), heights = c(1,3), respect = FALSE)
  image(1:ncol(Tmp.Data),1,t(Tmp.Data),xlab='',ylab='',axes=FALSE,cex=0.8,breaks=Scale,col=Colors)
  axis(1,at=1:ncol(Tmp.Data),labels=Labels,las=2,tick=FALSE)
  dev.off()
  #--------------------------------------------#

  #--------------------------------------------#
  nGroups <- length(Group.Files)
  N <- 0
  for (Group.File in Group.Files) {N<-N+length(as.vector(read.table(Group.File,header=FALSE,sep="\t",quote="")[,1])) }
  Blank.Data <- matrix(NA,ncol=ncol(Data),nrow=round(N/120))
    
  New.Data <- matrix(nrow=0,ncol=ncol(Data))
  colnames(New.Data) <- colnames(Data)
  for (i in 1:nGroups) {
     Group.File        <- Group.Files[i]
     Group.Individuals <- as.vector(read.table(Group.File,header=FALSE,sep="\t",quote="")[,1])

     if (length(Group.Individuals)>1) {
        Group.Data        <- Data[Group.Individuals,]
   	z                 <- (1+apply(Group.Data, 1, which.max)) %/% 2
   	#ddr               <- reorder( as.dendrogram( hclust(dist(Group.Data), method='ward.D2') ), z )
   	ddr               <- reorder( as.dendrogram( hclust(dist(Group.Data), method='ward') ), z )
   	rowInd            <- order.dendrogram(ddr)
   	Group.Data        <- Group.Data[rowInd,]
     } else {
        Group.Data        <- matrix(Data[Group.Individuals,],nrow=1)
	rownames(Group.Data) <- Group.Individuals
	colnames(Group.Data) <- colnames(New.Data)
     }	
     New.Data          <- rbind(New.Data , Group.Data)
     if (i < nGroups) {
        New.Data          <- rbind(New.Data , Blank.Data)
     }
  }
  
  #ddc    <- reorder( as.dendrogram( hclust(dist(t(Data)))), 1:dim(t(Data))[2])
  #colInd <- order.dendrogram(ddc)
  colInd  <- khorder(t(Data),200)
	  	
  New.Data <- New.Data[,colInd]
  New.Data <- New.Data[nrow(New.Data):1,]
  #--------------------------------------------#
    		
  #--------------------------------------------#
  pdf(file=OUT.File, width=7, height=7, title="HeatMAP")
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  layout(matrix(c(1,0,0,0), 2, 2, byrow = TRUE), widths  = c(7,1), heights = c(7,1), respect = FALSE)
  image(1:ncol(New.Data),1:nrow(New.Data),t(New.Data),xlab='',ylab='',axes=FALSE,cex=0.8,col=Colors,breaks=Scale)
  axis(1, at=1:ncol(New.Data), labels=colnames(New.Data), las=3,tick=FALSE,cex.axis=0.5)
  axis(4, at=1:nrow(New.Data), labels=rownames(New.Data), las=2,tick=FALSE,cex.axis=0.5)
  dev.off()
  #--------------------------------------------#
 

  #--------------------------------------------#
  png(filename = sprintf("%s.png",OUT.File), width = 5000, height = 5000)
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  layout(matrix(c(1,0,0,0), 2, 2, byrow = TRUE), widths  = c(7,1), heights = c(7,1), respect = FALSE)
  image(1:ncol(New.Data),1:nrow(New.Data),t(New.Data),xlab='',ylab='',axes=FALSE,cex=0.8,col=Colors,breaks=Scale)
  axis(1, at=1:ncol(New.Data), labels=colnames(New.Data), las=3,tick=FALSE,cex.axis=0.5)
  axis(4, at=1:nrow(New.Data), labels=rownames(New.Data), las=2,tick=FALSE,cex.axis=0.5)
  dev.off()
  #--------------------------------------------#

  return(OUT.File)
}
#--------------------------------------------#

khorder <- function(Data,n) {
   if (nrow(Data) < n) {
      return(rownames(Data)[hclust(dist(Data), method='ward')$order])
   } 
   if (sum(apply(Data,2,sd)) <= 0) {
      return(rownames(Data))
   }  

   k <- 2
   cat("Split dataset into",k,"clusters (k-means)\n")
   k.res <- kmeans(Data, k, iter.max = 1000 , nstart = 1 ,algorithm='Hartigan-Wong')

   indexes <- c()
   for (ki in 1:k) {
       k.indexes   <- rownames(Data)[which(k.res$cluster == ki)]
	  
       if (length(k.indexes) > n) {
          cat("Sort cluster",ki,"->",length(k.indexes),"... run khorder again!\n")
       	  indexes <- c(indexes , khorder(Data[k.indexes,],n))
       } else if (length(k.indexes) > 2) {      
          cat("Sort cluster",ki,"->",length(k.indexes),"... run hclust!\n")
          indexes <- c(indexes , k.indexes[hclust(dist(Data[k.indexes,]), method='ward')$order])
       } else {
          indexes <- c(indexes , k.indexes)
       }
    }
    return(indexes)
}
