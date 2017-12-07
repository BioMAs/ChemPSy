#--------------------------------------------#
# DrugMatrix downloaded at ftp://anonftp.niehs.nih.gov/drugmatrix
# folder organization is similar
#--------------------------------------------#

#--------------------------------------------#
setwd("/home/genouest/irset/archives/projects/ChemPSy")
.libPaths("scripts/RLIBRARIES")
library(irr)
#--------------------------------------------#

#--------------------------------------------#
PCA.RData.File <- commandArgs(TRUE)[1]
Output.Dir <- commandArgs(TRUE)[2]

if (0==1) {
   PCA.RData.File   <- "processed_data/Rattus_norvegicus/ALL/4_pca_analysis/PCA.RData"
   Output.Dir       <- "test"
}
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
Condordance.coefficient <- function(x,y,weights) {
   n11 <- sum(weights[which(x==y & x!=0)])
   return(n11)
}

Similarity.coefficients <- function(x,y,weights) {
   # Adapted to -1, 0, 1 matrix
   # if x=1 and y=-1, n10 is incremented
   # if y=1 and x=-1, n01 is incremented
   # n11, is incremented when x=y=1 or x=y=-1
   # n00, is incremented when x=y=0
   
   n11 <- sum(weights[which(x==y & x!=0)])
   n00 <- sum(weights[which(x==y & x==0)])
   n10 <- sum(weights[which((x==1 & x>y) | (x==-1 & y==0))])
   n01 <- sum(weights[which((y==1 & x<y) | (y==-1 & x==0))])

   r <- n11
   R <- sum(weights[which(x!=0)])
   n <- sum(weights[which(y!=0)])
   N <- sum(weights)

   ZS = (r-n*R/N) / sqrt( (n*R/N) * (1-R/N) * (1-(n-1)/(N-1)))
   HT = phyper(r-1,n,N-n,R,lower.tail=FALSE)

   CT = cor.test(x,y,method="pearson",alternative ="greater")
   CC = CT$estimate
   CP = CT$p.value
   i  = which( x!=0 | y!=0)
   if (length(i) > 2) {
      CT2= cor.test(x[i],y[i],method="pearson",alternative ="greater")
      CC2 = CT2$estimate
      CP2 = CT2$p.value
   } else {
      CT2=CT
      CC2=CC
      CP2=CP
   }
   ED = sqrt(sum((x[i]-y[i])^2))
   OC = n11 / min( sum(weights[which(x!=0)]) , sum(weights[which(y!=0)]))
   OR = (n11*n00) / (n01*n10)
   JI = n11 / (n01+n10+n11)
   JD = 1-JI
   HS = (n11*n00 - n10*n01) / sqrt((n11+n10)*(n01+n00)*(n11+n01)*(n10+n00))
  
   m <- matrix(c(x,y),ncol=2,byrow=FALSE)
   KU = kappa2(m,"unweighted")$value
   KE = kappa2(m,"equal")$value
   KS = kappa2(m,"squared")$value


  return( list(Condordance=n11, Overlap.coefficient=OC, Odds.ratio=OR, Jaccard.index=JI, Jaccard.distance=JD, HubertGamma.score=HS, Euclidean.distance=ED, Correlation.coefficient=CC, Correlation.pvalue=CP, Correlation.coefficient2=CC2, Correlation.pvalue2=CP2, Z.score=ZS, Hypergeometric.test=HT, Kappa.unweighted=KU, Kappa.equal=KE, Kappa.squared=KS) )
}
#--------------------------------------------#

#--------------------------------------------#
Data <- t(Reduced.Data.bin)
weights1 <- rep(1,ncol(Data))
weights2 <- 1/apply(abs(Data),2,sum)
weights2[which(apply(abs(Data),2,sum) == 0)] <- 0
#--------------------------------------------#

#--------------------------------------------#
system.time({
  N <- nrow(Data)
  IDs <- rownames(Data)

   Similarities.res     <- matrix(NA,ncol=2+16,nrow=K)
   colnames(Similarities.res) <- c("Chemical1","Chemical2","Concordance","Overlap.coefficient","Odds.ratio","Jaccard.index","Jaccard.distance","HuberGamma.score","Euclidean.distance","Correlation.coefficient","Correlation.pvalue","Correlation.coefficient.no00","Correlation.pvalue.no00","Z.score","Hypergeometric.test","Kappa.unweighted","Kappa.equal","Kappa.squared")
   adjusted.Similarities.res <- Similarities.res

   k <- 0
   for (i in 1:(N-1)) {
       for (j in (i+1):N) {
       	   k <- k+1
	   
	   CC1 <- Condordance.coefficient(Data[i,],Data[j,],weights1)
	   if (CC1 < 10) {
	      Similarities.res[k,1:3]          <- c(IDs[i],IDs[j],CC1)
	      adjusted.Similarities.res[k,1:3] <- c(IDs[i],IDs[j],Condordance.coefficient(Data[i,],Data[j,],weights2))
	   } else {
	      Similarities.res[k,]          <- c(IDs[i],IDs[j],unlist(Similarity.coefficients(Data[i,],Data[j,],weights1)))
              adjusted.Similarities.res[k,] <- c(IDs[i],IDs[j],unlist(Similarity.coefficients(Data[i,],Data[j,],weights2)))
	   }
       }
   }
})
indexes1 <- which(as.double(adjusted.Similarities.res[,"HuberGamma.score"]) != "NA")
indexes2 <- which(as.double(adjusted.Similarities.res[,"HuberGamma.score"]) == "NA")
indexes1 <- indexes1[sort(as.double(adjusted.Similarities.res[indexes1,"HuberGamma.score"]),decreasing=TRUE,index.return=TRUE)$ix]
indexes  <- c(indexes1,indexes2)
Similarities.res                          <- Similarities.res[indexes,]
adjusted.Similarities.res                 <- adjusted.Similarities.res[indexes,]

indexes                                   <- which(as.double(Similarities.res[,3]) >= 10)
Reduced.Similarities.res                  <- Similarities.res[indexes,]
Reduced.adjusted.Similarities.res         <- adjusted.Similarities.res[indexes,]

indexes                                   <- 1:round(nrow(Reduced.Similarities.res)*0.05)
Top5.Reduced.Similarities.res             <- Reduced.Similarities.res[indexes,]
Top5.Reduced.adjusted.Similarities.res    <- Reduced.adjusted.Similarities.res[indexes,]
indexes                                   <- 1:round(nrow(Reduced.Similarities.res)*0.1)
Top10.Reduced.Similarities.res            <- Reduced.Similarities.res[indexes,]
Top10.Reduced.adjusted.Similarities.res   <- Reduced.adjusted.Similarities.res[indexes,]

write.table(Similarities.res                       ,file = sprintf("%s/similarities_table.tsv"                       ,Output.Dir), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)
write.table(adjusted.Similarities.res              ,file = sprintf("%s/similarities_table.adjusted.tsv"              ,Output.Dir), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)
write.table(Reduced.Similarities.res               ,file = sprintf("%s/similarities_table.reduced.tsv"               ,Output.Dir), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)
write.table(Reduced.adjusted.Similarities.res      ,file = sprintf("%s/similarities_table.adjusted.reduced.tsv"      ,Output.Dir), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)
write.table(Top5.Reduced.Similarities.res          ,file = sprintf("%s/similarities_table.reduced.top5.tsv"          ,Output.Dir), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)
write.table(Top5.Reduced.adjusted.Similarities.res ,file = sprintf("%s/similarities_table.adjusted.reduced.top5.tsv" ,Output.Dir), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)
write.table(Top10.Reduced.Similarities.res         ,file = sprintf("%s/similarities_table.reduced.top10.tsv"         ,Output.Dir), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)
write.table(Top10.Reduced.adjusted.Similarities.res,file = sprintf("%s/similarities_table.adjusted.reduced.top10.tsv",Output.Dir), sep = '\t' , col.names = TRUE , row.names = FALSE , quote=FALSE)

#--------------------------------------------#
save(list=ls()[grep(".Similarities.res",ls())], file = sprintf("%s/SimilaritiesCoeff4Topology.RData",Output.Dir))
#--------------------------------------------#

quit()
