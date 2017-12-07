print("******************* Chargement de l'envirronnement ******************")

setwd("/home/genouest/irset/archives/projects/ChemPSy/processed_data/Rattus_norvegicus/ALL/4_pca_analysis")
load("/home/genouest/irset/archives/projects/ChemPSy/processed_data/Rattus_norvegicus/ALL/4_pca_analysis/PCA.RData")

#setwd("Z:/projects/scripts/ChemPsy/ALL/")
#load("Z:/projects/scripts/ChemPsy/ALL/smallPCA.RData")
ls()


######################
# lecture 
######################

print("******************* Lecture Fichier ******************")
M = t(Reduced.Data.bin.pca.res$ind$coord)
C = colnames(M)
N = ncol(M)
head(C)

######################
# calcul de la corrélation 
######################

print("******************* Calcul de la Corrélation ******************")
corXY     <- function(x) cor(x,M)
head(corXY(M[,1]))
Mat_corXY <- apply(FUN=corXY, M, MARGIN=2)
Mat_corXY[1:5,1:5]



######################
# Nombre d'arêtes du graphe avec seuil choisi 
######################

print("******************* Nombre d'arêtes du graphe ******************")
quantile(Mat_corXY, seq(0.8,1,by=0.005))

S = quantile(Mat_corXY, 0.995) # top 0.5 percent = 0.5049595  
print((sum(Mat_corXY>S)-N)/2)

d = ((sum(Mat_corXY>S)-N)/2) / (N*(N-1)/2) 

######################
# Ecriture du graphe 
######################

print("******************* Ecriture du graphe ******************")
all_lines   <- c((sum(Mat_corXY>S)-N)/2)

for (i in 1:(N-1)){
  if (i %% 100 == 1) print(i)
  
  correlated   <- which(Mat_corXY[(i+1):N,i]>S)+i
  if (length(correlated)>0){
    compound   <- C[i]
    lines      <- function(x) paste(compound,"\t",C[x],"\t",Mat_corXY[i,x],sep="")
    all_lines  <- c(all_lines,lines(correlated))
  }
}
print (all_lines[1:50])
write.table(all_lines,file="20161018_AllOrgans_corXY.gr",quote=F,row.names=F,col.names=F,append=F)

rm(compound)
rm(lines)
rm(correlated)
rm(all_lines)
rm(corXY)
rm(S)
rm(i)
rm(d)




######################
# calcul de distance euclidienne
######################

print("******************* Calcul de la Distance Euclidienne (log(1+de)) ******************")
Mat_eucXY     <- log(1+as.matrix(dist(t(M),method="euclidean",diag=F,upper=F)))
Mat_eucXY[1:50,1:2]



######################
# Nombre d'arêtes du graphe avec seuil choisi 
######################

print("******************* Nombre d'arêtes du graphe ******************")

quantile(Mat_eucXY, seq(0,0.2,by=0.005))
S = quantile(Mat_eucXY, 0.005) # top 0.5 percent =  1.814912
print((sum(Mat_eucXY<S)-N)/2)

UPPER.LIMIT <- min(Mat_eucXY[Mat_eucXY != 0]) + 1
UPPER.LIMIT

######################
# Ecriture du graphe 
######################

print("******************* Ecriture du graphe ******************")
all_lines   <- c((sum(Mat_eucXY<S)-N)/2)


for (i in 1:(N-1)){
  if (i %% 100 == 1) print(i)
  
  close        <- which(Mat_eucXY[(i+1):N,i]<S)+i
  if (length(close)>0){
    compound   <- C[i]
    lines      <- function(x) paste(compound,"\t",C[x],"\t",UPPER.LIMIT-Mat_eucXY[i,x],sep="")
    all_lines  <- c(all_lines,lines(close))
  }
}
print (all_lines[1:50])
write.table(all_lines,file="20161018_AllOrgans_eucXY.gr",quote=F,row.names=F,col.names=F,append=F)

rm(compound)
rm(lines)
rm(close)
rm(all_lines)
rm(UPPER.LIMIT)
rm(S)
rm(i)





######################
# uniformisation des durées d'exposition
######################

print("******************* Uniformisation des durées d'exposition ******************")

splitted_all_cns = unlist(strsplit(C,"[+]"))
MC = matrix(data=splitted_all_cns, ncol=6, byrow=T)
colnames(MC) = c("Bank","Tissue","Molecule","Gen","Dose","exp")
head(MC)


expoTimeNorm <- function(x) { switch(x, 
                                     "1_d"    = {"h24"},
                                     "24_hr"  = {"h24"},
                                     "24h"    = {"h24"},
                                     "2d"     = {"h48"},
                                     "3_d"    = {"h72"},
                                     "4_day"  = {"h96"},
                                     "4_d"    = {"h96"},
                                     "5_d"    = {"h120"},
                                     "7_d"    = {"h168"},
                                     "8_day"  = {"h192"},
                                     "10d"    = {"h240"},
                                     "14_d"   = {"h336"},
                                     "15_day" = {"h360"},
                                     "29_day" = {"h696"},
                                     "2_hr"   = {"h2"},
                                     "3_hr"   = {"h3"},
                                     "6_hr"   = {"h6"},
                                     "6h"     = {"h6"},
                                     ".25_d"  = {"h6"},
                                     "8_hr"   = {"h8"},
                                     "9_hr"   = {"h9"},
                                     "NA"     = {"NA"},                                     
                                      {"NA"})
}


Exp=sapply(MC[,6], FUN=expoTimeNorm)
head(Exp)

rm(splitted_all_cns)
rm(MC)
rm(expoTimeNorm)

IDXT = which(Exp=="h24" | Exp=="h48"| Exp=="h72")
Ni   = length(IDXT)



######################
# Nombre d'arêtes du graphe avec seuil choisi 
######################

print("******************* Nombre d'arêtes du graphe ******************")
quantile(Mat_corXY[IDXT,IDXT], seq(0.8,1,by=0.005))

S = quantile(Mat_corXY[IDXT,IDXT], 0.99) # top 0.5 percent = 0.4964722, top 1 percent = 0.4198243  
print((sum(Mat_corXY[IDXT,IDXT]>S)-Ni)/2)

d = ((sum(Mat_corXY[IDXT,IDXT]>S)-Ni)/2) / (Ni*(Ni-1)/2) 


######################
# Ecriture du graphe 
######################

print("******************* Ecriture du graphe ******************")
all_lines   <- c((sum(Mat_corXY[IDXT,IDXT]>S)-Ni)/2)

for (i in 1:(N-1)){
  if (i %% 100 == 1) print(i)
  
  if (i %in% IDXT){
    correlatedi <- which(Mat_corXY[(i+1):N,i]>S)+i
    correlated   <- IDXT[IDXT %in% correlatedi]

    if (length(correlated)>0){
      compound   <- C[i]
      lines      <- function(x) paste(compound,"\t",C[x],"\t",Mat_corXY[i,x],sep="")
      all_lines  <- c(all_lines,lines(correlated))
    }
  }
}
print (all_lines[1:50])
write.table(all_lines,file="20161018_AllOrgans123_corXY.gr",quote=F,row.names=F,col.names=F,append=F)

rm(compound)
rm(lines)
rm(correlated)
rm(correlatedi)
rm(all_lines)
rm(S)
rm(i)




######################
# Nombre d'arêtes du graphe avec seuil choisi 
######################

print("******************* Nombre d'arêtes du graphe ******************")

quantile(Mat_eucXY[IDXT,IDXT], seq(0,0.2,by=0.005))
S = quantile(Mat_eucXY[IDXT,IDXT], 0.01) # top 0.5 percent =  1.874096, top 1 percent = 2.020678
print((sum(Mat_eucXY[IDXT,IDXT]<S)-Ni)/2)

d = ((sum(Mat_eucXY[IDXT,IDXT]<S)-Ni)/2) / (Ni*(Ni-1)/2)

UPPER.LIMIT <- min(Mat_eucXY[Mat_eucXY != 0]) + 1
UPPER.LIMIT

 



######################
# Ecriture du graphe 
######################

print("******************* Ecriture du graphe ******************")
all_lines   <- c((sum(Mat_eucXY[IDXT,IDXT]<S)-Ni)/2)


for (i in 1:(N-1)){
  if (i %% 100 == 1) print(i)
  if (i %in% IDXT){
    closei    <- which(Mat_eucXY[(i+1):N,i]<S)+i
    close     <- IDXT[IDXT %in% closei]  
    if (length(close)>0){
      compound   <- C[i]
      lines      <- function(x) paste(compound,"\t",C[x],"\t",UPPER.LIMIT-Mat_eucXY[i,x],sep="")
      all_lines  <- c(all_lines,lines(close))
    }
  }
}
print (all_lines[1:50])
write.table(all_lines,file="20161018_AllOrgans123_eucXY.gr",quote=F,row.names=F,col.names=F,append=F)

rm(compound)
rm(lines)
rm(close)
rm(closei)
rm(all_lines)
rm(UPPER.LIMIT)
rm(S)
rm(i)


