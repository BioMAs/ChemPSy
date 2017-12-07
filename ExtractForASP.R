setwd("/home/genouest/irset/archives/projects/ChemPSy/processed_data/Rattus_norvegicus/ALL/4_pca_analysis")
load("/home/genouest/irset/archives/projects/ChemPSy/processed_data/Rattus_norvegicus/ALL/4_pca_analysis/PCA.RData")
#setwd("Z:/projects/scripts/ChemPsy/ALL/3_reduce_data")
#load("Z:/projects/scripts/ChemPsy/ALL/3_reduce_data/data.reduced.RData")
ls()


######################
# lecture et ménage pour garder de l'espace mémoire
######################
rm(Data.log2fc)
rm(Data.pBH)
rm(Reduced.Data.log2fc)
rm(Reduced.Data.pBH)

all_cns = colnames(Reduced.Data.bin)
head(all_cns)

splitted_all_cns = unlist(strsplit(all_cns,"[+]"))
mat_cns = matrix(data=splitted_all_cns, ncol=6, byrow=T)
colnames(mat_cns) = c("Bank","Tissue","Molecule","Gen","Dose","exp")

head(mat_cns)


######################
# nombre de relations à exporter
######################
sum(Reduced.Data.bin == 1) + sum(Reduced.Data.bin == -1)


######################
# uniformisation des durées d'exposition
######################
ExpoTime = unique(mat_cns[,6])
unique(ExpoTime)

ExpoTime = ifelse((mat_cns[,6]=="1_d" | mat_cns[,6]=="24_hr" | mat_cns[,6]=="24h"),"h24",
                        ifelse(mat_cns[,6]=="3_d","h72", ifelse(mat_cns[,6]=="5_d","h120",
                                      ifelse(mat_cns[,6]=="7_d","h168", ifelse(mat_cns[,6]=="3_hr","h3",
                                                    ifelse(mat_cns[,6]=="15_day","h360", ifelse(mat_cns[,6]=="29_day","h696",
                                                                  ifelse((mat_cns[,6]=="4_day" | mat_cns[,6]=="4_d"),"h96",
                                                                         ifelse((mat_cns[,6]=="6_hr" | mat_cns[,6]=="6h" | mat_cns[,6]==".25_d"),"h6",
                                                                                ifelse(mat_cns[,6]=="8_day","h192", ifelse(mat_cns[,6]=="9_hr","h9",
                                                                                              ifelse(mat_cns[,6]=="14_d","h336", ifelse(mat_cns[,6]=="2_hr","h2",
                                                                                                            ifelse(mat_cns[,6]=="8_hr","h8", ifelse(mat_cns[,6]=="10d","h240",
                                                                                                                          ifelse(mat_cns[,6]=="2d","h48","na"))))))))))))))))


mat_cns2 = cbind(mat_cns[,-6],ExpoTime)
head(mat_cns2)


######################
# création des concepts
######################
mypaste  <- function (x,y) paste(x,y,sep='","')
redpaste <- function (x)   paste('c("',Reduce(mypaste,x),'")',sep="")
names <- apply(FUN=redpaste, mat_cns2, MARGIN=1)
head(names)

genesnames <- rownames(Reduced.Data.bin)
head(genesnames)


######################
# export format ASP
######################
for (i in 1:ncol(Reduced.Data.bin)){
  if (i %% 100 == 1) print(i)
  genesup   <- which(Reduced.Data.bin[,i]==1)
  genesdo   <- which(Reduced.Data.bin[,i]==-1)
  if ((length(genesup) > 0) & (length(genesdo) > 0)){
    linesup    <- function(x) paste('rel(g',genesnames[x],',',names[i],',1)',sep="")
    linesdo    <- function(x) paste('rel(g',genesnames[x],',',names[i],',-1)',sep="")
    lines = c(linesup(genesup),linesdo(genesdo))
    write.table(lines,file="20161013_AllOrgans_ASP.txt",quote=F,row.names=F,col.names=F,append=T)
  }
  else if ((length(genesup) > 0)){
    linesup    <- function(x) paste('rel(g',genesnames[x],',',names[i],',1)',sep="")
    write.table(linesup(genesup),file="20161013_AllOrgans_ASP.txt",quote=F,row.names=F,col.names=F,append=T)
  }
  else if ((length(genesdo) > 0)){
    linesdo    <- function(x) paste('rel(g',genesnames[x],',',names[i],',-1)',sep="")
    write.table(linesdo(genesdo),file="20161013_AllOrgans_ASP.txt",quote=F,row.names=F,col.names=F,append=T)
  }
}



