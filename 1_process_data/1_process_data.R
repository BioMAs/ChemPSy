#--------------------------------------------#
# DrugMatrix downloaded at ftp://anonftp.niehs.nih.gov/drugmatrix
# folder organization is similar
#--------------------------------------------#

#--------------------------------------------#
RLIBS<-"/home/genouest/irset/archives/softs/R/3.1.0/libs"
setwd("/home/genouest/irset/archives/projects/ChemPSy")
.libPaths(RLIBS)
#--------------------------------------------#

#--------------------------------------------#
# Custom CDF file linked to Entrez Gene IDs must be downloaded from the BRAINARRAY website: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
Association.File <- commandArgs(TRUE)[1]
RawData.Dir      <- commandArgs(TRUE)[2]
CDF.RData.File   <- commandArgs(TRUE)[3]
Output.Dir       <- commandArgs(TRUE)[4]

#Association.File <- "Affymetrix_data/Rattus_norvegicus/HEART/GSE57800/Experimental_conditions/GSE57800+HEART+tacrolimus+F0+134_mgkg+5_d/treatment.info"
#RawData.Dir      <- "Affymetrix_data/Rattus_norvegicus/HEART/GSE57800/Individual_experiments"
#CDF.RData.File   <- "scripts/CDF/Rat2302_Rn_ENTREZG.cdf.RData"
#Output.Dir       <- "Testerror/GSE57800+HEART+tacrolimus+F0+134_mgkg+5_d"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
if("R.utils" %in% rownames(installed.packages()) == FALSE) {install.packages("R.utils", lib=RLIBS)}
if("preprocessCore" %in% rownames(installed.packages()) == FALSE) {install.packages("preprocessCore", lib=RLIBS)}
if("affy" %in% rownames(installed.packages()) == FALSE) { biocLite("affy")}
if("limma" %in% rownames(installed.packages()) == FALSE) { biocLite("limma")}
library(R.utils)
library(affy)
library(limma)
load(CDF.RData.File)
#--------------------------------------------#

#--------------------------------------------#
Association.Data <- as.matrix(read.table(Association.File,header=FALSE,sep="\t",quote=""))
Treated.GZ.CEL.Filenames <- unique(Association.Data[,1])
Control.GZ.CEL.Filenames <- unique(Association.Data[,2])
GZ.CEL.Filenames <- c(Treated.GZ.CEL.Filenames,Control.GZ.CEL.Filenames)
CEL.Filenames    <- gsub(".gz$"    ,"",GZ.CEL.Filenames)
Filenames        <- gsub(".CEL.gz$","",GZ.CEL.Filenames)
#--------------------------------------------#

#--- Uncompress CEL files -------------------#
for (GZ.CEL.Filename in GZ.CEL.Filenames) {
	GZ.CEL.File <- sprintf("%s/%s",RawData.Dir,GZ.CEL.Filename)
	CEL.File    <- gsub(".gz$","",sprintf("%s/%s",Output.Dir,GZ.CEL.Filename))
	gunzip(filename=GZ.CEL.File,destname=CEL.File,remove=FALSE,overwrite=TRUE,ext="gz",FUN=gzfile)
}
#--------------------------------------------#

#--- Read raw data --------------------------#
# Read in the CEL files in the directory, then normalize the data
Affy.Batch    <- ReadAffy(filenames=sprintf("%s/%s",Output.Dir,CEL.Filenames),sampleNames=GZ.CEL.Filenames)
Affy.Batch@cdfName <- 'CDF.Env'
#--------------------------------------------#

#--- RMA normalization ----------------------#
RMA.Data <- exprs(rma(Affy.Batch))
rownames(RMA.Data) <- gsub("_at$","",rownames(RMA.Data))
RMA.Data <- RMA.Data[!grepl("^AFFX",rownames(RMA.Data)),]
# Save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.table(RMA.Data, file=sprintf("%s/normdata.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
# Finally, remove the CEL files
file.remove(sprintf("%s/%s",Output.Dir,CEL.Filenames))
#--------------------------------------------#

#--- Quality control before pre-processing --#
SampleColors=c(rep("green",length(Treated.GZ.CEL.Filenames)),rep("grey",length(Control.GZ.CEL.Filenames)))
pdf(file=sprintf("%s/qc_boxplot_beforenormalization.pdf",Output.Dir), width=7, height=7, title='Boxplot_Insentities')
par(mar = c(10, 5, 5, 1))
boxplot(Affy.Batch,names=FALSE,col=SampleColors,ylab='Log2 signal intensities')
title(main=sprintf('%s\nMedian:%s, mean:%s',basename(dirname(Association.File)),round(median(log2(exprs(Affy.Batch))),digit=3),round(mean(log2(exprs(Affy.Batch))),digit=3)),cex.main =1.0)
axis(1, at=1:length(GZ.CEL.Filenames), labels=GZ.CEL.Filenames, las=3, cex.axis = 0.8)
dev.off()

Expr.Colors <- rainbow(20)
for (i in 1:length(GZ.CEL.Filenames)) {
	GZ.CEL.Filename <- GZ.CEL.Filenames[i]
	Tmp.Data <- log2(matrix(exprs(Affy.Batch[,i]),nrow=nrow(Affy.Batch),ncol=ncol(Affy.Batch),byrow=TRUE))
	png(file=sprintf("%s/qc_image_%s.png",Output.Dir,GZ.CEL.Filename), width=600, height=600)
	image(1:nrow(Affy.Batch), 1:ncol(Affy.Batch), Tmp.Data, axes = FALSE, xlab = '', ylab = '',col=Expr.Colors)
	title(main=sprintf('Surface Intensity Distribution of the %s sample',GZ.CEL.Filename),cex.main = 0.8)
	dev.off()
}
#--------------------------------------------#

#--- Quality control after pre-processing ---#
SampleColors=c(rep("green",length(Treated.GZ.CEL.Filenames)),rep("grey",length(Control.GZ.CEL.Filenames)))
pdf(file=sprintf("%s/qc_boxplot_afternormalization.pdf",Output.Dir), width=7, height=7, title='Boxplot_Insentities')
par(mar = c(10, 5, 5, 1))
boxplot(as.data.frame(RMA.Data),names=FALSE,col=SampleColors,ylab='Log2 signal intensities')
title(main=sprintf('%s\nMedian:%s, mean:%s',basename(dirname(Association.File)),round(median(RMA.Data),digit=3),round(mean(RMA.Data),digit=3)),cex.main =1.0)
axis(1, at=1:ncol(RMA.Data), labels=colnames(RMA.Data), las=3, cex.axis = 0.8)
dev.off()

Distance.matrix    <- as.dist(1-cor(RMA.Data))
HC.Data            <- hclust(Distance.matrix)
Dendo.Data         <- reorder(as.dendrogram(HC.Data), 1:dim(RMA.Data)[2])
Column.Index       <- order.dendrogram(Dendo.Data)
Tmp.Data           <- RMA.Data[,Column.Index]
Correlation.matrix <- cor(Tmp.Data)

SampleColors       <- c(rep("green",length(Treated.GZ.CEL.Filenames)),rep("black",length(Control.GZ.CEL.Filenames)))[Column.Index]
pdf(file=sprintf("%s/qc_corrmatrix_afternormalization.pdf",Output.Dir), width=7, height=7, title='Correlation_Matrix')
layout(matrix(c(0,2,0,3,1,0,0,0,0), 3, 3, byrow = TRUE), widths  = c(1,7,2), heights = c(1,7,2), respect = FALSE)
par(mar = c(0, 0, 0, 0))
image(1:nrow(Correlation.matrix),1:ncol(Correlation.matrix),Correlation.matrix,xlab='',ylab='',col=gray(0:100/100),axes=FALSE,main=sprintf('%s\nCorrelation matrix',basename(dirname(Association.File))),cex=0.8)
Map(function(x,y,z) 
  axis(1,at=x,col.axis=y,labels=z,lwd=0,las=3,tick=FALSE),
  1:ncol(Correlation.matrix),
  SampleColors,
  colnames(Correlation.matrix)
)
Map(function(x,y,z) 
  axis(4,at=x,col.axis=y,labels=z,lwd=0,las=2,tick=FALSE),
  1:nrow(Correlation.matrix),
  SampleColors,
  rownames(Correlation.matrix)
)
plot(Dendo.Data, horiz = FALSE, axes = FALSE, xaxs = 'i', leaflab = 'none')
plot(Dendo.Data, horiz = TRUE , axes = FALSE, yaxs = 'i', leaflab = 'none')
dev.off()
#--------------------------------------------#

#--- Median data ----------------------------#
Med.RMA.Data <- matrix(ncol=2,nrow=nrow(RMA.Data))
colnames(Med.RMA.Data) <- c("Treated","Control")
rownames(Med.RMA.Data) <- rownames(RMA.Data)
Med.RMA.Data[,1] <- apply(RMA.Data[,Treated.GZ.CEL.Filenames,drop=FALSE],1,median)
Med.RMA.Data[,2] <- apply(RMA.Data[,Control.GZ.CEL.Filenames,drop=FALSE],1,median)
write.table(Med.RMA.Data, file=sprintf("%s/mednormdata.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
#--------------------------------------------#

#--- FC data ----------------------------#
FC.Data <- matrix(round(Med.RMA.Data[,1]-Med.RMA.Data[,2],digit=5),ncol=1,nrow=nrow(Med.RMA.Data))
colnames(FC.Data) <- c("FC")
rownames(FC.Data) <- rownames(Med.RMA.Data)
write.table(FC.Data, file=sprintf("%s/log2fcchangedata.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
#--------------------------------------------#

#--- Design matrix --------------------------#
Design.matrix <- matrix(0,ncol=2,nrow=ncol(RMA.Data))
colnames(Design.matrix) <- c("Treated","Control")
rownames(Design.matrix) <- colnames(RMA.Data)
Design.matrix[Treated.GZ.CEL.Filenames,1] <- 1
Design.matrix[Control.GZ.CEL.Filenames,2] <- 1
write.table(Design.matrix, file=sprintf("%s/designmatrix.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
#--------------------------------------------#

#--- Contrast matrix ------------------------#
Contrast.matrix <- matrix(c(1,-1),nrow=2,ncol=1)
colnames(Contrast.matrix) <- "Treated vs Control"
rownames(Contrast.matrix) <- c("Treated","Control")
write.table(Contrast.matrix, file=sprintf("%s/contrastmatrix.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
#--------------------------------------------#

#--------------------------------------------#
BG.Cutoff <- median(RMA.Data)
FC.Cutoff <- 1.5
PV.Cutoff <- 0.05
#--------------------------------------------#

#--------------------------------------------#
step1a      <- rownames(Med.RMA.Data)[which(apply(Med.RMA.Data,1,max)>=BG.Cutoff)]
step1b.up   <- rownames(FC.Data)[which(FC.Data[,1] >= log2(FC.Cutoff))]
step1b.down <- rownames(FC.Data)[which(FC.Data[,1] <= -1*log2(FC.Cutoff))]
step1b      <- union(step1b.up,step1b.down)
step2.up    <- intersect(step1a,step1b.up)
step2.down  <- intersect(step1a,step1b.down)
step2       <- union(step2.up,step2.down)
#--------------------------------------------#

#--------------------------------------------#
if (length(step2) > 0) {
   step2.data <- RMA.Data[step2,]	
   if (length(step2) == 1) {
      step2.data <- matrix(step2.data,nrow=1)
      colnames(step2.data) <- colnames(RMA.Data)
      rownames(step2.data) <- step2
   }

   # Fitting the linear model using the design matrix
   fit  <- lmFit(step2.data,design=Design.matrix)
   # Fitting the linear model for the contrasts
   fit2 <- contrasts.fit(fit,contrasts=Contrast.matrix)
   # Empirical Bayes statistics
   eb   <- eBayes(fit2)
   # F-values 
   F    <- eb$F.p.value
   FBH  <- p.adjust(F,method='BH')
   step3      <- step2[which(FBH <= PV.Cutoff)]
   step3.up   <- intersect(step3,step2.up)
   step3.down <- intersect(step3,step2.down)
} else {
   step3      <- c()
   step3.up   <- c()
   step3.down <- c()
   FBH        <- c()
   F          <- c()
}
#--------------------------------------------#

#--------------------------------------------#
res <- matrix(NA,ncol=7,nrow=nrow(RMA.Data))
rownames(res)                  <- rownames(RMA.Data)
colnames(res)                  <- c("step1a","step1b","step2","step3","log2FC","P","PBH")

res[      ,"step1a"]           <-  0
res[step1a,"step1a"]           <-  1

res[           ,"step1b"]      <-  0
res[step1b.up  ,"step1b"]      <-  1
res[step1b.down,"step1b"]      <- -1

res[step1a     ,"step2"]       <-  0
res[step2.up   ,"step2"]       <-  1
res[step2.down ,"step2"]       <- -1

res[step2      ,"step3"]       <-  0
res[step3.up   ,"step3"]       <-  1
res[step3.down ,"step3"]       <- -1

res[           ,"log2FC"   ] <-  FC.Data[,1]
res[step2      ,"P"    ] <-  F
res[step2      ,"PBH"] <-  FBH
write.table(res, file=sprintf("%s/filtration.txt",Output.Dir),quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)
#--------------------------------------------#

quit()
