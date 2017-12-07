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
# Custom CDF file linked to Entrez Gene IDs must be downloaded from the BRAINARRAY website: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
CTDChemicalIDs2Diseases.File <- commandArgs(TRUE)[1]
CTDHGIDs2Diseases.File       <- commandArgs(TRUE)[2]
HomoloGene4R.File            <- commandArgs(TRUE)[3]
Chemicals2CTD.File           <- commandArgs(TRUE)[4]
Output.Dir                   <- commandArgs(TRUE)[5]

if (0 == 1) {
  CTDChemicalIDs2Diseases.File <- "CTD_data/CTD_chemical2CTD_diseases.txt"
  CTDHGIDs2Diseases.File       <- "CTD_data/CTD_HomoloGene2CTD_diseases.txt"
  HomoloGene4R.File            <- "CTD_data/HomoloGene4R.txt"
  Chemicals2CTD.File           <- "CTD_data/ChemPSy_MESH.tsv"
  Output.Dir                   <- "CTD_data"
}
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
CTD.Chemical2Disease.Data      <- read.table(CTDChemicalIDs2Diseases.File,header=TRUE,row.names=1,sep="\t",quote="")
CTD.HomoloGene2Disease.Data    <- read.table(CTDHGIDs2Diseases.File,header=TRUE,row.names=1,sep="\t",quote="")
Gene2HomoloGene.Data           <- read.table(HomoloGene4R.File,header=FALSE,row.names=1,sep="\t",quote="")
CTD.Chemicals.Data             <- read.table(Chemicals2CTD.File,header=FALSE,row.names=1,sep="\t",quote="")
#--------------------------------------------#

save(CTD.Chemical2Disease.Data, CTD.HomoloGene2Disease.Data, Gene2HomoloGene.Data, CTD.Chemicals.Data, file = sprintf("%s/CTD_diseases_data.RData",Output.Dir))

quit()
