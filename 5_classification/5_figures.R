#--------------------------------------------#
# DrugMatrix downloaded at ftp://anonftp.niehs.nih.gov/drugmatrix
# folder organization is similar
#--------------------------------------------#

#--------------------------------------------#
RLIBS<-"/home/genouest/irset/archives/softs/R/3.2.3/libs"
setwd("/home/genouest/irset/archives/projects/ChemPSy")
.libPaths(RLIBS)
source("scripts/heatmap.R")
#--------------------------------------------#

#--------------------------------------------#
Input.Dir      <- commandArgs(TRUE)[1]
PCA.RData.File <- commandArgs(TRUE)[2]
Output.Dir     <- commandArgs(TRUE)[3]

#Input.Dir       <- "test_all_recmclust"
#PCA.RData.File  <- "processed_data/Rattus_norvegicus/ALL/3_reduce_data/4_pca_analysis/PCA.RData"
#Output.Dir      <- "test"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
load(PCA.RData.File)
#--------------------------------------------#

#--------------------------------------------#
All.Group.Files   <- list.files(path=Input.Dir,pattern=".txt")
All.Group.Files   <- All.Group.Files[grep("_C[0-9]+\\.[0-9]+\\.txt$",All.Group.Files)]
All.Group.Methods <- unique(gsub("[0-9]+\\.txt$","",All.Group.Files))
#--------------------------------------------#

#--------------------------------------------#
for (Group.Method in All.Group.Methods) {
    Group.Files <- list.files(path=Input.Dir , pattern=Group.Method)
    Group.Files <- Group.Files[grep(sprintf("^%s[0-9]+\\.txt$",Group.Method),Group.Files)]
    Group.Files <- sprintf("%s/%s",Input.Dir,Group.Files)

    cat(Group.Method,"-> Heatmap -> bin\n")
    Out.File    <- sprintf("%s/%sHeatmap_bin.pdf",Output.Dir,Group.Method )
    try(HeatMapFC(t(Reduced.Data.bin),Group.Files,Out.File),silent=TRUE)

    cat(Group.Method,"-> Heatmap -> log2fc\n")
    Out.File    <- sprintf("%s/%sHeatmap_log2fc.pdf",Output.Dir,Group.Method )
    try(HeatMapFC(t(Reduced.Data.log2fc),Group.Files,Out.File),silent=TRUE)
}
#--------------------------------------------#

quit()
