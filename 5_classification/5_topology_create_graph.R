#--------------------------------------------#
# DrugMatrix downloaded at ftp://anonftp.niehs.nih.gov/drugmatrix
# folder organization is similar
#--------------------------------------------#

#--------------------------------------------#
setwd("/home/genouest/irset/archives/projects/ChemPSy")
.libPaths("scripts/RLIBRARIES")
#--------------------------------------------#

#--------------------------------------------#
SimilaritiesCoeff.RData.File <- commandArgs(TRUE)[1]
Output.Dir                   <- commandArgs(TRUE)[2]

#SimilaritiesCoeff.RData.File   <- "test/SimilaritiesCoeff4Topology.RData"
#Output.Dir       <- "test"
#--------------------------------------------#

#--------------------------------------------#
dir.create(c(Output.Dir),recursive=TRUE,showWarnings=FALSE)
#--------------------------------------------#

#--------------------------------------------#
source("http://bioconductor.org/biocLite.R")
if("igraph" %in% rownames(installed.packages()) == FALSE) {install.packages("igraph", lib="scripts/RLIBRARIES",dependencies=TRUE)}
library(igraph)
load(SimilaritiesCoeff.RData.File)
#--------------------------------------------#


names <- as.character(unique(c(data\[,1\],data\[,2\])))

ids <- 1:length(names)"
names(ids) <- as.character(names)"

#--- GRAPH CREATION ------------------------------------------
#--- ADD VERTICES --------------------------------------------
g <- graph.empty(directed=FALSE)
g <- add.vertices(g, length(ids), name=as.character(names))
#--- ADD EDGES -----------------------------------------------
from <- as.character(data[,1])
to <- as.character(data[,2])
weight <- as.integer(data[,3])
edges <- matrix(c(ids[from],ids[to]),ncol=2)
g <- add.edges(g, t(edges))
E(g)$weight <- weight
#-------------------------------------------------------------
quit()
