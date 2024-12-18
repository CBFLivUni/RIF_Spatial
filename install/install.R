options(Ncpus = 6)


#function to install libararies from bioconductor
installPackage <- function(libName) {
  if(libName %in% rownames(installed.packages()) == FALSE){
    BiocManager::install(libName,ask = FALSE)
  }}

install.packages("BiocManager")

#read the libraries needed
packagesToInstall <- read.delim("install/packagesToInstall.txt",header=FALSE, 
                                stringsAsFactors = FALSE)

#install all the libraries
sapply(packagesToInstall[,1],installPackage)

devtools::install_github("Nanostring-Biostats/SpatialOmicsOverlay")
devtools::install_github("mojaveazure/seurat-disk")
devtools::install_github("knncolle/BiocNeighbors")
devtools::install_github('saezlab/liana')

library(ExperimentHub)
eh <- ExperimentHub()
db <- eh[["EH3226"]]
eh[["EH3226"]]
