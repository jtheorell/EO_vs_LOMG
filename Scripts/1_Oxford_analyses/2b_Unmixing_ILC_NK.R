library(flowCore)
library(flowSpecs)

refSet <- read.flowSet(path = "../External/Oxford/FCS/Raw/ILC_NK_beads",
                       transformation = FALSE)

#Now, we are going to create one crazy unmixing matrix, with all the 42 single-stains
#that we have, and that is then going to be picked apart into two matrices outside
#of R. 
#We have thymus cells too for autofluorescence correction, if needed in the 
#future, but as the main focus is on the PBMC, we will use them primarily. 
rawMat <- specMatCalc(refSet, groupNames = c("Beads_", "Dead_"), 
                      autoFluoName = "PBMC_unstained.../External/Oxford/FCS")
dir.create("Data/Oxford/Unmixing_matrices")
write.csv(rawMat, "Data/Oxford/Unmixing_matrices/Raw_mat_ILC_NK.csv")
#After export, the order of the rows and the naming is changed a bit. 

usedMat <- read.csv("Data/Oxford/Unmixing_matrices/Used_mat_ILC_NK.csv", row.names = 1,
                  check.names = FALSE)

#Now, we will run through this file-by file, as each file is so big. 
dir.create("../External/Oxford/FCS/Unmixed/ILC_NK", recursive = TRUE)
fileList <- list.files("../External/Oxford/FCS/Raw/ILC_NK_panel", full.names = TRUE)
for(i in fileList){
  print(i)
  locFF <- read.../External/Oxford/FCS(i, transformation = FALSE)
  unmixFF <- specUnmix(locFF, usedMat)
  #Here, we are making a change to the internals of the file, to avoid tautologies
  #in the naming 
  unmixFF@parameters@data$name <- gsub(".+_|", "", colnames(unmixFF@exprs))
  unmixFF@parameters@data$desc <- gsub("|_.+", "", colnames(unmixFF@exprs))
  write.../External/Oxford/FCS(unmixFF, paste0("../External/Oxford/FCS/Unmixed/ILC_NK/Unmixed_", 
                            gsub("../External/Oxford/FCS/Raw/ILC_NK_panel/", "", i)))
}


