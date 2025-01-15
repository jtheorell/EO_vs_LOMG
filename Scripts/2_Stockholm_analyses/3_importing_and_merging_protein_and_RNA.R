#Here, we take the 24 iterations from the totalVI run and average them, 
#whereafter they are imported into the common singleCellExperiment.

#We will import the first columnof each CSV file and average them, whereafter
#we will create one final frame. 
library(SingleCellExperiment)
library(readr)
library(Matrix)
library(scater)
library(DepecheR)

#First, we do this for the common proteins, that are extracted from the common
#protein model, to maximise the similarity between the observed and the approximated
#expression. Then, the same is done for all the other proteins. 
protNameList <- list(colnames(read.csv("../External/Stockholm/Data/Anndata_input/adt_counts_common.csv", row.names = 1)),
                     colnames(read.csv("../External/Stockholm/Data/Anndata_input/adt_counts.csv", row.names = 1)))


modelNameList <- list(list.files("../External/Stockholm/Data/Normalised_protein_models_common_proteins", 
                                 full.names = TRUE, 
                                 pattern = "Model_"),
                      list.files("../External/Stockholm/Data/Normalised_protein_models", 
                                 full.names = TRUE, 
                                 pattern = "Model_"))

nColList <- list(ncol(read_csv("../External/Stockholm/Data/Anndata_input/adt_counts_common.csv")),
                 ncol(read_csv("../External/Stockholm/Data/Anndata_input/adt_counts.csv")))

colNamesForAll <- row.names(read.csv("../External/Stockholm/Data/Anndata_input/adt_counts_common.csv", row.names = 1))

protDatList <- lapply(1:2, function(i){
  #Here, we have to remove the first column, as it contains the rownames
  #which are not supported as a separate entity in read_csv
  sumDat <- read_csv(modelNameList[[i]][[1]], col_select = 2:nColList[[i]],
                     show_col_types = FALSE)
  for(j in 2:length(modelNameList[[i]])){
    locDir <- modelNameList[[i]][[j]]
    locDat <- read_csv(locDir, col_select = 2:nColList[[i]],
                       show_col_types = FALSE)
    sumDat <- sumDat+locDat
  }
  averageProtData <- t(sumDat)/length(modelNameList[[i]])
  rownames(averageProtData) <- protNameList[[i]]
  colnames(averageProtData) <- colNamesForAll
  averageProtData
})

#Now, these are saved.  
write.csv(protDatList[[1]], "../External/Stockholm/Data/Normalised_protein_models_common_proteins/Final_average.csv")

write.csv(protDatList[[2]], "../External/Stockholm/Data/Normalised_protein_models/Final_average.csv")

#And now, we unify them, so that the approximations for the common proteins come
#from the common protein analysis. 
allProt <- as.data.frame(t(protDatList[[2]]))
commonProt <- as.data.frame(t(protDatList[[1]]))
allProt[,which(colnames(allProt) %in% colnames(commonProt))] <- commonProt

#And this is saved

write.csv(allProt, "../External/Stockholm/Data/Combined_average_proteins.csv")

#allProt <- read.csv("../External/Stockholm/Data/Combined_average_proteins.csv", row.names = 1)

#And now, we can combine it all. 
gexCounts <- t(readMM("../External/Stockholm/Data/Anndata_input/gex_counts.mtx"))

rowData <- read.csv("../External/Stockholm/Data/Anndata_input/geneData.csv")
colData <- read.csv("../External/Stockholm/Data/Anndata_input/cellData.csv")
colnames(gexCounts) <- colData$X
rownames(gexCounts) <- rowData$X

fullSce <- SingleCellExperiment(assays = list("counts" = gexCounts), 
                                colData = colData, 
                                rowData = rowData)
altExp(fullSce, "NormADT") <- 
    SingleCellExperiment(assays = list("counts" = t(allProt)))

#We also import the latent representation of the combination of the common transcripts and 
#common proteins. 
latRep <- read.csv("../External/Stockholm/Data/Normalised_protein_models_common_proteins/Latent_representation_common_proteins.csv", row.names = 1)
row.names(latRep) <- colnames(gexCounts)
reducedDim(fullSce, "totalVI_lat_rep") <- latRep
set.seed(1002)
fullSce <- runUMAP(fullSce, dimred = "totalVI_lat_rep")

dDensityPlot(reducedDim(fullSce, "UMAP"), idsVector = fullSce$batch, plotName = "Batch", 
             plotDir = "Diagnostics/Stockholm/Batch_densities")
#This looks very good. 

#No errors or warnings, mean that the cells are in the same order. 
dir.create("../External/Stockholm/Data/SCE_versions")
saveRDS(fullSce, "../External/Stockholm/Data/SCE_versions/1_full_SCE.rds")

