#Here, we are going to check if the three populations identified in this
#analysis overlap with the ones produced on the other side. 
#Interestingly, all three populations have the same name in both analyses. 

signOxStoc <- read.csv("Results_per_celltype/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv")

#Now, for these three populations, what proportion overlaps between the per-lymphocyte
#and per-cell type analyses?
stockholmOverlap <- lapply(signOxStoc[,1], function(i){
  cellType <- gsub("|_.OMG.+", "", i)
  stockPerCellType <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation_per_cell_type/", cellType, "_Stockholm_data_with_Ox_neighbors.rds"))
  stockPerLymphocytes <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/", cellType, "_Stockholm_data_with_Ox_neighbors.rds"))
  focRowsCellType <- which(stockPerCellType$Ox_hit_clust == i)
  otherSignRowsCellType <- which(stockPerCellType$Ox_hit_clust %in% c(i, "None") == FALSE)
  resVecCellType <- rep("allOtherCells", nrow(stockPerCellType))
  resVecCellType[focRowsCellType] <- "cellTypeDefinedCluster"
  resVecCellType[otherSignRowsCellType] <- "otherCellTypeDefinedClusters"
  
  focRowsLymphocyte <- which(stockPerLymphocytes$Ox_hit_clust == i)
  otherSignRowsLymphocyte <- which(stockPerLymphocytes$Ox_hit_clust %in% c(i, "None") == FALSE)
  resVecLymphocyte <- rep("allOtherCells", nrow(stockPerLymphocytes))
  resVecLymphocyte[focRowsLymphocyte] <- "lymphocyteDefinedCluster"
  resVecLymphocyte[otherSignRowsLymphocyte] <- "otherLymphocyteDefinedClusters"
  list("StockholmTable" = table(resVecCellType, resVecLymphocyte), 
       "StockholmFisher" = fisher.test(resVecCellType, resVecLymphocyte, 
                                       simulate.p.value = TRUE))
})

#And the same for the Oxfor data
oxfordOverlap <- lapply(signOxStoc[,1], function(i){
  cellType <- gsub("|_.OMG.+", "", i)
  subset <- gsub("^.{2,4}_|", "", i)
  oxPerCellType <- readRDS(paste0("../External/Oxford/Resulting_data_analysis_per_celltype/", cellType, "/", cellType, 
                                  "_full_file_post_Euclid_including_subclusts.rds"))
  
  oxPerLymphocytes <- readRDS(paste0("../External/Oxford/Resulting_data/", cellType, "/", cellType, 
                                     "_full_file_post_Euclid_including_subclusts.rds"))

  focRowsCellType <- which(oxPerCellType$smoothGroupDeep == subset)
  otherSignRowsCellType <- which(oxPerCellType$smoothGroupDeep %in% c(subset, "None") == FALSE)
  resVecCellType <- rep("allOtherCells", nrow(oxPerCellType))
  resVecCellType[focRowsCellType] <- "cellTypeDefinedCluster"
  resVecCellType[otherSignRowsCellType] <- "otherCellTypeDefinedClusters"
  
  focRowsLymphocyte <- which(oxPerLymphocytes$smoothGroupDeep == subset)
  otherSignRowsLymphocyte <- which(oxPerLymphocytes$smoothGroupDeep %in% c(subset, "None") == FALSE)
  resVecLymphocyte <- rep("allOtherCells", nrow(oxPerLymphocytes))
  resVecLymphocyte[focRowsLymphocyte] <- "lymphocyteDefinedCluster"
  resVecLymphocyte[otherSignRowsLymphocyte] <- "otherLymphocyteDefinedClusters"
  list("StockholmTable" = table(resVecCellType, resVecLymphocyte), 
       "StockholmFisher" = fisher.test(resVecCellType, resVecLymphocyte, 
                                       simulate.p.value = TRUE))
})

#And this is saved.
saveRDS(oxfordOverlap, "Results_per_celltype/Ox_lymph_vs_cell_type_defined_hits.rds")
saveRDS(stockholmOverlap, "Results_per_celltype/Stock_lymph_vs_cell_type_defined_hits.rds")

#The results indicate that the CD4 population is strongly overlapping, that 
#the CD8 population does not overlap at all and that the NK population is mixed. 
#With this, we put this way of analysing downl, at least for now. 
