#And now, we can go into the large loop, where we use these normalised values
#to identify the nearest neighbors in each cell population. We will not in this
#case restrict the number of neighbors for each donor, as it is beyond the point: 
#if one cell is the closest for the whole dataset - so be it. 
library(DepecheR)
library(harmony)
library(gmodels)
library(uwot)
library(FNN)
library(viridis)
#Here, are the proteins of focus for each cell population. In this case, 
#We will recombine the ILC and NK cells, as they are not separated in the original
#file

markerList <- list("B" = c("CD38", "CCR7", "CD95", "CD45", "CCR4", "IgD", 
                           "CD11b", "CD27", "CD20", "CD45RA", "CD31", 
                           "CD24", "IgM"),
                   "CD4T" = c("CD38", "KLRG1", "CCR7", "CD95", "CD25", "CD4", "CCR4", "CD71", 
                              "CD11b", "CD27", "CD64", "CD20", "CD45RA", "CD2", "CD31", "CD161", "CD7",
                              "CD24", "CD152"),
                   "CD8T" = c("CD38", "KLRG1", "CCR7", "CD95", "CD25", "CCR4", "CD71", "CD8",
                              "CD11b", "CD27", "CD64", "CD20", "CD45RA", "CD2", "CD31", "CD161", "CD7",
                              "CD24", "CD152"),
                   "TCRgd" = c("CD38", "KLRG1", "CCR7", "CD95", "CD25", "CCR4", "CD71", "CD8",
                               "CD11b", "CD27", "CD64", "CD20", "CD45RA", "CD2", "CD31", "CD161", "CD7",
                               "CD24", "CD152"),
                   "NK" = c("NKp30","HLADR","CCR7", 
                            "CD183","CD161","CD16",
                            "CD34","CD56", "CD8",
                            "CD127", "CD27", "CD57",
                            "NKG2A", "CD2", "NKG2C",
                            "CRTH2", "CD7","CD117","CD218a"),
                   "ILC" = c("NKp30", "HLADR",
                             "CCR7", "CD183",
                             "CD161", "CD56", 
                             "CD127", "CD27", 
                             "NKG2A", "CD2", 
                             "CRTH2", "CD7",
                             "CD117", "CD218a"))

#How long is actually the list of unique markers? 
length(unique(unlist(lapply(markerList, function(x) x))))
#36. As CD3, CD11c and CD19 are not above, that means a total of 32+7 measured markers per dataset. 

listOfPopVecs <- list("BT" = c("B", "CD4T", "CD8T", "TCRgd"),
                         "ILC_NK" = c("ILC", "NK"))
dir.create("Results/Graphics/Oxford_and_Stockholm/Harmonisation")
for(z in names(listOfPopVecs)){
  print(z)
  stockRef <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/Stock_harmonised_to_", z, ".rds"))
  
  #This is included for the dScaling
  oxRef <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/Ox_reference_", z, ".rds"))
  for(j in listOfPopVecs[[z]]){
    print(j)
    allDat <- readRDS(paste0("../External/Oxford/Resulting_data/", j, "/", j, "_full_file_post_Euclid_including_subclusts.rds"))
    #We start by zooming in on the cells of relevance, i.e. the pre-thymectomy
    #PBMC samples from patients. 
    locDat <- allDat[which(allDat$tissue == "pre" &
                             allDat$group != "Ctrl"),]
    
    if(j == "B"){
      #We start by making one change here: the CD3 positive cells are turned negative
      #in the IgD/TCRgd channel, and then this is renamed to IgD. To avoid making
      #the distribution strange by adding only one number here, we will give 
      #all the T cells a non-IgD positive value.
      IgDDat <- locDat$IgD.TCRgd
      potentiallyTCRgdPosRows <- which(locDat$CD3 > 375)
      IgDDat[potentiallyTCRgdPosRows] <- IgDDat[sample(which(IgDDat < 300), 
                                                       length(potentiallyTCRgdPosRows))]
      locDat$IgD <- IgDDat
    } 
    
    #Now, we identify the markers among the local population defining ones that
    #are present in both datasets. 
    markerOverlap <- Reduce(intersect, list(markerList[[j]], colnames(stockRef)))
    #We now reduce the three datasets to this and dScale our dataset to the
    #Oxford reference, to put it in the right place. 
    locDatRed <- locDat[,which(colnames(locDat) %in% markerOverlap)]
    stockRefRed <- stockRef[,which(colnames(stockRef) %in% markerOverlap)]
    oxRefRed <- oxRef[,which(colnames(oxRef) %in% markerOverlap)]
    locDatOrd <- locDatRed[,match(colnames(oxRefRed), colnames(locDatRed))]
    stockRefOrd <- stockRefRed[,match(colnames(oxRefRed), colnames(stockRefRed))]
    all(unlist(lapply(list(locDatOrd, stockRefOrd), function(x){
      identical(colnames(oxRefRed), colnames(x))
    }))) #It is TRUE in all cases. 
    #Now, we scale the local data
    locDatScaled <- dScale(locDatOrd, control = oxRefRed)
    #Here, the stockRef is now reduced to the relevant cell type. As the ILC
    #compartment was too small to be detected in the Stockholm cohort, we will
    #here have to stick to the NK compartment, which should in theory also include
    #these cells. 
    
    if(j == "ILC"){
      stockRefRows <- which(stockRef$cellType == "NK")
    } else {
      stockRefRows <- which(stockRef$cellType == j)
    }
    
    stockRefLoc <- stockRefOrd[stockRefRows,]
    oxHitVec <- locDat$smoothGroupDeep
    oxHitVec[which(oxHitVec != "None")] <- paste0(j, "_", oxHitVec[which(oxHitVec != "None")])
    #Now, we combine the datasets, create a PCA and then harmonise this PCA, 
    #followed by an umap creation on a reduced dataset. 
    locDatScaled$Set <- "Ox"
    stockRefLoc$Set <- "Stock"
    locFullDat <- rbind(locDatScaled, stockRefLoc)
    
    #set.seed(111)
    #oxStockPCA <- fast.prcomp(locFullDat[,1:(ncol(locFullDat)-1)])$x
    #if(ncol(oxStockPCA) > 10){
    #  oxStockPCA <- oxStockPCA[,1:10]
    #}
    #saveRDS(oxStockPCA, paste0("../External/Oxford_and_Stockholm/Harmonisation/", j, "_raw_PCA.rds"))
    oxStockPCA <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/", j, "_raw_PCA.rds"))
    
    #set.seed(1010)
    #harmonyPCA <- 
    #  RunHarmony(oxStockPCA, meta_data = locFullDat$Set)
    #saveRDS(harmonyPCA, paste0("../External/Oxford_and_Stockholm/Harmonisation/", j, "_harmony_PCA.rds"))
    harmonyPCA <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/", j, "_harmony_PCA.rds"))
    
    #And now comes the key thing - we identify the nearest neighbors for all
    #the Oxford datapoints in the Stockholm dataset.
    oxRows <- which(locFullDat$Set == "Ox")
    stockRows <- which(locFullDat$Set == "Stock")
    
    timing <- Sys.time()
    oxNeigh <- knnx.index(harmonyPCA[oxRows,], 
                             harmonyPCA[stockRows,], k = 1)
    Sys.time()-timing
    #This is now converted to the cell type distribution. 
    oxType <- oxHitVec[oxNeigh[,1]]
    
    #And this is now saved. 
    neighborised_Stock <- cbind(stockRef[stockRefRows,], "Ox_neigh" = oxNeigh[,1],
                        "Ox_hit_clust" = oxType)
    saveRDS(neighborised_Stock, paste0("../External/Oxford_and_Stockholm/Harmonisation/", j, "_Stockholm_data_with_Ox_neighbors.rds"))
  }
}


