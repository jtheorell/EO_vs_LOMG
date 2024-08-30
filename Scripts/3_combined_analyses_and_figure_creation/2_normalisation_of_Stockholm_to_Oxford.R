#Here, we are going to import the original screen data and see if the there 
#selected significant cells are reproducibly different between the groups in the
#second screen. 

#To this end, we will have to scale the data in some reasonable way. The most 
#seemingly sensible is to use peak normalisation, even if the data is less
#similar in nature than one would like it to be, with the screen data to a large degree being 
#approximations based on transcript characteristics, but we will make do with
#what is at hand, as this will be a very straight forward analysis to interpret
#if it turns out to hold. 

library(flowSpecs)
library(flowCore)
library(DepecheR)
library(SingleCellExperiment)
library(ggplot2)
sthlmProtDat <- readRDS("../External/Stockholm/Data/SCE_versions/Stockholm_protein_data_including_colData.rds")

dir.create("Results/Data/Oxford_and_Stockholm/Harmonisation")
#Now, the row names are changed to match the ones in the second experiment
sthlmProts <- gsub("_PROT", "", row.names(sthlmProtDat))
sthlmProts <- gsub("CD00", "CD", sthlmProts)
sthlmProts <- gsub("CD0", "CD", sthlmProts)
row.names(sthlmProtDat) <- sthlmProts

sthlmProtDf <- as.data.frame(t(logcounts(sthlmProtDat)))

#Now, we will take the data from the ILC_NK_all and BT_all panels and use these
#to normalise the common protein markers in the sthlmProts dataset. After this, 
#we will run the single-cell nearest neighbor analysis to identify common hit cells.
protsToChangeList <- list("BT" = c("CD2", "CD7", "CD8", 
                                   "CD19", "CD20", "CD24",
                                   "CD27", "CD38", "CD45RA", "IgD"),
                          "ILC_NK" = c("CD2", "CD3", "CD7", "CD8", 
                                       "CD16",  "CD19", "CD27", "CD38", "IgD"))
protChangeValList <- list("BT" = c("CD2" = 0.5, "CD7" = 10, 
                                   "CD8" = 10, "CD19" = 0.5, 
                                   "CD20" = 1, "CD24" = 1, 
                                   "CD27" = 3, "CD38" = 30, 
                                   "CD45RA" = 20, "IgD" = 1),
                         "ILC_NK" = c("CD2" = 0.5,"CD3" = 1, "CD7" = 10, 
                            "CD8" = 10, "CD16" = 1, "CD19" = 0.5, 
                            "CD27" = 3, "CD38" = 30, "IgD" = 1))
overarchGroups <- list("BT","ILC_NK")
for(z in overarchGroups){
  loc_file <- readRDS(paste0("Results/Data/", z, "_all/", z, "_all_full_file.rds"))
  #We exclude all controls, as we believe that the patients are better matched
  #to the patient cells. 
  loc_file <- loc_file[-which(loc_file$group == "Ctrl"),]
  #We start by making one change here: the CD3 positive cells are turned negative
  #in the IgD/TCRgd channel, and then this is renamed to IgD. To avoid making
  #the distribution strange by adding only one number here, we will give 
  #all the T cells a non-IgD positive value.
  IgDDat <- loc_file$IgD.TCRgd
  potentiallyTCRgdPosRows <- which(loc_file$CD3 > 375)
  IgDDat[potentiallyTCRgdPosRows] <- IgDDat[sample(which(IgDDat < 300), 
                                                   length(potentiallyTCRgdPosRows))]
  loc_file$IgD <- IgDDat
  locStockMarkers <- Reduce(intersect, list(colnames(loc_file), colnames(sthlmProtDf)))
  
  loc_fileRed <- loc_file[,which(colnames(loc_file) %in% locStockMarkers)]
  sthlmProtDatRedloc <- sthlmProtDf[,which(colnames(sthlmProtDf) %in% locStockMarkers)]
  sthlmProtDatlocOrd <- sthlmProtDatRedloc[,match(colnames(loc_fileRed), 
                                                  colnames(sthlmProtDatRedloc))]
  identical(colnames(loc_fileRed), colnames(sthlmProtDatlocOrd))
  
  #It is clear from investigations of the histograms of individual parameters in
  #the first rounds of analyses below that the negative population of a few of the
  #parameters could be optimised in the Stockholm dataset to better reflect the
  #Oxford shapes. THis is especially clear for CD38 which has been inadvertendly 
  #split. This will be changed here. 

  
  for(i in protsToChangeList[[z]]){
    print(i)
    locCol <- which(colnames(sthlmProtDatlocOrd) == i)
    locRaw <- sthlmProtDatlocOrd[,locCol]
    locUntransed <- sinh(locRaw)*5
    locRetransed <- asinh(locUntransed/protChangeValList[[z]][i])
    sthlmProtDatlocOrd[,locCol] <- locRetransed
  }
  
  
  #We will also, just to make these values plottable, scale the values using the 
  #dScale function in DepecheR.
  loc_fileRedScaled <- dScale(loc_fileRed)
  sthlmProtDatlocScaled <- dScale(sthlmProtDatlocOrd)
  
  stockholm_norm_to_loc <- as.data.frame(exprs(peakNorm(flowSet(flowFrame(as.matrix(sthlmProtDatlocScaled))), 1,
                                                        flowFrame(as.matrix(loc_fileRedScaled)))[[1]]))
  
  #This works well for all markers apart from CD31, that has a more even distributio
  #in the stockholm data, and is therefore not being picked up in this analysis. 
  #We will therefore manually adjust it here. 
  if("CD31" %in% colnames(stockholm_norm_to_loc)){
    CD31Ref <- quantile(loc_fileRedScaled$CD31, c(0.1, 0.9))
    CD31RawQuantiles <- quantile(sthlmProtDatlocScaled$CD31, c(0.1, 0.9))
    firstCD31Step <- (sthlmProtDatlocScaled$CD31-CD31RawQuantiles[1])/
      (CD31RawQuantiles[2]-CD31RawQuantiles[1])
    CD31moved <- (firstCD31Step*(CD31Ref[2]-CD31Ref[1]))+CD31Ref[1]
    sthlmProtDatlocScaled$CD31 <- CD31moved
  }
  
  dir.create(paste0("Results/Graphics/Oxford_and_Stockholm/OxStock_norm/", z), recursive = TRUE)
  sthlmProtDatlocScaled$Set <- "Stockholm_raw"
  stockholm_norm_to_loc$Set <- "Stockholm_norm"
  loc_fileRedScaled$Set <- "Oxford"
  
  largeStocklocDf <- rbind(sthlmProtDatlocScaled, stockholm_norm_to_loc, loc_fileRedScaled)
  
  for(i in seq_len(ncol(largeStocklocDf)-1)){
    locName <- colnames(largeStocklocDf)[i]
    locDf <- largeStocklocDf[,c(i, which(colnames(largeStocklocDf) == "Set"))]
    colnames(locDf)[1] <- "var"
    ggplot(locDf, aes(x=var, fill=Set)) +
      geom_density(alpha = 0.8) + xlab(locName) +ylab("")
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/OxStock_norm/", z, "/", locName, ".pdf"))
  }
  
  dir.create("Results/Data/Oxford_and_Stockholm/Harmonisation")
  largeStockNorm <- largeStocklocDf[which(largeStocklocDf$Set == "Stockholm_norm"),]
  largeStockNorm$cellType <- sthlmProtDat$cellType
  largeStockNorm$hitClust <- sthlmProtDat$smoothSubClust
  #One more change is needed for the integration
  largeStockNorm$cellType[which(largeStockNorm$cellType == "gdT")] <- "TCRgd"
  row.names(largeStockNorm) <- colnames(sthlmProtDat)
  saveRDS(largeStockNorm, 
          paste0("Results/Data/Oxford_and_Stockholm/Harmonisation/Stock_harmonised_to_", z, ".rds"))
  saveRDS(loc_fileRed, 
          paste0("Results/Data/Oxford_and_Stockholm/Harmonisation/Ox_reference_", z, ".rds"))
}

