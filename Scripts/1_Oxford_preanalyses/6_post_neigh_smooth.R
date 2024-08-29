library(mclust)
library(DepecheR)
library(BiocParallel)
#Here, we will investigate all the remaining populations. We start by checking
#if the populations should be subclustered and then subcluster the ones that fulfill
#the criteria set up for this. 

popVec <- c("B", "CD4T", "CD8T", "TCRgd", "NK", "ILC")

markerList <- list("B" = c("CD38", "CCR7", "CD95", "CD45", "CCR4", "IgD.TCRgd", 
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

dir.create("Results/Data/Combined_analysis_post_smooth", recursive = TRUE)
dir.create("Results/Graphics/BIC_model_plots", recursive = TRUE)
subclusterability <- as.list(1:length(popVec))
names(subclusterability) <- popVec
hitPops <- c("EOMGhigh", "EOMGlow", "LOMGhigh", "LOMGlow")
set.seed(3412)
for(x in popVec){
  print(x)
  locFile <- readRDS(paste0("../External/Oxford/Resulting_data/", x, "/", x, "_full_file_post_Euclid.rds"))
  locRes <- bplapply(hitPops, function(y){
    if(length(which(locFile$SmoothGroup == y)) > 1){
      mclustData <- as.matrix(dScale(locFile[which(locFile$SmoothGroup == y),markerList[[x]]]))
      if(nrow(mclustData) > 5000){
        mclustData <- mclustData[sample(1:nrow(mclustData), 5000),]
      }
      res <- mclustBIC(mclustData)
      pdf(paste0("Results/Graphics/BIC_model_plots/", x, "_", y, "_model.pdf"))
      plot(res)
      dev.off()
      #Now, our criterion is that we subcluster if the improvement in the statistic
      #is more than 5% with any number of clusters and any method. We also
      #have the additional criterion that if the number of cells are fewer than 100, 
      #we will not subcluster them. 
      
      #In some cases, there are many NAs. For this reason, we will split the data
      #if need be
      
      if(length(which(is.na(as.vector(res[2:nrow(res),])))) > 1){
        allButNoClustRes <- as.vector(res[2:nrow(res),])
        allButNoClustRes <- allButNoClustRes[-which(is.na(allButNoClustRes))]
        maxImprove <- max(res[1,])/max(allButNoClustRes)
      } else {
        maxImprove <- max(res[1,])/max(res[2:nrow(res),])
      }
      
      if(maxImprove > 1.05 && nrow(mclustData) > 1000){
        TRUE
      } else {
        FALSE
      }
    } else {
      FALSE
    }
    
  })
  names(locRes) <- hitPops
  subclusterability[[x]] <- locRes
}

saveRDS(subclusterability, "Results/Data/Combined_analysis_post_smooth/Subclusterability.rds")
#subclusterability <- readRDS("Results/Data/Combined_analysis_post_smooth/Subclusterability.rds")

#Individual penalties have been tested in a previous round, and to optimise the 
#ARIs and keep the number of clusters reasonable, the following values have been
#established
standardPenalties <- list("EOMGhigh" = 2^seq(7, 9, by = 0.5),
                          "EOMGlow" = 2^seq(7, 9, by = 0.5),
                          "LOMGhigh" = 2^seq(7, 9, by = 0.5),
                          "LOMGlow" = 2^seq(7, 9, by = 0.5))
penaltyList <- list("B" = standardPenalties,
                  "CD4T" = list("EOMGhigh" = 2^seq(8, 9, by = 0.5),
                                "EOMGlow" = 2^seq(6, 8, by = 0.5),
                                "LOMGhigh" = 2^seq(7, 9, by = 0.5),
                                "LOMGlow" = 2^seq(7, 9, by = 0.5)),
                  "CD8T" = standardPenalties,
                  "TCRgd" = standardPenalties,
                  "NK" = list("EOMGhigh" = 2^seq(7, 9, by = 0.5),
                              "EOMGlow" = 2^seq(7, 9, by = 0.5),
                              "LOMGhigh" = 2^seq(7, 9, by = 0.5),
                              "LOMGlow" = 2^seq(7, 9, by = 0.5)),
                  "ILC" = standardPenalties)

standardPenalties <- list("EOMGhigh" = 2^seq(5, 9, by = 0.5),
                          "EOMGlow" = 2^seq(5, 9, by = 0.5),
                          "LOMGhigh" = 2^seq(5, 9, by = 0.5),
                          "LOMGlow" = 2^seq(5, 9, by = 0.5))
penaltyList <- list("B" = standardPenalties,
                    "CD4T" = standardPenalties,
                    "CD8T" = list("EOMGhigh" = 2^seq(5, 9, by = 0.5),
                                  "EOMGlow" = 2^seq(5, 9, by = 0.5),
                                  "LOMGhigh" = 2^seq(5, 9, by = 0.5),
                                  "LOMGlow" = 2^seq(4, 6.5, by = 0.5)),
                    "TCRgd" = standardPenalties,
                    "NK" = standardPenalties,
                    "ILC" = standardPenalties)
#############
#FIC CD8 LOMG low, that needs fewer high penalties. The three lowest are enough. 


#Now for the big loop
freqResList <- list()
for(i in popVec){
  print(i)
  datDir <- paste0("../External/Oxford/Resulting_data/", i, "/")
  graphDir <- paste0("Results/Graphics/", i, "/")
  locFile<- readRDS(paste0(datDir, i, "_full_file_post_Euclid.rds"))

  modelDf <- locFile[,markerList[[i]]]

  ########
  #Secondary clustering
  locFile$smoothGroupDeep <- locFile$SmoothGroup
  for(x in hitPops){
    if(subclusterability[[i]][[x]]){
      innerDir <- paste0(datDir, "Smooth_model_hit_subclustering_DepecheR/", x)
      dir.create(innerDir, recursive = TRUE)
      locPenalties <- penaltyList[[i]][[x]]
      locRows <- which(locFile$SmoothGroup == x)
      locDepModel <- depeche(modelDf[locRows,], 
                              penalties = locPenalties,
                              plotDir = innerDir)
      saveRDS(locDepModel, 
              paste0(innerDir, "/DepecheR_model.rds"))
      locFile$smoothGroupDeep[locRows] <- paste0(x, "_", locDepModel$clusterVector)
    } else if(length(which(locFile$SmoothGroup == x)) > 0){
      locRows <- which(locFile$SmoothGroup == x)
      locFile$smoothGroupDeep[locRows] <- paste0(x, "_", locFile$SmoothGroup[locRows])
    }
  }
  saveRDS(locFile, paste0(datDir, i, "_full_file_post_Euclid_including_subclusts.rds"))
  #locFile <- readRDS(paste0(datDir, i, "_full_file_post_Euclid_including_subclusts.rds"))
  
  #Now over to the stats. 
  locLongNameVec <- paste0(locFile$group, "_age_",
                           locFile$age, "_", 
                           locFile$idTime)
  
  countTab <- table(locLongNameVec, locFile$smoothGroupDeep)
  countTab <- countTab[,-which(colnames(countTab) %in% c("None", "Ctrl"))]
  ageVec <- as.numeric(gsub(".+_age_|_.+_.+", "", row.names(countTab)))
  
  lymphVec <- unlist(lapply(row.names(countTab), function(x){
    locFile$lymphCount[which(locLongNameVec == x)[1]]
  }))
  
  freqTab <- as.data.frame(do.call("rbind", lapply(row.names(countTab), function(x){
    locDat <- countTab[which(row.names(countTab) == x),]
    locLymph <- lymphVec[which(row.names(countTab) == x)]
    locDat/locLymph
  })))
  row.names(freqTab) <- row.names(countTab)
  #Here, the just over the threshold LOMG patient is also excluded:
  lomgDat <- freqTab[which(grepl("LOMG_", row.names(freqTab)) &
                             grepl("pre", row.names(freqTab))),]
  eomgDat <- freqTab[which(grepl("EOMG_", row.names(freqTab)) &
                             grepl("pre", row.names(freqTab))),]
  youngCtrlDat <- freqTab[which(grepl("Ctrl_", row.names(freqTab)) &
                                  ageVec < 50 & 
                                  grepl("pre", row.names(freqTab))),]
  oldCtrlDat <- freqTab[which(grepl("Ctrl_", row.names(freqTab)) &
                                ageVec >= 50 & 
                                grepl("pre", row.names(freqTab))),]
  
  #Now, we run statistics on this 
  
  pRes <- do.call("rbind", lapply(colnames(freqTab), function(x){
    if(grepl("high", x)){
      side <- "greater"
    } else {
      side <- "less"
    }
    if(grepl("EOMG", x)){
      eoCtrlP <- wilcox.test(eomgDat[,x], youngCtrlDat[,x], 
                             alternative = side, exact = FALSE)$p.value
      loCtrlP <- 1
      eoLoP <- wilcox.test(eomgDat[,x], lomgDat[,x], 
                           alternative = side,  exact = FALSE)$p.value
    } else {
      eoCtrlP <- 1
      loCtrlP <- wilcox.test(lomgDat[,x], oldCtrlDat[,x], 
                             alternative = side,  exact = FALSE)$p.value
      eoLoP <- wilcox.test(lomgDat[,x], eomgDat[,x], 
                           alternative = side,  exact = FALSE)$p.value
    }
    
    
    c(eoCtrlP, loCtrlP, eoLoP)
  }))
  rownames(pRes) <- colnames(freqTab)
  colnames(pRes) <- c("eomgCtrlP", "lomgCtrlP", "eomgLomgP")
  
  freqTab$Group <- gsub("|_age_.+", "", row.names(freqTab))
  freqTab$Group[which(row.names(freqTab) %in% row.names(youngCtrlDat))] <- "Young_ctrl"
  freqTab$Group[which(row.names(freqTab) %in% row.names(oldCtrlDat))] <- "Old_ctrl"
  
  freqResList[[i]] <- list("pVals" = pRes, "freqs" = freqTab)
}

pValDf <- do.call("rbind", lapply(popVec, function(x){
  locDat <- freqResList[[x]][[1]]
  rownames(locDat) <- paste0(x, "_", rownames(locDat))
  locDat
}))
write.csv(pValDf, "Results/Data/Combined_analysis_post_smooth/pValTable.csv")

#We now export a full frequency table. 
freqTabBig <- do.call("cbind", lapply(popVec, function(x){
  locDat <- freqResList[[x]][[2]]
  colnames(locDat) <- paste0(x, "_", colnames(locDat))
  locDat
}))
#Now, we reduce the group columns to one. 
freqTabBig$Group <- freqTabBig$B_Group

freqTabBig <- freqTabBig[,-grep("_Group", colnames(freqTabBig))]
write.csv(freqTabBig, "Results/Data/Combined_analysis_post_smooth/freqTableAllPops.csv")
