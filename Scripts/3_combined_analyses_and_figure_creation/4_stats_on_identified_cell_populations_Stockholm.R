library(SingleCellExperiment)
#Now stats. 
cellPops <- list("B", "CD4T", "CD8T", "TCRgd", "ILC", "NK")
sthlmProtDat <- readRDS("../External/Stockholm/Data/SCE_versions/Stockholm_protein_data_including_colData.rds")
stockColDat <- colData(sthlmProtDat)

freqResList <- list()
for(i in cellPops){
  print(i)
  locDatPat <- readRDS(paste0("Results/Data/Oxford_and_Stockholm/Harmonisation/", i, "_Stockholm_data_with_Ox_neighbors.rds"))
  if(i == "ILC"){
    locMetaDatPat <- stockColDat[which(stockColDat$cellType == "NK"),]
  } else if(i == "TCRgd"){
    locMetaDatPat <- stockColDat[which(stockColDat$cellType == "gdT"),]
  } else {
    locMetaDatPat <- stockColDat[which(stockColDat$cellType == i),]
  }
  
  if(identical(row.names(locMetaDatPat), row.names(locDatPat))){
    locLongNameVec <- paste0(locMetaDatPat$MG_type, "_age_", 
                             locMetaDatPat$Age, "_", 
                             locMetaDatPat$sample)

    countTab <- table(locLongNameVec, locDatPat$Ox_hit_clust)
    countTab <- countTab[,-which(colnames(countTab) == "None")]
    #Here, we introduce an arbitrary criterion, namely that only clusters with #
    #more than 50 cells are considered. 
    countTabFilt <- apply(as.matrix(countTab), 2, function(x){
      if(sum(x) < 50){
        rep(0, length(x))
      } else {
        x
      }
    })
    row.names(countTabFilt) <- row.names(countTab)
    #Here, we do not go on if there is not a single cluster that has passed. 
    if(sum(countTabFilt) > 50){
      ageVec <- as.numeric(gsub(".+_age_|_.+_.+", "", row.names(countTabFilt)))
      lymphVec <- unlist(lapply(row.names(countTabFilt), function(x){
        length(which(sthlmProtDat$sample == gsub(".+_age_.._|", "", x)))
      }))
      freqTab <- as.data.frame(do.call("rbind", lapply(row.names(countTabFilt), function(x){
        locDat <- countTabFilt[which(row.names(countTabFilt) == x),]
        locLymph <- lymphVec[which(row.names(countTabFilt) == x)]
        locDat/locLymph
      })))
      row.names(freqTab) <- row.names(countTabFilt)
      
      lomgDat <- freqTab[which(grepl("LOMG_", row.names(freqTab))),]
      eomgDat <- freqTab[which(grepl("EOMG_", row.names(freqTab))),]
      #Now, we run statistics on this. 
      #We use one-tailed tests, as we do know the direction we expect each time. 
      
      pRes <- do.call("rbind", lapply(colnames(freqTab), function(x){
        if(grepl("EOMGhigh|LOMGlow", x)){
          side <- "greater"
        } else {
          side <- "less"
        }
        eoLoP <- wilcox.test(eomgDat[,x], lomgDat[,x], 
                             alternative = side, 
                             exact = FALSE)$p.value
      }))
      rownames(pRes) <- colnames(freqTab)
      colnames(pRes) <- "eomgLomgP"
      freqTab$Group <- gsub("|_age_.+", "", row.names(freqTab))
      freqResList[[i]] <- list("pVals" = pRes, "freqs" = freqTab)
    } else {
      list("pVals" = NULL, "freqs" = NULL)
    }
    
  } else {
    stop("Disarray reigns among the row names. Check the reason.")
  }
}

#There were no ILCs. 
pValDf <- do.call("rbind", lapply(cellPops[-which(cellPops == "ILC")], function(x){
  locDat <- freqResList[[x]][[1]]
}))
write.csv(pValDf, "Results/Data/Oxford_and_Stockholm/Harmonisation/pValTable_Stockholm.csv")

#We now exclude the ones that are not significant. 
pValDfSig <- pValDf[which(pValDf[,1] < 0.05),]
round(pValDfSig, 4)

#   B_LOMGlow_2 CD4T_EOMGlow_3 CD8T_LOMGlow_2   NK_EOMGlow_1   NK_EOMGlow_2 
#        0.0306         0.0227         0.0119         0.0085         0.0037 
#  NK_EOMGlow_3 
#        0.0010 

write.csv(pValDfSig, "Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stockholm.csv")

#Nice. We also save the frequency table. 
freqTabBig <- do.call("cbind", lapply(cellPops[-which(cellPops == "ILC")], function(x){
  locDat <- freqResList[[x]][[2]]
}))
write.csv(freqTabBig, "Results/Data/Oxford_and_Stockholm/Harmonisation/freqTable_Stockholm.csv")




