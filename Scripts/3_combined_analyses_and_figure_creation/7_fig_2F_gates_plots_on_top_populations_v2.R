#Here, we will run a two-step optimisation. First, we will see which levels of gates
#for all the markers of interest for each cell type that can be used to create an 
#as oure replica of our populations as possible. Then, we will take the markers
#and see which the smallest set of markers is that can be used to do this. 
#As this analysis is heavy, the marker sets are slightly reduced compared to
#previous sets. 

library(ggplot2)
library(dfoptim)
library(reshape2)
library(tidysdm)
library(flowCore)

#The 
pValFoc <- 
  suppressWarnings(read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1))
#This creates an irrelevant warning, hence the suppression. 

#Previously, we had populations from B and CD4 in addition to CD8 and NK, but no longer. 
#popVec <- c("B", "CD4T", "CD8T", "NK")
popVec <- c("CD8T", "NK")

markerList <- list("B" = c("CD38", "CCR7", "CCR4", "IgD.TCRgd", 
                           "CD11b", "CD27", "CD20", "CD45RA", "CD31", 
                           "CD24", "IgM"),
                   "CD4T" = c("CD38", "CCR7", "CD25", "CCR4", 
                              "CD11b", "CD27", "CD20", "CD45RA", 
                              "CD2", "CD31", "CD161", "CD7",
                              "CD24", "CD152"),
                   "CD8T" = c("CD38", "KLRG1", "CCR7", "CCR4", 
                              "CD11b", "CD27", "CD20", "CD45RA", "CD2", 
                              "CD31", "CD161", "CD7",
                              "CD24", "CD152"),
                   "NK" = c("NKp30","HLADR","CCR7", 
                            "CD183","CD161","CD16",
                            "CD56", "CD57",
                            "NKG2A", "CD2", "NKG2C"))

#After a number of iterations, it has been realised that the best thing is to 
#focus on the markers that are needed for separation as well as a few phenotypic
#markers. 
focMarkerlist <- list("B" = c("CD38", "IgD.TCRgd", "CD27", 
                              "CD24", "IgM", "CD20"),
                      "CD4T" = c("CCR4", "CD45RA", "CCR7", "CD27",
                                 "CD31", "CD25"),
                      "CD8T" = c("KLRG1", "CCR7", "CD161",
                                 "CD45RA", "CD7", "CD27"),
                      "NK" = c("NKp30","CD16", "NKG2C", "CD2",
                               "CD56", "CD57"))

#Now, we will do something reasonably simple. We will first identify the 99th 
#and the 1st percentile for each marker for the population of interest. We
#will then rank all these markers based on their "uniqueness" compared to the
#rest of the population, i.e. where the ratio of our population to the 
#rest of the population is the best. Then we take these, one-by-one, and run a 
#Wilcoxon rank sum test, and we do this in a loop that breaks when we reach 
#significance.

bestPvals <- list()
selectedMarkers <- list()
freqList <- list()
gateValList <- list()
for(i in popVec){
  datDir <- paste0("../External/Oxford/Resulting_data/", i, "/")
  locFile <- readRDS(paste0(datDir, i, "_full_file_post_Euclid_including_subclusts.rds"))
  locFilePre <- locFile[which(locFile$tissue == "pre"),]
  for(j in 1:length(grep(i, row.names(pValFoc)))){
    locFocClust <- row.names(pValFoc)[grep(i, row.names(pValFoc))][j]
    print(i)
    print(locFocClust)
    focRows <- which(locFilePre$smoothGroupDeep %in% gsub("^.{1,4}_|", "", locFocClust))
    focPop <- locFilePre[focRows,markerList[[i]]]
    focPopMeta <- locFilePre[focRows,c("idTime", "group", "age", "lymphCount")]
    altPop <- locFilePre[-focRows,markerList[[i]]]
    fullPop <- locFilePre[,markerList[[i]]]
    #Now, to construct the gates, we will use the madFilter function. For this, 
    #we do however need to clarify which side of the peak that the gate should be 
    #present, which will be done through median value calvulations for both
    #populations. 
    focMedianVec <- apply(focPop, 2, median)
    altMedianvec <- apply(altPop, 2, median)
    
    gateSideVec <- sapply(seq_along(focMedianVec), function(x){
      if(focMedianVec[x] > altMedianvec[x]){
        "high"
      } else {
        "low"
      }
    })
    names(gateSideVec) <- names(focMedianVec)
    
    #Now, we check how many peaks we get for each, and all the stats we might need
    #to create the most logical gate there is for each one. 
    gateDataFull <- lapply(colnames(fullPop), function(x){
      locPeakRes <- flowSpecs:::peakIdenti(fullPop[,x], returnStats = TRUE)
      if(length(locPeakRes$PeakPos) == 2){
        c("gateVal" = locPeakRes$Width[[1]][2], 
          "nPeakDist" = locPeakRes$PeakPos[2]-locPeakRes$PeakPos[1]) 
      } else {
        #If there is only one peak, we define the middle between the peak of the 
        #focus and non-focus populations. 
        focPeakRes <- flowSpecs:::peakIdenti(focPop[,x], returnStats = TRUE)
        altPeakRes <- flowSpecs:::peakIdenti(altPop[,x], returnStats = TRUE)
        peakDistance <- abs(focPeakRes$PeakPos-altPeakRes$PeakPos)
        
        c("gateVal" = min(c(focPeakRes$PeakPos, altPeakRes$PeakPos))+(peakDistance/2), 
          "nPeakDist" = peakDistance) 
      }
    })
    names(gateDataFull) <- colnames(fullPop)
    gateData <- unlist(lapply(gateDataFull, "[[", 1))
    peakDist <- unlist(lapply(gateDataFull, "[[", 2))
    if(i == "CD4T"){
      gateData["CCR7"] <-  375
    } else if(i == "CD8T"){
      gateData["CD45RA"] <- 500
      
    } else if(i == "NK"){
      gateData["CD57"] <- 625
      gateData["NKp30"] <- 270
      gateData["CD16"] <- 575
      gateData["NKG2C"] <- 250
      gateData["CD56"] <- 375
    }
    
    
    #Now, we calculate the fraction of the non-focus cells within this region of 
    #interest compared to the fraction of the focus cells and then rank the markers 
    #based on this. 
    altPopFracs <- sapply(colnames(altPop), function(x){
      if(gateSideVec[x] == "high"){
        length(which(altPop[,x] > gateData[x]))/nrow(altPop)
      } else {
        length(which(altPop[,x] < gateData[x]))/nrow(altPop)
      }
    })
    popFracs <- sapply(colnames(focPop), function(x){
      if(gateSideVec[x] == "high"){
        length(which(focPop[,x] > gateData[x]))/nrow(focPop)
      } else {
        length(which(focPop[,x] < gateData[x]))/nrow(focPop)
      }
    })
    
    fracLeverage <- popFracs/altPopFracs
    #We also weigth these, so that those with more separated peaks are priorited. 
    #We do however scale this effect down somewhat by taking the square root of the
    #distance between the peaks. 
    fracLeverageWeighted <- fracLeverage*sqrt(peakDist)
    
    #So now, the markers are ranked. 
    fracLeveragesOrd <- fracLeverageWeighted[order(fracLeverageWeighted, decreasing = TRUE)]
    
    #And here, we restrict the markers to the ones we focus on. 
    fracLeveragesOrd <- fracLeveragesOrd[which(names(fracLeveragesOrd) %in% focMarkerlist[[i]])]
    #And now for the while loop
    pPat <- 1
    pCtrl <- 1
    #We restrict this, so that it never uses less than two markers
    numOfMarkers <- 1
    while(numOfMarkers < length(fracLeveragesOrd) && (pPat > 0.05 || pCtrl > 0.05)){
      numOfMarkers <- numOfMarkers+1
      usedMarkers <- names(fracLeveragesOrd)[1:numOfMarkers]
      gatedRows <- Reduce(intersect, lapply(usedMarkers, function(x){
        if(gateSideVec[x] == "high"){
          which(locFilePre[,x] > gateData[x])
        } else {
          which(locFilePre[,x] < gateData[x])
        }
      }))
      
      locFocMeta <- locFilePre[gatedRows,]
      eomgCounts <- as.matrix(table(locFocMeta$idTime[which(locFocMeta$group == "EOMG")]))
      eomgFreqs <- sapply(row.names(eomgCounts), function(x){
        eomgCounts[which(row.names(eomgCounts) == x),1]/
          locFocMeta$lymphCount[which(locFocMeta$idTime == x)][1]
      })
      lomgCounts <- as.matrix(table(locFocMeta$idTime[which(locFocMeta$group == "LOMG")]))
      lomgFreqs <- sapply(row.names(lomgCounts), function(x){
        lomgCounts[which(row.names(lomgCounts) == x),1]/
          locFocMeta$lymphCount[which(locFocMeta$idTime == x)][1]
      })
      #Now, we know that we only have populations that are lower for both LOMG
      #and EOMG. For that reason, we can focus entirely on these one-tailed analyses.
      if(grepl("EOMG", locFocClust)){
        ctrlCounts <- 
          as.matrix(table(locFocMeta$idTime[which(locFocMeta$group == "Ctrl" &
                                                    locFocMeta$age < 50)]))
      } else {
        ctrlCounts <- 
          as.matrix(table(locFocMeta$idTime[which(locFocMeta$group == "Ctrl" &
                                                    locFocMeta$age > 50)]))
      }
      ctrlFreqs <- sapply(row.names(ctrlCounts), function(x){
        ctrlCounts[which(row.names(ctrlCounts) == x),1]/
          locFocMeta$lymphCount[which(locFocMeta$idTime == x)][1]
      })
      #And now for the stats
      if(grepl("EOMG", locFocClust)){
        pPat <- wilcox.test(eomgFreqs, lomgFreqs, alternative = "less")$p.value
        pCtrl <- wilcox.test(eomgFreqs, ctrlFreqs, alternative = "less")$p.value
      } else {
        pPat <- wilcox.test(lomgFreqs, eomgFreqs, alternative = "less")$p.value
        pCtrl <- wilcox.test(lomgFreqs, ctrlFreqs, alternative = "less")$p.value
      }
    }
    freqList[[locFocClust]] <- list("EOMG" = eomgFreqs,
                          "LOMG" = lomgFreqs,
                          "Ctrl" = ctrlFreqs)
    bestPvals[[locFocClust]] <- c("patPat" = pPat, "patCtrl" = pCtrl)
    if(all(bestPvals[[locFocClust]] < 0.05)){
      selectedMarkers[[locFocClust]] <- usedMarkers
    } else {
      selectedMarkers[[locFocClust]] <- "None"
    }
    
    #Now, we are going to select the markers of interest for each population. 
    #We will do this by taking the most different markers, based on our gating
    #and ordering analysis above. We will make a supplementary showing all markers
    #but focus on the top ones in the figure. 
    
    focMarks <- names(fracLeveragesOrd)
    set.seed(198765432)
    focPopFoc <- focPop[sample(1:nrow(focPop), 5000),
                        which(colnames(focPop) %in% focMarks)]
    altPopFoc <- altPop[sample(1:nrow(altPop), 5000),
                        which(colnames(altPop) %in% focMarks)]
    focPopFoc$Population <- "focus"
    altPopFoc$Population <- "other"
    locDatFull <- rbind(focPopFoc, altPopFoc)
    
    locDatLong <- melt(locDatFull, id.vars = "Population")
    locDatLong$variable <- factor(locDatLong$variable, 
                                  levels = focMarks)
    
    locDatLong$minVals <- locDatLong$maxVals <- NA
    if(selectedMarkers[[locFocClust]][1] != "None"){
      for(j in selectedMarkers[[locFocClust]]){
        if(gateSideVec[j] == "high"){
          locDatLong$minVals[which(locDatLong$variable == j)] <- 
            gateData[j]
          locDatLong$maxVals[which(locDatLong$variable == j)] <- 
            max(locDatLong$value)
        } else {
          locDatLong$minVals[which(locDatLong$variable == j)] <- 
            0
          locDatLong$maxVals[which(locDatLong$variable == j)] <- 
            gateData[j]
        }
        
      }
    }
    
    locDatLong$numVariable <- NA
    #We know that the number of rows is the same for each. 
    numberOfRows <- nrow(locDatLong)/length(unique(locDatLong$variable))
    for(j in 1:length(unique(locDatLong$variable))){
      whichRows <- which(as.numeric(locDatLong$variable) == j)
      locDatLong$numVariable[whichRows] <- seq(j-0.5, j+0.5, 
                                               length.out = numberOfRows)
    }
    #locDatLong$numVariable <- locDatLong$numVariable[order(locDatLong$variable)]
    
    if(grepl("EOMG", locFocClust)){
      colorForPlot <- "red4"
    } else {
      colorForPlot <- "blue4"
    }
    
    p <- ggplot(locDatLong, aes(
      x = variable,
      y = value,
      fill = Population)) +
      geom_split_violin(scale = "width") +
      theme_bw() + scale_fill_manual(values = c(colorForPlot, "#DDDDDD")) +
      scale_y_continuous(limits = c(0, 1050), expand = c(0,0)) 
    if(selectedMarkers[[locFocClust]][1] != "None"){
      p <- p + geom_ribbon(locDatLong, mapping = aes(x = numVariable,
                                                     ymin = minVals, 
                                                     ymax = maxVals),
                           alpha = 0.3, inherit.aes = FALSE)
    }
    
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                  locFocClust, ".pdf"), plot = p)
    p <- p+theme(legend.position = "None", 
                 axis.title = element_blank())
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                  locFocClust, "_stripped.pdf"), plot = p, 
           width = 4, height = 2)
    #As this does not work as intended for the "gates", these being corrupted
    #repeatedly when importing into illustrator, I will also create a png output 
    #for this purpose. 
    p <- p+theme_void() + theme(legend.position = "none")
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                  locFocClust, "_stripped_bare.png"), plot = p, 
           width = 5, height = 2)
    
    #And now at the end, we also plot the full thing. 
    set.seed(18)
    focPopFoc <- focPop[sample(1:nrow(focPop), 5000),]
    altPopFoc <- altPop[sample(1:nrow(altPop), 5000),]
    focPopFoc$Population <- "focus"
    altPopFoc$Population <- "other"
    locDatFull <- rbind(focPopFoc, altPopFoc)
    
    allNames <- colnames(focPop)
    nameOrd <- allNames[order(allNames)]
    CDMols <- grep("CD", nameOrd)
    nameOrd[CDMols] <- 
      nameOrd[CDMols][order(as.numeric(gsub("CD|RA|a|b", "", nameOrd[CDMols])))]
    
    locDatLong <- melt(locDatFull, id.vars = "Population")
    locDatLong$variable <- factor(locDatLong$variable, 
                                  levels = nameOrd)
    p <- ggplot(locDatLong, aes(
      x = variable,
      y = value,
      fill = Population)) +
      geom_split_violin(scale = "width") +
      theme_bw() + scale_fill_manual(values = c(colorForPlot, "#DDDDDD")) +
      scale_y_continuous(limits = c(0, 1050), expand = c(0,0)) 
    
    
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                  locFocClust, "_all_markers.pdf"), plot = p)
    p <- p+theme(legend.position = "None", 
                 axis.title = element_blank())
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                  locFocClust, "_all_markers_stripped.pdf"), plot = p, 
           width = 7, height = 2)
    p <- p+theme_void() + theme(legend.position = "none")
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                  locFocClust, "_all_markers_stripped_bare.png"), plot = p, 
           width = 5, height = 2)
    gateValList[[locFocClust]] <- list("Gates" = gateData, 
                             "Interesting_side_of_gates" = gateSideVec)
  }

}

bestPvals
#$CD8T_LOMGlow_1
#     patPat     patCtrl 
#0.001560824 0.030223270 
#
#$CD8T_LOMGlow_2
#     patPat     patCtrl 
#0.036770002 0.00138297
#
#$NK_EOMGlow_1
#     patPat     patCtrl 
#0.021034815 0.005576467

selectedMarkers
#$CD8T_LOMGlow_1
#[1] "CCR7"   "CD45RA"
#
#$CD8T_LOMGlow_2
#[1] "CD45RA" "CD7"   
#
#$NK_EOMGlow_1
#[1] "CD57"  "CD2"   "NKG2C" "NKp30" "CD16" 

#Now, we save the gate values and the selected markers
dir.create("Results/Data/Oxford_and_Stockholm/Gates_finding_clusters")
saveRDS(selectedMarkers, "Results/Data/Oxford_and_Stockholm/Gates_finding_clusters/selectedMarkers.rds")
saveRDS(freqList, "Results/Data/Oxford_and_Stockholm/Gates_finding_clusters/freqList.rds")
saveRDS(gateValList, "Results/Data/Oxford_and_Stockholm/Gates_finding_clusters/gateValList.rds")

#selectedMarkers <- readRDS("Results/Data/Oxford_and_Stockholm/Gates_finding_clusters/selectedMarkers.rds")
#gateValList <- readRDS("Results/Data/Oxford_and_Stockholm/Gates_finding_clusters/gateValList.rds")
#Now, to convince Dr Ruffin, I will also create conventional gates, for show, that 
#will go in hte supplement. It is either 2 or 4 parameters that are needed
#to identify the population of interest. 
dir.create("Results/Graphics/Oxford_and_Stockholm/Gates_finding_significant/")
for(i in popVec){
  datDir <- paste0("../External/Oxford/Resulting_data/", i, "/")
  locFile <- readRDS(paste0(datDir, i, "_full_file_post_Euclid_including_subclusts.rds"))
  locFilePre <- locFile[which(locFile$tissue == "pre"),]
  locFocClusts <- row.names(pValFoc)[grep(i, row.names(pValFoc))]
  for(locFocClust in locFocClusts){
    print(i)
    print(locFocClust)
    focRows <- which(locFilePre$smoothGroupDeep %in% gsub("^.{1,4}_|", "", locFocClust))
    focPop <- locFilePre[focRows,selectedMarkers[[locFocClust]]]
    altPop <- locFilePre[-focRows,selectedMarkers[[locFocClust]]]
    set.seed(10101101)
    if(nrow(focPop) > 5000){
      focPop <- focPop[sample(1:nrow(focPop), 5000),]
    }
    if(nrow(altPop) > 10000){
      altPop <- altPop[sample(1:nrow(altPop), 10000),]
    }
    locSelMark <- selectedMarkers[[locFocClust]]
    #Now, in the case where the number of markers defining a population is uneven, 
    #we will reuse the most separating marker also for this last one. 
    if(length(locSelMark) %% 2 != 0){
      locSelMark <- c(locSelMark, locSelMark[1])
    }
    for(j in 1:(length(locSelMark)/2)){
      xMark <- locSelMark[((j*2)-1)]
      yMark <- locSelMark[((j*2))]
      xMarkFoc <- focPop[,xMark]
      xMarkAlt <- altPop[,xMark]
      yMarkFoc <- focPop[,yMark]
      yMarkAlt <- altPop[,yMark]
      ggDat <- data.frame("xMark" = c(xMarkAlt, xMarkFoc),
                          "yMark" = c(yMarkAlt, yMarkFoc),
                          "Group" = c(rep("Ref", nrow(altPop)),
                                      rep("Focus", nrow(focPop)))
      )
      #Now, the gate is constructed
      xGate <- gateValList[[locFocClust]]$Gates[xMark]
      xSide <- gateValList[[locFocClust]]$Interesting_side_of_gates[xMark]
      yGate <- gateValList[[locFocClust]]$Gates[yMark]
      ySide <- gateValList[[locFocClust]]$Interesting_side_of_gates[yMark]
      rect_coords <- data.frame(
        xmin = ifelse(xSide == "high", xGate, min(ggDat$xMark)),   
        xmax = ifelse(xSide == "high", max(ggDat$xMark), xGate),    
        ymin = ifelse(ySide == "high", yGate, min(ggDat$yMark)),   
        ymax = ifelse(ySide == "high", max(ggDat$yMark), yGate)     
      )
      p <- ggplot(ggDat, aes(x = xMark, y = yMark, color = Group)) +
        geom_point() + theme_bw() +
        scale_color_manual(values = c(ifelse(i %in% c("CD4T", "NK"), "#FF000011", "#0000FF44"), "#33333377")) +
        geom_rect(data = rect_coords, inherit.aes = FALSE,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  color = "black", fill = NA, alpha = 1) + 
        theme(legend.position = "none") + ylab(yMark) + xlab(xMark) + 
        scale_x_continuous(limits = c(0, 1050), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1050), expand = c(0, 0))
      p
      ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Gates_finding_significant/", locFocClust, "_", xMark, "_vs_", yMark, "gates.pdf"), width = 5, height = 5)
      #And we also save a super-stripped, png version
      p <- p+theme_void() + theme(legend.position = "none")
      p
      ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Gates_finding_significant/", locFocClust, "_", xMark, "_vs_", yMark, "gates.png"), width = 3, height = 3)
      
  }
      
  }
}
