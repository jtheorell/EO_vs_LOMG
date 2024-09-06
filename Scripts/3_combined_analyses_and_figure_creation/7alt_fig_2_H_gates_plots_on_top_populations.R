#Here, we will run a two-step optimisation. First, we will see which levels of gates
#for all the markers of interest for each cell type that can be used to create an 
#as oure replica of our populations as possible. Then, we will take the markers
#and see which the smallest set of markers is that can be used to do this. 
#As this analysis is heavy, the marker sets are slightly reduced compared to
#previous sets. 

library(ggplot2)
library(dfoptim)
library(reshape2)

#The 
pValFoc <- 
  read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1)

popVec <- c("B", "CD4T", "CD8T", "NK")

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
                              "CD24", "IgM", "CD31"),
                      "CD4T" = c("CCR7", "CCR4","CD27", "CD45RA", 
                                 "CD31", "CD2"),
                      "CD8T" = c("KLRG1", "CCR7", 
                                 "CD27",  "CD45RA", "CD7", "CD161"),
                      "NK" = c("NKp30","CD16", "HLADR", "CD2",
                               "CD56", "CD57", "CD183",
                               "NKG2C"))

#Now, we will do something reasonably simple. We will first identify the 99th 
#and the 1st percentile for each marker for the population of interest. We
#will then rank all these markers based on their "uniqueness" compared to the
#rest of the population, i.e. where the ratio of our population to the 
#rest of the population is the best. Then we take these, one-by-one, and run a 
#Wilcoxon rank sum test, and we do this in a loop that breaks when we reach 
#significance. 

bestPvals <- list()
selectedMarkers <- list()
for(i in popVec){
  datDir <- paste0("../External/Oxford/Resulting_data/", i, "/")
  locFile <- readRDS(paste0(datDir, i, "_full_file_post_Euclid_including_subclusts.rds"))
  locFilePre <- locFile[which(locFile$tissue == "pre"),]
  locFocClust <- row.names(pValFoc)[grep(i, row.names(pValFoc))]
  print(i)
  print(locFocClust)
  focRows <- which(locFilePre$smoothGroupDeep %in% gsub("^.{1,4}_|", "", locFocClust))
  focPop <- locFilePre[focRows,markerList[[i]]]
  focPopMeta <- locFilePre[focRows,c("idTime", "group", "age", "lymphCount")]
  altPop <- locFilePre[-focRows,markerList[[i]]]
  #Now, to construct the gates, we will sweep from inclusion of the full focus
  #population mid two quartiles
  lowSeq <- seq(0,0.49, by = 0.01)
  lowerBoundDf <- do.call("rbind", lapply(colnames(focPop), function(x){
    
    locRatios <- unlist(lapply(lowSeq, function(y){
      lowerBound <- quantile(focPop[,x], y)
      upperBound <- max(focPop[,x])
      fracOfAlt <- length(which(altPop[,x] >= lowerBound &
                                  altPop[,x] < upperBound))/nrow(altPop)
      (1-y)/fracOfAlt
    }))
    bestQuantile <- lowSeq[which.max(locRatios)]
    c("quantGate" = unname(quantile(focPop[,x], bestQuantile)), 
      "fracRatio" = max(locRatios), 
      "quantile" = bestQuantile)
  }))
  highSeq <- seq(0.5,1, by = 0.01)
  upperBoundDf <- do.call("rbind", lapply(colnames(focPop), function(x){
    
    locRatios <- unlist(lapply(highSeq, function(y){
      lowerBound <- min(focPop[,x])
      upperBound <- quantile(focPop[,x], y)
      fracOfAlt <- length(which(altPop[,x] >= lowerBound &
                                  altPop[,x] < upperBound))/nrow(altPop)
      (y)/fracOfAlt
    }))
    bestQuantile <- highSeq[which.max(locRatios)]
    c("quantGate" = unname(quantile(focPop[,x], bestQuantile)), 
      "fracRatio" = max(locRatios), 
      "quantile" = bestQuantile)
  }))
  
  lowerBoundVec <- lowerBoundDf[,1]
  upperBoundVec <- upperBoundDf[,1]
  names(lowerBoundVec) <- names(upperBoundVec) <- colnames(focPop)
  #Now, we calculate the fraction of the non-focus cells within this region of 
  #interest compared to the fraction of the focus cells and then rank the markers 
  #based on this. 
  altPopFracs <- sapply(colnames(altPop), function(x){
    length(which(altPop[,x] > lowerBoundVec[x] &
                   altPop[,x] < upperBoundVec[x]))/nrow(altPop)
  })
  popFracs <- sapply(colnames(focPop), function(x){
    length(which(focPop[,x] > lowerBoundVec[x] &
                   focPop[,x] < upperBoundVec[x]))/nrow(focPop)
  })
  fracLeverage <- popFracs/altPopFracs
  #So now, the markers are ranked. 
  fracLeveragesOrd <- fracLeverage[order(fracLeverage, decreasing = TRUE)]
  
  #And now for the while loop
  pPat <- 1
  pCtrl <- 1
  numOfMarkers <- 0
  while(numOfMarkers < length(fracLeveragesOrd) && (pPat > 0.05 || pCtrl > 0.05)){
    numOfMarkers <- numOfMarkers+1
    usedMarkers <- names(fracLeveragesOrd)[1:numOfMarkers]
    gatedRows <- Reduce(intersect, lapply(usedMarkers, function(x){
          which(locFilePre[,x] >= lowerBoundVec[x] & locFilePre[,x] < upperBoundVec[x])
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
    if(grepl("EOMG", row.names(pValFoc)[grep(i, row.names(pValFoc))])){
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
    if(grepl("EOMG", row.names(pValFoc)[grep(i, row.names(pValFoc))])){
      pPat <- wilcox.test(eomgFreqs, lomgFreqs, alternative = "less")$p.value
      pCtrl <- wilcox.test(eomgFreqs, ctrlFreqs, alternative = "less")$p.value
    } else {
      pPat <- wilcox.test(lomgFreqs, eomgFreqs, alternative = "less")$p.value
      pCtrl <- wilcox.test(lomgFreqs, ctrlFreqs, alternative = "less")$p.value
    }
  }
  bestPvals[[i]] <- c("patPat" = pPat, "patCtrl" = pCtrl)
  if(all(bestPvals[[i]] < 0.05)){
    selectedMarkers[[i]] <- usedMarkers
  } else {
    selectedMarkers[[i]] <- "None"
  }
  
  #Now, we are going to select the markers of interest for each population. 
  #We will do this by taking the most different markers, based on our gating
  #and ordering analysis above. We will make a supplementary showing all markers
  #but focus on the top ones in the figure. 
  
  focMarks <- names(fracLeveragesOrd)[which(names(fracLeveragesOrd) %in% focMarkerlist[[i]])]
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
  if(selectedMarkers[[i]][1] != "None"){
    for(j in selectedMarkers[[i]]){
      locDatLong$minVals[which(locDatLong$variable == j)] <- 
        lowerBoundVec[j]
    }
    for(j in selectedMarkers[[i]]){
      locDatLong$maxVals[which(locDatLong$variable == j)] <- 
        upperBoundVec[j]
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
  
  if(grepl("EOMG", row.names(pValFoc)[grep(i, row.names(pValFoc))])){
    colorForPlot <- "red4"
  } else {
    colorForPlot <- "blue4"
  }

  p <- ggplot(locDatLong, aes(
    x = variable,
    y = value,
    fill = Population)) +
    geom_split_violin(scale = "width") +
    theme_bw() + scale_fill_manual(values = c(colorForPlot, "grey")) +
    scale_y_continuous(limits = c(0, 1050), expand = c(0,0)) 
  if(selectedMarkers[[i]][1] != "None"){
    p <- p + geom_ribbon(locDatLong, mapping = aes(x = numVariable,
                                                  ymin = minVals, 
                                                  ymax = maxVals),
                        alpha = 0.3, inherit.aes = FALSE)
  }
    
  ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                i, ".pdf"), plot = p)
  p <- p+theme(legend.position = "None", 
               axis.title = element_blank())
  ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                i, "_stripped.pdf"), plot = p, 
         width = 4, height = 2)
  #As this does not work as intended for the "gates", these being corrupted
  #repeatedly when importing into illustrator, I will also create a png output 
  #for this purpose. 
  p <- p+theme_void() + theme(legend.position = "none")
  ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                i, "_stripped_bare.png"), plot = p, 
         width = 4, height = 2)
  
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
    theme_bw() + scale_fill_manual(values = c(colorForPlot, "grey")) +
    scale_y_continuous(limits = c(0, 1050), expand = c(0,0)) 
  
  
  ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                i, "_all_markers.pdf"), plot = p)
  p <- p+theme(legend.position = "None", 
               axis.title = element_blank())
  ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                i, "_all_markers_stripped.pdf"), plot = p, 
         width = 7, height = 2)
  p <- p+theme_void() + theme(legend.position = "none")
  ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant/", 
                i, "_all_markers_stripped_bare.png"), plot = p, 
         width = 4, height = 2)
}

bestPvals
#$B
#     patPat     patCtrl 
#0.037730199 0.003026506 
#
#$CD4T
#    patPat    patCtrl 
#0.02239038 0.01275888 
#
#$CD8T
#     patPat     patCtrl 
#0.041615062 0.002067874 
#
#$NK
#     patPat     patCtrl 
#0.003317116 0.049308792

selectedMarkers
#$B
#[1] "CD38" "IgM"  "CD24" "CD31"
#
#$CD4T
#[1] "CCR4" "CD2" 
#
#$CD8T
#[1] "CD45RA" "CD7"    "KLRG1"  "CCR7"  
#
#$NK
#[1] "CD57"  "HLADR" "CD2"   "CD16"  "CD56"  "CD183" "NKp30"

