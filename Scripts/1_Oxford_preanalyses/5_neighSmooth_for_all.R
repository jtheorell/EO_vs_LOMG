popVec <- c("B", "CD4T", "CD8T", "TCRgd", "NK", "ILC", "BT_all", "ILC_NK_all")

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
                            "CD56", "CD8",
                            "CD127", "CD27", "CD57",
                            "NKG2A", "CD2", "NKG2C",
                            "CD7","CD218a"), 
                   "ILC" = c("NKp30", "HLADR",
                             "CCR7", "CD183",
                             "CD161", "CD56", 
                             "CD127", "CD27", 
                             "NKG2A", "CD2", 
                             "CRTH2", "CD7",
                             "CD117", "CD218a"),
                   "BT_all" = c("CD38", "KLRG1", "CCR7", "CD95", "CD25", "CD45", 
                                "CD4", "CCR4", "IgD.TCRgd", "CD71","CD8", "CD19",
                                "CD3", "CD11b", "CD27", "CD20", "CD45RA","CD2",
                                "CD31", "CD161", "CD7", "CD24", "IgM", "CD152"),
                   "ILC_NK_all" = c("NKp30", "HLADR", "CCR7", "CD183", "CD161", "CD45", 
                                    "IgD.TCRgd", "CD56", "CD8", "CD19", "CD3", "CD127", 
                                    "CD27", "CD57", "NKG2A", "CD2", "NKG2C", "CRTH2", 
                                    "CD7", "CD117", "CD218a"))
library(gmodels)
library(ClusterR)
library(DepecheR)
for(a in popVec){
  print(a)
  dataDir <- paste0("../External/Oxford/Resulting_data/", a, "/")
  big_all_file <- readRDS(paste0(dataDir, a, "_full_file.rds"))
  
  #Here, we importantly change the categorisation of the 50-year old LOMG patient
  #to unclear and rename preLOMG to unclear as well, as we have no good way
  #of knowing these samples´ category. 
  big_all_file$group[which(big_all_file$group == "preLOMG")] <- "unclear"
  big_all_file$group[which(big_all_file$group == "LOMG" &
                             big_all_file$age == 50)] <- "unclear"

  #Now, we create the group vectors that are going to be used to generate the 
  #datasets below. 
  fullEomgSet <- unique(big_all_file$idTime[which(big_all_file$group == "EOMG" &
                                                    big_all_file$tissue == "pre")])
  
  fullCtrlSet <- unique(big_all_file$idTime[which(big_all_file$group == "Ctrl")])
  fullLomgSet <- unique(big_all_file$idTime[which(big_all_file$group == "LOMG" &
                                                    big_all_file$tissue == "pre")])
  fullSet <- c(fullEomgSet, fullCtrlSet, fullLomgSet)
  
  #Here, we turn the unclear into controls, but their rows will never be used
  #as neighbors anyway
  simpleGroup <- big_all_file$group
  simpleGroup[which(simpleGroup == "unclear")] <- "LOMG"
  groupDummy <- as.data.frame(matrix(0, nrow = length(simpleGroup), 
                                     ncol = 3))
  colnames(groupDummy) <- c("Ctrl", "EOMG", "LOMG")
  
  
  for(i in colnames(groupDummy)){
    groupDummy[which(simpleGroup == i),
               which(colnames(groupDummy) == i)] <- 1
  }
  
  modelDf <- big_all_file[,markerList[[a]]]
  
  #Here comes the run vector, where no samples are excluded:
  runVec <- sort(unique(big_all_file$idTime))
  
  big_all_file$pop_perc <- NA
  for(i in runVec){
    locRowVec <- which(big_all_file$idTime == i)
    big_all_file$pop_perc[locRowVec] <- 
      length(locRowVec)/big_all_file$lymphCount[locRowVec][1]
  }
  eomgRows <- which(big_all_file$idTime %in% fullEomgSet)
  ctrlRows <- which(big_all_file$idTime %in% fullCtrlSet)
  lomgRows <- which(big_all_file$idTime %in% fullLomgSet)
  
  eomgDonMinSize <- min(big_all_file$lymphCount[eomgRows])
  ctrlDonMinSize <- min(big_all_file$lymphCount[ctrlRows])
  lomgDonMinSize <- min(big_all_file$lymphCount[lomgRows])
  #Here, we calculate the size of the dataset if each donor would get the same
  #representation in each group. 
  eomgNum <- length(fullEomgSet)
  ctrlNum <- length(fullCtrlSet)
  lomgNum <- length(fullLomgSet)
  minEomgDatSize <- eomgDonMinSize*eomgNum
  minCtrlDatSize <- ctrlDonMinSize*ctrlNum
  minLomgDatSize <- lomgDonMinSize*lomgNum
  
  #So, now the smallest of these is identified
  minDatSize <- min(c(minEomgDatSize, minCtrlDatSize, minLomgDatSize))
  #206622
  
  #This value is now used to identify how many lymphocytes that should be retrieved
  #from each donor. 
  eomgDonSizeFinal <- round(minDatSize/eomgNum)
  ctrlDonSizeFinal <- round(minDatSize/ctrlNum)
  lomgDonSizeFinal <- round(minDatSize/lomgNum)
  
  #There is a case where one donor has fewer cells than this, 
  #if the number of donors in one of the groups would be much larger
  #than the others. In such a case, we will
  #resample to make up for the difference. 
  
  sampSizeList <- list(list(fullEomgSet, eomgDonSizeFinal), 
                       list(fullCtrlSet, ctrlDonSizeFinal),
                       list(fullLomgSet, lomgDonSizeFinal))
  #Importantly, we will here for each individual relate the percentage to the
  #total size of the lymphocyte population and then downsize the specific cell population
  #accordingly. The FALSE and TRUE at the end of the generated lists are for
  #replacvements below. 
  
  indivSampSizes <- unlist(lapply(sampSizeList, function(y){
    locList <- lapply(y[[1]], function(z){
      locRows <- which(big_all_file$idTime == z)
      locLymphSize <- big_all_file$lymphCount[locRows][1]
      locPopPerc <- big_all_file$pop_perc[locRows][1]
      if(locLymphSize >= y[[2]]){
        list("Cell_number" = round(y[[2]]*locPopPerc), 
             "Replace" = FALSE)
      } else {
        list("Cell_number" = round(y[[2]]*locPopPerc), 
             "Replace" = TRUE)
      }
    })
    names(locList) <- y[[1]]
    locList
    
  }), recursive = FALSE)
  
  #None get replacements, generally. 
  
  #Next thing is: how large is this neighbor set? If larger than 100000, we
  #will limit it to around this value. 
  if(sum(unlist(lapply(indivSampSizes, "[[", 1))) > 100000){
    limSize <- 100000/sum(unlist(lapply(indivSampSizes, "[[", 1)))
    indivSampSizes <- lapply(indivSampSizes, function(x){
      x$Cell_number <- round(x$Cell_number*limSize)
      x
    })
  }
  
  set.seed(119)
  modelPca <- fast.prcomp(modelDf)$x[,1:10]
  
  kMeansK <- max(1, round(nrow(modelPca)/5000))
  #Here, we set a cap on 1000, for this not to be ridiculous
  #if(kMeansK > 1000){
  #  kMeansK <- 1000
  #}
  #bef <- Sys.time()
  #modelClusts <- KMeans_rcpp(modelPca, kMeansK)
  #Sys.time()-bef
  #
  #saveRDS(modelClusts, paste0(dataDir, "modelClusts.rds"))
  #To avoid re-running this by mistake, we have commented the code above out, 
  #and instead activated the code below. Change this during the first run, of course
  
  modelClusts <- readRDS(paste0(dataDir, "modelClusts.rds"))
  
  dir.create(paste0(dataDir, "NeighSmoothModels"))
  for(i in 1:100){
    print(paste0("Here comes cycle number ", i))
    selectedRows <- unlist(lapply(names(indivSampSizes), function(y){
      locDat <- indivSampSizes[[which(names(indivSampSizes) == y)]]
      set.seed(i+100)
      sample(which(big_all_file$idTime == y), 
             size = locDat[[1]], 
             replace = locDat[[2]])
    }))
    set.seed(i+3)
    neighSmoothRes <- neighSmooth(focusData = groupDummy, 
                                  euclidSpaceData = modelPca, 
                                  neighRows = selectedRows,
                                  kNeighK = 15,
                                  kMeansCenters = modelClusts$centroids,
                                  kMeansClusters = modelClusts$clusters)
    #And this is saved, to avoid getting the computer clogged with bizarre amounts 
    #of data. 
    saveRDS(neighSmoothRes, 
            paste0(dataDir, "NeighSmoothModels/neighSmoothRes_", i, ".rds"))
    if(i > 2){
      smoothFileList <- list.files(paste0(dataDir, "NeighSmoothModels"), 
                                   full.names = TRUE)
      
      summedSmooth <- readRDS(smoothFileList[1])
      smoothResidualSumList <- list(1)
      smoothResDerivativeList <- list(1)
      for(j in 2:length(smoothFileList)){
        locFileName <- smoothFileList[j]
        locFile <- readRDS(locFileName)
        #Here, 
        smoothResidualSumList[[j]] <- sqrt(sum(abs((summedSmooth/(j-1))-
                                                     ((summedSmooth+locFile)/j))))
        if(j == 2){
          smoothResDerivativeList[[j]] <- 1
        } else {
          distFromLast <- abs(smoothResidualSumList[[j]]-smoothResidualSumList[[j-1]])
          resSumRange <- max(unlist(smoothResidualSumList))-
            min(unlist(smoothResidualSumList))
          smoothResDerivativeList[[j]] <- distFromLast/resSumRange
        }
        summedSmooth <- summedSmooth+locFile
      }
      if(smoothResDerivativeList[[i]] >= 0.01){
        print(paste0("The additional square root improvement of the residuals by",
                     "running the last round was ", 
                     round(smoothResDerivativeList[[i]]*100, 1), "%, so we continue"))
      } else {
        message("This cycle was the first to lessen the square root of the residuals", 
             "less than 1%, and therefore we stop here")
        break
      }
    }
  }
  
  #Now, we import these one by one and sum them, and then divide by 100. 
  smoothFileList <- list.files(paste0(dataDir, "NeighSmoothModels"), 
                               full.names = TRUE)
  
  summedSmooth <- readRDS(smoothFileList[1])
  for(i in smoothFileList[2:length(smoothFileList)]){
    print(i)
    locFile <- readRDS(i)
    summedSmooth <- summedSmooth+locFile
  }
  
  summedSmoothMean <- summedSmooth/length(smoothFileList)
  summedSmoothDf <- as.data.frame(summedSmoothMean)
  #And here we are!
  big_all_file$SmoothGroup <- "None"
  #We now call the high and the low events. In the case of the low, we
  #require the other two groups to be above their expected probability
  big_all_file$SmoothGroup[which(summedSmoothDf$EOMG > 0.825)] <- "EOMGhigh"
  big_all_file$SmoothGroup[which(summedSmoothDf$LOMG > 0.825)] <- "LOMGhigh"
  
  big_all_file$SmoothGroup[which( big_all_file$SmoothGroup == "None" &
                             summedSmoothDf$EOMG < 0.132 &
                             summedSmoothDf$LOMG > 0.33 &
                             summedSmoothDf$Ctrl > 0.33)] <- "EOMGlow"
  big_all_file$SmoothGroup[which( big_all_file$SmoothGroup == "None" &
                                    summedSmoothDf$LOMG < 0.132 &
                                    summedSmoothDf$EOMG > 0.33 &
                                    summedSmoothDf$Ctrl > 0.33)] <- "LOMGlow"
  
  saveRDS(big_all_file, paste0(dataDir, a, "_full_file_post_Euclid.rds"))
}