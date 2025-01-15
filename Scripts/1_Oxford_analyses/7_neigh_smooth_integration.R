library(BiocParallel)
#Here, we will take the results from the three different neighbor numbers and
#integrate them, so that only the ones that are classified as significant
#by all three are considered truly different. 
popVec <- c("B", "CD4T", "CD8T", "TCRgd", "NK", "ILC")
for(x in popVec){
  print(x)
  locFile_5 <- readRDS(paste0("../External/Oxford/Resulting_data_5_neighbors/", x, "/", x, "_full_file_post_Euclid_5_neighbors.rds"))
  locFile_15 <- readRDS(paste0("../External/Oxford/Resulting_data/", x, "/", x, "_full_file_post_Euclid.rds"))
  locFile_45 <- readRDS(paste0("../External/Oxford/Resulting_data_45_neighbors/", x, "/", x, "_full_file_post_Euclid_45_neighbors.rds"))
  if(all(identical(locFile_5[,1:10], locFile_15[,1:10]), identical(locFile_5[,1:10], locFile_45[,1:10]))){
    consensusSmoothGroup <- unlist(bplapply(1:nrow(locFile_5), function(y){
      if(all(identical(locFile_5$SmoothGroup[y],
                       locFile_15$SmoothGroup[y]),
             identical(locFile_5$SmoothGroup[y],
                       locFile_45$SmoothGroup[y]))){
        locFile_5$SmoothGroup[y]
      } else {
        "None"
      }
    }))
    locFileCommon <- locFile_5
    locFileCommon$SmoothGroup <- consensusSmoothGroup
    print(identical(locFile_5$SmoothGroup, consensusSmoothGroup))
    saveRDS(locFileCommon, paste0( "../External/Oxford/Resulting_data/", x, "/", x, "_full_file_post_Euclid_consensus.rds"))
  } else {
    stop("Something is up with the order of the events, check and get back")
  }
}