library(viridis)
library(DepecheR)
library(uwot)
#Here, we create plots for the first figure, which is an extended flow chart. 
neighborised_Stock <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/CD8T_Stockholm_data_with_Ox_neighbors.rds"))

oxDat <- readRDS(paste0("../External/Oxford/Resulting_data/CD8T/CD8T_full_file_post_Euclid_including_subclusts.rds"))
#We start by zooming in on the cells of relevance, i.e. the pre-thymectomy
#PBMC samples from patients. 
locDat <- oxDat[which(oxDat$tissue == "pre" &
                         oxDat$group != "Ctrl"),]
rawPCA <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/CD8T_raw_PCA.rds"))
harmonyPCA <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/CD8T_harmony_PCA.rds"))

oxRows <- 1:nrow(locDat)
stockRows <- (nrow(locDat)+1):(nrow(locDat)+nrow(neighborised_Stock))

#We will downsample the Oxford data here, but still
#keep a considerable difference between them. This is only done for the CD8
#as it takes some time, and is only used for the example in figure 1. 

set.seed(11)
focOxRows <- sample(oxRows, length(stockRows)*20)
umapRows <- c(focOxRows, stockRows)
#rawUmap <- umap(rawPCA[umapRows,])
#saveRDS(rawUmap, "Results/Data/Oxford_and_Stockholm/Harmonisation/CD8T_UMAP_raw.rds")
rawUmap <- readRDS("Results/Data/Oxford_and_Stockholm/Harmonisation/CD8T_UMAP_raw.rds")

#harmonyUmap <- umap(harmonyPCA[umapRows,])
#saveRDS(harmonyUmap, "Results/Data/Oxford_and_Stockholm/Harmonisation/CD8T_UMAP.rds")
harmonyUmap <- readRDS("Results/Data/Oxford_and_Stockholm/Harmonisation/CD8T_UMAP.rds")
#We create two simplified plotting vectors here, where all EOMG and LOMG
#defining cells are given one colour each.
oxHitVec <- locDat$smoothGroupDeep
plotEuclidClusterList <- list("Ox" = oxHitVec[focOxRows], 
                              "Stock" = neighborised_Stock$Ox_hit_clust)

oxPlotRows <- 1:length(focOxRows)
stockPlotRows <- (length(focOxRows)+1):(length(focOxRows)+length(stockRows))

#Now, we are going to change the order, to get the None plotted first, at the
#bottom.
oxPlotOrder <- order(plotEuclidClusterList[[1]], decreasing = TRUE)
stockPlotOrder <- order(plotEuclidClusterList[[2]], decreasing = TRUE)

plotEuclidClusterListOrd <- plotEuclidClusterList
plotEuclidClusterListOrd[[1]] <- plotEuclidClusterList[[1]][oxPlotOrder]
plotEuclidClusterListOrd[[2]] <- plotEuclidClusterList[[2]][stockPlotOrder]

#Here, the names are changed for the stockholm data, to remove the prefix "CD8T"
plotEuclidClusterListOrd[[2]] <- gsub("CD8T_", "", plotEuclidClusterListOrd[[2]])
#Here, we create a color vector for the subclusters. 
subclustColVecShort <- c(dColorVector(1:5, colors = c("white",  "#440154FF", "black"))[2:4],
                         dColorVector(1:4, colors = c("white",  "#35B779FF", "black"))[2:3],
                         dColorVector(1:5, colors = c("white",  "#FDE725FF", "black"))[2:4],
                         dColorVector(1:4, colors = c("white",  "#31688EFF", "black"))[2:3],
                         "grey")
names(subclustColVecShort) <- sort(unique(plotEuclidClusterList[[1]]))

subclustColListLong <- lapply(plotEuclidClusterListOrd, function(x){
  unlist(sapply(x, function(y){
    unname(subclustColVecShort[which(names(subclustColVecShort) == y)])
  }))
})

#We change the grey color for the stockholm data somewhat
subclustColListLong[[2]][which(subclustColListLong[[2]] == "grey")] <- "#777777"

plotEuclidSimpleList <- lapply(plotEuclidClusterListOrd, function(x){
  gsub("^.{1,5}_|_.$", "", x)
})
rawUmapOrd <- rawUmap
rawUmapOrd[oxPlotRows,] <- rawUmapOrd[oxPlotRows,][oxPlotOrder,]
rawUmapOrd[stockPlotRows,] <- rawUmapOrd[stockPlotRows,][stockPlotOrder,]

harmonyUmapOrd <- harmonyUmap
harmonyUmapOrd[oxPlotRows,] <- harmonyUmapOrd[oxPlotRows,][oxPlotOrder,]
harmonyUmapOrd[stockPlotRows,] <- harmonyUmapOrd[stockPlotRows,][stockPlotOrder,]

dir.create(paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j))
#Step one is showing the neighbor vectors, so we have to import those. 
smoothFileList <- list.files("../External/Oxford/Resulting_data/CD8T/NeighSmoothModels", 
                             full.names = TRUE)

summedSmooth <- readRDS(smoothFileList[1])
for(i in smoothFileList[2:length(smoothFileList)]){
  print(i)
  locFile <- readRDS(i)
  summedSmooth <- summedSmooth+locFile
}

summedSmoothMean <- summedSmooth/length(smoothFileList)
summedSmoothDf <- as.data.frame(summedSmoothMean)

locSmoothDf <- summedSmoothDf[which(oxDat$tissue == "pre" &
                                      oxDat$group != "Ctrl"),]
eomgProb <- locSmoothDf$EOMG[focOxRows][oxPlotOrder]
lomgProb <- locSmoothDf$LOMG[focOxRows][oxPlotOrder]
ctrlProb <- locSmoothDf$Ctrl[focOxRows][oxPlotOrder]
j <- "CD8T"
dColorPlot(eomgProb, 
           xYData = rawUmapOrd[oxPlotRows,], 
           colors = c("#35B779FF","grey", "grey", "#440154FF"), dotSize = 3,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, 
                             "/", j, "_1a_Oxford_umap_EOMG_smooth"))
dColorPlot(lomgProb, 
           xYData = rawUmapOrd[oxPlotRows,], 
           colors = c("#31688EFF","grey", "grey", "#FDE725FF"), dotSize = 3,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, 
                             "/", j, "_1b_Oxford_umap_LOMG_smooth"))

dColorPlot(ctrlProb, 
           xYData = rawUmapOrd[oxPlotRows,], 
           colors = c("orange","grey", "grey", "black"), dotSize = 3,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, 
                             "/", j, "_1c_Oxford_umap_ctrl_smooth"))


dColorPlot(plotEuclidSimpleList[[1]], 
           xYData = rawUmapOrd[oxPlotRows,], 
           colors = c("#440154FF", "#35B779FF", "#FDE725FF", "#31688EFF", "grey"), dotSize = 3,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, "/", j, "_2_Oxford_umap_raw"))

dColorPlot(plotEuclidClusterListOrd[[1]], 
           xYData = rawUmapOrd[oxPlotRows,], 
           colors = subclustColListLong[[1]], dotSize = 3,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, "/", j, "_3a_Oxford_umap_raw_subclusts"))

dColorPlot(plotEuclidClusterListOrd[[1]], 
           xYData = harmonyUmapOrd[oxPlotRows,], 
           colors = subclustColListLong[[1]], dotSize = 3,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, "/", j, "_3_Oxford_umap_harmony_subclusts"))

dColorPlot(plotEuclidClusterListOrd[[2]], 
           xYData = harmonyUmapOrd[stockPlotRows,], 
           colors = subclustColListLong[[2]], dotSize = 6,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, "/", j, "_4_Stockholm_umap_harmony_subclusts"))

#Now, we also create two grey plots, that will be combined in the figure. 
dColorPlot(rep("Ox", length(plotEuclidClusterListOrd[[1]])), 
           xYData = harmonyUmapOrd[oxPlotRows,], 
           colors = c("grey", "black"), dotSize = 3,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, "/", j, "_5_Oxford_umap_harmony_no_colors"))

dColorPlot(rep("Stock", length(plotEuclidClusterListOrd[[2]])), 
           xYData = harmonyUmapOrd[stockPlotRows,], 
           colors = c("#777777", "black"), dotSize = 6,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, "/", j, "_6_Stockholm_umap_harmony_no_colors"))

#Finally, we only show the one population that is significant. This is in
#a way illogical, as that should not be known at this stage of the analysis,
#but this change is introduced post-hoc. 
hitEuclidClusterList <- plotEuclidClusterListOrd
hitEuclidClusterList[[1]][-grep("LOMGlow", hitEuclidClusterList[[1]])] <- "None"
hitEuclidClusterList[[2]][-grep("LOMGlow", hitEuclidClusterList[[2]])] <- "None"

dColorPlot(hitEuclidClusterList[[1]], 
           xYData = harmonyUmapOrd[oxPlotRows,], 
           colors = c(subclustColVecShort[grep("LOMGlow", names(subclustColVecShort))], "grey"), dotSize = 3,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, "/", j, "_7_Oxford_umap_harmony_hit_pop"))

dColorPlot(hitEuclidClusterList[[2]], 
           xYData = harmonyUmapOrd[stockPlotRows,], 
           colors = c(subclustColVecShort[which(names(subclustColVecShort) == "LOMGlow_2")], "#777777"), dotSize = 6,
           plotName = paste0("Results/Graphics/Oxford_and_Stockholm/Harmonisation/", j, "/", j, "_8_Stockholm_umap_harmony_hit_pop"))

