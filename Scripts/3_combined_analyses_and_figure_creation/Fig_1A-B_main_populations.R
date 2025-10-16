#Here, we are simply going to compare the main populations with the same kind
#of dot plots that we have for all the other populations downstream. 
library(ggplot2)
library(ggforce)
stockMeta <- read.csv("Data/Stockholm/Metadata/Metadata_combined_used.csv")

stockholmColDat <- 
  colData(readRDS("../External/Stockholm/Data/SCE_versions/5_post_exclusions.rds"))

stockholmCountDat <- 
  as.data.frame.matrix(table(stockholmColDat$sample, stockholmColDat$cellType))

stockholmFreqDat <-stockholmCountDat/rowSums(stockholmCountDat) 

#Now, we need to set the groups here. 
stockholmFreqDat$idTime <- row.names(stockholmFreqDat)

stockholmFreqDat$Group <- sapply(stockholmFreqDat$idTime, function(x){
  locRow <- which(stockMeta$Sample_ID == x)
  if(stockMeta$MG_type[locRow] == "LOMG" &
     (stockMeta$Thymectomy[locRow] | stockMeta$Age[locRow] == 50)){
    "uncertain"
  } else if(stockMeta$MG_type[locRow] == "EOMG" & 
            stockMeta$Years_since_thymectomy[locRow] > 0){
    "EOMG_post_thym"
  } else {
    stockMeta$MG_type[locRow]
  }
})

stockholmFreqDat$Sex <- sapply(stockholmFreqDat$idTime, function(x){
  locRow <- which(stockMeta$Sample_ID == x)
  stockMeta$Sex[locRow]
})

#Now, we make it long
stockholmFreqDatLong <- reshape2::melt(stockholmFreqDat, value.name = "Freq")
colnames(stockholmFreqDatLong)[which(colnames(stockholmFreqDatLong) == "variable")] <- "CellType"


#Here, we remove the DPT and DNT, as they are so very few. 
stockholmFreqDatLong <- stockholmFreqDatLong[-which(stockholmFreqDatLong$CellType %in% c("DNT", "DPT")),]

#And finally, before going to Oxford, we change the name o the gdT pop
stockholmFreqDatLong$CellType <- as.character(stockholmFreqDatLong$CellType)
stockholmFreqDatLong$CellType[which(stockholmFreqDatLong$CellType == "gdT")] <- "TCRgd"

#And now to the slightly more convoluted creation of the Oxford freqs.
popVec <- list("B", "CD4T", "CD8T", "TCRgd",
                   "NK", "ILC")

oxfordFreqDatLong <- do.call("rbind", lapply(popVec, function(x){
  print(x)
  big_all_file <- readRDS(paste0("../External/Oxford/Resulting_data/", x, "/", x,  "_full_file_post_Euclid_consensus.rds"))
  
  big_all_file$group[which(big_all_file$group == "preLOMG")] <- "uncertain"
  big_all_file$group[which(big_all_file$group == "LOMG" &
                             big_all_file$age == 50)] <- "uncertain"
  #We also divide the controls here. 
  big_all_file$group[which(big_all_file$group == "Ctrl" &
                             big_all_file$age < 50)] <- "Young_ctrl"
  big_all_file$group[which(big_all_file$group == "Ctrl" &
                             big_all_file$age > 50)] <- "Old_ctrl"
  
  #And we also introdce the post thym category to make this complete. 
  big_all_file$group[grep("post", big_all_file$idTime)] <- "EOMG_post_thym"
  #Ad here, we remove hte thymus samples
  big_all_file <- big_all_file[-which(big_all_file$tissue == "thy"),]
  
  freqDf <- do.call("rbind", lapply(unique(big_all_file$group), function(y){
    locSampList <- unique(big_all_file$idTime[which(big_all_file$group == y)])
    locRes <- unlist(lapply(locSampList, function(z){
        iPoss <- which(big_all_file$idTime == z)
        locCount <- length(iPoss)
        locFreq <- locCount/big_all_file$lymphCount[iPoss][1]
    }))
    locResDf <- data.frame("CellType" = x, 
                           "Freq" = locRes, 
                           "Group" = y,
                           "idTime" = locSampList)
    locResDf$Sex <- sapply(locResDf$idTime, function(x){
      big_all_file$sex[which(big_all_file$idTime == x)][1]
    })
    locResDf
  }))
}))
dir.create("Results/Data/For_figure_file")

write.csv(oxfordFreqDatLong, "Results/Data/For_figure_file/Figure_1A_UK.csv", row.names = FALSE)
write.csv(stockholmFreqDatLong, "Results/Data/For_figure_file/Figure_1B_SE.csv", row.names = FALSE)

#And now, we are ready to create the figures. 
freqTabList <- list("Stock" = stockholmFreqDatLong, "Ox" = oxfordFreqDatLong)
dir.create("Results/Graphics/Oxford_and_Stockholm/Freq_of_main_pops")
for(i in names(freqTabList)){
  freqTabLong <- freqTabList[[i]]
  freqTabLong$Group <- factor(freqTabLong$Group, 
                              levels = c("EOMG",
                                         "EOMG_post_thym",
                                         "Young_ctrl",
                                         "LOMG",
                                         "Old_ctrl", 
                                         "uncertain"))
  freqTabLong$CellType <- factor(freqTabLong$CellType, 
                              levels = c("B",
                                         "CD4T",
                                         "CD8T",
                                         "TCRgd",
                                         "NK", 
                                         "ILC"))
  freqTabLong$log10Percent <- log10(freqTabLong$Freq*100)
  
  #If any population frequencies were 0 before, they will be infinite now, 
  #so they are given a value here.
  if(length(which(is.infinite(freqTabLong$log10Percent))) > 0){
    freqTabLong$log10Percent[which(is.infinite(freqTabLong$log10Percent))] <- 
      -2.5
  } 
  
  if(i == "Stock"){
    widthForGraph <- 7
  } else {
    widthForGraph <- 14
  }
  
  sexShape <- sapply(freqTabLong$Sex, function(x) {
    switch(x, "F" =16, "M" = 17) })
  freqTabLong$sexShape <- sexShape[order(freqTabLong$CellType, freqTabLong$Group)]
  set.seed(10)
  p <- ggplot(freqTabLong, aes(x = CellType,
                               y = log10Percent,
                               fill = Group, 
                               color = Group)) +
    geom_sina(jitter_y = FALSE, 
              scale = "width", 
              maxwidth = 0.85, size = 4,
              shape = freqTabLong$sexShape) +
    geom_violin(alpha = 0, scale = "width", width = 0.85, linewidth = 1, color = "black") +
    theme_bw() + scale_color_manual(values = c("red", "red4", "grey", "blue", "#555555", "#AA00FF"),
                                    drop = FALSE)+
    scale_y_continuous(limits = c(-2.1, 2), expand = c(0,0))
  p
  ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freq_of_main_pops/", i, "main_pop_freqs.pdf"))
  
  p <- p + theme(legend.position = "none",
                 axis.text = element_blank(),
                 axis.title = element_blank())
  
  p
  ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freq_of_main_pops/", i, "main_pop_freqs_stripped.pdf"), 
         width = widthForGraph, height = 4)
}

#And now for the stats: 
pValList <- unlist(lapply(names(freqTabList), function(i){
  freqTabLong <- freqTabList[[i]]
  #Here, we merge the pre- and post thymectomy samples for the Stockholm dataset
  if(length(unique(freqTabLong$idTime[which(freqTabLong$Group == "EOMG")])) < 4){
    freqTabLong$Group[which(freqTabLong$Group == "EOMG_post_thym")] <- "EOMG"
    pVals <- sapply(unique(freqTabLong$CellType), function(x){
      locType <- freqTabLong[which(freqTabLong$CellType == x),]
      wilcox.test(locType$Freq[which(locType$Group == "EOMG")], 
                  locType$Freq[which(locType$Group == "LOMG")], 
                  exact = FALSE)$p.value
    })
    list("Stockholm" = pVals)
  } else {
    pValsPat <- sapply(unique(freqTabLong$CellType), function(x){
      locType <- freqTabLong[which(freqTabLong$CellType == x),]
      wilcox.test(locType$Freq[which(locType$Group == "EOMG")], 
                  locType$Freq[which(locType$Group == "LOMG")], 
                  exact = FALSE)$p.value
    })
    pValsEOMG <- sapply(unique(freqTabLong$CellType), function(x){
      print(x)
      locType <- freqTabLong[which(freqTabLong$CellType == x),]
      wilcox.test(locType$Freq[which(locType$Group == "EOMG")], 
                  locType$Freq[which(locType$Group == "Young_ctrl")], 
                  exact = FALSE)$p.value
    })
    pValsLOMG <- sapply(unique(freqTabLong$CellType), function(x){
      locType <- freqTabLong[which(freqTabLong$CellType == x),]
      wilcox.test(locType$Freq[which(locType$Group == "EOMG")], 
                  locType$Freq[which(locType$Group == "Old_ctrl")], 
                  exact = FALSE)$p.value
    })
    list("EOMG_LOMG" = pValsPat, 
         "EOMG_ctrl" = pValsEOMG, 
         "LOMG_ctrl" = pValsLOMG)
  }
}), recursive = FALSE)

lapply(pValList, p.adjust, method = "fdr", n = length(unlist(pValList)))