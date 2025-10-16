library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggmosaic)
library(ggplot2)
#Here, we are going to classify the selected hit cells into cell types using
#SingleR again. This time with another dataset, to increase the resolution. 
signPops <- row.names(read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1))

sceAllRaw <- readRDS("../External/Stockholm/Data/SCE_versions/5_post_exclusions.rds")

#Already here, we are going to sync the dataset with the expression dataset
tSinglerData <- MonacoImmuneData()
nkSinglerData <- HumanPrimaryCellAtlasData()

datasetList <- list(tSinglerData, tSinglerData, nkSinglerData)


#Here, we have curated two CD8T lists and one NK list for our purposes. 
NKSinglerVector <- c("NK_cell", "NK_cell:IL2", "NK_cell:CD56hiCD62L+")
CD8TSinglerVector <- c("Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells", "MAIT cells")
singlerFineLabelsList <- list(CD8TSinglerVector, CD8TSinglerVector, NKSinglerVector)

labelsNK <- readRDS("../External/Oxford_and_Stockholm/Harmonisation/NK_Stockholm_data_with_Ox_neighbors.rds")
labelsCD8T <- readRDS("../External/Oxford_and_Stockholm/Harmonisation/CD8T_Stockholm_data_with_Ox_neighbors.rds")

labelList <- list(labelsCD8T, labelsCD8T, labelsNK)


colorList <- list(rev(magma(5)), rev(magma(5)), mako(4)[2:4])

names(datasetList) <- names(labelList) <- names(singlerFineLabelsList) <- names(colorList) <- signPops
#And now we start looping. We will find MAITs, hence the name. 
dir.create("Results/Graphics/Stockholm/SingleR_and_MAIT")
singlerList <- lapply(signPops, function(x){
  print(x)
  locMother <- gsub("|_.OMGlow_.", "", x)
  locLabels <- labelList[[x]]
  locSCE <- sceAllRaw[,which(sceAllRaw$cellType == locMother)]
  
  locSingler <- datasetList[[x]]
  commonTranscripts <- Reduce(intersect, list(row.names(locSingler), rowData(locSCE)$hgnc_symbol))
  
  commonSinglerData <- locSingler[which(row.names(locSingler) %in% commonTranscripts),]
  
  sceAllCommon <- locSCE[which(rowData(locSCE)$hgnc_symbol %in% commonTranscripts),]
  row.names(sceAllCommon) <- rowData(sceAllCommon)$hgnc_symbol
  
  sceAllOrd <- sceAllCommon[match(row.names(commonSinglerData), row.names(sceAllCommon)),]
  identical(row.names(sceAllOrd), row.names(commonSinglerData)) #TRUE of course. 
  
  if(identical(row.names(locLabels), colnames(sceAllOrd))){
    locFocSCE <- sceAllOrd[,which(locLabels$Ox_hit_clust == x)]
    focSingler <- commonSinglerData[,which(commonSinglerData$label.fine %in% singlerFineLabelsList[[x]])]
    innerRes <- SingleR(locFocSCE, focSingler, 
                        labels = focSingler$label.fine, 
                        num.threads = 1)
    ggDat <- as.data.frame(innerRes)
    ggDat$labels <- factor(ggDat$labels, levels = singlerFineLabelsList[[x]])
    #Now, we are plotting
    p <- ggplot(ggDat, aes(x = 1, fill = labels)) +
      geom_bar(position="stack") + theme_void() +
      scale_fill_manual(values = colorList[[x]], drop = FALSE)
    p
    ggsave(paste0("Results/Graphics/Stockholm/SingleR_and_MAIT/", x, "_Singler_categories.pdf"))
    p <- p+theme(legend.position = "none")
    ggsave(paste0("Results/Graphics/Stockholm/SingleR_and_MAIT/", x, "_Singler_categories_stripped.pdf"), width = 1, height = 5)
    #And we save the figure data here. 
    write.csv(ggDat, paste0("Results/Data/For_figure_file/Figure_3B_SingleR_", x, ".csv"))
    innerRes
  } else {
    stop("Get your house in order before trying to analyse")
  }
  
  
})

lapply(singlerList, function(x) table(x$labels))
#Extremely nice in a way: we confirm the findings from others but in such a refined way. 
#   Central memory CD8 T cells                    MAIT cells             Naive CD8 T cells Terminal effector CD8 T cells 
#                           84                             1                          1156                             1 
#   Central memory CD8 T cells   Effector memory CD8 T cells                    MAIT cells Terminal effector CD8 T cells 
#                            9                             6                           127                             4 
#             NK_cell NK_cell:CD56hiCD62L+          NK_cell:IL2 
#                 288                    1                   27 

#No 

#So the first population is naihve, as we knew, and the third is NK. But the second
#one is MAIT. It makes so much sense! Now, what about the chains? 

singler2 <- singlerList[[2]]

signPop2ColData <- colData(sceAllRaw[,which(colnames(sceAllRaw) %in% row.names(singler2))])
table(signPop2ColData$TCR_Japha33, signPop2ColData$TCR_Vapha7.2)
#FALSE TRUE
#FALSE    79    6
#TRUE      1   30

#I might not have caught the right associated beta chains, because there is no overlap
#with the alpha chains. But they also seem very over-represented.
table(signPop2ColData$TCR_Vbeta13, signPop2ColData$TCR_Vbeta2)
#        FALSE TRUE
#  FALSE    76   17
#  TRUE     23    0

#This is the background population. 
table(sceAllRaw$TCR_Japha33, sceAllRaw$TCR_Vapha7.2)

#        FALSE  TRUE
#  FALSE 25726   137
#  TRUE    184   110
110/(25726+137+184)
#0.004223135
#And for the beta chains: 
table(sceAllRaw$TCR_Vbeta13, sceAllRaw$TCR_Vbeta2)
#        FALSE  TRUE
#  FALSE 22624  1683
#  TRUE   1850     0

#So in fact, this is extremely significant: 
fisher.test(matrix(c(30, ncol(signPop2ColData), c(110, nrow(sceAllRaw))), 2,2))
#p-value smaller than 2.2e-16. 

#As we do know from another set of graphs that
#there are signficant differences also in terms of the TCRVbeta family distribution, 
#but that clonality is stable, we are going to integrate all these factors here and
#include the first population as well, for good measure. 

CD8ColDat <- colData(sceAllRaw[,which(sceAllRaw$TCR & sceAllRaw$cellType == "CD8T")])
labelsCD8TTCR <- labelsCD8T[which(rownames(labelsCD8T) %in% rownames(CD8ColDat)),]
#Just checking
identical(rownames(CD8ColDat), rownames(labelsCD8TTCR)) #TRUE
freqTabList <- list()
for(i in c(signPops[1:2], "Ref")){
  if(i == "Ref"){
    locColDat <- CD8ColDat
  } else {
    locColDat <- CD8ColDat[which(labelsCD8TTCR$Ox_hit_clust == i),]
  }
  longTab <- as.data.frame(table(locColDat$TCR_Vapha7.2, locColDat$TCR_Japha33, locColDat$TCR_Clonal, locColDat$TCR_vFamB))
  freqTab <- longTab
  freqTab$Freq <- longTab$Freq/sum(longTab$Freq)
  colnames(freqTab)[1:4] <- c("Valpha7.2", "Jalpha33", "TCR_Clonal", "TCRVB_fam")
  freqTab[,1:3] <- apply(freqTab[,1:3], 2, as.logical)
  #Now, we combine the first three into one. 
  freqTab$clonalAndAlpha <- "NonClonal_None"
  freqTab$clonalAndAlpha[which(freqTab$Valpha7.2 == FALSE &
                                 freqTab$Jalpha33 == FALSE &
                                 freqTab$TCR_Clonal)] <- "Clonal_None"
  freqTab$clonalAndAlpha[which(freqTab$Valpha7.2 == TRUE &
                                 freqTab$Jalpha33 == FALSE &
                                 freqTab$TCR_Clonal == FALSE)] <- "NonClonal_Valpha7.2"
  freqTab$clonalAndAlpha[which(freqTab$Valpha7.2 == TRUE &
                                 freqTab$Jalpha33 == FALSE &
                                 freqTab$TCR_Clonal == TRUE)] <- "Clonal_Valpha7.2"
  freqTab$clonalAndAlpha[which(freqTab$Valpha7.2 == FALSE &
                                 freqTab$Jalpha33 == TRUE &
                                 freqTab$TCR_Clonal == FALSE)] <- "NonClonal_Jalpha33"
  freqTab$clonalAndAlpha[which(freqTab$Valpha7.2 == FALSE &
                                 freqTab$Jalpha33 == TRUE &
                                 freqTab$TCR_Clonal == TRUE)] <- "Clonal_Jalpha33"
  freqTab$clonalAndAlpha[which(freqTab$Valpha7.2 == TRUE &
                                 freqTab$Jalpha33 == TRUE &
                                 freqTab$TCR_Clonal == FALSE)] <- "NonClonal_MAIT"
  freqTab$clonalAndAlpha[which(freqTab$Valpha7.2 == TRUE &
                                 freqTab$Jalpha33 == TRUE &
                                 freqTab$TCR_Clonal == TRUE)] <- "Clonal_MAIT"
  freqTab$clonalAndAlpha <- factor(freqTab$clonalAndAlpha, 
                                   levels = c("NonClonal_None", 
                                              "Clonal_None", 
                                              "NonClonal_Valpha7.2", 
                                              "Clonal_Valpha7.2", 
                                              "NonClonal_Jalpha33",
                                              "Clonal_Jalpha33", 
                                              "NonClonal_MAIT", 
                                              "Clonal_MAIT"))
  p <- ggplot(freqTab, aes(x=clonalAndAlpha, y=Freq, fill=TCRVB_fam)) + 
    geom_bar(position="stack", stat="identity", , show.legend = T) + theme_bw() +
    scale_fill_manual(values = cividis(8), drop = FALSE) +
    scale_y_continuous(limits = c(0,0.55), expand = c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p
  ggsave(paste0("Results/Graphics/Stockholm/SingleR_and_MAIT/", i, "_clonal_vs_MAIT_vs_BVFam.pdf"))
  p <- p+theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.text.y = element_blank())
  p
  ggsave(paste0("Results/Graphics/Stockholm/SingleR_and_MAIT/", i, "_clonal_vs_MAIT_vs_BVFam_stripped.pdf"), width = 5, height = 5)
  #And we save the csv 
  write.csv(freqTab, paste0("Results/Data/For_figure_file/Figure_3F_AIRR_", i, ".csv"), row.names = FALSE)
  
  freqTabList[[i]] <- c("Freqs" = freqTab, "Counts" = longTab)
}

#Now, a few tests. 
fisherMat <- do.call("rbind", lapply(freqTabList[c(2,3)], function(x){
  c(sum(x$Counts.Freq), 
    x$Counts.Freq[which(x$Freqs.clonalAndAlpha == "NonClonal_MAIT" & x$Freqs.TCRVB_fam == "TRBV6")])
}))

#               [,1] [,2]
#CD8T_LOMGlow_2  111   12
#Ref            6534   3
12/123
#0.09756098
3/6537
#0.0004589261

fisher.test(fisherMat)
#p-value = 1.987e-11

fisherMat <- do.call("rbind", lapply(freqTabList[c(2,3)], function(x){
  c(sum(x$Counts.Freq[which(x$Freqs.clonalAndAlpha == "NonClonal_MAIT" & x$Freqs.TCRVB_fam != "TRBV6")]), 
    x$Counts.Freq[which(x$Freqs.clonalAndAlpha == "NonClonal_MAIT" & x$Freqs.TCRVB_fam == "TRBV6")])
}))

fisher.test(fisherMat)

#p-value = 0.5718, the cell numbers are simply too low to allow for true statistical comparisons. 
