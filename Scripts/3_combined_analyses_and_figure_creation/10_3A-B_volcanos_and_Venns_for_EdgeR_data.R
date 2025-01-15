logSCE <- readRDS("../External/Stockholm/Data/SCE_versions/6_with_Ox_neigh_clust.rds")
signRNA <- read.csv("Results/Data/Stockholm/Top_transcripts_EdgeR.csv", row.names = 1)

edgerFullOutcomes <- readRDS("Results/Data/Stockholm/EdgeR_all_transcriptomes.rds")
edgerSignOutcomes <- readRDS("Results/Data/Stockholm/EdgeR_significant_transcriptomes.rds")

#How many genes do we have? 
lapply(edgerSignOutcomes, function(x) length(which(x$logFC < 0)))
#$CD8T_LOMGlow_1
#[1] 356
#
#$CD8T_LOMGlow_2
#[1] 75
#
#$NK_EOMGlow_1
#[1] 25

lapply(edgerSignOutcomes, function(x) length(which(x$logFC > 0)))
#$CD8T_LOMGlow_1
#[1] 400
#
#$CD8T_LOMGlow_2
#[1] 40
#
#$NK_EOMGlow_1
#[1] 13

library(ggplot2)
library(viridis)
dir.create("Results/Graphics/Stockholm/Volcano")
plasmaCols <- plasma(4)
for(i in names(edgerFullOutcomes)){
  locDat <- edgerFullOutcomes[[i]]
  locSign <- edgerSignOutcomes[[i]]
  locDat$Group <- "None"
  locDat$Group[which(locDat$logFC > 0 &
                        row.names(locDat) %in% row.names(locSign))] <- "Up"
  locDat$Group[which(locDat$logFC < 0 &
                       row.names(locDat) %in% row.names(locSign))] <- "Down"
  #if(grepl("LOMG", i)){
  #  locDat$Group[which(row.names(locDat) %in% eomgLowGeneList$Gene &
  #                       locDat$Group == "Up")] <- "LOMGlow_common_up"
  #  locDat$Group[which(row.names(locDat) %in% eomgLowGeneList$Gene &
  #                       locDat$Group == "Down")] <- "LOMGlow_common_down"
  #  #And now, we reorder them, so that the common ones are on top. 
  #  locDat <- rbind(locDat[-grep("OMG", locDat$Group),],
  #                  locDat[grep("OMG", locDat$Group),])
  #} #else {
    #locDat$Group[which(row.names(locDat) %in% lomgLowGeneList$Gene &
    #                     locDat$Group == "Up")] <- "LOMGlow_common_up"
    #locDat$Group[which(row.names(locDat) %in% lomgLowGeneList$Gene &
    #                     locDat$Group == "Down")] <- "LOMGlow_common_down"
  #}
  
  locDat$Group <- factor(locDat$Group, levels = c("Up", "Down", "None"))
  print(max(-log10(locDat$PValue)))
  print(max(locDat$logFC))
  print(min(locDat$logFC))
  p <- ggplot(locDat, aes(x = logFC, y = -log10(PValue), colour = Group)) +
    geom_point(size = 4) + scale_color_manual(values = c("#901E76", "#EDB953",
                                                 "grey"), drop = FALSE) +
    theme_bw() +
    scale_x_continuous(limits=c(-14, 14), expand = c(0, 0)) +
    scale_y_continuous(limits=c(-0, 20), expand = c(0, 0))
  p
  ggsave(paste0("Results/Graphics/Stockholm/Volcano/", i, ".pdf"))
  p <- p+theme(legend.position = "None")
  ggsave(paste0("Results/Graphics/Stockholm/Volcano/", i, "_no_legend_.pdf"),
         width = 5, height = 5)
  
}