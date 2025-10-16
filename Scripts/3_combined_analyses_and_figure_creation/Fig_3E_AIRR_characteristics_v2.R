library(SingleCellExperiment)
library(ggplot2)
library(viridis)
#Here, we investigate the BCR and TCR expression differences between
#the significant and non-significant populations. 
logSCE <- readRDS("../External/Stockholm/Data/SCE_versions/6_with_Ox_neigh_clust.rds")
significantClusters <- read.csv("Results/Data/Oxford_and_Stockholm/harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1)

logSCE$signOxClust <- logSCE$oxNeighClust
logSCE$signOxClust[-which(logSCE$oxNeighClust %in% row.names(significantClusters))] <- "None"

#Now, CD8
cd8Sce <- logSCE[,which(logSCE$TCR & logSCE$cellType == "CD8T")]


barGraphDat <- colData(cd8Sce)[,c("signOxClust","TCR_Clonal", "TCR_vFamB")]

#We do need to make a few changes here. 
barGraphDat$signOxClust <- factor(barGraphDat$signOxClust, levels = c("None","CD8T_LOMGlow_1", "CD8T_LOMGlow_2"))
barGraphDat$CD8T_Clonal <- factor(barGraphDat$TCR_Clonal, c(TRUE, FALSE))
barGraphDat$CD8T_V_gene_fam <- factor(barGraphDat$TCR_vFamB, levels = c("TRBV1", "TRBV2", "TRBV3", "TRBV4", "TRBV5", "TRBV6", "TRBV7", "TRBV9"))

colsAndNames <- list(
    list("CD8T_Clonal", c("grey", "black")),
    list("CD8T_V_gene_fam", cividis(8))
)

for(i in colsAndNames){
    print(i[[1]])
    locCounts <- table(barGraphDat[,i[[1]]], barGraphDat$signOxClust)
    #Here, we are seeing all the cells, including the clusters of interest, as the reference, 
    #as especially the naive cluster otherwise takes out such a chunk of the data. 
    locCounts[,1] <- locCounts[,1]+rowSums(locCounts[,2:3])
    locDat <- as.data.frame(locCounts)
    colnames(locDat) <- c(i[[1]], "Cluster", "value")
    p <- ggplot(locDat,                         
                aes(x = Cluster,
                    y = value,
                    fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        scale_fill_manual(values=i[[2]], drop = FALSE)
    p
    ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/", i[[1]], "_no_legend.pdf"),
           width = 9, height = 6)
    
    #And the data is saved
    write.csv(locDat, paste0("Results/Data/For_figure_file/Figure_3E_", i[[1]], ".csv"))
    #The stats are also included here. 
    print(round(fisher.test(locCounts[,1:2], simulate.p.value = TRUE)$p.value, 8))
    print(round(fisher.test(locCounts[,c(1,3)], simulate.p.value = TRUE)$p.value, 8))
}

#[1] "CD8T_Clonal"
#Saving 8.69 x 6.68 in image
#[1] 0
#[1] 0.3468512
#[1] "CD8T_V_gene_fam"
#Saving 8.69 x 6.68 in image
#[1] 0.015992
#[1] 0.00049975

#So this makes sense: there is no difference between the overall population and
#the LOMG_2 cluster in terms of clonality, but all other categories differ. 

#So we need to investigate the CD8T v-family finding further, checking which V families that are contributing
#to this difference, especially for the LOMG_1 population, where there is no apparrent difference. 
CD8VFamDistrP_LOMG1 <- unlist(lapply(row.names(locCounts), function(x){
    focRow <- locCounts[which(row.names(locCounts) == x),1:2]
    otherRows <- colSums(locCounts[-which(row.names(locCounts) == x),1:2])
    fisher.test(rbind(focRow, otherRows))$p.value
}))

names(CD8VFamDistrP_LOMG1) <- row.names(locCounts)
round(p.adjust(CD8VFamDistrP_LOMG1, method = "fdr"), 5)
#  TRBV1   TRBV2   TRBV3   TRBV4   TRBV5   TRBV6   TRBV7   TRBV8   TRBV9 
#0.50541 0.50541 1.00000 1.00000 0.01882 0.50541 0.10345 1.00000 0.50541 

#So here it is only TRBV5 that differs slightly. 

CD8VFamDistrP_LOMG2 <- unlist(lapply(row.names(locCounts), function(x){
  focRow <- locCounts[which(row.names(locCounts) == x),c(1,3)]
  otherRows <- colSums(locCounts[-which(row.names(locCounts) == x),c(1,3)])
  fisher.test(rbind(focRow, otherRows))$p.value
}))

names(CD8VFamDistrP_LOMG2) <- row.names(locCounts)
round(p.adjust(CD8VFamDistrP_LOMG2, method = "fdr"), 5)

#  TRBV1   TRBV2   TRBV3   TRBV4   TRBV5   TRBV6   TRBV7   TRBV8   TRBV9 
#0.00002 0.79149 0.57783 0.26909 0.79149 0.00000 0.00003 1.00000 0.02280 

#So this is really interesting. It is possible that it is the very significant
#aberrations in TCRVB6 and 7 that drive the rest, but it does not necessarily matter. 
#There are no TRBV8 genes here, hence the 1. All explained by the MAIT finding. 

