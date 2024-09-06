library(pheatmap)
library(viridis)
library(robustbase)

logSCE <- readRDS("../External/Stockholm/Data/SCE_versions/6_with_Ox_neigh_clust.rds")
signRNA <- read.csv("Results/Data/Stockholm/Top_transcripts_EdgeR.csv", row.names = 1)

edgerOutcomes <- readRDS("Results/Data/Stockholm/EdgeR_significant_transcriptomes.rds")
#In this analysis, we focus on the genes that are co-regulated by all the 
#selected cell subsets in each of the EOMGlow and LOMGlow groups. To make these
#analyses meaningful, we will further split them into over- and underexpressed genes. 

EOMGlowLong <- Reduce(intersect, lapply(edgerOutcomes[grep("EOMG", names(edgerOutcomes))], row.names))
LOMGlowLong <- Reduce(intersect, lapply(edgerOutcomes[grep("LOMG", names(edgerOutcomes))], row.names))

#These now need to be further filtered, to make sure that all populations have
#the genes regulated in the same direction. 
EOMGlowDf <- do.call("rbind", lapply(EOMGlowLong, function(x){
    locVec <- unlist(lapply(edgerOutcomes[grep("EOMG", names(edgerOutcomes))], function(y){
        y$logFC[which(row.names(y) == x)]
    }))
    if(all(locVec > 0)){
        c(x, "high")
    } else if(all(locVec < 0)){
        c(x, "low")
    } else {
        NULL
    }
}))
colnames(EOMGlowDf) <- c("Gene", "Expression_compared_to_ref")
write.csv(EOMGlowDf, "Results/Data/Stockholm/EOMGlow_gene_list.csv", row.names = FALSE)

#This reduced the number by 20%.

LOMGlowDf <- do.call("rbind", lapply(LOMGlowLong, function(x){
  locVec <- unlist(lapply(edgerOutcomes[grep("LOMG", names(edgerOutcomes))], function(y){
    y$logFC[which(row.names(y) == x)]
  }))
  if(all(locVec > 0)){
    c(x, "high")
  } else if(all(locVec < 0)){
    c(x, "low")
  } else {
    NULL
  }
}))
#Here, we have both high and low genes. 
colnames(LOMGlowDf) <- c("Gene", "Expression_compared_to_ref")
write.csv(LOMGlowDf[order(LOMGlowDf[,2]),], "Results/Data/Stockholm/LOMGlow_gene_list.csv", row.names = FALSE)

#Here it reduced the number from 247 to 176. 

#Now, we select the populations that should be shown. 
focPops <- colnames(signRNA)
popList <- list("Mothers" = unique(gsub("|_.+_.", "\\1", focPops)),
                "Daughters" = focPops)

#Now, we calculate the median value for each gene in the dataset for each
#of these populations. 
focSCE <- logSCE[which(rowData(logSCE)$hgnc_symbol %in% c(EOMGlowLong, LOMGlowLong)),]
row.names(focSCE) <- rowData(focSCE)$hgnc_symbol

#We first identify the max value for each of the selected proteins. 
maxVals <- rowMax(as.matrix(logcounts(focSCE)))

meanMat <- do.call("cbind", lapply(names(popList), function(x){
    locPops <- popList[[x]]
    locRes <- do.call("cbind", lapply(locPops, function(y){
        if(x == "Mothers"){
            rowMeans(as.matrix(logcounts(focSCE)[,which(focSCE$cellType == y)]))
        } else {
            rowMeans(as.matrix(logcounts(focSCE)[,which(focSCE$oxNeighClust == y)]))
        }
    }))
    colnames(locRes) <- locPops
    locRes
}))

#Here, we for visualisation purposes set each row to the same scale. 
meanMatNorm <- t(apply(meanMat, 1, function(x){
    x/max(x)
}))

#And now, at the end, we create two different datasets: one for the EOMGlow 
#and one for the LOMGlow, both just showing the genes that are co-regulated in 
#all the cell subsets. 

meanMatEomgLow <- meanMatNorm[which(row.names(meanMatNorm) %in% EOMGlowLong),
                              colnames(meanMatNorm)[c(grep("CD4T", colnames(meanMatNorm)),
                                                      grep("NK", colnames(meanMatNorm)))]]

meanMatLomgLow <- meanMatNorm[which(row.names(meanMatNorm) %in% LOMGlowLong), 
                              colnames(meanMatNorm)[c(grep("B", colnames(meanMatNorm)),
                                                      grep("CD8T", colnames(meanMatNorm)))]]

dir.create("Results/Graphics/Stockholm")
pdf("Results/Graphics/Stockholm/RNA_heatmap_EOMGlow.pdf")
pheatmap(meanMatEomgLow, cluster_rows = TRUE, cluster_cols = FALSE, 
         border_color = NA, show_rownames = FALSE, treeheight_row = 0,
         color = viridis(20))
dev.off()

pdf("Results/Graphics/Stockholm/RNA_heatmap_LOMGlow.pdf")
pheatmap(meanMatLomgLow, cluster_rows = TRUE, cluster_cols = FALSE, 
         border_color = NA, show_rownames = FALSE, treeheight_row = 0,
         color = viridis(20))
dev.off()

