#Here, we are going to look at some central transcripts to confirm that the
#populations we have identified actually have observed and not just imputed data
#confirming their phenotype. We will simply check if some of the most highly
#regulated proteins are also present and different on the transcript level
library(pheatmap)
library(viridis)
signPops <- suppressWarnings(read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1))


#Here, we also include the inhibtitory KIR genes (the activating ones are not present), to see if they are expressed or not, given the lack of NKG2C. 
transcriptList <- list("CD8T" = list(list("SignPop" = rownames(signPops)[grep("CD8T", rownames(signPops))][1],
                                          "Transcripts" = c("CCR7", "KLRG1", "CD27"),
                                          "Direction" = c("Up", "Down", "Up")),
                                     list("SignPop" = rownames(signPops)[grep("CD8T", rownames(signPops))][2],
                                          "Transcripts" = c("CD7", "CCR7", "KLRB1"),
                                          "Direction" = c("Down", "Down", "Up"))),
                       "NK" = list(list("SignPop" = rownames(signPops)[grep("NK", rownames(signPops))],
                                        "Transcripts" = c("B3GAT1", "CD2", "KLRC2"),
                                        "Direction" = c("Up", "Down", "Down")),
                                   list("SignPop" = rownames(signPops)[grep("NK", rownames(signPops))],
                                        "Transcripts" = rowData(completeRNA)$hgnc_symbol[grep("KIR2|KIR3", rowData(completeRNA)$hgnc_symbol)],
                                        "Direction" = rep("Unknown", length(grep("KIR2|KIR3", rowData(completeRNA)$hgnc_symbol))))))

completeRNA <- readRDS("../External/Stockholm/Data/SCE_versions/5_post_exclusions.rds")
#Starting with the 
freqExpressedList <- list()
for(a in 1:length(names(transcriptList))){
  i <- names(transcriptList)[a]
    print(i)
    protDat <- readRDS(paste0("../External/Oxford_and_Stockholm/Harmonisation/", i, "_Stockholm_data_with_Ox_neighbors.rds"))
    locList <- transcriptList[[i]]
    for(j in 1:length(locList)){
      transcriptDat <- completeRNA[which(rowData(completeRNA)$hgnc_symbol %in% locList[[j]]$Transcripts),
                                   which(completeRNA$cellType == i)]
      row.names(transcriptDat) <- rowData(transcriptDat)$hgnc_symbol
      if(identical(rownames(protDat), colnames(transcriptDat))){
        freqExpressed <- sapply(locList[[j]]$Transcripts, function(z){
          x <- counts(transcriptDat)[which(row.names(transcriptDat) == z),]
          hitDat <- x[which(protDat$Ox_hit_clust == locList[[j]]$SignPop)]
          nonHitDat <- x[-which(protDat$Ox_hit_clust == locList[[j]]$SignPop)]
          hitFreq <- length(which(hitDat > 0))/length(hitDat)
          nonHitFreq <- length(which(nonHitDat > 0))/length(nonHitDat)
          c(nonHitFreq, hitFreq, locList[[j]]$Direction[which(locList[[j]]$Transcripts == z)])
        })
        colnames(freqExpressed) <- locList[[j]]$Transcripts
        rownames(freqExpressed) <- c("Ref", "Hit", "Expected_direction_compared_to_ref")
        freqExpressedList[[(a*(a-1))+j]] <- freqExpressed
      } else {
        stop("Cells are in disarray. Order them better between the files!")
      }
    }
    
}
names(freqExpressedList) <- c(row.names(signPops), "NK_EOMGlow_1_KIR")

#We now remove the last row and turn the matrices numeric and then use that data for
#a few 

dir.create("Results/Graphics/Stockholm/Transcript_verifications")
lapply(names(freqExpressedList), function(x){
  locDatComplete <- freqExpressedList[[x]]
  locDatNumeric <- apply(locDatComplete[1:2,], 2, as.numeric)
  pdf(paste0("Results/Graphics/Stockholm/Transcript_verifications/", x, "_transcript_heatmap.pdf"))
  pheatmap(locDatNumeric, color = viridis(100), cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0,1, length = 101))
  dev.off()
})