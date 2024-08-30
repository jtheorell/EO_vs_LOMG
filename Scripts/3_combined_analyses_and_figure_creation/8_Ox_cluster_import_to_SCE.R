#Here, we reimport the data from the Oxford-Stockholm comparison and use that 
#to identify the top transcripts for each cluster. 

#Step one is to create a vector with all the clusters from the Oxford neighbor
#analysis. 
patSCE <- readRDS("../External/Stockholm/Data/SCE_versions/5_post_exclusions.rds")

#We will not include the ILC for two reasons: first, no such clusters were found
#to be significant and second, no such cells are defined in hte Stockholm data, 
#all being confined to the NK space. 
popVec <- c("B", "CD4T", "CD8T", "TCRgd", "NK")
patSCE$oxNeighClust <- NA
for(i in popVec){
    print(i)
    locDat <- readRDS(paste0("Results/Data/Oxford_and_Stockholm/Harmonisation/", i, "_Stockholm_data_with_Ox_neighbors.rds"))
    if(i == "TCRgd"){
        locRows <- which(patSCE$cellType == "gdT")
    } else {
        locRows <- which(patSCE$cellType == i)
    }
    
    if(identical(row.names(locDat), colnames(patSCE)[locRows])){
        patSCE$oxNeighClust[locRows] <- locDat$Ox_hit_clust
    } else {
        print("something is up, change this")
    }
}


table(patSCE$oxNeighClust, patSCE$cellType, useNA = "always")
#There are a few cells that are not defined and they are of course double positive
#and double negative T cells, as these are too few to be possible to analyse.

saveRDS(patSCE, "../External/Stockholm/Data/SCE_versions/6_with_Ox_neigh_clust.rds")

