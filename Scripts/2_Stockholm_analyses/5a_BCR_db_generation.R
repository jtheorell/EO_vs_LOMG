library(shazam)
library(alkazam)
library(SingleCellExperiment)
library(BiocParallel)
#We start by importing and combining all files. 
impList <- list.files("../External/Stockholm/Data/Immcantation_output_BCR", full.names = TRUE)
impList <- unname(sapply(impList, paste0, "/filtered_contig_db-pass.tsv"))
allBCRList <- lapply(impList, read.delim)

#Now, we add a correct cell name to each of them, before integrating. There are 13
#cells with the same barcode from different experiments, so this is an important step. 
allBCRList <- lapply(seq_along(allBCRList), function(x){
    locNum <- gsub("../External/Stockholm/Data/Immcantation_output_BCR/VDJ_21_0|_ENRB/filtered_contig_db-pass.tsv", 
                   "", impList[[x]])
    locDat <- allBCRList[[x]]
    locDat$cell_id <- paste0(locNum, "_", locDat$cell_id)
    locDat
})
fullBCR <- do.call("rbind", allBCRList)

uniqueBCRCells <- unique(fullBCR$cell_id)

bcrH <- fullBCR[which(fullBCR$locus == "IGH"),]

#At this stage it is of course pertinent to investigate the cells with two heavy chains. 
bcrDupH <- bcrH[c(which(duplicated(bcrH$cell_id)), which(duplicated(bcrH$cell_id, fromLast = TRUE))),]
bcrDupHOrd <- bcrDupH[order(bcrDupH$cell_id),]

#When looking at UMI count, junction length and sequence, there simply seems to
#be too many real differences for this to be anything else than a group of doublets.
#They are however kept in for now, but we add a new vector to the dataframe with 
#the information about these, so that they can be excluded downstream in the SCE.
bcrH$duplicated <- FALSE
bcrH$duplicated[c(which(duplicated(bcrH$cell_id)), which(duplicated(bcrH$cell_id, fromLast = TRUE)))] <- TRUE
fullBCR$duplicated <- FALSE
fullBCR$duplicated[which(fullBCR$cell_id %in% bcrDupH$cell_id)] <- TRUE

#Now, we cluster the data with a conservative, Hamming-based method, which 
#has been compared to other methods and been identified as giving identical results.
dist_ham <- distToNearest(bcrH, sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=8)

#Among the threshold selection methods, it is the simpler density method that
#works best for this data, considering both dist methods. 
output <- findThreshold(dist_ham$dist_nearest, method="density")
#This method does for clear reasons not find a threshold. But we will introduce one.
output@threshold <- 0.08
plot(output, binwidth=0.02, title="Density Method")

#This shows that the model does not set a reasonable threshold for clonal definitions. 
#So I will set a manual one to 0.08
changeoClusteringFile <- bcrH[,c("cell_id", "sequence_id", "v_call", "j_call",
                                        "junction")]
dir.create("../External/Stockholm/Data/BCR_dbs")
write.table(changeoClusteringFile, "../External/Stockholm/Data/BCR_dbs/1_changeoClusteringFile.tsv", 
            sep = "\t",row.names = FALSE)

#THis is run in the terminal
#cd /Users/jakob.theorell/Labbet/2023/230824_MG_10X_analysis_3_with_scvi-tools/External/Data/BCR_dbs
#DefineClones.py -d 1_changeoClusteringFile.tsv -o 2_changeoClusteringFile_Ham_result.tsv --act set --model ham --norm len --format changeo --sf JUNCTION --vf V_CALL --jf J_CALL --dist 0.08


#Now, we re-import these. 
hamClones <- read.table("../External/Stockholm/Data/BCR_dbs/2_changeoClusteringFile_Ham_result.tsv", header = TRUE)

#In two cases, there are cells that are classified as beloning to two clones. These
#are surely among the duplicated cells. Therefore, we add the extra requirement of only
#taking the first clone that is associated to the cell in question. 
fullBCR$HamClone <- unlist(lapply(fullBCR$cell_id, function(x){
    locVal <- hamClones$CLONE[which(hamClones$CELL_ID == x)][1]
}))

bcrH <- fullBCR[which(fullBCR$locus == "IGH"),]

#Now, we unname all clones that consist of one cell only. 
cloneSizeTab <- table(bcrH$HamClone)
realClones <- as.numeric(names(cloneSizeTab)[which(cloneSizeTab > 1)])
cloneSizeTab[which(cloneSizeTab > 1)]

fullBCR$HamClone[-which(fullBCR$HamClone %in% realClones)] <- NA

#Over to checking the intraclonal distance with Levenshtein distances. 
fullBCR$intraClonalDistance <- NA
cloneNames <- unique(fullBCR$HamClone)
cloneNames <- cloneNames[-which(is.na(cloneNames))]
for(i in cloneNames){
    locRows <- which(fullBCR$HamClone == i & fullBCR$locus == "IGH")
    locNames <- fullBCR$cell_id[locRows]
    uniqueLocNames <- unique(locNames)
    locCloneDist <- adist(fullBCR[locRows,"junction"])
    #To ensure that we do not get zero everywhere, we remove the diagonal
    diag(locCloneDist) <- 1000
    colnames(locCloneDist) <- locNames
    if(length(uniqueLocNames) > 1){
        for(j in uniqueLocNames){
            if(length(which(locNames == j)) > 1){
                innerCloneDist <- locCloneDist[,which(locNames == j)]
                minVal <- min(innerCloneDist[,which.min(apply(innerCloneDist,2,min))])
            } else {
                minVal <- min(locCloneDist[,j])
            }
            fullBCR$intraClonalDistance[which(fullBCR$cell_id == j)] <- 
                minVal
        }
    } else {
        NA
    }
    
}

table(fullBCR$intraClonalDistance)

#   0    1    2    3    4 
#1314   12   10   14   11 
(11+14+10+12)/(1314+11+14+10+12)
#0.03453343
#Which shows that this is a conservative but sensible threshold. 
#And the same is done for the light chain
fullBCR$intraClonalDistanceLight <- NA
cloneNames <- unique(fullBCR$HamClone)
cloneNames <- cloneNames[-which(is.na(cloneNames))]
for(i in cloneNames){
    locRows <- which(fullBCR$HamClone == i & fullBCR$locus != "IGH")
    locNames <- fullBCR$cell_id[locRows]
    uniqueLocNames <- unique(locNames)
    locCloneDist <- adist(fullBCR[locRows,"junction"])
    #To ensure that we do not get zero everywhere, we remove the diagonal
    diag(locCloneDist) <- 1000
    colnames(locCloneDist) <- locNames
    if(length(uniqueLocNames) > 1){
        for(j in uniqueLocNames){
            if(length(which(locNames == j)) > 1){
                innerCloneDist <- locCloneDist[,which(locNames == j)]
                minVal <- min(innerCloneDist[,which.min(apply(innerCloneDist,2,min))])
            } else {
                minVal <- min(locCloneDist[,j])
            }
            fullBCR$intraClonalDistanceLight[which(fullBCR$cell_id == j)] <- 
                minVal
        }
    } else {
        NA
    }
    
}

table(fullBCR$intraClonalDistanceLight)
#   0    1    2    6    8   14 
#1293    6   18    8    2    4
(14+8+6+2+1)/(1293+14+8+6+2+1)
#0.0234139
#Here there is a bit more laxity which is not strange considering that the clones
#have been defined on theheavy chain only, but it also shows that for the vast
#Majority of the cells, the clonal definitions make sense. Arguably, even a 
#Levenshtein distance of 14 is possible to accept as a true intraclonal distance. 

#Now, we add light chain information
fullBCR$light_type <- NA
for(i in unique(fullBCR$cell_id)){
    locRows <- which(fullBCR$cell_id == i)
    locLoc <- fullBCR$locus[locRows]
    lightLoc <- unique(locLoc[which(locLoc != "IGH")])
    if(length(lightLoc) == 2){
        fullBCR$light_type[locRows] <-  "double"
    } else if(length(lightLoc) == 1){
        fullBCR$light_type[locRows] <-  lightLoc
    } else {
        fullBCR$light_type[locRows] <-  "none"
    }
}

table(fullBCR$light_type)
#double    IGK    IGL   none 
#   174   4648   2728     33

#Now, we add mutational information
fullBCR <- observedMutations(fullBCR, sequenceColumn = "sequence_alignment",
                                   germlineColumn = "germline_alignment",
                                   regionDefinition = IMGT_V)

#We also add a column with the total number of mutations
fullBCR$All_mutations <- rowSums(fullBCR[,c("mu_count_cdr_r",
                                                        "mu_count_cdr_s",
                                                        "mu_count_fwr_r",
                                                        "mu_count_fwr_s"),])

fullBCR$Non_silent_mutations <- rowSums(fullBCR[,c("mu_count_cdr_r",
                                                               "mu_count_fwr_r"),])

#Now, we make some further mutational counts that will be included in the 
#SCE. 
fullBCR$All_mutations_H <- unlist(bplapply(seq_len(nrow(fullBCR)), function(x){
    locDat <- fullBCR[which(fullBCR$cell_id == fullBCR$cell_id[x]),]
    if(any(locDat$locus == "IGH")){
        #To make sure that we pick the most expressed one in case there 
        #are multiple, we here focus on the top expressor
        hMut <- locDat[which(locDat$locus == "IGH"),]
        hMut$All_mutations[which.max(hMut$umi_count)]
    } else {
        NA
    }
}))

fullBCR$All_mutations_L <- unlist(bplapply(seq_len(nrow(fullBCR)), function(x){
    locDat <- fullBCR[which(fullBCR$cell_id == fullBCR$cell_id[x]),]
    if(any(locDat$locus %in% c("IGK", "IGL"))){
        #To make sure that we pick the most expressed one in case there 
        #are multiple, we here focus on the top expressor
        lMut <- locDat[which(locDat$locus %in% c("IGK", "IGL")),]
        lMut$All_mutations[which.max(hMut$umi_count)]
    } else {
        NA
    }
}))

fullBCR$All_mutations_HL <- fullBCR$All_mutations_H+fullBCR$All_mutations_L

#And a simple column denoting if the cell is clonal or not. 
fullBCR$Clonal <- FALSE
fullBCR$Clonal[-which(is.na(fullBCR$HamClone))] <- TRUE

#Now, the V-family
fullBCR$vFam <- substr(fullBCR$v_call, 1, 5)

fullBCR$vFamH <- sapply(seq_len(nrow(fullBCR)), function(x){
    locDat <- fullBCR[which(fullBCR$cell_id == fullBCR$cell_id[x]),]
    #Here, we grab the dominant heavy chain if there are more than one. 
    if(any(grepl("H", locDat$vFam))){
        innerLocDat <- locDat[grep("H", locDat$vFam),]
        innerLocDat$vFam[which.max(innerLocDat$umi_count)]
    } else {
        NA
    }
})

#And here, we add a column that shows if the cell has a complete, productive
#BCR
fullBCR$complete <- unlist(bplapply(seq_len(nrow(fullBCR)), function(x){
    locDat <- fullBCR[which(fullBCR$cell_id == fullBCR$cell_id[x]),]
    if(any(locDat$locus == "IGH") &&
       (any(locDat$locus == "IGK") || 
        any(locDat$locus == "IGL"))){
        TRUE
    } else {
        FALSE
    }
}))


#We of course add back all the 
write.csv(fullBCR, "../External/Stockholm/Data/BCR_dbs/3_clonality_and_mutations_added.csv", 
          row.names = FALSE)
