library(shazam)
library(alakazam)
library(SingleCellExperiment)
library(BiocParallel)
#We start by importing and combining all files. 
impList <- list.files("../External/Stockholm/Data/Immcantation_output_TCR", full.names = TRUE, pattern = "VDJ")
impList <- unname(sapply(impList, paste0, "/filtered_contig_db-pass.tsv"))
allTCRList <- lapply(impList, read.delim)

#Now, we add a correct cell name to each of them, before integrating. There are a number
#of cells that share bar code, so this is necessary. 
allTCRList <- lapply(seq_along(allTCRList), function(x){
    locNum <- gsub("../External/Stockholm/Data/Immcantation_output_TCR/VDJ_21_0|_ENRT/filtered_contig_db-pass.tsv", 
                   "", impList[[x]])
    locDat <- allTCRList[[x]]
    locDat$cell_id <- paste0(locNum, "_", locDat$cell_id)
    locDat
})
fullTCR <- do.call("rbind", allTCRList)

uniqueTCRCells <- unique(fullTCR$cell_id)
#This then nicely gives us 33544 TCR containing cells. 

#Now, we are going to look into clonality among these remaining cells
#Here, we use the Beta chain and a very low cutoff value, to allow for some
#errors in the reading, in accordance with the bottom lines at
#https://changeo.readthedocs.io/en/stable/examples/cloning.html
TCR_B <- fullTCR[which(fullTCR$locus == "TRB"),]

changeoClusteringFile <- TCR_B[,c("cell_id", "sequence_id", "v_call", "j_call",
                                 "junction")]
dir.create("../External/Stockholm/Data/TCR_dbs")
write.table(changeoClusteringFile, "../External/Stockholm/Data/TCR_dbs/1_changeoClusteringFile.tsv", 
            sep = "\t",row.names = FALSE)

#THis is run in the terminal
#cd /Users/jakob.theorell/Labbet/2023/230824_MG_10X_analysis_3_with_scvi-tools/External/Data/TCR_dbs
#DefineClones.py -d 1_changeoClusteringFile.tsv -o 2_changeoClusteringFile_Ham_result.tsv --act set --model ham --norm len --format changeo --sf JUNCTION --vf V_CALL --jf J_CALL --dist 0.03

#Now, we re-import these. 
hamClones <- read.table("../External/Stockholm/Data/TCR_dbs/2_changeoClusteringFile_Ham_result.tsv", header = TRUE)
library(BiocParallel)
fullTCR$HamClone <- unlist(bplapply(fullTCR$cell_id, function(x){
    locVal <- hamClones$CLONE[which(hamClones$CELL_ID == x)][1]
}))

tcrB <- fullTCR[which(fullTCR$locus == "TRB"),]

#Now, we unname all clones that consist of one cell only. 
cloneSizeTab <- table(tcrB$HamClone)
realClones <- as.numeric(names(cloneSizeTab)[which(cloneSizeTab > 1)])

fullTCR$HamClone[-which(fullTCR$HamClone %in% realClones)] <- NA

#Over to checking the intraclonal distance with Levenshtein distances. 
fullTCR$intraClonalDistance <- NA
cloneNames <- unique(fullTCR$HamClone)
cloneNames <- cloneNames[-which(is.na(cloneNames))]
for(i in cloneNames){
    locRows <- which(fullTCR$HamClone == i & fullTCR$locus == "TRB")
    locNames <- fullTCR$cell_id[locRows]
    uniqueLocNames <- unique(locNames)
    locCloneDist <- adist(fullTCR[locRows,"junction"])
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
            fullTCR$intraClonalDistance[which(fullTCR$cell_id == j)] <- 
                minVal
        }
    } else {
        NA
    }
    
}

table(fullTCR$intraClonalDistance)


#    0     1 
#13664    90
#This seems like an entirely reasonable result: a few cells are accepted as clonal
#members despite a very minor difference in their junction sequece, as these
#can be assumed to be related to sequencing errors. 

#Over to checking the intraclonal distance in the TCRA
fullTCR$intraClonalDistanceTCRA <- NA

for(i in cloneNames){
    locRows <- which(fullTCR$HamClone == i & fullTCR$locus == "TRA")
    locNames <- fullTCR$cell_id[locRows]
    uniqueLocNames <- unique(locNames)
    locCloneDist <- adist(fullTCR[locRows,"junction"])
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
            fullTCR$intraClonalDistanceTCRA[which(fullTCR$cell_id == j)] <- 
                minVal
        }
    } else {
        NA
    }
    
}

table(fullTCR$intraClonalDistanceTCRA)
#    0     1     2     3     5     8     9    10    11    12    13    14    15 
#10585    11    18     4     8    10     5     4     4    11    15    19    52 
#   16    17    18    19    20    21    22    23    24    25    26    29    30 
#   76    60    73    96    82    44    44    49    25     5     6     8     2 
#So this shows an interesting distribution, with a peak Levenshtein difference of
#19, apart from that 
length(which(fullTCR$intraClonalDistanceTCRA == 0))/length(fullTCR$intraClonalDistanceTCRA[-which(is.na(fullTCR$intraClonalDistanceTCRA))])
#or 94%, have a distance of 0. 
hist(fullTCR$intraClonalDistanceTCRA[which(fullTCR$intraClonalDistanceTCRA > 0)])
#So a small group seems to have changed the alpha chain after the initiation of the
#clonal expansion. 
#Here, a column of clonality is also added
fullTCR$Clonal <- FALSE
fullTCR$Clonal[-which(is.na(fullTCR$HamClone))] <- TRUE

#Now, we also check how many that have more than one functional alpha or beta chain
dual <- as.data.frame(do.call("rbind", bplapply(unique(fullTCR$cell_id), function(x){
    locDat <- fullTCR[which(fullTCR$cell_id == x),]
    if(nrow(locDat) > 1){
        locDatA <- locDat[which(locDat$locus == "TRA"),]
        locDatB <- locDat[which(locDat$locus == "TRB"),]
        if((length(which(locDatA$productive)) == 2) ||
           (length(which(locDatB$productive)) == 2)){
            c("cell_id" = x, "dual" = TRUE)
        } else {
            c("cell_id" = x, "dual" = FALSE)
        }
    } else {
        c("cell_id" = x, "dual" = FALSE)
    }
})))

table(dual)
#FALSE  TRUE 
#29789  3755

fullTCR$dual <- FALSE
fullTCR$dual[which(fullTCR$cell_id %in% 
                       dual$cell_id[which(dual$dual == "TRUE")])] <- TRUE
table(fullTCR$dual)
#FALSE  TRUE 
#53800 11870


multiplet <- as.data.frame(do.call("rbind", bplapply(unique(fullTCR$cell_id), function(x){
    locDat <- fullTCR[which(fullTCR$cell_id == x),]
    if(nrow(locDat) > 3){
        locDatA <- locDat[which(locDat$locus == "TRA"),]
        locDatB <- locDat[which(locDat$locus == "TRB"),]
        if((length(which(locDatA$productive)) > 2) ||
           (length(which(locDatB$productive)) > 2)){
            c("cell_id" = x, "multiplet" = TRUE)
        } else {
            c("cell_id" = x, "multiplet" = FALSE)
        }
    } else {
        c("cell_id" = x, "multiplet" = FALSE)
    }
})))

table(multiplet)
#FALSE 
#33544

#So, all the multiplets have already been excluded. Very good. 

#Now, the V-family
fullTCR$vFam <- substr(fullTCR$v_call, 1, 5)

fullTCR$vFamB <- unlist(bplapply(seq_len(nrow(fullTCR)), function(x){
    locDat <- fullTCR[which(fullTCR$cell_id == fullTCR$cell_id[x]),]
    #Here, we grab the dominant heavy chain if there are more than one. 
    if(any(grepl("B", locDat$vFam))){
        innerLocDat <- locDat[grep("B", locDat$vFam),]
        innerLocDat$vFam[which.max(innerLocDat$umi_count)]
    } else {
        NA
    }
}))

#And here, we add a column that shows if the cell has a complete, productive
#TCR
fullTCR$complete <- unlist(bplapply(seq_len(nrow(fullTCR)), function(x){
    locDat <- fullTCR[which(fullTCR$cell_id == fullTCR$cell_id[x]),]
    if(any(locDat$locus == "TRB") &&
       any(locDat$locus == "TRA")){
        TRUE
    } else {
        FALSE
    }
}))

#These are exported
write.csv(fullTCR, "../External/Stockholm/Data/TCR_DBs/3_clonal_TCR_db.csv", 
          row.names = FALSE)


