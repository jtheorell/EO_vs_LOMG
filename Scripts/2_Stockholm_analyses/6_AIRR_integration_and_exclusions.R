library(SingleCellExperiment)
library(BiocParallel)
fullSce <- readRDS("../External/Stockholm/Data/SCE_versions/2_normalised_SCE_with_full_protein.rds")

#We need to start by renaming the cell id column to cell id, from x. 
colnames(colData(fullSce))[which(colnames(colData(fullSce)) == "X")] <- "Cell_id"
bcrDb <- read.csv("../External/Stockholm/Data/BCR_dbs/3_clonality_and_mutations_added.csv")
tcrDb <- read.csv("../External/Stockholm/Data/TCR_DBs/3_clonal_TCR_db.csv")

#Now, we start by defining cells as having a TCR or a BCR. 
fullSce$BCR <- FALSE
fullSce$BCR[which(fullSce$Cell_id %in% bcrDb$cell_id)] <- TRUE
fullSce$TCR <- FALSE
fullSce$TCR[which(fullSce$Cell_id %in% tcrDb$cell_id)] <- TRUE

table(fullSce$BCR, fullSce$TCR)
#        FALSE  TRUE
#  FALSE 70779 26398
#  TRUE   2649    71

#Now, we start with the inclusions of new columns. 
#Common columns first. 
bcrCols <- c("cell_id", "productive",
             "complete", "HamClone", 
             "Clonal", "intraClonalDistance",
             "intraClonalDistanceLight",
             "light_type", "vFamH",
             "All_mutations_H", 
             "All_mutations_L", 
             "All_mutations_HL", 
             "duplicated", "c_call")
colDatBCR <- bcrDb[,which(colnames(bcrDb) 
                          %in% bcrCols)]
colnames(colDatBCR) <- paste0("BCR_", colnames(colDatBCR))
#And now, as all these columns are per cell and not per entry, we can collate
#this to the unique cells. 
bcrCells <- unique(colDatBCR$BCR_cell_id)
colDatBCRReduced <- do.call("rbind", bplapply(bcrCells, function(x){
    
  colDatBCR[which(colDatBCR$BCR_cell_id == x)[1],]
}))

#Here, we exclude the BCRs that do not have an associated transcriptome. 
colDatBCRPresent <- colDatBCRReduced[which(colDatBCRReduced$BCR_cell_id %in%
                                             fullSce$Cell_id),]
#This leads to the devastating loss of 
(3827-2720)/3827
#29% of the BCRs. But what can we do?

bcrDummyDf <- as.data.frame(
  matrix(NA, (ncol(fullSce)-nrow(colDatBCRPresent)), ncol(colDatBCRPresent)))
colnames(bcrDummyDf) <- colnames(colDatBCR)
bcrDummyDf$BCR_cell_id <- 
  fullSce$Cell_id[-which(fullSce$Cell_id %in% colDatBCRPresent$BCR_cell_id)]
fullBcrColDat <- rbind(colDatBCRPresent, bcrDummyDf)
fullBcrColDatOrdered <- fullBcrColDat[
  match(fullSce$Cell_id, fullBcrColDat$BCR_cell_id),]
identical(fullBcrColDatOrdered$BCR_cell_id, fullSce$Cell_id) #TRUE

#So now, over to the TCR data. 
tcrCols <- c("cell_id", "productive",
             "complete", "HamClone", 
             "Clonal", "intraClonalDistance",
             "intraClonalDistanceTCRA",
              "vFamB", "dual", "Vapha7.2", "Japha33", "Vbeta13", "Vbeta2")

colDatTCR <- tcrDb[,which(colnames(tcrDb) 
                          %in% tcrCols)]
colnames(colDatTCR) <- paste0("TCR_", colnames(colDatTCR))
#And now, as all these columns are per cell and not per entry, we can collate
#this to the unique cells. 
tcrCells <- unique(colDatTCR$TCR_cell_id)
colDatTCRReduced <- do.call("rbind", bplapply(tcrCells, function(x){
  colDatTCR[which(colDatTCR$TCR_cell_id == x)[1],]
}))

#Here, we exclude the BCRs that do not have an associated transcriptome. 
colDatTCRPresent <- colDatTCRReduced[which(colDatTCRReduced$TCR_cell_id %in%
                                             fullSce$Cell_id),]

#This leads to the devastating loss of  (33544-26469)/33544 or 
#21% of the TCRs, so slightly better anyway. 
tcrDummyDf <- as.data.frame(
  matrix(NA, (ncol(fullSce)-nrow(colDatTCRPresent)), ncol(colDatTCRPresent)))
colnames(tcrDummyDf) <- colnames(colDatTCR)
tcrDummyDf$TCR_cell_id <- 
  fullSce$Cell_id[-which(fullSce$Cell_id %in% colDatTCRPresent$TCR_cell_id)]
fullTcrColDat <- rbind(colDatTCRPresent, tcrDummyDf)
fullTcrColDatOrdered <- fullTcrColDat[
  match(fullSce$Cell_id, fullTcrColDat$TCR_cell_id),]
identical(fullTcrColDatOrdered$TCR_cell_id, fullSce$Cell_id) #TRUE

#And so, this data can be incorporated
colDatFullSce <- colData(fullSce)
colDatPlusAirr <- cbind(colDatFullSce, fullBcrColDatOrdered, fullTcrColDatOrdered)
colData(fullSce) <- colDatPlusAirr

#And there we are, integration-wise. 
saveRDS(fullSce, "../External/Stockholm/Data/SCE_versions/3_including_AIRR.rds")

