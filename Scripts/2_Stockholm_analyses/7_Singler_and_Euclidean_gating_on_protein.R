#Now, we are creating clusters based on the protein expression among the shared 
#proteins in both datasets. They way we will go about this is to conventionally 
#gate out the main populations, but leaving populations in the middle, where it is
#uncerain if the cell is positive or not, and then we will use the euclidean distance
#to the centroids of each of the identified populations to associate each cell to
#its right cluster. We will identify CD3+CD4+, CD3+CD8+, CD3+DN, CD3+DP, NK, myeloid
#and B in this way. We will additionally use a singler classification to identify
#highly likely gamma-delta T cells and then integrate this information as well. 
library(SingleCellExperiment)
library(DepecheR)
library(FNN)
library(SingleR)
library(viridis)
library(BiocParallel)
fullSce <- readRDS("../External/Stockholm/Data/SCE_versions/3_including_AIRR.rds")

commonProts <- c("CD003_PROT", "CD004_PROT", "CD008_PROT", "CD011c_PROT", "CD016_PROT", 
                 "CD019_PROT", "CD056_PROT")

#We start here by defining a color vector function that will be used downstream
plasmaColors <- plasma(8)
colorSwitch <- function(cellTypeVector){
  sapply(cellTypeVector, switch, 
         "B" = plasmaColors[1], 
         "CD4T" = plasmaColors[6], 
         "CD8T" = plasmaColors[2],
         "DPT" = plasmaColors[8],
         "DNT" = plasmaColors[8],
         "gdT" = plasmaColors[5], 
         "ILC" = plasmaColors[3], 
         "NK" = plasmaColors[4],
         "Myeloid" = plasmaColors[7],
         "doublet" = "#333333",
         "none" = "grey")
}

protDat <- as.data.frame(t(logcounts(altExp(fullSce)))[,commonProts])
#We start with CD3 and CD19

dir.create("Diagnostics/Stockholm/Euclidean_gating")
png("Diagnostics/Stockholm/Euclidean_gating/1_CD3_vs_CD19.png", width = 300, height = 300)
par(mar = c(0,0,0,0), lwd = 3)
plot(protDat$CD003_PROT, protDat$CD019_PROT, 
     col = "#00000055", 
     xaxt = "n", yaxt = "n",
     ann=FALSE,
     pch = 19)
segments(x0 = -10, 
         y0 = 0.6, 
         x1 = 0.6, 
         y1 = 0.6,
         col = "red")
segments(x0 = 0.6, 
         y0 = 1, 
         x1 = 5, 
         y1 = 1,
         col = "red")
segments(x0 = 0.6, 
         y0 = 0.6, 
         x1 = 0.6, 
         y1 = 1,
         col = "red")
segments(x0 = 0.75, 
         y0 = 1, 
         x1 = 0.75, 
         y1 = 4,
         col = "red")
dev.off()

cellType <- rep("none", nrow(protDat))
cellType[which(protDat$CD003_PROT > 0.75 &
                 protDat$CD019_PROT > 1)] <- "doublet"
cellType[which(cellType == "none" & 
                 (protDat$CD003_PROT < 0.6 &
                    protDat$CD019_PROT > 0.6) |
                 (protDat$CD003_PROT < 0.75 &
                    protDat$CD019_PROT > 1))] <- "B"


pdf("Diagnostics/Stockholm/Euclidean_gating/1_CD3_vs_CD19_gated.pdf", width = 5, height = 5)
plot(protDat$CD003_PROT, protDat$CD019_PROT, col = colorSwitch(cellType))
segments(x0 = -10, 
         y0 = 0.6, 
         x1 = 0.6, 
         y1 = 0.6,
         col = "red")
segments(x0 = 0.6, 
         y0 = 1, 
         x1 = 5, 
         y1 = 1,
         col = "red")
segments(x0 = 0.6, 
         y0 = 0.6, 
         x1 = 0.6, 
         y1 = 1,
         col = "red")
segments(x0 = 0.75, 
         y0 = 1, 
         x1 = 0.75, 
         y1 = 4,
         col = "red")
dev.off()

#Now, we go on with the non-doublet, non-B. 
nonDoubNonB <- protDat[which(cellType == "none"),]

png("Diagnostics/Stockholm/Euclidean_gating/2_CD3_vs_CD56.png", width = 300, height = 300)
par(mar = c(0,0,0,0), lwd = 3)
plot(nonDoubNonB$CD003_PROT, nonDoubNonB$CD056, 
     col = "#00000055", 
     xaxt = "n", yaxt = "n",
     ann=FALSE,
     pch = 19)
segments(x0 = -5, 
         y0 = 0.3, 
         x1 = 0.4, 
         y1 = 0.3,
         col = "red")
segments(x0 = 0.4, 
         y0 = 0.3, 
         x1 = 0.4, 
         y1 = 2,
         col = "red")
dev.off()

cellType[which(cellType == "none" & 
                 (protDat$CD003_PROT < 0.4 &
                    protDat$CD056_PROT > 0.3))] <- "NK"
pdf("Diagnostics/Stockholm/Euclidean_gating/2_CD3_vs_CD56_gated.pdf", width = 5, height = 5)
plot(protDat$CD003_PROT, protDat$CD056, col = colorSwitch(cellType))
segments(x0 = -5, 
         y0 = 0.3, 
         x1 = 0.4, 
         y1 = 0.3,
         col = "red")
segments(x0 = 0.4, 
         y0 = 0.3, 
         x1 = 0.4, 
         y1 = 2,
         col = "red")
dev.off()

#Now, myeloid
nonDoubNonBNonNK <- protDat[which(cellType == "none"),]
png("Diagnostics/Stockholm/Euclidean_gating/3a_CD3_vs_CD11c.png", width = 300, height = 300)
par(mar = c(0,0,0,0), lwd = 3)
plot(nonDoubNonBNonNK$CD003_PROT, 
     nonDoubNonBNonNK$CD011c_PROT, 
     col = "#00000055", 
     xaxt = "n", yaxt = "n",
     ann=FALSE,
     pch = 19)

#Myeloid vs T
segments(x0 = -5, 
         y0 = 1.8, 
         x1 = 0.7, 
         y1 = 1.8,
         col = "red")
segments(x0 = 0.7, 
         y0 = 1.8, 
         x1 = 0.7, 
         y1 = 5,
         col = "red")
#T
segments(x0 = 1, 
         y0 = -1, 
         x1 = 1, 
         y1 = 1,
         col = "red")
segments(x0 = 1, 
         y0 = 1, 
         x1 = 3, 
         y1 = 1,
         col = "red")
dev.off()

png("Diagnostics/Stockholm/Euclidean_gating/3b_CD56_vs_CD11c.png", width = 300, height = 300)
par(mar = c(0,0,0,0), lwd = 3)
plot(nonDoubNonBNonNK$CD056_PROT, 
     nonDoubNonBNonNK$CD011c_PROT, 
     col = "#00000055", 
     xaxt = "n", yaxt = "n",
     ann=FALSE,
     pch = 19)

#Myeloid vs NK
segments(x0 = -0.5, 
         y0 = 1.8, 
         x1 = 0.15, 
         y1 = 1.8,
         col = "red")
segments(x0 = 0.15, 
         y0 = 1.8, 
         x1 = 0.15, 
         y1 = 5,
         col = "red")
dev.off()


cellType[which(cellType == "none" &
                 protDat$CD003_PROT > 1 &
                 protDat$CD011c_PROT < 1)] <- "T"
cellType[which(cellType == "none" &
                 protDat$CD003_PROT < 0.7 &
                 protDat$CD011c_PROT > 1.8 &
                 protDat$CD056_PROT < 0.15)] <- "Myeloid"

pdf("Diagnostics/Stockholm/Euclidean_gating/3_CD3_vs_CD11c_gated.pdf", width = 5, height = 5)
plot(protDat$CD003_PROT, protDat$CD011c, col = colorSwitch(cellType))
segments(x0 = -5, 
         y0 = 1.8, 
         x1 = 0.7, 
         y1 = 1.8,
         col = "red")
segments(x0 = 0.7, 
         y0 = 1.8, 
         x1 = 0.7, 
         y1 = 5,
         col = "red")
segments(x0 = 1, 
         y0 = -1, 
         x1 = 1, 
         y1 = 5,
         col = "red")
dev.off()

#Now, we pause and classify all the cells into these major clusters. These are
#only used during this though, and they will not be kept. 

#And now, we will define the nearest neighbor for each cell and thereby define
#its cluster. 
latDat <- reducedDim(fullSce, "totalVI_lat_rep")

#uniqueCellTypes <- unique(cellType)
#uniqueCellTypes <- uniqueCellTypes[-which(uniqueCellTypes == "none")]
#centroidTab <- do.call("rbind", lapply(uniqueCellTypes, function(x){
#    colMeans(latDat[which(cellType == x),])
#}))
#rownames(centroidTab) <- uniqueCellTypes

#And now, we use this to associate each cell to its closest centroid
noneCellRows <- which(cellType == "none")
latDatNone <- latDat[noneCellRows,]
latDatDefined <- latDat[-noneCellRows,]

nearestCellType <- knnx.index(latDatDefined, latDatNone, k = 1)[,1]
cellTypesOriginal <- cellType[-noneCellRows]
nearestCellTypeNamed <- nearestCellType
for(i in unique(cellTypesOriginal)){
    locRows <- which(cellTypesOriginal == i)
    nearestCellTypeNamed[which(nearestCellType %in% locRows)] <- i
}

cellTypeIntermediate <- cellType
cellTypeIntermediate[noneCellRows] <- nearestCellTypeNamed
table(cellTypeIntermediate)
#      B doublet Myeloid      NK       T 
#   8754     194   18200    8680   64069 

#At this stage, before going into the T cell compartment, we will use SingleR
#to sub-classify all the T cells as gamma-delta or non-gamma-delta
fullSce$cellType <- cellTypeIntermediate

monDat <- MonacoImmuneData()

#Now, these are synced. 
possibleGDCols <- which(fullSce$cellType == "T" &
                            (is.na(fullSce$TCR_productive 
                                   | fullSce$TCR_productive != TRUE)))
namedSce <- fullSce[,possibleGDCols]
row.names(namedSce) <- rowData(fullSce)$hgnc_symbol

monDatUsed <- monDat[which(row.names(monDat) %in% row.names(namedSce)),]
namedSceUsed <- namedSce[which(row.names(namedSce) %in% row.names(monDatUsed)),]
monDatOrd <- monDatUsed[match(row.names(namedSceUsed), row.names(monDatUsed)),]

identical(row.names(monDatOrd), row.names(namedSceUsed)) #TRUE

#This is the data we will use here: 
TMon <- monDatOrd[,grep("T cells", monDatOrd$label.main)]

locSceColVec <- ceiling(seq_len(ncol(namedSce))/1000)

#Now, we introduce an extra harsh criterion here, only accepting the cell type
#if the delta.next is below 0.1, which is the case for between 80 and 90% of the
#cells in this dataset. We also only care about the gamma-delta. 
monacoIdList <- bplapply(unique(locSceColVec), function(y){
    print(y)
    innerSce <- namedSce[,which(locSceColVec == y)]
    innerRes <- SingleR(innerSce, TMon, 
                        labels = TMon$label.fine, 
                        num.threads = 1)
    innerLabel <- innerRes$labels
    innerLabel[-grep("gd T", innerLabel)] <- "none"
    innerLabel[grep("gd T", innerLabel)] <- "gdT"
    #Here, we also define the highest correlation seen to any of the two cell
    #types, and we only allow for the ones over a correlation of 0.25 to be included
    dgDat <- as.data.frame(innerRes$scores[,grep("gd T", colnames(innerRes$scores))])
    dgDat$maxVal <- apply(dgDat, 1, max)
    innerLabel[which(dgDat$maxVal < 0.25)] <- "none"
    innerLabel[which(innerRes$delta.next > 0.1)] <- "none"
    data.frame("Cell_id" = innerSce$Cell_id, 
               "Monaco" = innerLabel)
})

monacoId <- do.call("rbind", monacoIdList)
monacoIdOrdered <- monacoId[match(namedSce$Cell_id, monacoId$Cell_id),]
identical(monacoIdOrdered$Cell_id, namedSce$Cell_id) #TRUE

#This is now transferred to the cell type vector
fullNoneTVector <- rep("none", length(cellTypeIntermediate))
fullNoneTVector[possibleGDCols] <- monacoIdOrdered$Monaco
nonNonePoss <- which(fullNoneTVector != "none")

cellType[nonNonePoss] <- fullNoneTVector[nonNonePoss]
cellTypeIntermediate[nonNonePoss] <- fullNoneTVector[nonNonePoss]
#And finally, we define the T-cell compartment for the remaining cells. 

TDat <- protDat[which(cellTypeIntermediate == "T"),]

png("Diagnostics/Stockholm/Euclidean_gating/4_CD4_vs_CD8_on_T.png", width = 300, height = 300)
par(mar = c(0,0,0,0), lwd = 3)
plot(TDat$CD004_PROT, TDat$CD008, 
     col = "#00000055", 
     xaxt = "n", yaxt = "n",
     ann=FALSE,
     pch = 19)
#DN
segments(x0 = -2, 
         y0 = 1, 
         x1 = 1, 
         y1 = 1,
         col = "red")
segments(x0 = 1, 
         y0 = -2, 
         x1 = 1, 
         y1 = 1,
         col = "red")
#CD8
segments(x0 = -2, 
         y0 = 2, 
         x1 = 1.5, 
         y1 = 2,
         col = "red")
segments(x0 = 1.5, 
         y0 = 2, 
         x1 = 1.5, 
         y1 = 6,
         col = "red")
#CD4
segments(x0 = 1.5, 
         y0 = 2.5, 
         x1 = 7, 
         y1 = 2.5,
         col = "red")
segments(x0 = 1.5, 
         y0 = -2, 
         x1 = 1.5, 
         y1 = 2.5,
         col = "red")
#DP
segments(x0 = 2, 
         y0 = 3, 
         x1 = 7, 
         y1 = 3,
         col = "red")
segments(x0 = 2, 
         y0 = 3, 
         x1 = 2, 
         y1 = 7,
         col = "red")
dev.off()

cellType[which(cellTypeIntermediate == "T" &
                 protDat$CD004_PROT < 1 &
                 protDat$CD008_PROT < 1)] <- "DNT"
cellType[which(cellTypeIntermediate == "T" &
                 protDat$CD004_PROT < 1.5 &
                 protDat$CD008_PROT > 2)] <- "CD8T"
cellType[which(cellTypeIntermediate == "T" &
                 protDat$CD004_PROT > 1.5 &
                 protDat$CD008_PROT < 2.5)] <- "CD4T"
cellType[which(cellTypeIntermediate == "T" &
                 protDat$CD004_PROT > 2 &
                 protDat$CD008_PROT > 3)] <- "DPT"

#And now, we project these defined cell types on the umap. 
plotCellType <- cellType
plotCellType[which(plotCellType == "T")] <- "none"
dColorPlot(plotCellType, xYData = reducedDim(fullSce, "UMAP"), colors = colorSwitch(plotCellType), 
           plotName = "Manually_defined_cells", plotDir = "Diagnostics/Stockholm/Euclidean_gating")

table(cellType)
#cellType
#      B    CD4T    CD8T     DNT doublet     DPT     gdT Myeloid      NK    none       T 
#   8648   41212   17340     515     189      80    3219   16862    5046    5510    1276 


#And now, we use this to associate each cell to its closest neighbor
noneCellRows <- which(cellType == "none")
latDatNone <- latDat[noneCellRows,]
latDatDefined <- latDat[-noneCellRows,]

nearestCellType <- knnx.index(latDatDefined, latDatNone, k = 1)[,1]
cellTypesOriginal <- cellType[-noneCellRows]
nearestCellTypeNamed <- nearestCellType
for(i in unique(cellTypesOriginal)){
    locRows <- which(cellTypesOriginal == i)
    nearestCellTypeNamed[which(nearestCellType %in% locRows)] <- i
}

cellType[noneCellRows] <- nearestCellTypeNamed
table(cellType)


#      B    CD4T    CD8T     DNT doublet     DPT     gdT Myeloid      NK       T 
#   8677   41608   17527     646     193      85    3360   17914    8560    1327 

#ANd now, we take one step further, and define the remaining unspecified T-cells as
#belonging to one of the five T-cell populations. We do this in two steps: first, 
#all the cells are defined as belonging to one of the subsets
otherTRows <- which(grepl("T", cellType) & cellType != "T")
latDatOtherT <- latDat[otherTRows,]
tCellRows <- which(cellType == "T")
latDatT <- latDat[tCellRows,]

nearestCellType <- knnx.index(latDatOtherT, latDatT, k = 1)[,1]
cellTypesOriginal <- cellType[otherTRows]
nearestCellTypeNamed <- nearestCellType
for(i in unique(cellTypesOriginal)){
    locRows <- which(cellTypesOriginal == i)
    nearestCellTypeNamed[which(nearestCellType %in% locRows)] <- i
}

cellType[tCellRows] <- nearestCellTypeNamed

#Then the alpha-beta T cells in the gamma-delta category are re-defined. 

aphaBetaRows <- which(cellType == "gdT" & fullSce$TCR_productive == TRUE)
length(aphaBetaRows) #58
latDatNonGgT <- latDat[aphaBetaRows,]

otherTRows <- grep("CD4|CD8|DN|DP", cellType)
latDatOtherT <- latDat[otherTRows,]

nearestCellType <- knnx.index(latDatOtherT, latDatNonGgT, k = 1)[,1]
cellTypesOriginal <- cellType[otherTRows]
nearestCellTypeNamed <- nearestCellType
for(i in unique(cellTypesOriginal)){
    locRows <- which(cellTypesOriginal == i)
    nearestCellTypeNamed[which(nearestCellType %in% locRows)] <- i
}


cellType[aphaBetaRows] <- nearestCellTypeNamed
table(cellType)
#      B    CD4T    CD8T     DNT doublet     DPT     gdT Myeloid      NK 
#   8677   41829   18229     741     193      93    3661   17914    8560
length(which(cellType == "gdT" & fullSce$TCR_productive)) 
#0

dColorPlot(cellType, xYData = reducedDim(fullSce, "UMAP"), 
           colors = colorSwitch(cellType),
           plotName = "Cells_defined_by_Euclidean_gating", 
           plotDir = "Diagnostics/Stockholm/Euclidean_gating")

dDensityPlot(xYData = reducedDim(fullSce, "UMAP"), idsVector = cellType,
           plotDir = "Diagnostics/Stockholm/Euclidean_gating/Cluster_density")

#This seems clean and meaningful. 
fullSce$cellType <- cellType
#Before saving, we also make one small correction: the cell name column is named X
saveRDS(fullSce, "../External/Stockholm/Data/SCE_versions/4_including_cell_type.rds")

#Later addition of new coldata to this file: 
#fullSce <- readRDS("../External/Stockholm/Data/SCE_versions/4_including_cell_type.rds")
#fullSceNewColData <- readRDS("../External/Stockholm/Data/SCE_versions/3_including_AIRR.rds")
#identical(colnames(fullSce), colnames(fullSceNewColData)) #TRUE
#identical(colData(fullSce), colData(fullSceNewColData)) #FALSE, as we had planned
#newCols <- colData(fullSceNewColData)[,-which(colnames(colData(fullSceNewColData)) %in% colnames(colData(fullSce)))]
#colData(fullSce) <- cbind(colData(fullSce), newCols)
#saveRDS(fullSce, "../External/Stockholm/Data/SCE_versions/4_including_cell_type.rds")





    