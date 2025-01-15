#Now, we are also going to test the size of two categories, namely all Tfh and all MAIT
library(SingleCellExperiment)
library(celldex)
library(SingleR)
tSinglerData <- MonacoImmuneData()

sceAllRaw <- readRDS("../External/Stockholm/Data/SCE_versions/5_post_exclusions.rds")

sceT <- sceAllRaw[,grep("T", sceAllRaw$cellType)]
sceABT <- sceT[,-grep("gd", sceT$cellType)]

commonTranscripts <- Reduce(intersect, list(row.names(tSinglerData), rowData(sceABT)$hgnc_symbol))
commonSinglerData <- tSinglerData[which(row.names(tSinglerData) %in% commonTranscripts),grep("T", tSinglerData$label.fine)]

abTSingler <- commonSinglerData[,-grep("gd", commonSinglerData$label.fine)]
sceABTCommon <- sceABT[which(rowData(sceABT)$hgnc_symbol %in% commonTranscripts),]
row.names(sceABTCommon) <- rowData(sceABTCommon)$hgnc_symbol

sceABTOrd <- sceABTCommon[match(row.names(commonSinglerData), row.names(sceABTCommon)),]
identical(row.names(sceABTOrd), row.names(commonSinglerData)) #TRUE of course. 

singlerRes <- SingleR(sceABTOrd, abTSingler , 
                    labels = abTSingler $label.fine)

#Now, what we are after is the statistics for the two populations.
tfhNumber <- table(as.factor(sceABTOrd$sample[which(singlerRes$labels == "Follicular helper T cells")]))
maitNumber <- table(as.factor(sceABTOrd$sample[which(singlerRes$labels == "MAIT cells")]))

allNumber <- table(sceAllRaw$sample)

all(identical(names(allNumber), names(tfhNumber)), identical(names(allNumber), names(maitNumber)))
#TRUE

tfhFreq <- tfhNumber/allNumber

maitFreq <- maitNumber/allNumber

#Now, we divide them per group. 

metaDat <- read.csv("Data/Stockholm/Metadata/Metadata_combined_used.csv")

groupVec <- sapply(names(tfhFreq), function(x){
  sceAllRaw$MG_type[which(sceAllRaw$sample == x)][1]
})

wilcox.test(tfhFreq[which(groupVec == "EOMG")], 
            tfhFreq[which(groupVec == "LOMG")])
#p-value = 0.9497

wilcox.test(maitFreq[which(groupVec == "EOMG")], 
            maitFreq[which(groupVec == "LOMG")])
#p-value = 0.11

#At this time, we also check the ones that have a full MAIT alpha TCR
realMaitCount <- table(sceABTOrd$sample[which(sceABTOrd$TCR_Japha33 & sceABTOrd$TCR_Vapha7.2)])
realMaitFreq <- realMaitCount/allNumber

wilcox.test(realMaitFreq[which(groupVec == "EOMG")], 
            realMaitFreq[which(groupVec == "LOMG")])

#And the p-value here is 0.345

#Now, to the most important control: what about T cells with a MAIT phenotype 
#and the two chains? 

pureMaitCount <- table(sceABTOrd$sample[which(sceABTOrd$TCR_Japha33 & sceABTOrd$TCR_Vapha7.2 &
                                                singlerRes$labels == "MAIT cells")])
pureMaitFreq <- pureMaitCount/allNumber

wilcox.test(pureMaitFreq[which(groupVec == "EOMG")], 
            pureMaitFreq[which(groupVec == "LOMG")])

#p-value = 0.2284
  