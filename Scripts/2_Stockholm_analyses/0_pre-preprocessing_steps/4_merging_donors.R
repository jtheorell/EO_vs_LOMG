#Here, we will first merge the three runs from each batch and then merge the 
#batches. All information about the original pool/batch will of course be kept
library(SingleCellExperiment)

sceList <- list.files("../External/Stockholm/Data/SCE", full.names = TRUE, pattern = "Run_")
#Here we divide the runs into the two batches.
sceListList <- list(sceList[1:3], sceList[4:6])

#Here, we individualize the cell names and the sample names for the two batches, 
#and add batch and runId number to each
fullSCE <- do.call("cbind", lapply(sceListList, function(x){
    do.call("cbind", lapply(x, function(y){
        locSce <- readRDS(y)
        locSce$runId <- as.numeric(gsub("../External/Stockholm/Data/SCE/Run_|\\.rds", "", y))
        colnames(locSce) <- paste0(locSce$runId[1], "_", colnames(locSce))
        locSce$batch <- if(locSce$runId[1] %in% c(19, 20, 21)){1}else{2}
        locSce$sample <- paste0(locSce$batch[1], "_", locSce$sample)
        locSce
    }))
}))

#Now, we just check the distributions of counts and mitochondrial percentages. 
plot(log(fullSCE$sum), fullSCE$subsets_Mito_percent,
     xlab="Total count", ylab='Mitochondrial %')
#It all makes sense.
saveRDS(fullSCE, "../External/Stockholm/Data/SCE/Raw_myasthenia.rds")

