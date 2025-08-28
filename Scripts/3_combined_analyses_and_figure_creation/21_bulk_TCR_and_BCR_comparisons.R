library(SingleCellExperiment)

sceAllRaw <- 
  readRDS("../External/Stockholm/Data/SCE_versions/5_post_exclusions.rds")

colDatAllRaw <- colData(sceAllRaw)
#We start with the most straightforward comparisons of BCR and TCR clonality

lapply(c("BCR_Clonal", "TCR_Clonal"), function(x){
  airrClonTablePerDon <- as.data.frame.matrix(table(colDatAllRaw$sample, 
                                                    colDatAllRaw[,which(colnames(colDatAllRaw) == x)]))
  
  airrClonTablePerDonFrac <- airrClonTablePerDon/rowSums(airrClonTablePerDon)
  
  airrClonTablePerDonFrac$Group <- sapply(row.names(airrClonTablePerDonFrac), function(y){
    sceAllRaw$MG_type[which(sceAllRaw$sample == y)][1]
  })
  
  wilcox.test(airrClonTablePerDonFrac$`TRUE`[which(airrClonTablePerDonFrac$Group == "EOMG")],
              airrClonTablePerDonFrac$`TRUE`[which(airrClonTablePerDonFrac$Group == "LOMG")])$p.value
})

#[[1]]
#[1] 0.1418581
#
#[[2]]
#[1] 0.6620047

#So that is very clearly not different. 

#Now, the BCR mutations. 
airrClonTablePerDon <- as.data.frame.matrix(table(colDatAllRaw$sample, 
                                                  colDatAllRaw$BCR_All_mutations_H))

airrClonTablePerDonFrac <- airrClonTablePerDon/rowSums(airrClonTablePerDon)

airrClonTablePerDonFrac$Group <- sapply(row.names(airrClonTablePerDonFrac), function(y){
  sceAllRaw$MG_type[which(sceAllRaw$sample == y)][1]
})

#In this case, as this data is reasonably complex, we will reduce it to sums
#means for each group and then run a two-group Kolmogorov-Smirnov-test. I have
#done these analyses up and down, and I have never found anything. Here is another
#example of this:
eomgFracs <- rowMeans(airrClonTablePerDonFrac[which(airrClonTablePerDonFrac$Group == "EOMG"),
                                              1:ncol(airrClonTablePerDonFrac)-1])
lomgFracs <- colMeans(airrClonTablePerDonFrac[which(airrClonTablePerDonFrac$Group == "LOMG"),
                                              1:ncol(airrClonTablePerDonFrac)-1])

ks.test(eomgFracs, lomgFracs)
#p-value = 0.1184

#Now, finally: what about the V gene families?

lapply(c("BCR_vFamH", "TCR_vFamB"), function(x){
  airrClonTablePerDon <- as.data.frame.matrix(table(colDatAllRaw$sample, 
                                                    colDatAllRaw[,which(colnames(colDatAllRaw) == x)]))
  
  airrClonTablePerDonFrac <- airrClonTablePerDon/rowSums(airrClonTablePerDon)
  
  groupVec <- sapply(row.names(airrClonTablePerDonFrac), function(y){
    sceAllRaw$MG_type[which(sceAllRaw$sample == y)][1]
  })
  pVals <- sapply(colnames(airrClonTablePerDonFrac), function(y){
    wilcox.test(airrClonTablePerDonFrac[which(groupVec == "EOMG"),
                                        which(colnames(airrClonTablePerDonFrac) == y)],
                airrClonTablePerDonFrac[which(groupVec == "LOMG"),
                                        which(colnames(airrClonTablePerDonFrac) == y)],
                exact = FALSE)$p.value
  })
  p.adjust(pVals, method = "fdr")
})

#IGHV1     IGHV2     IGHV3     IGHV4     IGHV5     IGHV6     IGHV7 
#0.9485325 0.3678307 0.9485325 0.9485325 0.9485325 0.9485325 0.9485325 
#
#TRBV1     TRBV2     TRBV3     TRBV4     TRBV5     TRBV6     TRBV7     TRBV9 
#0.1353975 0.9673721 1.0000000 0.3254445 0.9673721 0.5867503 0.8980413 0.8980413 