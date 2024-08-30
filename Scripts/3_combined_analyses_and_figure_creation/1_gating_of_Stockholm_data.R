library(SingleCellExperiment)
sthlmProtDat <- readRDS("../External/Stockholm/Data/SCE_versions/Stockholm_protein_data_including_colData.rds")

#Now, we will gate the data here in R, focusing on the populations identified 
#as potentially significant in the Oxford data. 

protDf <- as.data.frame(t(logcounts(sthlmProtDat)))

#The analysis needs to be performed by hand, on a population by population basis
pResLow <- read.csv("Results/Data/Oxford/Gating/Interesting_gated_populations_p_vals_Oxford.csv",
                    row.names = 1)

#We start with the B-cells, where we look for these populatsions: 
#B_IgDposCD27negCD24posCD38posTrans
#B_IgDnegCD27neg  
BDf <- protDf[which(sthlmProtDat$cellType == "B"),]
plot(BDf$IgD_PROT, BDf$CD027_PROT)
abline(v = 0.5)
abline(h = 0.6)
B_IgDnegCD27neg_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "B" &
                                                 protDf$IgD_PROT < 0.5 &
                                                 protDf$CD027_PROT < 0.6)]
                   

plot(BDf$CD024_PROT, BDf$CD038_PROT)
abline(v = 0.34)
abline(h = 2)
B_IgDposCD27negCD24posCD38posTrans_Counts <- 
  sthlmProtDat$sample[which(sthlmProtDat$cellType == "B" &
                                    protDf$IgD_PROT >= 0.5 &
                                    protDf$CD027_PROT < 0.6 &
                                    protDf$CD024_PROT > 0.35 &
                                    protDf$CD038_PROT > 2)]
                   
#Now CD4. 
#CD4T_naive_CD25                             
#CD4T_EM_CD161  
CD4Df <- protDf[which(sthlmProtDat$cellType == "CD4T"),]
#We annoyingly lack CCR7 in this panel, so we will need to improvise. We will mainly
#include CD62L as a proxy for CCR7 here. 
hist(CD4Df$CD045RA_PROT, breaks = 50)
abline(v = 0.75)
hist(CD4Df$CD045RO_PROT, breaks = 50)
abline(v = 0.5)
hist(CD4Df$CD062L_PROT, breaks = 50)
abline(v = 0.7)
hist(CD4Df$CD028_PROT, breaks = 50)
abline(v = 0.9)
hist(CD4Df$CD057_PROT, breaks = 50)
abline(v = 0.4)
hist(CD4Df$CD025_PROT, breaks = 50)
abline(v = 0.25)
hist(CD4Df$CD161_PROT, breaks = 50)
abline(v = 0.5)
CD4T_naive_CD25_Counts <-  sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD4T" &
                                                              protDf$CD045RA_PROT > 0.75 &
                                                              protDf$CD045RO_PROT < 0.5 &
                                                              protDf$CD028_PROT > 0.9 & 
                                                              protDf$CD057_PROT < 0.4 &
                                                              protDf$CD062L_PROT > 0.7 &
                                                              protDf$CD025_PROT > 0.5)]
                                              
CD4T_EM_CD161_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD4T" &
                                                      protDf$CD045RA_PROT < 0.75 &
                                                      protDf$CD045RO_PROT > 0.5 &
                                                      protDf$CD062L_PROT < 0.7 &
                                                      protDf$CD161_PROT > 0.5)]        
  
#And now CD8.
#We annoyingly lack CCR7 in this panel, so we will need to improvise. We will mainly
#include CD62L as a proxy for CCR7 here. 

#"CD8T_CM_CD25"                       "CD8T_CM_CD38"                      
#"CD8T_EM_CD25"                       "CD8T_EM_CD38"                      
#"CD8T_EM_CD161"                      "CD8T_naive"                        
#"CD8T_naive_CD7"                     "CD8T_naive_CD27"                   
#"CD8T_naive_CD31"                    "CD8T_naive_CD38"                   
#"CD8T_naive_CD161"
CD8Df <- protDf[which(sthlmProtDat$cellType == "CD8T"),]
hist(CD8Df$CD045RA_PROT, breaks = 50)
abline(v = 0.75)
hist(CD8Df$CD045RO_PROT, breaks = 50)
abline(v = 0.75)
hist(CD8Df$CD062L_PROT, breaks = 50)
abline(v = 0.75)
hist(CD8Df$CD028_PROT, breaks = 50)
abline(v = 0.4)
hist(CD8Df$CD057_PROT, breaks = 50)
abline(v = 0.4)

hist(protDf$CD007_PROT, breaks = 50)
abline(v = 1)
hist(protDf$CD027_PROT, breaks = 50)
abline(v = 0.7)
hist(CD8Df$CD031_PROT, breaks = 50)
abline(v = 0.75)
hist(CD8Df$CD038_PROT, breaks = 50)
abline(v = 1)
hist(CD8Df$CD161_PROT, breaks = 50)
abline(v = 0.7)
hist(CD8Df$CD025_PROT, breaks = 50)
abline(v = 0.2)


CD8T_naive_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                      protDf$CD045RA_PROT > 0.75 &
                                                      protDf$CD045RO_PROT < 0.75 &
                                                      protDf$CD028_PROT > 0.4 & 
                                                      protDf$CD057_PROT < 0.4 &
                                                      protDf$CD062L_PROT > 0.75)]

CD8T_naive_CD7_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                            protDf$CD045RA_PROT > 0.75 &
                                                            protDf$CD045RO_PROT < 0.75 &
                                                            protDf$CD028_PROT > 0.4 & 
                                                            protDf$CD057_PROT < 0.4 &
                                                            protDf$CD062L_PROT > 0.75 &
                                                            protDf$CD007_PROT > 1)]

CD8T_naive_CD27_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                            protDf$CD045RA_PROT > 0.75 &
                                                            protDf$CD045RO_PROT < 0.75 &
                                                            protDf$CD028_PROT > 0.4 & 
                                                            protDf$CD057_PROT < 0.4 &
                                                            protDf$CD062L_PROT > 0.75 &
                                                            protDf$CD027_PROT > 0.7)]

CD8T_naive_CD31_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                             protDf$CD045RA_PROT > 0.75 &
                                                             protDf$CD045RO_PROT < 0.75 &
                                                             protDf$CD028_PROT > 0.4 & 
                                                             protDf$CD057_PROT < 0.4 &
                                                             protDf$CD062L_PROT > 0.75 &
                                                             protDf$CD031_PROT > 0.75)]

CD8T_naive_CD38_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                             protDf$CD045RA_PROT > 0.75 &
                                                             protDf$CD045RO_PROT < 0.75 &
                                                             protDf$CD028_PROT > 0.4 & 
                                                             protDf$CD057_PROT < 0.4 &
                                                             protDf$CD062L_PROT > 0.75 &
                                                             protDf$CD038_PROT > 1)]

CD8T_naive_CD161_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                             protDf$CD045RA_PROT > 0.75 &
                                                             protDf$CD045RO_PROT < 0.75 &
                                                             protDf$CD028_PROT > 0.4 & 
                                                             protDf$CD057_PROT < 0.4 &
                                                             protDf$CD062L_PROT > 0.75 &
                                                             protDf$CD161_PROT > 0.7)]

CD8T_CM_CD25_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                   protDf$CD045RA_PROT < 0.75 &
                                                   protDf$CD045RO_PROT > 0.75 &
                                                   protDf$CD057_PROT < 0.4 &
                                                   protDf$CD062L_PROT > 0.75 &
                                                   protDf$CD025_PROT > 0.2)]

CD8T_EM_CD25_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                   protDf$CD045RA_PROT < 0.75 &
                                                   protDf$CD045RO_PROT > 0.75 &
                                                   protDf$CD062L_PROT < 0.75 &
                                                   protDf$CD025_PROT > 0.2)]

CD8T_CM_CD38_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                   protDf$CD045RA_PROT < 0.75 &
                                                   protDf$CD045RO_PROT > 0.75 &
                                                   protDf$CD057_PROT < 0.4 &
                                                   protDf$CD062L_PROT > 0.75 &
                                                   protDf$CD038_PROT > 1)]

CD8T_EM_CD38_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                   protDf$CD045RA_PROT < 0.75 &
                                                   protDf$CD045RO_PROT > 0.75 &
                                                   protDf$CD062L_PROT < 0.75 &
                                                   protDf$CD038_PROT > 1)]

CD8T_EM_CD161_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "CD8T" &
                                                    protDf$CD045RA_PROT < 0.75 &
                                                    protDf$CD045RO_PROT > 0.75 &
                                                    protDf$CD062L_PROT < 0.75 &
                                                    protDf$CD161_PROT > 0.7)]

#Now the gamma-delta
#TCRgd_CM_CD161
#TCRgd_EM_CD161
#TCRgd_naive_CD25
TCRgdDf <- protDf[which(sthlmProtDat$cellType == "gdT"),]
hist(TCRgdDf$CD045RA_PROT, breaks = 50)
abline(v = 1)
hist(TCRgdDf$CD045RO_PROT, breaks = 50)
abline(v = 0.3)
hist(TCRgdDf$CD062L_PROT, breaks = 50)
abline(v = 0.75)
hist(TCRgdDf$CD028_PROT, breaks = 50)
abline(v = 0.4)
hist(TCRgdDf$CD057_PROT, breaks = 50)
abline(v = 0.4)
hist(TCRgdDf$CD025_PROT, breaks = 50)
abline(v = 0.12)
hist(TCRgdDf$CD161_PROT, breaks = 50)
abline(v = 0.5)

TCRgd_naive_CD25_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "gdT" &
                                                       protDf$CD045RA_PROT > 1 &
                                                       protDf$CD045RO_PROT < 0.3 &
                                                       protDf$CD028_PROT > 0.4 & 
                                                       protDf$CD057_PROT < 0.4 &
                                                       protDf$CD062L_PROT > 0.75 &
                                                       protDf$CD025_PROT > 0.12)]

TCRgd_CM_CD161_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "gdT" &
                                                     protDf$CD045RA_PROT < 1 &
                                                     protDf$CD045RO_PROT > 0.3 &
                                                     protDf$CD057_PROT < 0.4 &
                                                     protDf$CD062L_PROT < 0.75 &
                                                     protDf$CD161_PROT > 0.5)]

TCRgd_EM_CD161_Counts <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "gdT" &
                                                     protDf$CD045RA_PROT < 1 &
                                                     protDf$CD045RO_PROT > 0.3 &
                                                     protDf$CD062L_PROT < 0.75 &
                                                     protDf$CD161_PROT > 0.5)]

#Abd finally, the double negative population
#"CD4negCD8negT"  
DNT_Count <- sthlmProtDat$sample[which(sthlmProtDat$cellType == "DNT")]

#Now, these are all put together.
allGatedCountList <- list("B_IgDnegCD27neg" = B_IgDnegCD27neg_Counts, 
                          "B_IgDposCD27negCD24posCD38posTrans" = B_IgDposCD27negCD24posCD38posTrans_Counts,
                          "CD4T_naive_CD25" = CD4T_naive_CD25_Counts,
                          "CD4T_EM_CD161" = CD4T_EM_CD161_Counts,
                          "CD8T_naive" = CD8T_naive_Counts,
                          "CD8T_naive_CD7" = CD8T_naive_CD7_Counts,
                          "CD8T_naive_CD27" = CD8T_naive_CD27_Counts,
                          "CD8T_naive_CD31" = CD8T_naive_CD31_Counts,
                          "CD8T_naive_CD38" = CD8T_naive_CD38_Counts,
                          "CD8T_naive_CD161" = CD8T_naive_CD161_Counts,
                          "CD8T_CM_CD25" = CD8T_CM_CD25_Counts,
                          "CD8T_CM_CD38" = CD8T_CM_CD38_Counts,
                          "CD8T_EM_CD25" = CD8T_EM_CD25_Counts,
                          "CD8T_EM_CD38" = CD8T_EM_CD38_Counts,
                          "CD8T_EM_CD161" = CD8T_EM_CD161_Counts,
                          "TCRgd_naive_CD25" = TCRgd_naive_CD25_Counts,
                          "TCRgd_CM_CD161" = TCRgd_CM_CD161_Counts,
                          "TCRgd_EM_CD161" = TCRgd_EM_CD161_Counts,
                          "DNT" = DNT_Count)

allGatedCounts <- as.data.frame(do.call("cbind", lapply(allGatedCountList, function(x){
  locDatFac <- factor(x, levels = c("1_1", "1_2", "1_3", "1_4", "1_5", "1_6",
                                    "1_7", "1_8", "2_1", "2_2", "2_3", "2_4",
                                    "2_5", "2_6", "2_7", "2_8"))
  locTab <- as.matrix(table(locDatFac))
  
})))

colnames(allGatedCounts) <- names(allGatedCountList)

#This is now turned into fractions of the total number per donor
perDonTotCounts <- sapply(row.names(allGatedCounts), function(x){
  length(which(sthlmProtDat$sample == x))
})

allGatedFracs <- apply(allGatedCounts, 2, function(x){
  x/perDonTotCounts
})
#Now, we are splitting these into LOMG and EOMG. 
#First, exclusion of the uncertain. 
groupVec <- sapply(row.names(allGatedFracs), function(x){
  sthlmProtDat$MG_type[which(sthlmProtDat$sample == x)][1]
})

eomgDat <- allGatedFracs[which(groupVec == "EOMG"),]
lomgDat <- allGatedFracs[which(groupVec == "LOMG"),]

#Now, to the loop. However, as we do know what direction we are expecting
#Of the p-values again, we will run a set of one-tailers. 
eomgtoLomgSize <- pResLow$EOMG_to_LOMG
eomgtoLomgSizeOrd <- eomgtoLomgSize[match(colnames(allGatedFracs), rownames(pResLow))]
names(eomgtoLomgSizeOrd) <- colnames(allGatedFracs)

pRes <- do.call("rbind", lapply(colnames(allGatedFracs), function(x){

  eoLoP <- wilcox.test(eomgDat[,x], lomgDat[,x], 
                       alternative = eomgtoLomgSizeOrd[x], 
                       exact = FALSE)$p.value
}))
rownames(pRes) <- colnames(allGatedFracs)
colnames(pRes) <- "eomgLomgP"

#Now, we zoom in on the ones that survive: 
pResLowStock <- pRes[which(pRes[,1] < 0.05),]
pResLowStock
#     CD8T_naive  CD8T_naive_CD7 CD8T_naive_CD27 CD8T_naive_CD31 CD8T_naive_CD38 
#     0.03060732      0.03060732      0.03060732      0.02269396      0.01162696 
#  CD8T_EM_CD161 
#     0.01193422

dir.create("Results/Data/Oxford_and_Stockholm/Gating")
#These are saved. 
write.csv(pResLowStock, "Results/Data/Oxford_and_Stockholm/Gating/Sign_in_Stock_and_Ox.csv")

#And we also save the frequencies
allGatedFracs <- as.data.frame(allGatedFracs)
allGatedFracs$Group <- groupVec
write.csv(allGatedFracs, "Results/Data/Oxford_and_Stockholm/Gating/freqTable_Stockholm.csv")



