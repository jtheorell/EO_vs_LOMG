#Here, we are going to start by selecting the top populations by sparse partial
#least squares analysis, and then we are going to hone in on the ones that are 
#significant 
#We start by creating one master frame. In this analysis, part of the work, classifying
#the B and T cells, has been conducted by Dr Chiara Sorini. The second analysis 
#of the BT panelby myself concerns only the populations identified through 
#gating, targeting the clusters of signicicance. 
library(mixOmics)
library(ggplot2)

soriniBT <- read.csv("Data/Oxford/Gating/2024-08-23_BT_Sorini_analysis_JT_changes.csv", row.names = 1)
ILCNK <- read.csv("Data/Oxford/Gating/240819_NK_ILC_cell_numbers.csv", row.names = 1)

#For the soriniBT file, this is a very straightforward matter, as all the daughters
#are called as their parents. 

#We have to start by adding two columns to the ILCNK file. 
ILCNK$ILC <- ILCNK$ILC1+ILCNK$ILC2+ILCNK$ILCpILC3
ILCNK$NK <- ILCNK$CD56bright+ILCNK$Adaptive+ILCNK$CD56dimNonAdapt

#Now, we calculate the fraction of the main population for each individual separately
#and then combine them.
soriniMainPops <- c("B", "CD4T", "CD8T",  "TCRgd")
ILCNKMainPops <- c("ILC", "NK")
soriniList <- list(soriniBT, soriniMainPops)
ILCNKList <- list(ILCNK, ILCNKMainPops)

fracOfSubsetDf <- do.call("cbind", lapply(list(soriniList,  ILCNKList), function(x){
  do.call("cbind", lapply(x[[2]], function(y){
    if(y == "NK"){
      locDat <- x[[1]][,grep("CD56|Adaptive", colnames(x[[1]]))]
    } else {
      locDat <- x[[1]][,grep(y, colnames(x[[1]]))]
    }
    fracDat <- locDat/x[[1]][,which(colnames(x[[1]]) == y)]
  }))
}))

#Now, we also add the label and tissue to this.
fracOfSubsetDf$Label <- soriniBT$Label
fracOfSubsetDf$Tissue <- soriniBT$Tissue

#Time to introduce the other metadata. 
metaData <- read.csv("Data/Oxford/Documentation/Metadata_per_id.csv")

fracOfSubsetDf$Sex <- unlist(lapply(fracOfSubsetDf$Label, function(x){
  metaData$Sex[which(metaData$Label == x)]
}))

fracOfSubsetDf$Age <- unlist(lapply(fracOfSubsetDf$Label, function(x){
  metaData$Age[which(metaData$Label == x)]
}))

fracOfSubsetDf$Group <- unlist(lapply(fracOfSubsetDf$Label, function(x){
  metaData$Group[which(metaData$Label == x)]
}))

#We also make the change to the label of the 50-year old LOMG pat and we
#introduce the two control groups here already
fracOfSubsetDf$Group[which(fracOfSubsetDf$Group == "preLOMG")] <- "uncertain"
fracOfSubsetDf$Group[which(fracOfSubsetDf$Age == 50)] <- "uncertain"
fracOfSubsetDf$Group[which(fracOfSubsetDf$Group == "Ctrl" &
                            fracOfSubsetDf$Age < 50)] <- "Young_ctrl"
fracOfSubsetDf$Group[which(fracOfSubsetDf$Group == "Ctrl" &
                            fracOfSubsetDf$Age > 50)] <- "Old_ctrl"
table(fracOfSubsetDf$Group, fracOfSubsetDf$Tissue)
#             post pre thy
#  EOMG          3  13   7
#  LOMG          0  16   1
#  Old_ctrl      0  10   0
#  uncertain     0   3   2
#  Young_ctrl    0  10   0

#This looks right. So now, lets run an sPLS-DA!
#Before doing so, we need to run some exclusions. 
freqDfPat <- fracOfSubsetDf[-grep("ctrl", fracOfSubsetDf$Group),]

freqDfPre <- freqDfPat [grep("pre", freqDfPat$Tissue),]

nonPreLomg <- freqDfPre[-grep("uncertain", freqDfPre$Group),]
preLomg <- freqDfPre[grep("uncertain", freqDfPre$Group),]

metaDatCols <- which(colnames(nonPreLomg) %in% c("Tissue", colnames(metaData)))

#Now, we run Wilcoxon analyses for all of these.  

lomgDat <- nonPreLomg[which(nonPreLomg$Group == "LOMG"),]
eomgDat <- nonPreLomg[which(nonPreLomg$Group == "EOMG"),]
youngCtrlDat <- fracOfSubsetDf[which(fracOfSubsetDf$Group == "Young_ctrl"),]
oldCtrlDat <- fracOfSubsetDf[which(fracOfSubsetDf$Group == "Old_ctrl"),]

#Now, we do this in a two-step fashion, starting with significance for the 
#EOMC/LOMG comparison
pResMG <- do.call("rbind", lapply(colnames(lomgDat)[-metaDatCols], function(x){
  eoLoP <- wilcox.test(eomgDat[,x], lomgDat[,x], exact = FALSE)$p.value
}))
dimnames(pResMG) <- list(colnames(lomgDat)[-metaDatCols], c("eomgLomgP"))

#Now, given the number of comparisons, we of course need to adjust for
#multiple copmarisons in some way. However, given that we are making three
#comparisons (Pat-pat in Oxford, Pat-ctrl in Oxford, Pat-pat in Stockholm), 
#and we in each case will apply a threshold of 0.05, it corresponds to 
#a p-value of 0.05*0.05*0.05 = 0.000125. Adjusting for 211 comparisons in the
#most conservative way would mean applying the Bonferroni correction, i.e. dividing
#0.05 by 211, which is 0.00023. This means that with our triple-step method, we
#are in fact more conservative even without multiple comparison correction, 
#in the traditional sense. 

pResMGLowP <- pResMG[which(pResMG[,1] < 0.05),]

#These populations are now categorised into such where EOMG is higher or lower. 
EOMGlowHighVec <- sapply(names(pResMGLowP), function(x){
  if(median(eomgDat[,x]) > median(lomgDat[,x])){
    "greater"
  } else if(median(eomgDat[,x]) < median(lomgDat[,x])){
    "less"
  }
})

#This leads us on to the next step, where we now know what to expect of the controls, 
#and we can therefore run one-sided tests. 
pResCtrl <- t(sapply(names(pResMGLowP), function(x){
  if(EOMGlowHighVec[x] == "greater"){
    alternativeVec <- c("greater", "less")
  } else {
    alternativeVec <- c("less", "greater")
  }
  eoCtrlP <- wilcox.test(eomgDat[,x], youngCtrlDat[,x], 
                         alternative = alternativeVec[1], 
                         exact = FALSE)$p.value
  loCtrlP <- wilcox.test(lomgDat[,x], oldCtrlDat[,x], 
                         alternative = alternativeVec[2], 
                         exact = FALSE)$p.value
  c(eoCtrlP, loCtrlP)
}))
pRes <- data.frame(pResCtrl, pResMGLowP, EOMGlowHighVec)
colnames(pRes) <- c("eomgCtrlP", "lomgCtrlP", "eomgLomgP", "EOMG_to_LOMG")

#Now, we further restrict the set by zooming in on the ones that have at least
#one low p-value here. 
lowValVec <- apply(pRes[,1:2], 1, min)

pResLow <- pRes[which(lowValVec < 0.05),]

#And this is saved. 
dir.create("Results_per_celltype/Data/Oxford/Gating", recursive = TRUE)

#In this case, we now include the NK population that has been identified 
#as significant in the other analysis. 
write.csv(pResLow, 
          "Results_per_celltype/Data/Oxford/Gating/Interesting_gated_populations_p_vals_Oxford.csv")

#And we save the patient frequencies 
row.names(fracOfSubsetDf) <- sapply(row.names(fracOfSubsetDf), function(x){
  label <- soriniBT$Label[which(row.names(soriniBT) == x)]
  tissue <- soriniBT$Tissue[which(row.names(soriniBT) == x)]
  paste0(label, "_", tissue)
})
write.csv(fracOfSubsetDf, "Results_per_celltype/Data/Oxford/Gating/freqTableAllPops.csv")


