#Here, we are first going to create an PLS-DA analysis for the Oxford patients
#based on the percentages of the four identified clusters. We will create a 2D PLS-DA, 
#as we will include the controls as well. Then, we are going to run a ROC analysi, 
#identifying the optimal separation between the patient groups here. Finally, we will
#use the same PLS-DA grid to separate the patients in the Stockholm data. 

library(mixOmics)
library(caret)
library(plotROC)
oxEuclidAllFreqs <- read.csv("Results/Data/Oxford/Analysis_post_smooth/freqTableAllPops.csv",
                             row.names = 1)
stockEuclidAllFreqs <- read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/freqTable_Stockholm.csv",
                                row.names = 1)
signPops <- suppressWarnings(read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", 
                                      row.names = 1))
signPopNames <- row.names(signPops)
#We start by reducing the data to the populations and donors of interest. 
oxEuclidFocFreqs <- oxEuclidAllFreqs[,which(colnames(oxEuclidAllFreqs) %in% signPopNames),]
stockEuclidFocFreqs <- stockEuclidAllFreqs[,which(colnames(stockEuclidAllFreqs) %in% signPopNames)]

#Already here, we log10-transforming the data has no effect on the outcome. sqrt
#transformation does not have any effect on the ROC, even if the data looks more
#tidy with it. So just to simplify things, we do not transform the data here 
#at all.
#oxEuclidFocFreqs <- sqrt(oxEuclidFocFreqs)
#stockEuclidFocFreqs <- sqrt(stockEuclidFocFreqs)

#And here, we separate the data into the groups we are interested in. 
oxEOMG <- oxEuclidFocFreqs[grep("EOMG.+pre", row.names(oxEuclidAllFreqs)),]
oxLOMG <- oxEuclidFocFreqs[grep("LOMG........+pre", row.names(oxEuclidAllFreqs)),]
oxCtrl <- oxEuclidFocFreqs[grep("Ctrl", row.names(oxEuclidAllFreqs)),]

oxUncertain <- oxEuclidFocFreqs[grep("nclear........+pre", row.names(oxEuclidAllFreqs)),]
oxPost <- oxEuclidFocFreqs[grep("post", row.names(oxEuclidAllFreqs)),]

stockEOMG <- stockEuclidFocFreqs[grep("EOMG", row.names(stockEuclidAllFreqs)),]
stockLOMG <- stockEuclidFocFreqs[grep("LOMG", row.names(stockEuclidAllFreqs)),]
stockUncertain <- stockEuclidFocFreqs[grep("uncertain", row.names(stockEuclidAllFreqs)),]

#Now, we create the PLS-DA
plsDaRes <- mixOmics::plsda(rbind(oxEOMG, oxLOMG), 
                  Y = c(rep("EOMG", nrow(oxEOMG)), 
                        rep("LOMG", nrow(oxLOMG))),
                  ncomp = 1)
#Now, we make a few predictions. First, the uncertain cases as well as the 
#post-thymectomy ones in the Oxford data. 
plsDaPredOx <- predict(plsDaRes, rbind(oxUncertain, oxPost, oxCtrl))

plsDaPred <- predict(plsDaRes, rbind(stockEOMG, stockLOMG, stockUncertain))

#We will not use these predictions directly, as they make no sense when the
#input is 2D, but we will use the PLS vectors. 

#Now, this is plotted, just to get an overview first. 
oxMeta <- read.csv("Data/Oxford/Documentation/Metadata_per_id.csv")

oxPlotDat <- data.frame("Sample" = row.names(rbind(oxEOMG, oxLOMG, oxUncertain, oxPost, oxCtrl)),
                        "PLSDA1" = c(plsDaRes$variates$X[,1], plsDaPredOx$variates[,1]),
                        "Group" = c(rep("EOMG", nrow(oxEOMG)),
                                    rep("LOMG", nrow(oxLOMG)),
                                    rep("Uncertain", nrow(oxUncertain)),
                                    rep("EOMG_post_thym", nrow(oxPost)),
                                    rep("Ctrl", nrow(oxCtrl))))

oxPlotDat$Age <- as.numeric(gsub(".+_age_|_.+_p.+", "", oxPlotDat$Sample))
oxPlotDat$Sex <- sapply(oxPlotDat$Sample, function(x){
  oxMeta$Sex[which(oxMeta$Label == gsub(".+_age_.._|_p.+", "", x))]
})

#And we refine the controls here: 
oxPlotDat$Group[which(oxPlotDat$Group == "Ctrl" & oxPlotDat$Age < 50)] <- "Young_ctrl"
oxPlotDat$Group[which(oxPlotDat$Group == "Ctrl" & oxPlotDat$Age > 50)] <- "Old_ctrl"

oxPlotDat$Group <- factor(oxPlotDat$Group, levels = c("EOMG", "EOMG_post_thym",
                                                      "Young_ctrl", "LOMG",
                                                      "Old_ctrl", "Uncertain"))
p <- ggplot(oxPlotDat, aes(x = PLSDA1, 
                           y = Age, 
                           color = Group,
                           shape = Sex)) +
  geom_point(, size = 6) + theme_bw() +
  scale_color_manual(values = c("red", "red4", "grey", "blue", "#555555", "#AA00FF")) 
p


#The separation is of course along this axis
oxROCDat <- oxPlotDat[which(grepl("OMG", oxPlotDat$Group) &
                              grepl("thym", oxPlotDat$Group)==FALSE),]

scoreVals <- seq(min(oxROCDat$PLSDA1), 
                 max(oxROCDat$PLSDA1), length.out = 50)

conMatList <- lapply(scoreVals, function(z){
  locPred <- rep("EOMG", nrow(oxROCDat))
  locPred[which(oxROCDat$PLSDA1 > z)] <- "LOMG"
  confusionMatrix(factor(locPred, levels = c("EOMG", "LOMG")), 
                  factor(oxROCDat$Group, levels = c("EOMG", "LOMG")),
                  positive = "EOMG")
})
balancedAccList <- lapply(conMatList, function(x){
  x$byClass["Balanced Accuracy"]
})

#If there are multiple values with the maximal accuracy, we choose the middle one.
maxVal <- max(unlist(balancedAccList))
maxValPoss <- which(balancedAccList == maxVal)
#midPoss <- maxValPoss[round(length(maxValPoss)/2)]
bestScoreThreshold <- scoreVals[maxValPoss[1]]
conMatList[[max(which(balancedAccList == maxVal))]]

#So this works well: 
#          Reference
#Prediction EOMG LOMG
#      EOMG   10    0
#      LOMG    2   16
#            Sensitivity : 0.8333           
#            Specificity : 1.0000
#Balanced Accuracy : 0.9167 

#And if that threshold is tried on the Stockholm data: 
trueLength <- nrow(plsDaPred$variates)-2
stockPredY <- rep("EOMG", nrow(plsDaPred$variates))
stockPredY[which(plsDaPred$variates[,1] > bestScoreThreshold)] <- "LOMG"
stockPredY <- stockPredY[1:trueLength]

confusionMatrix(factor(stockPredY, 
                       levels = c("EOMG", "LOMG")), 
                factor(c(rep("EOMG", nrow(stockEOMG)),
                         rep("LOMG", nrow(stockLOMG))), levels = c("EOMG", "LOMG")),
                positive = "EOMG")

#          Reference
#Prediction EOMG LOMG
#EOMG    8    1
#LOMG    0    5

#And now, if the same splsDA threshold is used to separate the young from the 
#old controls in the Oxford data:
oxCtrlDat <- oxPlotDat[grep("ctrl", oxPlotDat$Group),]

oxCtrlPredY <- rep("Young_ctrl", nrow(oxCtrlDat))
oxCtrlPredY[which(oxCtrlDat$PLSDA1 > bestScoreThreshold)] <- "Old_ctrl"

confusionMatrix(factor(oxCtrlPredY, 
                       levels = c("Young_ctrl", "Old_ctrl")), 
                factor(oxCtrlDat$Group, 
                       levels = c("Young_ctrl", "Old_ctrl")),
                positive = "Young_ctrl")
#            Reference
#Prediction   Young_ctrl Old_ctrl
#  Young_ctrl         10        4
#  Old_ctrl            0        6
#Sensitivity : 1.0000          
#Specificity : 0.6000

#Now, we are going to greate a ROC curve for the
#Oxford data, which will show why we end up in the place we do end up in for the 
#Stockholm cohort. As there is a perfect separation in the Stockholm cohort, 
#it is not meaningful to show that ROC curve. 
oxROCDat$GroupNumeric <- 0
oxROCDat$GroupNumeric[which(oxROCDat$Group == "LOMG")] <- 1

basicPlot <- ggplot(oxROCDat, aes(d = GroupNumeric, m = PLSDA1)) + 
  geom_roc(n.cuts = 0)

basicPlot + style_roc(xlab = "1 - Specificity") +
  annotate("text", x = .75, y = .25, 
           label = paste("AUC =", round(calc_auc(basicPlot)$AUC, 2)))
ggsave("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/ROC_PLSDA1_Oxford.pdf", width = 5, height = 5)

#Now, this information is used to plot everything for both the Oxford and the
#Stockholm data. 
#Oxford first
p <- ggplot(oxPlotDat, aes(x = PLSDA1, 
                           y = Age, 
                           color = Group,
                           shape = Sex)) +
  geom_point(, size = 6) + theme_bw() +
  scale_color_manual(values = c("red", "red4", "grey", "blue", "#555555", "#AA00FF")) +
  geom_vline(xintercept = bestScoreThreshold, linetype="dashed",
             color = "black", linewidth=1) +
  scale_x_continuous(limits=c(-4.5, 2.2), expand = c(0, 0))
p
ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/Used_Oxford_Euclid_PLSDA.pdf"))

p + theme(legend.position = "none")
ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/Used_Oxford_Euclid_PLSDA_stripped.pdf"), 
       width = 5, height = 5)

#And here the Stockholm data. 

stockMeta <- read.csv("Data/Stockholm/Metadata/Metadata_combined_used.csv")

stockPlotDat <- data.frame("Sample" = c(row.names(stockEuclidFocFreqs)[grep("EOMG", row.names(stockEuclidAllFreqs))],
                                     row.names(stockEuclidFocFreqs)[grep("LOMG", row.names(stockEuclidAllFreqs))],
                                     row.names(stockEuclidFocFreqs)[grep("uncertain", row.names(stockEuclidAllFreqs))]),
                        "PLSDA1" = c(plsDaPred$variates[,1]), 
                        "Group" = c(rep("EOMG", nrow(stockEOMG)),
                                    rep("LOMG", nrow(stockLOMG)),
                                    rep("uncertain", nrow(stockUncertain))))
stockPlotDat$Age <- as.numeric(gsub(".+_age_|_._.", "", stockPlotDat$Sample))

stockMetDat <- as.data.frame(do.call("cbind", lapply(c("Years_since_thymectomy", "Sex"), function(x){
  sapply(stockPlotDat$Sample, function(y){
    stockMeta[grep(gsub(".+_age_.._|", "", y), stockMeta$Sample_ID), which(colnames(stockMeta) == x)]
  })
})))
colnames(stockMetDat) <- c("Years_since_thymectomy", "Sex")

stockPlotDat$Sex <- stockMetDat$Sex

#And we refine the groups with the info re the thymectomy. 
stockPlotDat$Group[which(stockMetDat$Years_since_thymectomy > 0)] <- "EOMG_post_thym"

stockPlotDat$Group <- factor(stockPlotDat$Group, 
                             levels = c("EOMG",
                                        "EOMG_post_thym",
                                        "LOMG",
                                        "uncertain"))

p <- ggplot(stockPlotDat, aes(x = PLSDA1, 
                           y = Age, 
                           color = Group,
                           shape = Sex)) +
  geom_point(, size = 6) + theme_bw() +
  scale_color_manual(values = c("red", "red4", "blue", "#AA00FF")) +
  geom_vline(xintercept = bestScoreThreshold, linetype="dashed",
             color = "black", linewidth=1) +
  scale_x_continuous(limits=c(-4.5, 2.2), expand = c(0, 0))
p
ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/Used_Stockholm_Euclid_PLSDA_on_Oxford_grid.pdf"))

p + theme(legend.position = "none")
ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/Used_Stockholm_Euclid_PLSDA_on_Oxford_grid_stripped.pdf"), 
       width = 5, height = 5)

#Here is the best score threshold: 
bestScoreThreshold
#-0.129528

#And the data is saved. 
oxPlotDat$scoreThreshold <- stockPlotDat$scoreThreshold <- bestScoreThreshold
write.csv(oxPlotDat, "Results/Data/For_figure_file/Figure_4D_PLSDA_data_UK.csv")
write.csv(oxROCDat, "Results/Data/For_figure_file/Figure_4E_ROC_PLSDA_data_UK.csv")
write.csv(stockPlotDat, "Results/Data/For_figure_file/Figure_4F_PLSDA_data_SE.csv")



