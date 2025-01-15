library(ggplot2)
library(reshape2)
#Here, we correlate the populations to the AChR concentration and QMG
stockMeta <- read.csv("Data/Stockholm/Metadata/Metadata_combined_used.csv")

#First, we create a new group variable here. 
stockMeta$newGroup <- stockMeta$MG_type
stockMeta$newGroup[which(stockMeta$MG_type == "LOMG" & stockMeta$Age == 50)] <- "uncertain"
stockMeta$newGroup[which(stockMeta$MG_type == "LOMG" & is.na(stockMeta$Years_since_thymectomy) == FALSE)] <- "uncertain"
stockMeta$newGroup[which(stockMeta$newGroup == "EOMG" & stockMeta$Years_since_thymectomy > 0)] <- "EOMG_post_thym"

allFreq <- read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/freqTable_Stockholm.csv", row.names = 1)

signPops <- row.names(read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1))

signFreq <- allFreq[,which(colnames(allFreq) %in% signPops)]

signFreq$Id <- gsub(".+_age_.._|", "", row.names(signFreq))

signFreq$QMG <- sapply(signFreq$Id, function(x){
  stockMeta$QMG[which(stockMeta$Sample_ID == x)]
})

signFreq$AChR_Abs <- sapply(signFreq$Id, function(x){
  stockMeta$AChR_Abs[which(stockMeta$Sample_ID == x)]
})

signFreq$Group <- sapply(signFreq$Id, function(x){
  stockMeta$newGroup[which(stockMeta$Sample_ID == x)]
})

signFreq$Sex <- sapply(signFreq$Id, function(x){
  stockMeta$Sex[which(stockMeta$Sample_ID == x)]
})

#And now we make the six correlations
apply(signFreq[,1:3], 2, function(x){
  locRes <- cor.test(x, signFreq$AChR_Abs, method = "spearman", exact = FALSE)
  c(locRes$estimate, "p.value" = locRes$p.value)
})
#        CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#rho         -0.0500000     0.09705882  -0.54117647
#p.value      0.8564875     0.72135801   0.03270894

apply(signFreq[,1:3], 2, function(x){
  locRes <- cor.test(x, signFreq$QMG, method = "spearman", exact = FALSE)
  c(locRes$estimate, "p.value" = locRes$p.value)
})
#        CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#rho        -0.43311740    -0.61789445    0.3311205
#p.value     0.09377491     0.01074743    0.2102955

p.adjust(c(0.856487, 0.72135801, 0.03270894, 
           0.09377491, 0.01074743,0.2102955), 
         method = "fdr")

#0.85648700 0.85648700 0.09812682 0.18754982 0.06448458 0.31544325
#So nothing is significant, but it is reasonably close at least. 

#We need to refine these analyses, so that they show what happens within the
#group of interest, as we otherwise get artificial correlations, due to the higher
#QMG titres in the EOMG group. 

eomgFreq <- signFreq[which(signFreq$Id %in% stockMeta$Sample_ID[which(stockMeta$MG_type == "EOMG" & stockMeta$Thymus == "Hyperplasia")]),]

apply(eomgFreq[,1:3], 2, function(x){
  locRes <- cor.test(x, eomgFreq$AChR_Abs, method = "spearman", exact = FALSE)
  c(locRes$estimate, "p.value" = locRes$p.value)
})
#        CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#rho         -0.3571429     0.04761905   -0.1904762
#p.value      0.3851206     0.91084917    0.6514015

apply(eomgFreq[,1:3], 2, function(x){
  locRes <- cor.test(x, eomgFreq$QMG, method = "spearman", exact = FALSE)
  c(locRes$estimate, "p.value" = locRes$p.value)
})
#        CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#rho         -0.1796439     -0.4550980    0.1916202
#p.value      0.6703443      0.2571921    0.6494102
#This is cool: even within the EOMG population, there is a correlation between disease severity and this subset.

lomgFreq <- signFreq[which(signFreq$Id %in% stockMeta$Sample_ID[which(stockMeta$MG_type == "LOMG" & stockMeta$Age > 50)]),]

apply(lomgFreq[,1:3], 2, function(x){
  locRes <- cor.test(x, lomgFreq$AChR_Abs, method = "spearman", exact = FALSE)
  c(locRes$estimate, "p.value" = locRes$p.value)
})
#        CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#rho         -0.5428571     0.02857143   -0.6571429
#p.value      0.2657026     0.95715452    0.1561749

apply(lomgFreq[,1:3], 2, function(x){
  locRes <- cor.test(x, lomgFreq$QMG, method = "spearman", exact = FALSE)
  c(locRes$estimate, "p.value" = locRes$p.value)
})
#        CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#rho         -0.2648204              0   -0.7944613
#p.value      0.6120552              1    0.0590276

#So now, we are going to plot all of these for a supplementary figure. 
allFreqLong <- melt(signFreq, id.vars = c("Id", "QMG", "AChR_Abs", "Group", "Sex"), value.name = "Frequency")

allFreqLong$Group <- factor(allFreqLong$Group, c("EOMG", "EOMG_post_thym", "LOMG", "uncertain"))
dir.create("Results/Graphics/Stockholm/Correlations_to_clinical")
for(i in signPops){
  locDat <- allFreqLong[which(allFreqLong$variable == i),]
  #The x value range needs to be defined
  maxX <- max(locDat$Frequency)*1.05
  p <- ggplot(locDat, aes(x = Frequency, y = AChR_Abs, color = Group, shape = Sex)) + 
    geom_point(size = 3) + theme_bw() +scale_color_manual(values = c("red","red4", "blue", "#AA00FF")) +
    scale_y_continuous(limits = c(0, 310), expand = c(0,0)) + 
    scale_x_continuous(limits = c(0, maxX), expand = c(0,0)) + 
    theme(legend.position = "none")
  ggsave(paste0("Results/Graphics/Stockholm/Correlations_to_clinical/AChR_agains_freq_of", i, ".pdf"), 
         plot= p)
  p <- p + theme(axis.text.x=element_blank(), 
               axis.ticks.x=element_blank(), 
               axis.text.y=element_blank(), 
               axis.ticks.y=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank()) 
  ggsave(paste0("Results/Graphics/Stockholm/Correlations_to_clinical/AChR_agains_freq_of", i, "_stripped.png"), 
         plot= p, width = 3, height = 3)
  p <- ggplot(locDat, aes(x = Frequency, y = QMG, color = Group, shape = Sex)) + 
    geom_point(size = 3) + theme_bw() +scale_color_manual(values = c("red","red4", "blue", "#AA00FF")) +
    scale_y_continuous(limits = c(0, 20), expand = c(0,0)) + theme(legend.position = "none") +
    scale_x_continuous(limits = c(0, maxX), expand = c(0,0))
  ggsave(paste0("Results/Graphics/Stockholm/Correlations_to_clinical/QMG_agains_freq_of", i, ".pdf"), 
         plot= p)
  p <- p + theme(axis.text.x=element_blank(), 
                 axis.ticks.x=element_blank(), 
                 axis.text.y=element_blank(), 
                 axis.ticks.y=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank()) 
  ggsave(paste0("Results/Graphics/Stockholm/Correlations_to_clinical/QMG_agains_freq_of", i, "_stripped.png"),
         width = 3, height = 3, plot= p)
}

