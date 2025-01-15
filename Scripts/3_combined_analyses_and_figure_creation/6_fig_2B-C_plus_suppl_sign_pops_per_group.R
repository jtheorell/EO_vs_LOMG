#Now for some plotting of individuals in groups with the gated and the Euclidean
#data, also including the filtered data . 
library(ggplot2)
library(ggforce)
library(mixOmics)
library(pls)

runVec <- list("gated" = data.frame("focPops" = "Results/Data/Oxford_and_Stockholm/Gating/Sign_in_Stock_and_Ox.csv",
                           "ox" = "Results/Data/Oxford/Gating/freqTableAllPops.csv",
                           "stock" = "Results/Data/Oxford_and_Stockholm/Gating/freqTable_Stockholm.csv"), 
               "neighSmoothed" = data.frame("focPops" = "Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv",
                                   "ox" = "Results/Data/Oxford/Analysis_post_smooth/freqTableAllPops.csv",
                                   "stock" = "Results/Data/Oxford_and_Stockholm/Harmonisation/freqTable_Stockholm.csv"))
oxMeta <- read.csv("Data/Oxford/Documentation/Metadata_per_id.csv")
stockMeta <- read.csv("Data/Stockholm/Metadata/Metadata_combined_used.csv")

dir.create("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA")
for(i in runVec){
  print(i[1])
  focusPopulations <- read.csv(i$focPops)$X
  for(j in c("ox", "stock")){
    print(j)
    locDat <- read.csv(i[,j], row.names = 1)
    if(j == "ox"){
      locNonRows <- grep("thy", row.names(locDat))
      locDatFoc <- locDat[-locNonRows,which(colnames(locDat) %in% focusPopulations)]
      locDatFoc$Group <- locDat$Group[-locNonRows]
      locDatFoc$Group[grep("post", row.names(locDatFoc))] <- "EOMG_post_thym"
      if(grepl("Gating", i$ox)){
        locDatFoc$Age <- locDat$Age[-locNonRows]
        locDatFoc$Sex <- locDat$Sex[-locNonRows]
        locDatFoc$Donor <- gsub("|_.+", "", row.names(locDatFoc))
        plotName <- "Gated"
        widthForGraph <- 14
        ##############################
        #This is a slightly risky piece of hard code that might not work if we change
        #the analysis again, but it needs to be here to create some consistency
        locDatFoc <- locDatFoc[,c(colnames(locDatFoc)[grep("naive", colnames(locDatFoc))],
                                  colnames(locDatFoc)[-grep("naive", colnames(locDatFoc))])]
        
      } else {
        locDatFoc$Age <- sapply(row.names(locDatFoc), function(x){
          oxMeta$Age[which(oxMeta$Label == gsub(".+_age_.._|_.+", "", x))]
        })
        locDatFoc$Sex <- 
          sapply(row.names(locDatFoc), function(x){
            oxMeta$Sex[which(oxMeta$Label == gsub(".+_age_.._|_.+", "", x))]
          })
        locDatFoc$Donor <- gsub(".+_age_.._|_.+", "", row.names(locDatFoc))
        plotName <- "Euclid_cluster"
        widthForGraph <- 14
      }
      locDatFoc$Group[which(locDatFoc$Group == "unclear")] <- "uncertain"
      freqTabLong <- reshape2::melt(locDatFoc, 
                                    id.vars = c("Group", "Donor", "Age", "Sex"))
      
      freqTabLong$Group <- factor(freqTabLong$Group, 
                                  levels = c("EOMG",
                                             "EOMG_post_thym",
                                             "Young_ctrl",
                                             "LOMG",
                                             "Old_ctrl", 
                                             "uncertain"))
      freqTabLong$log10Percent <- log10(freqTabLong$value*100)
      
      if(length(which(is.infinite(freqTabLong$log10Percent))) > 0){
          freqTabLong$log10Percent[which(is.infinite(freqTabLong$log10Percent))] <- 
          -3
      } 
      #There are also two values that are slighly below 0.0001, and these are
      #rounded to this value. 
      if(any(freqTabLong$log10Percent < -3)){
        freqTabLong$log10Percent[which(freqTabLong$log10Percent < -3)] <- -3
      }
     
      sexShape <- sapply(freqTabLong$Sex, function(x) {
        switch(x, "F" =16, "M" = 17) })
      freqTabLong$sexShape <- sexShape[order(freqTabLong$variable, freqTabLong$Group)]
      set.seed(10)
      p <- ggplot(freqTabLong, aes(x = variable,
                                   y = log10Percent,
                                   fill = Group, 
                                   color = Group)) +
        geom_sina(jitter_y = FALSE, 
                  scale = "width", 
                  maxwidth = 0.85, size = 4,
                  shape = freqTabLong$sexShape) +
        geom_violin(alpha = 0, scale = "width", width = 0.85, linewidth = 1, color = "black") +
        theme_bw() + scale_color_manual(values = c("red", "red4", "grey", "blue", "#555555", "#AA00FF"))+
        scale_y_continuous(limits = c(-3, 2), expand = c(0,0))
      p
      ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/", plotName, "_freqs_", j, ".pdf"))
      
      p <- p + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.title = element_blank())
      
      p
      ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/", plotName, "_freqs_", j, "_stripped.pdf"), 
             width = widthForGraph, height = 4)
      
      #Now, for the PLS-DA
      mgDat <- locDatFoc[-grep("ctrl", locDatFoc$Group),]

    } else {
      locDatFoc <- locDat[,which(colnames(locDat) %in% focusPopulations)]
      locDatFoc$Group <- locDat$Group
      if(grepl("Gating", i$stock)){
        locNames <- row.names(locDatFoc)
        widthForGraph <- 8
      } else {
        locNames <- gsub(".+_age_.._|", "", row.names(locDat))
        widthForGraph <- 8
      }
      locDatFoc$Age <- sapply(locNames, function(x){
        stockMeta$Age[which(stockMeta$Sample_ID == x)]
      })
      locDatFoc$Sex <- sapply(locNames, function(x){
        stockMeta$Sex[which(stockMeta$Sample_ID == x)]
      })
      yearsSinceThym <- sapply(locNames, function(x){
        stockMeta$Years_since_thymectomy[which(stockMeta$Sample_ID == x)]
      })
      locDatFoc$Group[which(yearsSinceThym > 0)] <- "EOMG_post_thym"
      
      freqTabLong <- reshape2::melt(locDatFoc, 
                                    id.vars = c("Group", "Age", "Sex"))
      
      freqTabLong$Group <- factor(freqTabLong$Group, 
                                  levels = c("EOMG",
                                             "EOMG_post_thym",
                                             "LOMG",
                                             "uncertain"))
      freqTabLong$log10Percent <- log10(freqTabLong$value*100)
      if(length(which(is.infinite(freqTabLong$log10Percent))) > 0){
        freqTabLong$log10Percent[which(is.infinite(freqTabLong$log10Percent))] <- 
          -3
      }
      #There are also two values that are slighly below 0.0001, and these are
      #rounded to this value. 
      if(any(freqTabLong$log10Percent < -3)){
        freqTabLong$log10Percent[which(freqTabLong$log10Percent < -3)] <- -3
      }
      
      sexShape <- sapply(freqTabLong$Sex, function(x) {
        switch(x, "F" =16, "M" = 17) })
      freqTabLong$sexShape <- sexShape[order(freqTabLong$variable, freqTabLong$Group)]
      
      set.seed(10)
      p <- ggplot(freqTabLong, aes(x = variable,
                                   y = log10Percent,
                                   fill = Group, 
                                   color = Group)) +
        geom_sina(jitter_y = FALSE, 
                  scale = "width", 
                  maxwidth = 0.85, size = 4,
                  shape = freqTabLong$sexShape) +
        geom_violin(alpha = 0, scale = "width", width = 0.85, linewidth = 1, color = "black") +
        theme_bw() + scale_color_manual(values = c("red", "red4", "blue", "#AA00FF"))+
        scale_y_continuous(limits = c(-3, 2), expand = c(0,0))
      p
      ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/", plotName, "_freqs_", j, ".pdf"))
      
      p <- p + theme(legend.position = "none",
                     axis.text = element_blank(),
                     axis.title = element_blank())
      
      p
      ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/", plotName, "_freqs_", j, "_stripped.pdf"), 
             width = widthForGraph, height = 4)
      
      mgDat <- locDatFoc
    }
    #And here comes the rest of the PLS-DA analysis. 
    eoLoDat <- mgDat[-grep("thym|uncertain", mgDat$Group),]
    eoLoDatPLS <- eoLoDat[,-which(colnames(eoLoDat) %in% c("Group", "Donor", "Age", "Sex"))]
    plsDaRes <- plsda(X = eoLoDatPLS, 
                      Y = eoLoDat$Group,
                      ncomp = 1)
    
    #And now we use this for prediction of the two unknown samples and the
    #post thymectomy samples too. 
    otherDat <- mgDat[grep("thym|uncertain", mgDat$Group),]
    otherDatPLS <- otherDat[,-which(colnames(otherDat) %in% 
                                      c("Group", "Donor", "Age", "Sex"))]
    pslDaPred <- predict(plsDaRes, otherDatPLS)
    
    allvariates <- rbind(plsDaRes$variates$X, pslDaPred$variates)
    
    allvariatesOrdered <- allvariates[match(row.names(mgDat), row.names(allvariates)),]
    identical(names(allvariatesOrdered), row.names(mgDat)) #TRUE
    
    mgDat$plsDaVec <- allvariatesOrdered
    mgDat$Group <- factor(mgDat$Group, 
                          levels = c("EOMG",
                                     "EOMG_post_thym",
                                     "LOMG",
                                     "uncertain"))
    p <- ggplot(mgDat, aes(x = plsDaVec, 
                           y = Age, 
                           color = Group,
                           shape = Sex, size = 3)) +
      geom_point() + theme_bw() +
      scale_color_manual(values = c("red","red4", "blue", "#AA00FF"))
    p
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/", plotName, "_sPLSDA_", j, ".pdf"))
    
    p + theme(legend.position = "none")
    ggsave(paste0("Results/Graphics/Oxford_and_Stockholm/Freqs_and_PLSDA/", plotName, "_sPLSDA_", j, "_stripped.pdf"), 
           width = 5, height = 5)
    
  }
  
}


