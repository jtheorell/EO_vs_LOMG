#Here, as a final part of the final figure, we ask the simple question: is there
#a correlation between the size of the thymic hyperplasia and any of the identified
#EOMGlow populations?
library(ggplot2)
#We will create a framework for comparison, which is the major cell populations. 
popVec <- c("B", "CD4T", "CD8T", "TCRgd", "NK", "ILC")

freqTabMajorPops <- do.call("rbind", lapply(popVec, function(x){
  print(x)
  big_all_file_with_pJor <- readRDS(paste0("../External/Oxford/Resulting_data/", x, "/", x,  "_full_file.rds"))
  #EXCLUSION OF PJOR
  #As it is clear that the donor PJor is more or less dead (95% of the cells are
  #gated out at the dump/dead gate stage), it will be excluded here, before going
  #in to the real analyses. 
  big_all_file <- big_all_file_with_pJor[-which(big_all_file_with_pJor$id == "PJor"),]
  freqPerDon <- unlist(lapply(unique(big_all_file$idTime), function(y){
    locCount <- length(which(big_all_file$idTime == y))
    locFreq <- locCount/big_all_file$lymphCount[which(big_all_file$idTime == y)][1]
  }))
  names(freqPerDon) <- unique(big_all_file$idTime)
  freqPerDon
}))
row.names(freqTabMajorPops) <- popVec

freqTabBigIncomplete <- read.csv("Results/Data/Oxford/Analysis_post_smooth/freqTableAllPops.csv",
                       row.names = 1)

#And now, we add in the major populations. 
freqTabMajorPopsT <- t(freqTabMajorPops)

#The rows are now renamed, to match between them
rownames(freqTabMajorPopsT) <- sapply(rownames(freqTabMajorPopsT), function(x) rownames(freqTabBigIncomplete)[grep(x, rownames(freqTabBigIncomplete))])
freqTabMajorPopsTOrd <- freqTabMajorPopsT[match(rownames(freqTabBigIncomplete), rownames(freqTabMajorPopsT)),]

#Are they identical?
identical(rownames(freqTabMajorPopsTOrd), rownames(freqTabBigIncomplete)) #TRUE

#So then, they can be combined. 
freqTabBig <- cbind(freqTabMajorPopsTOrd, freqTabBigIncomplete)

signPops <- read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1)

freqTabBig$Donor <- gsub(".+_age_.._|_.+", "", row.names(freqTabBig))

metaData <- read.csv("Data/Oxford/Documentation/Metadata_per_id.csv")

freqTabBig$hyperplasia <- unlist(lapply(freqTabBig$Donor, function(x){
  metaData$Thymic_hyperplasia[which(metaData$Label == x)]
}))

freqTabPre <- freqTabBig[grep("_pre", row.names(freqTabBig)),]

freqTabThy <- freqTabPre[-which(is.na(freqTabPre$hyperplasia)),]
freqTabThy$hyperNumeric <- sapply(freqTabThy$hyperplasia, switch, 
                                  "None" = 0,
                                  "Minimal" = 1, 
                                  "Moderate" = 2, 
                                  "Major" = 3)

#Now, are there any correlations to be found here? 
#corTestList <- lapply(row.names(signPops), function(x){
#  cor.test(freqTabThy$hyperNumeric, freqTabThy[,which(colnames(freqTabThy) == x)])
#})
#
#names(corTestList) <- row.names(signPops)
#
#round(p.adjust(unlist(lapply(corTestList, function(x) x$p.value)), method = "fdr"), 3)
#CD8T_LOMGlow_1 CD8T_LOMGlow_2   NK_EOMGlow_1 
#         0.673          0.175          0.099

#round(unlist(lapply(corTestList, function(x) x$estimate)), 3)

#CD8T_LOMGlow_1.cor CD8T_LOMGlow_2.cor   NK_EOMGlow_1.cor 
#             0.144              0.441             -0.643

#After closer consideration, it is only meaningful to do this for the EOMG pop. 


#We also need to investigate whether this is the case for the major cell populations
corTestList <- lapply(popVec, function(x){
  cor.test(freqTabThy$hyperNumeric, freqTabThy[,which(colnames(freqTabThy) == x)])
})
names(corTestList) <- popVec
round(p.adjust(unlist(lapply(corTestList, function(x) x$p.value)), method = "fdr"), 3)
#    B  CD4T  CD8T TCRgd    NK   ILC 
#0.868 0.868 0.480 0.251 0.266 0.266
#So a tendency towards 

round(unlist(lapply(corTestList, function(x) x$estimate)), 3)
#    B.cor  CD4T.cor  CD8T.cor TCRgd.cor    NK.cor   ILC.cor 
#    0.057    -0.071     0.331     0.620    -0.482     0.499

#Nothing significant here. 

#ow, what about the thymic samples?
freqTabThySamps <- freqTabBig[grep("_thy", row.names(freqTabBig)),]

freqTabThySamps <- freqTabThySamps[-which(is.na(freqTabThySamps$hyperplasia)),]
freqTabThySamps$hyperNumeric <- sapply(freqTabThySamps$hyperplasia, switch, 
                                  "None" = 0,
                                  "Minimal" = 1, 
                                  "Moderate" = 2, 
                                  "Major" = 3)

#And now for the correlations in the thymus. 
#corTestListThy <- lapply(row.names(signPops), function(x){
#  cor.test(freqTabThySamps$hyperNumeric, freqTabThySamps[,which(colnames(freqTabThySamps) == x)])
#})
#
#names(corTestListThy) <- row.names(signPops)
#names(corTestListThy) <- row.names(signPops)
#round(p.adjust(unlist(lapply(corTestListThy, function(x) x$p.value)), method = "fdr"), 3)
#CD8T_LOMGlow_1 CD8T_LOMGlow_2   NK_EOMGlow_1 
#         0.346          0.364          0.070
#round(unlist(lapply(corTestListThy, function(x) x$estimate)), 3)
#CD8T_LOMGlow_1.cor CD8T_LOMGlow_2.cor   NK_EOMGlow_1.cor 
#             0.444             -0.345             -0.738 

#It is only really meaningful to look at the EOMG cluster. 

cortTestRes <- cor.test(freqTabThySamps$hyperNumeric, freqTabThySamps[,which(colnames(freqTabThySamps) == "NK_EOMGlow_1")])


#Now, what about the major cell populations? 
corTestListThy <- lapply(popVec, function(x){
  cor.test(freqTabThySamps$hyperNumeric, freqTabThySamps[,which(colnames(freqTabThySamps) == x)])
})
names(corTestListThy) <- popVec

round(p.adjust(unlist(lapply(corTestListThy, function(x) x$p.value)), method = "fdr"), 3)
#    B  CD4T  CD8T TCRgd    NK   ILC 
#0.293 0.293 0.293 0.293 0.293 0.293

round(unlist(lapply(corTestListThy, function(x) x$estimate)), 3)
#    B.cor  CD4T.cor  CD8T.cor TCRgd.cor    NK.cor   ILC.cor 
#    0.426     0.395    -0.449    -0.639    -0.486    -0.518

#So, overall, there is a small trend that NK cells in general (as well as gamma-deltas)
#follow the same pattern as the NK_EOMGlow_1 population, but the trends are 
#much weaker, and it is not explained by a general increase in the T-cell subsets, 
#as the CD4 go up a little but the CD8 go down too. 


freqTabNKFoc <- data.frame("Reference" = log10(freqTabThy[,which(colnames(freqTabThy) == "NK")]), 
                           "NK_cluster" = log10(freqTabThy[,which(colnames(freqTabThy) == "NK_EOMGlow_1")]),
                           "HypeprlasiaNumeric" = freqTabThy$hyperNumeric)

freqTabNKLong <- reshape2::melt(freqTabNKFoc, id.vars = "HypeprlasiaNumeric")
#To show all the data in the y direction, we will introduce some noise.
set.seed(11)
freqTabNKLong$HypeprlasiaNumeric <- freqTabNKLong$HypeprlasiaNumeric+rnorm(length(freqTabNKLong$value), sd = 0.03)
dir.create("Results/Graphics/Oxford/Thymus")
ggplot(freqTabNKLong, aes(x = value, y = HypeprlasiaNumeric, colour = variable)) +
  geom_point(size = 4) + theme_bw() + scale_color_manual(values = c("grey", "black")) +
  theme(legend.position = "none") + scale_y_continuous(limits = c(-0.1,3.1))
ggsave("Results/Graphics/Oxford/Thymus/NK_EOMGlow_1_hyperplasia_correlation.pdf",
       width = 5, height = 5)


#The same is of course done for the thymus cells.
freqTabNKFoc <- data.frame("Reference" = log10(freqTabThySamps[,which(colnames(freqTabThySamps) == "NK")]), 
                           "NK_cluster" = log10(freqTabThySamps[,which(colnames(freqTabThySamps) == "NK_EOMGlow_1")]),
                           "HypeprlasiaNumeric" = freqTabThySamps$hyperNumeric)


freqTabNKLongThymus <- reshape2::melt(freqTabNKFoc, id.vars = "HypeprlasiaNumeric")
#To show all the data in the y direction, we will introduce some noise.
set.seed(110)
freqTabNKLongThymus$HypeprlasiaNumeric <- freqTabNKLongThymus$HypeprlasiaNumeric+rnorm(length(freqTabNKLongThymus$value), sd = 0.03)
ggplot(freqTabNKLongThymus, aes(x = value, y = HypeprlasiaNumeric, colour = variable)) +
  geom_point(size = 4) + theme_bw() + scale_color_manual(values = c("grey", "black")) +
  theme(legend.position = "none")  + scale_y_continuous(limits = c(-0.1,3.1))
ggsave("Results/Graphics/Oxford/Thymus/NK_EOMGlow_1_hyperplasia_correlation_thymus_cells.pdf",
       width = 5, height = 5)

#At this stage, we should also obviously correlate thymus cell numbers to blood cell numbers
freqTabThySamps <- freqTabBig[grep("_thy", row.names(freqTabBig)),]

freqTabPreThySampCorr <- freqTabPre[which(freqTabPre$Donor %in% freqTabThySamps$Donor),]

identical(freqTabThySamps$Donor, freqTabPreThySampCorr$Donor) #TRUE
identical(colnames(freqTabThySamps), colnames(freqTabPreThySampCorr)) #TRUE

#Now, does cluster 49 correlate? 
corTestList <- lapply("NK_EOMGlow_1", function(x){
    cor.test(freqTabThySamps[,which(colnames(freqTabThySamps) == x)], 
             freqTabPreThySampCorr[,which(colnames(freqTabPreThySampCorr) == x)])
})

round(unlist(lapply(corTestList, function(x) x$p.value)), 5)

lapply(corTestList, function(x) x$estimate)
