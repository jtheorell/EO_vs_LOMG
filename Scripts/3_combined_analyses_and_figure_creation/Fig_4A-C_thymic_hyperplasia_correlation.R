#Here, as a final part of the final figure, we ask the simple question: is there
#a correlation between the size of the thymic hyperplasia and any of the identified
#EOMGlow populations?
library(ggplot2)
#We will create a framework for comparison, which is the major cell populations. 
popVec <- c("B", "CD4T", "CD8T", "TCRgd", "NK", "ILC")

freqTabMajorPops <- do.call("rbind", lapply(popVec, function(x){
  print(x)
  big_all_file_with_pJor <- readRDS(paste0("../External/Oxford/Resulting_data_15_neighbors/", x, "/", x,  "_full_file.rds"))
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

#After closer consideration, it is only meaningful to do this for the EOMG pop as 
#we lack data for the others.. 
cortTestRes <- cor.test(freqTabThy$hyperNumeric, log(freqTabThy[,which(colnames(freqTabThy) 
                                                                       == "NK_EOMGlow_1")]),
                        method = "pearson")
p.adjust(cortTestRes$p.value, method = "fdr", n = 2) #0.1119545
cortTestRes$estimate #-0.5901353

#We also need to investigate whether this is the case for the major cell populations
corTestList <- lapply(popVec, function(x){
  cor.test(freqTabThy$hyperNumeric, log(freqTabThy[,which(colnames(freqTabThy) == x)]),
           method = "pearson")
})
names(corTestList) <- popVec
round(p.adjust(unlist(lapply(corTestList, function(x) x$p.value)), method = "fdr"), 3)
#    B  CD4T  CD8T TCRgd    NK   ILC 
#0.905 0.905 0.322 0.035 0.322 0.251
#So a seeming significance for gamma-deltas. 

round(unlist(lapply(corTestList, function(x) x$estimate)), 3)
#    B.rho  CD4T.rho  CD8T.rho TCRgd.rho    NK.rho   ILC.rho
#    0.041    -0.062     0.407     0.768    -0.418     0.544

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
#It is only really meaningful to look at the EOMG cluster. 

cortTestRes <- cor.test(freqTabThySamps$hyperNumeric, 
                        log(freqTabThySamps[,which(colnames(freqTabThySamps) == "NK_EOMGlow_1")]),
                        method = "pearson")
p.adjust(cortTestRes$p.value, method = "fdr", n = 2) #0.1016018
cortTestRes$estimate #-0.6646862

#Now, what about the major cell populations? 
corTestListThy <- lapply(popVec, function(x){
  cor.test(freqTabThySamps$hyperNumeric, log(freqTabThySamps[,which(colnames(freqTabThySamps) == x)]),
           method = "pearson")
})
names(corTestListThy) <- popVec

round(p.adjust(unlist(lapply(corTestListThy, function(x) x$p.value)), method = "fdr"), 3)
#    B  CD4T  CD8T TCRgd    NK   ILC 
#0.305 0.305 0.305 0.298 0.406 0.305

round(unlist(lapply(corTestListThy, function(x) x$estimate)), 3)
#    B.cor  CD4T.cor  CD8T.cor TCRgd.cor    NK.cor   ILC.cor 
#    0.505     0.425    -0.478    -0.667    -0.317    -0.441

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

#And this data is saved. 
write.csv(freqTabNKLong, "Results/Data/For_figure_file/Figure_4C_PBMC_thymic_hyperplasia_corr.csv", row.names = FALSE)

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

write.csv(freqTabNKLongThymus, "Results/Data/For_figure_file/Figure_4C_thymus_thymic_hyperplasia_corr.csv", row.names = FALSE)

#At this stage, we should also obviously correlate thymus cell numbers to blood cell numbers
freqTabThySamps <- freqTabBig[grep("_thy", row.names(freqTabBig)),]

freqTabPreThySampCorr <- freqTabPre[which(freqTabPre$Donor %in% freqTabThySamps$Donor),]

identical(freqTabThySamps$Donor, freqTabPreThySampCorr$Donor) #TRUE
identical(colnames(freqTabThySamps), colnames(freqTabPreThySampCorr)) #TRUE

#Now, does cluster 49 correlate? 
corTestList <- lapply(c("NK_EOMGlow_1", "NK"), function(x){
    cor.test(log(freqTabThySamps[,which(colnames(freqTabThySamps) == x)]), 
             log(freqTabPreThySampCorr[,which(colnames(freqTabPreThySampCorr) == x)]),
             method = "spearman")
})

round(unlist(lapply(corTestList, function(x) x$p.value)), 5) 
#NK_EOMGlow_1: 0.005557 
#NK: 0.3488

lapply(corTestList, function(x) x$estimate) 
#NK_EOMGlow_1: 0.830303
#NK: 0.3333333

#What about the other subsets?
corTestList <- lapply(c("B", "CD4T", "CD8T", "TCRgd", "ILC"), function(x){
  cor.test(freqTabThySamps[,which(colnames(freqTabThySamps) == x)], 
           freqTabPreThySampCorr[,which(colnames(freqTabPreThySampCorr) == x)],
           method = "spearman")
})

round(unlist(lapply(corTestList, function(x) x$p.value)), 5) 
#0.58354 0.06647 0.86475 0.34885 0.40695 without correction. 

round(unlist(lapply(corTestList, function(x) x$estimate)), 5) 
#    cor      cor      cor      cor      cor 
# 0.20000  0.61212 -0.06667 -0.33333 -0.29697

#And now, we plot the NK/NK and NK_EOMGlow_1 PBMC/thymus correlations. 
eomgNKPlotDat <- data.frame("Thymus" = freqTabThySamps[,which(colnames(freqTabThySamps) == "NK_EOMGlow_1")],
                            "PBMC" = freqTabPreThySampCorr[,which(colnames(freqTabPreThySampCorr) == "NK_EOMGlow_1")])
ggplot(eomgNKPlotDat, aes(x = log10(Thymus), y = log10(PBMC))) +
  geom_point(size = 4, color = "black") + theme_bw() + theme(legend.position = "none")
  
ggsave("Results/Graphics/Oxford/Thymus/NK_EOMGlow_1_PBMC_thymus_corr.pdf",
       width = 5, height = 5)

write.csv(eomgNKPlotDat, "Results/Data/For_figure_file/Figure_4A_PBMC_thymus_corr_cluster_49.csv", row.names = FALSE)


eomgNKPlotDat <- data.frame("Thymus" = freqTabThySamps[,which(colnames(freqTabThySamps) == "NK")],
                            "PBMC" = freqTabPreThySampCorr[,which(colnames(freqTabPreThySampCorr) == "NK")])
ggplot(eomgNKPlotDat, aes(x = log10(Thymus), y = log10(PBMC))) +
  geom_point(size = 4, color = "grey") + theme_bw() + theme(legend.position = "none")
ggsave("Results/Graphics/Oxford/Thymus/NK_all_PBMC_thymus_corr.pdf",
       width = 5, height = 5)

write.csv(eomgNKPlotDat, "Results/Data/For_figure_file/Figure_4B_PBMC_thymus_corr_all_NK.csv", row.names = FALSE)

