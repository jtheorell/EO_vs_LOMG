#Here, as a final part of the final figure, we ask the simple question: is there
#a correlation between the size of the thymic hyperplasia and any of the identified
#EOMGlow populations?
library(ggplot2)
freqTabBig <- read.csv("Results/Data/Oxford/Analysis_post_smooth/freqTableAllPops.csv",
                       row.names = 1)
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

#Now, are there any correlations to be found here? THis is question one. The
#seond one is obvious, but it comes later. 
corTestList <- lapply(row.names(signPops), function(x){
  cor.test(freqTabThy$hyperNumeric, freqTabThy[,which(colnames(freqTabThy) == x)])
})

names(corTestList) <- row.names(signPops)

round(unlist(lapply(corTestList, function(x) x$p.value)), 3)
#B_LOMGlow_2 CD4T_EOMGlow_3 CD8T_LOMGlow_2   NK_EOMGlow_1 
#0.313          0.398          0.200          0.022

round(unlist(lapply(corTestList, function(x) x$estimate)), 3)

#   B_LOMGlow_2.cor CD4T_EOMGlow_3.cor CD8T_LOMGlow_2.cor   NK_EOMGlow_1.cor 
#            -0.319              0.269              0.398             -0.649

#So there is a considerable negative correlation between the size of the NK 
#populations and the grade of hyperplasia. 

#Very beautiful. Now, what about the THYMIC SAMPLES THAT WE HAVE?
freqTabThySamps <- freqTabBig[grep("_thy", row.names(freqTabBig)),]

freqTabThySamps <- freqTabThySamps[-which(is.na(freqTabThySamps$hyperplasia)),]
freqTabThySamps$hyperNumeric <- sapply(freqTabThySamps$hyperplasia, switch, 
                                  "None" = 0,
                                  "Minimal" = 1, 
                                  "Moderate" = 2, 
                                  "Major" = 3)

#And now for the correlations. 
corTestListThy <- lapply(row.names(signPops), function(x){
  cor.test(freqTabThySamps$hyperNumeric, freqTabThySamps[,which(colnames(freqTabThySamps) == x)])
})

names(corTestListThy) <- row.names(signPops)

round(unlist(lapply(corTestListThy, function(x) x$p.value)), 3)
#   B_LOMGlow_2 CD4T_EOMGlow_3 CD8T_LOMGlow_2   NK_EOMGlow_1 
#         0.243          0.101          0.484          0.024

round(unlist(lapply(corTestListThy, function(x) x$estimate)), 3)
#   B_LOMGlow_2.cor CD4T_EOMGlow_3.cor CD8T_LOMGlow_2.cor   NK_EOMGlow_1.cor 
#             0.434              0.581              0.269             -0.736 

plot(freqTabThySamps$NK_EOMGlow_1, freqTabThySamps$hyperNumeric)

#And it holds true. 
#Now, we will plot the most significant populatin which is NK_EOMGlow_1
dir.create("Results/Graphics/Oxford/Thymus")
ggplot(freqTabThy, aes(x = NK_EOMGlow_1, y = hyperNumeric)) +
  geom_point(size = 4) + theme_bw()
ggsave("Results/Graphics/Oxford/Thymus/NK_EOMGlow_1_hyperplasia_correlation.pdf",
       width = 2.5, height = 5)

ggplot(freqTabThySamps, aes(x = NK_EOMGlow_1, y = hyperNumeric)) +
  geom_point(size = 4) + theme_bw()
ggsave("Results/Graphics/Oxford/Thymus/NK_EOMGlow_1_in_thymus_hyperplasia_correlation.pdf",
       width = 2.5, height = 5)


