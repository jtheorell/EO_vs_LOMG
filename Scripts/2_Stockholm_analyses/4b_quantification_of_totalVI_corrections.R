library(DepecheR)
#Here, we are quantifying and more outrightly visualizing the errors introduced by
#totalVI and how these are solved. 
adtCountsCtrl <- read.csv("Data/Stockholm/totalVI_protein_corrections/adtCountsCtrl.csv", row.names = 1)
normAdtCountsCtrl <- read.csv("Data/Stockholm/totalVI_protein_corrections/normAdtCountsCtrl.csv", row.names = 1)
transformed_x <- read.csv("Data/Stockholm/totalVI_protein_corrections/transformed_x.csv", row.names = 1)

#Now, as we in reality only will use a reduced set, we will focus on that set here
markerList <- unique(unlist(list("B" = c("CD38", "CCR7", "CD95", "CD45", "CCR4", "IgD.TCRgd", 
                           "CD11b", "CD27", "CD20", "CD45RA", "CD31", 
                           "CD24", "IgM"),
                   "CD4T" = c("CD38", "KLRG1", "CCR7", "CD95", "CD25", "CD4", "CCR4", "CD71", 
                              "CD11b", "CD27", "CD64", "CD20", "CD45RA", "CD2", "CD31", "CD161", "CD7",
                              "CD24", "CD152"),
                   "CD8T" = c("CD38", "KLRG1", "CCR7", "CD95", "CD25", "CCR4", "CD71", "CD8",
                              "CD11b", "CD27", "CD64", "CD20", "CD45RA", "CD2", "CD31", "CD161", "CD7",
                              "CD24", "CD152"),
                   "TCRgd" = c("CD38", "KLRG1", "CCR7", "CD95", "CD25", "CCR4", "CD71", "CD8",
                               "CD11b", "CD27", "CD64", "CD20", "CD45RA", "CD2", "CD31", "CD161", "CD7",
                               "CD24", "CD152"),
                   "NK" = c("NKp30","HLADR","CCR7", 
                            "CD183","CD161","CD16",
                            "CD56", "CD8",
                            "CD127", "CD27", "CD57",
                            "NKG2A", "CD2", "NKG2C",
                            "CD7","CD218a"), 
                   "ILC" = c("NKp30", "HLADR",
                             "CCR7", "CD183",
                             "CD161", "CD56", 
                             "CD127", "CD27", 
                             "NKG2A", "CD2", 
                             "CRTH2", "CD7",
                             "CD117", "CD218a"),
                   "BT_all" = c("CD38", "KLRG1", "CCR7", "CD95", "CD25", "CD45", 
                                "CD4", "CCR4", "IgD.TCRgd", "CD71","CD8", "CD19",
                                "CD3", "CD11b", "CD27", "CD20", "CD45RA","CD2",
                                "CD31", "CD161", "CD7", "CD24", "IgM", "CD152"),
                   "ILC_NK_all" = c("NKp30", "HLADR", "CCR7", "CD183", "CD161", "CD45", 
                                    "IgD.TCRgd", "CD56", "CD8", "CD19", "CD3", "CD127", 
                                    "CD27", "CD57", "NKG2A", "CD2", "NKG2C", "CRTH2", 
                                    "CD7", "CD117", "CD218a"))))
colnamesInOxFormat <- colnames(adtCountsCtrl)
colnamesInOxFormat <- gsub("_PROT", "", colnamesInOxFormat)
colnamesInOxFormat <- gsub("CD0", "CD", colnamesInOxFormat)
colnamesInOxFormat <- gsub("CD0", "CD", colnamesInOxFormat)

redColList <- colnames(adtCountsCtrl)[which(colnamesInOxFormat %in% markerList)]

adtCountsCtrlRed <- adtCountsCtrl[,redColList]
normAdtCountsCtrlRed <- normAdtCountsCtrl[,redColList]
transformed_xRed <- transformed_x[,redColList]

#If we quantify the differences by comparing the correlation matrices, and rank
#the markers based on this, we get the following list: 
before <- apply(adtCountsCtrlRed, 2, function(x){
  apply(adtCountsCtrlRed, 2, cor, y = x, method = "spearman")
})

afterTotalVi <- apply(normAdtCountsCtrlRed, 2, function(x){
  apply(normAdtCountsCtrlRed, 2, cor, y = x, method = "spearman")
})

afterCorr <- apply(transformed_xRed, 2, function(x){
  apply(transformed_xRed, 2, cor, y = x, method = "spearman")
})

#If we now look at the correlation difference, it looks like this: 
corDif <- afterTotalVi-before
hist(as.vector(corDif), breaks = 50)
summary(as.vector(corDif))

#So these correlations have gone up very significantly. The top five in this list are: 
corDifOrd <- order(corDif, decreasing = TRUE)
topCorDif <- corDif[corDifOrd[1:5]]
topCols <- apply(corDif, 2, function(x) any(x %in% topCorDif))
topRows <- apply(corDif, 1, function(x) any(x %in% topCorDif))

colnames(corDif)[topCols]
#"CD003_PROT"  "CD011b_PROT" "CD027_PROT"  "CD127_PROT" 

#From the perspective of us already having applied the correction, and even more
#interesting metric might be: which are the ones where the summed effect of the
#correction makes the largest difference: 
correctionEffect <- corDif+(afterTotalVi-afterCorr)

#And the max candidates here are: 
correctionEffectOrd <- order(correctionEffect, decreasing = TRUE)
topCorrectionEffect <- correctionEffect[correctionEffectOrd[1:5]]
topCols <- apply(correctionEffect, 2, function(x) any(x %in% topCorrectionEffect))
colnames(correctionEffect)[topCols]
#"CD003_PROT"  "CD011b_PROT" "CD027_PROT"  "CD127_PROT" 
#Whith is the same. So now, these are investigated. 
plot(adtCountsCtrl[1:10000,"CD003_PROT"], adtCountsCtrl[1:10000,"CD011b_PROT"])
plot(normAdtCountsCtrl[1:10000,"CD003_PROT"], normAdtCountsCtrl[1:10000,"CD011b_PROT"])
plot(transformed_x[1:10000,"CD003_PROT"], transformed_x[1:10000,"CD011b_PROT"])

plot(adtCountsCtrl[1:10000,"CD003_PROT"], adtCountsCtrl[1:10000,"CD027_PROT"])
plot(normAdtCountsCtrl[1:10000,"CD003_PROT"], normAdtCountsCtrl[1:10000,"CD027_PROT"])
plot(transformed_x[1:10000,"CD003_PROT"], transformed_x[1:10000,"CD027_PROT"])

plot(adtCountsCtrl[1:10000,"CD003_PROT"], adtCountsCtrl[1:10000,"CD127_PROT"])
plot(normAdtCountsCtrl[1:10000,"CD003_PROT"], normAdtCountsCtrl[1:10000,"CD127_PROT"])
plot(transformed_x[1:10000,"CD003_PROT"], transformed_x[1:10000,"CD127_PROT"])

plot(adtCountsCtrl[1:10000,"CD011b_PROT"], adtCountsCtrl[1:10000,"CD027_PROT"])
plot(normAdtCountsCtrl[1:10000,"CD011b_PROT"], normAdtCountsCtrl[1:10000,"CD027_PROT"])
plot(transformed_x[1:10000,"CD011b_PROT"], transformed_x[1:10000,"CD027_PROT"])

#THis is a candidate. 
plot(adtCountsCtrl[1:10000,"CD011b_PROT"], adtCountsCtrl[1:10000,"CD127_PROT"])
cor(adtCountsCtrl[,"CD011b_PROT"], adtCountsCtrl[,"CD127_PROT"], method = "spearman")
plot(normAdtCountsCtrl[1:10000,"CD011b_PROT"], normAdtCountsCtrl[1:10000,"CD127_PROT"])
cor(normAdtCountsCtrl[,"CD011b_PROT"], normAdtCountsCtrl[,"CD127_PROT"], method = "spearman")
plot(transformed_x[1:10000,"CD011b_PROT"], transformed_x[1:10000,"CD127_PROT"])
cor(transformed_x[,"CD011b_PROT"], transformed_x[,"CD127_PROT"], method = "spearman")


plot(adtCountsCtrl[1:10000,"CD027_PROT"], adtCountsCtrl[1:10000,"CD127_PROT"])
plot(normAdtCountsCtrl[1:10000,"CD027_PROT"], normAdtCountsCtrl[1:10000,"CD127_PROT"])
plot(transformed_x[1:10000,"CD027_PROT"], transformed_x[1:10000,"CD127_PROT"])

#So now, we go to real plotting. 
dir.create("Diagnostics/Estimation_controls/Supplementary_figure_plots")
dDensityPlot(xYData = adtCountsCtrl[,c("CD011b_PROT", "CD127_PROT")],
             plotName = "Diagnostics/Estimation_controls/Supplementary_figure_plots/1_pre_totalVI_11b_vs_127")
dDensityPlot(xYData = normAdtCountsCtrl[,c("CD011b_PROT", "CD127_PROT")],
             plotName = "Diagnostics/Estimation_controls/Supplementary_figure_plots/2_post_totalVI_11b_vs_127")
dDensityPlot(xYData = transformed_x[,c("CD011b_PROT", "CD127_PROT")],
             plotName = "Diagnostics/Estimation_controls/Supplementary_figure_plots/3_post_correction_11b_vs_127")


