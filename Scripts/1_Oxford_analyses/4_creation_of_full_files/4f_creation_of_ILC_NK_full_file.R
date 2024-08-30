#We start by importing the files, that after export have gone through a round 
#of NameChanger name changes, to make sure that the timing (if PBMC pre if before thymectomy
#and post if after) or cell type is clear. Importantly, for the JR numbered 
#controls, there are two
#formats in the original names: JRXXXX_XXXX, and JRXXXX. The former is, very 
#counterintuitively, showing the number of the individual sample first and then
#the number of the patient. For this reason, I have removed the first number
#in nameChanger, before going into this analysis. 
#The cells are importantly exported as channel values, as this corresponds to the
#transformed values. 
library(FNN)
filenames <- list.files(path = "../External/Oxford/Exported_csv/ILC_NK_pre-gated", 
                        full.names = TRUE)

big_all_file <- do.call("rbind", lapply(filenames, function(x){
  ret <- read.csv(x)
  id <- gsub(pattern="../External/Oxford/Exported_csv/ILC_NK_pre-gated/ILCNK_|_.+\\.csv", "\\1", x)
  ret$id <- id
  tissue <- gsub(pattern="../External/Oxford/Exported_csv/ILC_NK_pre-gated/ILCNK_.+_|\\.csv", "\\1", x)
  ret$tissue <- tissue 
  ret 
}))

#Now, we start adding in metadata. 
#There are only three samples that are post-thymectomy, but we know how long 
#post thymectomy they are, so we are adding that data here:
big_all_file$postTxDays <- NA
big_all_file$postTxDays[which(big_all_file$id == "EJon" &
                                big_all_file$tissue == "post")] <- 12
big_all_file$postTxDays[which(big_all_file$id == "EHal" &
                                big_all_file$tissue == "post")] <- 51
big_all_file$postTxDays[which(big_all_file$id == "PJor" &
                                big_all_file$tissue == "post")] <- 32
table(big_all_file$postTxDays, big_all_file$tissue, useNA = "ifany")
#          post     pre     thy
#  12     98950       0       0
#  32     19261       0       0
#  51     42470       0       0
#  <NA>       0 6440786   71137

#All other metadata is not dependent on the tissue, and can thus be imported on a
#per-id basis. 
metadata <- read.csv("Data/Oxford/Documentation/Metadata_per_id.csv")

big_all_file$group <- big_all_file$age <- 
  big_all_file$sex <- big_all_file$hyperplasia <- big_all_file$sampleYear <- NA

for(i in metadata$Label){
  print(i)
  locMetRow <- which(metadata$Label == i)
  locDatRows <- which(big_all_file$id == i)
  big_all_file$group[locDatRows] <- metadata$Group[locMetRow]
  big_all_file$age[locDatRows] <- metadata$Age[locMetRow]
  big_all_file$sex[locDatRows] <- metadata$Sex[locMetRow]
  big_all_file$hyperplasia[locDatRows] <- metadata$Thymic_hyperplasia[locMetRow]
  big_all_file$sampleYear[locDatRows] <- metadata$Sample_year[locMetRow]
}

#Now, checking that it worked: 
table(big_all_file$group, useNA = "always")
#   Ctrl    EOMG    LOMG preLOMG    <NA> 
#2955745 1570421 1294682  851756       0 
table(big_all_file$sex, useNA = "always")
#      F       M    <NA> 
#3113273 3559331       0 

#We also import the lymphocyte counts that will be used downstream for
#normalisation purposes
big_all_file$idTime <- paste0(big_all_file$id, "_", big_all_file$tissue)

lymphCounts <- read.csv("Data/Oxford/Lymphocyte_counts/ILC_NK_lymphocyte_counts.csv")
lymphCounts$idTime <- paste0(lymphCounts$Label, "_", lymphCounts$Tissue)

identical(sort(unique(lymphCounts$idTime)), sort(unique(big_all_file$idTime)))
#TRUE
big_all_file$lymphCount <- NA
for(i in unique(big_all_file$idTime)){
  big_all_file$lymphCount[which(big_all_file$idTime == i)] <- 
    lymphCounts$Lymphocyte_count[which(lymphCounts$idTime == i)]
}
table(big_all_file$lymphCount, useNA = "always")
#No NA left. 

#Now, it turns out that it is very hard to work with this massive amount of 
#data, and further, it makes any analyses of individual ILC subsets impossible, 
#as they drown in the NK cell populations. Therefore, we will gate out the
#bona fide NK cells and ILC from this pool of mixed cells.

#ILC will in this context be defined as CD127+NKG2A-NKG2C-. It is possible
#That this leads to the inclusion of some atypical cells, especially for the 
#thymic samples, but we will have to live with that. We will use Euclidean
#gates to make sure that we catch all the relevant cells. 
testSamp <- sample(1:nrow(big_all_file), 10000)
hist(big_all_file$NKG2A[testSamp], breaks = 50)
abline(v = 300, col = "red")
hist(big_all_file$NKG2C[testSamp], breaks = 50)
abline(v = 200, col = "red")
hist(big_all_file$CD127[testSamp], breaks = 50)
abline(v = 500, col = "red")


big_all_file$cellType <- NA
big_all_file$cellType[which(big_all_file$NKG2A < 300 &
                              big_all_file$NKG2C < 200 &
                              big_all_file$CD127 > 500)] <- "ILC"
hist(big_all_file$CD56[testSamp], breaks = 50)
abline(v = 300, col = "red")
hist(big_all_file$CD16[testSamp], breaks = 50)
abline(v = 400, col = "red")
hist(big_all_file$CD7[testSamp], breaks = 50)
abline(v = 400, col = "red")
hist(big_all_file$NKG2A[testSamp], breaks = 50)
abline(v = 400, col = "red")

#For the NK cells
big_all_file$cellType[which(is.na(big_all_file$cellType) &
                              big_all_file$CD56 > 300 |
                              (big_all_file$CD16 > 400 &
                              big_all_file$CD7 > 400))] <- "NK"

#We will also take this opportunity to divide the CD56bright and dim, just to avoid
#categorising any borderline ILC/CD56bright into the wrong category. 
big_all_file$cellType[which(big_all_file$cellType == "NK" &
                              big_all_file$CD56 > 720 )] <- "CD56brightNK"

#And finally, we define a definite dummy population, negative for all the lineage
big_all_file$cellType[which(is.na(big_all_file$cellType) &
                              big_all_file$CD7 < 200 &
                              big_all_file$CD16 < 200 &
                              big_all_file$CD56 < 100 &
                              big_all_file$CD127 < 300)] <- "Something_else"



table(big_all_file$cellType, useNA = "always")

#  CD56brightNK            ILC             NK Something_else           <NA> 
#        175650          28990        6313495              5         154464 

#Now, we will investigate the remaining cells, to make sure that that
#fraction do not contain any obvious canditates for any of the cell types. 
naFrac <- sample(which(is.na(big_all_file$cellType)), 10000)
plot(big_all_file$CD7[naFrac],
     big_all_file$CD56[naFrac])
plot(big_all_file$CD7[naFrac],
     big_all_file$CD16[naFrac])
plot(big_all_file$CD7[naFrac],
     big_all_file$CD57[naFrac])

#It does not seem that way. 

#So, that means that it seems like we have a reasonably clean setup. Then
#we will use these four gates to define the remaining cells. We will for 
#this purpose use a reduced set of markers: 
modelDf <- big_all_file[,c("NKp30","HLADR","CD183","CD161",
                           "CD16","CD34","CD56","CD8",
                           "CD127","CD27","CD57","NKG2A",
                           "CD2","NKG2C","CRTH2",
                           "CD7","CD117","CD218a")]

#Now, we pre-cluster everything, to make things go a bit faster. 
realCellTypes <- unique(big_all_file$cellType)
realCellTypes <- realCellTypes[-which(is.na(realCellTypes))]
centroidMat <- do.call("rbind", lapply(realCellTypes, function(x){
  print(x)
  locRows <- which(big_all_file$cellType == x)
  k <- max(1, round(sqrt(length(locRows)/100)))
  set.seed(1111)
  locRes <- kmeans(modelDf[locRows,], k)
  clustCentroids <- locRes$centers
  rownames(clustCentroids) <- 
    paste0(x, "_", rownames(clustCentroids))
  clustCentroids
}))

#And now, each of the remaining cells are associated to one of those clusters.
set.seed(111)
cellIndexForEach <- 
  knnx.index(centroidMat,
             modelDf[which(is.na(big_all_file$cellType)),], 
             k = 1)

#Now, just to make sure that we get the numbers right, we will associate
#each cluster to its own cluster number. 
clusterIndex <- knnx.index(centroidMat, centroidMat, k = 1)
length(unique(clusterIndex[,1])) #300, as expected, and they are all 
#correctly numbered: if plotted, a diagonal line is created. 
#This is now used to create a key used for the other data: 
row.names(clusterIndex) <- row.names(centroidMat)
cellTypeForEach <- rep(NA, length(cellIndexForEach))
for(i in unique(cellIndexForEach)){
  cellTypeForEach[which(cellIndexForEach == i)] <- 
    row.names(clusterIndex)[which(clusterIndex[,1] == i)]
}
table(cellTypeForEach, useNA = "always")
#There are no NA left. 
length(unique(cellTypeForEach)) #177

#And with that, we are ready to collapse these back to the original groups
cellTypeForEachCollapsed <- gsub("|_.+", "", cellTypeForEach)

#And now, we can reintroduce this data into the form.
big_all_file$cellType[which(is.na(big_all_file$cellType))] <- 
  cellTypeForEachCollapsed

#Unintededly, all the ones previously called Something_else are now called
#Something as there is an underscore in between. 
big_all_file$cellType[which(big_all_file$cellType == "Something")] <- "Something_else"
table(big_all_file$cellType, useNA = "always")
          
#  CD56brightNK            ILC             NK Something_else           <NA> 
#        175664          73434        6421875           1631              0 

#So this landed nicely. We will of course collapse the two NK populations 
#downstream, but we will break out the ILC. 
#However, as the ILC population has grown by 300%, we will make sure that it has
#not been contaminated by any significant populations. 
ilcRows <- which(big_all_file$cellType == "ILC")
plot(big_all_file$NKG2A[ilcRows], big_all_file$CD127[ilcRows])
plot(big_all_file$CD56[ilcRows], big_all_file$NKG2A[ilcRows])
plot(big_all_file$CD16[ilcRows], big_all_file$NKG2A[ilcRows])
plot(big_all_file$CD7[ilcRows], big_all_file$NKG2A[ilcRows])
#There is astrange, small fraction of these cells that are NKG2A+, CD127+ and CD7+
#but negative for CD56 and low/negative for CD16. We will for this analysis include
#them rather in the NK cell population, as that is such a big population that a few
#cells of uncertain population status will not harm the analysis, whereas
#fheir faulty inclusion in the ILC analysis might screw things up, given that the
#total number of cells here is only 1% of the NK cell population. 
big_all_file$cellType[which(big_all_file$cellType == "ILC" &
                              big_all_file$NKG2A > 500)] <- "NK"
ilcRows <- which(big_all_file$cellType == "ILC")
plot(big_all_file$NKG2A[ilcRows], big_all_file$CD127[ilcRows])
plot(big_all_file$CD56[ilcRows], big_all_file$CD127[ilcRows])
plot(big_all_file$CD16[ilcRows], big_all_file$CD127[ilcRows])
plot(big_all_file$CD161[ilcRows], big_all_file$CD127[ilcRows])
plot(big_all_file$CD161[ilcRows], big_all_file$CD127[ilcRows])
plot(big_all_file$CRTH2[ilcRows], big_all_file$CD127[ilcRows])
plot(big_all_file$CD117[ilcRows], big_all_file$CD127[ilcRows])
#There is further a small population of cells that have been included here 
#that do not express CD127 and that lack all ILC markers apart from CD117. 
#They are for definition reasons also excluded. 
hist(big_all_file$CD127[ilcRows], breaks = 50)
abline(v = 280, col = "red")
big_all_file$cellType[which(big_all_file$cellType == "ILC" &
                              big_all_file$CD127 < 280)] <- "Something_else"

ilcRows <- which(big_all_file$cellType == "ILC")
plot(big_all_file$NKG2A[ilcRows], big_all_file$CD127[ilcRows])

#One more subpopulation that has been previously shown not to be an ILC, but rather 
#a non-surface expressing T-cell population is CD127+CD8+ cells. These will also
#be excluded. 
plot(big_all_file$CD8[ilcRows], big_all_file$CD127[ilcRows])
line_eq <- function(x) 7 * x - 2500

# Add the line to the plot
curve(line_eq, add = TRUE, col = "red")

# Identify points on one side of the line
points_on_one_side <- which(big_all_file$CD127[ilcRows] > 
                              line_eq(big_all_file$CD8[ilcRows]))

# Highlight points on one side of the line
points(big_all_file$CD8[ilcRows][points_on_one_side], 
       big_all_file$CD127[ilcRows][points_on_one_side], col = "blue")

#And we are ready for the exclusion
big_all_file$cellType[which(big_all_file$cellType == "ILC" &
                              big_all_file$CD127 < 
                              line_eq(big_all_file$CD8))] <- "Something_else"

ilcRows <- which(big_all_file$cellType == "ILC")
plot(big_all_file$CD8[ilcRows], big_all_file$CD127[ilcRows])

#Now, we investigate the NK cells. 
brigthRows <- which(big_all_file$cellType == "CD56brightNK")
plot(big_all_file$CD56[brigthRows], big_all_file$CD16[brigthRows])
#These all look like bona fide CD56bright NK cells
dimRows <- sample(which(big_all_file$cellType == "NK"), 10000)
plot(big_all_file$CD56[dimRows], big_all_file$CD16[dimRows])
plot(big_all_file$CD56[dimRows], big_all_file$CD7[dimRows])
#This shows that almost all the cells are CD56 positive, and the few that are not
#are CD7+CD16+. So we are happy here. 
table(big_all_file$cellType, useNA = "always")
#  CD56brightNK            ILC             NK Something_else           <NA> 
#        175664          62804        6422932          11204              0 

#Beautiful. 

#And with that, we are ready to save the file. 
dir.create("../External/Oxford/Resulting_data/ILC_NK")
saveRDS(big_all_file, "../External/Oxford/Resulting_data/ILC_NK/ILC_NK_full_file.rds")

#Now, we focus on the separate compartments. 
ilcDf <- big_all_file[which(big_all_file$cellType == "ILC"),]
dir.create("../External/Oxford/Resulting_data/ILC")
saveRDS(ilcDf, "../External/Oxford/Resulting_data/ILC/ILC_full_file.rds")

nkDf <- big_all_file[grep("NK", big_all_file$cellType),]
dir.create("../External/Oxford/Resulting_data/NK")
saveRDS(nkDf, "../External/Oxford/Resulting_data/NK/NK_full_file.rds")


