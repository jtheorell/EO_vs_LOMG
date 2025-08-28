library(SingleCellExperiment)
popVec <- c("B", "CD4T", "CD8T", "TCRgd", "NK", "ILC")

for(i in popVec){
  print(i)
  datDir <- paste0("../External/Oxford/Resulting_data/", i, "/")
  locFile <- readRDS(paste0(datDir, i, "_full_file_post_Euclid_consensus.rds"))
  print(length(which(locFile$SmoothGroup != "None"))/nrow(locFile))
}
#[1] "B"
#[1] 0.07505516
#[1] "CD4T"
#[1] 0.07671225
#[1] "CD8T"
#[1] 0.2006817
#[1] "TCRgd"
#[1] 0.2989536
#[1] "NK"
#[1] 0.1023844
#[1] "ILC"
#[1] 0.04986927

#Now, instead, we identify the numbers of cells per group
groupMat <- do.call("rbind", lapply(popVec, function(x){
  print(x)
  datDir <- paste0("../External/Oxford/Resulting_data/", x, "/")
  locFile <- readRDS(paste0(datDir, x, "_full_file_post_Euclid_consensus.rds"))
  locGroup <- locFile$group
  locGroup[which(locGroup == "Ctrl")] <- "youngCtrl"
  locGroup[which(locGroup == "youngCtrl" & locFile$age > 50)] <- "oldCtrl"
  table(locGroup)
}))
colSums(groupMat)
#    Ctrl     EOMG     LOMG  unclear 
#10011513  8474790  3280095  2058739
sum(groupMat)
#23825137

#Now, furthermore, how many are the cells in each of the three hit clusters and
#what percentage do they make up?
CD8T <- readRDS(paste0("../External/Oxford/Resulting_data/CD8T/CD8T_full_file_post_Euclid_including_subclusts.rds"))
signOxStock <- read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1)
numCellsCD8T_LOMGlow_1 <- length(which(CD8T$smoothGroupDeep == gsub("CD8T_|", "", row.names(signOxStock)[1])))
#365185
numCellsCD8T_LOMGlow_2 <- length(which(CD8T$smoothGroupDeep == gsub("CD8T_|", "", row.names(signOxStock)[2])))
#128071
fracOfAllCD8_LOMGlow_1 <- numCellsCD8T_LOMGlow_1/nrow(CD8T)
#0.08929242
fracOfAllCD8_LOMGlow_2 <- numCellsCD8T_LOMGlow_2/nrow(CD8T)
#0.031315

NK <- readRDS(paste0("../External/Oxford/Resulting_data/NK/NK_full_file_post_Euclid_including_subclusts.rds"))

numCellsNK_EOMGlow_1 <- length(which(NK$smoothGroupDeep == gsub("NK_|", "", row.names(signOxStock)[3])))
#395043
fracOfAllNK_EOMGlow_1 <- numCellsNK_EOMGlow_1/nrow(NK)
#0.06007305

#And now for the Stockholm cohort: 

stockDatCD8T <- readRDS("../External/Oxford_and_Stockholm/Harmonisation/CD8T_Stockholm_data_with_Ox_neighbors.rds")
numCellsCD8T_LOMGlow_1_Stock <- length(which(stockDatCD8T$Ox_hit_clust == row.names(signOxStock)[1]))
#1242
fracOfAllCD8_LOMGlow_1_Stock <- numCellsCD8T_LOMGlow_1_Stock/nrow(stockDatCD8T)
#0.1680877
numCellsCD8T_LOMGlow_2_Stock <- length(which(stockDatCD8T$Ox_hit_clust == row.names(signOxStock)[2]))
#146
fracOfAllCD8_LOMGlow_2_Stock <- numCellsCD8T_LOMGlow_2_Stock/nrow(stockDatCD8T)
#0.0197591

stockDatNK <- readRDS("../External/Oxford_and_Stockholm/Harmonisation/NK_Stockholm_data_with_Ox_neighbors.rds")
numCells_NK_EOMGlow_1_Stock <- length(which(stockDatNK$Ox_hit_clust == row.names(signOxStock)[3]))
#316
fracOfAllNK_EOMGlow_1_Stock <- numCells_NK_EOMGlow_1_Stock/nrow(stockDatNK)
#0.08499193


#Ans if all the cell types are analysed?

sthlmProtDat <- readRDS(paste0("../External/Stockholm/Data/SCE_versions/Stockholm_protein_data_including_colData.rds"))
table(sthlmProtDat$MG_type)
#     EOMG      LOMG uncertain 
#    17625     12294      5815
ncol(sthlmProtDat)
#35734

