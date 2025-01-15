
oxPVals <- read.csv("Results_per_celltype/Data/Oxford/Analysis_post_smooth/pValTable.csv", row.names = 1)

oxSign <- oxPVals[which((oxPVals[,1] < 0.05 | oxPVals[,2] < 0.05) &
                          oxPVals[,3] < 0.05),]
stockSign <- read.csv("Results_per_celltype/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stockholm.csv", row.names = 1)

signOxStockNames <- Reduce(intersect, list(row.names(oxSign), row.names(stockSign)))


signOxStock <- stockSign[which(row.names(stockSign) %in% signOxStockNames),]
names(signOxStock) <- row.names(stockSign)[which(row.names(stockSign) %in% signOxStockNames)]

write.csv(signOxStock, "Results_per_celltype/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv")