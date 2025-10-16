library(fmsb)

#First, we import a manually annotated file for this purpose.
disc_clust_char <- read.csv("Data/Discriminant_cluster_characteristics.csv", row.names = 1)
#This throws a warning. Not to worry. 

#First, we need to normalise the values. E.g. for the MAITs in UK, we need to relate
#the sizes of this cluster from EOMG, LOMG and relevant controls
updatedDat <- disc_clust_char
updatedDat[,grep("in_..MG_UK|ctrls", colnames(updatedDat))] <- 
  t(apply(disc_clust_char[,grep("in_..MG_UK|ctrls", colnames(updatedDat))], 1, function(x){
    x/max(x)
  }))
updatedDat[,grep("in_..MG_SE", colnames(updatedDat))] <- 
  t(apply(disc_clust_char[,grep("in_..MG_SE", colnames(updatedDat))], 1, function(x){
    x/max(x)
  }))

#We need to add max and min values. 

fullData <- rbind(c(rep(16,2),rep(1,5),56), 
                  rep(0, ncol(updatedDat)),
                  updatedDat)
write.csv(fullData, "Results/Data/For_figure_file/Figure_4G_radar_graph.csv")
dirName <- "Results/Graphics/Oxford_and_Stockholm/Discriminant_cluster_characteristics/"
dir.create(dirName)
pdf(paste0(dirName, "Discriminant_cluster_characteristics.pdf"), width = 10, height = 10)
radarchartcirc(fullData, pfcol = c("#34158855", "#FFBF0088", "#88888888"),
               pcol = c("#341588", "#FFBF0077", "#888888"))
dev.off()

