focPops <- read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv")
allPercentPops <- read.csv("Results/Data/Oxford/Analysis_post_smooth/freqTableAllPops.csv", row.names = 1)

focPercentPops <- allPercentPops[,which(colnames(allPercentPops) %in% focPops$X)]

#Now, the specific question here concerns the LOMGs and the fact that they
#are quite spread in terms of their NK cell percentage. 
LOMGPercentPops <- focPercentPops[which(allPercentPops$Group == "LOMG"),]

apply(LOMGPercentPops, 2, function(x){
  apply(LOMGPercentPops, 2, function(y){
    cor(x, y, method = "spearman")
  })
})

#               CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#CD8T_LOMGlow_1     1.00000000      0.6004902  -0.09313725
#CD8T_LOMGlow_2     0.60049020      1.0000000   0.16666667
#NK_EOMGlow_1      -0.09313725      0.1666667   1.00000000

#Now, we calculate the median percentages for all groups in this dataset
apply(focPercentPops, 2, function(x){
  splitDat <- split(x, allPercentPops$Group)
  unlist(lapply(splitDat, median))
  
})

#           CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#EOMG         0.0091476846   0.0069115044  0.000920339
#LOMG         0.0005861561   0.0004399757  0.004168808
#Old_ctrl     0.0052346776   0.0035998972  0.009191378
#unclear      0.0097426471   0.0030162162  0.007061940
#Young_ctrl   0.0279417993   0.0038808156  0.004413771
#THis is directly transferred to the file used to construct the final figure

#Now, we do the same for the Stockholm data. 

stockDat <- read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/freqTable_Stockholm.csv")

focStockDat <- stockDat[,which(colnames(stockDat) %in% focPops$X)]
apply(focStockDat, 2, function(x){
  splitDat <- split(x, stockDat$Group)
  unlist(lapply(splitDat, median))
  
})

#          CD8T_LOMGlow_1 CD8T_LOMGlow_2 NK_EOMGlow_1
#EOMG          0.03938030   0.0042678878  0.003231110
#LOMG          0.01544051   0.0008957187  0.023078709
#uncertain     0.05752310   0.0037397083  0.002868537

