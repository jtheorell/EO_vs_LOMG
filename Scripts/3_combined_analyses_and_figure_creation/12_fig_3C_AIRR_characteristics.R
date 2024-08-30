library(SingleCellExperiment)
library(ggplot2)
library(viridis)
#Here, we investigate the BCR and TCR expression differences between
#the significant and non-significant populations. 
logSCE <- readRDS("../External/Stockholm/Data/SCE_versions/6_with_Ox_neigh_clust.rds")
significantClusters <- read.csv("../External/Stockholm/Data/SCE_versions/6_with_Ox_neigh_clust.csv", row.names = 1)

logSCE$signOxClust <- logSCE$oxNeighClust
logSCE$signOxClust[-which(logSCE$oxNeighClust %in% row.names(significantClusters))] <- "None"

#Now, we are going to make one B and one T cell analysis. 

#We start with the B-cells. 

bcrSce <- logSCE[,which(logSCE$BCR)]


barGraphDat <- colData(bcrSce)[,c("signOxClust","BCR_Clonal", "BCR_light_type", 
                                  "BCR_c_call", "BCR_vFamH", 
                               "BCR_All_mutations_HL")]
#For the categorical data, we start. 
#We do need to make a few changes here. 
barGraphDat$signOxClust <- factor(barGraphDat$signOxClust, levels = c("None", "B_LOMGlow_2"))
barGraphDat$Clonal <- factor(barGraphDat$BCR_Clonal, c(TRUE, FALSE))
barGraphDat$light_type <- factor(barGraphDat$BCR_light_type, levels = c("IGK", "IGL"))
#Here, we remove all non-IgG. 
barGraphDat$Isotype <- factor(barGraphDat$BCR_c_call, levels = c("IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2"))
barGraphDat$V_gene_fam <- factor(barGraphDat$BCR_vFamH, levels = c("IGHV1", "IGHV2", "IGHV3", "IGHV4", "IGHV5", "IGHV6", "IGHV7"))

colsAndNames <- list(
    list("Clonal", c("grey", "black")),
    list("V_gene_fam", gray.colors(12, start =1, end = 0)[3:9]),
    list("light_type", c("orange", "darkblue")),
    list("Isotype", c("#B1CFFF", "#AE6552", "#904952", "#722E52", "#541352", "#238A8D", "#106060"))
)
"#8D8D8D"
dir.create("Results/Graphics/Stockholm/CloneBarGraphs", recursive = TRUE)
for(i in colsAndNames){
    print(i[[1]])
    locCounts <- table(barGraphDat[,i[[1]]], barGraphDat$signOxClust)
    locDat <- as.data.frame(locCounts)
    colnames(locDat) <- c(i[[1]], "Cluster", "value")
    p <- ggplot(locDat,                         
                aes(x = Cluster,
                    y = value,
                    fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        scale_fill_manual(values=i[[2]], drop = FALSE)
    p
    ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/", i[[1]], "_no_legend.pdf"),
           width = 9, height = 6)
    #The stats are also included here. 
    print(round(fisher.test(locCounts)$p.value, 2))
}
#The non-adjusted Fisher p-values are 0.34, 0.33, 0.62, 0.69 

#Finally, we will look at number of mutations. Here, we will have to take a slightly
#different approach, given that we cannot use fisher tests. 

mutationTable <- do.call("rbind", lapply(unique(bcrSce$sample), function(x){
    locSignMuts <- colData(bcrSce)$BCR_All_mutations_HL[which(bcrSce$sample == x &
                                                                  bcrSce$oxNeighClust != "None")]
    locSignMuts <- locSignMuts[-which(is.na(locSignMuts))]
    if(length(locSignMuts) > 5){
        locOtherMuts <- colData(bcrSce)$BCR_All_mutations_HL[which(bcrSce$sample == x &
                                                                       bcrSce$oxNeighClust == "None")]
        locOtherMuts <- locOtherMuts[-which(is.na(locOtherMuts))]
        c("OxClust" = median(locSignMuts),
          "Other" = median(locOtherMuts))
    } else {
        c("OxClust" = NA,
          "Other" = NA)
    }
    
}))
row.names(mutationTable) <- unique(bcrSce$sample)

mutationTable <- mutationTable[-which(is.na(mutationTable[,1])),]
mutationTable
#    OxClust Other
#1_1    14.5    15
#1_8    28.0    21
#2_2    27.0    21
#2_4    24.0    18
#2_6    17.0    18

#And this is now tested. 
wilcox.test(mutationTable[,1], mutationTable[,2], paired = TRUE, exact = FALSE)
#p-value = 0.2785


#We also produce a figure for this purpose. 
locDat <- as.data.frame(table(barGraphDat$BCR_All_mutations_HL, barGraphDat$signOxClust))
locDat[,1] <- as.numeric(as.character(locDat[,1]))
colnames(locDat) <- c("Mutations", "Cluster", "value")
p <- ggplot(locDat,                         
            aes(x = Cluster,
                y = value,
                fill = locDat[,1])) + 
    geom_bar(stat = "identity",
             position = "fill") + theme_bw() +
    scale_fill_viridis(discrete = FALSE, 
                       option = "magma")
p
ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/Mutations.pdf"))
p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/Mutations_no_legend.pdf"),
       width = 9, height = 6)


#ANd now over to the T cells. Here, we will run one CD4 and one CD8 analysis. 
cd4Sce <- logSCE[,which(logSCE$TCR & logSCE$cellType == "CD4T")]


barGraphDat <- colData(cd4Sce)[,c("signOxClust","TCR_Clonal", "TCR_vFamB")]

#We do need to make a few changes here. 
barGraphDat$signOxClust <- factor(barGraphDat$signOxClust, levels = c("None", "CD4T_EOMGlow_3"))
barGraphDat$CD4T_Clonal <- factor(barGraphDat$TCR_Clonal, c(TRUE, FALSE))
barGraphDat$CD4T_V_gene_fam <- factor(barGraphDat$TCR_vFamB, levels = c("TRBV1", "TRBV2", "TRBV3", "TRBV4", "TRBV5", "TRBV6", "TRBV7", "TRBV8", "TRBV9"))

colsAndNames <- list(
    list("CD4T_Clonal", c("grey", "black")),
    list("CD4T_V_gene_fam", gray.colors(12, start =1, end = 0)[2:11])
)

for(i in colsAndNames){
    print(i[[1]])
    locCounts <- table(barGraphDat[,i[[1]]], barGraphDat$signOxClust)
    locDat <- as.data.frame(locCounts)
    colnames(locDat) <- c(i[[1]], "Cluster", "value")
    p <- ggplot(locDat,                         
                aes(x = Cluster,
                    y = value,
                    fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        scale_fill_manual(values=i[[2]], drop = FALSE)
    p
    ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/", i[[1]], "_no_legend.pdf"),
           width = 9, height = 6)
    #The stats are also included here. 
    print(round(fisher.test(locCounts, simulate.p.value = TRUE)$p.value, 2))
}

#The non-adjusted p-values are 0.02 and 0.58

#Now, CD8
cd8Sce <- logSCE[,which(logSCE$TCR & logSCE$cellType == "CD8T")]


barGraphDat <- colData(cd8Sce)[,c("signOxClust","TCR_Clonal", "TCR_vFamB")]

#We do need to make a few changes here. 
barGraphDat$signOxClust <- factor(barGraphDat$signOxClust, levels = c("None", "CD8T_LOMGlow_2"))
barGraphDat$CD8T_Clonal <- factor(barGraphDat$TCR_Clonal, c(TRUE, FALSE))
barGraphDat$CD8T_V_gene_fam <- factor(barGraphDat$TCR_vFamB, levels = c("TRBV1", "TRBV2", "TRBV3", "TRBV4", "TRBV5", "TRBV6", "TRBV7", "TRBV8", "TRBV9"))

colsAndNames <- list(
    list("CD8T_Clonal", c("grey", "black")),
    list("CD8T_V_gene_fam", gray.colors(12, start =1, end = 0)[2:11])
)

for(i in colsAndNames){
    print(i[[1]])
    locCounts <- table(barGraphDat[,i[[1]]], barGraphDat$signOxClust)
    locDat <- as.data.frame(locCounts)
    colnames(locDat) <- c(i[[1]], "Cluster", "value")
    p <- ggplot(locDat,                         
                aes(x = Cluster,
                    y = value,
                    fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        scale_fill_manual(values=i[[2]], drop = FALSE)
    p
    ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/Graphics/Stockholm/CloneBarGraphs/", i[[1]], "_no_legend.pdf"),
           width = 9, height = 6)
    #The stats are also included here. 
    print(round(fisher.test(locCounts, simulate.p.value = TRUE)$p.value, 2))
}
#The p-values are 0.36 and 0. 


#So we need to investigate the CD8T v-family finding further, checking which V families that are contributing
#to this difference for the CD8T. 
CD8VFamDistrP <- unlist(lapply(row.names(locCounts), function(x){
    focRow <- locCounts[which(row.names(locCounts) == x),]
    otherRows <- colSums(locCounts[-which(row.names(locCounts) == x),])
    fisher.test(rbind(focRow, otherRows))$p.value
}))

names(CD8VFamDistrP) <- row.names(locCounts)
round(p.adjust(CD8VFamDistrP, method = "fdr"), 5)
#  TRBV1   TRBV2   TRBV3   TRBV4   TRBV5   TRBV6   TRBV7   TRBV8   TRBV9 
#0.00008 0.56869 0.90165 0.31402 0.02335 0.00000 0.00000 1.00000 0.00285

#So this is really interesting. It is possible that it is the very significant
#aberrations in TCRVB6 and 7 that drive the rest, but it does not necessarily matter. 
#There are no TRBV8 genes here, hence the 1.  

