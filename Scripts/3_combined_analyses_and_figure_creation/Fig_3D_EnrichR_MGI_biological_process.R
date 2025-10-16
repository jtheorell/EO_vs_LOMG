library(enrichR)
library(BiocParallel)
library(ggplot2)
library(stringr)
#Here, we are going to make this as sensical as possible, trying to identify underlying transcription factor
#profiles identifiable with these up- and downregulated gene sets. We will also
#create a range of random gene sets for each cell type and filter out genes with
#a high likelihood of just representing the cell population more in general. 
#dbs <- listEnrichrDbs()

edgerSignOutcomes <- readRDS("Results/Data/Stockholm/EdgeR_significant_transcriptomes.rds")
enrichedGeneList <- lapply(edgerSignOutcomes, function(i){
  upDat <- i[which(i$logFC > 0),]
  enrichedUp <- enrichr(row.names(upDat), "MGI_Mammalian_Phenotype_Level_4_2024")
  Sys.sleep(5)
  downDat <- i[which(i$logFC < 0),]
  enrichedDown <- enrichr(row.names(downDat), "MGI_Mammalian_Phenotype_Level_4_2024")
  Sys.sleep(5)
  list("Up" = enrichedUp, "Down" = enrichedDown)
})

#Now, the same is done for all the random populations for the 
edgerRandomOutcomes <- readRDS("../External/Stockholm/Data/RandomPopEdgeR/All_together.rds")

#Here, we take any term associated to anything. Something strange happens when
#the enrichr function is put inside too many loops, so I avoid that here. 
randomGeneList <- lapply(edgerRandomOutcomes, function(x){
  locTFList <- lapply(x, function(y){
    print(length(which(y$logFC > 0)))
    upDat <- y[which(y$logFC > 0),]
    if(nrow(upDat) > 0){
      enrichedUp <- enrichr(row.names(upDat), "MGI_Mammalian_Phenotype_Level_4_2024")
      longGeneNamesUp <- enrichedUp$MGI_Mammalian_Phenotype_Level_4_2024$Term
      Sys.sleep(5)      
    } else {
      longGeneNamesUp <- NA
    }
    downDat <- y[which(y$logFC < 0),]
    if(nrow(downDat) > 0){
      enrichedDown <- enrichr(row.names(downDat), "MGI_Mammalian_Phenotype_Level_4_2024")
      longGeneNamesDown <- enrichedDown$MGI_Mammalian_Phenotype_Level_4_2024$Term
      Sys.sleep(5)
    } else {
      longGeneNamesDown <- NA
    }
    list("Up" = longGeneNamesUp, "Down" = longGeneNamesDown)
  })
})
#Takes about 40 minutes to run on my computer, as a sleep break needs to be included. 
#Otherwise the function is lazy enough to just retrieve the results from the previous
#round. 

saveRDS(randomGeneList, "Results/Data/Stockholm/Random_gene_list.rds")
#randomGeneList <- readRDS("Results/Data/Stockholm/Random_gene_list.rds")
#Now, we are going to order these genes depending on how common they are per cell type, per direction
#and then we can filter out the ones that are not of interest, i.e. that are just 
#generally over-represented when working with this data. First, we reorganise the
#data per dataset. 
randomGeneListPerCelltype <- lapply(1:3, function(x) lapply(randomGeneList, "[[", x))

randomGeneSummaries <- lapply(randomGeneListPerCelltype, function(x){
  #Now, we restructure and count all genes. 
  innerUpTable <- table(unlist(lapply(x, "[[", 1)))
  innerUpTableOrdered <- innerUpTable[order(innerUpTable, decreasing = TRUE)]
  innerDownTable <- table(unlist(lapply(x, "[[", 2)))
  innerDownTableOrdered <- innerDownTable[order(innerDownTable, decreasing = TRUE)]
  #Now, we introduce a filter here: any transcription factor that occurs in more than 5% of the
  #random runs should not be considered of relevance
  innerUpTableFinal <- innerUpTableOrdered[which(innerUpTableOrdered > 4)]
  innerDownTableFinal <- innerDownTableOrdered[which(innerDownTableOrdered > 4)]
  list("Up" = innerUpTableFinal, "Down" = innerDownTableFinal)
})
names(randomGeneSummaries) <- row.names(read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1))

#So now, we use this to filter out the random genes from the ones we have found: 
filteredGeneList <- lapply(names(enrichedGeneList), function(x){
  print(x)
  logSignEdgeRList <- edgerSignOutcomes[[x]]
  upN <- length(which(logSignEdgeRList$logFC > 0))
  downN <- length(which(logSignEdgeRList$logFC < 0))
  numOfSignGenesList <- list("Up" = upN, "Down" = downN)
  locSignEnricheRList <- enrichedGeneList[[x]]
  locRandomList <- randomGeneSummaries[[x]]
  locFilteredList <- lapply(c("Up", "Down"), function(y){
    innerSign <- locSignEnricheRList[[y]]
    innerSignGeneNames <- innerSign$MGI_Mammalian_Phenotype_Level_4_2024$Term
    innerSign$MGI_Mammalian_Phenotype_Level_4_2024$Gene_names <- innerSignGeneNames
    #Here, we count the number of genes associated to the TF in question, and further
    #calculate what fraction of the total number of genes that this corresponds to
    nAssGenes <- str_count(innerSign$MGI_Mammalian_Phenotype_Level_4_2024$Genes, ";")+1
    innerSign$MGI_Mammalian_Phenotype_Level_4_2024$fracAssGenes <- nAssGenes/numOfSignGenesList[[y]]
    innerRandom <- locRandomList[[y]]
    if(any(innerSignGeneNames %in% row.names(innerRandom))){
      innerSignNonRandom <- innerSign$MGI_Mammalian_Phenotype_Level_4_2024[-which(innerSignGeneNames %in% row.names(innerRandom)),]
    } else {
      innerSignNonRandom <- innerSign$MGI_Mammalian_Phenotype_Level_4_2024
    }
    #And here, we take the top 5 most regulated. In the scenario where there are 
    #more than five with the same level of significance, all are exported
    innerTrueSign <- innerSignNonRandom[c(1:5),]
  })
  names(locFilteredList) <- c("Up", "Down")
  locFilteredList
})
names(filteredGeneList) <- names(enrichedGeneList)

#Now, we plot these. First, we identify the max value for the fraction of associated genes
maxFrac <- max(unlist(lapply(filteredGeneList, function(x) lapply(x, function(z) z$fracAssGenes))))
dir.create("Results/Graphics/Stockholm/EnrichR_MGI")
for(i in names(filteredGeneList)){
  for(j in c("Up", "Down")){
    locDat <- filteredGeneList[[i]][[j]]
    locColor <- ifelse(j == "Up", "#901E76", "#EDB953")
    #We will arrange these so that they go in different directions depending 
    #on if they are associated with up- or downregulated genes. 
    
    locDat$Names <- gsub("| MP:.......", "", locDat$Gene_names)
    if(j == "Up"){
      nameOrder <- locDat$Names[order(locDat$fracAssGenes)]
    } else {
      nameOrder <- locDat$Names[order(locDat$fracAssGenes, decreasing = TRUE)]
    }
    locDat$Names_factor <- factor(locDat$Names, 
                                       levels = nameOrder)
    p <- ggplot(locDat, aes(x = Names_factor, y = fracAssGenes, label = Names_factor)) +
      geom_bar(stat = "identity", fill=locColor ) +
      scale_y_reverse(limits = c(maxFrac*1.05, 0), expand = c(0,0)) +
      theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
      #geom_text(angle = 90, hjust=1)
    p
    ggsave(paste0("Results/Graphics/Stockholm/EnrichR_MGI/", i, "_", j, ".pdf"))
    p <- p+theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                 axis.text.y = element_blank(), axis.title.y = element_blank())
    ggsave(paste0("Results/Graphics/Stockholm/EnrichR_MGI/", i, "_", j, "_stripped.pdf"), width = 5.5, height = 5)
    
    #And we save the data
    write.csv(locDat, paste0("Results/Data/For_figure_file/Figure_3D_MGI_terms_", i, "_", j, ".csv"))
  }
}
#And we also export the gene list to the same directory. 
geneDf <- data.frame(do.call("cbind", lapply(names(filteredGeneList), function(x){
  locRes <- do.call("cbind", lapply(c("Up", "Down"), function(y){
    filteredGeneList[[x]][[y]]$Gene_names
  }))
  locRes
})))
colnames(geneDf) <- paste0(rep(names(filteredGeneList), each = 2), "_", c("Up", "Down"))
#These are duplicated: 
as.matrix(geneDf)[which(duplicated(as.vector(as.matrix(geneDf))))]
#None, which is very reassuring. 

geneDf
#CD8T_LOMGlow_1_Up
#1                         Increased T Cell Proliferation MP:0005348
#2                             Decreased B-1a Cell Number MP:0008168
#3                      Abnormal Thymus Cortex Morphology MP:0002371
#4                        Abnormal T Cell Differentiation MP:0002145
#5 Decreased CD8-positive, Naive Alpha-Beta T Cell Number MP:0013435
#CD8T_LOMGlow_1_Down                                    CD8T_LOMGlow_2_Up
#1                                                 Abnormal T Cell Physiology MP:0002444 Abnormal Extraembryonic Tissue Physiology MP:0004264
#2                                              Abnormal Thymocyte Activation MP:0003850                            Spleen Atrophy MP:0003643
#3 Increased CD4-positive, CD25-positive, Alpha-Beta Regulatory T Cell Number MP:0010168           Abnormal B Cell Differentiation MP:0002144
#4                                          Decreased Interleukin-2 Secretion MP:0008688                      Decreased IgG1 Level MP:0008495
#5                                           Increased T-helper 1 Cell Number MP:0008086                                    Anemia MP:0001577
#CD8T_LOMGlow_2_Down                                          NK_EOMGlow_1_Up
#1 Abnormal Extraembryonic Tissue Physiology MP:0004264                     Absent Endometrial Glands MP:0009097
#2                            Spleen Atrophy MP:0003643                         Increased Embryo Size MP:0001699
#3           Abnormal B Cell Differentiation MP:0002144 Increased Heart Left Ventricle Wall Thickness MP:0021157
#4                Decreased Thymocyte Number MP:0000715                            Genetic Imprinting MP:0003121
#5   Increased Double-Positive T Cell Number MP:0005091    Abnormal Type IV Hypersensitivity Reaction MP:0002534
#NK_EOMGlow_1_Down
#1 Decreased Circulating Parathyroid Hormone Level MP:0002905
#2              Abnormal Bone Trabecula Morphology MP:0010867
#3                   Abnormal Outer Ear Morphology MP:0002177
#4                    Increased Plasma Cell Number MP:0008097
#5                 Arrested B Cell Differentiation MP:0001802