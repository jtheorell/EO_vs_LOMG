
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
  enrichedUp <- enrichr(row.names(upDat), "ARCHS4_TFs_Coexp")
  Sys.sleep(5)
  downDat <- i[which(i$logFC < 0),]
  enrichedDown <- enrichr(row.names(downDat), "ARCHS4_TFs_Coexp")
  Sys.sleep(5)
  list("Up" = enrichedUp, "Down" = enrichedDown)
})

#Now, the same is done for all the random populations for the 
edgerRandomOutcomes <- readRDS("../External/Stockholm/Data/RandomPopEdgeR/All_together.rds")

#Here, we take any term associated to anything. It turns out that the query 
#is only done about once every two seconds to the database, so when the algorithm
#is run in parallel mode, or for that matter in sequential mode, but without 
#sleep, it returns faulty results, that are just long ranges of exactly the same. 
#hence this algorithm, where one second of sleep is sufficient on my computer. 

randomGeneList <- list()
for(i in seq_along(edgerRandomOutcomes)){
  print(i)
  x <- edgerRandomOutcomes[[i]]
  randomGeneList[[i]] <- lapply(x, function(y){
    print(length(which(y$logFC > 0)))
    upDat <- y[which(y$logFC > 0),]
    if(nrow(upDat) > 0){
      enrichedUp <- enrichr(row.names(upDat), "ARCHS4_TFs_Coexp")
      longGeneNamesUp <- enrichedUp$ARCHS4_TFs_Coexp$Term[which(enrichedUp$ARCHS4_TFs_Coexp$Adjusted.P.value < 0.05)]
      Sys.sleep(1)      
    } else {
      longGeneNamesUp <- NA
    }
    downDat <- y[which(y$logFC < 0),]
    if(nrow(downDat) > 0){
      enrichedDown <- enrichr(row.names(downDat), "ARCHS4_TFs_Coexp")
      longGeneNamesDown <- enrichedDown$ARCHS4_TFs_Coexp$Term[which(enrichedUp$ARCHS4_TFs_Coexp$Adjusted.P.value < 0.05)]
      Sys.sleep(1)
    } else {
      longGeneNamesDown <- NA
    }
    list("Up" = longGeneNamesUp, "Down" = longGeneNamesDown)
  })
}
#Takes about 40 minutes to run on my computer, as a sleep break needs to be included. 
#Otherwise the function is lazy enough to just retrieve the results from the previous
#round. 

saveRDS(randomGeneList, "Results/Data/Stockholm/Random_gene_list_TFs.rds")
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
    innerSignGeneNames <- innerSign$ARCHS4_TFs_Coexp$Term
    innerSign$ARCHS4_TFs_Coexp$Gene_names <- innerSignGeneNames
    #Here, we count the number of genes associated to the TF in question, and further
    #calculate what fraction of the total number of genes that this corresponds to
    nAssGenes <- str_count(innerSign$ARCHS4_TFs_Coexp$Genes, ";")+1
    innerSign$ARCHS4_TFs_Coexp$fracAssGenes <- nAssGenes/numOfSignGenesList[[y]]
    innerRandom <- locRandomList[[y]]
    if(any(innerSignGeneNames %in% row.names(innerRandom))){
      innerSignNonRandom <- innerSign$ARCHS4_TFs_Coexp[-which(innerSignGeneNames %in% row.names(innerRandom)),]
    } else {
      innerSignNonRandom <- innerSign$ARCHS4_TFs_Coexp
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
dir.create("Results/Graphics/Stockholm/EnrichR_TF")
for(i in names(filteredGeneList)){
  for(j in c("Up", "Down")){
    locDat <- filteredGeneList[[i]][[j]]
    locColor <- ifelse(j == "Up", "#901E76", "#EDB953")
    #We will arrange these so that they go in different directions depending 
    #on if they are associated with up- or downregulated genes. 
    
    #locDat$Names <- gsub("| MP:.......", "", locDat$Gene_names)
    locDat$Names <- locDat$Gene_names
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
    ggsave(paste0("Results/Graphics/Stockholm/EnrichR_TF/", i, "_", j, ".pdf"))
    p <- p+theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                 axis.text.y = element_blank(), axis.title.y = element_blank())
    ggsave(paste0("Results/Graphics/Stockholm/EnrichR_TF/", i, "_", j, "_stripped.pdf"), width = 5.5, height = 5)
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
#[1] "ZNF831 human tf ARCHS4 coexpression" "SCML4 human tf ARCHS4 coexpression"  "BCL11B human tf ARCHS4 coexpression"
#[4] "ETS1 human tf ARCHS4 coexpression" 

geneDf
#                    CD8T_LOMGlow_1_Up                 CD8T_LOMGlow_1_Down                    CD8T_LOMGlow_2_Up
#1 BCL11B human tf ARCHS4 coexpression  DRAP1 human tf ARCHS4 coexpression TSC22D3 human tf ARCHS4 coexpression
#2   TCF7 human tf ARCHS4 coexpression  NR1H2 human tf ARCHS4 coexpression     FOS human tf ARCHS4 coexpression
#3   LEF1 human tf ARCHS4 coexpression  STAT4 human tf ARCHS4 coexpression   SCML4 human tf ARCHS4 coexpression
#4 ZNF831 human tf ARCHS4 coexpression ZNF524 human tf ARCHS4 coexpression    RORA human tf ARCHS4 coexpression
#5   ETS1 human tf ARCHS4 coexpression  IFI16 human tf ARCHS4 coexpression ZFP36L1 human tf ARCHS4 coexpression
#                  CD8T_LOMGlow_2_Down                     NK_EOMGlow_1_Up                   NK_EOMGlow_1_Down
#1 ZNF831 human tf ARCHS4 coexpression RNF166 human tf ARCHS4 coexpression NFATC2 human tf ARCHS4 coexpression
#2  SCML4 human tf ARCHS4 coexpression  ASCL2 human tf ARCHS4 coexpression THAP11 human tf ARCHS4 coexpression
#3  IKZF3 human tf ARCHS4 coexpression   KLF2 human tf ARCHS4 coexpression  ZFP42 human tf ARCHS4 coexpression
#4 BCL11B human tf ARCHS4 coexpression ZNF710 human tf ARCHS4 coexpression   ZIC2 human tf ARCHS4 coexpression
#5   ETS1 human tf ARCHS4 coexpression ZNF276 human tf ARCHS4 coexpression  TBX21 human tf ARCHS4 coexpression

#So all the NK genes are unique, which is reassuring. 