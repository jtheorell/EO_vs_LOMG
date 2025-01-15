library(enrichR)
library(BiocParallel)
library(ggplot2)
#Here, we are going to make this as sensical as possible, trying to identify underlying transcription factor
#profiles identifiable with these up- and downregulated gene sets. We will also
#create a range of random gene sets for each cell type and filter out genes with
#a high likelihood of just representing the cell population more in general. 
dbs <- listEnrichrDbs()



edgerSignOutcomes <- readRDS("Results/Data/Stockholm/EdgeR_significant_transcriptomes.rds")
enrichedGeneList <- lapply(edgerSignOutcomes, function(i){
  upDat <- i[which(i$logFC > 0),]
  downDat <- i[which(i$logFC < 0),]
  enrichedGenes <- lapply(list(upDat, downDat), function(x){
    enriched <- enrichr(row.names(x), "WikiPathways_2024_Human")
  })
  names(enrichedGenes) <- c("Up", "Down")
  enrichedGenes
})

#Now, the same is done for all the random populations for the 
edgerRandomOutcomes <- readRDS("../External/Stockholm/Data/RandomPopEdgeR/All_together.rds")

randomGeneList <- bplapply(edgerRandomOutcomes, function(x){
  locTFList <- lapply(x, function(y){
    upDat <- y[which(y$logFC > 0),]
    downDat <- y[which(y$logFC < 0),]
    enrichedGenes <- lapply(list(upDat, downDat), function(z){
      enriched <- enrichr(row.names(z), "WikiPathways_2024_Human")
      longGeneNames <- enriched$WikiPathways_2024_Human$Term[which(enriched$WikiPathways_2024_Human$Adjusted.P.value < 0.05)]
    })
    names(enrichedGenes) <- c("Up", "Down")
    enrichedGenes
  })
})

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
    innerSignGeneNames <- innerSign$WikiPathways_2024_Human$Term
    innerSign$WikiPathways_2024_Human$Gene_names <- innerSignGeneNames
    #Here, we count the number of genes associated to the TF in question, and further
    #calculate what fraction of the total number of genes that this corresponds to
    nAssGenes <- str_count(innerSign$WikiPathways_2024_Human$Genes, ";")+1
    innerSign$WikiPathways_2024_Human$fracAssGenes <- nAssGenes/numOfSignGenesList[[y]]
    innerRandom <- locRandomList[[y]]
    if(any(innerSignGeneNames %in% row.names(innerRandom))){
      innerSignNonRandom <- innerSign$WikiPathways_2024_Human[-which(innerSignGeneNames %in% row.names(innerRandom)),]
    } else {
      innerSignNonRandom <- innerSign$WikiPathways_2024_Human
    }
    #And here, we take the top 5 most regulated.
    innerTrueSign <- innerSignNonRandom[1:5,]
  })
  names(locFilteredList) <- c("Up", "Down")
  locFilteredList
})
names(filteredGeneList) <- names(enrichedGeneList)

#Now, we plot these. First, we identify the max value for the fraction of associated genes
maxFrac <- max(unlist(lapply(filteredGeneList, function(x) lapply(x, function(z) z$fracAssGenes))))
dir.create("Results/Graphics/Stockholm/EnrichR")
for(i in names(filteredGeneList)){
  for(j in c("Up", "Down")){
    locDat <- filteredGeneList[[i]][[j]]
    locColor <- ifelse(j == "Up", "#901E76", "#EDB953")
    #We will arrange these so that they go in different directions depending 
    #on if they are associated with up- or downregulated genes. 
    if(j == "Up"){
      nameOrder <- rev(locDat$Gene_names)
    } else {
      nameOrder <- locDat$Gene_names
    }
    locDat$Gene_names_factor <- factor(locDat$Gene_names, 
                                       levels = nameOrder)
    ggplot(locDat, aes(x = Gene_names_factor, y = fracAssGenes)) +
      geom_bar(stat = "identity", fill=locColor ) +
      scale_y_reverse(limits = c(maxFrac+0.02, 0), expand = c(0,0)) +
      theme_bw()
    ggsave(paste0("Results/Graphics/Stockholm/EnrichR/", i, "_", j, ".pdf"), width = 5, height = 5)
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


geneDf
´