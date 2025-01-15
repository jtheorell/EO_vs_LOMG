library(SingleCellExperiment)
library(scran)
library(scater)
library(edgeR)
library(BiocParallel)

#Now, we are going to construct a few random gene sets. 
logSCE <- readRDS("../External/Stockholm/Data/SCE_versions/6_with_Ox_neigh_clust.rds")
strippedSce <- SingleCellExperiment(assays = list("counts" = counts(logSCE)),
                                    colData = colData(logSCE),
                                    rowData = rowData(logSCE))

#We then do the same pre-selections as for the edgeR analysis
proteinSce <- strippedSce[which(rowData(strippedSce)$gene_biotype == "protein_coding"),]
noMitoSce <- proteinSce[-grep("MT|GL", rowData(proteinSce)$chromosome_name),]
library(org.Hs.eg.db)
library(GO.db)
riboResults <- AnnotationDbi::select(org.Hs.eg.db, keys=c("GO:0005840"), columns = c('SYMBOL'), keytype = "GOALL")
ribo_gene_symbols <- unique(riboResults$SYMBOL)
length(ribo_gene_symbols)
#246, which seems realistic. 
noRiboSce <- noMitoSce[-which(rowData(noMitoSce)$hgnc_symbol %in% ribo_gene_symbols),]
#We will also make one further pre-exclusion: only if a transcript above 10, 
#i.e. that it either is expressed quite well by a few cells or that 10 cells 
#express it lowly, will it be accepted. 
sceSums <- rowSums(counts(noRiboSce))
expressedSce <- noRiboSce[which(sceSums > 10),]

nrow(strippedSce)-nrow(expressedSce)
#So a total number of 6677 genes are removed in all these steps. 
#Now, we remove row.names that are non-unique. About 1000 transcripts will be gone here. 
allSymbols <- rowData(expressedSce)$hgnc_symbol
duplicatedNames <- unique(allSymbols[which(duplicated(allSymbols))])
duplicatedNames #""

namedSce <- expressedSce[-which(allSymbols %in% duplicatedNames),]

rownames(namedSce) <- rowData(namedSce)$hgnc_symbol
#Now, we are down to 13502 genes. 

#Here, we are going to generate populations that are the same size as our hit populations. 
signPops <- row.names(read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1))

signPopSizes <- lapply(signPops, function(x){
  length(which(logSCE$oxNeighClust == x))
})
names(signPopSizes) <- signPops

#Now, we generate gene lists based on this. 
patGroups <- c("EOMG", "LOMG")
randomEdgeROutcomes <- list()
dir.create("../External/Stockholm/Data/RandomPopEdgeR")
for(i in 1:100){
  print(i)
  set.seed(i)
  randomEdgeROutcomes[[i]] <- bplapply(signPops, function(x){
    print(x)
    
    locMother <- gsub("|_.OMGlow_.", "\\1", x)
    #Now, we construct a random cluster which is the same size as the significant cluster
    
    locDaughterCols <- sample(which(logSCE$cellType == locMother), signPopSizes[[x]])
    print(length(locDaughterCols))
    
    locMotherCols <- which(logSCE$cellType == locMother)
    locOtherCols <- locMotherCols[-which(locMotherCols %in% locDaughterCols)]
    #Here, we reduce the data if it is unnecessarily extensive
    if(length(locOtherCols) > 5000){
      locOtherCols <- sample(locOtherCols, 5000)
    }
    
    locClusterCols <- c(locDaughterCols, locOtherCols)
    
    hitClustFull <- rep(FALSE, ncol(logSCE))
    
    hitClustFull[locDaughterCols] <- TRUE
    hitClust <- hitClustFull[locClusterCols]
    
    locSCE <- namedSce[,locClusterCols]
    locSCE$hitClust <- hitClust
    
    current <- aggregateAcrossCells(locSCE, 
                                    ids=colData(locSCE)[,c("sample", "hitClust")])
    
    # Creating up a DGEList object for use in edgeR:
    y <- DGEList(counts(current), samples=colData(current))
    
    #Here, we discard the small samples
    discarded <- current$ncells < 10
    y <- y[,!discarded]
    
    #Here, we need to be clever. In normal instances, there is often a reasonable
    #balance in cell number between the compared groups. In this case, however, 
    #the cells we are interested in are very, very much fewer. Therefore, we
    #cannot use the same criteria to define which genes that are of interest for
    #both groups. Instead, we will select genes that are either massively expressed
    #by the large cluster, or reasonably expressed by ours. 
    yFoc <- y[,y$samples$hitClust]
    
    keep <- suppressWarnings(filterByExpr(yFoc,min.count = 10, 
                                          min.total.count = 15, large.n = 20))
    
    yOther <- y[,which(y$samples$hitClust==FALSE)]
    
    keep2 <- suppressWarnings(filterByExpr(yOther,min.count = 100, 
                                           min.total.count = 150, large.n = 20))
    keep[keep2] <- TRUE
    
    print(length(which(keep)))
    y <- y[keep,]
    
    y <- calcNormFactors(y, method = "upperquartile")
    
    design <- model.matrix(~factor(hitClust), y$samples)
    
    
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit, coef=ncol(design))
    
    #Now, we pick out all the ones that are significantly upregulated. If there
    #are none, which is the case for the CD4T, that have too few cells, we will
    #take the top 10
    topRes <- topTags(res)
    resSigEdgeR <- topRes@.Data[[1]]
    
    #And these are ordered according to their fold change
    resSigEdgeROrd <- resSigEdgeR[order(abs(resSigEdgeR$logFC), decreasing = TRUE),]
    #Now, the results are sent back
    resSigEdgeROrd
    
  })
  saveRDS(randomEdgeROutcomes[[i]], paste0("../External/Stockholm/Data/RandomPopEdgeR/Run_number_", i, ".rds"))
}

#And now, the whole thing is saved
saveRDS(randomEdgeROutcomes, "../External/Stockholm/Data/RandomPopEdgeR/All_together.rds")
