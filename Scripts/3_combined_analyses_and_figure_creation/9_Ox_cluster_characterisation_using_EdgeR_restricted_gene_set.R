library(SingleCellExperiment)
library(scran)
library(scater)
library(edgeR)
library(BiocParallel)

#Now, we are going to analyse the interesting clusters here they are: 
significantClusters <- read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1)
logSCE <- readRDS("../External/Stockholm/Data/SCE_versions/6_with_Ox_neigh_clust.rds")

runVec <- row.names(significantClusters)

strippedSce <- SingleCellExperiment(assays = list("counts" = counts(logSCE)),
                                    colData = colData(logSCE),
                                    rowData = rowData(logSCE))

#Here, we make a few restrictions on the transcript side. First, we hone in 
#on protein coding genomic transcripts. 
proteinSce <- strippedSce[which(rowData(strippedSce)$gene_biotype == "protein_coding"),]
#Here, 3000 transcripts are lost. 
#Now, we remove the mitochondrial as well as the experimental transcripts, which 
#are less than 20 in total, but that still make up a cnsiderable portion of the 
#transcriptome. 

noMitoSce <- proteinSce[-grep("MT|GL", rowData(proteinSce)$chromosome_name),]

#Now, we go on with the hardest one: we will remove the ribosomal transcripts
#as these make up a very large portion of highly expressed transcripts, 
#and tend to distort the whole analysis, not the least by giving rise to 
#gene set enrichments as informative as "large ribosome", etc.
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

patGroups <- c("EOMG", "LOMG")
set.seed(110)
edgerOutcomes <- bplapply(runVec, function(x){
    print(x)
    if(grepl("TCRgd", x)){
        locMother == "gdT"
    } else {
        locMother <- gsub("|_.OMG.+_.", "\\1", x)
    }
    
    signGroup <- gsub("|_.$", "\\1", x)
    #Now, we compare the specific cluster to all the other cells within the 
    #mother cluster, apart from the clusters significant from the same patient
    #group and regulated in the same direction. 
    locDaughterCols <- which(logSCE$oxNeighClust == x)
    
    locOtherCols <- which(logSCE$cellType == locMother &
                              (grepl(signGroup, logSCE$oxNeighClust) == FALSE))
    #Here, we reduce the data if it is unnecessarily extensive
    if(length(locOtherCols) > 5000){
        locOtherCols <- sample(locOtherCols, 5000)
    }
    
    locClusterCols <- c(locDaughterCols, locOtherCols)
    
    discriminantFull <- rep(paste0(1, "_other"), ncol(logSCE))
    
    discriminantFull[locDaughterCols] <- paste0(2, "_", x)
    discriminant <- discriminantFull[locClusterCols]
    
    locSCE <- namedSce[,locClusterCols]
    locSCE$discriminant <- discriminant
    #Here, we are identrifying the normalised counts as counts, to satisfy the software. 
    #but we will not ruin things by doubly transforming the data below. 
    current <- aggregateAcrossCells(locSCE, 
                                    ids=colData(locSCE)[,c("sample", "discriminant")])
    
    # Creating up a DGEList object for use in edgeR:
    y <- DGEList(counts(current), samples=colData(current))
    
    #Here, we discard the small samples
    discarded <- current$ncells < 10
    y <- y[,!discarded]
    
    y <- calcNormFactors(y)
    design <- model.matrix(~factor(discriminant), y$samples)
    
    keep <- filterByExpr(y, design = design)
    y <- y[keep,]
    
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit, coef=ncol(design))$table
    
    #Now, we pick out all the ones that are significantly upregulated. 
    resSigEdgeR <- res[which(res$PValue < 0.05),]
    
    #And these are ordered according to their fold change
    resSigEdgeROrd <- resSigEdgeR[order(abs(resSigEdgeR$logFC), decreasing = TRUE),]
    
})
names(edgerOutcomes) <- runVec
dir.create("Results/Data/Stockholm")
saveRDS(edgerOutcomes, "Results/Data/Stockholm/EdgeR_significant_transcriptomes.rds")
#Now, we will create a considerably reduced set for plotting, where we show
#only the ones with a logFC above 1. 

range(unlist(lapply(edgerOutcomes, function(x){
    length(which(abs(x$logFC) > 1))
})))
#64 802

#Now, we create a big one based on this, only collecting the logFC data for all above logFC 1.

allUsedGeneNames <- unique(unlist(lapply(edgerOutcomes, function(x){
    row.names(x)[which(abs(x$logFC) > 1)]
} )))

length(allUsedGeneNames) #2475

allUsedGeneNamesOrdered <- allUsedGeneNames[order(allUsedGeneNames)]

allMat <- do.call("cbind", lapply(edgerOutcomes, function(x){
    locRes <- unlist(lapply(allUsedGeneNamesOrdered, function(y){
        if(y %in% row.names(x)){
            x$logFC[which(row.names(x) == y)]
        } else {
            NA
        }
    }))
    names(locRes) <- allUsedGeneNamesOrdered
    locRes
}))

write.csv(allMat, "Results/Data/Stockholm/Top_transcripts_EdgeR.csv")

