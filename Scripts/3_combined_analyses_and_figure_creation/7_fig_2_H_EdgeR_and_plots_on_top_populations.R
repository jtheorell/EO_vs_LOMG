#Here, we will identify which markers that define the hit populations compared
#to the rest of the non-hit cells from the same cells. We will use EdgeR for this purpose, as has been done
#before in the flow world, see e.g. 
#https://www.researchgate.net/post/Flow_Cytometry_Data_Anlaysis_tools
library(edgeR)
library(scuttle)
library(SingleCellExperiment)
library(ggplot2)
library(tidysdm)
library(reshape2)

#The 
pValFoc <- 
  read.csv("Results/Data/Oxford_and_Stockholm/Harmonisation/Sign_in_Stock_and_Ox.csv", row.names = 1)

popVec <- c("B", "CD4T", "CD8T", "NK")

markerList <- list("B" = c("CD38", "CCR7", "CD45", "CCR4", "IgD.TCRgd", 
                           "CD11b", "CD27", "CD20", "CD45RA", "CD31", 
                           "CD24", "IgM"),
                   "CD4T" = c("CD38", "KLRG1", "CCR7", "CD25", "CD4", "CCR4", "CD71", 
                              "CD11b", "CD27", "CD64", "CD20", "CD45RA", "CD2", "CD31", "CD161", "CD7",
                              "CD24", "CD152"),
                   "CD8T" = c("CD38", "KLRG1", "CCR7", "CD25", "CCR4", "CD71", "CD8",
                              "CD11b", "CD27", "CD64", "CD20", "CD45RA", "CD2", "CD31", "CD161", "CD7",
                              "CD24", "CD152"),
                   "NK" = c("NKp30","HLADR","CCR7", 
                            "CD183","CD161","CD16",
                            "CD34","CD56", "CD8",
                            "CD127", "CD27", "CD57",
                            "NKG2A", "CD2", "NKG2C",
                            "CRTH2", "CD7","CD117","CD218a"))

dir.create("Results/Graphics/Oxford_and_Stockholm/Violins_on_significant")
locEdgeResultAll <- list()
for(i in popVec){
  print(i)
  datDir <- paste0("../External/Oxford/Resulting_data/", i, "/")
  locFile <- readRDS(paste0(datDir, i, "_full_file_post_Euclid_including_subclusts.rds"))
  locHitPops <- gsub("[BCD4T8NKTCRgd]+_|", "", row.names(pValFoc)[grep(i, row.names(pValFoc))])
  popMarks <- markerList[[i]]
  locEdgeResult <- 
    lapply(locHitPops, function(x){
      focPop <- locFile[which(locFile$smoothGroupDeep == x),c(popMarks, "id")]
      otherPop <- locFile[-which(locFile$smoothGroupDeep == x),c(popMarks, "id")]
      #Here, we downsample the data, to speed things up. 
      set.seed(111)
      downSampPops <- lapply(list(focPop, otherPop), function(y){
        if(nrow(y) > 100000){
          y[sample(1:nrow(y), 100000),]
        } else {y}
      })
      reshuffleDat <- do.call("rbind", downSampPops)
      discriminant <- c(rep("Focus", nrow(downSampPops[[1]])),
                        rep("Referece", nrow(downSampPops[[2]])))
      sampleVec <- reshuffleDat$id
      locDatProt <- reshuffleDat[,-which(colnames(reshuffleDat) == "id")]
      #Here, we are identrifying the normalised counts as counts, to satisfy the software. 
      #but we will not ruin things by doubly transforming the data below. 
      current <- aggregateAcrossCells(SingleCellExperiment(assays = list("counts" = t(locDatProt))), 
                                      ids=DataFrame(sampleVec, discriminant))
      
      # Creating up a DGEList object for use in edgeR:
      y <- DGEList(counts(current), samples=colData(current))
      
      #Here, we discard the small samples
      discarded <- current$ncells < 10
      table(discarded)
      y <- y[,!discarded]
      
      design <- model.matrix(~factor(discriminant), y$samples)
      y <- estimateDisp(y, design)
      fit <- glmQLFit(y, design, robust=TRUE)
      res <- glmQLFTest(fit, coef=ncol(design))$table
      #Now, we add one more column, namely adjusted p-values
      res$p.adjust <- p.adjust(res$PValue, method = "fdr")
      #And now, we combine the new protein file
      locDatProt$Group <- discriminant
      list("edgerRes" = res, "protDat" = locDatProt)
    })
  names(locEdgeResult) <- locHitPops
  locEdgeResultAll[[i]] <- locEdgeResult
  #Now, we zoom in on the markers relevant to any population under the cell type
  focMarks <- unique(unlist(lapply(lapply(locEdgeResult, "[[", 1), function(x){
    row.names(x)[which(x$p.adjust < 0.05 & abs(x$logFC) > 0.2)]
  })))
  nameOrd <- focMarks[order(focMarks)]
  CDMols <- grep("CD", nameOrd)
  nameOrd[CDMols] <- 
    nameOrd[CDMols][order(as.numeric(gsub("CD|RA|a|b", "", nameOrd[CDMols])))]
  #And now for the plotting
  lapply(names(locEdgeResult), function(x){
    locSubset <- locEdgeResult[[x]]
    locDat <- locSubset$protDat[,c(nameOrd, "Group")]
    locDatLong <- melt(locDat, id.vars = "Group")
    locDatLong$variable <- factor(locDatLong$variable, 
                                     levels = nameOrd)
    if(grepl("EOMG", x)){
      colorForPlot <- "red4"
    } else {
      colorForPlot <- "blue4"
    }
    p <- ggplot(locDatLong, aes(
      x = variable,
      y = value,
      fill = Group
    )) +
      geom_split_violin(scale = "width") +
      theme_bw() + scale_fill_manual(values = c(colorForPlot, "grey")) +
      scale_y_continuous(limits = c(0, 1050), position = "right")  
    ggsave(paste0("Results/Graphics/Violins_on_significant/", 
                  i, "_", x, ".pdf"), plot = p)
    p <- p+theme(legend.position = "None", 
                 axis.title = element_blank())
    ggsave(paste0("Results/Graphics/Violins_on_significant/", 
                  i, "_", x, "_stripped.pdf"), plot = p, 
           width = 7, height = 2)
  })
  
}