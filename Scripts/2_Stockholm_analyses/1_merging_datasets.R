#Here, we start by combining the datasets into one singleCellExperiment. 
#We here deviate from the previous analysis in three ways when it comes to
#the protein data: first of all, we only combine the truly identical markers, 
#i.e., we do not construct any CD14/CD33 marker and we also do not construct a
#CD45ish marker, as we are hopeful that we will be able to get properly assigned
#CD45 and CD33/CD14 expression down the line for both datasets. 
#The data has been preprocessed to this state. See 

library(SingleCellExperiment)
library(singleCellTK)
#To make it possible to store tehse files on github which as a 100 MB
#size limit, I have divided the patient file into three. 
rawPatSCE1 <- readRDS("Data/Stockholm/Raw_data/Raw_myasthenia1.rds")
rawPatSCE2 <- readRDS("Data/Stockholm/Raw_data/Raw_myasthenia2.rds")
rawPatSCE3 <- readRDS("Data/Stockholm/Raw_data/Raw_myasthenia3.rds")
rawPatSCE <- cbind(rawPatSCE1, rawPatSCE2, rawPatSCE3)
rawCtrlSCE <- readRDS("Data/Stockholm/Raw_data/Raw_Kotliarov.rds") 

#Now, we will sync the metadata and merge the datasets.
rowNamesPat <- row.names(rawPatSCE)
rowNamesCtrl <- row.names(rawCtrlSCE)

commonNamesPat <- rowNamesPat[which(rowNamesPat %in% rowNamesCtrl)]
commonNamesCtrl <- rowNamesCtrl[which(rowNamesCtrl %in% rowNamesPat)]

identical(commonNamesPat[order(commonNamesPat)], commonNamesCtrl[order(commonNamesCtrl)])
#TRUE. 

#Before reducing, we need to add the total transcriptome count sum to the control dataset
#It is already present for the other one, and is obdained in the following way: 
rawCtrlSCE$sum <- colSums(counts(rawCtrlSCE))
#So then we can reduce both datasets to this. 
rawPatSCECommon <- rawPatSCE[which(row.names(rawPatSCE) %in% commonNamesPat),]
rawCtrlSCECommon <- rawCtrlSCE[which(row.names(rawCtrlSCE) %in% commonNamesCtrl),]
rawPatSCECommonOrdered <- rawPatSCECommon[order(row.names(rawPatSCECommon)),]
rawCtrlSCECommonOrdered <- rawCtrlSCECommon[order(row.names(rawCtrlSCECommon)),]

identical(row.names(rawPatSCECommonOrdered), row.names(rawCtrlSCECommonOrdered))
#TRUE

#Now over to the much simpler task of stramlining the colData
colnames(colData(rawPatSCECommonOrdered))
geneSums<- rowSums(counts(rawPatSCECommonOrdered))
#[1] "sum"                   "detected"              "subsets_Mito_sum"     
#[4] "subsets_Mito_detected" "subsets_Mito_percent"  "total"                
#[7] "sample"                "runId"                 "batch"    

colnames(colData(rawCtrlSCECommonOrdered))
# [1] "nGene"                       "nUMI"                        "orig.ident"                 
# [4] "pctMT"                       "barcode_check"               "tenx_lane"                  
# [7] "cohort"                      "batch"                       "hash_maxID"                 
#[10] "hash_secondID"               "hash_margin"                 "hto_classification"         
#[13] "hto_classification_global"   "hash_ID"                     "adjmfc.time"                
#[16] "DEMUXLET.RD.PASS"            "DEMUXLET.N.SNP"              "DMX_GLOBAL_BEST"            
#[19] "DEMUXLET.BARCODE"            "sample"                      "sampleid"                   
#[22] "joint_classification_global" "dmx_hto_match"               "timepoint"                  
#[25] "sum"                         "detected"                    "altexps_adt_sum"            
#[28] "altexps_adt_detected"        "altexps_adt_percent"         "total"  
#So it seems possible to sync these.

colData(rawPatSCECommonOrdered) <- colData(rawPatSCECommonOrdered)[,c(1,7,9,5)]
colnames(colData(rawPatSCECommonOrdered))[4] <- "Mito_percent"
rawPatSCECommonOrdered$group <- "Pat"

altExpNames(rawPatSCECommonOrdered) <- altExpNames(rawCtrlSCECommonOrdered)
colnames(rowData(rawPatSCECommonOrdered)) <- "originalName"
#And of course there are differences in the altExps also. 

colData(rawCtrlSCECommonOrdered) <- colData(rawCtrlSCECommonOrdered)[,c(25,21,8,4)]
colnames(colData(rawCtrlSCECommonOrdered)) <- c("sum", "sample", 
                                                "batch", "Mito_percent")
rawCtrlSCECommonOrdered$group <- "Ctrl"

#Now, let us add some more colData
colMetaData <- read.csv("Data/Stockholm/Metadata/Metadata_combined_used.csv")

patMetaData <- colMetaData[-which(is.na(colMetaData$Database_ID)),]

#Now, we create a large frame that we populate with data. 
dummyPatDf <- data.frame(matrix("error", nrow = nrow(colData(rawPatSCECommonOrdered)),
                      ncol = 14))
length(which(as.vector(as.matrix(dummyPatDf)) == "error"))
#628726
colnames(dummyPatDf) <- colnames(patMetaData)[2:15]

for(i in patMetaData$Sample_ID){
    dummyPatDf[which(colData(rawPatSCECommonOrdered)$sample == i),] <- 
        patMetaData[which(patMetaData$Sample_ID == i),c(2:15)]
}
length(which(as.vector(as.matrix(dummyPatDf)) == "error"))
#0, so that worked. 
colData(rawPatSCECommonOrdered) <- cbind(colData(rawPatSCECommonOrdered), dummyPatDf)

#Now for the controls
ctrlMetaData <- colMetaData[which(is.na(colMetaData$Database_ID)),]

#Now, we create a large frame that we populate with data. 
dummyCtrlDf <- data.frame(matrix("error", nrow = nrow(colData(rawCtrlSCECommonOrdered)),
                                ncol = 14))
length(which(as.vector(as.matrix(dummyCtrlDf)) == "error"))
#769832
colnames(dummyCtrlDf) <- colnames(ctrlMetaData)[2:15]

for(i in ctrlMetaData$Sample_ID){
    dummyCtrlDf[which(colData(rawCtrlSCECommonOrdered)$sample == i),] <- 
        ctrlMetaData[which(ctrlMetaData$Sample_ID == i),c(2:15)]
}
length(which(as.vector(as.matrix(dummyCtrlDf)) == "error"))
#0, so that worked. 
colData(rawCtrlSCECommonOrdered) <- cbind(colData(rawCtrlSCECommonOrdered), dummyCtrlDf)

#Now, for the ADT data. 

#Control data cleanup
#To investigate this, the cite-seq portion of the Kotliarov data has been used in a fcs-transformed
#format, found here: https://github.com/jtheorell/Cite_seq_FCS
#The files have been investigated in FlowJo for visual simplicity.  
#There are a number of markers with insufficient or aberrant expression in the 
#control dataset. This can be either because they actually do not work, or because
#the cells have not been properly activated to induce their expression, but regardless
#we lack the transcriptomic-surface correlation needed to be able to sufficiently well
#transfer the information from one dataset to the other.

#CD14 (annoyingly too weak) 

#CD133: too widely expressed and no transcript correlation.

#CD134 (OX40) and CD137: should be expressed on activated T cells, but neither 
#HLA-DR nor CD69+ T-cells express them above background. 

#CD138 (syndecan-1): should be expressed on plasmablasts, here defined as
#CD19+CD20-CD38+(CD27+), but it is not systematically expressed there, and 
#expressed in other compartments, non-systematically.

#CD152 (CTLA-4): shoudl be expressed on T-reg (here deined as CD3+CD4+CD25+) and
#activated T cells (Here defined as CD3+CD69 or HLA-DR+), but correlates to nothing

#CD183/CXCR3, CD184/CXCR4, CD197/CCR7: signal is too weak. 

#CD206 (Mannose receptor, no expression on CD33+ cells, which it should be)

#CD223 (LAG3) no expression on T cells where it should be. 

#CD273 and CD274 (PD-L1 and PD-L2), should be expressed on myeloid APCs, but are not

#CD275 (ICOSLG) does not show any systematic expression. Possibly, this is
#because B-cells and myeloid APC need to be activated to express it, but nonetheless, 
#as it is not there, we cannot use it. 

#CD294 (CRTH2, should be expressed on ILC (a very rough definition being 
#lineage negative, CD127+), which it is not)

#CD357, TNFRSF18, should be expressed on CD3+CD4+CD25+ which it is not

#CD366, Tim-3, should be expressed on Th1 cells (here defined as CD3+CD4+CCR5+), 
#but there is no expression above background. 

#IgA and IgM are not specific enough, being expressed on large groups of
#non B-linerage cells.

#TCR-gd is sadly also expressed too lowly, and the correlation to CD3 expression
#is not at all as 100% as it should be. 

#Annexin V is negative, as all the dead cells have
#been pre-excluded, so this is removed. 

#The isotypes are also removed as they
#by definition should not show any biological signal and therefore cannot be
#plausibly inferred based on the other markers. 

#Now, we are lacking CD45 from this dataset, and after hte cleanup, CD14 will be
#lacking as well. This will however be corrected later on in the multiVI run. 

ctrlCounts <- counts(altExp(rawCtrlSCECommonOrdered))

#And then all the others are removed
ctrlCounts <- ctrlCounts[-which(row.names(ctrlCounts) %in% 
                                    c("AnnexinV_PROT", 
                                      "CD14_PROT",
                                      "CD133_PROT",
                                      "CD134_PROT",
                                      "CD137_PROT",
                                      "CD138_PROT",
                                      "CD152_PROT", 
                                      "CD184_PROT", 
                                      "CD184_PROT", 
                                      "CD197_PROT", 
                                      "CD206_PROT",
                                      "CD223_PROT",
                                      "CD273_PROT",
                                      "CD274_PROT",
                                      "CD275_PROT",
                                      "CD294_PROT",
                                      "CD357_PROT",
                                      "CD366_PROT",
                                      "IgA_PROT",
                                      "IgM_PROT",
                                      "TCRgd_PROT")),]

ctrlCounts <- ctrlCounts[-grep("sotype", row.names(ctrlCounts)),]

#Naming
rowNamesLoose <- row.names(ctrlCounts)
#We start by removing two erroneous additions of spaces in the names
rowNamesLoose <- gsub(" ", "", rowNamesLoose)

#Now, for clarity, we will change the numbering of the CD molecules, so that they
#are all three-digit. 

#Now, we divide this into specific categories
CDlowPoss <- which(grepl("CD", rowNamesLoose) &
                                         ((nchar(rowNamesLoose) == 8) |
                       (nchar(rowNamesLoose) == 9 &
                            grepl("a|b|c|d", rowNamesLoose))))

CDmidPoss <- which(grepl("CD", rowNamesLoose) &
                           ((nchar(rowNamesLoose) == 9 &!
                           grepl("a|b|c|d", rowNamesLoose)) |
                       (nchar(rowNamesLoose) == 10 &
                            grepl("a|b|c|d|L", rowNamesLoose)) |
                           (nchar(rowNamesLoose) == 11 &
                                grepl("RA|RO", rowNamesLoose))))

rowNamesLoose[CDlowPoss] <- 
    paste0("CD00", substring(rowNamesLoose[CDlowPoss], first = 3))
rowNamesLoose[CDmidPoss] <- 
    paste0("CD0", substring(rowNamesLoose[CDmidPoss], first = 3))

row.names(ctrlCounts) <- rowNamesLoose

#And now, they are ordered. 
ctrlCounts <- ctrlCounts[order(row.names(ctrlCounts)),]

#Now, we will make the patient data comparable to this. 
patCounts <- counts(altExp(rawPatSCECommonOrdered))
#We start by removing the bar codes. 
patCounts <- patCounts[1:9,]

#Now, the names are changed to reflect the naming of the Kotliarov data. 
patCountNames <- row.names(patCounts)
patCountNames <- gsub("Prot_", "", patCountNames)

patCountNames[-c(5,7,9)] <- paste0("CD0", substring(patCountNames[-c(5,7,9)], 3))
patCountNames[c(5,7,9)] <- paste0("CD00", substring(patCountNames[c(5,7,9)], 3))
patCountNames <- paste0(patCountNames, "_PROT")

row.names(patCounts) <- patCountNames
patCounts <- patCounts[order(row.names(patCounts)),]

#For some reason, the python-internal algorithms do not like the formatting of 
#these, so they do not accept their combination. For this reason, we will 
#construct the whole, integrated dataset here in four files, and then export 
#that to python. 

#We start by creating two integrated, comparable datasets. 

#Step one is to make the ADT portions fully comparable. 
nonCtrlNames <- row.names(patCounts)[-which(row.names(patCounts) %in% 
                                                row.names(ctrlCounts))]
dummyCtrlMat <- matrix(0, nrow = length(nonCtrlNames),
                                        ncol = ncol(ctrlCounts))
rownames(dummyCtrlMat) <- nonCtrlNames
ctrlCountsFull <- rbind(ctrlCounts, dummyCtrlMat)
ctrlCountsFull <- ctrlCountsFull[order(row.names(ctrlCountsFull)),]

nonPatNames <- row.names(ctrlCounts)[-which(row.names(ctrlCounts) %in% 
                                                row.names(patCounts))]
dummyPatMat <- matrix(0, nrow = length(nonPatNames),
                       ncol = ncol(patCounts))
rownames(dummyPatMat) <- nonPatNames
patCountsFull <- rbind(patCounts, dummyPatMat)
patCountsFull <- patCountsFull[order(row.names(patCountsFull)),]

identical(row.names(ctrlCountsFull), row.names(patCountsFull))
#TRUE

#So now, we are ready to integrate. 
allAdtCounts <- cbind(ctrlCountsFull, patCountsFull)

#Now, integration of the RNA
allGexCounts <- cbind(counts(rawCtrlSCECommonOrdered), counts(rawPatSCECommonOrdered))

#So, over to the col and row data. 
allCellData <- rbind(colData(rawCtrlSCECommonOrdered), colData(rawPatSCECommonOrdered))

#We will make one change here. The batches are namely confusingly called 1 and 2
#for both experiments. We will change them to something more useful. 
allCellData$batch <- paste0(allCellData$group, "_", allCellData$batch)

#To make sure that we can use the ctrl data for both, we just check that the
#row names are identical, again. 
identical(row.names(rawCtrlSCECommonOrdered), row.names(rawPatSCECommonOrdered))
#TRUE
#So on we go with one of them. 
gexRowData <- rowData(rawCtrlSCECommonOrdered)

#So now, we have two datasets that should be comparable. We will now export them
#and use them as input into python. 
dir.create("../External/Stockholm/Data/Anndata_input", recursive = TRUE)

#We need to transpose the data to get this right in the AnnData format. 
writeMM(t(allGexCounts), 
        "../External/Stockholm/Data/Anndata_input/gex_counts.mtx")
write.csv(t(as.matrix(allAdtCounts)), 
        "../External/Stockholm/Data/Anndata_input/adt_counts.csv")
write.csv(gexRowData, 
        "../External/Stockholm/Data/Anndata_input/geneData.csv")
write.csv(allCellData, 
          "../External/Stockholm/Data/Anndata_input/cellData.csv")

#Here, we add another version of the ADT counts, that only includes the
#common proteins, which will be used to generate the low-dimensional representation
#used for the logistic regression analysis, all UMAPs and the normalised version of the
#common proteins themselves. 

allAdtCounts <- read.csv("../External/Stockholm/Data/Anndata_input/adt_counts.csv", 
                         row.names = 1)
allAdtCountsCommon <- allAdtCounts[,c("CD003_PROT",
                                      "CD004_PROT",
                                      "CD008_PROT",
                                      "CD011c_PROT",
                                      "CD016_PROT",
                                      "CD019_PROT",
                                      "CD056_PROT")]

write.csv(allAdtCountsCommon, 
          "../External/Stockholm/Data/Anndata_input/adt_counts_common.csv")



