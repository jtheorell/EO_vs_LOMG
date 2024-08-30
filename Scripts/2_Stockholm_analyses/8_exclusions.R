#Now, we are going to exclude cells, both if they display an incompatible pattern
#of AIRR expression and phenotype, if they are defined as doublets or if they simply
#display surface phenotype features of doublets (mainly co-expression of CD3 and CD19
#in this dataset)
library(SingleCellExperiment)
fullSce <- readRDS("../External/Stockholm/Data/SCE_versions/4_including_cell_type.rds")

dim(fullSce)
#20182 99897

#First, we exclude all cells defined as doublets
realSce <- fullSce[,-which(fullSce$cellType == "doublet")]
dim(realSce) #20182 99704, so here we lose a mere 193 cells. 

#Now, we remove TCR-expressing non-T-cells
realSce <- realSce[,-which(realSce$TCR &! grepl("T", realSce$cellType))]
dim(realSce) #20182 99451, so this means the loss of another 253 cells. 

#Here, we instead remove all BCR-expressing non-B cells. 
realSce <- realSce[,-which(realSce$BCR &! realSce$cellType == "B")]
dim(realSce) #20182 99434, losing 17 cells here. 

#Now, we exclude cells with more than one BCR heavy chain. 
realSce <- realSce[,-which(realSce$BCR_duplicated)]
dim(realSce) #20182 99403 so another 31 cells taken away. 

#With these techniques, a total of 494 cells, or 0.6% of the cells were lost. 

#At this stage, we exclude the myeloid cells, as they are highly problematic, 
#given that the EOMG/LOMG cohorts differ somewhat in the freezing protocol, 
#and that, for this reason, it is known that much fewer myeloid cells have
#survived in the EOMG cohort. 
nonMyeloidSCE <- realSce[,-which(realSce$cellType == "Myeloid")]

dim(nonMyeloidSCE) #20182 81695, this on the otherside,is a devastating, but 
#unavoidable loss of 20% of the data. 

#We here also make another large exclusion: we do not need the controls for any
#purposes further downstream, and they are therefore removed here. 
patSCE <- nonMyeloidSCE[,-which(nonMyeloidSCE$group == "Ctrl")]
dim(patSCE)
#20182 35734
#So this is it. 

#Now, we are adding one more piece of information here, and that is that two
#of the individuals have an uncertain MG subtype. These are defined here: 
patSCE$MG_type[which(patSCE$MG_type == "LOMG" & 
                         patSCE$Age == 41)] <- "uncertain" 
patSCE$MG_type[which(patSCE$Age == 50)] <- "uncertain" 

saveRDS(patSCE, "../External/Stockholm/Data/SCE_versions/5_post_exclusions.rds")

#We also export a protein and coldata object here
sthlmProtDat <- altExp(patSCE)
colData(sthlmProtDat) <- colData(patSCE)
saveRDS(sthlmProtDat, "../External/Stockholm/Data/SCE_versions/Stockholm_protein_data_including_colData.rds")
