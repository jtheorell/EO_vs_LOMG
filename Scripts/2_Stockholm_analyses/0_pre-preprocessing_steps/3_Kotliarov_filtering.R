library(SingleCellExperiment)
library(biomaRt)
library(scRNAseq)
library(scater)
library(scran)
library(scuttle)

commonName <- "Kotliarov"
#Here, we let the Kotliarov cells go throught he same filtering as our cells. 
rawCtrlSCE <- KotliarovPBMCData(mode = c("rna", "adt"), ensembl = TRUE, location = TRUE)

rawCtrlSCE   <- getBMFeatureAnnos(rawCtrlSCE, filters = "ensembl_gene_id", 
                                attributes = c("ensembl_gene_id", "hgnc_symbol",
                                               "chromosome_name", "gene_biotype", 
                                               "start_position", "end_position"), 
                                dataset = "hsapiens_gene_ensembl")



rawCtrlSCE <- addPerCellQCMetrics(rawCtrlSCE)
colnames(colData(rawCtrlSCE))

#Now, we divide the data into the two batches and run the QC separately, as we 
#have done that (and even been further granular, actually, runing each cell separately)
#for the other experiment
batch1 <- rawCtrlSCE[,which(rawCtrlSCE$batch == 1)]

reasons1 <- perCellQCFilters(batch1, 
                            sub.fields=c("pctMT"), nmads = 5)
colSums(as.matrix(reasons1))
#low_lib_size low_n_features     high_pctMT        discard 
#445            790           1395           1550 

dir.create(paste0("Diagnostics/Stockholm/", commonName))
plot(log(batch1$sum), batch1$pctMT,
     xlab="Log total count", ylab='Mitochondrial %')
abline(h=attr(reasons1$high_pctMT, "thresholds")["higher"], col="red")
abline(v=log(attr(reasons1$low_lib_size, "thresholds"))["lower"], col="blue")
dev.copy(pdf, paste0("Diagnostics/Stockholm/", commonName, "/Mito_vs_lib_size_w_thresholds_batch_1.pdf"))
dev.off()

hist(log(batch1$detected), breaks = 200)
abline(v=log(attr(reasons1$low_n_features, "thresholds"))["lower"], col="blue")
dev.copy(pdf, paste0("Diagnostics/Stockholm/", commonName, "/N_features_w_threshold.pdf"))
dev.off()

batch2 <- rawCtrlSCE[,which(rawCtrlSCE$batch == 2)]

reasons2 <- perCellQCFilters(batch2, 
                            sub.fields=c("pctMT"), nmads = 5)
colSums(as.matrix(reasons2))
#low_lib_size low_n_features     high_pctMT        discard 
#534           1008           1882           2116 

dir.create(paste0("Diagnostics/Stockholm/", commonName))
plot(log(batch2$sum), batch2$pctMT,
     xlab="Log total count", ylab='Mitochondrial %')
abline(h=attr(reasons2$high_pctMT, "thresholds")["higher"], col="red")
abline(v=log(attr(reason2$low_lib_size, "thresholds"))["lower"], col="blue")
dev.copy(pdf, paste0("Diagnostics/Stockholm/", commonName, "/Mito_vs_lib_size_w_thresholds_batch_2.pdf"))
dev.off()

hist(log(batch2$detected), breaks = 200)
abline(v=log(attr(reasons2$low_n_features, "thresholds"))["lower"], col="blue")
dev.copy(pdf, paste0("Diagnostics/Stockholm/", commonName, "/N_features_w_threshold_batch_2.pdf"))
dev.off()

batch1 <- batch1[,-which(reasons1$discard)]
batch2 <- batch2[,-which(reasons2$discard)]

filtCtrlSCE <- cbind(batch1, batch2)

dir.create("../External/Stockholm/Data/SCE")
saveRDS(filtCtrlSCE, "../External/Stockholm/Data/SCE/Raw_Kotliarov.rds")
