#We will, as always, follow the OSCA work very closely. As it is constantly
#expanding, it is hard to provided good, long-lived links, but here is an attempt. 
#OSCA advanced: http://bioconductor.org/books/3.14/OSCA.advanced/
#OSCA basic: http://bioconductor.org/books/3.14/OSCA.basic/

#To reduce the work load, we will start by running the six runs separately, 
#despite run 1-3 (19-21) being from the same 8 individuals and the same being
#true for the last three. This script will look identical for all six runs. 
#We start by importing the gene expression ../External/Stockholm/Data, and sorting out the 
#empty droplets. 
#There was an error in the naming, so all the files have lost the two 
#first letters and have "an" there instead. This means that anw means raw. I
#Have renamed them all for clarity.
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dropletUtils)
library(scuttle)
#The difference between the runs is only that they are called 019, 020, etc, 
#in all places below. Here the 019 run is shown. 
commonName <- "GE_019"

#This is taken from 
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices#r-load-mat
#In this case we have pre-selected the filtered results from cellRanger, as that
#produces very similar results to doing the droplet exclusion myself. 
matrix_dir = "../External/Stockholm/Data/MG_10X_GE_cellRanger7/GE_019/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
raw_mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(raw_mat) = barcode.names$V1
rownames(raw_mat) = feature.names$V1

#So now step one is to kill off the empty barcodes, before we start loading anything
#else onto here. 
sce.pbmc <- SingleCellExperiment(assays = list(counts =raw_mat))

#Here, we add meta../External/Stockholm/Data for the genes.
sce.pbmc   <- getBMFeatureAnnos(sce.pbmc, filters = "ensembl_gene_id", 
                             attributes = c("ensembl_gene_id", "hgnc_symbol",
                                            "chromosome_name", "gene_biotype", 
                                            "start_position", "end_position"), 
                             ../External/Stockholm/Dataset = "hsapiens_gene_ensembl")


row../External/Stockholm/Data(sce.pbmc)$Symbol <- feature.names$V2
is.mito <- grep("MT", row../External/Stockholm/Data(sce.pbmc)$chromosome_name)
sce.pbmc <- addPerCellQCMetrics(sce.pbmc, subsets=list(Mito=is.mito))
colnames(col../External/Stockholm/Data(sce.pbmc))

reasons <- perCellQCFilters(sce.pbmc, 
                            sub.fields=c("subsets_Mito_percent"), nmads = 5)
colSums(as.matrix(reasons))
#low_lib_size            low_n_features high_subsets_Mito_percent 
#27                        16                       295 
#discard 
#303

#Makes sense to keep the thresholds as high as this, looking at the plots below.
dir.create(paste0("Diagnostics/Stockholm/", commonName))
plot(log(sce.pbmc$sum), sce.pbmc$subsets_Mito_percent,
     xlab="Total count", ylab='Mitochondrial %')
abline(h=attr(reasons$high_subsets_Mito_percent, "thresholds")["higher"], col="red")
abline(v=log(attr(reasons$low_lib_size, "thresholds"))["lower"], col="blue")
dev.copy(pdf, paste0("Diagnostics/Stockholm/", commonName, "/Mito_vs_lib_size_w_thresholds.pdf"))
dev.off()

hist(log(sce.pbmc$detected), breaks = 200)
abline(v=log(attr(reasons$low_n_features, "thresholds"))["lower"], col="blue")
dev.copy(pdf, paste0("Diagnostics/Stockholm/", commonName, "/N_features_w_thresholds.pdf"))
dev.off()

sce.pbmc <- sce.pbmc[,-which(reasons$discard)]


#And with that, we are done with filtering out cells. Now it is time to bring
#in the surface ../External/Stockholm/Data
matrix_dir = "../External/Stockholm/Data/Surface/VDJ_21_019_surface/"
barcode.path <- paste0(matrix_dir, "raw_feature_bc_matrix_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "raw_feature_bc_matrix_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "raw_feature_bc_matrix_matrix.mtx.gz")
raw_mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(raw_mat) = barcode.names$V1
rownames(raw_mat) = feature.names$V2

#Now we exclude all that did not have a transcriptome.
reduced_mat <- raw_mat[,which(colnames(raw_mat) %in% colnames(sce.pbmc))]

#Checking that the cells have the same order
identical(colnames(reduced_mat), colnames(sce.pbmc))
#TRUE

surface../External/Stockholm/Data <- SingleCellExperiment(assays = list(counts =reduced_mat))

altExp(sce.pbmc) <- surface../External/Stockholm/Data
altExp(sce.pbmc) <- addPerCellQCMetrics(altExp(sce.pbmc))

#Now, we should not have to run a second purge of empty droplets, as we already
#selected the cells with large and diverse transcriptomes, but here comes a check
plot(sce.pbmc$total, altExp(sce.pbmc)$total,log="xy",
     xlab="Total gene count", ylab="Total HTO count")
dev.copy(pdf, paste0("Diagnostics/Stockholm/", commonName, "/Gene_count_vs_surface_count.pdf"))
dev.off()
#
#This showed not surprisingly that the vast majority of cells have considerable 
#surface counts, with only a very minor number of exceptions, which anyway
#will be excluded in the next step, as these cells will not be associated 
#to any donor.

hto.mat <- counts(altExp(sce.pbmc))[grep("B2M_", row.names(altExp(sce.pbmc))),]
hash.stats <- hashedDrops(hto.mat)

#This below looks better than in the setup ../External/Stockholm/Data
hist(hash.stats$LogFC, xlab="Log fold-change from best to second HTO", main="", breaks = 100)
table(hash.stats$Best[hash.stats$Confident])
#1   2   3   4   5   6   7   8 
#638 865 954 668 770 509 869 683 

#At this stage, we will just believe that the current strategy is good enough: 
#we will not attempt to improve it by our own analysis strategies. 
sce.pbmc$sample <- hash.stats$Best
sce.pbmc <- sce.pbmc[,hash.stats$Confident]

dir.create("../External/Stockholm/Data/SCE")
saveRDS(sce.pbmc, "../External/Stockholm/Data/SCE/Run_019.rds")
