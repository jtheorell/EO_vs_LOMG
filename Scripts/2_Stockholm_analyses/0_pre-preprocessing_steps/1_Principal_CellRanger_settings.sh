#These settings and versions were used to generate the count matrix files. 
#For Gene expression data: 
module load cellranger/7.0.1
cellranger count --id [AN ID] \
--transcriptome=$CELLRANGER_DATA/refdata-gex-GRCh38-2020-A \
--fastqs [folder where the resulting files were kept]

#For VDJ and cell surface data, cellranger version 6.1.2 were used. Principal
#library data for the cell surface information can be found in 
#Data/Stockholm/Input_files_for_Bianca
