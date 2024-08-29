#We start by importing the files, that after export have gone through a round 
#of NameChanger name changes, to make sure that the timing (if PBMC pre if before thymectomy
#and post if after) or cell type is clear. Importantly, for the JR numbered 
#controls, there are two
#formats in the original names: JRXXXX_XXXX, and JRXXXX. The former is, very 
#counterintuitively, showing the number of the individual sample first and then
#the number of the patient. For this reason, I have removed the first number
#in nameChanger, before going into this analysis. 
#The cells are importantly exported as channel values, as this corresponds to the
#transformed values. 
library(FNN)
filenames <- list.files(path = "../External/Oxford/Exported_csv/BT_all", 
                        full.names = TRUE)

#Here, we import and subsample at the same time.
set.seed(1100)
big_all_file <- do.call("rbind", lapply(filenames, function(x){
  ret <- read.csv(x)
  if(nrow(ret) > 5000){
    ret <- ret[sample(1:nrow(ret), 5000),]
  }
  id <- gsub(pattern="../External/Oxford/Exported_csv/BT_all/BT_lymphocytes_|_.+\\.csv", "\\1", x)
  ret$id <- id
  tissue <- gsub(pattern="../External/Oxford/Exported_csv/BT_all/BT_lymphocytes_.+_|\\.csv", "\\1", x)
  ret$tissue <- tissue 
  ret 
}))

#Now, we start adding in metadata. 
#There are only three samples that are post-thymectomy, but we know how long 
#post thymectomy they are, so we are adding that data here:
big_all_file$postTxDays <- NA
big_all_file$postTxDays[which(big_all_file$id == "EJon" &
                                big_all_file$tissue == "post")] <- 12
big_all_file$postTxDays[which(big_all_file$id == "EHal" &
                                big_all_file$tissue == "post")] <- 51
big_all_file$postTxDays[which(big_all_file$id == "PJor" &
                                big_all_file$tissue == "post")] <- 32
table(big_all_file$postTxDays, big_all_file$tissue, useNA = "ifany")
#          post     pre     thy
#  12     25000       0       0
#  32     17719       0       0
#  51     25000       0       0
#  <NA>       0 1290894  250000

#All other metadata is not dependent on the tissue, and can thus be imported on a
#per-id basis. 
metadata <- read.csv("Data/Oxford/Documentation/Metadata_per_id.csv")

big_all_file$group <- big_all_file$age <- 
  big_all_file$sex <- big_all_file$hyperplasia <- big_all_file$sampleYear <- NA

for(i in metadata$Label){
  print(i)
  locMetRow <- which(metadata$Label == i)
  locDatRows <- which(big_all_file$id == i)
  big_all_file$group[locDatRows] <- metadata$Group[locMetRow]
  big_all_file$age[locDatRows] <- metadata$Age[locMetRow]
  big_all_file$sex[locDatRows] <- metadata$Sex[locMetRow]
  big_all_file$hyperplasia[locDatRows] <- metadata$Thymic_hyperplasia[locMetRow]
  big_all_file$sampleYear[locDatRows] <- metadata$Sample_year[locMetRow]
}

#Now, checking that it worked: 
table(big_all_file$group, useNA = "always")
#   Ctrl    EOMG    LOMG preLOMG    <NA> 
# 100000  115000   90000   20000       0
table(big_all_file$sex, useNA = "always")
#      F       M    <NA> 
#170000 155000      0

#We also import the lymphocyte counts that will be used downstream for
#normalisation purposes
big_all_file$idTime <- paste0(big_all_file$id, "_", big_all_file$tissue)

#For simplicity's sake, we stick with the calculations based on the lymphocyte
#counts, eventhough of course in this case, it is the same number as the number
#of cells from the donor. In this case, it is also always 50000 cells. 
lymphCounts <- as.data.frame(as.matrix(table(big_all_file$idTime)))
colnames(lymphCounts) <- "Lymphocyte_count"
#TRUE
big_all_file$lymphCount <- NA
for(i in unique(big_all_file$idTime)){
  big_all_file$lymphCount[which(big_all_file$idTime == i)] <- 
    lymphCounts$Lymphocyte_count[which(row.names(lymphCounts) == i)]
}
table(big_all_file$lymphCount, useNA = "always")
#No NA left. 

#And with that, we are ready to save the file. 
dir.create("../External/Oxford/Resulting_data/BT_all")
saveRDS(big_all_file, "../External/Oxford/Resulting_data/BT_all/BT_all_full_file.rds")
