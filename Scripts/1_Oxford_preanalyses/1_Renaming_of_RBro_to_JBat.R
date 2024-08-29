#Here, we are correcting the internal naming of the RBro files, to JBat, as
#this sample was wrongly taken from the freezer in the first place. I have
#already renamed the files, but the internal naming still seems to be 
#RBro, as FlowJo keeps mislabeling them. 

library(flowCore)

jBatFiles <- list.files("../External/Oxford/FCS/JBat_files/Pre_FIL_change", full.names = TRUE)

for(i in jBatFiles){
  print(i)
  locFF <- read.FCS(i, transformation = FALSE, truncate_max_range = FALSE)
  keyword(locFF)$TUBENAME <- "JBat"
  keyword(locFF)$`$FIL` <- gsub("RBro", "JBat", keyword(locFF)$`$FIL`)
  write.FCS(locFF, gsub("Pre", "Post", i))
}

#After this procedure had commenced, the two unmixed JBat files were moved back
#in with the other samples, to avoid future confusion or missing out of this sample. 



