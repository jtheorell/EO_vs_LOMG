metadata <- read.csv("Data/Oxford/Documentation/Metadata_per_id.csv")

metadata$ageChunks <- "<30"
metadata$ageChunks[which(metadata$Age >=30 & metadata$Age < 40)] <- "30-39"
metadata$ageChunks[which(metadata$Age >=40 & metadata$Age < 50)] <- "40-49"
metadata$ageChunks[which(metadata$Age >=50 & metadata$Age < 60)] <- "50-59"
metadata$ageChunks[which(metadata$Age >=60 & metadata$Age < 70)] <- "60-69"
metadata$ageChunks[which(metadata$Age >=70 )] <- ">70"

table(metadata$ageChunks, metadata$Group, metadata$Sex)

#F
#        Ctrl EOMG LOMG preLOMG
#  <30      2    1    0       0
#  30-39    4    4    0       0
#  40-49    2    2    0       0
#  50-59    1    0    2       0
#  60-69    2    0    1       0
#  >70      2    0    4       0

#M
#        Ctrl EOMG LOMG preLOMG
#  <30      0    2    0       0
#  30-39    1    4    0       2
#  40-49    1    0    0       0
#  50-59    1    0    1       0
#  60-69    3    0    5       0
#  >70      1    0    4       0

#This shows that we have two young male EOMG donors that sort of derail the balance,
#as we lack controls for them. We will however keep them in anyway. It also shows that we can
#call all controls below 50 young and all controls above 50 old. 
