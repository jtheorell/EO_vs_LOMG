library(ggplot2)
#Here, we are going to plot the death rate in the Oxford data against 
#the time the sample has spent in the freezer. 
deathDat <- read.csv("Data/Oxford/Gating/Live_dumpneg.csv")
metaDat <- read.csv("Data/Oxford/Documentation/Metadata_per_id.csv")

fullDeathDat <- as.data.frame(do.call("rbind", lapply(seq_along(deathDat$Label), function(x){
  locLab <- deathDat$Label[x]
  c("Label" = locLab, 
    "Live_dumpneg" = deathDat$Live_dumpneg[x], 
    "Sample_year" = metaDat$Sample_year[which(metaDat$Label == locLab)], 
    "Group" = metaDat$Group[which(metaDat$Label == locLab)],
    "Tissue" = deathDat$Tissue[x])
  
})))

#Now, the sample year is turned into a meaningful number
fullDeathDat$Time_in_freezer <- 2023-as.numeric(fullDeathDat$Sample_year)
fullDeathDat$Live_dumpneg <- as.numeric(fullDeathDat$Live_dumpneg)

#And now, we remove the thmyic samples, as they have very many more dead thymocytes
#and also are much older. 
redDeathDat <- fullDeathDat[-grep("thy|post", fullDeathDat$Tissue),]
dir.create("Results/Graphics/Oxford/Death_in_freezer")
ggplot(redDeathDat, aes(Time_in_freezer, Live_dumpneg, colour = Group)) +
  geom_point() + theme_bw() + scale_color_manual(values = c("grey", "red", "blue", "#AA00FF")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 50), expand = c(0,0))
ggsave("Results/Graphics/Oxford/Death_in_freezer/Death_rate.pdf", width = 6, height = 5)

