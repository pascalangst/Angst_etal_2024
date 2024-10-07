# compare each sample to each sample (but not to itself)

# load data
my_files <- list.files("")
my_data <- lapply(my_files, read.delim)

# results data.frame
SH_all.frame <- data.frame(NSH = vector(), Het_value = vector(), Intra_island = vector(), Intra_pond = vector(), sampleA = vector(), sampleB = vector())

for (i in 1:length(my_data)) {
  for (j in i:length(my_data)) {
    
    k <- j + 1
    
    if (k > length(my_data)) next
    
    # check if from same island or pond
    Intra_island <- gsub("-.*", "", names(my_data[i])) == gsub("-.*", "", names(my_data[k]))
    Intra_pond <- gsub("_.*", "", names(my_data[i])) == gsub("_.*", "", names(my_data[k]))
    
    data_A <- my_data[[i]]
    data_B <- my_data[[k]]
    
    # check if missing interpreted as Het (I did not see any case)
    data_A$Het_corrA <- data_A$nHet - data_A$nMiss
    data_B$Het_corrB <- data_B$nHet - data_B$nMiss
    
    Het.frame <- cbind(data_A, data_B)
    
    # get number of shared heterozygote positions
    NSH <- nrow(Het.frame[Het.frame$Het_corrA == 1 & Het.frame$Het_corrB == 1, ])
    
    # divide it by total number of heterozygote positions for each sample
    SHA <- NSH / nrow(Het.frame[Het.frame$Het_corrA == 1, ])
    SHB <- NSH / nrow(Het.frame[Het.frame$Het_corrB == 1, ])
    
    # take smaller proportion
    SH <- min(SHA, SHB)
    
    # add to results
    SH_all.frame <- rbind(SH_all.frame, cbind(NSH, SH, Intra_island, Intra_pond, names(my_data[i]), names(my_data[k])))
    
  }
}

# plot results

SH_all.frame$SH <- as.numeric(SH_all.frame$SH)
SH_all.frame$NSH <- as.numeric(SH_all.frame$NSH)
SH_all.frame$Intra_island <- as.factor(SH_all.frame$Intra_island)
SH_all.frame$Intra_pond <- as.factor(SH_all.frame$Intra_pond)

# load pre-run file for host
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
SH_all.frame_magna <- loadRData("SH_all.frame_magna_incl2023.R")

SH_all.frame_magna$SH <- as.numeric(SH_all.frame_magna$SH)
SH_all.frame_magna$NSH <- as.numeric(SH_all.frame_magna$NSH)
SH_all.frame_magna$Intra_island <- as.factor(SH_all.frame_magna$Intra_island)
SH_all.frame_magna$Intra_pond <- as.factor(SH_all.frame_magna$Intra_pond)
SH_all.frame_magna$condition <- "host"

SH_all.frame$condition <- "parasite"
SH_all.frame <- rbind(SH_all.frame,SH_all.frame_magna)

facet_names <- c(
  `host` = "D. magna",
  `parasite` = "H. tvaerminnensis"
)

ggplot(SH_all.frame, aes(x= condition, y=SH, fill = Intra_island, col = Intra_pond)) + facet_wrap(~factor(condition, levels=c("parasite","host")), scales = "free_x", labeller = as_labeller(facet_names)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        legend.position = "bottom", strip.text = element_text(face = "italic"), strip.text.x = element_text(size = 15)) +
  scale_color_manual(labels=c("Inter-Island", "Intra-Island"), values = c("#eb6223", "#522888")) +
  scale_fill_manual(labels=c("Inter-Pond", "Intra-Pond"), values = c(alpha("#eb6223", 0.5), alpha("#522888", 0.5))) +
  guides(
    col = guide_legend(order = 1, nrow=2),   # Color legend first
    fill = guide_legend(order = 2, nrow=2)   # Fill legend second
  ) +
  scale_y_continuous(name = "pairwise SH", limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1), minor_breaks = seq(0.125,0.875, 0.25)) +
  geom_hline(yintercept=0.9,linetype="dashed", color = "red4") +
  xlab(NULL)
