## Script to calculate pairwise FST between samples and plot it

library(poolfstat)  # FST calculation of pools
library(ggplot2)    # plotting
library(reshape2)   # plotting
library(geodist)    # measure geographic distance
library(ggpubr)     # plotting and stats
library(adespatial) # dbmem by RDA
library(vegan)      # dbmem by RDA
library(data.table) # transpose data
library(tibble) 	  # plotting
library(tidyr)      # transform data
library(ggside)     # for plotting
library(cowplot)    # for plotting
library(dplyr)      # for stats

# get samplenames
samplenames <- system("bcftools query -l VCF.vcf", intern = TRUE)

# import data
pooldata<-vcf2pooldata("VCF.vcf", poolsizes=rep(100,length(samplenames)), poolnames=samplenames)

# calculate pairwiseFST between all samples
PairwiseFST<-compute.pairwiseFST(pooldata)

# set negative values to 0
PairwiseFST@PairwiseFSTmatrix[PairwiseFST@PairwiseFSTmatrix < 0 & is.finite(PairwiseFST@PairwiseFSTmatrix)] <- 0

## geographic distance vs. FST ----

PairwiseFST@PairwiseFSTmatrix[,0] %>%
  as.data.frame %>%
  rownames_to_column("island_pond_seaseon") %>%
  separate(island_pond_seaseon,c("island", "pond", "season")) %>%
  group_by(season) %>% count()

# avoid pseudo-replication
samples_FST_spatial_smr_14_15_union <- paste0(union(gsub(pattern = "_s.r201.", "", grep("smr2014", colnames(PairwiseFST@PairwiseFSTmatrix), value = T)),
                                                    gsub(pattern = "_s.r201.", "", grep("smr2015", colnames(PairwiseFST@PairwiseFSTmatrix), value = T))), "_smr2015")
samples_FST_spatial_smr_14_15_union <- samples_FST_spatial_smr_14_15_union[samples_FST_spatial_smr_14_15_union != "N-27A_smr2015"]
samples_FST_spatial_smr_15 <- colnames(PairwiseFST@PairwiseFSTmatrix[colnames(PairwiseFST@PairwiseFSTmatrix) %in% samples_FST_spatial_smr_14_15_union, colnames(PairwiseFST@PairwiseFSTmatrix) %in% samples_FST_spatial_smr_14_15_union])
FST_spatial <- PairwiseFST@PairwiseFSTmatrix[colnames(PairwiseFST@PairwiseFSTmatrix) %in% c(samples_FST_spatial_smr_15, gsub("smr2015", "smr2014", setdiff(samples_FST_spatial_smr_14_15_union, samples_FST_spatial_smr_15))), colnames(PairwiseFST@PairwiseFSTmatrix) %in% c(samples_FST_spatial_smr_15, gsub("smr2015", "smr2014", setdiff(samples_FST_spatial_smr_14_15_union, samples_FST_spatial_smr_15)))]

# load coordinates
all.covariates <- read.csv("all.covariates.csv", sep=",")
all.covariates$season <- ifelse(all.covariates$sample == "1","spr", ifelse(all.covariates$sample == "2","smr",NA))
all.covariates$samplename <- paste0(all.covariates$poolname, "_", all.covariates$season, all.covariates$year)

# subset covariates for available data
covariates<-all.covariates[all.covariates$samplename %in% colnames(FST_spatial),]

# get coordinates
coordinates <- covariates[c("samplename","latitude_corr", "longitude_corr")]
distance <- geodist(coordinates)
distance[upper.tri(distance, diag=TRUE)] <- NA
FST_spatial[upper.tri(FST_spatial, diag=TRUE)] <- NA
combined <- data.frame(distance = as.vector(as.matrix(distance)), FST = as.vector(FST_spatial/(1-FST_spatial)))

# dbMEM by RDA
FST.pc <- prcomp(as.dist(FST_spatial))
xy.dbmem <- dbmem(as.dist(distance))
geodist.rda<-rda(FST.pc$x, xy.dbmem)
RsquareAdj(geodist.rda)

## intra-island
islands <- unique(gsub(pattern = "-[0-9A]*_s.r201.", "-",colnames(FST_spatial))) # include "-" to differentiate islands when grep
data_islands <- vector(mode = "list")
for (i in islands) data_islands[[i]] <- FST_spatial[grep(paste0("^",i), colnames(FST_spatial)),grep(paste0("^",i), colnames(FST_spatial))]
mean <- mean(unlist(data_islands), na.rm = TRUE)

colnames(distance)<-colnames(FST_spatial)
rownames(distance)<-colnames(FST_spatial)
distance_islands <- vector(mode = "list")
for (i in islands) distance_islands[[i]] <- distance[grep(paste0("^",i), colnames(distance)),grep(paste0("^",i), colnames(distance))]
combined_islands <- data.frame(distance = as.vector(as.matrix(unlist(distance_islands))), FST = as.vector(unlist(data_islands)/(1-unlist(data_islands))))

# dbMEM by RDA for island N
FST.pc <- prcomp(as.dist(data_islands$`N-`)) # r-stats
xy.dbmem <- dbmem(as.dist(distance_islands$`N-`))
geodist.rda<-rda(FST.pc$x, xy.dbmem)
RsquareAdj(geodist.rda)

## temporal analysis ----

# get pond names
ponds <- unique(gsub("_s.r20..", "", colnames(PairwiseFST@PairwiseFSTmatrix)))

# aggregate samples from each pond
FSTdata_ponds <- vector(mode = "list")
for (i in ponds) FSTdata_ponds[[i]] <- as.data.frame(PairwiseFST@PairwiseFSTmatrix[grep(paste0("^",i,"_"), colnames(PairwiseFST@PairwiseFSTmatrix)),grep(paste0("^",i,"_"), colnames(PairwiseFST@PairwiseFSTmatrix))])

# relative age from 2014
for (i in 1:length(FSTdata_ponds)) {
  if (nrow(FSTdata_ponds[[i]])>1){
    FSTdata_ponds[[i]]["year"] <- gsub("^.*_s.r", "", rownames(FSTdata_ponds[[i]])) 
    FSTdata_ponds[[i]]["season"] <- gsub("201.", "", gsub("^.*_", "", rownames(FSTdata_ponds[[i]])))
    FSTdata_ponds[[i]]["year"] <- as.numeric(as.data.frame(FSTdata_ponds[[i]])$year) - 2014
    FSTdata_ponds[[i]]["season"] <- ifelse(FSTdata_ponds[[i]]["season"]=="spr", 0, 0.2)
  }
}

dict <- list("0" = "spr2014", "0.2" = "smr2014",
             "1" = "spr2015", "1.2" = "smr2015",
             "2" = "spr2016", "2.2" = "smr2016",
             "3" = "spr2017", "3.2" = "smr2017",
             "4" = "spr2018", "4.2" = "smr2018")

# function to get absolute differences between samples
custom_fun <- function(x, y) {
  z <- abs(x - y)
  return(z)
}

# data frame for absolute differences between samples
delta.frame <- data.frame()

# calculate delta between all samples
for (i in 1:length(FSTdata_ponds)) {
  if (nrow(FSTdata_ponds[[i]])>1){
    # delta age calculation
    age_matrix <- outer(as.numeric(unlist(FSTdata_ponds[[i]]["year"])) + as.numeric(unlist(FSTdata_ponds[[i]]["season"])),as.numeric(unlist(FSTdata_ponds[[i]]["year"])) + as.numeric(unlist(FSTdata_ponds[[i]]["season"])),custom_fun)
    
    # add to data frame
    delta.frame <- rbind(delta.frame, cbind(age_matrix[lower.tri(age_matrix)], FSTdata_ponds[[i]][1:(length(FSTdata_ponds[[i]])-2)][lower.tri(FSTdata_ponds[[i]][1:(length(FSTdata_ponds[[i]])-2)])]))
  }
}

## Prepare host data for comparison ----

# load data
samplenames_magna <- system("bcftools query -l Host_VCF.vcf", intern = TRUE)
pooldata_magna <-vcf2pooldata("Host_VCF.vcf", poolsizes=rep(100,length(samplenames_magna)), poolnames=samplenames_magna)

# calculate Pairwise FST
FST_data_magna <- compute.pairwiseFST(pooldata_magna)
FST_data_magna@PairwiseFSTmatrix[FST_data_magna@PairwiseFSTmatrix < 0 & is.finite(FST_data_magna@PairwiseFSTmatrix)] <- 0
FST_data_magna <- FST_data_magna@PairwiseFSTmatrix

# focus on samples with parasite data
FST_spatial_magna <- FST_data_magna[colnames(FST_data_magna) %in% colnames(FST_spatial), colnames(FST_data_magna) %in% colnames(FST_spatial)]

## parasite FST vs. host FST

# focus on samples with host and parasite data
FST_spatial <- FST_spatial[colnames(FST_spatial_magna), colnames(FST_spatial_magna)]

# aggregate FST from both
FST_spatial_magna[upper.tri(FST_spatial_magna, diag=TRUE)] <- NA
FST_combined <- data.frame(parasite = as.vector(FST_spatial), host= unlist(as.vector(FST_spatial_magna)))
FST_combined$host <- FST_combined$host/(1-FST_combined$host)
FST_combined$parasite <- FST_combined$parasite/(1-FST_combined$parasite)

# initialize an empty matrix (to be filled where intra_island = T) for plotting
n <- nrow(FST_spatial)

# all F; col/row names are island ID
result_matrix <- matrix(FALSE, nrow = n, ncol = n, dimnames = list(gsub("-.*", "", rownames(FST_spatial)), gsub("-.*", "", colnames(FST_spatial))))

# fill in the matrix with TRUE where row name matches column name
for (i in 1:n) {
  for (j in 1:n) {
    if (rownames(result_matrix)[i] == colnames(result_matrix)[j]) {
      result_matrix[i, j] <- TRUE
    }
  }
}

# data frame with parasite and host FST and (intra)island info
FST_combined <- data.frame(parasite = as.vector(FST_spatial), host= unlist(as.vector(FST_spatial_magna)), intra_island=unlist(as.vector(result_matrix)))
FST_combined$host <- FST_combined$host/(1-FST_combined$host)
FST_combined$parasite <- FST_combined$parasite/(1-FST_combined$parasite)

wilcox.test(FST_combined[FST_combined$intra_island == "TRUE",]$parasite, FST_combined[FST_combined$intra_island == "FALSE",]$parasite)
wilcox.test(FST_combined[FST_combined$intra_island == "TRUE",]$host, FST_combined[FST_combined$intra_island == "FALSE",]$host)

# transform data
FST_combined_long <- gather(FST_combined, condition, measurement, host:parasite, factor_key=TRUE)

# plot Inter- vs. Intra-Island data for both
facet_names <- c(
  `host` = "D. magna",
  `parasite` = "H. tvaerminnensis"
)

Fig2A <- ggplot(FST_combined, aes(parasite, host)) + geom_point(aes(col = intra_island)) + theme_bw() +
  #stat_smooth(method='loess',se=T) 
  stat_smooth(method='lm',se=T, col ="black") + ylab(expression("Host [pairwise"~paste(italic(F) [ST])~"/(1 - pairwise"~paste(italic(F) [ST])~")]")) + xlab(expression("Parasite [pairwise"~paste(italic(F) [ST])~"/(1 - pairwise"~paste(italic(F) [ST])~")]")) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(labels=c("Inter-Island", "Intra-Island"), values = c("#eb6223", "#522888")) +
  geom_xsideboxplot(aes(y = intra_island, col = intra_island), orientation = "y", notch = T, show.legend = F) +
  geom_ysideboxplot(aes(x = intra_island, col = intra_island), orientation = "x", notch = T, show.legend = F) +
  scale_xsidey_discrete() +
  scale_ysidex_discrete() +
  theme(ggside.axis.ticks = element_blank(), ggside.axis.text = element_blank(), ggside.panel.scale = 0.2) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5)))

# dbMEM by RDA
FST.pc <- prcomp(as.dist(FST_spatial))
FST.dbmem <- dbmem(as.dist(FST_spatial_magna))
FST.rda<-rda(FST.pc$x, FST.dbmem)
RsquareAdj(FST.rda)

## Host temporal analysis

# get pond names
ponds_magna <- unique(gsub("_s.r201.", "", colnames(FST_data_magna)))

# aggregate samples from each pond
FSTdata_ponds_magna <- vector(mode = "list")
for (i in ponds_magna) FSTdata_ponds_magna[[i]] <- as.data.frame(FST_data_magna[grep(paste0("^",i,"_"), colnames(FST_data_magna)),grep(paste0("^",i,"_"), colnames(FST_data_magna))])

# relative age from 2014
for (i in 1:length(FSTdata_ponds_magna)) {
  if (nrow(FSTdata_ponds_magna[[i]])>1){
    FSTdata_ponds_magna[[i]]["year"] <- gsub("^.*_s.r", "", rownames(FSTdata_ponds_magna[[i]])) 
    FSTdata_ponds_magna[[i]]["season"] <- gsub("201.", "", gsub("^.*_", "", rownames(FSTdata_ponds_magna[[i]])))
    FSTdata_ponds_magna[[i]]["year"] <- as.numeric(as.data.frame(FSTdata_ponds_magna[[i]])$year) - 2014
    FSTdata_ponds_magna[[i]]["season"] <- ifelse(FSTdata_ponds_magna[[i]]["season"]=="spr", 0, 0.2)
  }
}

# data frame for absolute differences between samples
delta.frame_magna <- data.frame()

# calculate delta between all samples
for (i in 1:length(FSTdata_ponds_magna)) {
  if (nrow(FSTdata_ponds_magna[[i]])>1){
    # delta age calculation
    age_matrix <- outer(as.numeric(unlist(FSTdata_ponds_magna[[i]]["year"])) + as.numeric(unlist(FSTdata_ponds_magna[[i]]["season"])),as.numeric(unlist(FSTdata_ponds_magna[[i]]["year"])) + as.numeric(unlist(FSTdata_ponds_magna[[i]]["season"])),custom_fun)
    
    # add to data frame
    delta.frame_magna <- rbind(delta.frame_magna, cbind(age_matrix[lower.tri(age_matrix)], FSTdata_ponds_magna[[i]][1:(length(FSTdata_ponds_magna[[i]])-2)][lower.tri(FSTdata_ponds_magna[[i]][1:(length(FSTdata_ponds_magna[[i]])-2)])]))
  }
}

# plot pairwise FST vs. delta years
delta.frame$species <- "parasite"
delta.frame_magna$species <- "host"
delta.frame.comb <- rbind(delta.frame, delta.frame_magna)
delta.frame.comb$species <- factor(delta.frame.comb$species, levels = c("parasite", "host"))
delta.frame.comb$V2 <- as.numeric(delta.frame.comb$V2)
delta.frame.comb[delta.frame.comb$V2 < 0,]$V2 <- 0
delta.frame.comb$V2 <- delta.frame.comb$V2 + 0.001

Fig2B <- delta.frame.comb %>%
  ggplot(aes(x = V1, y = log(V2/(1-V2)), color = species, size = species)) +
  #ggplot(aes(x = V1, y = V2, color = species, size = species)) +
  geom_point(alpha = 0.3) + theme_bw() + 
  geom_smooth(data = delta.frame.comb[delta.frame.comb$species == "host",], method='nls', formula = y ~ a*b^(1/x), method.args = list(start = list(a = 1, b = 1)), se = F, lty = 2, size = 1) +
  geom_smooth(data = delta.frame.comb[delta.frame.comb$species == "parasite",], method='nls', formula = y ~ a*b^(1/x), method.args = list(start = list(a = 1, b = 1)), se = F, size = 1) +
  labs(x=expression("Delta years (="~paste(delta)~"t)"), y= expression("Pairwise"~paste(italic(F) [ST])~"/(1 - pairwise"~paste(italic(F) [ST])~") (ln)")) +
  #labs(x=expression("Delta years (="~paste(delta)~"t)"), y= expression("Pairwise"~paste(italic(F) [ST]))) +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), legend.position = "bottom") +
  scale_color_manual(name="", values=c("steelblue4", "#ADFF2F"), labels=c(expression(italic("H. tvaerminnensis")), expression(italic("D. magna")))) +
  scale_size_manual(name = "", values = c(2,1), labels=c(expression(italic("H. tvaerminnensis")), expression(italic("D. magna")))) + 
  scale_x_continuous(breaks = seq(0, 10, by = 2))

ggsave(filename = "FST_host_parasite_box_temporal_trafo_col.pdf", plot_grid(Fig2A, Fig2B, labels = c('A', 'B')), width = 10, height = 5)
