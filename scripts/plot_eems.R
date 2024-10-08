# Script for plotting output of eems software

library(reemsplots2)
library(ggplot2)
library(sf)
library(ggsn)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(rnaturalearth)
library(ggmap)
library(cowplot)

# plot data from eems
plots.Dm <- make_eems_plots("Magna-run/", longlat = TRUE, add_demes = F, eems_colors = rev(c("#ddb877", "#e3c692", "#eedbbb", "#f8f0e3", "#FFFFFF", "#d2dfe7", "#90b0c3", "#4e819f", "#226288")))
plots.Ht <- make_eems_plots("Ht-run/", longlat = TRUE, add_demes = F, eems_colors = rev(c("#ddb877", "#e3c692", "#eedbbb", "#f8f0e3", "#FFFFFF", "#d2dfe7", "#90b0c3", "#4e819f", "#226288")))

# shapefile from Finish environment institute
uso <- st_read("data/metadata_metapopulation/ranta10meret/meri10.shp")

# transform Finnland to international coordinates
uso_tr <- st_transform(uso, "EPSG:4326")
#uso_tr <- st_transform(uso, "EPSG:3395")

# get map
map_tr <- ggplot() +
  geom_sf(data = uso_tr[7], aes(color = "red")) +
  xlim(c(23.242983, 23.263650)) +
  ylim(c(59.813667, 59.833717))

# transform sf object to spatial (sp) object to be able to map it as a layer in ggplot
uso_tr_sp <- as(uso_tr, "Spatial")

# load coordinates of populations
daphnia <- read.table("Magna_5.coord", quote="\"", comment.char="")
hami <- read.table("Ht_5.coord", quote="\"", comment.char="")

# island labels
coord_islands <- data.frame(island=c("K","N","LON","LONA", "G", "LA", "M", "FSS", "FS", "SK", "LG"), lat=c(59.82377, 59.823, 59.81994, 59.81987, 59.8154, 59.8265, 59.8220, 59.828, 59.828, 59.8305, 59.818), lon=c(23.25273, 23.26072, 23.25067,23.25090, 23.248, 23.242, 23.246, 23.249, 23.2505, 23.255, 23.243))

# Host figure
FigDm <- plots.Dm$mrates01 +
  theme_bw() +
  ggspatial::layer_spatial(data = uso_tr_sp[7], color = "black", fill = NA, lwd = 0.5) +
  scale_x_continuous(limits = c(23.24, 23.2675), label = abs) +
  scale_y_continuous(limits = c(59.814, 59.8335), label = abs) +
  geom_point(aes(x=V1, y=V2), data=daphnia, shape=16, size=2, col="chartreuse3") +
  theme(axis.text = element_text(size=12,colour = "black"), axis.title = element_text(size=14), legend.key.size=unit(1.5, "cm"), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.key = element_rect(fill = "transparent"))+theme(legend.key=element_rect(color="transparent"), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scalebar(dist = 250, transform = T, dist_unit = "m", model = "WGS84", y.min = 59.8146, y.max = 59.8165, x.min = 23.252, x.max = 23.265, st.size = 4, location = "bottomleft",box.fill = c("black", "black"), st.dist = 0.3) +
  labs(x="\n Longitude East", y="Latitude North \n") +
  geom_text_repel(data=coord_islands, aes(x=lon, y=lat, label=island), size = 4, max.overlaps = Inf,
                  #box.padding = 0.35,
                  #segment.curvature = -0.1,
                  #segment.ncp = 5,
                  #segment.angle = 20,
                  #segment.color = "gray24",
                  segment.alpha=.7,
                  fontface = "bold",
                  nudge_x      = 0.0015,
                  min.segment.length = Inf)

# Parasite figure
FigHt <- plots.Ht$mrates01 +
  theme_bw() +
  ggspatial::layer_spatial(data = uso_tr_sp[7], color = "black", fill = NA, lwd = 0.5) +
  scale_x_continuous(limits = c(23.24, 23.2675), label = abs) +
  scale_y_continuous(limits = c(59.814, 59.8335), label = abs) +
  geom_point(aes(x=V1, y=V2), data=hami, shape=16, size=2, col="chartreuse3") +
  theme(axis.text = element_text(size=12,colour = "black"), axis.title = element_text(size=14), legend.key.size=unit(1.5, "cm"), legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  theme(legend.key = element_rect(fill = "transparent"))+theme(legend.key=element_rect(color="transparent"), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scalebar(dist = 250, transform = T, dist_unit = "m", model = "WGS84", y.min = 59.8146, y.max = 59.8165, x.min = 23.252, x.max = 23.265, st.size = 4, location = "bottomleft",box.fill = c("black", "black"), st.dist = 0.3) +
  labs(x="\n Longitude East", y="Latitude North \n") +
  geom_text_repel(data=coord_islands, aes(x=lon, y=lat, label=island), size = 4, max.overlaps = Inf,
                  #box.padding = 0.35,
                  #segment.curvature = -0.1,
                  #segment.ncp = 5,
                  #segment.angle = 20,
                  #segment.color = "gray24",
                  segment.alpha=.7,
                  fontface = "bold",
                  nudge_x      = 0.0015,
                  min.segment.length = Inf)

# Finland map for inset 
worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf')
finland <- worldmap[worldmap$name == 'Finland',]
FigFin <-ggplot() + geom_sf(data = finland) + theme_classic()+
  geom_point(aes(x=23.252, y=59.8285), alpha=.5,color="black", fill = "red", shape=22, size=3)+
  xlab("") + ylab("") + theme_inset() +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks=element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), plot.background = element_blank())

# Combine host and parasite figures
FigPatchWork <- FigDm + plot_layout(guides = 'keep') + 
  ggtitle(expression(italic("Daphnia magna"))) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="\n Longitude East", y="Latitude North \n") | 
  FigHt + guides(fill="none") + 
  ggtitle(expression(italic("Hamiltosporidium tvaerminnensis"))) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x="\n Longitude East", y = "")

# add Finland
ggdraw() +
  draw_plot(FigPatchWork) +
  draw_plot(FigFin, x = 0.835, y = 0.07, width = 0.3, height = 0.3) +
  draw_plot(FigFin, x = 0.3, y = 0.07, width = 0.3, height = 0.3)
