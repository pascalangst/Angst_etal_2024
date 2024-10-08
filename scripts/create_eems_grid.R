## create population grid (demes and edges) for eems

library(readxl)
library(ggplot2)

# get demes
coordinates <- read_excel("Pools_coordinates_2017_vers7b.xlsx")

# function to calculate Euclidean distances between two points
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}

# calculate distances between coordinates
dist_matrix <- matrix(0, nrow = nrow(coordinates), ncol = nrow(coordinates))
for (i in 1:nrow(coordinates)) {
  for (j in 1:nrow(coordinates)) {
    if (i != j) {
      dist_matrix[i, j] <- euclidean_distance(coordinates$latitude_corr[i], coordinates$longitude_corr[i], coordinates$latitude_corr[j], coordinates$longitude_corr[j])
    }
  }
}

# find the five closest neighbors for each deme (plus some random connections to immitate inter-island migration)
set.seed(12345)
closest_neighbors <- apply(dist_matrix, 1, function(x) order(x)[c(2:6, sample(c(7:(nrow(dist_matrix)/4), 7:(nrow(dist_matrix)/4), rep(NA, nrow(dist_matrix)/2)),1))])

# Create a data frame to store the connections
connections <- data.frame()
for (i in 1:nrow(coordinates)) {
  counter <- 0
  for (j in closest_neighbors[, i]) {
    counter <- counter + 1
    if (i > j & i %in% closest_neighbors[,closest_neighbors[counter, i]]) { # connections only once (avoid connections if there is one already)
      next
    }
    else {
    connections <- rbind(connections, data.frame(
      deme_1 = i,
      from_latitude_corr = coordinates$latitude_corr[i],
      from_longitude_corr = coordinates$longitude_corr[i],
      deme_2 = j,
      to_latitude_corr = coordinates$latitude_corr[j],
      to_longitude_corr = coordinates$longitude_corr[j]
    )) }
  }
}

# deme number for plot
dot_numbers <- data.frame(
  latitude_corr = coordinates$latitude_corr,
  longitude_corr = coordinates$longitude_corr,
  index = seq_len(nrow(coordinates))
)

# plot grid (demes and connections=edges)
p <- ggplot() +
  geom_point(data = coordinates, aes(x = longitude_corr, y = latitude_corr), color = "blue", size = 3) +
  geom_segment(data = connections, aes(x = from_longitude_corr, y = from_latitude_corr, xend = to_longitude_corr, yend = to_latitude_corr), color = "gray", size = 0.5) +
  labs(x = "longitude_corrgitude", y = "latitude_corritude", title = "Coordinate Grid with Connections") +
  theme_minimal() +
  geom_text(data = dot_numbers, aes(x = longitude_corr, y = latitude_corr, label = index), color = "black", size = 3, vjust = -0.5)

p

# save demes and edges (if valid connections)
connections <- connections[complete.cases(connections),]
write.table(coordinates[c("longitude_corr", "latitude_corr")], file="gridpath.demes", row.names=F, col.names=F)
write.table(connections[c("deme_1", "deme_2")], file="gridpath.edges", row.names=F, col.names=F)
