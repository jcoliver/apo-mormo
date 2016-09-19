# Mapping Apodemia localities
# Jeff Oliver
# 2016-01-25
# jcoliver@email.arizona.edu

setwd("~/Documents/Other Work/Apodemia/Apo-R")
localities <- read.delim(file = "data/Apo_localities.txt", header = TRUE)
#install.packages("maps")
#install.packages("mapdata")
library(maps)
library(mapdata)

point.cols <- rep(x = "black", times = nrow(localities))
point.cols[which(localities$pop.name == "langei")] <- "red"

map(database = "worldHires", xlim = c(-125, -105), ylim = c(30, 50), col = "lightgrey", fill = TRUE)
points(x = localities$longitude, y = localities$latitude, cex = 0.8, pch = 19, col = point.cols)

# Different colors based on structure K = 2 results
group2.names <- c("PointLoma", "WildhorseMeadows", "Dockweiler", "CampPendleton", "GrasslandsNPSK", "Borrego")
group2 <- which(localities$pop.name %in% group2.names)

point.cols <- rep(x = "#B21212FF", times = nrow(localities)) # x = "red"
point.cols[group2] <- "#12B2B2FF" # "cadetblue"

pdf(file = "output/Apo map.pdf", useDingbats = FALSE)
  map(database = "worldHires", xlim = c(-125, -105), ylim = c(30, 50), col = "lightgrey", fill = TRUE)
  points(x = localities$longitude, y = localities$latitude, cex = 0.8, pch = 21, bg = point.cols, col = "black")
dev.off()
