# Mapping Apodemia localities
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-22

# INCOMPLETE

################################################################################
# SUMMARY

################################################################################
# SETUP

#install.packages("maps")
library(maps)
#install.packages("mapdata")
library(mapdata)

localities.file <- "data/Apo_localities.txt"
date.filename <- date.filename <- format(Sys.Date(), "%Y-%m-%d")
map.out.file <- paste0("output/population-map-", date.filename, ".pdf")

################################################################################
# DATA PREP

localities <- read.delim(file = localities.file, header = TRUE)

# Remove the GrasslandsNPSK locality
localities <- localities[-(which(localities$pop.name == "GrasslandsNPSK")), ]

# Different colors based on structure K = 2 results
group2.names <- c("PointLoma", "WildhorseMeadows", "Dockweiler", "CampPendleton", "Borrego")
group2 <- which(localities$pop.name %in% group2.names)

point.cols <- rep(x = "#B21212FF", times = nrow(localities)) # a dark red
point.cols[group2] <- "#12B2B2FF" # "cadetblue"

pdf(file = map.out.file, useDingbats = FALSE)
  map(database = "worldHires", xlim = c(-125, -115), ylim = c(31, 41), col = "lightgrey", fill = TRUE)
  points(x = localities$longitude, y = localities$latitude, cex = 1.2, pch = 21, bg = point.cols, col = "black")
dev.off()
