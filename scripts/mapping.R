# Mapping Apodemia localities
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-22

################################################################################
# SUMMARY
# Reads in locality latitude / longitude coordinates and plots points on a map.
# Populations colored by north / south distinction

################################################################################
# SETUP
# Load dependancies & establish input and output filenames

# Load dependancies
#install.packages("maps")
library(maps)
#install.packages("mapdata")
library(mapdata)

# Establish input and output filenames
localities.file <- "data/Apo_localities.txt"
date.filename <- date.filename <- format(Sys.Date(), "%Y-%m-%d")
map.out.file <- paste0("output/population-map-", date.filename, ".pdf")

################################################################################
# DATA PREP
# Read in data, remove GrasslandsNPSK population, and establish color vector to 
# color the points based on which STRUCTURE cluster they belong to (north or 
# south)

# Read in data
localities <- read.delim(file = localities.file, header = TRUE)

# Remove the GrasslandsNPSK locality
localities <- localities[-(which(localities$pop.name == "GrasslandsNPSK")), ]

# Identify those populations in the south cluster
group2.names <- c("PointLoma", "WildhorseMeadows", "Dockweiler", "CampPendleton", "Borrego")
group2 <- which(localities$pop.name %in% group2.names)

# Color vector based on north / south
point.cols <- rep(x = "#B21212FF", times = nrow(localities)) # a dark red for north populations
point.cols[group2] <- "#12B2B2FF" # "cadetblue" for south populations


################################################################################
# PLOT

pdf(file = map.out.file, useDingbats = FALSE)
  map(database = "worldHires", 
      xlim = c(-125, -115), 
      ylim = c(31, 41), 
      col = "#E7E7E7", 
      fill = TRUE)
  points(x = localities$longitude, 
         y = localities$latitude, 
         cex = 1.2, 
         pch = 21, 
         bg = point.cols, 
         col = "black")
dev.off()
