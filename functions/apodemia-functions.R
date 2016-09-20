# Functions for use in Apodemia project
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-06

################################################################################
#' Haversine formula for great circle distances
#' 
#' @usage Haversine(lat1, long1, lat2, long2)
#' 
#' @param lat1 numeric vector of latitude in decimal degrees
#' @param long1 numeric vector of longitude in decimal degrees
#' @param lat2 numeric vector of latitude in decimal degrees
#' @param long2 numeric vector of longitude in decimal degrees
#' 
#' @details Calculates the distance between two latitude / longitude 
#' coordinates using the Haversine formula for great circle distance. See  
#' http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
#' 
#' @return numeric distance in kilometers
#' 
#' @examples 
#' distance <- Haversine(lat1 = 32.23, long1 = -110.95, lat2 = 39.00, long2 = -114.31)
#' distance # 6508.232
Haversine <- function(lat1, long1, lat2, long2) {
  # Convert to radians
  lat1.rad <- (lat1 * pi) / 180
  long1.rad <- (long1 * pi) / 180
  lat2.rad <- (lat2 * pi) / 180
  long2.rad <- (long2 * pi) / 180

  earth.rad <- 6371 # km
  delta.lat <- lat1.rad - lat2.rad
  delta.long <- long1.rad - long2.rad
  a <- sin(delta.lat/2)^2 + cos(lat1.rad) * cos(lat2.rad) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = earth.rad * c
  return(d) # Distance in km
}

################################################################################
#' Haversine formula for great circle distances
#' 
FormatLocalities <- function(file, sep = "\t", genind, omit = c()) {
  localities <- read.table(file = file, sep = sep, header = TRUE)
  if (is.null(localities$pop.name)) {
    stop(paste0("File ", file, " is missing required 'pop.name' vector."))
  }
  
  if (length(omit) > 0) {
    localities <- localities[-which(localities$pop.name %in% omit), ]
  }
  
  genind.localities <- data.frame(pop.name = rownames(genind@tab), 
                                  pop.number = genind@pop)
  genind.localities$pop.name <- gsub("_[0-9A-Za-z]+", "", genind.localities$pop.name)
  genind.localities <- genind.localities[order(genind.localities$pop.name), ]
  genind.localities <- genind.localities[match(unique(genind.localities$pop.name), genind.localities$pop.name), ]
  
  # Sanity check
  if (nrow(localities) != nrow(genind.localities)) {
    stop(paste0("Different number of localities in file (", nrow(localities), ") and genind object (", nrow(genind.localities), ")"))
  }
  
  localities <- localities[, -which(colnames(localities) == "pop.number")]
  localities <- merge(x = localities, y = genind.localities, by = "pop.name")
  
  localities <- localities[order(localities$pop.number), ]
  rownames(localities) <- NULL
  return(localities)
}

################################################################################
#' Calculate pairwise geographic distance matrix
#' 
GeoDistances <- function(localities) {
  if (is.null(localities$latitude) || is.null(localities$longitude)) {
    stop("Missing at least one required vector (latitude, longitude) from localities data.frame.")
  }
  if (is.null(localities$pop.number)) {
    stop("Missing required vector 'pop.number' from localities data.frame.")
  }
  
  geo.dist <- matrix(data = NA, nrow = nrow(localities), ncol = nrow(localities))
  colnames(geo.dist) <- localities$pop.number[1:length(localities$pop.number)]
  rownames(geo.dist) <- localities$pop.number[1:length(localities$pop.number)]
  
  for (i in 1:(length(localities$latitude) - 1)) {
    lat1 <- localities$latitude[i]
    long1 <- localities$longitude[i]
    for (j in (i + 1):length(localities$latitude)) {
      lat2 <- localities$latitude[j]
      long2 <- localities$longitude[j]
      d <- Haversine(lat1 = lat1, long1 = long1, lat2 = lat2, long2 = long2)
      geo.dist[j, i] <- d
    }
  }
  geo.dist <- as.dist(geo.dist)
  return(geo.dist)
}

################################################################################
#' Generate an indicator matrix based on list of populations
#' 
IndicatorMatrix <- function(pop.set, localities) {
  if (is.null(pop.set) || length(pop.set) < 1) {
    stop("Object 'pop.set' empty or null; pop.set must include at least one population.")
  }
  
  if (is.null(localities$pop.name)) {
    stop("Required 'pop.name' vector missing from localities object.")
  }
  
  if (is.null(localities$pop.number)) {
    stop("Required 'pop.number' vector missing from localities object.")
  }
  
  to.mark <- which(localities$pop.name %in% pop.set)
  indicator.matrix <- matrix(data = NA, nrow = nrow(localities), ncol = nrow(localities))
  rownames(indicator.matrix) <- colnames(indicator.matrix) <- localities$pop.number
  for (i in 1:nrow(localities) - 1) {
    for (j in 1:nrow(localities)) {
      d <- 0
      if ((i %in% to.mark && !(j %in% to.mark)) 
          || (j %in% to.mark && !(i %in% to.mark))) {
        d <- 1
      }
      indicator.matrix[j, i] <- d
    }
  }
  indicator.matrix <- as.dist(indicator.matrix)
  return(indicator.matrix)
}