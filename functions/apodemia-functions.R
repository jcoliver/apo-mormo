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

