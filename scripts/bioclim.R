# Exploratory analysis of climate variable predictors
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-11-02

rm(list = ls())
################################################################################
# Much of this comes from http://rspatial.org/sdm/rst/4_sdm_envdata.html

localities.file = "output/reconciled-localities.RData"

# install.packages("raster")
library("raster")

# Download climate variables into RasterStack
# For high resolution (res = 0.5), only parts of the data are downloaded; 
# unfortunately, the cutoff between data sets is longitude -120, so two 
# downloads are necessary. While this is a east/west division, it effectively 
# splits along south/north division in our populations (south = east, north = 
# west).
bioclim.west <- getData("worldclim", download = TRUE, path = "bioclim", var = "bio", res = 0.5, lat = 40, lon = -121)
bioclim.east <- getData("worldclim", download = TRUE, path = "bioclim", var = "bio", res = 0.5, lat = 40, lon = -119)

# Pull out values for our localities
load(file = localities.file)
coords <- localities[, c("longitude", "latitude")]
bioclim.east.extract <- extract(x = bioclim.east, y = coords)
bioclim.west.extract <- extract(x = bioclim.west, y = coords)

# Join into one data frame
bioclim.vars <- bioclim.east.extract
bioclim.vars[is.na(bioclim.vars)] <- bioclim.west.extract[is.na(bioclim.vars)]
rownames(bioclim.vars) <- localities$pop.name

# Check for colinearity among climate variables?
# Skip because PCA effectively does this?
# Could still think about this, as it would reduce the number of variables
# TODO: calculate r^2 among all bioclim var pairs; drop secondaries that have 
# r^2 > 0.7...
# raster::pairs(x = bioclim.subset[, 1:4])
fit <- prcomp(x = bioclim.vars, scale = TRUE)
summary(fit)
x.lab <- paste0("PC 1 (", round(fit$sdev[1]/sum(fit$sdev) * 100, 0), "%)")
y.lab <- paste0("PC 2 (", round(fit$sdev[2]/sum(fit$sdev) * 100, 0), "%)")

# Points
plot(x = fit$x[, 1], y = fit$x[, 2], xlab = x.lab, ylab = y.lab)

# Pop. names
plot(x = fit$x[, 1], y = fit$x[, 2], xlab = x.lab, ylab = y.lab, type = "n")
text(x = fit$x[, 1], y = fit$x[, 2], labels = localities$pop.name)
