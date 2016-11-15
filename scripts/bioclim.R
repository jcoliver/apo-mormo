# Exploratory analysis of climate variable predictors
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-11-02

rm(list = ls())
################################################################################
# Much of this comes from http://rspatial.org/sdm/rst/4_sdm_envdata.html
# Data are from http://www.worldclim.org/bioclim

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

# Check for colinearity among climate variables
bioclim.scaled <- scale(x = bioclim.vars)
r.matrix <- matrix(nrow = ncol(bioclim.scaled), ncol = ncol(bioclim.scaled))
colnames(r.matrix) <- rownames(r.matrix) <- colnames(bioclim.scaled)
for (i in 1:(ncol(bioclim.scaled) - 1)) {
  for (j in (i+1):ncol(bioclim.scaled)) {
    bioclim.model <- lm(bioclim.scaled[, i] ~ bioclim.scaled[, j])
    bioclim.summary <- summary(bioclim.model)
    bioclim.r.squared <- bioclim.summary$r.squared
    r.matrix[j, i] <- sqrt(bioclim.r.squared)
  }
}

# Skip because PCA effectively does this?
# Could still think about this, as it would reduce the number of variables
# TODO: calculate r^2 among all bioclim var pairs; drop secondaries that have 
# r^2 > 0.7...
fit <- prcomp(x = bioclim.vars, scale. = TRUE)
summary(fit)
fit$rotation[, 1:3]
x.lab <- paste0("PC 1 (", round(fit$sdev[1]^2/sum(fit$sdev^2) * 100, 0), "%)")
y.lab <- paste0("PC 2 (", round(fit$sdev[2]^2/sum(fit$sdev^2) * 100, 0), "%)")

# Points
plot(x = fit$x[, 1], y = fit$x[, 2], xlab = x.lab, ylab = y.lab)

# Pop. names
plot(x = fit$x[, 1], y = fit$x[, 2], xlab = x.lab, ylab = y.lab, type = "n")
text(x = fit$x[, 1], y = fit$x[, 2], labels = localities$pop.name)

# Just considering four variables:
# 1 Mean temperature
# 5 Max temperature
# 6 Min temperature
# 12 Mean precipitation
small.vars <- bioclim.vars[, c("bio1_12", "bio5_12", "bio6_12", "bio12_12")]
fit.small <- prcomp(x = small.vars, scale. = TRUE)
summary(fit.small)

x.lab <- paste0("PC 1 (", round(fit.small$sdev[1]^2/sum(fit.small$sdev^2) * 100, 0), "%)")
y.lab <- paste0("PC 2 (", round(fit.small$sdev[2]^2/sum(fit.small$sdev^2) * 100, 0), "%)")

# Points
plot(x = fit.small$x[, 1], y = fit.small$x[, 2], xlab = x.lab, ylab = y.lab)

# Pop. names
plot(x = fit.small$x[, 1], y = fit.small$x[, 2], xlab = x.lab, ylab = y.lab, type = "n")
text(x = fit.small$x[, 1], y = fit.small$x[, 2], labels = localities$pop.name)
