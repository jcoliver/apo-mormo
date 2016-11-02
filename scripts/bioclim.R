# Exploratory analysis of climate variable predictors
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-11-02

rm(list = ls())
################################################################################
# install.packages("raster")
library("raster")

# Download climate variables into RasterStack
bioclim <- getData("worldclim", download = TRUE, path = "bioclim", var = "bio", res = 2.5)

# Pull out values for our localities


# Check for colinearity among climate variables

# From raster guide help on getData
# If name='worldclim' you must also provide arguments var, and a resolution res. 
# Valid variables names are ’tmin’, ’tmax’, ’prec’ and ’bio’.  Valid resolutions 
# are 0.5, 2.5, 5, and 10 (minutes of a degree).  In the case of res=0.5 , you 
# must also provide a lon and lat argument for a tile; for the lower resolutions 
# global data will be downloaded. In all cases there are 12 (monthly) files for 
# each variable except for ’bio’ which contains 19 files.
# getData('worldclim', var='tmin', res=0.5, lon=5, lat=45)
# getData('worldclim', var='bio', res=10)