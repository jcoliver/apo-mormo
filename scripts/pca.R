# Principal Components Analysis
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-26

################################################################################
# SUMMARY
# * Reads in genetic data from genind object; also reads in locality information
# for plot coloring purposes
# * Performs Principal Components Analysis (PCA) on allele frequency data
# * Plots first three axes of PCA (1 & 2 in one plot, 2 & 3 in a second plot)

################################################################################
# SETUP

# Load dependancies
#install.packages("adegenet")
library("adegenet")

# Establish data files
genind.file = "output/genind-object.RData"
localities.file = "output/reconciled-localities.RData"
plot.out.file = "output/pca-plots.pdf"

# Read in genind file
load(file = genind.file)

# Extract allele frequency data, replacing NA values with mean allele 
# frequencies (mean across individuals of all populations)
allele.data <- tab(x = apo.str.genind, freq = TRUE, NA.method = "mean")

################################################################################
# PCA

# Run pca, skipping the scaling step (scale = FALSE), because we already scaled
# allele frequencies via tab call above.
allele.pca <- dudi.pca(df = allele.data, 
                       center = TRUE, 
                       scale = FALSE, 
                       scannf = FALSE, 
                       nf = 3)

# Create plots
# Load localities (for coloring north & south populations differently)
load(file = localities.file)
pop.nums <- pop(apo.str.genind)

south.pop.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
south.pop.nums <- localities$pop.number[which(localities$pop.name %in% south.pop.names)]


group2 <- which(localities$pop.name %in% group2.names)


pop.vector <- pop(apo.str.genind)


col.vector <- funky(15)

point.cols <- rep(x = "#B21212FF", times = nrow(localities)) # a dark red for north populations
point.cols[group2] <- "#12B2B2FF" # "cadetblue" for south populations


# Plot first two axes
xaxis <- 1
yaxis <- 2
s.class(dfxy = allele.pca$li, 
        fac = pop(apo.str.genind), 
        xax = xaxis, 
        yax = yaxis, 
        col = funky(15),
        clabel = 0.6)
add.scatter.eig(w = allele.pca$eig[1:50], 
                nf = 3, 
                xax = xaxis, 
                yax = yaxis, 
                ratio = 0.3)

# Plot second and third axes
xaxis <- 2
yaxis <- 3
s.class(dfxy = allele.pca$li, 
        fac = pop(apo.str.genind), 
        xax = xaxis, 
        yax = yaxis, 
        col = funky(15),
        clabel = 0.6)
add.scatter.eig(w = allele.pca$eig[1:50], 
                nf = 3, 
                xax = xaxis, 
                yax = yaxis, 
                ratio = 0.3)


################################################################################
# PCoA

# Produces results ~identical to PCA, above
allele.pco <- dudi.pco(d = dist(allele.data), scannf = FALSE, nf = 3)

# Plot first two axes
xaxis <- 1
yaxis <- 2
s.class(dfxy = allele.pco$li, 
        fac = pop(apo.str.genind), 
        xax = xaxis, 
        yax = yaxis, 
        col = funky(15),
        clabel = 0.6)
add.scatter.eig(w = allele.pco$eig[1:50], 
                nf = 3, 
                xax = xaxis, 
                yax = yaxis, 
                ratio = 0.3)

