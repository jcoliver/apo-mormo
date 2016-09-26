# Principal Components Analysis
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-26

################################################################################
# SUMMARY

################################################################################
# SETUP

# Load dependancies
#install.packages("adegenet")
library("adegenet")

# Establish data files
str.file = "data/Apodemia_0.9-noDockweiler-GNPSK.str"
localities.file = "output/reconciled-localities.RData"

# Read in STRUCTURE file
apo.str.genind <- read.structure(file = str.file,
                                 n.ind = 102,
                                 n.loc = 4057,
                                 onerowperind = TRUE,
                                 col.lab = 1,
                                 col.pop = 2,
                                 col.others = 0,
                                 row.marknames = 0,
                                 NA.char = "-9",
                                 sep = "\t")

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

