# Principal Components Analysis
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-26

# TODO: Add reference to tutorial pdfs for adegenet

################################################################################
# SUMMARY
# * Reads in genetic data from genind object; also reads in locality information
# for plot coloring purposes
# * Performs Principal Components Analysis (PCA) on allele frequency data
# * Plots first three axes of PCA (1 & 2 in one plot, 2 & 3 in a second plot)
# * Runs preliminary Principal Coordinates Analyses (PCoA), which appears to 
# produce results identical to PCA.

################################################################################
# SETUP
# Load dependancies, load data, and extract allele frequencies

# Load dependancies
#install.packages("adegenet")
library("adegenet")

# Establish data files
genind.file = "output/genind-object.RData"
localities.file = "output/reconciled-localities.RData"
plot.out.file = "output/pca-plots.pdf"

# Read in genind file
load(file = genind.file)

################################################################################
# DATA PREP
# Extract allele frequencies & set up color vector for coloring plots

# Extract allele frequency data, replacing NA values with mean allele 
# frequencies (mean across individuals of all populations)
allele.data <- tab(x = apo.str.genind, freq = TRUE, NA.method = "mean")

# Load localities (for coloring north & south populations differently)
load(file = localities.file)
pop.nums <- pop(apo.str.genind)

# Lots of work to give south populations range of red colors and north 
# populations a range of blue colors

# Identify north and south population numbers
south.pop.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
south.pop.nums <- localities$pop.number[which(localities$pop.name %in% south.pop.names)]
north.pop.nums <- localities$pop.number[-which(localities$pop.name %in% south.pop.names)]

# Generate palettes for north and south populations
south.color.endpoints <- c("#00c4c4", "#428282")
south.colors <- colorRampPalette(colors = south.color.endpoints)(length(south.pop.nums))
north.color.endpoints <- c("#DF0000", "#7C3737")
north.colors <- colorRampPalette(colors = north.color.endpoints)(length(north.pop.nums))
col.vector <- c(south.colors, north.colors)

# Use population number as element names in col.vector so we can re-order the 
# elements in col.vector based on levels in the pop.vector; necessary so the 
# populations are actually colored as we expect them to be
names(col.vector) <- c(south.pop.nums, north.pop.nums)
col.vector <- col.vector[levels(pop.nums)]

################################################################################
# PCA
# Run PCA, identify populations for plot coloring purposes, and plot first 
# three PCA axes (in two graphs)

# Run pca, skipping the scaling step (scale = FALSE), because we already scaled
# allele frequencies via tab call above.
allele.pca <- dudi.pca(df = allele.data, 
                       center = TRUE, 
                       scale = FALSE, 
                       scannf = FALSE, 
                       nf = 3)

# Create plots
# Open pdf device and send plots
pdf(file = plot.out.file, useDingbats = FALSE)
par(mfrow = c(2,1))

# Plot first two axes
xaxis <- 1
yaxis <- 2
s.class(dfxy = allele.pca$li, 
        fac = pop.nums, 
        xax = xaxis, 
        yax = yaxis, 
        col = col.vector,
        clabel = 0.6)
add.scatter.eig(w = allele.pca$eig[1:20], 
                nf = 3, 
                xax = xaxis, 
                yax = yaxis, 
                ratio = 0.3)

# Plot second and third axes
xaxis <- 2
yaxis <- 3
s.class(dfxy = allele.pca$li, 
        fac = pop.nums, 
        xax = xaxis, 
        yax = yaxis, 
        col = col.vector,
        clabel = 0.6)
add.scatter.eig(w = allele.pca$eig[1:20], 
                nf = 3, 
                xax = xaxis, 
                yax = yaxis, 
                ratio = 0.3, 
                posi = "topleft")
dev.off()

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

