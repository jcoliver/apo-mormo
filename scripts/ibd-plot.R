# Apodemia mormo isolation by distance plot
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-19

################################################################################
# SUMMARY
# * Reads in fst and geographic distance matrices, as well as localities object
# * Creates scatterplot of distance vs. Fst, coloring points based on the 
# comparison (north-north, south-south, and north-south)
# * Creates scatterplot of distance vs. Fst, coloring all comparisons with 
# langei population different from remaining comparisons

################################################################################
# SETUP
# Establish data files and constants
pairwise.fst.file = "output/pairwise-fst.RData"
localities.file = "output/reconciled-localities.RData"
geo.dist.file = "output/pairwise-geo-dist.RData"

################################################################################
# DATA PREP

# Load pairwise Fst matrix
load(file = pairwise.fst.file)
# Want Fst/(1 - Fst) as our differentiation matrix
p.fst <- pairwise.fst/(1 - pairwise.fst)

# Load localities
load(file = localities.file)

# Load geographic distance matrix
load(file = geo.dist.file)
# Want log-transformed geographic distances, per Rousset 1997
# Rousset does not explicitly mention the base of the log function, 
# but http://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-6-13
# say base-10, so we go with that
geo.dist <- log(x = geo.dist, base = 10)

################################################################################
# IBD PLOT

# Identify south populations for coloring each type of comparison
# dot.colors will ultimately have three possible values:
# "white" for north-north
# "black" for south-south
# "red" for north-south
#
# Identify comparions with langei population, using triangles for those 
# comparions
#
# Note there is some funkiness because the distance matrices are have 
# number of rows & columns equal to the number of localities minus one;
# these distance matrices have columns corresponding to populations 1 through 
# n - 1 (where n is the number of populations) and rows corresponding to 
# populations 2 through n.
south.pop.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego", "Dockweiler")
south.pop.numbers <- localities$pop.number[which(localities$pop.name %in% south.pop.names)]
dot.colors <- matrix(data = "purple", nrow = nrow(localities) - 1, ncol = nrow(localities) - 1)
rownames(dot.colors) <- localities$pop.number[2:length(localities$pop.number)]
colnames(dot.colors) <- localities$pop.number[1:(length(localities$pop.number) - 1)]

langei.pop.name <- c("langei")
langei.pop.number <- localities$pop.number[which(localities$pop.name %in% langei.pop.name)]
dot.shapes <- matrix(data = 0, nrow = nrow(localities) - 1, ncol = nrow(localities) - 1)
rownames(dot.shapes) <- localities$pop.number[2:length(localities$pop.number)]
colnames(dot.shapes) <- localities$pop.number[1:(length(localities$pop.number) - 1)]

# Loop over row/columns to categorize comparison and assign appropriate color and shape
for (i in 1:nrow(dot.colors)) {
  for (j in 1:i) {
    locality.i <- i + 1
    locality.j <- j

    locality.i.number <- localities$pop.number[locality.i]
    locality.j.number <- localities$pop.number[locality.j]
    
    d <- "red"
    if (locality.i.number %in% south.pop.numbers
        && locality.j.number %in% south.pop.numbers) { # both south
      d <- "black"
    } else if (!(locality.i.number %in% south.pop.numbers)
               && !(locality.j.number %in% south.pop.numbers)) { # both NOT south, so north-north
      d <- "white"
    }
    dot.colors[i, j] <- d
    
    s <- 21
    if (locality.i.number %in% langei.pop.number
        || locality.j.number %in% langei.pop.number) {
      s <- 25
    }
    dot.shapes[i, j] <- s
  }
}

# We only want the elements from the lower triangle of the matrix (diagonal included)
dot.tri <- lower.tri(x = dot.colors, diag = TRUE)
dot.colors.lower <- dot.colors[dot.tri]
dot.shapes.lower <- dot.shapes[dot.tri]

# Plot to a file
pdf(file = paste0("output/Apo-ibd-plot.pdf"), useDingbats = FALSE)
  plot(x = geo.dist, 
       y = p.fst, 
       col = "black", 
       # pch = 21, 
       pch = dot.shapes.lower,
       bg = dot.colors.lower, 
       xlab = "Log(distance)", 
       ylab = "Fst/1 - Fst")
  abline(lm(p.fst ~ geo.dist))
  legend("topleft", 
         legend = c("North-North", "South-South", "North-South"), 
         col = "black", 
         pt.bg = c("white", "black", "red"),
         pch = 21,
         cex = 0.8)
dev.off()

################################################################################
# Plot with langei as different color dot
# "white" for most
# "black" for langei-north & langei-south
langei.pop.name <- "langei"
langei.pop.number <- localities$pop.number[which(localities$pop.name %in% langei.pop.name)]
langei.colors <- matrix(data = "purple", nrow = nrow(localities) - 1, ncol = nrow(localities) - 1)
rownames(langei.colors) <- localities$pop.number[2:length(localities$pop.number)]
colnames(langei.colors) <- localities$pop.number[1:(length(localities$pop.num) - 1)]

# Loop over row/columns to categorize comparison and assign appropriate color
for (i in 1:nrow(langei.colors)) {
  for (j in 1:i) {
    locality.i <- i + 1
    locality.j <- j
    
    locality.i.number <- localities$pop.number[locality.i]
    locality.j.number <- localities$pop.number[locality.j]
    
    d <- "white"
    if ((locality.i.number %in% langei.pop.number
        && !(locality.j.number %in% langei.pop.number))
        || (locality.j.number %in% langei.pop.number
          && !(locality.i.number %in% langei.pop.number))){
      d <- "black"
    }
    langei.colors[i, j] <- d
  }
}

# Pull out lower triangle of the matrix (diagonal included)
langei.tri <- lower.tri(x = langei.colors, diag = TRUE)
langei.colors.lower <- langei.colors[langei.tri]

# Plot to a file
pdf(file = paste0("output/Apo-ibd-plot-langei.pdf"), useDingbats = FALSE)
  plot(x = geo.dist, 
       y = p.fst, 
       col = "black", 
       pch = 21, 
       bg = langei.colors.lower, 
       xlab = "Log(distance)", 
       ylab = "Fst/1 - Fst")
  abline(lm(p.fst ~ geo.dist))
dev.off()
