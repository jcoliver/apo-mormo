# Apodemia mormo isolation by distance plot
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-19

#install.packages("adegenet")
library("adegenet")
source(file = "functions/apodemia-functions.R")

################################################################################
# Load files & prep data
apo.str.genind <- read.structure(file = "data/Apodemia_0.9-noDockweiler-GNPSK.str",
                                 n.ind = 102,
                                 n.loc = 4057,
                                 onerowperind = TRUE,
                                 col.lab = 1,
                                 col.pop = 2,
                                 col.others = 0,
                                 row.marknames = 0,
                                 NA.char = "-9",
                                 sep = "\t")
# Pairwise Fst matrix, see ibd.R
load(file = "output/pairwise-fst-Apodemia_0.9-noDockweiler-GNPSK.RData")

localities <- FormatLocalities(file = "data/Apo_localities.txt",
                               genind = apo.str.genind, 
                               omit = c("Dockweiler", "GrasslandsNPSK"))

geo.dist <- GeoDistances(localities = localities)

# Transform distances for plot
p.fst <- pairwise.fst/(1 - pairwise.fst)
geo.dist <- log(x = geo.dist, base = 10)

################################################################################
# Identify south populations for coloring each type of comparison
# dot.colors will ultimately have three possible values:
# "white" for north-north
# "black" for south-south
# "red" for north-south
# Note there is some funkiness because the distance matrices are have 
# number of rows & columns equal to the number of localities minus one;
# these distance matrices have columns corresponding to populations 1 through 
# n - 1 (where n is the number of populations) and rows corresponding to 
# populations 2 through n.
south.pop.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
south.pop.numbers <- localities$pop.number[which(localities$pop.name %in% south.pop.names)]
dot.colors <- matrix(data = "purple", nrow = nrow(localities) - 1, ncol = nrow(localities) - 1)
rownames(dot.colors) <- localities$pop.number[2:length(localities$pop.number)]
colnames(dot.colors) <- localities$pop.number[1:(length(localities$pop.num) - 1)]

# Loop over row/columns to categorize comparison and assign appropriate color
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
  }
}

# We only want the elements from the lower triangle of the matrix (diagonal included)
dot.tri <- lower.tri(x = dot.colors, diag = TRUE)
dot.colors.lower <- dot.colors[dot.tri]

# Plot to a file
date.filename <- format(Sys.Date(), "%Y-%m-%d")
pdf(file = paste0("output/Apo-ibd-plot-", date.filename, ".pdf"), useDingbats = FALSE)
plot(x = geo.dist, y = p.fst, col = "black", pch = 21, bg = dot.colors.lower, xlab = "Log(distance)", ylab = "Fst/1 - Fst")
abline(lm(p.fst ~ geo.dist))
dev.off()
