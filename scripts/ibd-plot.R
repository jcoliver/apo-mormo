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

south.pop.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")
south.pop.numbers <- as.integer(localities$pop.number[which(localities$pop.name %in% south.pop.names)])
dot.colors <- matrix(data = 0, nrow = nrow(localities), ncol = nrow(localities))
rownames(dot.colors) <- colnames(dot.colors) <- localities$pop.number
dot.colors <- as.dist(dot.colors)


for (i in 1:(nrow(localities) - 1)) {
  for (j in 1: nrow(localities)) {
    d <- 1
    if (localities$pop.number[i] %in% south.pop.numbers
        && localities$pop.number[j] %in% south.pop.numbers) {
      
    }
  }
}

north <- c()
south <- c()




south.names <- c("PointLoma", "WildhorseMeadows", "CampPendleton", "Borrego")

# Plotting distance vs. Fst (both transformed)
# Coloring comparisons with North vs. South black
cols <- north.south
cols[cols == 0] <- "white"
cols[cols == 1] <- "black"
date.filename <- format(Sys.Date(), "%Y-%m-%d")
pdf(file = paste0("output/Apo-ibd-ns-graph-", date.filename, ".pdf"), useDingbats = FALSE)
plot(x = geo.dist, y = p.fst, col = "black", pch = 21, bg = cols, xlab = "Log(distance)", ylab = "Fst/1 - Fst")
abline(lm(p.fst ~ geo.dist))
dev.off()

