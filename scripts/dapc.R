# DAPC on Apodemia data
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-22

# INCOMPLETE

# install.packages("adegenet")
library("adegenet")

load(file = "output/genind-object.RData")


# Doing some clustering, we see two groups
# See https://nescent.github.io/popgenInfo/DifferentiationSNP.html

grp <- find.clusters(apo.str.genind, max.n.clust = 10, n.pca = 20, choose.n.clust = FALSE) 
dapc1 <- dapc(apo.str.genind, grp$grp, n.pca = 20, n.da = 6) 
scatter(dapc1)
