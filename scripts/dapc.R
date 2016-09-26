# DAPC on Apodemia data
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-22

# INCOMPLETE - considerable variation in results among runs of `find.clusters`

# install.packages("adegenet")
library("adegenet")

load(file = "output/genind-object.RData")


# Doing some clustering, we see two groups
# See https://nescent.github.io/popgenInfo/DifferentiationSNP.html

# Also see http://diversityinlife.weebly.com/uploads/1/2/8/6/12864099/pop_structure_lecture.html

grp20 <- find.clusters(apo.str.genind, max.n.clust = 15, n.pca = 20, choose.n.clust = FALSE) 
dapc20 <- dapc(apo.str.genind, grp20$grp, n.pca = 20, n.da = 6) 
scatter(dapc20)

# need to do scalegen?

# Looking at between-rep runs of the optimal number of clusters
num.reps <- 30
optimal.K <- integer(length = num.reps)
print.freq <- 10
apo.scaled.genind <- scaleGen(x = apo.str.genind, NA.method = "mean")
apo.scaled.genind <- scaleGen()
for (i in 1:num.reps) {
  suppressMessages( # to avoid 'The "ward" method has been renamed..." message with every call
    groups <- find.clusters(apo.scaled.genind, max.n.clust = 20, n.pca = 20, choose.n.clust = FALSE)
  )
  optimal.K[i] <- as.integer(gsub(pattern = "K=", replacement = "", x = names(groups$stat)))
  if (i %% print.freq == 0) {
    cat("Replicate ", i, " of ", num.reps, " complete.\n", sep = "")
  }
}
hist(optimal.K)

# Some plotting options
library(reshape2)
library(ggplot2)
dat <- melt(dapc_4pop$posterior)
ggplot(dat, aes(x = Var1, y = value, fill = as.factor(Var2))) +
  geom_bar(stat = "identity", width = 1)