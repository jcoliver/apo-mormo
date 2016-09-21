# Boxplot of pairwise Fst values
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-19

################################################################################
# SUMMARY

################################################################################
# SETUP
# Load dependancies & establish data and output files

# Load dependancies
#install.packages("tidyr")
library("tidyr")
#install.packages("ggplot2")
library("ggplot2")
#source(file = "functions/apodemia-functions.R")

# Establish data & output files
pairwise.fst.file = "output/pairwise-fst.RData"
localities.file = "output/reconciled-localities.RData"

plot.file = "output/"

################################################################################
# DATA PREP

# Read in Fst data and convert to long format
load(file = pairwise.fst.file)
fst.matrix <- as.matrix(pairwise.fst)
fst.df <- as.data.frame(fst.matrix)
fst.long <- gather(data = fst.df, key = pop.number, value = fst)
# Drop zeros
fst.long <- fst.long[which(fst.long$fst > 0), ]

# Read in localities files & get names of localities
load(file = localities.file)
fst.long <- merge(x = fst.long, y = localities, by = "pop.number")

# Re-level pop.name factor by increasing latitude
pop.name.levels <- localities$pop.name[order(localities$latitude, decreasing = TRUE)]
fst.long$pop.name <- factor(fst.long$pop.name, levels = pop.name.levels)

#boxplot(formula = fst ~ pop.number, data = fst.long)
ggplot(data = fst.long, aes(pop.name, fst)) +
  geom_boxplot() +
  xlab(label = "Pop.") +
  ylab(label = "Fst") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

