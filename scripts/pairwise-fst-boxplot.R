# Boxplot of pairwise Fst values
# Jeffrey C. Oliver
# jcoliver@email.arizona.edu
# 2016-09-19

################################################################################
# SUMMARY
# Reads in pairwise Fst matrix and creates a boxplot for each population

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

plot.file = paste0("output/fst-boxplot.pdf", sep = "")

################################################################################
# DATA PREP
# Read in data & convert to long, re-level pop.name so populations are 
# displayed in plot in decreasing latitude

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


################################################################################
# PLOT

pdf(file = plot.file, useDingbats = FALSE)
  pairwise.plot <- ggplot(data = fst.long, aes(pop.name, fst)) + 
    geom_boxplot() +
  #  geom_jitter(position = position_dodge(width = 1)) +
    xlab(label = "Population") +
    ylab(label = "Fst") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  print(pairwise.plot)
dev.off()
