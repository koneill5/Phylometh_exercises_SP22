library(geiger)
library(rotl)
library(ape)
library(phytools)

#Tree and data are coming from a mix of Brandley et al., 2008 and Luke J. Harmon :https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/squamate.phy
tree <- read.tree("squamate.phy")
plot(tree)

discrete.data.old <- read.csv("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/squamate.phy", stringsAsFactors = FALSE) #I don't think this will work for the HW, so don't use

discrete.data <- read.csv("brandley_table.csv")

#Getting discrete data by counting species with zero-length fore- and hind limbs as "limbess"
limbless <- as.numeric(discrete.data[,"FLL"]==0 & discrete.data[,"HLL"]==0)
sum(limbless)

#match names:
nn <- discrete.data[,1]
nn2 <- sub("","_",names)
names(limbless) <- nn2

cleandata <- treedata(tree,limbless)

model <- fitDiscrete(cleandata$phy, cleandata$data)
cleaned.discrete.phydata <- phangorn::phyDat(cleandata)



