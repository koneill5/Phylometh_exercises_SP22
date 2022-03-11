library(geiger)
library(rotl)
library(ape)
library(phytools)
library(phangorn)
library(tidyverse)

#Tree and data are coming from a mix of Brandley et al., 2008 and Luke J. Harmon :https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/squamate.phy
tree <- read.tree("squamate.phy")
plot(tree)

discrete.data <- read.csv("brandley_table.csv", stringsAsFactors = FALSE)

#Getting discrete data by counting species with zero-length fore- and hind limbs as "limbless"
limbless <- as.numeric(discrete.data[,"FLL"]==0 & discrete.data[,"HLL"]==0)
sum(limbless)

limbed.dat <- as.numeric(discrete.data[,"FLL"]>0 & discrete.data[,"HLL"]>0)
sum(limbed.dat)

#match names:
nn <- discrete.data[,1]
nn2 <- sub("","_",nn)
names(limbless) <- nn2

new.disc.dat <- cbind(discrete.data$Species,limbed.dat,limbless)
newdat2 <- gsub(" ","_",new.disc.dat[,1])

#need to have the species names as the rows ***There needs to be a space between the quotemarks 
rownames(discrete.data) <- gsub(" ","_", discrete.data$Species)

cleaned <- geiger::treedata(phy = tree, data = discrete.data, sort = TRUE)
print(ape::Ntip(cleaned$phy))
print(nrow(cleaned$data))

visualizedata <- cleaned

#Use parsimony to look at ancestral states:
cleanded.discrete.phy <- phangorn::phyDat(discrete.data, type = "USER", levels = all)
anc.test <- phangorn::phyDat(new.disc.dat, type = "USER", levels = c(all))

anc.p <- phangorn::ancestral.pars(tree,cleanded.discrete.phy)
plotAnc(tree, anc.p, 1)

#Uncertainty? What does it mean?

# Likelihood estimates:
anc.ml <- ancestral.pml(pml(tree, cleanded.discrete.phy), type = "ml")
plot(anc.ml)

#How do parsimony and likelihood differ?
# Likelihood seeks to find tree topology that confers on the highest probability on the observed traits of tip species 
# Parsimony seeks to find the tree topology that requires the fewest changes in trait states to produce the tip species

#Why do they differ? 

#Biological questions about the data using the package corHMM:
#Estimate the transitions between states

# How could you examine if transitions rates are equal?

# Go back to the Lewis, 2001 article and the MKV model (also, think about the Harmon example)
 # Are traits variable?
  #Will this make sense for your data?
  # Test using the different methods and do results change?
  #How could you test the order of state evolution?




#### Things to delete:
ape::read.tree(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/squamate.phy"))
read.csv(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/brandley_table.csv"))


discrete.data.old <- read.csv("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/squamate.phy", stringsAsFactors = FALSE) #I don't think this will work for the HW, so don't use

cleandata <- treedata(tree, limbless)

model <- fitDiscrete(cleandata$phy, cleandata$data)
cleaned.discrete.phydata <- phangorn::phyDat(cleandata)

#likelihood
td <- treedata(tree, discrete.data)

#Look at the data:
viewdata <- phytools::contMap(tree, discrete.data[1:258,])

#Check concordance between data file and tree (tree.not.data= taxa in tree but not data, Data.not.tree=taxa in data but not tree):
Check <- name.check(tree, discrete.data, data.names = NULL) 

## matching data to tree
newphytree <- drop.tip(discrete.data$phy)
name.check(newphytree, discrete.data$dat)

