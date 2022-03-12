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

new.disc.dat <- cbind(limbed.dat, limbless)
rownames(new.disc.dat) <- gsub(" ", "_", discrete.data$Species)
 
cleaned <- geiger::treedata(phy = tree,data = new.disc.dat,sort = TRUE)
print(ape::Ntip(cleaned$phy))
print(nrow(cleaned$data))

visualizedata <- cleaned

#Use parsimony to look at ancestral states:
cleaned.discrete.phy <- phangorn::phyDat(cleaned$data, type = "USER", levels = c(0,1))
anc.p <- phangorn::ancestral.pars(cleaned$phy, cleaned.discrete.phy)
plotAnc(tree,anc.p, 1)


# Likelihood estimates:
anc.ml <- ancestral.pml(pml(tree,cleaned.discrete.phy), type = "ml")
plotAnc(tree,anc.ml, 1)



#How do parsimony and likelihood differ?
# Likelihood seeks to find tree topology that confers on the highest probability on the observed traits of tip species 
# Parsimony seeks to find the tree topology that requires the fewest changes in trait states to produce the tip species

#Why do they differ? 
# Parsimony and likelihood differ due to the contributions each trait/character mkes to the length of the overall tree. In parsimony the number of steps a character has in a tree is indepdendent from all else. In a likelihood, the value of a character's contribution/steps depend on the information contributed from all other characters, thus is not indepdent.  

#What does uncertainty mean?
# I think uncertainty uncertainty is asking whether certain traits are capable of rapid evolutionary change and if those traits change in a biased manner over evolutionary time. 
# An example would be to examine the number of times a trait was gained or lost over the evolutionary time of taxa of interest. 


#Biological questions about the data using the package corHMM:
library(corHMM)

#1. How can you estimate the transitions between states?

# 2. How could you examine if transitions rates are equal?
# There are two main methods used to examine transition rates, a likelihood-based ancestral state reconstruction or a hidden Markov rates model. 
#Beaulieu, J. M., O’Meara, B. C., & Donoghue, M. J. (2013). Identifying Hidden Rate Changes in the Evolution of a Binary Morphological Character: The Evolution of Plant Habit in Campanulid Angiosperms. Systematic Biology, 62(5), 725–737. https://doi.org/10.1093/sysbio/syt034


# 3. Go back to the Lewis, 2001 article and the MKV model (also, think about the Harmon example)
 # 3a.Are traits variable?
  #3b. Will this make sense for your data?
  # 3c.Test using the different methods and do results change?
MKV.test. <- corHMM(phy = tree, data = new.disc.dat, rate.cat = 1)
  
#4. How could you test the order of state evolution?
# I think you can use a likelihood ratio test (going back to Pagel, 1994) to answer questions regarding the evolutionary order of change among two traits (i.e., does the second trait evolve from state 0 to state 1 in a similar manner as the first trait).
# You can also test for order of state evolution from a Bayesian framework. In this instance, the probability of the tree and the specific observed traits are conditional on the prior probability and the parameters set. 



#### Things to delete/didn't work:
#ape::read.tree(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/squamate.phy"))
#read.csv(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/brandley_table.csv"))

#oldnew.disc.dat <- cbind(discrete.data$Species,limbed.dat,limbless)
#old: newdat2 <- gsub(" ","_",new.disc.dat[,1])
#need to have the species names as the rows ***There needs to be a space between the quotemarks 
#rownames(discrete.data) <- gsub(" ","_", discrete.data$Species)
#cleaned <- geiger::treedata(phy = tree, data = discrete.data, sort = TRUE) : this didn't work

#discrete.data.old <- read.csv("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/squamate.phy", stringsAsFactors = FALSE) #I don't think this will work for the HW, so don't use
#old: cleanded.discrete.phy <- phangorn::phyDat(discrete.data, type = "USER", levels = all)

#anc.p <- phangorn::ancestral.pars(tree,cleanded.discrete.phy)
#plotAnc(tree, anc.p, 1)
#anc.test <- phangorn::phyDat(new.disc.dat, type = "USER", levels = c(all))
#cleandata <- treedata(tree, limbless)

#model <- fitDiscrete(cleandata$phy, cleandata$data)
#cleaned.discrete.phydata <- phangorn::phyDat(cleandata)
#anc.ml <- ancestral.pml(pml(tree, cleanded.discrete.phy), type = "ml")
#plot(anc.ml)
#likelihood
#td <- treedata(tree, discrete.data)

#Look at the data:
#viewdata <- phytools::contMap(tree, discrete.data[1:258,])

#Check concordance between data file and tree (tree.not.data= taxa in tree but not data, Data.not.tree=taxa in data but not tree):
#Check <- name.check(tree, discrete.data, data.names = NULL) 

## matching data to tree
#newphytree <- drop.tip(discrete.data$phy)
#name.check(newphytree, discrete.data$dat)

