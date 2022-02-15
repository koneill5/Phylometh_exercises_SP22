
devtools::install_github("bomeara/phybase") # this installs a modified package B. OMeara wrote for creating trees.

library(rotl)
library(ape)
phy <- get_study_tree("ot_485", "tree1")
plot(phy, cex=0.3)

#Simplify the tree
install.packages("geiger")
install.packages("phangorn")
library(geiger)
