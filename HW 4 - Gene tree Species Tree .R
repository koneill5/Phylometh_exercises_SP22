
devtools::install_github("bomeara/phybase") # this installs a modified package B. OMeara wrote for creating trees.

library(rotl)
library(ape)
phy <- get_study_tree("ot_485", "tree1")
plot(phy, cex=0.3)

#Simplify the tree
install.packages("geiger")
install.packages("phangorn")
library(geiger)

phy <- drop.random(phy, Ntip(phy) - 10)
plot(phy)
axisPhylo()

#Simulate gene trees
library(phybase)
gene.tree <- phybase::sim.coaltree.phylo(phy, pop.size = 1e-12)
plot(gene.tree)

#Comparing trees
library(phytools)
plot(cophylo(phy, gene.tree, cbind(sort(phy$tip.label), sort(gene.tree$tip.label))))

#Rotating nodes 
species.tree <- rcoal(7)
species.tree$edge.length <- species.tree$edge.length / (10*max(branching.times(species.tree)))
gene.tree <- phybase::sim.coaltree.phylo(species.tree)
plot(cophylo(species.tree, gene.tree, cbind(sort(species.tree$tip.label), sort(gene.tree$tip.label))))

#Lengthen the tips of the species tree 
tip.rows <- which(species.tree$edge[,2]<Ntip(species.tree))
species.tree2 <- species.tree
species.tree2$edge.length[tip.rows] <- 100 + species.tree2$edge.length[tip.rows]
gene.tree2 <- phybase::sim.coaltree.phylo(species.tree2)
plot(cophylo(species.tree2, gene.tree2, cbind(sort(species.tree2$tip.label), sort(gene.tree2$tip.label))))

# Plot cladogram and change branch lengths 
species.tree2.clado <- compute.brlen(species.tree2)
gene.tree2.clado <- compute.brlen(gene.tree2)
plot(cophylo(species.tree2.clado, gene.tree2.clado, cbind(sort(species.tree2.clado$tip.label), sort(gene.tree2.clado$tip.label))))
