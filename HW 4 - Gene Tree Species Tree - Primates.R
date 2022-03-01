
devtools::install_github("bomeara/phybase")

library(rotl)
library(ape)
phylo <- get_study_tree("ot_1553", "tree1") #These are primate NWM primate taxa 
plot(phylo, cex=0.3)

#Simplify the tree 
install.packages("geiger")
install.packages("phangorn")
library(geiger)

phylo.dropped <- drop.random(phylo, Ntip(phylo) - 10)
plot(phylo.dropped)
axisPhylo()

#Simulate gene tree
library(phybase)
gene.tree.dropped <- phybase::sim.coaltree.phylo(phylo.dropped, pop.size = 1e-12)
plot(gene.tree.dropped)

#Comparing trees
library(phytools)
plot(cophylo(phylo.dropped, gene.tree.dropped, cbind(sort(phylo.dropped$tip.label), sort(gene.tree.dropped$tip.label))))

#Rotating nodes
species.tree.dropped <- rcoal(7)
species.tree.dropped$edge.length <- species.tree.dropped$edge.length / (10*max(branching.times(species.tree.dropped)))
gene.tree.dropped <- phybase::sim.coaltree.phylo(species.tree.dropped)
plot(cophylo(species.tree.dropped, gene.tree.dropped, cbind(sort(species.tree.dropped$tip.label), sort(gene.tree.dropped$tip.label))))

#Lengthen the tips of the species tree
tip.rows.dropped <- which(species.tree.dropped$edge[,2]<Ntip(species.tree.dropped))
species.tree.dropped2 <- species.tree.dropped
species.tree.dropped2$edge.length[tip.rows.dropped] <- 100 + species.tree.dropped2$edge.length[tip.rows.dropped]
gene.tree.dropped2 <- phybase::sim.coaltree.phylo(species.tree.dropped2)
plot(cophylo(species.tree.dropped2, gene.tree.dropped2, cbind(sort(species.tree.dropped2$tip.label), sort(gene.tree.dropped2$tip.label))))

#Plot cladogram and change branch lengths 
species.tree.dropped2.clado <- compute.brlen(species.tree.dropped2)
gene.tree.dropped2.clado <- compute.brlen(gene.tree.dropped2)
plot(cophylo(species.tree.dropped2.clado, gene.tree.dropped2.clado, cbind(sort(species.tree.dropped2.clado$tip.label), sort(gene.tree.dropped2.clado$tip.label))))


