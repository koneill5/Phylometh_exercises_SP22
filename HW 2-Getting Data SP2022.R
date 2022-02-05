library(rotl)
library(ape)

#Open Tree of life 
musteloidea.trees <- studies_find_trees(property = "ot:ottTaxonName", value = "Carnivora", detailed = FALSE)
musteloidea.studies.id <- unlist(musteloidea.trees$study_ids)


musteloidea.study2017.metadata <- rotl::get_study_meta(musteloidea.studies.id[45])
print(rotl::get_publication(musteloidea.study2017.metadata))


#Getting the tree I'm interested in
musteloidea.study2017.tree1 <- get_study(musteloidea.studies.id[1])[[1]] 

#Plot the tree 
ape::plot.phylo(musteloidea.study2017.tree1, type = "phylogram", cex = 0.2) #has branch lengths 
ape::plot.phylo(musteloidea.study2017.tree1, type = "cladogram", cex = 0.2)
ape::plot.phylo(musteloidea.study2017.tree1, type = "fan", cex = 0.2)
ape::plot.phylo(musteloidea.study2017.tree1, type = "radial", cex = 0.2)
ape::plot.phylo(musteloidea.study2017.tree1, type = "unrooted", cex = 0.2)





















