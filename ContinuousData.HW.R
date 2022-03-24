library(geiger)
library(rotl)
library(ape)
library(phytools)
library(phangorn)
library(tidyverse)
library(OUwie)
eval=TRUE


tree<-read.tree(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/ilhabela/master/workingFiles/continuousModels/anolis.phy"))
plotTree(tree)

continuous.data <- read.csv("Anole.LJHarmon.csv", stringsAsFactors = FALSE)
colnames(continuous.data) 

row.names(continuous.data) <- continuous.data$Species
continuous.data.SVL <- continuous.data$SVL
names(continuous.data.SVL) <- rownames(continuous.data)
#SVL is the anole body size measurements

clead.cont.data <- geiger::treedata(phy = tree, data = continuous.data, sort = TRUE)
visualize.cont.data <- geiger::treedata(phy = tree, data = continuous.data)

tree.map <- contMap(tree, continuous.data.SVL)


#What is the rate of evolution of your trait on the tree? 
BM1 <- geiger::fitContinuous(tree, continuous.data.SVL, model = "BM")
print(paste("The rate of evolution is",sigsq=0.136160, "in unites of",))

OUI <- fitContinuous(tree, continuous.data.SVL, model = "OU")
plot(tree, show.tip.label=FALSE)
ou.tree <- rescale(tree, model="OU", alpha=0, sigsq=0.136160)
plot(ou.tree)
# talk about difference 

#Compare trees:
aic.scores <-c(13.400807, 15.400807)
names(aic.scores) <- c("BM","OU")
aicw(aic.scores)

AIC.BM1 <- 13.400807
AIC.OU1 <- 15.400807
delta.AIC.BM1 <- 0
delta.AIC.OU1 <- 2



## Run OWie 
one.discrete.char <- continuous.data.SVL
reconstruction.info <- ace(one.discrete.char, tree, type="discrete", method = "ML", CI=TRUE)
best.states <- colnames(reconstruction.info$lik.anc)[apply(reconstruction.info$lik.anc, 1, which.max)]

plotTree(tree, node.numbers=TRUE)

#Label tree 
label.tree
select.reg <- character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree)

continuous.data[1:10,]

nodelabels(pch=21, bg=select.reg)

OUwie(tree,continuous.data.SVL,model = c("OUMV"))


run1 <- OUwie(tree,continuous.data,model = c("OUMV"), root.station = FALSE, algorithm = "invert", opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="50", "ftol_abs"=0.001))
nodeBased.OUMV <- OUwie(tree, one.discrete.char, model = "OUMV", simmap.tree = FALSE, diagn = FALSE)
nodeBased.OUMV2 <- OUwie(tree, clead.cont.data, model = "OUMV", simmap.tree = FALSE, diagn = FALSE)  

#fit an OU model based on a clade of interest 
OUwie(tree,continuous.data.SVL,model = c("OUMV"), clade = c("ahli", "websteri"), algorithm="three.point")


#OUwie models:
models <- c("BM1", "BMS", "OUI","OUM", "OUMV", "OUMA", "OUMVA" )
results <- lapply(models, phy=tree, data =continuous.data.SVL)








