library(geiger)
library(rotl)
library(ape)
library(phytools)
library(phangorn)
library(tidyverse)
library(OUwie)
eval=TRUE


#data.contin <- ape::read.tree(text = RCurl::getURL("https://github.com/lukejharmon/ilhabela/blob/master/workingFiles/continuousModels/anolisDataAppended.csv"))

save(data, file = 'anole.data.rda')
load(file = 'anole.data.rda')
anole.data <- read_file("anole.data.rda")


tree<-read.tree(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/ilhabela/master/workingFiles/continuousModels/anolis.phy"))
plotTree(tree)

continuous.data <- read.csv("Anole.LJHarmon.csv", stringsAsFactors = FALSE)
#continuous.data<-read.csv(text=RCurl::getURL("https://github.com/lukejharmon/ilhabela/blob/master/workingFiles/continuousModels/anolisDataAppended.csv"),stringsAsFactors=FALSE)


