#install.packages(c("ape", "TreeSim", "geiger", "diversitree", "devtools"))
library(ape)
library(TreeSim)
library(geiger)
library(diversitree)
#devtools::install_github("thej022214/hisse")
library(hisse)

my.tree <- TreeSim::sim.bd.taxa(n=300, numbsim=1, lambda=0.1, mu=0)[[1]]
plot(my.tree)
ape::ltt.plot(my.tree)
ape::ltt.plot(my.tree, log="y")

yule.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.1, mu=0, complete=FALSE)
#stop("How to do a multiple ltt plot?")

bd.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=.9, complete=FALSE)
ape::mltt.plot(bd.trees, log="y", legend=FALSE)

depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
max.depth <- sum(abs(depth.range)) #ape rescales depths
plot(x=c(0, -1*max.depth), y=c(1, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)

#To do:
my.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=0.005, complete=FALSE)
ape::mltt.plot(my.trees, log="y", legend=FALSE)

#What happens if speciation rate is much higher than extinction? 
# lambda is the constant rate of speciation, mu is the constant rate of extinction 
#

#How does the simulation change with different values, but keeping their difference constant? 

#If their sum is constant?



##Tree plus trait model ***make sure to re-run libraries before running the tree.musse
speciation.rates <- c(0.1, 0.1, 0.1, 0.2) #0A, 1A, 0B, 1B
extinction.rates <- rep(0.03, 4)
transition.rates <- c(0.01,0.01,0, 0.01, 0, 0.01, 0.01,0,0.01, 0,0.01,0.01)
pars <- c(speciation.rates, extinction.rates, transition.rates)
phy <- tree.musse(pars, max.taxa=50, x0=1, include.extinct=FALSE)
sim.dat.true <- data.frame(names(phy$tip.state), phy$tip.state)
sim.dat <- sim.dat.true
# Now to hide the "hidden" state
sim.dat[sim.dat[,2]==3,2] = 1
sim.dat[sim.dat[,2]==4,2] = 2
# and convert states 1,2 to 0,1
sim.dat[,2] = sim.dat[,2] - 1

#Look:
plot(phy)
knitr::kable(cbind(sim.dat, true.char=sim.dat.true$phy.tip.state))

#index vectors: 
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)

#Separate turnover rates for states 0A and 1A
turnover.anc = c(1,2,0,0)

#Full hisse model 
turnover.anc = c(1,2,3,4)

#Model that corresponds to no extinction / Yule (pure birth) equivalent
eps.anc = c(0,0,0,0)

#Setting up the transition rate matrix
trans.rates = TransMatMaker.old(hidden.states = TRUE)
trans.rates

trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual

trans.rates.nodual.equal16 = ParEqual(trans.rates.nodual, c(1,6))
trans.rates.nodual.equal16

trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

trans.rates.nodual.allequal = trans.rates.nodual
trans.rates.nodual.allequal[!is.na(trans.rates.nodual.allequal) & !trans.rates.nodual.allequal == 0] = 1
trans.rates.nodual.allequal

#BiSSE model in a HiSSE:
trans.rates.bisse = TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse

pp = hisse.old(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc,
               eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)

#A common mistake 
#Test where a hidden state is associated with a single observed state
turnover.anc = c(1,2,0,3)
eps.anc = c(1,2,0,3)

trans.rates <- TransMatMaker.old(hidden.states=TRUE)
trans.rates.nodual.no0B <- ParDrop(trans.rates, c(2,3,5,7,8,9,10,12))
trans.rates.nodual.no0B

#Changing the output 
pp = hisse.old(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc,
               eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal, output.type="net.div")

#Setting up the 2-state character-independent (CID-2) model 
#"null-two" model:
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)

#Assume the usual 8 transitions
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#Assuming all rates are equal:
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

#third option:
# Now we want three specific rates:
trans.rates.nodual.threerates <- trans.rates.nodual
# Set all transitions from 0->1 to be governed by a single rate:
to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1
# Now set all transitions from 1->0 to be governed by a single rate:
to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2
# Finally, set all transitions between the hidden state to be a single rate (essentially giving
# you an estimate of the rate by which shifts in diversification occur:
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3
trans.rates.nodual.threerates

#HiSSE
pp = hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc,
           eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)

#Plotting HiSSE
load("testrecon1.rda")
class(pp.recon)
pp.recon

#Plot net diversification
plot.hisse.states(pp.recon, rate.param="net.div", show.tip.label=FALSE)

#plot with a vector with minimum and maximum rate
plot.hisse.states(pp.recon, rate.param="net.div", show.tip.label=FALSE, rate.range=c(0,0.072))
#The blue areas become black 

pp.recon$aic
#Use the AIC in the model
pp.recon = MarginRecon(phy, sim.dat, f=c(1,1), hidden.states=TRUE, pars=pp$solution,
                       aic=pp$aic, n.cores=2)

hisse.results.list = list()
load("testrecon1.rda")
hisse.results.list[[1]] = pp.recon
load("testrecon2.rda")
hisse.results.list[[2]] = pp.recon
load("testrecon3.rda")
hisse.results.list[[3]] = pp.recon
# Now supply the list the plotting function
plot.hisse.states(hisse.results.list, rate.param="net.div", show.tip.label=FALSE, rate.range=c(0,0.072))


# First, suck in all the files with .Rsave line ending in your working directory:
files = system("ls -1 | grep .rda", intern=TRUE)
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(files[i])
  hisse.results.list[[i]] = pp.recon
  rm(pp.recon)
}

















