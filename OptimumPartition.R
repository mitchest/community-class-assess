## Author: Mitchell Lyons (mitchell.lyons@gmail.com)
## No licence, though please attribute me and publish any derivitives open source =)

# model fitting/clustering
library(mvabund)
library(tweedie)
library(statmod)
library(vegan)
library(ordinal)
library(clustsig)
library(cluster)
library(labdsv)
#library(betareg) # using gamlss instead
library(gamlss)
library(gamlss.dist)
# titbits
library(Cairo)
library(plyr)
library(reshape)
library(ggplot2)

setwd("")
source("OptimumPartition_functions.R")
source("vegsim.R")

load("AbundanceData.RData")


# explore data a little ---------------------------------------------------
cover.unique = apply(cover, 1, unique)
sums = numeric(522)
for (i in 1:552) {
  sums[i] = sum(!cover.unique[[i]] %in% c(0:5, seq(10,100,5)))
}
sum(sums)

par(mfcol=c(2,2))
plot(as.matrix(scores),as.matrix(cover), xlab="scores", ylab="cover", main="cover ~ scores")
plot(as.matrix(scores),as.matrix(counts), xlab="scores", ylab="counts", main="counts ~ scores")
plot(as.matrix(cover),as.matrix(counts), xlab="cover", ylab="counts", main="counts ~ cover")



# define clusters ---------------------------------------------------------

# transform and felxible beta clustering
cover.log = log(cover+1)
cover.log.dist = dsvdis(cover.log, 'bray/curtis')
cover.flex = agnes(cover.log.dist, method='flexible', par.method=0.625)

# check for obvious outliers
for (i in 2:50){
  print(table(cutree(cover.flex,k=i)))
}



# Fisher's exact test (OptimClass) ----------------------------------------

# get p-values for each clustering solution
pvals.allocs.list = list()
for (i in 2:50) { # choose max k clusters here
  print(i)
  alloc = cutree(cover.flex,k=i)
  pvals.allocs.list[[i]] = fishers.pvals(pa, alloc)
}
save(pvals.allocs.list, file="pvals.allocs.list.RData")
load("pvals.allocs.list.RData")


CairoPDF(file="figures/OptimClass_results.pdf", width=10, height=10)
par(mfcol=c(3,2))

OptimClass1(pvals.allocs.list, 0.001)
OptimClass1(pvals.allocs.list, 0.000001)
OptimClass1(pvals.allocs.list, 0.000000001)

OptimClass2(pvals.allocs.list, 0.001, 3)
OptimClass2(pvals.allocs.list, 0.000001, 3)
OptimClass2(pvals.allocs.list, 0.000000001, 3)

dev.off()



# sum-of-aic for difference response types --------------------------------

# loop through different cutoff levels and model types to test sum-of-AIC on groupings
# set cluster numbers and cluster source
groups = 2:50
clust = cover.flex

# counts
data = mvabund(counts)
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manyglm(formula=data~as.factor(alloc), family="negative.binomial")$AICsum
}
counts.AICsum = data.frame(AICsum,groups)
save(counts.AICsum,file="counts.AICsum.RData")

# P/A
data = mvabund(pa)
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manyglm(formula=data~as.factor(alloc), family="binomial")$AICsum
}
pa.AICsum = data.frame(AICsum,groups)
save(pa.AICsum,file="pa.AICsum.RData")

# scores
data = data.frame(lapply(scores, as.factor))
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manyclm.AICsum(data, as.factor(alloc))
}
scores.AICsum = data.frame(AICsum,groups)
save(scores.AICsum,file="scores.AICsum.RData")

# cover
# use ordinal regression
data = data.frame(lapply(cover, as.factor))
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manyclm.AICsum(data, as.factor(alloc))
}
cover.AICsum.ord = data.frame(AICsum,groups)
save(cover.AICsum.ord,file="cover.AICsum.ord.RData")
#data = ((cover/100)*(nrow(cover)-1)+0.5)/nrow(cover) # this is a transformation for betareg
data = cover/100
# use zero/one inflated beta
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manybeta.AICsum(data, as.factor(alloc))
}
cover.AICsum.infbe = data.frame(AICsum,groups)
save(cover.AICsum.infbe,file="cover.AICsum.infbe.RData")
# use compund Poisson-gamma (i.e. tweedie with power= 1.5)
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  Xdata = data.frame(alloc=cutree(clust, groups[i]), dummy=integer(nrow(data)))
  AICsum[i] = sum(AIC(manyany("glm", data, data~alloc, data=Xdata,
                              family=tweedie(var.power=1.2, link.power=0), var.power=1.2, composition=FALSE)))
}
cover.AICsum.tw = data.frame(AICsum,groups)
save(cover.AICsum.tw,file="cover.AICsum.tw.RData")


# load results if starting fresh
load("counts.AICsum.RData")
load("pa.AICsum.RData")
load("scores.AICsum.RData")
load("cover.AICsum.infbe.RData")
load("cover.AICsum.tw.RData")
load("cover.AICsum.ord.RData")

# plot AIC
#CairoWin()
counts.AICsum$AICsum = minus.min(counts.AICsum)
pa.AICsum$AICsum = minus.min(pa.AICsum)
scores.AICsum$AICsum = minus.min(scores.AICsum)
cover.AICsum.infbe$AICsum = minus.min(cover.AICsum.infbe)
cover.AICsum.tw$AICsum = minus.min(cover.AICsum.tw)
cover.AICsum.ord$AICsum = minus.min(cover.AICsum.ord)

CairoPDF(file="figures/sumAIC_real.pdf", width=12, height=8)
par(mfrow=c(2,2), cex.lab=1.5, cex.main=1.7, oma=c(4,4,0,0)+0.1, mar=c(2,0,4,3)+0.1)

plot(counts.AICsum$AICsum~counts.AICsum$groups, main="counts - negative binonial", 
     ylab="", xlab="", type='l')
points(x=WhichMin(counts.AICsum)[2],y=WhichMin(counts.AICsum)[1], pch=10, cex=2)

plot(pa.AICsum$AICsum~pa.AICsum$groups, main="pres/abs - binomal", 
     ylab="", xlab="", type='l')
points(x=WhichMin(pa.AICsum)[2],y=WhichMin(pa.AICsum)[1], pch=10, cex=2)

plot(scores.AICsum$AICsum~scores.AICsum$groups, main="ordinal scores - prop. odds", 
     ylab="", xlab="", type='l')
points(x=WhichMin(scores.AICsum)[2],y=WhichMin(scores.AICsum)[1], pch=10, cex=2)

plot(cover.AICsum.ord$AICsum~cover.AICsum.ord$groups, main="% cover - prop. odds", 
     ylab="", xlab="", type='l')
points(x=WhichMin(cover.AICsum.ord)[2],y=WhichMin(cover.AICsum.ord)[1], pch=10, cex=2)

# plot(cover.AICsum.infbe$AICsum~cover.AICsum.infbe$groups, main="% cover - inf. beta", 
#      ylab="", xlab="", type='l')
# points(x=WhichMin(cover.AICsum.infbe)[2],y=WhichMin(cover.AICsum.infbe)[1], pch=10, cex=2)
# 
# plot(cover.AICsum.tw$AICsum~cover.AICsum.tw$groups, main="% cover - Poisson-gamma", 
#      ylab="", xlab="", type='l')
# points(x=WhichMin(cover.AICsum.tw)[2],y=WhichMin(cover.AICsum.tw)[1], pch=10, cex=2)

title(xlab = "number of clusters",
      ylab = "sum-of-AIC",
      outer = TRUE, line = 3)

dev.off()



# indicator species -------------------------------------------------------

alloc = cutree(cover.flex, 11)

# P/A
data = mvabund(pa)
pa.manyglm = manyglm(formula=data~as.factor(alloc), family="binomial")
pa.manyglm.null = manyglm(formula=data~1, family="binomial")
pa.manyglm.species = data.frame(sort(pa.manyglm.null$aic-pa.manyglm$aic, decreasing=T))
# clean
pa.manyglm.species = data.frame(species=row.names(pa.manyglm.species), 
                                dAIC=pa.manyglm.species$sort.pa.manyglm.null.aic)
save(pa.manyglm.species, file="pa.manyglm.species.RData")
load("pa.manyglm.species.RData")

# counts
data = mvabund(counts)
counts.manyglm = manyglm(formula=data~as.factor(alloc), family="negative.binomial")
counts.manyglm.null = manyglm(formula=data~1, family="negative.binomial")
counts.manyglm.species = data.frame(sort(counts.manyglm.null$aic-counts.manyglm$aic, decreasing=T))
# clean
counts.manyglm.species = data.frame(species=row.names(counts.manyglm.species), 
                                    dAIC=counts.manyglm.species$sort.counts.manyglm.null.aic)
save(counts.manyglm.species, file="counts.manyglm.species.RData")
load("counts.manyglm.species.RData")

# scores
data = data.frame(lapply(scores, as.factor))
scores.manyclm = manyclm.AIC(data,as.factor(alloc))
scores.manyclm.null = manyclm.AIC(data,1)
scores.manyclm.species = data.frame(sort(scores.manyclm.null-scores.manyclm, decreasing=T))
# clean
scores.manyclm.species = data.frame(species=row.names(scores.manyclm.species),
                                    dAIC=scores.manyclm.species$sort.scores.manyclm.null)
save(scores.manyclm.species, file="scores.manyclm.species.RData")
load("scores.manyclm.species.RData")

# cover
data = data.frame(lapply(cover, as.factor))
cover.manyclm = manyclm.AIC(data,as.factor(alloc))
cover.manyclm.null = manyclm.AIC(data,1)
cover.manyclm.species = data.frame(sort(cover.manyclm.null-cover.manyclm, decreasing=T))
# clean
cover.manyclm.species = data.frame(species=row.names(cover.manyclm.species),
                                   dAIC=cover.manyclm.species$sort.cover.manyclm.null)
save(cover.manyclm.species, file="cover.manyclm.species.RData")
load("cover.manyclm.species.RData")


# check top contributors
top20 = cbind(pa.manyglm.species[1:20,], scores.manyclm.species[1:20,])
sum(pa.manyglm.species[1:100,1] %in% scores.manyclm.species[1:100,1])
write.csv(top20, file="figures/top20_PAScores.csv", row.names=F)


# OptimClass species
# calculate shared species bewteen top 100 dAIC/fiathful species
faithful.species = OptimClass.species(pvals.allocs.list, 100)
faithful.11 = names(faithful.species[[11]])
sum(faithful.11[1:100] %in% pa.manyglm.species[1:100,1])
sum(faithful.11[1:100] %in% scores.manyclm.species[1:100,1])

# calculate shared species at different numbers faithful species and different Fisher's cuts
CairoPDF(file="figures/OptimClass_dAIC.pdf", width=6, height=6)

faithful.species = OptimClass.species(pvals.allocs.list, 712)
faithful.11 = names(faithful.species[[11]][faithful.species[[11]]<0.001])
n.shared = numeric(length(faithful.11)-20)
for (i in 20:length(faithful.11)){
  n.shared[i-19] = sum(unique(faithful.11[1:i]) %in% pa.manyglm.species[1:100,1])
}
plot(n.shared ~ c(20:length(faithful.11)), type="l",
     main=bquote(atop("OptimClass faithful species shared",
                      "with top 100 " ~ Delta ~ "AIC species")),
     ylab="no. shared species", xlab="no. most faithful species considered")
abline(a=0, b=1, lty=2, col="red")

faithful.11 = names(faithful.species[[11]][faithful.species[[11]]<0.000001])
n.shared = numeric(length(faithful.11)-20)
for (i in 20:length(faithful.11)){
  n.shared[i-19] = sum(unique(faithful.11[1:i]) %in% pa.manyglm.species[1:100,1])
}
abline(v=c(20:length(faithful.11))[length(n.shared)], lty=2)

faithful.11 = names(faithful.species[[11]][faithful.species[[11]]<0.000000001])
n.shared = numeric(length(faithful.11)-20)
for (i in 20:length(faithful.11)){
  n.shared[i-19] = sum(unique(faithful.11[1:i]) %in% pa.manyglm.species[1:100,1])
}
abline(v=c(20:length(faithful.11))[length(n.shared)], lty=2)

mtext(text=bquote(atop("cut =" ~10^{-~9})), side=1, line=-1.5, at=c(90))
mtext(text=bquote(atop("cut =" ~10^{-~6})), side=1, line=-1.5, at=c(185))
mtext(text=bquote(atop("cut =" ~10^{-~3})), side=1, line=-1.5, at=c(300))

dev.off()


# indicator values (fidelity and abundance)
indval.pa = indval(x=pa, clustering=alloc)
high.indval.pa = list()
for (i in 1:11) {
  vals = sort(indval.pa$indval[[i]], decreasing=T)[1:20]
  high.indval.pa[[i]] = names(pa)[which(indval.pa$indval[[i]] %in% vals[vals>0])]
  rm(vals)
}
high.indval.pa = unique(unlist(high.indval.pa))
sum(high.indval.pa %in% pa.manyglm.species[1:100,1])
sum(high.indval.pa %in% faithful.11[1:100])

n.shared.indval = numeric(30)
for (j in 1:30) {
  high.indval.pa = list()
  for (i in 1:11) {
    vals = sort(indval.pa$indval[[i]], decreasing=T)[1:j]
    high.indval.pa[[i]] = names(pa)[which(indval.pa$indval[[i]] %in% vals[vals>0])]
    rm(vals)
  }
  high.indval.pa = unique(unlist(high.indval.pa))
  n.shared.indval[i] = sum(high.indval.pa %in% pa.manyglm.species[1:100,1])
}


indval.scores = indval(x=scores, clustering=alloc)
high.indval.scores = list()
for (i in 1:11) {
  vals = sort(indval.scores$indval[[i]], decreasing=T)[1:20]
  high.indval.scores[[i]] = names(scores)[which(indval.scores$indval[[i]] %in% vals[vals>0])]
  rm(vals)
}
high.indval.scores = unique(unlist(high.indval.scores))
sum(high.indval.scores %in% pa.manyglm.species[1:100,1])
sum(high.indval.scores %in% faithful.11[1:100])



# sum-of-aic backward selection -------------------------------------------

start.k = 50 # choose initial k for merging
n.iter = 48 # choose number of merging iterations

alloc.test = cutree(cover.flex, start.k) 
for (i in 1:n.iter) {
  print(paste0("Starting merging iteration ",i,":"))
  pairs = GeneratePairwise(alloc.test)
  #PairwiseCombinedSpeciesCount(pa, alloc.test, pairs)
  alloc.test.pairwise = PairwiseManyglm(pa, alloc.test, pairs)
  #alloc.test.pairwise.sorted = alloc.test.pairwise[with(alloc.test.pairwise, order(deltaSumAIC)),]
  # calculate summary and merge clusters using this information
  PairwiseSummary(alloc.test, alloc.test.pairwise)
  merge.pair = c(as.numeric(alloc.test.pairwise$target[alloc.test.pairwise$deltaSumAIC==min(alloc.test.pairwise$deltaSumAIC)]),
                 as.numeric(alloc.test.pairwise$test[alloc.test.pairwise$deltaSumAIC==min(alloc.test.pairwise$deltaSumAIC)]))
  alloc.test[alloc.test==merge.pair[1] | alloc.test==merge.pair[2]] = 1000 + i
  print("#############")
  print(paste0("merged cluster: ",merge.pair[1]," and ",merge.pair[2]))
  print("#############")
  assign(x=paste0("alloc.test.",start.k-i), value=alloc.test)
}

# fit and get AIC

alloc.list = list()
for (i in 1:length(2:49)) {
  alloc.list[[i]] = get(paste0("alloc.test.",c(2:49)[i]))
}
save(alloc.list, file="alloc.list.RData")

data = mvabund(pa)
alloc.sumAIC = numeric(length(alloc.list))
for (i in 1:length(alloc.list)) {
  alloc.sumAIC[i] = manyglm(formula=data~as.factor(alloc.list[[i]]), family="binomial")$AICsum
}

alloc.pairwise.AIC = data.frame(AICsum=alloc.sumAIC, groups=c(2:49))
load("pa.AICsum.RData")
pa.AICsum.pairwise = rbind(pa.AICsum, alloc.pairwise.AIC)
save(pa.AICsum.pairwise, file="pa.AICsum.pairwise.RData")
load("pa.AICsum.pairwise.RData")


# plot pairwise modelling merges
AIC.hclust = pa.AICsum.pairwise[1:49,]
AIC.pairwise = pa.AICsum.pairwise[50:97,]

#CairoWin()
CairoPDF(file="figures/ClassRefine.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(AIC.hclust$AICsum~AIC.hclust$groups, type='l', 
     ylim=c(min(AIC.pairwise$AICsum),max(AIC.hclust$AICsum)),
     ylab="sum-of-AIC", xlab="number of clusters", main="")
points(x=WhichMin(AIC.hclust)[2],y=WhichMin(AIC.hclust)[1], pch=10, cex=2)
lines(y=AIC.pairwise$AICsum, x=AIC.pairwise$groups, lty=3)
points(x=WhichMin(AIC.pairwise)[2],y=WhichMin(AIC.pairwise)[1], pch=10, cex=2)
legend("topleft", legend=expression("dendrogram cutting", Delta~" sum-of-AIC pruning"), lty=c(1,3))
dev.off()



# sum-of-AIC for simulated abundance --------------------------------------

# counts

run.sims.AIC = function(uniformWidth, randomSample, nboot=100, groups=2:20) {
  # loop through different cutoff levels and simulated data sets to test sum-of-AIC on groupings
  # set cluster numbers and cluster source
  groups = groups
  # start nboot runs
  counts.sim.AICsum = list()
  nboot = nboot
  sim.clust.AICsum = numeric(nboot)
  for (sim in 1:nboot) {
    # source vegSim.R/run VegSim() to generate some synthetic data
    set.seed(sim+42*42)
    sampled = VegSim(C=9, S=10, Nmult=10, width=0.2, uniformWidth=uniformWidth, randomSample=randomSample, plotSim=F, plotOpt=F)
    counts.sim = sampled[,-1]
    #clusters.sim = sampled[,1]
    # cluster the simulated data
    clust = agnes(dsvdis(counts.sim, 'bray/curtis'), method='flexible', par.method=0.625)
    #clust = hclust(vegdist(counts.sim, "bray"), "average")
    # loop through grouping levels 
    data = mvabund(counts.sim)
    AICsum = numeric(length(groups))
    for (i in 1:length(groups)) {
      alloc = cutree(clust, groups[i])
      AICsum[i] = manyglm(formula=data~as.factor(alloc), family="negative.binomial")$AICsum
    }
    counts.sim.AICsum[[sim]] = data.frame(AICsum,groups)
    sim.clust.AICsum[sim] = manyglm(formula=data~as.factor(sampled[,1]), family="negative.binomial")$AICsum
  }
  #save for narrow niche
  return(list(counts.sim.AICsum,sim.clust.AICsum))
}

sim.results = run.sims.AIC(uniformWidth=T, randomSample=F)
AICsum.narrow.high = sim.results[[1]]
true.narrow.high = sim.results[[2]]

sim.results = run.sims.AIC(uniformWidth=T, randomSample=T)
AICsum.narrow.low = sim.results[[1]]
true.narrow.low = sim.results[[2]]

sim.results = run.sims.AIC(uniformWidth=F, randomSample=F)
AICsum.variable.high = sim.results[[1]]
true.variable.high = sim.results[[2]]

sim.results = run.sims.AIC(uniformWidth=F, randomSample=T)
AICsum.variable.low = sim.results[[1]]
true.variable.low = sim.results[[2]]

# plots
# sorry, this section is lazy - should functionalise the plotting...
#CairoWin()
CairoPDF(file="figures/sumAIC_simulated.pdf", width=8, height=8)
par(mfrow=c(2,2), mar=c(3,3,4,0.15), cex.lab=1.2, cex.main=1.2, mgp=c(2.1,1,0))

maxAIC = numeric(length(AICsum.narrow.high))
minAIC = numeric(length(AICsum.narrow.high))
for (i in 1:length(AICsum.narrow.high)) { maxAIC[i]=(max(c(AICsum.narrow.high[[i]]$AICsum,true.narrow.high)))}
for (i in 1:length(AICsum.narrow.high)) { minAIC[i]=(min(c(AICsum.narrow.high[[i]]$AICsum,true.narrow.high)))}
plot(AICsum.narrow.high[[1]]$AICsum~AICsum.narrow.high[[1]]$groups, type="n",
     ylim=c(min(minAIC),max(maxAIC)), ylab="sum-of-AIC", #yaxt='n', 
     xlab="", main="Narrow niches, high membership prob.")
mtext(text="(a)", side=1, line=-1.5, at=3)
for (i in 3:length(AICsum.narrow.high)) {
  lines(AICsum.narrow.high[[i]]$AICsum~AICsum.narrow.high[[i]]$groups, col="grey")
  points(x=WhichMin(AICsum.narrow.high[[i]])[2],y=WhichMin(AICsum.narrow.high[[i]])[1], col="grey", pch=10)
}
for (i in 1:2) {
  col=i+1
  lines(AICsum.narrow.high[[i]]$AICsum~AICsum.narrow.high[[i]]$groups, col=col)
  points(x=WhichMin(AICsum.narrow.high[[i]])[2],y=WhichMin(AICsum.narrow.high[[i]])[1], col=col, pch=10, cex=1.5)
  points(x=9, y=true.narrow.high[i], pch=9, col=col, cex=1.5)
}

maxAIC = numeric(length(AICsum.narrow.low))
minAIC = numeric(length(AICsum.narrow.low))
for (i in 1:length(AICsum.narrow.low)) { maxAIC[i]=(max(c(AICsum.narrow.low[[i]]$AICsum,true.narrow.low)))}
for (i in 1:length(AICsum.narrow.low)) { minAIC[i]=(min(c(AICsum.narrow.low[[i]]$AICsum,true.narrow.low)))}
plot(AICsum.narrow.low[[1]]$AICsum~AICsum.narrow.low[[1]]$groups, type="n",
     ylim=c(min(minAIC),max(maxAIC)), ylab="", #yaxt='n', 
     xlab="", main="Narrow niches, low membership prob.")
mtext(text="(b)", side=1, line=-1.5, at=3)
for (i in 3:length(AICsum.narrow.low)) {
  lines(AICsum.narrow.low[[i]]$AICsum~AICsum.narrow.low[[i]]$groups, col="grey")
  points(x=WhichMin(AICsum.narrow.low[[i]])[2],y=WhichMin(AICsum.narrow.low[[i]])[1], col="grey", pch=10)
}
for (i in 1:2) {
  col=i+1
  lines(AICsum.narrow.low[[i]]$AICsum~AICsum.narrow.low[[i]]$groups, col=col)
  points(x=WhichMin(AICsum.narrow.low[[i]])[2],y=WhichMin(AICsum.narrow.low[[i]])[1], col=col, pch=10, cex=1.5)
  points(x=9, y=true.narrow.low[i], pch=9, col=col, cex=1.5)
}

maxAIC = numeric(length(AICsum.variable.high))
minAIC = numeric(length(AICsum.variable.high))
for (i in 1:length(AICsum.variable.high)) { maxAIC[i]=(max(c(AICsum.variable.high[[i]]$AICsum,true.variable.high)))}
for (i in 1:length(AICsum.variable.high)) { minAIC[i]=(min(c(AICsum.variable.high[[i]]$AICsum,true.variable.high)))}
plot(AICsum.variable.high[[1]]$AICsum~AICsum.variable.high[[1]]$groups, type="n",
     ylim=c(min(minAIC),max(maxAIC)), ylab="sum-of-AIC", #yaxt='n', 
     xlab="number of clusters", main="Variable niches, high membership prob.")
mtext(text="(c)", side=1, line=-1.5, at=3)
for (i in 3:length(AICsum.variable.high)) {
  lines(AICsum.variable.high[[i]]$AICsum~AICsum.variable.high[[i]]$groups, col="grey")
  points(x=WhichMin(AICsum.variable.high[[i]])[2],y=WhichMin(AICsum.variable.high[[i]])[1], col="grey", pch=10)
}
for (i in 1:2) {
  col = i+1
  lines(AICsum.variable.high[[i]]$AICsum~AICsum.variable.high[[i]]$groups, col=col)
  points(x=WhichMin(AICsum.variable.high[[i]])[2],y=WhichMin(AICsum.variable.high[[i]])[1], col=col, pch=10, cex=1.5)
  points(x=9, y=true.variable.high[i], pch=9, col=col, cex=1.5)
}

maxAIC = numeric(length(AICsum.variable.low))
minAIC = numeric(length(AICsum.variable.low))
for (i in 1:length(AICsum.variable.low)) { maxAIC[i]=(max(c(AICsum.variable.low[[i]]$AICsum,true.variable.low)))}
for (i in 1:length(AICsum.variable.low)) { minAIC[i]=(min(c(AICsum.variable.low[[i]]$AICsum,true.variable.low)))}
plot(AICsum.variable.low[[1]]$AICsum~AICsum.variable.low[[1]]$groups, type="n",
     ylim=c(min(minAIC),max(maxAIC)), ylab="", #yaxt='n', 
     xlab="number of clusters", main="Variable niches, low membership prob.")
mtext(text="(d)", side=1, line=-1.5, at=3)
for (i in 3:length(AICsum.variable.low)) {
  lines(AICsum.variable.low[[i]]$AICsum~AICsum.variable.low[[i]]$groups, col="grey")
  points(x=WhichMin(AICsum.variable.low[[i]])[2],y=WhichMin(AICsum.variable.low[[i]])[1], col="grey", pch=10)
}
for (i in 1:2) {
  col = i+1
  lines(AICsum.variable.low[[i]]$AICsum~AICsum.variable.low[[i]]$groups, col=col)
  points(x=WhichMin(AICsum.variable.low[[i]])[2],y=WhichMin(AICsum.variable.low[[i]])[1], col=col, pch=10, cex=1.5)
  points(x=9, y=true.variable.low[i], pch=9, col=col, cex=1.5)
}

dev.off()



# residual plots for supp material ----------------------------------------

alloc = cutree(cover.flex, 11)

data = mvabund(counts)
# pois = manyglm(formula=data~as.factor(alloc), family="poisson")
nb = manyglm(formula=data~as.factor(alloc), family="negative.binomial")
plot(nb, which=1, caption="counts - neg bin regression", var.subset=sample(1:712,10))

data = mvabund(pa)
pa.fm = manyglm(formula=data~as.factor(alloc), family="binomial")
plot(pa.fm, which=1, caption="pres/abs - binomial regression", var.subset=sample(1:712,10))

clm.residuals = function(obj) {
  preds = predict(obj, type="cum.prob")
  residuals = runif(length(preds$cprob1))*(preds$cprob1-preds$cprob2) + preds$cprob2
  return(qnorm(residuals))
}
ord = clm(as.factor(scores$Abutoxyc)~as.factor(alloc)) # don't know how to define residual for ordinal regression
plot(clm.residuals(ord) ~ alloc, main="cover-abund - ordinal regression", xlab="cluster no.", ylab="Dunn-Smyth Residuals")
specs = sample(1:712,10)
for (i in 1:length(specs)){
  ord = clm(as.factor(scores[,i])~as.factor(alloc))
  points(clm.residuals(ord) ~ alloc, col=i)
}

Mean = colwise(function(x){mean(x)})(cover)
Variance = colwise(function(x){var(x)})(cover)
plot(as.numeric(log(Variance))~as.numeric(log(Mean)),
     main="mean-variance, no levels", xlab="log(Mean)", ylab="log(Variance)")
abline(a=0, b=1, lty=2)

