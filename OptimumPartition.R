## Author: Mitchell Lyons (mitchell.lyons@gmail.com)
## No licence, though please attribute me and publish any derivitives open source =)

# model fitting/clustering
library(mvabund)
library(tweedie)
library(statmod)
library(vegan)
library(ordinal)
library(clustsig)
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


# define clusters ---------------------------------------------------------

# create a dummy classification (bray curtis/hclust)
counts.bray = vegdist(counts, "bray")
counts.clust = hclust(counts.bray, "average")

# look at allocation in mean~var plot
counts.alloc = cutree(counts.clust, 8)
meanvar.plot(mvabund(counts)~as.factor(counts.alloc))



# sum-of-aic for difference response types --------------------------------

# loop through different cutoff levels and model types to test sum-of-AIC on groupings
# set cluster numbers and cluster source
groups = 2:20
clust = counts.clust

# counts
data = mvabund(counts)
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manyglm(formula=data~as.factor(alloc), family="negative.binomial")$AICsum
}
counts.AICsum = data.frame(AICsum,groups)

# P/A
data = mvabund(pa)
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manyglm(formula=data~as.factor(alloc), family="binomial")$AICsum
}
pa.AICsum = data.frame(AICsum,groups)

# scores
data = data.frame(lapply(scores, as.factor))
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manyclm.AICsum(data, as.factor(alloc))
}
scores.AICsum = data.frame(AICsum,groups)

# cover
#data = ((cover/100)*(nrow(cover)-1)+0.5)/nrow(cover) # this is a transformation for betareg
data = cover/100
# use zero/one inflated beta
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manybeta.AICsum(data, as.factor(alloc))
}
cover.AICsum.infbe = data.frame(AICsum,groups)

# use compund Poisson-gamma (i.e. tweedie with power= 1.5)
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  Xdata = data.frame(alloc=cutree(clust, groups[i]), dummy=integer(nrow(data)))
  AICsum[i] = sum(AIC(manyany("glm", data, data~alloc, data=Xdata,
                              family=tweedie(var.power=1.2, link.power=0), var.power=1.2, composition=FALSE)))
}
cover.AICsum.tw = data.frame(AICsum,groups)

# use ordinal regression
data = data.frame(lapply(cover, as.factor))
AICsum = numeric(length(groups))
for (i in 1:length(groups)) {
  alloc = cutree(clust, groups[i])
  AICsum[i] = manyclm.AICsum(data, as.factor(alloc))
}
cover.AICsum.ord = data.frame(AICsum,groups)

# plot AIC
#CairoWin()

CairoPDF(file="figures/sumAIC_real.pdf", width=12, height=8)
par(mfrow=c(2,3), cex.lab=1.5, cex.main=1.7, oma=c(4,4,0,0)+0.1, mar=c(2,0,4,3)+0.1)

plot(counts.AICsum$AICsum~counts.AICsum$groups, main="counts - negative binonial", 
     ylab="", xlab="", type='l', yaxt='n')
points(x=WhichMin(counts.AICsum)[2],y=WhichMin(counts.AICsum)[1], pch=10, cex=2)

plot(pa.AICsum$AICsum~pa.AICsum$groups, main="pres/abs - binomal", 
     ylab="", xlab="", type='l', yaxt='n')
points(x=WhichMin(pa.AICsum)[2],y=WhichMin(pa.AICsum)[1], pch=10, cex=2)

plot(scores.AICsum$AICsum~scores.AICsum$groups, main="ordinal scores - prop. odds", 
     ylab="", xlab="", type='l', yaxt='n')
points(x=WhichMin(scores.AICsum)[2],y=WhichMin(scores.AICsum)[1], pch=10, cex=2)

plot(cover.AICsum.infbe$AICsum~cover.AICsum.infbe$groups, main="% cover - inf. beta", 
     ylab="", xlab="", type='l', yaxt='n')
points(x=WhichMin(cover.AICsum.infbe)[2],y=WhichMin(cover.AICsum.infbe)[1], pch=10, cex=2)

plot(cover.AICsum.tw$AICsum~cover.AICsum.tw$groups, main="% cover - Poisson-gamma", 
     ylab="", xlab="", type='l', yaxt='n')
points(x=WhichMin(cover.AICsum.tw)[2],y=WhichMin(cover.AICsum.tw)[1], pch=10, cex=2)

plot(cover.AICsum.ord$AICsum~cover.AICsum.ord$groups, main="% cover - prop. odds", 
     ylab="", xlab="", type='l', yaxt='n')
points(x=WhichMin(cover.AICsum.ord)[2],y=WhichMin(cover.AICsum.ord)[1], pch=10, cex=2)

title(xlab = "number of clusters",
      ylab = "sum-of-AIC",
      outer = TRUE, line = 3)

dev.off()



# characteristic species -------------------------------------------------------

alloc = cutree(counts.clust, 8)

# counts
data = mvabund(counts)
counts.manyglm = manyglm(formula=data~as.factor(alloc), family="negative.binomial")
counts.manyglm.null = manyglm(formula=data~1, family="negative.binomial")
counts.manyglm.species = data.frame(sort(counts.manyglm.null$aic-counts.manyglm$aic, decreasing=T))
# clean
counts.manyglm.species = data.frame(species=row.names(counts.manyglm.species), 
                                    dAIC=counts.manyglm.species$sort.counts.manyglm.null.aic)

# scores
data = data.frame(lapply(scores, as.factor))
scores.manyclm = manyclm.AIC(data,as.factor(alloc))
scores.manyclm.null = manyclm.AIC(data,1)
scores.manyclm.species = data.frame(sort(scores.manyclm.null-scores.manyclm, decreasing=T))
# clean
scores.manyclm.species = data.frame(species=row.names(scores.manyclm.species),
                                    dAIC=scores.manyclm.species$sort.scores.manyclm.null)

# check top contributors
top20 = cbind(counts.manyglm.species[1:20,], scores.manyclm.species[1:20,])
sum(counts.manyglm.species[1:100,1] %in% scores.manyclm.species[1:100,1])
write.csv(top20, file="figures/top20_CountsScores.csv", row.names=F)


# indicator values (fidelity and abundance)
indicatorValue = indval(x=counts, clustering=alloc)

high.indval = list()
for (i in 1:8) {
  vals = sort(indicatorValue$relfrq[[i]], decreasing=T)[1:20]
  high.indval[[i]] = names(counts)[which(indicatorValue$relfrq[[i]] %in% vals[vals>0])]
  rm(vals)
}
high.indval = unlist(high.indval)
high.indval = unique(high.indval)

# check match-up between top dAIC species and high indval species
counts.shared = numeric(15)
scores.shared = numeric(15)
for (i in 1:15) {
  counts.shared[i] = sum(counts.manyglm.species[1:(i*10),1] %in% high.indval)
  scores.shared[i] = sum(scores.manyclm.species[1:(i*10),1] %in% high.indval)
}

# indval top20
counts.shared.20 = counts.shared
scores.shared.20 = scores.shared
# indval top10
counts.shared.10 = counts.shared
scores.shared.10 = scores.shared
# relfrq top20
counts.shared.rf = counts.shared
scores.shared.rf = scores.shared
# relabund top20
counts.shared.ra = counts.shared
scores.shared.ra = scores.shared


#CairoWin()
CairoPDF(file="figures/AIC_indval.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(4,4,2,0.2), cex.axis=0.7, cex.main=0.9, cex.sub=0.7, cex.lab=1)
plot(y=scores.shared.20, x=1:15*10, ylab="shared species", 
     xlab=expression(paste("top ",Delta,"  AIC species considered")), 
     main=expression(paste(Delta,"  AIC and 'indicator value' match-up")), pch=1)
points(y=counts.shared.20, x=1:15*10, pch=16)
points(y=scores.shared.10, x=1:15*10, pch=1)
points(y=counts.shared.10, x=1:15*10, pch=16)

plot(y=scores.shared.rf, x=1:15*10, ylab="", 
     xlab=expression(paste("top ",Delta,"  AIC species considered")), 
     main=expression(paste(Delta,"  AIC and relative freq./abund. match-up")), pch=2)
points(y=counts.shared.rf, x=1:15*10, pch=17)
points(y=scores.shared.ra, x=1:15*10, pch=5)
points(y=counts.shared.ra, x=1:15*10, pch=18)
dev.off()



# sum-of-aic backward selection -------------------------------------------

# clusters
alloc.1 = cutree(counts.clust, 20)

# pairwise comparison based merging
# pairwise manyglm based on alloc.1 (clustering with 20 groups)
# we do this manually here to describe the process - but for operational use, just automate merging the lowest delta AIC pairs
pairs = GeneratePairwise(alloc.1)
PairwiseCombinedSpeciesCount(counts, alloc.1, pairs)
alloc.1.pairwise = PairwiseManyglm(counts, alloc.1, pairs)
alloc.1.pairwise.sorted = alloc.1.pairwise[with(alloc.1.pairwise, order(-deltaSumAIC)),]
# calculate summary and merge clusters using this information
PairwiseSummary(alloc.1, alloc.1.pairwise)
alloc.1.1 = alloc.1
alloc.1.1[alloc.1.1==7 | alloc.1.1==17] = 21

# repeat
pairs = GeneratePairwise(alloc.1.1)
PairwiseCombinedSpeciesCount(counts, alloc.1.1, pairs)
alloc.1.1.pairwise = PairwiseManyglm(counts, alloc.1.1, pairs)
PairwiseSummary(alloc.1.1, alloc.1.1.pairwise)
alloc.1.2 = alloc.1.1
alloc.1.2[alloc.1.2==21 | alloc.1.2==14] = 22

# repeat
pairs = GeneratePairwise(alloc.1.2)
PairwiseCombinedSpeciesCount(counts, alloc.1.2, pairs)
alloc.1.2.pairwise = PairwiseManyglm(counts, alloc.1.2, pairs)
PairwiseSummary(alloc.1.2, alloc.1.2.pairwise)
alloc.1.3 = alloc.1.2
alloc.1.3[alloc.1.3==22 | alloc.1.3==18] = 23

# repeat
pairs = GeneratePairwise(alloc.1.3)
PairwiseCombinedSpeciesCount(counts, alloc.1.3, pairs)
alloc.1.3.pairwise = PairwiseManyglm(counts, alloc.1.3, pairs)
PairwiseSummary(alloc.1.3, alloc.1.3.pairwise)
alloc.1.4 = alloc.1.3
alloc.1.4[alloc.1.4==5 | alloc.1.4==10] = 24

# repeat
pairs = GeneratePairwise(alloc.1.4)
PairwiseCombinedSpeciesCount(counts, alloc.1.4, pairs)
alloc.1.4.pairwise = PairwiseManyglm(counts, alloc.1.4, pairs)
PairwiseSummary(alloc.1.4, alloc.1.4.pairwise)
alloc.1.5 = alloc.1.4
alloc.1.5[alloc.1.5==23 | alloc.1.5==15] = 25

# repeat
pairs = GeneratePairwise(alloc.1.5)
PairwiseCombinedSpeciesCount(counts, alloc.1.5, pairs)
alloc.1.5.pairwise = PairwiseManyglm(counts, alloc.1.5, pairs)
PairwiseSummary(alloc.1.5, alloc.1.5.pairwise)
alloc.1.6 = alloc.1.5
alloc.1.6[alloc.1.6==1 | alloc.1.6==4] = 26

# repeat
pairs = GeneratePairwise(alloc.1.6)
PairwiseCombinedSpeciesCount(counts, alloc.1.6, pairs)
alloc.1.6.pairwise = PairwiseManyglm(counts, alloc.1.6, pairs)
PairwiseSummary(alloc.1.6, alloc.1.6.pairwise)
alloc.1.7 = alloc.1.6
alloc.1.7[alloc.1.7==25 | alloc.1.7==16] = 27

# repeat
pairs = GeneratePairwise(alloc.1.7)
PairwiseCombinedSpeciesCount(counts, alloc.1.7, pairs)
alloc.1.7.pairwise = PairwiseManyglm(counts, alloc.1.7, pairs)
PairwiseSummary(alloc.1.7, alloc.1.7.pairwise)
alloc.1.8 = alloc.1.7
alloc.1.8[alloc.1.8==9 | alloc.1.8==20] = 28

# repeat
pairs = GeneratePairwise(alloc.1.8)
PairwiseCombinedSpeciesCount(counts, alloc.1.8, pairs)
alloc.1.8.pairwise = PairwiseManyglm(counts, alloc.1.8, pairs)
PairwiseSummary(alloc.1.8, alloc.1.8.pairwise)
alloc.1.9 = alloc.1.8
alloc.1.9[alloc.1.9==28 | alloc.1.9==19] = 29

# repeat
pairs = GeneratePairwise(alloc.1.9)
PairwiseCombinedSpeciesCount(counts, alloc.1.9, pairs)
alloc.1.9.pairwise = PairwiseManyglm(counts, alloc.1.9, pairs)
PairwiseSummary(alloc.1.9, alloc.1.9.pairwise)
alloc.1.10 = alloc.1.9
alloc.1.10[alloc.1.10==27 | alloc.1.10==13] = 30

# repeat
pairs = GeneratePairwise(alloc.1.10)
PairwiseCombinedSpeciesCount(counts, alloc.1.10, pairs)
alloc.1.10.pairwise = PairwiseManyglm(counts, alloc.1.10, pairs)
PairwiseSummary(alloc.1.10, alloc.1.10.pairwise)
alloc.1.11 = alloc.1.10
alloc.1.11[alloc.1.11==2 | alloc.1.11==11] = 31

# repeat
pairs = GeneratePairwise(alloc.1.11)
PairwiseCombinedSpeciesCount(counts, alloc.1.11, pairs)
alloc.1.11.pairwise = PairwiseManyglm(counts, alloc.1.11, pairs)
PairwiseSummary(alloc.1.11, alloc.1.11.pairwise)
alloc.1.12 = alloc.1.11
alloc.1.12[alloc.1.12==30 | alloc.1.12==29] = 32

# fit and get AIC
data = mvabund(counts)
alloc.pairwise.AIC = data.frame(alloc=c("alloc.1.1","alloc.1.2","alloc.1.3","alloc.1.4","alloc.1.5","alloc.1.6",
                                        "alloc.1.7","alloc.1.8","alloc.1.9","alloc.1.10","alloc.1.11","alloc.1.12"),
                                groups=c(19:8),
                                AICsum=c(manyglm(formula=data~as.factor(alloc.1.1), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.2), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.3), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.4), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.5), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.6), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.7), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.8), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.9), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.10), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.11), family="negative.binomial")$AICsum,
                                         manyglm(formula=data~as.factor(alloc.1.12), family="negative.binomial")$AICsum))

alloc.pairwise.AIC$AICzero = alloc.pairwise.AIC$AICsum-min(alloc.pairwise.AIC$AICsum)

counts.AICsum.pairwise = rbind(counts.AICsum, alloc.pairwise.AIC[,c(3,2)])


# plot hclust vs. pairwise modelling merges
AIC.hclust = counts.AICsum.pairwise[1:19,]
AIC.paiwise = counts.AICsum.pairwise[19:31,]

CairoPDF(file="figures/ClassRefine.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(AIC.hclust$AICsum~AIC.hclust$groups, type='l', 
     ylim=c(min(AIC.paiwise$AICsum),max(AIC.hclust$AICsum)),
     ylab="sum-of-AIC", yaxt='n', xlab="number of clusters", main="")
points(x=WhichMin(AIC.hclust)[2],y=WhichMin(AIC.hclust)[1], pch=10, cex=2)
lines(y=AIC.paiwise$AICsum, x=AIC.paiwise$groups, lty=3)
points(x=WhichMin(AIC.paiwise)[2],y=WhichMin(AIC.paiwise)[1], pch=10, cex=2)
legend('bottomright', legend=c("hclust", "model-based"), lty=c(1,3))
dev.off()



# sum-of-AIC for simulated abundance --------------------------------------

# counts
# loop through different cutoff levels and simulated data sets to test sum-of-AIC on groupings
# set cluster numbers and cluster source
groups = 2:20
# start nboot runs
counts.sim.AICsum = list()
nboot = 3
sim.clust.AICsum = numeric(nboot)
for (sim in 1:nboot) {
  # source vegSim.R/run VegSim() to generate some synthetic data
  sampled = VegSim(C=9, S=50, Nmult=50, width=0.25, uniformWidth=F, randomSample=F, plotSim=F, plotOpt=F)
  counts.sim = sampled[,-1]
  #clusters.sim = sampled[,1]
  # cluster the simulated data
  clust = hclust(vegdist(counts.sim, "bray"), "average")
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

# save for narrow niche
AICsum.narrow = counts.sim.AICsum
true.narrow = sim.clust.AICsum
#### rerun sims above before continuing - bad coding, forgive me please =(
# save for wide niche
AICsum.wide = counts.sim.AICsum
true.wide = sim.clust.AICsum

# plots
#CairoWin()
CairoPDF(file="figures/sumAIC_simulated.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(4,4,4,0.5), cex.lab=1.05, cex.main=0.95)

maxAIC = numeric(length(AICsum.narrow))
minAIC = numeric(length(AICsum.narrow))
for (i in 1:length(AICsum.narrow)) { maxAIC[i]=(max(c(AICsum.narrow[[i]]$AICsum,true.narrow)))}
for (i in 1:length(AICsum.narrow)) { minAIC[i]=(min(c(AICsum.narrow[[i]]$AICsum,true.narrow)))}
plot(AICsum.narrow[[1]]$AICsum~AICsum.narrow[[1]]$groups, type="n",
     ylim=c(min(minAIC),max(maxAIC)), ylab="sum-of-AIC", yaxt='n', 
     xlab="number of clusters", main="Narrow niches, high membership prob.")
for (i in 1:length(AICsum.narrow)) {
  col=i+1
  lines(AICsum.narrow[[i]]$AICsum~AICsum.narrow[[i]]$groups, col=col)
  points(x=WhichMin(AICsum.narrow[[i]])[2],y=WhichMin(AICsum.narrow[[i]])[1], pch=10)
  points(x=9, y=true.narrow[i], pch=9, col=col, cex=1.5)
}

maxAIC = numeric(length(AICsum.wide))
minAIC = numeric(length(AICsum.wide))
for (i in 1:length(AICsum.wide)) { maxAIC[i]=(max(c(AICsum.wide[[i]]$AICsum,true.wide)))}
for (i in 1:length(AICsum.wide)) { minAIC[i]=(min(c(AICsum.wide[[i]]$AICsum,true.wide)))}
plot(AICsum.wide[[1]]$AICsum~AICsum.wide[[1]]$groups, type="n",
     ylim=c(min(minAIC),max(maxAIC)), ylab="", yaxt='n', 
     xlab="number of clusters", main="Wide niches, low membership prob.")
for (i in 1:length(AICsum.wide)) {
  col = i+1
  lines(AICsum.wide[[i]]$AICsum~AICsum.wide[[i]]$groups, col=col)
  points(x=WhichMin(AICsum.wide[[i]])[2],y=WhichMin(AICsum.wide[[i]])[1], pch=10)
  points(x=9, y=true.wide[i], pch=9, col=col, cex=1.5)
}
dev.off()


# dAIC and indicator values
# create simulated data
sampled = VegSim(C=9, S=20, Nmult=500, width=1, uniformWidth=F, customMultiplier=c(0.1,2.5), 
                 randomSample=F, plotSim=F, plotsToPrint=NULL, plotOpt=T)
counts.sim = data.frame(sampled[,-1])
alloc = sampled[,1]

# calculate dAIC
data = mvabund(counts.sim)
#meanvar.plot(data~factor(alloc))
counts.manyglm = manyglm(formula=data~as.factor(alloc), family="negative.binomial")
counts.manyglm.null = manyglm(formula=data~1, family="negative.binomial")
counts.manyglm.species = data.frame(sort(counts.manyglm.null$aic-counts.manyglm$aic, decreasing=T))
# clean
counts.manyglm.species = data.frame(species=row.names(counts.manyglm.species), 
                                    dAIC=counts.manyglm.species$sort.counts.manyglm.null.aic)

dAIC = join(x=counts.manyglm.species, y=data.frame(species=colnames(counts.sim), VarTol=attr(sampled,"VarTol")), by="species")

# calculate indvals
indicatorValue = indval(x=counts.sim, clustering=alloc)

indvals = list()
for (i in 1:attr(sampled,"C")) {
  vals = indicatorValue$indval[paste0(i)]
  indvals.df = data.frame(species=row.names(vals), indval=vals[,1])
  indvals[[i]] = indvals.df[indvals.df$indval>0.1,]
}
indvals = do.call(rbind, indvals)
indvals = join(x=indvals, y=data.frame(species=colnames(counts.sim), VarTol=attr(sampled,"VarTol")), by="species")


# plot - compare high dAIC or indval to niche width 
library(gridExtra)
niches = as.character(unique(round(attr(sampled,"VarTol"),1)))
niches[seq(2,length(niches),2)] = ""

dAIC$VarTol = factor(dAIC$VarTol)
indvals$VarTol = factor(indvals$VarTol)

CairoPDF(file="figures/niche_simulations.pdf", width=5, height=8)
#CairoWin()
niche.dAIC = ggplot(aes(y=dAIC, x=VarTol), data=dAIC) + geom_boxplot(lwd=0.3, outlier.size=0.5) +
  scale_x_discrete(labels=niches) + xlab("Species niche width") + ylab(expression(paste(Delta,"  AIC"))) +
  theme_classic() + theme(axis.line=element_line(size=0.4)) +
  geom_hline(yintercept=0, linetype="dashed", lwd=0.3)

niche.indval = ggplot(aes(y=indval, x=VarTol), data=indvals) + geom_boxplot(lwd=0.3, outlier.size=0.5) + 
  scale_x_discrete(labels=niches) + xlab("Species niche width") + ylab("'indval'") + 
  theme_classic() + theme(axis.line=element_line(size=0.4)) +
  geom_hline(yintercept=0.1, linetype="dashed", lwd=0.3)

grid.arrange(niche.dAIC, niche.indval, ncol=1)

dev.off()



# residual plots for supp material ----------------------------------------

data = mvabund(counts)
pois = manyglm(formula=data~as.factor(counts.alloc), family="poisson")
nb = manyglm(formula=data~as.factor(counts.alloc), family="negative.binomial")

data = mvabund(pa)
pa = manyglm(formula=data~as.factor(counts.alloc), family="binomial")

# ord = clm(as.factor(scores$Pomaumbe)~as.factor(counts.alloc)) # don't know how to define residuals for ordinal regression

data = cover/100
Xdata = data.frame(alloc=counts.alloc, dummy=integer(nrow(data)))
tw = manyany("glm", data, data~as.factor(alloc), data=Xdata, family=tweedie(var.power=1.2, link.power=0), var.power=1.2, composition=FALSE)
infbeta = gamlss(data$Eucaalbe~as.factor(counts.alloc), family=BEZI())

Mean = colwise(function(x){mean(x)})(cover)
Mean = as.numeric(t(mv.m))
Variance = colwise(function(x){var(x)})(cover)
Variance = as.numeric(t(mv.v))

CairoPDF(file="figures/supp_resids.pdf", width=8, height=6)
par(mfrow=c(2,3))

plot(log(Variance)~log(Mean), main="mean-variance, no levels")
abline(a=0, b=1, lty=2)
plot(pois$residuals[,'Eucaalbe']~jitter(pois$linear.predictor[,'Eucaalbe'],amount=0.1), 
     main="counts - Poisson", xlab="fitted", ylab="residual")
plot(nb$residuals[,'Eucaalbe']~jitter(pois$linear.predictor[,'Eucaalbe'],amount=0.1), 
     main="counts - negative binomial", xlab="fitted", ylab="residual")
plot(pa$residuals[,'Eucaalbe']~jitter(pois$linear.predictor[,'Eucaalbe'],amount=0.1), 
     main="pres/abs - binomial", xlab="fitted", ylab="residual")
plot(pois$residuals[,'Eucaalbe']~jitter(pois$linear.predictor[,'Eucaalbe'],amount=0.1), 
     main="cover - Poisson-gamma", xlab="fitted", ylab="residual")
plot(infbeta$residuals~jitter(fitted(infbeta),amount=0.005),
     main="cover - inflated beta", xlab="fitted", ylab="residual")

dev.off()

