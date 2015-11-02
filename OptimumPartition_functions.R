## Author: Mitchell Lyons (mitchell.lyons@gmail.com)
## No licence, though please attribute me and publish any derivitives open source =)


# function definitions ----------------------------------------------------

# function for calcualting minimum values for plotting, from sum-of-AIC elements
WhichMin = function(x) {
  idx = which(x==min(x$AICsum))
  return(c(x$AICsum[idx],x$groups[idx]))
}

minus.min = function(x) {
  x$AICsum - min(x$AICsum)
}

manyclm.AICsum = function(response, predictor) {
  AIC = numeric(nrow(response))
  for (i in 1:ncol(response)) {
    AIC[i] = AIC(clm(response[,i]~predictor))
  }
  return(sum(AIC))
}

manyclm.AIC = function(response, predictor) {
  if (length(predictor)==1) { # null model
    AIC = numeric(nrow(response))
    species = character(nrow(response))
    for (i in 1:ncol(response)) {
      AIC[i] = AIC(clm(response[,i]~1))
      species[i] = names(response)[i]
    }
  } else {
    AIC = numeric(nrow(response))
    species = character(nrow(response))
    for (i in 1:ncol(response)) {
      AIC[i] = AIC(clm(response[,i]~predictor))
      species[i] = names(response)[i]
    }
  }
  names(AIC) = species
  return(AIC)
}

manybeta.AICsum = function(response, predictor) {
  AIC = numeric(nrow(response))
  for (i in 1:ncol(response)) {
    AIC[i] = gamlss(data[,i]~predictor, family=BEZI())$aic
    print(cat("------", "\n", i, "\n","------"))
  }
  return(sum(AIC))
}

GeneratePairwise = function(communities) {
  # generate pairwise tests
  communities = as.character(unique(communities))
  target.com = character(round(length(communities)^2/2))
  test.com = character(round(length(communities)^2/2))
  # loop across communities to create pairwise comparisons, without reverse tests
  ticker=0
  for (target in 1:(length(communities)-1)){
    for (test in target:(length(communities)-1)){
      ticker = ticker+1
      target.com[ticker] = communities[target]
      test.com[ticker] = communities[test+1]
    }
  }
  target.com=target.com[!target.com==""]
  test.com=test.com[!test.com==""]
  return(data.frame(target=target.com, test=test.com, stringsAsFactors=F))
}

## function that calculates the species count after removing species with 0 occurance in a pariwise combination
PairwiseCombinedSpeciesCount = function(multivar, predictor, pairs) {
  target = as.character(pairs$target)
  test = as.character(pairs$test)
  for (i in 1:length(target)){
    # subset to target sites
    target.rows = which(predictor==target[i])
    target.multivar = multivar[target.rows,]
    target.pred = predictor[target.rows]
    # subset to test sites
    test.rows = which(predictor==test[i])
    test.multivar = multivar[test.rows,]
    test.pred = predictor[test.rows]
    # create data objects for the pariwise test
    combined.multivar = rbind(target.multivar, test.multivar)
    combined.multivar = combined.multivar[,colSums(combined.multivar)>0] # remove species that don't occur in either community
    print(ncol(combined.multivar))
  }
}

PairwiseManyglm = function(multivar, predictor, pairs) {
  target = as.character(pairs$target)
  test = as.character(pairs$test)
  rank = numeric(nrow(pairs))
  nspecies = numeric(nrow(pairs))
  deltaSumAIC = numeric(nrow(pairs))
  for (i in 1:length(target)){
    # subset to target sites
    target.rows = which(predictor==target[i])
    target.multivar = multivar[target.rows,]
    target.pred = predictor[target.rows]
    # subset to test sites
    test.rows = which(predictor==test[i])
    test.multivar = multivar[test.rows,]
    test.pred = predictor[test.rows]
    # create data objects for the pariwise test
    combined.multivar = rbind(target.multivar, test.multivar)
    combined.multivar = mvabund(data.matrix(combined.multivar))
    #combined.multivar = mvabund(data.matrix(combined.multivar[,colSums(combined.multivar)>0])) # remove double absences
    combined.pred = as.factor(c(as.character(target.pred), as.character(test.pred)))
    # fit/test the model
    fit = manyglm(combined.multivar ~ combined.pred, family="negative.binomial")
    fit.null = manyglm(combined.multivar ~ 1, family="negative.binomial")
    # calculate % of species for which dAIC (from null) is >n
    dAIC = fit.null$aic-fit$aic
    rank[i] = sum(dAIC>4)/length(dAIC)
    nspecies[i] = ncol(combined.multivar)
    deltaSumAIC[i] = fit.null$AICsum - fit$AICsum
  }
  return(data.frame(pairs, rank, nspecies, deltaSumAIC, stringsAsFactors=F))
}

## function that calculates lowest ranked pairwise comparison and lowest ranked cluster
# alloc==allocation used, pairwisemanyglm=object from PairwiseManyglm() using alloc
PairwiseSummary = function(alloc, pairwisemanyglm) {
  # print cluster numbers
  print(table(alloc))
  # find lowest ranked pairwise comparison
  print(pairwisemanyglm[pairwisemanyglm$deltaSumAIC==min(pairwisemanyglm$deltaSumAIC),])
  # find cluster that has lowest mean ranks
  #   for (i in unique(alloc)) {
  #     print(paste0("Cluster ",i," mean delta sum-of-AIC:"))
  #     print(mean(pairwisemanyglm$deltaSumAIC[pairwisemanyglm$target==i | pairwisemanyglm$test==i]))
  #   }
}


## OptimClass functions

fishers.pval = function(x, alloc, cluster.n) {
  fisher.test(table(x==0,alloc!=cluster.n), alternative="greater")$p.value
}

fishers.pvals = function(species.data, alloc) {
  pvals.list = list()
  for (i in unique(alloc)) {
    cluster.n = i
    pvals = apply(X=species.data, MARGIN=2, FUN=fishers.pval, alloc=alloc, cluster.n=cluster.n)
    names(pvals) = names(species.data)
    pvals.list[[i]] = pvals
  }
  return(pvals.list)
}

extract.faithful = function(x, alpha) {
  names(x[x < alpha])
}

OptimClass1 = function(pvals.allocs.list, alpha) {
  k.clusters = 2:length(pvals.allocs.list)
  n.faithful.species = numeric(length(k.clusters))
  for (i in k.clusters) {
    n.faithful.species[i-1] = length(unique(unlist(
      lapply(X=pvals.allocs.list[[i]], FUN=extract.faithful, alpha=alpha)
    )))
  }
  print(plot(n.faithful.species~k.clusters, ylab="No. faithful species", xlab="No. clusters", 
             main=paste0("OptimClass1 - cut ",alpha), type="l"))
  return(cbind(k.clusters,n.faithful.species))
}

length.faithful = function(x, alpha) {
  length(x[x < alpha])
}

OptimClass2 = function(pvals.allocs.list, alpha, n.faithful) {
  k.clusters = 2:length(pvals.allocs.list)
  clusters.gt.n = numeric(length(k.clusters))
  for (i in k.clusters) {
    clusters.gt.n[i-1] = sum(unlist(
      lapply(X=pvals.allocs.list[[i]], FUN=length.faithful, alpha=alpha)
    )>n.faithful)
  }
  print(plot(clusters.gt.n~k.clusters, ylab=paste0("No. clusters with > ",n.faithful," faithful species"), 
             xlab="No. clusters", main=paste0("OptimClass2 - cut ",alpha), type="l"))
  print(abline(a=0, b=1, lty=2, col="red"))
  return(cbind(k.clusters,clusters.gt.n))
}

OptimClass.species = function(pvals.allocs.list, n.species) {
  k.clusters = 2:length(pvals.allocs.list)
  faithful.species = list()#numeric(length(k.clusters))
  for (i in k.clusters) {
    faithful.species[[i]] = sort(unlist(pvals.allocs.list[[i]]), decreasing=F)[1:n.species]
  }
  return(faithful.species)
}
