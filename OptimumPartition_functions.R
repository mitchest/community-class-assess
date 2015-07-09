## Author: Mitchell Lyons (mitchell.lyons@gmail.com)
## No licence, though please attribute me and publish any derivitives open source =)


# function definitions ----------------------------------------------------


# if you want to use these function more stand-alone, then add require() or library() calls in as needed


# function for calcualting minimum values for plotting, from sum-of-AIC elements
WhichMin = function(x) {
  idx = which(x==min(x$AICsum))
  return(c(x$AICsum[idx],x$groups[idx]))
}

manyclm.AICsum = function(response, predictor) {
  AIC = numeric(nrow(response))
  for (i in 1:ncol(response)) {
    AIC[i] = AIC(clm(response[,i]~predictor))
  }
  return(sum(AIC))
}

manyclm.AIC = function(response, predictor) {
  if (predictor==1) { # null model
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
  return(data.frame(target=target.com, test=test.com))
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
    nspecies[i] = ncol(combined.multivar)
    deltaSumAIC[i] = fit.null$AICsum - fit$AICsum
  }
  return(cbind(pairs, rank, nspecies, deltaSumAIC))
}

## function that calculates lowest ranked pairwise comparison and lowest ranked cluster
# alloc==allocation used, pairwisemanyglm=object from PairwiseManyglm() using alloc
PairwiseSummary = function(alloc, pairwisemanyglm) {
  # print cluster numbers
  print(table(alloc))
  # find lowest ranked pairwise comparison
  print(pairwisemanyglm[pairwisemanyglm$deltaSumAIC==min(pairwisemanyglm$deltaSumAIC),])
  # find cluster that has lowest mean ranks
  for (i in unique(alloc)) {
    print(paste0("Cluster ",i," mean delta sum-of-AIC:"))
    print(mean(pairwisemanyglm$deltaSumAIC[pairwisemanyglm$target==i | pairwisemanyglm$test==i]))
  }
}