# library(coenocliner)
# library(class)


# bivariate community gradient --------------------------------------------

## function to simulate veg communities, where:
# C = 9 --> number of communities - a square number
# S = 3 --> number of species per community
# Nmult = 10 --> multiplier for how many samples in each community
# stddev = 0.01 --> standard deviation for species optimums
# width = 0.25 --> width of the species response curves (tolerence on gradient)
# uniformWidth = T --> if true, species gradient tolerences will be uniform, if false they will vary around "width" (more realistic?)
# customMultiplier = F --> if true
# randomSample = F --> whether to sample close to species optima or randomly across gradients
# plotSim = FALSE --> plot species/community distribution
# plotsToPrint = NULL --> which species numbers to plot response/sample plots for
# plotOpt = TRUE --> plot sample optima
VegSim = function(C=9, S=3, Nmult=10, width=0.25, uniformWidth=T, customMultiplier=NULL, randomSample=F, plotSim=F, plotsToPrint=NULL, plotOpt=T) {
  require(coenocliner)
  require(class)
  # check if parameterised correctly
  if (!sqrt(C)%%1==0) stop("number of communities [C] must be a square number")
  
  # total number of species
  M = S * C
  # number of eventual samples
  N = C * Nmult
  
  # variance of species gradient means
  stddev = 0.01
  # standard deviation for sample optimums to choose community "spread"/"tightness"
  sampleSpread=0.15
  # tolerances if not using uniform tol
  if (is.null(customMultiplier)) {multiplier = rep(seq(0.5,1.5,length.out=S), C)}
  if (!is.null(customMultiplier)) {multiplier = rep(seq(customMultiplier[1],customMultiplier[2],length.out=S), C)}
  tol = rep(width,M) * multiplier
  
  ## First (x) gradient
  minx = 0 # x gradient minimum...
  maxx = sqrt(C)+1 # ...and maximum
  gradx = seq(minx, maxx, length=N) # gradient x locations
  # create species optima on x gradient clustered into communities
  optx = rnorm(M, rep(1:sqrt(C), each=M/sqrt(C)), sd=stddev)
  if (uniformWidth==T) {
    tolx = rep(width, M) # uniform species tolerances
  } else {
    tolx = tol # varying tolerences/width as defined above
  }
  h = ceiling(rlnorm(M, meanlog=3)) # max abundances
  paramx = cbind(opt=optx, tol=tolx, h=h) # put in a matrix
  
  ## Second (y) gradient
  miny = 0 # y gradient minimum...
  maxy = sqrt(C)+1 # ...and maximum
  grady = seq(miny, maxy, length=N) # gradient y locations
  opty = rnorm(M, rep(1:sqrt(C), each=S), sd=stddev) # species optima on y gradient
  if (uniformWidth==T) {
    toly = rep(width, M) # uniform species tolerances
  } else {
    toly = tol # varying tolerences/width as defined above
  }
  paramy = cbind(opt=opty, tol=toly) # put in a matrix
  # combine x and y gradients
  params = list(px=paramx, py=paramy) # put parameters into a list
  
  # generate high-res species/community distribution plots, if sim'd dataset is small enough 
  if (plotSim==T & M<30) {
    # gaussian responses
    grad = expand.grid(x=gradx, y=grady) # put gradient locations together for a large perfect Gaussian model
    responses = coenocline(grad, responseModel = "gaussian",
                       params = params, extraParams = list(corr = 0.5),
                       expectation = TRUE)
    # plot them
    par(mar=c(1,0.5,1,0.5))
    for (i in plotsToPrint) {
      persp(gradx, grady, matrix(responses[, i], ncol=length(grady)),
            ticktype="detailed", zlab="", ylab="", xlab="",
            theta=45, phi=30)
    }
    
    # simulate from gaussian respnse with negbin error
    simulated = coenocline(grad, responseModel = "gaussian",
                        params = params, extraParams = list(corr = 0.5),
                        countModel = "negbin", countParams = list(alpha = 1))
    # plot them
    par(mar=c(1,0.5,1,0.5))
    for (i in plotsToPrint) {
      persp(gradx, grady, matrix(simulated[, i], ncol=length(grady)),
            ticktype="detailed", zlab="", ylab="", xlab="",
            theta=45, phi=30)
    }
  }
  
  if (randomSample==F) {
    # create cluster label vector
    clusters = rep(1:C, each=N/C)
    # sample around the original optimum values to get known communities
    sample.optx = rnorm(N, rep(1:sqrt(C), each=N/sqrt(C)), sd=sampleSpread)
    sample.opty = rnorm(N, rep(rep(1:sqrt(C), each=N/C), sqrt(C)), sd=sampleSpread)
    sample.opt = cbind(x=sample.optx, y=sample.opty)
  } else {
    # sample randomly across the gradient landscape
    sample.optx = runif(n=N,min=minx+width,max=sqrt(C)+width)
    sample.opty = runif(n=N,min=minx+width,max=sqrt(C)+width)
    sample.opt = cbind(x=sample.optx, y=sample.opty)
    # create knn clusters
    centres = cbind(rep(1:sqrt(C), each=sqrt(C)),rep(1:sqrt(C), sqrt(C)))
    class = factor(1:C)
    clusters = knn(train=centres, test=sample.opt, cl=class, k=1)
  }
    
  if (plotOpt == T) {
    plot(x=sample.optx, y=sample.opty, col=clusters, pch=16, xlab="covariate 1", ylab="covariate 2")
    points(x=optx, y=opty, pch=1)
  }
  
  sampled = coenocline(sample.opt, responseModel = "gaussian",
                         params = params, extraParams = list(corr = 0.5),
                         countModel = "negbin", countParams = list(alpha = 1))
  
  sampled = cbind(clusters, sampled)
  # remove rows with no species
  sampled = sampled[rowSums(sampled[,-1])>0,]
  attr(sampled,"S") = S
  attr(sampled,"C") = C
  attr(sampled,"VarTol") = tol
  attr(sampled,"SampleOpts") = sample.opt 
  return(sampled)
}
