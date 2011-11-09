div.param.native <- function(spSite, phy, fun, nrep=100, null.model = c("taxa.labels", "richness",
    "frequency", "sample.pool", "phylogeny.pool", "independentswap",
    "trialswap","constrained"),suit=NULL,taxa=NULL){
    # to test: spSite=output$native$communities; phy=Cophi; fun=Cofun; nrep=100; null.model = "taxa.labels"
    # get observed values
    obs <- sapply(div.param.native.obs(spSite, phy, fun), function(z) z)

    # prepare randomizations for null models
    null.model <- match.arg(null.model)
    phyNULL <- lapply(1:nrep, function(x) taxaShuffle(phy))
    funNULL <- lapply(1:nrep, function(x) taxaShuffle(fun))
    
    spSiteNULL <- switch(null.model, 
      		taxa.labels = replicate(nrep, spSite),
      		richness = replicate(nrep, randomizeMatrix(spSite,null.model = "richness")),
      		frequency = replicate(nrep, randomizeMatrix(spSite,null.model = "frequency")),
      		independentswap = replicate(nrep, randomizeMatrix(spSite,null.model = "independentswap")),
      		trialswap = replicate(nrep, randomizeMatrix(spSite,null.model = "trialswap")),
          constrained = replicate(nrep,constrained.NM(spSite,taxa,suit)))
    
    # calculate null model observations            
    tmp <- switch(null.model,
      		taxa.labels = lapply(1:nrep, function(x) div.param.native.obs(spSite, phyNULL[[x]], funNULL[[x]])),
      		richness = lapply(1:nrep, function(x) div.param.native.obs(spSiteNULL[[x]], phy, fun)),
      		frequency = lapply(1:nrep, function(x) div.param.native.obs(spSiteNULL[[x]], phy, fun)),
      		independentswap = lapply(1:nrep, function(x) div.param.native.obs(spSiteNULL[[x]], phy, fun)),
      		trialswap = lapply(1:nrep, function(x) div.param.native.obs(spSiteNULL[[x]], phy, fun))
          constrained = lapply(1:nrep, function(x) div.param.native.obs(spSiteNULL[[x]], phy, fun)))
      
    # get output for null model observations
    tmp1 <- lapply(1:nrow(obs), function(z) lapply(1:length(colnames(obs)), function(y) sapply(1:nrep, function(x) tmp[[x]][[y]][[z]])))
    meanNULL <- t(sapply(tmp1, function(z) sapply(z, mean, na.rm = TRUE)))
    sdNULL <- t(sapply(tmp1, function(z) sapply(z, sd, na.rm = TRUE)))
    rankNULL <- t(sapply(1:nrow(obs), function(z) sapply(1:length(colnames(obs)), function(y) rank(c(obs[z,y],tmp1[[z]][[y]]))[1]/(nrep+1))))
    zNULL <- (obs - meanNULL) / sdNULL		# zNULL = NaN when sdNULL=0...
    colnames(zNULL) <- colnames(rankNULL) <- colnames(meanNULL) <- colnames(sdNULL) <- colnames(obs)
    return(list(obs=obs, zNULL=zNULL, rankNULL=rankNULL, meanNULL=meanNULL, sdNULL=sdNULL))
}
