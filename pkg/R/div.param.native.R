div.param.native <- function(spSite, niche.opt, tree, phy, fun, nrep=100, 
					null.model = "taxa.labels", indX.nat){
						
    # to test: spSite=all.abundances2; niche.opt=niche.optima.nat; tree=tree.nat; phy=dist.phy.nat; fun=dist.fun.nat; nrep=n.rep.null.model; null.model = null.model; indX.nat=indX.nat; spSite=all.abundances2; niche.opt=niche.optima.nat; tree=tree.nat; phy=dist.phy.nat; fun=dist.fun.nat; nrep=n.rep.null.model; null.model = null.model; indX.nat=indX.nat

    # get observed values
    obs <- sapply(div.param.native.obs(spSite=spSite, phy=phy, fun=fun, niche.optima=niche.opt, tree=tree, indX.nat=indX.nat), function(z) z)
    if (is.null(nrow(obs))){obs<-t(as.data.frame(obs))}
	
    # prepare randomizations for null models
	if(!is.null(null.model)){
	    #  phyNULL <- lapply(1:nrep, function(x) taxaShuffle(phy))
	    #  funNULL <- lapply(1:nrep, function(x) taxaShuffle(fun))
	
	    spSiteNULL <- switch(null.model, 
	    		taxa.labels = replicate(nrep, as.matrix(spSite[,sample(1:ncol(spSite))])),
	      		richness = replicate(nrep, randomizeMatrix(spSite, null.model = "richness")),
	      		frequency = replicate(nrep, randomizeMatrix(spSite, null.model = "frequency")),
	      		independentswap = replicate(nrep, randomizeMatrix(spSite, null.model = "independentswap")),
	      		trialswap = replicate(nrep, randomizeMatrix(spSite, null.model = "trialswap")))
	    
	    # calculate null model observations 
	    tmp <- lapply(1:nrep, function(x) div.param.native.obs(spSiteNULL[,,x], phy, fun, niche.optima=niche.opt, tree=tree, indX.nat=indX.nat))           
	    
	    # Testing imbalance of local trees
	    # colless.indices <- apply(spSite, 1, function(x){
	    #  tree.red <- drop.tip(tree, as.character(tree$tip.label[x==0]))
	    #  col.test <- colless.test.no.print(as.treeshape(tree.red), alternative="greater", n.mc=nrep)$p.value
	    # })
	      
	    # get output for null model observations
	    tmp1 <- lapply(1:nrow(obs), function(z) lapply(1:length(colnames(obs)), function(y) sapply(1:nrep, function(x) tmp[[x]][[y]][[z]])))
	    meanNULL <- t(sapply(tmp1, function(z) sapply(z, mean, na.rm = TRUE)))
	    sdNULL <- t(sapply(tmp1, function(z) sapply(z, sd, na.rm = TRUE)))
	    rankNULL <- t(sapply(1:nrow(obs), function(z) sapply(1:length(colnames(obs)), function(y) rank(c(obs[z,y],tmp1[[z]][[y]]))[1]/(nrep+1))))
	    zNULL <- (obs - meanNULL) / sdNULL		# zNULL = NaN when sdNULL=0...
	    colnames(zNULL) <- colnames(rankNULL) <- colnames(meanNULL) <- colnames(sdNULL) <- colnames(obs)
	} else {
    	zNULL <- rankNULL <- meanNULL <- sdNULL <- NULL
    }
    return(list(obs=obs, zNULL=zNULL, rankNULL=rankNULL, meanNULL=meanNULL, sdNULL=sdNULL))
}
