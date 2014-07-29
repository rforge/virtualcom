div.param.native <- function(spSite, niche.opt, tree, phy, fun, nrep = 100, null.model = "taxa.labels", indX.nat) {
    # get observed values
    obs <- sapply(div.param.native.obs(spSite = spSite, phy = phy, fun = fun, niche.optima = niche.opt, tree = tree, indX.nat = indX.nat), function(z) z)
    if (is.null(nrow(obs))) {
        obs <- t(as.data.frame(obs))
    }
    # prepare randomizations for null models
    if (!is.null(null.model)) {
        spSiteNULL <- switch(null.model, taxa.labels = replicate(nrep, as.matrix(spSite[, sample(1:ncol(spSite))])), richness = replicate(nrep, randomizeMatrix(spSite, 
            null.model = "richness")), frequency = replicate(nrep, randomizeMatrix(spSite, null.model = "frequency")), independentswap = replicate(nrep, randomizeMatrix(spSite, 
            null.model = "independentswap")), trialswap = replicate(nrep, randomizeMatrix(spSite, null.model = "trialswap")))
        # calculate null model observations
        tmp <- lapply(1:nrep, function(x) div.param.native.obs(spSiteNULL[, , x], phy, fun, niche.optima = niche.opt, tree = tree, indX.nat = indX.nat))
        # get output for null model observations
        tmp1 <- lapply(1:nrow(obs), function(z) lapply(1:length(colnames(obs)), function(y) sapply(1:nrep, function(x) tmp[[x]][[y]][[z]])))
        meanNULL <- t(sapply(tmp1, function(z) sapply(z, mean, na.rm = TRUE)))
        sdNULL <- t(sapply(tmp1, function(z) sapply(z, sd, na.rm = TRUE)))
        rankNULL <- t(sapply(1:nrow(obs), function(z) sapply(1:length(colnames(obs)), function(y) rank(c(obs[z, y], tmp1[[z]][[y]]))[1]/(nrep + 1))))
        zNULL <- (obs - meanNULL)/sdNULL  # zNULL = NaN when sdNULL=0...
        colnames(zNULL) <- colnames(rankNULL) <- colnames(meanNULL) <- colnames(sdNULL) <- colnames(obs)
    } else {
        zNULL <- rankNULL <- meanNULL <- sdNULL <- NULL
    }
    return(list(obs = obs, zNULL = zNULL, rankNULL = rankNULL, meanNULL = meanNULL, sdNULL = sdNULL))
} 