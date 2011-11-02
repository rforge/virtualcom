div.param.invasion <- function(spSite, phy, fun, nrep=100, invad, null.model = c("native_inv")){
    # to test: spSite=mysamp; phy=Cophi; fun=Cofun; nrep=100; null.model = "native_inv" ; invad= myinva

    # get observed values
    obs <- div.param.invasion.obs(spSite, phy, fun, invad)

    # prepare randomizations for null models
    null.model <- match.arg(null.model)
    phyNULL <- lapply(1:nrep, function(x) taxaShuffle(phy))
    funNULL <- lapply(1:nrep, function(x) taxaShuffle(fun)) 
 
    spSiteNULL <- switch(null.model, 
        	native_inv = lapply(1:nrep, function(x) invasion_randomization(samp=spSite, inva=invad))        )
    
    # calculate null model observations            
    tmp <- switch(null.model,
      		native_inv = lapply(spSiteNULL, function(x) div.param.invasion.obs(x$samp, phy, fun, x$inv)))
        
    # get output for null model observations
    # so far only for one invader
    obs_1 <- obs[[1]]
    tmpp <- lapply(tmp, function(x) x[[1]])
    tmp1 <- lapply(1:nrow(obs_1), function(z) lapply(1:length(colnames(obs_1)), function(y) sapply(1:nrep, function(x) tmpp[[x]][z,y])))
    meanNULL <- t(sapply(tmp1, function(z) sapply(z, mean, na.rm = TRUE)))
    sdNULL <- t(sapply(tmp1, function(z) sapply(z, sd, na.rm = TRUE)))
    rankNULL <- t(sapply(1:nrow(obs_1), function(z) sapply(1:length(colnames(obs_1)), function(y) ifelse(is.na(obs_1[z,y]),NA,rank(c(obs_1[z,y],tmp1[[z]][[y]]))[1]/(nrep+1)))))	# good luck if you need to understand this line!!
    zNULL <- (obs_1 - meanNULL) / sdNULL
    colnames(zNULL) <- colnames(rankNULL) <- colnames(meanNULL) <- colnames(sdNULL) <- colnames(obs_1)
    return(list(obs=obs_1, zNULL=zNULL, rankNULL=rankNULL, meanNULL=meanNULL, sdNULL=sdNULL))
}
