div.param.invasion <- function(spSite, phy, fun, nrep=100, invad, null.model = c("native_inv")){
	# check the number of simulated communities, and adjust the input format accordingly
    Ncom <- nrow(spSite)
	if(Ncom==1) {spSite <- rbind(spSite, spSite) ; row.names(spSite) <- c(1,2) }						
    # get observed values
    obs <- div.param.invasion.obs(spSite, phy, fun, invad)
    # prepare randomizations for null models
    if(!is.null(null.model)){   	
		# Randomise the species position in the phylogeny 
	    phyNULL <- lapply(1:nrep, function(x) taxaShuffle(phy))
	    # Randomise the species functional trait values
	    funNULL <- lapply(1:nrep, function(x) taxaShuffle(fun)) 	 
	 	# Do not change the site-species matrix
	    spSiteNULL <- list(samp=spSite, inv=invad)
	    
	    # Loop for each invader
	    AllInv <- sapply(invad, function(INV){	
		    # Calculate the null distribution of invaders indices	
		    tmp <- lapply(1:nrep, function(x) { 
	    			div.param.invasion.obs(spSite, phyNULL[[x]], funNULL[[x]], INV) } )  								 
	    	# get output for null model observations
		    obs_1 <- obs[[INV]]
	    		tmpp <- lapply(tmp, function(x) x[[1]])
	   		tmp1 <- lapply(1:nrow(obs_1), function(z) { 
	   					lapply(1:length(colnames(obs_1)), function(y) { 
	   						sapply(1:nrep, function(x) tmpp[[x]][z,y])} )} )
	   	 	meanNULL <- t(sapply(tmp1, function(z) sapply(z, mean, na.rm = TRUE)))
	   		sdNULL <- t(sapply(tmp1, function(z) sapply(z, sd, na.rm = TRUE)))
	   		rankNULL <- t(sapply(1:nrow(obs_1), function(z) {
	   						sapply(1:length(colnames(obs_1)), function(y) {
	    						ifelse(is.na(obs_1[z,y]), NA, rank(c(obs_1[z,y], tmp1[[z]][[y]]))[1] / (nrep+1)) }) }) )
		    zNULL <- (obs_1 - meanNULL) / sdNULL
		    colnames(zNULL) <- colnames(rankNULL) <- colnames(meanNULL) <- colnames(sdNULL) <- colnames(obs_1)
			output <- list(obs=obs_1, zNULL=zNULL, rankNULL=rankNULL, meanNULL=meanNULL, sdNULL=sdNULL)
			# check the number of simulated communities, and adjust the input format accordingly
			if(Ncom==1){ 
				output <- sapply(names(output), function(x){
					matrix( output[[x]][1,], nrow=1, dimnames=list(1, colnames( output[[x]] )) ) }, simplify = F, USE.NAMES=T) 
			}
		    return(output) 
		},  simplify = F, USE.NAMES=T )   
    }  
    if(is.null(null.model)){ 
 		AllInv <- sapply(invad, function(INV){
			if(Ncom>=1){ obs_1 <- obs[[INV]] }
			if(Ncom==1){ obs_1 <- matrix( obs[[INV]][1,], nrow=1, dimnames=list(1, colnames(obs[[1]])) ) }
		   	zNULL <- rankNULL <- meanNULL <- sdNULL <- NULL 
		    return(list(obs=obs_1, zNULL=zNULL, rankNULL=rankNULL, meanNULL=meanNULL, sdNULL=sdNULL) )
 		},  simplify = F, USE.NAMES=T ) 
 	}	
    return(AllInv)
}
	    
	