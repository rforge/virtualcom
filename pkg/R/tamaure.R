
tamaure <- function(niche.breadth=5, niche.optima, env, beta.env=0, beta.comp=0, beta.abun=0, years=20, K=20, community.in=NA, species.pool.abundance=NA, plot=FALSE, ...){
  # to test: niche.breadth=5; niche.optima=niche.optima; beta.env=0; beta.comp=0; beta.abun=1; years=20; K=10; community.in=NA; plot=FALSE; env=100; species.pool.abundance=NA

  #-------
  # 1. getting input
  #-------
	species.count <- length(niche.optima)	# nb of species in the species pool
	if(is.na(species.pool.abundance))   species.pool.abundance <- rep(1, species.count) 
	
  #-------
  # 2. Creating the distance matrix from the niche overlap
  #-------
	species.niche.dist <- as.matrix(dist(niche.optima, diag=TRUE, upper=TRUE))   # Makes a matrix of species distances between their optima
    # Make a matrix for putting the results for biotic interactions	
	species.niche.overlap <- matrix(0, ncol=species.count, nrow=species.count)
	colnames(species.niche.overlap) <- rownames(species.niche.overlap) <- names(niche.optima)
 
  species.niche.overlap <- outer(niche.optima, niche.optima, 
                                function(x, y) {2*pnorm(-abs((x-y))/2, 
                                                        mean=0, sd=niche.breadth) })

  #-------
  # 3. The loop
  #-------
  out.p.comp <- matrix(NA, nrow=years, ncol=species.count)	# getting the infomation about the strength of competition after the loop
	out.community <- matrix(NA, nrow=years, ncol=K)	# getting the information about the species ID at each run after the loop
  if(plot == TRUE) mpd.now <- sapply(1:round(years/4), 
                                function(x) {community.null <- sample(names(niche.optima), K, replace=TRUE)
                                          mean(species.niche.dist[community.null,community.null])})
  
  	# make random community from all the species with a carrying capacity or take community.in
	if (is.na(community.in[1])) {
      community <- out.community[1, ] <- sample(names(niche.optima), K, replace=TRUE)
    } else {
		  community <- out.community[1, ] <- community.in		
	}  
	
  	# Calculate the standardized environmental preference (Penv) for all species (this is fix for all runs)	
	p.env <- dnorm(niche.optima, env, niche.breadth) / dnorm(env, env, niche.breadth) 

  for (year in 2:years) {
    abundance <- as.numeric(table(community)[names(niche.optima)])
    abundance <- ifelse(is.na(abundance), 0, abundance / max(abundance, na.rm = TRUE))

    for (i in 1:K) {
      p.comp <- 1 - colSums(species.niche.overlap[community,]) / K
  	  p.all <- exp(beta.env * log(p.env) + beta.comp * log(p.comp) + log(species.pool.abundance + beta.abun * abundance)) 
      out <- sample(seq(community), 1)	# Choose randomly one individual to exclude 
		  community[out] <- sample(names(niche.optima), 1, prob=p.all) # Replace this excluded species by lottery competition      
    }    
  	out.community[year,] <- community 
    
 	  if(plot == TRUE) {  
        n.rand <- round(years/4) + 1
        filters = cbind(comp=round(p.comp, 2), env=round(p.env, 2), abun=round(abundance, 2), all=round(p.all,2))
        mpd.now <- c(mpd.now, mean(species.niche.dist[community,community]))
        mypar <- par(mfcol=c(3,2))
           barplot(filters[,1], xlab="Species ID", ylab="p.comp"); barplot(filters[,2], xlab="Species ID", ylab="p.env"); barplot(filters[,3], xlab="Species ID", ylab="abundance")  
           barplot(filters[,4], xlab="Species ID", ylab="p.all"); barplot(sort(table(community), decreasing = TRUE), xlab="Species ID", ylab="Abundance"); plot((1:length(mpd.now) - n.rand), mpd.now, type="l", xlab="Time", ylab="mpd"); abline(v=0, col=2); abline(h=mean(mpd.now[1:n.rand]), col=3); if(year>2) abline(h=mean(mpd.now[(n.rand+1):length(mpd.now)]), col=3) 
        par(mypar)
    }
  
  }
        
  # prepare output
	abundance.vector <- as.numeric(table(out.community[years,])[names(niche.optima)])
	names(abundance.vector) <- names(niche.optima)
	return(list(community=out.community[years,], abundances=abundance.vector, all.communities=out.community))
}   

