tamaure <- function(niche.breadth=5, niche.optima, env, beta.env=0, beta.comp=0, beta.abun=0, years=20, K=20, community.in=NA, species.pool.abundance=NA, plot=FALSE, ...){
  # to test: niche.breadth=5; niche.optima=niche.optima; beta.env=0; beta.comp=0; beta.abun=1; years=20; K=10; community.in=NA; plot=TRUE; env=100; species.pool.abundance=species.pool.abundance=NA

  #-------
  # 1. getting input
  #-------
	species.count <- length(niche.optima)	# nb of species in the species pool
	if(is.na(species.pool.abundance[1]))   species.pool.abundance <- rep(1, species.count)

  #-------
  # 2. Creating the distance matrix from the niche overlap
  #-------
	species.niche.dist <- as.matrix(dist(niche.optima, diag=TRUE, upper=TRUE)) # to calculate mpd
  species.niche.overlap <- outer(niche.optima, niche.optima,
                                function(x, y) {2*pnorm(-abs((x-y))/2,
                                                        mean=0, sd=niche.breadth) })
  #-------
  # 3. The loop
  #-------
	out.community <- as.data.frame(matrix(NA, nrow=years+1, ncol=K)) # to store community composition over time
  if(plot == TRUE) mpd.now <- sapply(1:round(years/4),
                                function(x) {community.null <- sample(names(niche.optima), K, replace=TRUE)
                                          mean(species.niche.dist[community.null,community.null])})

 	# initialization: random or community.in
	if (is.na(community.in[1])) {
      community <- out.community[1, ] <- sample(names(niche.optima), K, replace=TRUE)
    } else {
		  community <- out.community[1, ] <- community.in
	}

	# environmental filter (fix through time)
	log.p.env <- dnorm(x=niche.optima, mean=env, sd=niche.breadth, log=TRUE) - dnorm(x=env, mean=env, sd=niche.breadth, log=TRUE)
  # loop - years
  for (year in 1:years) {
    # loop - asynchroneous updating within years
    for (i in 1:K) {
      abundance <- as.numeric(table(community)[names(niche.optima)])
      abundance <- ifelse(is.na(abundance), 0, abundance)
      p.comp <- 1 - colSums(species.niche.overlap[community,]) / K
  	  p.all <- exp(beta.env * log.p.env + beta.comp * log(p.comp) + log(species.pool.abundance + beta.abun * abundance))
      p.all <- ifelse(is.na(p.all), min(p.all, na.rm=TRUE), p.all)
      if(sd(p.all, na.rm=TRUE) == 0) p.all=NULL
      out <- sample(seq(community), 1)	# Choose randomly one individual to exclude
		  community[out] <- sample(names(niche.optima), 1, prob=p.all) # Replace this excluded species by lottery competition
    }
    out.community[year + 1,] <- community
  #-------
  # 4. Plot (if requested)
  #-------
    if(plot == TRUE) {
      n.rand <- round(years/4) + 1
      if(is.null(p.all))  p.all <- rep(0, species.count)
      filters = cbind(comp=round(p.comp, 2), env=round(exp(log.p.env), 2), abun=round(abundance, 2), all=round(p.all, 2))
      mpd.now <- c(mpd.now, mean(species.niche.dist[community,community]))
      mypar <- par(mfcol=c(3,2))
         barplot(filters[,1], xlab="Species ID", ylab="p.comp"); barplot(filters[,2], xlab="Species ID", ylab="p.env"); barplot(filters[,3], xlab="Species ID", ylab="abundance")
         barplot(filters[,4], xlab="Species ID", ylab="p.all"); hist(table(community), xlab="Abundance", ylab="Frequency", main=""); plot((1:length(mpd.now) - n.rand), mpd.now, type="l", xlab="Time", ylab="mpd"); abline(v=0, col=2); abline(h=mean(mpd.now[1:n.rand]), col=3); if(year>1) abline(h=mean(mpd.now[(n.rand+1):length(mpd.now)]), col=3)
      par(mypar)
    }
  }

  # prepare output
	abundance.vector <- as.numeric(table(community)[names(niche.optima)])
	abundance.vector <- ifelse(is.na(abundance.vector), 0, abundance.vector)
	names(abundance.vector) <- names(niche.optima)
	return(list(community=community, abundances=abundance.vector, communities.over.time=out.community))
}
   

