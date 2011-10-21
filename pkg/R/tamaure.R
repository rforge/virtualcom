tamaure <- function(niche.breath=5, opt, env=21, cstE=1, cstB=-0.1, cstSR = 1, years=10, kcapa=100, communityIn=NA, plot=FALSE, mypar=c(1,1), comp="NO", NN=NA, ...){
  #	env <-  21      	 # envir value
  #	cstE <- 1			 # const for environment
  #	cstB <- -0.1  		 # cst for the biotic term
  #	years <- 100			 # number of stabilization runs
  #	niche.breath <- 5	# niche breath
  #	kcapa <- 100		# Community carring capacity
  #	communityIn <- NA # do we start with a given community or with a randomly assembled one?
  # comp  # 1 is niche overlap competition and 2 is nearest neighbour competition
  # to test: niche.breath=myniche.breath; opt=myopt; env=50; cstE=mycstE; cstB=mycstB; years=timesteps; kcapa=mycapa; communityIn=NA; plot=FALSE; mypar=c(2,5); comp = "NO"; NN=myNN
	
  #-------
  # 1. getting input
  #-------
	nbspc <- length(opt)	# nb of species in the species pool
	stdN <- rep(niche.breath,nbspc) 			# std of the niche for all species (identical for all species)
	run <- kcapa * years # how often to repeat the replacement of individuals

  #-------
  # 2. Creating the distance matrix from the niche overlap
  #-------
	ourdist <- as.matrix(dist(opt,diag=TRUE,upper=TRUE))   # Makes a matrix of species distances between their optima
    # Make a matrix for putting the results for biotic interactions	
	didi <- matrix(0,ncol=nbspc,nrow=nbspc)
	colnames(didi) <- rownames(didi) <- names(opt)
 	if(!is.na(NN)) {
    	nearestN <- t(apply(ourdist, 1, function(x) order(x)[1:NN]))
    	furthestN <- t(apply(ourdist, 1, function(x) order(x)[NN:length(x)]))
    }			       
	for (i in 1:nbspc){
		# biotic interaction due to niche overlap
		if(comp!="NN"){ 	#CHANGED TAMI 22/08 (I just changed the order here)
   			for (j in 1:nbspc) didi[j,i] <- pnorm(opt[j] - ourdist[j,i]/2, opt[j], unlist(stdN[j]))*2 # pnorm = the smallest area (probability) of the a normal distribution at a certain value
    	}
    	# biotic interaction only with ourdist nearest neighbours
    	if(comp=="NN") didi[nearestN[i,],i] <- 1 
    	# biotic interaction due to niche overlap but only with ourdist nearest neighbours
        if(comp=="NO_cut") didi[furthestN[i,],i] <- 0 
	}
	# diag(didi) <- 0 #need to decide how much competition with individuals of own species

  #-------
  # 3. The loop
  #-------
  	outPbio <- list()	# getting the infomation about the strength of competition after the loop
	outspcC <- list()	# getting the information about the species ID at each run after the loop
  
  	# make random community from all the species with a carrying capacity or take communityIn
	if(is.na(communityIn[1])) {spcC <- outspcC[[1]] <- sample(names(opt),kcapa,replace=TRUE)} else {
		spcC <- outspcC[[1]] <- communityIn		
	}  
	
  	# Calculate the standardized environmental preference (Penv) for all species (this is fix for all runs)	
	Penv <- (dnorm(opt, env, stdN) - mean(dnorm(opt, env, stdN)))/sd(dnorm(opt, env, stdN))  

	j=2					# To start recording the loop outputs after the starting values (recorded in the first position of the list)
	for(i in 1:run){  
	  	# Calculating a standardized competition stength (Pbio)
		if(sd(colSums(didi[spcC,]))!= 0) Pbio <- colSums(didi[spcC,]) - mean(colSums(didi[spcC,]))/sd(colSums(didi[spcC,])) else Pbio <-  rep(1, length(opt)) 
    	# Calculate the overall likelyhood to enter for each species as a weighted sum of Penv and Pbio
		
		if(i %% kcapa == 0 | i == 1) {
			frQ <- as.matrix(table(spcC))   
			Abu <- rep(0,nbspc)
			names(Abu) <- names(opt)
			Abu[rownames(frQ)] <- frQ 			 # To do: Maybe make a baseline for the seed rain from outside the community?
			Abu <- (Abu - mean(Abu))/sd(Abu)
		}
		
		Pcomb <- transf(cstE * Penv + cstB * Pbio + cstSR * Abu) # transf shifts all values to zero and above
		out <- sample(1:length(spcC),1)	# Choose randomly one individual to exclude 
		spcC[out] <- sample(names(opt),1,prob=Pcomb) # Replace this excluded species by lottery competition	   
		# this is for recording information on development of communities and Pbio over time
    	if(i %% kcapa == 0){ 						    	# Saving the informations (%% is "modulo"): for instance this records the data every run number dividible by kcapa...
			outspcC[[j]] = spcC
			outPbio[[j]] = Pbio
			j = j+1
		}
	}
	
	# this is for plotting the community structure
	if (plot==TRUE){
		my_par <- par(mfrow=mypar)
   		for(i in 1:(years+1)){		
			hist(opt[outspcC[[i]]],xlim=c(0,100), main=paste("year =", i-1), xlab="Environmental optima",breaks=nbspc)
			abline(v=env,col="red")
    	}
    	par(my_par)
	} 
	
	# prepare output
	tt <- table(outspcC[[years+1]])
	vv <- as.vector(tt)
	names(vv) <- names(tt)

	#return(list(communities=outspcC, bio=outPbio, abio=Penv, spcXs=t(as.matrix(vv))))
	return(list(communities=outspcC[[years+1]], spcXs=t(as.matrix(vv))))
}   
