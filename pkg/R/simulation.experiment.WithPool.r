simulation.experiment.WithPool <- function(parameters, ...){
  	# to test: parameters=as.numeric(param.com[1,]); names(parameters) <- names(param.com[1,]); 
  
    #-------------------------------------
    # get input parameters
    #-------------------------------------  
    my.names <- names(parameters)
  	parameters <- as.numeric(parameters) 
    names(parameters) <- my.names 
     
    # general
  	years <- parameters["years"]		
  	invasion.time <- parameters["invasion.time"]	
  	niche.breadth <- parameters["niche.breadth"]
    InvDistrib <- parameters["InvDistrib"]
    
    # communities
  	n.communities <- parameters["n.communities"]			
  	K <- parameters["K"] 			
  	env <- parameters["env"]	
  	beta.env <- parameters["beta.env"]			
  	beta.comp <- parameters["beta.comp"]				
  	beta.abun <- parameters["beta.abun"]
    species.pool.abundance <- parameters["species.pool.abundance"]

  
    #-------------------------------------
    # Species pool: phylogenetic tree, trait values and invaders in the tree
    #-------------------------------------     
    pool <- pools[[parameters["poolID"]]]
    
    niche.optima <- pool$func$niche_evol
    names(niche.optima) <- pool$func$SpeciesID
    
    #-------------------------------------
    # Native community assembly 
    #-------------------------------------
    # Get the optima for the native species only 
      # depends on the number of invaders
    if(pool$parameters["n.invader.pool"] <=1){
      niche.optima.nat <- niche.optima[pool$func$Inv==0]	
  	  names(niche.optima.nat) <- pool$func[pool$func$Inv==0,"SpeciesID"]  
    }
    
    if(pool$parameters["n.invader.pool"] >1){
      if(InvDistrib == 1){
        niche.optima.nat <- niche.optima[pool$func$UnderInv==0]  
  	    names(niche.optima.nat) <- pool$func[pool$func$UnderInv==0,"SpeciesID"] 
      }
       if(InvDistrib == 2){
        niche.optima.nat <- niche.optima[pool$func$RandInv==0]  
        names(niche.optima.nat) <- pool$func[pool$func$RandInv==0,"SpeciesID"] 
      }
       if(InvDistrib == 3){
        niche.optima.nat <- niche.optima[pool$func$OverInv==0]  
        names(niche.optima.nat) <- pool$func[pool$func$OverInv==0,"SpeciesID"] 
      }
    }
    
    
  	# Initialization: creates the siteXspecies matrix
  	all.communities <- matrix(0, nrow=n.communities, ncol=K, dimnames=list(1:n.communities,1:K))    
  	all.abundances <- matrix(0, nrow=n.communities, ncol=length(niche.optima.nat), dimnames=list(1:n.communities,names(niche.optima.nat)))
  	
    for (i in 1:n.communities){
      one.community <- tamaure(niche.breadth=niche.breadth, niche.optima=niche.optima.nat, env=env, beta.env=beta.env, beta.comp=beta.comp, beta.abun=beta.abun, years=years, K=K, community.in=NA, plot=FALSE, species.pool.abundance=species.pool.abundance)
    	all.communities[i,] <- one.community$community 
      all.abundances[i,] <- one.community$abundances     
    }
    all.abundances <- ifelse(is.na(all.abundances), 0, all.abundances)
  	
    #Preparation of data
   
    # collect results
    output <- list()
    output$parameter <- parameters
    output$pool <- pool
    output$natives$all.abundances <- all.abundances

    #-------------------------------------
    # Invaded community assembly
    #-------------------------------------
    if(invasion.time > 0) {
       	all.communities.invaded <- matrix(0, nrow=n.communities, ncol=K, dimnames=list(1:n.communities,1:K))    
      	all.abundances.invaded <- matrix(0, nrow=n.communities, ncol=length(niche.optima), dimnames=list(1:n.communities,names(niche.optima)))
        invasion.dynamics <- list()
        
      	# Invasion assembly: with abiotic filters and biotic interactions and reproduction/seed rain	
        for (i in 1:n.communities){
          one.community <- tamaure(community.in=all.communities[i,], niche.optima=niche.optima, years=invasion.time, niche.breadth=niche.breadth, env=env, beta.env=beta.env, beta.comp=beta.comp, beta.abun=beta.abun, K=K, plot=FALSE, species.pool.abundance=species.pool.abundance)
        	all.communities.invaded[i,] <- one.community$community 
          all.abundances.invaded[i,] <- one.community$abundances
          invasion.dynamics[[i]] <- one.community$communities.over.time[-1,]
        }  
        all.abundances.invaded <- ifelse(is.na(all.abundances.invaded), 0, all.abundances.invaded)
       
        # Get indices and null models for diversity:
        #invader.ID <- paste("sp", pool$invader$ID, sep="")
      
        # collect results
        output$invaders$all.abundances <- all.abundances.invaded
        }  
    return(output)
}   
