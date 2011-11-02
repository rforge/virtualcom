simulation.experiment <- function(parameters, ...){
  	# to test: parameters=as.numeric(param[1,]); names(parameters) <- names(param[1,])
  	  	
    #-------------------------------------
    # get input parameters
    #-------------------------------------  
    my.names <- names(parameters)
  	parameters <- as.numeric(parameters) 
    names(parameters) <- my.names   
    # general
  	years <- parameters["years"]		
  	invasion.time <- parameters["invasion.time"]	
 	  n.rep.null.model <- parameters["n.rep.null.model"]			
    # species pool
  	n.species.pool <- parameters["n.species.pool"] 						
  	n.invader.pool <- parameters["n.invader.pool"] 			
  	niche.breadth <- parameters["niche.breadth"]     
 	  evol.model <- ifelse(parameters["evol.model"]==1, "BM", ifelse(parameters["evol.model"]==2, "OU", "delta")) 
  	min.phyl.signal <- parameters["min.phyl.signal"] 			
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
    pool <- create.pool(n.species.pool, n.invader.pool, evol.model, min.phyl.signal)
    niche.optima <- pool$func$niche_evol
    names(niche.optima) <- pool$func$SpeciesID
    
    #-------------------------------------
    # Native community assembly 
    #-------------------------------------
    # Get the optima for the native species only 
    niche.optima.nat <- niche.optima[pool$func$Inv==0]	
  	names(niche.optima.nat) <- pool$func[pool$func$Inv==0,"SpeciesID"]  
  	
  	# Initialization: creates the siteXspecies matrix
  	all.communities <- matrix(0, nrow=n.communities, ncol=K, dimnames=list(1:n.communities,1:K))    
  	all.abundances <- matrix(0, nrow=n.communities, ncol=length(niche.optima.nat), dimnames=list(1:n.communities,names(niche.optima.nat)))
  	
    for (i in 1:n.communities){
      one.community <- tamaure(niche.breadth=niche.breadth, niche.optima=niche.optima.nat, env=env, beta.env=beta.env, beta.comp=beta.comp, beta.abun=beta.abun, years=years, K=K, community.in=NA, plot=FALSE, species.pool.abundance=species.pool.abundance)
    	all.communities[i,] <- one.community$community 
      all.abundances[i,] <- one.community$abundances     
    }
    all.abundances <- ifelse(is.na(all.abundances), 0, all.abundances)
  	
    # Get indices and null models for diversity:
 	  dist.phy <- cophenetic(pool$phylo)
 	  dist.fun <- as.matrix(dist(niche.optima,diag=TRUE,upper=TRUE))
    indices.nat <- div.param.native(spSite=all.abundances, phy=dist.phy, fun=dist.fun, nrep=n.rep.null.model, null.model = "taxa.labels") # zNULL = NaN when sdNULL=0				  
    
    # collect results
    output <- list()
    output$parameter <- parameters
    output$pool <- pool
    output$natives <- list()
    output$native$all.abundances <- all.abundances
    output$native$indices <- indices.nat

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
        invader.ID <- paste("sp", pool$invader$ID, sep="")
        
       	colnames(all.abundances.invaded) <- paste("sp",colnames(all.abundances.invaded),sep="")
    	  colnames(dist.phy) <- rownames(dist.phy) <- paste("sp",colnames(dist.phy),sep="")
    	  colnames(dist.fun) <- rownames(dist.fun) <- paste("sp",colnames(dist.fun),sep="")
        dist.phy <- dist.phy[colnames(all.abundances.invaded),colnames(all.abundances.invaded)]
        dist.fun <- dist.fun[colnames(all.abundances.invaded),colnames(all.abundances.invaded)]
        
        if (sum(all.abundances.invaded[, invader.ID], na.rm=TRUE) !=0 ){ 
          indices.inv <- div.param.invasion(spSite=all.abundances.invaded, phy=dist.phy, fun=dist.fun, nrep=n.rep.null.model, inva=invader.ID, null.model="native_inv")
        }
         
        # collect results
        output$invaders <- list()
        output$invaders$all.abundances <- all.abundances.invaded
        output$invaders$indices <- indices.inv  
    }   
         
    return(output)
}   