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
    evol.model <- switch(parameters["evol.model"], 
      		"1" = "BM",
      		"2" = "deltaTree",
      		"3" = "kappaTree",
      		"4" = "ouTree") 
    evol.model.param <- parameters["evol.model.param"]
  	min.phyl.signal <- parameters["min.phyl.signal"] 			
  	
    # Null model characteristics    
    # null.model <- as.character(parameters$null.model)
    null.model <- switch(parameters["null.model"], 
      		"1" = "taxa.labels",
      		"2" = "constrained",
      		"3" = "richness",
      		"4" = "frequency",
      		"5" = "independentswap",
      		"6" = "trialswap")
    suitability <- ifelse(null.model=="constrained", parameters$suitability, FALSE) # TRUE or FALSE to the use of suitability
    taxa.level <- ifelse(null.model=="constrained", parameters$taxa.level, NA) # Either nothing, a numeric between 0 and 1 (height-type of cutting), an integer ("number of groups"-type of cutting)
    
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
    pool <- create.pool(n.species.pool, n.invader.pool, evol.model, min.phyl.signal, evol.model.param, n.rep.null.model)
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
  	
    #Preparation of data
    tree.nat<-drop.tip(pool$phylo,pool$phylo$tip.label[pool$func$Inv==1])
    niche.optima.nat<-niche.optima.nat[as.character(tree.nat$tip.label)]
    all.abundances2<-all.abundances[,as.character(tree.nat$tip.label)]
    
    # Get indices and null models for diversity:
 	  dist.phy.nat <- cophenetic(tree.nat)
 	  dist.fun.nat <- as.matrix(dist(niche.optima.nat,diag=TRUE,upper=TRUE))
    dist.phy <- cophenetic(pool$phylo)
 	  dist.fun<- as.matrix(dist(niche.optima,diag=TRUE,upper=TRUE))
    taxa<-sp.suit<-NULL
    if (suitability){sp.suit<-dnorm(x=niche.optima, mean=env, sd=niche.breadth)/dnorm(x=env, mean=env, sd=niche.breadth)}

    if(is.numeric(taxa.level)){
      taxa<-create.taxa(tree.nat,taxa.level)
      taxa<-taxa[[1]]#Maybe later, this will be corrected to allow the use of different null models or constraint to evaluate the same communities.
    }

    # Choose the indices you want to put in div.param.native:
	  indX.nat <-c("TD_pa_simpson", "TD_pa_shannon", "TD_ab_simpson", "TD_ab_shannon", "FD_pa_mpd", "FD_pa_mntd", "FD_pa_CWM", "FD_pa_CSD", 
        "FD_ab_mpd", "FD_ab_mntd", "FD_ab_CWM", "FD_ab_CSD", "PD_pa_mpd", "PD_pa_mntd", "PD_pa_faith", "PD_ab_mpd", "PD_ab_mntd", 
        "PD_pa_colless", "FD_pa_FEve", "FD_pa_FDis", "FD_pa_faith", "FD_ab_FEve", "FD_ab_FDis", "PD_pa_FEve", "PD_pa_FDis",
        "PD_ab_FEve", "PD_ab_FDis") # never put only 1 index (at least 2), it would change the output format!
	
    indices.nat <- div.param.native(spSite=all.abundances2, niche.opt=niche.optima.nat, tree=tree.nat, phy=dist.phy.nat, fun=dist.fun.nat, nrep=n.rep.null.model, null.model = null.model, suit=sp.suit, taxa=taxa, indX.nat=indX.nat) # zNULL = NaN when sdNULL=0				  
    
    # collect results
    output <- list()
    output$parameter <- parameters
    output$pool <- pool
    output$natives$all.abundances <- all.abundances
    output$natives$indices <- indices.nat

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
        
        # collect results
        output$invaders$all.abundances <- all.abundances.invaded
        if (sum(all.abundances.invaded[, invader.ID], na.rm=TRUE) !=0 ){ 
          indices.inv <- div.param.invasion(spSite=all.abundances.invaded, phy=dist.phy, fun=dist.fun, nrep=n.rep.null.model, inva=invader.ID, null.model="native_inv")
          # collect results
          output$invaders$indices <- indices.inv 
        }  
    }            
    return(output)
}   
