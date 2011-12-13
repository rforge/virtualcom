create.pool.invasion <- function(paramPool,...){
    # n.species.pool=n.species.pool; n.invader.pool=n.invader.pool; evol.model=evol.model; env=env;min.phyl.signal=min.phyl.signal; evol.model.param=evol.model.param; n.rep.null.model=n.rep.null.model
  
      
    #-------------------------------------
    # get input parameters
    #-------------------------------------  
    my.names <- names(paramPool)
  	parameters <- as.numeric(paramPool) 
    names(parameters) <- my.names 
     
    n.rep.null.model <- parameters["n.rep.null.model"]
    n.species.pool <- parameters["n.species.pool"] 						
  	n.invader.pool <- parameters["n.invader.pool"]
    evol.model <- switch(parameters["evol.model"], 
      		"1" = "BM",
      		"2" = "deltaTree",
      		"3" = "kappaTree",
      		"4" = "ouTree") 
    evol.model.param <- parameters["evol.model.param"]
  	min.phyl.signal <- parameters["min.phyl.signal"] 			
  	env <- parameters["env"]	
  
    pool <- trait.evolution(branchingRate=0.1, Nleaves=n.species.pool, Ninv=n.invader.pool, which.evolution.model=evol.model, extraTreeParam=evol.model.param)

  	niche.optima <- pool[[1]]$niche_evol
  	names(niche.optima) <- pool[[1]]$SpeciesID

  	# Redo if phylogenetic signal should not be below a certain value:
  	if(!is.na(min.phyl.signal)){
  		## If K < myphySig recreate the tree and traits again until K >= myphySig:
  		while(phylosignal(niche.optima, pool[[2]])$K < min.phyl.signal){
   	    pool <- trait.evolution(branchingRate=0.1, Nleaves=n.species.pool, Ninv=n.invader.pool, which.evolution.model=evol.model)
 	      niche.optima <- pool[[1]]$niche_evol
        names(niche.optima) <- pool[[1]]$SpeciesID
	    }
  	}
    names(pool) <-  c("func", "phylo")
    
    # calculate phylogenetic signal and tree imbalance
    phy.sign <- phylosignal(niche.optima, pool[[2]], reps = nrep)
    pool$indices$K <- phy.sign$K
    pool$indices$p_PIC <- phy.sign$PIC.variance.P
     
    col.test <- colless.test.no.print(as.treeshape(pool[[2]]),alternative="greater", n.mc = nrep)
    pool$indices$Colless <- colless(as.treeshape(pool[[2]]), norm="yule")
    pool$indices$p_Colless <- col.test$p.value
  	
    if(n.invader.pool==1){
      pool$invader <- list(ID=pool[[1]][pool[[1]]$Inv==1,"SpeciesID"], opt=niche.optima[pool[[1]][pool[[1]]$Inv==1,"SpeciesID"]], distToOpt= abs(env - niche.optima[pool[[1]][pool[[1]]$Inv==1,"SpeciesID"]]))
    }
    
    if(n.invader.pool > 1){
      pool$invader <- list(ID=list(OverInv=pool[[1]][pool[[1]]$OverInv==1,"SpeciesID"],RandInv=pool[[1]][pool[[1]]$RandInv==1,"SpeciesID"],UnderInv=pool[[1]][pool[[1]]$UnderInv==1,"SpeciesID"]),
                                  opt=list(OverInv=niche.optima[pool[[1]][pool[[1]]$OverInv==1,"SpeciesID"]],RandInv=niche.optima[pool[[1]][pool[[1]]$RandInv==1,"SpeciesID"]],UnderInv=niche.optima[pool[[1]][pool[[1]]$UnderInv==1,"SpeciesID"]]), 
                                  distToOpt= NA)
    }
    pool$parameters <- parameters
    
    return(pool)
}
