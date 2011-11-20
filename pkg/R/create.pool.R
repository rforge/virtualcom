create.pool <- function(n.species.pool, n.invader.pool, evol.model, min.phyl.signal, evol.model.param, nrep=499){
   	pool <- trait.evolution(branchingRate=0.1, Nleaves=n.species.pool, Ninv=n.invader.pool, which.evolution.model=evol.model, extraTreeParam=evol.model.param)
   	niche.optima <- pool[[1]]$niche_evol
  	names(niche.optima) <- pool[[1]]$SpeciesID

  	# Redo if phylogenetic signal should not be below a certain value:
  	if(!is.na(min.phyl.signal)){
  		## If K < myphySig recreate the tree and traits again until K >= myphySig:
  		while(phylosignal(niche.optima, pool[[2]])$K < min.phyl.signal){
   	    pool <- trait.evolution(branchingRate=0.1, Nleaves=n.species.pool, Ninv=n.invader.pool, which.evolution.model=evol.model, mytheta=1)
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
  	
    pool$invader <- list(ID=pool[[1]][pool[[1]]$Inv==1,"SpeciesID"], opt=niche.optima[pool[[1]][pool[[1]]$Inv==1,"SpeciesID"]], distToOpt= abs(env - niche.optima[pool[[1]][pool[[1]]$Inv==1,"SpeciesID"]]))
    return(pool)
}
