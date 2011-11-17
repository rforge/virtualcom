create.pool <- function(n.species.pool, n.invader.pool, evol.model, min.phyl.signal){
   	pool <- trait.evolution(branchingRate=0.1, Nleaves=n.species.pool, Ninv=n.invader.pool, which.evolution.model=evol.model, mytheta=1)
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
    phy.sign<-phylosignal(niche.optima, pool[[2]])
    pool$K$value<-phy.sign$K
    pool$K$pv<-phy.sign$PIC.variance.P
     
    col.test<-colless.test(as.treeshape(pool[[2]]),alternative="greater")
    pool$colless.rank<-col.test$p.value
  	names(pool) <-  c("func", "phylo")
    # pool$K <- phylosignal(niche.optima, pool[[2]])$K
    pool$invader <- list(ID=pool[[1]][pool[[1]]$Inv==1,"SpeciesID"], opt=niche.optima[pool[[1]][pool[[1]]$Inv==1,"SpeciesID"]], distToOpt= abs(env - niche.optima[pool[[1]][pool[[1]]$Inv==1,"SpeciesID"]]))
    return(pool)
}
