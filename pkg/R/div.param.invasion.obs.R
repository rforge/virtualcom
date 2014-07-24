div.param.invasion.obs <- function(ispSite, iphy, ifun, imyinva){
	# ispSite = all.abundances.invaded ; iphy=dist.phy ; ifun=dist.fun ; imyinva=invader.ID_Pres 
	
	# Check the number of simulated communites, and adjust the INtput format
	Ncom <- nrow(ispSite)
	if(Ncom==1) {ispSite <- rbind(ispSite, ispSite) ; row.names(ispSite) <- c(1,2) }	
	
	# Run the calculations	
     phyloRes <- MDNC_MDNN_Invasive_new(samp=ispSite, dis=iphy, inva=imyinva)
     funRes <- MDNC_MDNN_Invasive_new(samp=ispSite, dis=ifun, inva=imyinva)

     iobs <- sapply(imyinva, function(x) {
     			tm <- cbind(sapply(names(phyloRes), function(y) unlist(phyloRes[[y]][x])), 
     						sapply(names(funRes), function(y)  unlist(funRes[[y]][x])) ) 
              	rownames(tm) <- 1: nrow(ispSite) 
              	colnames(tm) <- c(paste("phy_", names(phyloRes), sep=""), 
              						paste("fun_", names(phyloRes), sep=""))
              	tm}, simplify=F, USE.NAMES=T)
	# Check the number of simulated communites, and adjust the OUTtput format
     if(Ncom==1){ for(i in 1:length(iobs)){ iobs[[1]] <- iobs[[1]][1,] }  }
     iobs
}

