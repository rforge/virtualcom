# shuffle species identities

invasion_randomization <- function(samp, inva){
	# samp=spSite ; inva=invad
	
	# identify the non successful natives
	new.samp <- samp[,sample(1:ncol(samp))] 
	colnames(new.samp) <- colnames(samp)
	new.samp
	return(list(samp=new.samp, inv=inva))
}
