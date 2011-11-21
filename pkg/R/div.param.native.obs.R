div.param.native.obs <- function(spSite, phy, fun, niche.optima, tree, indX.nat){
	#  spSite=all.abundances2; niche.optima=niche.optima.nat; tree=tree.nat; phy=dist.phy.nat; fun=dist.fun.nat; nrep=n.rep.null.model; null.model = null.model; suit=sp.suit; taxa=taxa; indX.nat=indX.nat
	
    spSite <- as.data.frame(spSite)
    if (ncol(spSite)==1){spSite <- as.data.frame(t(spSite))}
    spSite.dummy <- replace(spSite, spSite > 0, 1)
    
    spSite.rao <- as.data.frame(t(spSite))
    spSite.rao.dummy <- replace(spSite.rao,spSite.rao > 0, 1)
    results <- list(
				if("TD_pa_simpson" %in% indX.nat) {TD_pa_simpson = diversity(spSite.dummy, index="simpson")},
				if("TD_pa_shannon" %in% indX.nat) {TD_pa_shannon = diversity(spSite.dummy, index="shannon")},
				if("TD_ab_simpson" %in% indX.nat) {TD_ab_simpson = diversity(spSite, index="simpson")},
				if("TD_ab_shannon" %in% indX.nat) {TD_ab_shannon = diversity(spSite, index="shannon")},
				
				if("FD_pa_mpd" %in% indX.nat) {FD_pa_mpd= mpd(spSite, fun, abundance.weighted=FALSE)},
				if("FD_pa_mntd" %in% indX.nat) {FD_pa_mntd= mntd(spSite, fun, abundance.weighted=FALSE)},
        #FD_pa_mntd.med = mntd_median(spSite, fun),
				if("FD_pa_CWM" %in% indX.nat) {FD_pa_CWM = apply(spSite.dummy, 1, function(x){weighted.mean(niche.optima, w=x)})},
				if("FD_pa_CSD" %in% indX.nat) {FD_pa_CSD = apply(spSite.dummy, 1, function(x){sqrt(wtd.var(niche.optima, weights=x))})},
				
				if("FD_ab_mpd" %in% indX.nat) {FD_ab_mpd = mpd(spSite, fun, abundance.weighted=TRUE)},
				if("FD_ab_mntd" %in% indX.nat) {FD_ab_mntd = mntd(spSite, fun, abundance.weighted=TRUE)},
				if("FD_ab_CWM" %in% indX.nat) {FD_ab_CWM = apply(spSite, 1, function(x){weighted.mean(niche.optima,w=x)})},
				if("FD_ab_CSD" %in% indX.nat) {FD_ab_CSD = apply(spSite, 1, function(x){sqrt(wtd.var(niche.optima,weights=x))})},
   		    	
				if("PD_pa_mpd" %in% indX.nat) {PD_pa_mpd = mpd(spSite, phy, abundance.weighted=FALSE)},
				if("PD_pa_mntd" %in% indX.nat) {PD_pa_mntd = mntd(spSite, phy, abundance.weighted=FALSE)},
        #PD_pa_mntd.med = mntd_median(spSite, phy),
				if("PD_pa_faith" %in% indX.nat) {PD_pa_faith = apply(spSite,1,function(x){
									tree.red <- drop.tip(tree,as.character(tree$tip.label[x==0]))
									y <- sum(tree.red$edge.length)
									return(y)})},
				if("PD_pa_colless" %in% indX.nat) {PD_pa_colless = apply(spSite, 1, function(x){
            tree.red <- drop.tip(tree, as.character(tree$tip.label[x==0]))
            colless(as.treeshape(tree.red), norm="yule")})},
				
				if("PD_ab_mpd" %in% indX.nat) {PD_ab_mpd = mpd(spSite, phy, abundance.weighted=TRUE)},
				if("PD_ab_mntd" %in% indX.nat) {PD_ab_mntd = mntd(spSite, phy, abundance.weighted=TRUE)}
				)


		tt <- sapply(results,function(X){!is.null(X)})
		results2 <- results[tt]
		names(results2) <- indX.nat				
		return(results2)
}
