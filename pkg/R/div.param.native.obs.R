div.param.native.obs <- function(spSite, phy, fun, niche.optima, tree){
    spSite <- as.data.frame(spSite)
    if (ncol(spSite)==1){spSite <- as.data.frame(t(spSite))}
    spSite.dummy <- replace(spSite, spSite > 0, 1)
    
    spSite.rao <- as.data.frame(t(spSite))
    spSite.rao.dummy <- replace(spSite.rao,spSite.rao > 0, 1)
    return(list(
				TD_pa_simpson = diversity(spSite.dummy, index="simpson"),
				TD_pa_shannon = diversity(spSite.dummy, index="shannon"),
				TD_ab_simpson = diversity(spSite, index="simpson"),
				TD_ab_shannon = diversity(spSite, index="shannon"),
				
				FD_pa_mpd= mpd(spSite, fun, abundance.weighted=FALSE),
				FD_pa_mntd= mntd(spSite, fun, abundance.weighted=FALSE),
        #FD_pa_mntd.med = mntd_median(spSite, fun),
				FD_pa_rao = t(divc(spSite.rao.dummy, as.dist(fun))),
				FD_pa_CWM = apply(spSite.dummy, 1, function(x){weighted.mean(niche.optima, w=x)}),
				FD_pa_CSD = apply(spSite.dummy, 1, function(x){sqrt(wtd.var(niche.optima, weights=x))}),
				
				FD_ab_mpd = mpd(spSite, fun, abundance.weighted=TRUE),
				FD_ab_mntd = mntd(spSite, fun, abundance.weighted=TRUE),
				FD_ab_rao = t(divc(spSite.rao, as.dist(fun))),
				FD_ab_CWM = apply(spSite, 1, function(x){weighted.mean(niche.optima,w=x)}),
				FD_ab_CSD = apply(spSite, 1, function(x){sqrt(wtd.var(niche.optima,weights=x))}),
   		    	
				PD_pa_mpd = mpd(spSite, phy, abundance.weighted=FALSE),
				PD_pa_mntd = mntd(spSite, phy, abundance.weighted=FALSE),
        #PD_pa_mntd.med = mntd_median(spSite, phy),
				PD_pa_rao = t(divc(spSite.rao.dummy, as.dist(phy))),
				PD_pa_faith = apply(spSite,1,function(x){
									tree.red <- drop.tip(tree,as.character(tree$tip.label[x==0]))
									y <- sum(tree.red$edge.length)
									return(y)}),
				PD_pa_colless = apply(spSite, 1, function(x){
            tree.red <- drop.tip(tree, as.character(tree$tip.label[x==0]))
            colless(as.treeshape(tree.red), norm="yule")}),
				
				PD_ab_mpd = mpd(spSite, phy, abundance.weighted=TRUE),
				PD_ab_mntd = mntd(spSite, phy, abundance.weighted=TRUE),
				PD_ab_rao = t(divc(spSite.rao, as.dist(phy)))
				
				))
}
