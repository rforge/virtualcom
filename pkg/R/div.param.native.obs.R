div.param.native.obs <- function(spSite, phy, fun,niche.opt,tree){
    spSite.dummy<-replace(spSite,spSite>0,1)
    return(list(PD_pa_mpd = mpd(spSite, phy, abundance.weighted=FALSE),
    PD_pa_mntd = mntd(spSite, phy, abundance.weighted=FALSE),
    PD_pa_mntd.med = mntd_median(spSite, phy),
    FD_pa_mpd= mpd(spSite, fun, abundance.weighted=FALSE),
    FD_pa_mntd= mntd(spSite, fun, abundance.weighted=FALSE),
    FD_pa_mntd = mntd_median(spSite, fun),
    PD_ab_mpd = mpd(spSite, phy, abundance.weighted=TRUE),
    PD_ab_mntd = mntd(spSite, phy, abundance.weighted=TRUE),
    PD_rao.a = disc(as.data.frame(t(spSite)), as.dist(phy)),
    FD_ab_mpd = mpd(spSite, fun, abundance.weighted=TRUE),
    FD_ab_mntd = mntd(spSite, fun, abundance.weighted=TRUE),
    FD_ab_rao.a = disc(as.data.frame(t(spSite)), as.dist(fun)),
    TD_ab_simpson = diversity(spSite,index="simpson"),
    TD_ab_shannon = diversity(spSite,index="shannon"),
    TD_pa_simpson = diversity(spSite.dummy,index="simpson"),
    TD_pa_shannon = diversity(spSite.dummy,index="shannon"),
    FD_ab_CWM = apply(spSite,1,function(x){weighted.mean(niche.optima,w=x)}),
    FD_pa_CWM = apply(spSite.dummy,1,function(x){weighted.mean(niche.optima,w=x)}),
    ))
}
