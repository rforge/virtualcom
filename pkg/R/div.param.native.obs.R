div.param.native.obs <- function(spSite, phy, fun){
    return(list(mpd_phy = mpd(spSite, phy, abundance.weighted=FALSE),
    mntd_phy = mntd(spSite, phy, abundance.weighted=FALSE),
    mntd_phy_med = mntd_median(spSite, phy),
    mpd_fun = mpd(spSite, fun, abundance.weighted=FALSE),
    mntd_fun = mntd(spSite, fun, abundance.weighted=FALSE),
    mntd_fun_med = mntd_median(spSite, fun),
    mpd_phy_ab = mpd(spSite, phy, abundance.weighted=TRUE),
    mntd_phy_ab = mntd(spSite, phy, abundance.weighted=TRUE),
    mpd_fun_ab = mpd(spSite, fun, abundance.weighted=TRUE),
    mntd_fun_ab = mntd(spSite, fun, abundance.weighted=TRUE)))
}
