div.param.native.obs <- function(spSite, phy, fun, niche.optima, tree, indX.nat) {
    spSite <- as.data.frame(spSite)
    if (ncol(spSite) == 1) {
        spSite <- as.data.frame(t(spSite))
    }
    spSite.dummy <- replace(spSite, spSite > 0, 1)
    spIn <- colnames(spSite)[colSums(spSite) != 0]
    fun.NoZero <- fun[spIn, spIn]
    phy.NoZero <- phy[spIn, spIn]
    spSite.NoZero <- spSite[, spIn]
    if ("FD_pa_FEve" %in% indX.nat | "FD_pa_FDis" %in% indX.nat) {
        FD_pa_vill = dbFD(fun.NoZero, spSite.NoZero, w.abun = FALSE, calc.FRic = FALSE, calc.CWM = FALSE, calc.FDiv = FALSE, messages = FALSE)
    }
    if ("FD_ab_FEve" %in% indX.nat | "FD_ab_FDis" %in% indX.nat) {
        FD_ab_vill = dbFD(fun.NoZero, spSite.NoZero, w.abun = FALSE, calc.FRic = FALSE, calc.CWM = FALSE, calc.FDiv = FALSE, messages = FALSE)
    }
    if ("PD_pa_FEve" %in% indX.nat | "PD_pa_FDis" %in% indX.nat) {
        PD_pa_vill = dbFD(phy.NoZero, spSite.NoZero, w.abun = FALSE, calc.FRic = FALSE, calc.CWM = FALSE, calc.FDiv = FALSE, messages = FALSE)
    }
    if ("PD_ab_FEve" %in% indX.nat | "PD_ab_FDis" %in% indX.nat) {
        PD_ab_vill = dbFD(phy.NoZero, spSite.NoZero, w.abun = FALSE, calc.FRic = FALSE, calc.CWM = FALSE, calc.FDiv = FALSE, messages = FALSE)
    }
    
    results <- list()
    if ("TD_pa_simpson" %in% indX.nat) {
        results$TD_pa_simpson = rowSums(spSite.dummy)
    }
    if ("TD_pa_shannon" %in% indX.nat) {
        results$TD_pa_shannon = rowSums(spSite.dummy)
    }
    if ("TD_ab_simpson" %in% indX.nat) {
        results$TD_ab_simpson = diversity(spSite, index = "simpson")
    }
    if ("TD_ab_shannon" %in% indX.nat) {
        results$TD_ab_shannon = diversity(spSite, index = "shannon")
    }
    
    if ("FD_pa_mpd" %in% indX.nat) {
        results$FD_pa_mpd = mpd(spSite, fun, abundance.weighted = FALSE)
    }
    if ("FD_pa_mntd" %in% indX.nat) {
        results$FD_pa_mntd = mntd(spSite, fun, abundance.weighted = FALSE)
    }
    if ("FD_pa_CWM" %in% indX.nat) {
        results$FD_pa_CWM = apply(spSite.dummy, 1, function(x) {
            weighted.mean(niche.optima, w = x)
        })
    }
    if ("FD_pa_CSD" %in% indX.nat) {
        results$FD_pa_CSD = apply(spSite.dummy, 1, function(x) {
            sqrt(wtd.var(niche.optima, weights = x))
        })
    }
    if ("FD_pa_FEve" %in% indX.nat) {
        results$FD_pa_FEve = FD_pa_vill$FEve
    }
    if ("FD_pa_FDis" %in% indX.nat) {
        results$FD_pa_FDis = FD_pa_vill$FDis
    }
    if ("FD_pa_faith" %in% indX.nat) {
        results$FD_pa_faith = treedive(spSite, hclust(as.dist(fun)))
    }
    
    if ("FD_ab_mpd" %in% indX.nat) {
        results$FD_ab_mpd = mpd(spSite, fun, abundance.weighted = TRUE)
    }
    if ("FD_ab_mntd" %in% indX.nat) {
        results$FD_ab_mntd = mntd(spSite, fun, abundance.weighted = TRUE)
    }
    if ("FD_ab_CWM" %in% indX.nat) {
        results$FD_ab_CWM = apply(spSite, 1, function(x) {
            weighted.mean(niche.optima, w = x)
        })
    }
    if ("FD_ab_CSD" %in% indX.nat) {
        results$FD_ab_CSD = apply(spSite, 1, function(x) {
            sqrt(wtd.var(niche.optima, weights = x))
        })
    }
    if ("FD_ab_FEve" %in% indX.nat) {
        results$FD_ab_FEve = FD_ab_vill$FEve
    }
    if ("FD_ab_FDis" %in% indX.nat) {
        results$FD_ab_FDis = FD_ab_vill$FDis
    }
    if ("PD_pa_mpd" %in% indX.nat) {
        results$PD_pa_mpd = mpd(spSite, phy, abundance.weighted = FALSE)
    }
    if ("PD_pa_mntd" %in% indX.nat) {
        results$PD_pa_mntd = mntd(spSite, phy, abundance.weighted = FALSE)
    }
    if ("PD_pa_faith" %in% indX.nat) {
        results$PD_pa_faith = apply(spSite, 1, function(x) {
            tree.red <- drop.tip(tree, as.character(tree$tip.label[x == 0]))
            sum(tree.red$edge.length)
        })
    }
    if ("PD_pa_betasplit" %in% indX.nat) {
        results$PD_pa_betasplit = apply(spSite, 1, function(x) {
            tree.red <- drop.tip(tree, as.character(tree$tip.label[x == 0]))
            maxlik.betasplit(as.treeshape(tree.red))$max_lik
        })
    }
    if ("PD_pa_FEve" %in% indX.nat) {
        results$PD_pa_FEve = PD_pa_vill$FEve
    }
    if ("PD_pa_FDis" %in% indX.nat) {
        results$PD_pa_FDis = PD_pa_vill$FDis
    }
    
    if ("PD_ab_mpd" %in% indX.nat) {
        results$PD_ab_mpd = mpd(spSite, phy, abundance.weighted = TRUE)
    }
    if ("PD_ab_mntd" %in% indX.nat) {
        results$PD_ab_mntd = mntd(spSite, phy, abundance.weighted = TRUE)
    }
    if ("PD_ab_FEve" %in% indX.nat) {
        results$PD_ab_FEve = PD_ab_vill$FEve
    }
    if ("PD_ab_FDis" %in% indX.nat) {
        results$PD_ab_FDis = PD_ab_vill$FDis
    }
    tt <- sapply(results, function(X) {
        !is.null(X)
    })
    results2 <- results[tt]
    return(results2)
}