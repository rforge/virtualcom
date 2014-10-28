#' Simulates a species pool
#' 
#' Simulates a species pool containing a phylogeny, species trait values and information on invasion status.
#' 
#' @details
#' The different evolutionary models (argument: evol.mod) are based on 
#' tree transformations. These are done using the function rescale in 
#' the geiger package and are coded as follows: "1" = BM (Brownian motion),
#' "2" = deltaTree, "3" = kappaTree, "4" = lambdaTree, "5" = OUwie.sim (Ornstein-Uhlenbeck process). 
#' The tree transformations can be parameterized using the argument 
#' evol.model.param.
#' 
#' @param n.species.pool number of species (natives and invaders)
#' @param n.invader.pool number of invaders
#' @param evol.model choice of evolutionary model which determines phylogenetic signal, see details
#' @param min.phyl.signal minimum level of phylogenetic signal accepetd, if min.phyl.signal=NA no minimum limit
#' @param rescale if TRUE the trait values will be rescaled between 0 and 100
#' @param evol.model.param applies only if evol.model is deltaTree, kappaTree or lambdaTree, see details
#' @OUwie.sigma.sq applies only if evol.model is OUwie.sim, a numeric vector giving the values of sigma^2 for each selective regime 
#' @OUwie.theta applies only if evol.model is OUwie.sim, a numeric vector giving the values of theta for each selective regime 
#' @OUwie.alpha applies only if evol.model is OUwie.sim, a numeric vector giving the values of alpha for each selective regime 
#' @param nrep number of repetition for phylogenetic signal and tree imbalance tests
#' 
#' @return
#'   List of objects: 
#'   \describe{
#'   \item{func}{ evolved functional traits for each species of the pool }
#'   \item{phylo}{ phylogenetic tree of the species pool  }
#'   \item{indices}{ results of phylogenetic signal and colless test for phylogenetic tree shape } 
#'   \item{invader}{ list of the invader's identities and trait values } 
#'   \item{parameters}{ the parameters used for creating the species pool } 
#'   }
#' 
#' @seealso
#' \code{\link{simulation.experiment}} for running simulation 
#' experiments with an integrated simulation of species pools and 
#' diversity analyses of final community structures.
#' 
#' @examples
#' pool <- create.pool(n.species.pool=500, n.invader.pool=1, evol.model="deltaTree", min.phyl.signal=NA, evol.model.param=0.1, nrep=1)
#' # plot tree
#' plot(pool$phy)
#' # calculate phylogenetic signal
#' phylosignal(pool$func$niche_evol, pool$phy, checkdata=FALSE)
#' @export

create.pool <- function(n.species.pool, n.invader.pool, evol.model, rescale=TRUE, min.phyl.signal, evol.model.param=NA, OUwie.sigma.sq=NA, OUwie.theta=NA, OUwie.alpha=NA, nrep = 499) {
    # Create the species pool
    pool <- trait.evolution(branchingRate = 0.1, Nleaves = n.species.pool, Ninv = n.invader.pool, which.evolution.model = evol.model, rescale=rescale, tree.shape.param = evol.model.param, OUwie.sigma.sq=OUwie.sigma.sq, OUwie.theta=OUwie.theta, OUwie.alpha=OUwie.alpha)
    niche.optima <- pool[[1]]$niche_evol
    names(niche.optima) <- pool[[1]]$SpeciesID
    
    # Redo if phylogenetic signal has a minimum value:
    if (!is.na(min.phyl.signal)) {
        while (phylosignal(niche.optima, pool[[2]])$K < min.phyl.signal) {
            pool <- trait.evolution(branchingRate = 0.1, Nleaves = n.species.pool, Ninv = n.invader.pool, which.evolution.model = evol.model, rescale=rescale, tree.shape.param = evol.model.param, OUwie.sigma.sq=OUwie.sigma.sq, OUwie.theta=OUwie.theta, OUwie.alpha=OUwie.alpha)

            niche.optima <- pool[[1]]$niche_evol
            names(niche.optima) <- pool[[1]]$SpeciesID
        }
    }
    names(pool) <- c("func", "phylo")
    
    # calculate phylogenetic signal and tree imbalance
    phy.sign <- phylosignal(niche.optima, pool[[2]], reps = nrep)
    pool$indices$Blom_K <- phy.sign$K
    pool$indices$p_PIC <- phy.sign$PIC.variance.P
    col.test <- colless.test.no.print(as.treeshape(pool[[2]]), alternative = "greater", n.mc = nrep)
    pool$indices$Colless <- colless(as.treeshape(pool[[2]]), norm = "yule")
    pool$indices$p_Colless <- col.test$p.value
    
    # define invaders
    if (n.invader.pool == 1) {
        pool$invader <- list(ID = pool[[1]][pool[[1]]$RandInv == 1, "SpeciesID"], opt = niche.optima[pool[[1]][pool[[1]]$RandInv == 1, "SpeciesID"]])
    }
    if (n.invader.pool > 1) {
        pool$invader <- list(ID = list(UnderInv = pool[[1]][pool[[1]]$UnderInv == 1, "SpeciesID"], RandInv = pool[[1]][pool[[1]]$RandInv == 1, "SpeciesID"], OverInv = pool[[1]][pool[[1]]$OverInv == 
            1, "SpeciesID"]), opt = list(UnderInv = niche.optima[pool[[1]][pool[[1]]$UnderInv == 1, "SpeciesID"]], RandInv = niche.optima[pool[[1]][pool[[1]]$RandInv == 
            1, "SpeciesID"]], OverInv = niche.optima[pool[[1]][pool[[1]]$OverInv == 1, "SpeciesID"]]), distToOpt = NA)
    }
    pool$parameters <- data.frame(n.species.pool = n.species.pool, n.invader.pool = n.invader.pool, evol.model = evol.model, min.phyl.signal = min.phyl.signal, evol.model.param = evol.model.param, 
        nrep = nrep)
    return(pool)
} 