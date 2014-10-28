#' Runs simulation experiments with VirtualCom
#' 
#' This function creates (1) a species pool (see \code{\link{create.pool}}), (2) assembles communities (see \code{\link{tamaure}}) 
#' and (3) computes functional and phylogenetic diversity indices. 
#' After natural communities are assembled invasion can occur.
#' 
#' @param parameters The only argument is a named vector of input parameter values, the entries are given in details
#' @return 
#' Information on input parameters, species pool, native assembly structure 
#' and invader assembly structure are returned as a list:
#' \describe{
#' \item{parameter}{ a vector of the input parameters }
#' \item{pool}{ a list of different objects that define the species pool. "func" is a dataframe of trait values and invader identity. "phylo" is a phylo object. "invader" is a list of the invader's ID, their niche optima and their performances given the environment in the community }
#' \item{natives}{a list of the native commuity structure and diversity measures. "all.abundances" is a site-by-species matrix of abundances. "indices" is a list of matrices containing functional and phylohgenetic diversity information (see details) }
#' \item{invaders}{ a list of community structure and diversity measures after the invasion processes, see "natives" for the structure }
#' }
#' @details
#' \describe{
#' \item{n.species.pool}{ number of species in the species pool (including native and invasive species) }
#' \item{evol.model}{ choice of the trait evolution model (see also \code{\link{create.pool}}) }
#' \item{rescale}{ if TRUE the trait values will be rescaled between 0 and 100 }
#' \item{OUwie.sigma.sq}{ applies only if evol.model is OUwie.sim, a numeric vector giving the values of sigma^2 for each selective regime  }
#' \item{OUwie.theta}{ applies only if evol.model is OUwie.sim, a numeric vector giving the values of theta for each selective regime } 
#' \item{OUwie.alpha}{ applies only if evol.model is OUwie.sim, a numeric vector giving the values of alpha for each selective regime  }
#' \item{evol.model.param}{ parameterization of the trait evolution model (see also \code{\link{create.pool}}) }
#' \item{species.pool.abundance}{ vector of abundances in the species pool; if NA then all species are equally abundant }  
#' \item{min.phyl.signal}{ minimum value of phylogenetic signal in species pool (see also \code{\link{create.pool}}) }   
#' \item{years}{ number of simulated timesteps; If plots=TRUE equilibrium conditions can be checked by eye  }  
#' \item{n.communities}{ number of simulated communities}
#' \item{env}{ value of the environment in the simulated community (e.g. temperature) }
#' \item{K}{ value of carrying capacity, i.e. number of individuals in the community  }  
#' \item{niche.breadth}{ value of standard deviation of the Gaussian distribution that describes the niche traits (identical for all species)  } 
#' \item{beta.env}{ value of the strength of the environmental filter (see details) }
#' \item{beta.comp}{ value of the strength of the competition filter (see details) }
#' \item{beta.abun}{ value of the strength of the recruitment filter, i.e. the advantage of already being present in the community (see also \code{\link{tamaure}}) }   
#' \item{competition}{ choice between symmetric, asymmetric and hierarchical competition; in the two latter cases species with higher trait values put pressure on species with lower trait values (species with lower trait values do not influence species with higher trait values), under asymmetic competition, competitive strength depends on niche overlap (i.e. very different species do not compete), under hierarchical competition, niche overlap is unimportant}
#' \item{intra.sp.com}{ assigns the strength of intraspecific competition; the value should range between 0 (no intraspecific competition) and 1 (intraspecific competition always higher than interspecific competition) }
#' \item{invasion.time}{ number of time-steps to simulate invasion (after native community has established); if 0 than there is no invasion }   
#' \item{n.invader.pool}{ number of invaders in species pool  }                                                                                                                                 
#' \item{InvDistrib}{ The invader's distribution in the phylogney can be either clusterd ("1"), random ("2") or  overdispersed ("3)  }                                                                                                                                 
#' \item{n.rep.null.model}{ number of repetitions for null models to test the diversity indices describing community composition  }
#' \item{null.model}{ different null models can be choosen, see details  }   
#' }
#'
#' The following null models can be choosen to test for the structure of the native community: "1" = taxa.labels,"2" = richness, "3" = frequency, "4" = independentswap, "5" = trialswap. These models are explained in more detail in the function ses.mpd in the package picante. 
#' The output lists for natives and invaders consist of 5 matrices describing the diversity structure. 
#' They contain information on (1) observed diversity, and based on the chosen null models: (2) z-scores, (3) ranks, (4) mean of the null distribution and (5) standard deviation of the null distribution. Each matrix gives this information for all sites and for a set of diversity indices (mpd, mntd, ...).  
#' @seealso \code{\link{tamaure}} for the community assembly function
#' @examples
#'   # load pre-prepared parameter table
#'   data(simple_param)	
#'   # run the first line of the parameter table
#'   single.run <- simulation.experiment(simple_param[1,]) # run the first line
#'   str(single.run)
#'   
#'   # setting up a full experiment (this may take a few minutes) 
#'   wrapper <- function(a){
#'   library(VirtualCom)
#'   data(simple_param)
#'   print(a)
#'   return(try(simulation.experiment(simple_param[a,])))	
#'   }  
#'   
#'   # running the experiment either with lapply
#'   output <- lapply(1:nrow(simple_param), wrapper) 
#'   
#'   # # or running the experiment with sfLapply using snowfall to allow for parallel computing
#'   # require(snowfall)
#'   # sfInit(parallel=TRUE, cpus=11)
#'   # output <- sfLapply(1:nrow(simple_param), wrapper)
#'   # sfStop()
#'   
#'   # extract results from list
#'   result.table <- get.results(output=output, myvar="obs", invader="FALSE")
#'   result.table$process <- ifelse(result.table$beta.env==0 & result.table$beta.comp==0, "Random", ifelse(result.table$beta.env!=0 & result.table$beta.comp==0, "Env", ifelse(result.table$beta.env==0 & result.table$beta.comp!=0, "Comp", "Both") ))
#'   
#'   # plot the funtional diversity (mpd) of communities in dependence on assembly processes
#'   require(ggplot2)
#'   ggplot(data=result.table, aes(x=process, y=FD_ab_mpd)) + geom_boxplot(width=0.8) +  xlab("Assembly rules") +  ylab("Functional diversity")
#'   
#'   @keywords assembly invasion competition habitat filtering
#' @export


simulation.experiment <- function(parameters) {
    # ------------------------------------- get input parameters -------------------------------------
    my.names <- names(parameters)
    parameters <- as.numeric(parameters)
    names(parameters) <- my.names
    
    # general
    years <- parameters["years"]
    invasion.time <- parameters["invasion.time"]
    n.rep.null.model <- parameters["n.rep.null.model"]
    
    # species pool
    n.species.pool <- parameters["n.species.pool"]
    n.invader.pool <- parameters["n.invader.pool"]
    niche.breadth <- parameters["niche.breadth"]
    evol.model <- switch(parameters["evol.model"], `1` = "BM", `2` = "deltaTree", `3` = "kappaTree", `4` = "lambdaTree", `5` = "OUwie.sim")
    evol.model.param <- parameters["evol.model.param"]
    rescale <- parameters["rescale"]
    OUwie.sigma.sq <- parameters["OUwie.sigma.sq"]
    OUwie.theta <- parameters["theta"]
    OUwie.alpha <- parameters["alpha"]
    min.phyl.signal <- parameters["min.phyl.signal"]
    InvDistrib <- parameters["InvDistrib"]
    
    # Null model characteristics null.model <- as.character(parameters$null.model)
    null.model <- switch(as.character(parameters["null.model"]), `0` = NULL, `1` = "taxa.labels", `2` = "richness", `3` = "frequency", `4` = "independentswap", `5` = "trialswap")
    
    # communities
    n.communities <- parameters["n.communities"]
    K <- parameters["K"]
    env <- parameters["env"]
    beta.env <- parameters["beta.env"]
    beta.comp <- parameters["beta.comp"]
    beta.abun <- parameters["beta.abun"]
    species.pool.abundance <- parameters["species.pool.abundance"]
    competition <- switch(parameters["competition"], `1` = "symmetric", `2` = "asymmetric", `3` = "hierarchical")
    intra.sp.com <- parameters["intra.sp.com"]
    
    # ------------------------------------- Species pool: phylogenetic tree, trait values and invaders in the tree -------------------------------------
    pool <- create.pool(n.species.pool=n.species.pool, n.invader.pool=n.invader.pool, evol.model=evol.model, rescale=rescale, min.phyl.signal=min.phyl.signal, evol.model.param=evol.model.param, OUwie.sigma.sq=OUwie.sigma.sq, OUwie.theta=OUwie.theta, OUwie.alpha=OUwie.alpha, nrep=n.rep.null.model)
    niche.optima <- pool$func$niche_evol
    names(niche.optima) <- pool$func$SpeciesID
    
    

    # ------------------------------------- Native community assembly ------------------------------------- 
    if (n.invader.pool == 0) {
        niche.optima.nat <- niche.optima
        names(niche.optima.nat) <- names(niche.optima)
    }
    if (n.invader.pool == 1) {
        niche.optima.nat <- niche.optima[pool$func$RandInv == 0]
        names(niche.optima.nat) <- pool$func[pool$func$RandInv == 0, "SpeciesID"]
    }
    if (n.invader.pool > 1) {
        if (InvDistrib == 1) {
            niche.optima.nat <- niche.optima[pool$func$UnderInv == 0]
            names(niche.optima.nat) <- pool$func[pool$func$UnderInv == 0, "SpeciesID"]
        }
        if (InvDistrib == 2) {
            niche.optima.nat <- niche.optima[pool$func$RandInv == 0]
            names(niche.optima.nat) <- pool$func[pool$func$RandInv == 0, "SpeciesID"]
        }
        if (InvDistrib == 3) {
            niche.optima.nat <- niche.optima[pool$func$OverInv == 0]
            names(niche.optima.nat) <- pool$func[pool$func$OverInv == 0, "SpeciesID"]
        }
    }
    
    # Initialization: creates the siteXspecies matrix
    all.communities <- matrix(0, nrow = n.communities, ncol = K, dimnames = list(1:n.communities, 1:K))
    all.abundances <- matrix(0, nrow = n.communities, ncol = length(niche.optima.nat), dimnames = list(1:n.communities, names(niche.optima.nat)))
    
    for (i in 1:n.communities) {
        one.community <- tamaure(niche.breadth = niche.breadth, niche.optima = niche.optima.nat, env = env, beta.env = beta.env, beta.comp = beta.comp, beta.abun = beta.abun, 
            years = years, K = K, community.in = NA, plot = FALSE, species.pool.abundance = species.pool.abundance, competition=competition, intra.sp.com =intra.sp.com )
        all.communities[i, ] <- one.community$community
        all.abundances[i, ] <- one.community$abundances
    }
    all.abundances <- ifelse(is.na(all.abundances), 0, all.abundances)
    
    # Preparation of data
    if (n.invader.pool <= 1) {
        InvasiveID <- as.character(pool$invader$ID)
    } else {
        InvasiveID <- as.character(pool$invader$ID[[InvDistrib]])
    }
    tree.nat <- drop.tip(pool$phylo, InvasiveID)
    NativeID <- as.character(tree.nat$tip.label)
    niche.optima.nat <- niche.optima.nat[NativeID]
    if (n.communities > 1) {
        all.abundances2 <- all.abundances[, NativeID]
    } else {
        all.abundances2 <- matrix(all.abundances[, NativeID], nrow = 1, dimnames = list(1, NativeID))
    }
    
    # Get indices and null models for diversity:
    dist.phy.nat <- cophenetic(tree.nat)
    dist.fun.nat <- as.matrix(dist(niche.optima.nat, diag = TRUE, upper = TRUE))
    dist.phy <- cophenetic(pool$phylo)
    dist.fun <- as.matrix(dist(niche.optima, diag = TRUE, upper = TRUE))
    
    # Choose the indices you want to put in div.param.native:
    indX.nat <- c("TD_pa_simpson", "TD_pa_shannon", "TD_ab_simpson", "TD_ab_shannon", "FD_pa_mpd", "FD_pa_mntd", "FD_pa_CWM", "FD_pa_CSD", "FD_ab_mpd", "FD_ab_mntd", 
        "FD_ab_CWM", "FD_ab_CSD", "PD_pa_mpd", "PD_pa_mntd", "PD_ab_mpd", "PD_ab_mntd", "FD_pa_FEve", "FD_pa_FDis", "FD_pa_faith", "FD_ab_FEve", "FD_ab_FDis", "PD_pa_FEve", 
        "PD_pa_FDis", "PD_pa_faith", "PD_pa_colless", "PD_ab_FEve", "PD_ab_FDis", "PD_pa_betasplit")  # never put only 1 index (at least 2), it would change the output format!
    
    if(n.rep.null.model==0) {
    	    indices.nat <- NA
    	} else {
    		indices.nat <- div.param.native(spSite = all.abundances2, niche.opt = niche.optima.nat, tree = tree.nat, phy = dist.phy.nat, fun = dist.fun.nat, nrep = n.rep.null.model, 
        null.model = null.model, indX.nat = indX.nat)  # zNULL = NaN when sdNULL=0\t\t\t\t 
    } 
    
    # collect results
    output <- list()
    output$parameter <- parameters
    output$pool <- pool
    output$natives$all.abundances <- all.abundances
    output$natives$indices <- indices.nat
    
    # ------------------------------------- Invaded community assembly -------------------------------------
    if (invasion.time > 0) {
        all.communities.invaded <- matrix(0, nrow = n.communities, ncol = K, dimnames = list(1:n.communities, 1:K))
        all.abundances.invaded <- matrix(0, nrow = n.communities, ncol = length(niche.optima), dimnames = list(1:n.communities, names(niche.optima)))
        invasion.dynamics <- list()
        
        # Invasion assembly: with abiotic filters, biotic interactions and recruitment
        for (i in 1:n.communities) {
            one.community <- tamaure(community.in = all.communities[i, ], niche.optima = niche.optima, years = invasion.time, niche.breadth = niche.breadth, env = env, 
                beta.env = beta.env, beta.comp = beta.comp, beta.abun = beta.abun, K = K, plot = FALSE, species.pool.abundance = species.pool.abundance, competition=competition, intra.sp.com =intra.sp.com )
            all.communities.invaded[i, ] <- one.community$community
            all.abundances.invaded[i, ] <- one.community$abundances
            invasion.dynamics[[i]] <- one.community$communities.over.time[-1, ]
        }
        all.abundances.invaded <- ifelse(is.na(all.abundances.invaded), 0, all.abundances.invaded)
        
        # Get indices and null models for diversity:
        colnames(all.abundances.invaded) <- paste("sp", colnames(all.abundances.invaded), sep = "")
        AllSpeciesID <- colnames(all.abundances.invaded)
        colnames(dist.phy) <- rownames(dist.phy) <- paste("sp", colnames(dist.phy), sep = "")
        colnames(dist.fun) <- rownames(dist.fun) <- paste("sp", colnames(dist.fun), sep = "")
        dist.phy <- dist.phy[AllSpeciesID, AllSpeciesID]
        dist.fun <- dist.fun[AllSpeciesID, AllSpeciesID]
        
        # collect results
        output$invaders$all.abundances <- all.abundances.invaded
        
        # Calculate the indices for the invaded communities
        InvasiveIDsp <- paste("sp", InvasiveID, sep = "")
        if (n.invader.pool <= 1) {
            invader.ID_Pres <- vector()
            if (sum(all.abundances.invaded[, InvasiveIDsp]) > 0) {
                invader.ID_Pres <- InvasiveIDsp
            }
        } else {
            if (n.communities > 1) {
                invader.ID_Pres <- InvasiveIDsp[colSums(all.abundances.invaded[, InvasiveIDsp]) > 0]
            }
            if (n.communities == 1) {
                invader.ID_Pres <- InvasiveIDsp[all.abundances.invaded[, InvasiveIDsp] > 0]
            }
        }
        
        if (length(invader.ID_Pres) > 0) {
            null.model <- ifelse(is.null(null.model), NULL, "native_inv")
            if(n.rep.null.model==0) {
    	       indices.inv <- NA
    	    } else {
    		    indices.inv <- div.param.invasion(spSite = all.abundances.invaded, phy = dist.phy, fun = dist.fun, nrep = n.rep.null.model, invad = invader.ID_Pres, null.model = null.model)
            } 

            # collect results
            output$invaders$indices <- indices.inv
        } else {
            output$invaders$indices <- "Aliens have never been successful!"
        }
    }
    return(output)
}