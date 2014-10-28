#' Community assembly from a given species pool
#' 
#' This function simulates the spatially implicit assembly of individuals 
#' into a community over time. Individuals come from a given pool of 
#' species. In the pool each species is described by its trait's niche 
#' (mean and sd of a Gaussian distribution) and its frequency of 
#' occurrence in the studied region. 
#' 
#' @param niche.breadth value of standard deviation of the Gaussian distributions that describe the niches (identical for all species)
#' @param niche.optima vector of niche optima (means of the Gaussian distributions that describe the niches) in species-pool
#' @param env value of the environment in the simulated community (e.g. temperature)
#' @param beta.env value of the strength of the environmental filter (see details)
#' @param beta.comp value of the strength of the competition filter (see details)
#' @param beta.abun value of the strength of the recruitment filter (see details)
#' @param years number of simulated time-steps; if plots == TRUE equilibrium conditions can be checked by eye
#' @param K value of carrying capacity, i.e. number of individuals in the community
#' @param community.in vector of individuals (described by species names) already in the community; if NA community is initialized randomly from species pool
#' @param species.pool.abundance vector of species frequency of occurrence in species pool; if NA than all are equally frequent
#' @param plot if TRUE then community composition is plotted
#' @param competition choice between symmetric, asymmetric and hierarchical competition; in the two latter cases species with higher trait values put pressure on species with lower trait values (species with lower trait values do not influence species with higher trait values), under asymmetic competition, competitive strength depends on niche overlap (i.e. very different species do not compete), under hierarchical competition, niche overlap is unimportant
#' @param intra.sp.com assigns the strength of intraspecific competition; the value should range between 0 (no intraspecific competition) and 1 (intraspecific competition always higher than interspecific competition) 
#' @return
#'   List of objects: 
#'   \describe{
#'   \item{community}{ vector of individuals (names of species in community) in the final year) }
#'   \item{abundances}{ species abundance matrix in the final year }
#'   \item{communities.over.time}{ matrix of all individuals (identified by their species name) over all time-steps } 
#'   }
#'   
#'   The plots on the left side show the influence of each of the three filters (y-axes) on each species (x-axes) over time. The plots on the right side show the influence of the overall filter effects (top), the abundances distribution and the dynamics of mpd (mean functional pariwise distance from picante) over time (green horizontal lines indicate the range of random expectations).   
#' @details
#'   Community assembly is simulated by asynchroneous updating of K individuals per year. For each update a random individual is removed and replaced by an individual from the species pool. The choice of that individual follows a multinomial distribution, with probabilities being driven by an environmental filter, a competition filter and a recruitment filter. 
#'   
#'   Environmental filter, p.env = dnorm(niche.optimum[i], mean=env, sd=niche.breadth) / dnorm(env, mean=env, sd=niche.breadth), with i being the species identity  
#'   
#'   Competition filter, p.comp= 1 - sum(alpha_ij * N_j)/sum(N_j), with j being the individuals already in the community and alpha_ij the niche-overlap 
#'   
#'   Recruitment filter, p.abun = species.pool.abundance_i + beta.abun * abun_i, with abun beeing the species i abundance in community 
#'   
#'   The probability of a species to enter the community is proportional to: p.all = exp(beta.env * log(p.env) + beta.comp * log(p.comp) + log(p.abun))
#' @seealso \code{\link{simulation.experiment}} for running simulation experiments with an integrated simulation of species pools and diversity analyses of final community structures
#' @examples
#' # create species pool
#' pool <- create.pool(n.species.pool=500, n.invader.pool=1, evol.model="deltaTree", min.phyl.signal=NA, evol.model.param=0.1, nrep=1)   
#' niche.optima <- pool$func$niche_evol
#' names(niche.optima) <- pool$func$SpeciesID
#' 
#' # run community assembly either with only an environmental filter or with only a competition filter
#' env.filter.community <- tamaure(niche.breadth=2, niche.optima=niche.optima, env=50, beta.env=1, beta.comp=0, beta.abun=1, years=10, K=100, plot=TRUE)
#' competition.community <- tamaure(niche.breadth=2, niche.optima=niche.optima, env=50, beta.env=0, beta.comp=10, beta.abun=1, years=10, K=100, plot=TRUE)
#' @keywords assembly invasion competition habitat filtering
#' @export

tamaure <- function(niche.breadth = 5, niche.optima, env, beta.env = 0, beta.comp = 0, beta.abun = 0, years = 20, K = 20, community.in = NA, species.pool.abundance = NA, plot = FALSE, competititon="symmetric", intra.sp.com=1) {
    # ------- 1. getting input -------
    species.count <- length(niche.optima)  # nb of species in the species pool
    if (is.na(species.pool.abundance[1])) 
        species.pool.abundance <- rep(1, species.count)
    
    # ------- 2. Creating the distance matrix from the niche overlap -------
    species.niche.dist <- as.matrix(dist(niche.optima, diag = TRUE, upper = TRUE))  # to calculate mpd
    # symmetric competition
    species.niche.overlap.sym <- outer(niche.optima, niche.optima, function(x, y) {
        2 * pnorm(-abs((x - y))/2, mean = 0, sd = niche.breadth)
    })
    # asymmetric competition
    species.niche.overlap.asym <- outer(niche.optima, niche.optima, function(x, y) {
        sign <- ifelse(x > y, 1, 0)
        overlap <- 2 * pnorm(-abs((x - y))/2, mean = 0, sd = niche.breadth)
        sign * overlap
    }) 
    # hierarchical competition
    species.comp.hierarchy <- outer(niche.optima, niche.optima, function(x, y) {
        ifelse(x > y, 1, 0)
    }) 

    diag(species.niche.overlap.sym) <- intra.sp.com
    diag(species.niche.overlap.asym) <- intra.sp.com 
    diag(species.comp.hierarchy) <- intra.sp.com 
 
    # ------- 3. The loop -------
    out.community <- as.data.frame(matrix(NA, nrow = years + 1, ncol = K))  # to store community composition over time
    if (plot == TRUE) {
        community.null <- sample(names(niche.optima), K, replace = TRUE)
        mpd.now <- mean(species.niche.dist[community.null, community.null])
    }
    
    # initialization: random or community.in
    if (is.na(community.in[1])) {
        community <- out.community[1, ] <- sample(names(niche.optima), K, replace = TRUE)
    } else {
        community <- out.community[1, ] <- community.in
    }
    
    # environmental filter (fix through time)
    log.p.env <- dnorm(x = niche.optima, mean = env, sd = niche.breadth, log = TRUE) - dnorm(x = env, mean = env, sd = niche.breadth, log = TRUE)
    
    # loop - years
    for (year in 1:years) {
        
        # loop - asynchroneous updating within years
        for (i in 1:K) {
            abundance <- as.numeric(table(community)[names(niche.optima)])
            abundance <- ifelse(is.na(abundance), 0, abundance)
            if(competititon=="symmetric") p.comp <- 1 - colSums(species.niche.overlap.sym[community, ])/K
            if(competititon=="asymmetric") p.comp <- 1 - colSums(species.niche.overlap.asym[community, ])/K
            if(competititon=="hierarchical") p.comp <- 1 - colSums(species.comp.hierarchy[community, ])/K
            p.all <- exp(beta.env * log.p.env + beta.comp * log(p.comp) + log(species.pool.abundance + beta.abun * abundance))
            p.all <- ifelse(is.na(p.all), min(p.all, na.rm = TRUE), p.all)
            if (sd(p.all, na.rm = TRUE) == 0) 
                p.all = NULL
            out <- sample(seq(community), 1)  # Choose randomly one individual to exclude
            community[out] <- sample(names(niche.optima), 1, prob = p.all)  # Replace this excluded species by lottery competition
        }
        out.community[year + 1, ] <- community
        if (plot == TRUE) 
            mpd.now <- c(mpd.now, mean(species.niche.dist[community, community]))
        # ------- 4. Plot (if requested) -------
        if (plot == TRUE & year == years) {
            n.rand <- round(years/4) + 1
            if (is.null(p.all)) 
                p.all <- rep(0, species.count)
            filters = cbind(comp = round(p.comp, 2), env = round(exp(log.p.env), 2), abun = round(abundance, 2), all = round(p.all, 2))
            mpd.now <- c(mpd.now, mean(species.niche.dist[community, community]))
            mypar <- par(mfcol = c(3, 2))
            index <- order(niche.optima)
            barplot(filters[, 1][index], xlab = "Species ID", ylab = "p.comp")
            barplot(filters[, 2][index], xlab = "Species ID", ylab = "p.env")
            barplot(filters[, 3][index], xlab = "Species ID", ylab = "abundance")
            barplot(filters[, 4][index], xlab = "Species ID", ylab = "p.all")
            hist(table(community), xlab = "Abundance", ylab = "Frequency", main = "")
            mpd.null <- sapply(1:1000, function(x) {
                community.null <- sample(names(niche.optima), K, replace = TRUE)
                mean(species.niche.dist[community.null, community.null])
            })
            low.border <- sort(mpd.null)[50]
            high.border <- sort(mpd.null)[950]
            plot((1:length(mpd.now)), mpd.now, type = "l", xlab = "Time", ylab = "mpd", ylim = c(min(c(low.border, mpd.now)), max(c(high.border, mpd.now))))
            abline(h = c(low.border, high.border), col = 3)
            par(mypar)
        }
    }
    
    # ------- 5. Prepare the output -------
    abundance.vector <- as.numeric(table(community)[names(niche.optima)])
    abundance.vector <- ifelse(is.na(abundance.vector), 0, abundance.vector)
    names(abundance.vector) <- names(niche.optima)
    return(list(community = community, abundances = abundance.vector, communities.over.time = out.community))
}
