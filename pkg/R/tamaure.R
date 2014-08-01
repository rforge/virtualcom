tamaure <- function(niche.breadth = 5, niche.optima, env, beta.env = 0, beta.comp = 0, beta.abun = 0, years = 20, K = 20, community.in = NA, species.pool.abundance = NA, 
    plot = FALSE) {
    # ------- 1. getting input -------
    species.count <- length(niche.optima)  # nb of species in the species pool
    if (is.na(species.pool.abundance[1])) 
        species.pool.abundance <- rep(1, species.count)
    
    # ------- 2. Creating the distance matrix from the niche overlap -------
    species.niche.dist <- as.matrix(dist(niche.optima, diag = TRUE, upper = TRUE))  # to calculate mpd
    species.niche.overlap <- outer(niche.optima, niche.optima, function(x, y) {
        2 * pnorm(-abs((x - y))/2, mean = 0, sd = niche.breadth)
    })
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
            p.comp <- 1 - colSums(species.niche.overlap[community, ])/K
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
