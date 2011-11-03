\name{tamaure}
\alias{tamaure}

\title{ Assembly of communities from species pool with filters (competition, environment, abundance) }
\description{
  This function simulates the spatially implicit assembly of individuals into a community over time. Individuals come from a pool of species. In the pool each species is described by its niche (mean and sd of a Gaussian distribution) and its frequency. 
}
\usage{
tamaure(niche.breadth=5, niche.optima, env, beta.env=0, beta.comp=0, beta.abun=0, years=20, K=20, community.in=NA, species.pool.abundance=NA, plot=FALSE, ...)

}

\arguments{
  \item{niche.optima}{ vector of niche optima (means of the Gaussian distribution) in species-pool }
  \item{niche.breadth}{ value of standard deviation of the Gaussian distribution that describes the niches (for now the same for all species) }
  \item{env}{ value of the environmental optimum of the habitat the community is occupying  }
  \item{beta.env}{ value of the strength of the environmental filter (see details)  }
  \item{beta.comp}{ value of the strength of the competition filter (see details)  }
  \item{beta.abun}{ value of the strength of the local abundance filter, i.e. the advantage of already being present in the community (see details)  }
  \item{years}{ number of simulated years; If plots == TRUE equilibrium conditions can be checked by eye  }
  \item{K}{ value of carrying capacity, i.e. number of individuals in the community  }
  \item{community.in}{ vector of individuals (described by species names) already in the community; if NA than community is initialized randomly from species pool  }
  \item{species.pool.abundance}{ vector of abundances in species pool; if NA than all are set to 1  }
  \item{plot}{ if TRUE than development of filters and community composition is plotted over time  }
  }


\value{
  List of objects: 
  \item{community}{ vector of individuals (names of species in community in final year) }
  \item{abundances}{ matrix of species X abundances  }
  \item{all.communities}{ matrix of individuals (names of species in community) X years } 
  
  The plots on the left side show how the distribution of the three filters (besides the environmental filter which is kept constant) change over time. The plots on the right side show the dynamics of the overall filter (proportional to the probabilities in the multinomial distribution), the abundances distribution dynamics and the dynamics of mpd (mean pariwise distance from picante) over time (note that the negative values on the x-axis refer to random community assembly to fascilitate comparison).   
}

\details{
  Community assembly is simulated by asynchroneous updating of K individuals per year. For each update a random individual is removed and replaced by an individual from the species pool. The choice of that individual follows a multinomial distribution, with probabilities being driven by an environmental filter, a competition filter and a local abundance weight. 
  
  Environmental filter of species i, p.env: dnorm(niche.optimum[i], mean=env, sd=niche.breadth) / dnorm(env, mean=env, sd=niche.breadth)  
  Competition filter of species i, p.comp: 1 - sum(alpha_ij * N_j)/sum(N_j), with j beeing the individuals already in the community and alpha_ij the niche-overlap
  Abundance weight of species i, p.abun: species.pool.abundance_i + beta.abun * abun_i, with abun beeing abundance in community
  
  The probability of a species to enter the community is proportional to:
  p.all = exp(beta.env * log(p.env) + beta.comp * log(p.comp) + log(p.abun))
}

\author{ Laure Gallien, Tamara Munkemuller }
\seealso{ \code{\link{trait.evolution}}, \code{\link{get.attractors}} }

\examples{
  niche.optima <- rnorm(200, 100, 20) # create traits of species in species pool
  niche.optima <- sort(niche.optima) # nicer for the plots
  names(niche.optima) <- seq(niche.optima)
  
  # only environmental filter
  tamaure(beta.env=10, beta.comp=0, beta.abun=0, years=100, niche.breadth=5, niche.optima=niche.optima, env=100, K=50, communityIn=NA, plot=TRUE)
  # only competition filter; there is a border effect!
  tamaure(beta.env=0, beta.comp=10, beta.abun=0, years=100, niche.breadth=5, niche.optima=niche.optima, env=100, K=50, communityIn=NA, plot=TRUE)
  # only abundance weight; increased importance of local abundance changes variance of mpd!; very slow convergence
  tamaure(beta.env=0, beta.comp=0, beta.abun=100, years=100, niche.breadth=5, niche.optima=niche.optima, env=100, K=50, communityIn=NA, plot=TRUE)


}

\keyword{ assembly }
\keyword{ invasion }
\keyword{ competition }
\keyword{ habitat filtering }