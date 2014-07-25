\name{tamaure}
\alias{tamaure}

\title{ Community assembly from species pool  }
\description{
  This function simulates the spatially implicit assembly of individuals into a community over time. Individuals come from a pool of species. In the pool each species is described by a niche trait (mean and sd of a Gaussian distribution) and its frequency. 
}
\usage{
tamaure(niche.breadth=5, niche.optima, env, beta.env=0, beta.comp=0, beta.abun=0, years=20, K=20, community.in=NA, species.pool.abundance=NA, plot=FALSE, ...)

}

\arguments{
  \item{niche.optima}{ vector of niche optima (means of the Gaussian distributions that describe the niches) in species-pool }
  \item{niche.breadth}{ value of standard deviation of the Gaussian distributions that describe the niches (for now the same for all species) }
  \item{env}{ value of the environmental optimum of the plot the community is occupying  }
  \item{beta.env}{ value of the strength of the environmental filter (see details)  }
  \item{beta.comp}{ value of the strength of the competition filter (see details)  }
  \item{beta.abun}{ value of the strength of the local abundance filter, i.e. the advantage of already being present in the community (see details)  }
  \item{years}{ number of simulated time-steps; if plots == TRUE equilibrium conditions can be checked by eye  }
  \item{K}{ value of carrying capacity, i.e. number of individuals in the community  }
  \item{community.in}{ vector of individuals (described by species names) already in the community; if NA community is initialized randomly from species pool  }
  \item{species.pool.abundance}{ vector of abundances in species pool; if NA than all are equally abundant  }
  \item{plot}{ if TRUE than community composition is plotted }
  }


\value{
  List of objects: 
  \item{community}{ vector of individuals (names of species in community) in final year) }
  \item{abundances}{ matrix of species X abundances  }
  \item{communities.over.time}{ matrix of individuals (names of species in community) over time-steps } 
  
  The plots on the left side show the influence of each of the three filters (y-axes) on each species (x-axes) over time. The plots on the right side show the influence of the overall       filter (top), the abundances distribution and the dynamics of mpd (mean pariwise distance from picante) over time (green horizontal lines indicate the range of random expectations).   
}

\details{
  Community assembly is simulated by asynchroneous updating of K individuals per year. For each update a random individual is removed and replaced by an individual from the species   pool. The choice of that individual follows a multinomial distribution, with probabilities being driven by an environmental filter, a competition filter and a local abundance weight. 

  Environmental filter, p.env: dnorm(niche.optimum[i], mean=env, sd=niche.breadth) / dnorm(env, mean=env, sd=niche.breadth), with i  the species index  
  
  Competition filter, p.comp: 1 - sum(alpha_ij * N_j)/sum(N_j), with j beeing the individuals already in the community and alpha_ij the niche-overlap 
  
  Abundance weight, p.abun: species.pool.abundance_i + beta.abun * abun_i, with abun beeing abundance in community 
  
  The probability of a species to enter the community is proportional to: p.all = exp(beta.env * log(p.env) + beta.comp * log(p.comp) + log(p.abun))
}

\author{ Laure Gallien, Tamara Münkemüller }
\seealso{ \code{\link{simulation.experiment}} for running simulation experiments with an integrated simulation of species pools and diversity analyses of final community structures}

\examples{
# create species pool
pool <- create.pool(n.species.pool=500, n.invader.pool=1, evol.model="deltaTree", min.phyl.signal=NA, evol.model.param=0.1, nrep=1)   
niche.optima <- pool$func$niche_evol
names(niche.optima) <- pool$func$SpeciesID

# run community assembly either with only an environmental filter or with only a competition filter
env.filter.community <- tamaure(niche.breadth=2, niche.optima=niche.optima, env=50, beta.env=1, beta.comp=0, beta.abun=1, years=10, K=100, plot=TRUE)
competition.community <- tamaure(niche.breadth=2, niche.optima=niche.optima, env=50, beta.env=0, beta.comp=10, beta.abun=1, years=10, K=100, plot=TRUE)
}

\keyword{ assembly }
\keyword{ invasion }
\keyword{ competition }
\keyword{ habitat filtering }