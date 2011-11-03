\name{simulation.experiment}
\alias{simulation.experiment}

\title{ Runs simulation experiments for community assembly }
\description{
  This function creates (1) a species pool (with phylogenetic and functional attributes), (2) assembles communities (see \code{\link{tamaure}}) and (3) computes functional and phylogenetic diversity indices. After natural communities are assembled invasion can occur.
}
\usage{
  simulation.experiment(parameters)
}

\arguments{
  \item{parameters}{ vector of input parameter values (see details and examples)}
  }


\value{
  Several objects are created and returned as a list of elements: 
  
  \item{parameter}{ a vector of the parameters used in the simulation }
  \item{pool}{ a list of different objects that define the species pool. "func" is a dataframe of the trait values for each species and the identity of the invaders. "phylo" is a phylo object of the phylogenetic relationship between the species of the species pool. "invader" is a list that contains the ID of the invader, its niche optimum and its distance to the environmental optimum. }
  \item{natives}{a list of two elements "all.abundances" and "indices". all.abundances is a sites X species matrix of abundances. indices is a list of matrices containing functional and phylohgenetic diversity information (see details) }
  \item{invaders}{ like the 'natives' object but for invasion }
 
  
}

\details{
  The parameter vector needs to contain the following (see also \code{\link{tamaure}}):
  \item{niche.breadth}{ value of standard deviation of the Gaussian distribution that describes the niches (for now the same for all species) }
  \item{env}{ value of the environmental optimum of the habitat the community is occupying  }
  \item{beta.env}{ value of the strength of the environmental filter (see details)  }
  \item{beta.comp}{ value of the strength of the competition filter (see details)  }
  \item{beta.abun}{ value of the strength of the local abundance filter, i.e. the advantage of already being present in the community (see details)  }
  \item{years}{ number of simulated years; If plots == TRUE equilibrium conditions can be checked by eye  }
  \item{K}{ value of carrying capacity, i.e. number of individuals in the community  }
  \item{species.pool.abundance}{ vector of abundances in species pool; if NA than all are set to 1  }
  \item{invasion.time}{ number of years simulated for invasion (after native community has established); if 0 than there is no invasion }
  \item{n.rep.null.model}{ number of repetitios for null model  }
  \item{n.species.pool}{ size of species pool  }
  \item{n.invader.pool}{ number of invaders in species pool  }
  \item{evol.model}{ choice of evolutionary model which determines phylogenetic signal; 1: BM, 2: OU, 3: delta  }
  \item{min.phyl.signal}{ minimum value of phylogenetic signal in species pool (NA means that the first evolved value is taken)  }
  \item{n.communities}{ number of simulated communities  }
  
  There are 5 matrices describing the diversity structur (5 for natives and 5 for invaders). They contain information on (1) observed diversity, (2) z-scores, (3) ranks, (4) mean of the null distribution and (5) standard deviation of the null distribution. Each matrix gives these information for all sites and for a list of diversity measures (mpd, mntd, ...).  
  
  
}

\author{ Laure Gallien, Tamara Munkemuller }
\seealso{ \code{\link{tamaure}}, \code{\link{trait.evolution}}, \code{\link{get.attractors}} }

\examples{
  #################################################################################################
  # parameters (needed for all examples!)
  #################################################################################################
  # general
  years <- 10					
  invasion.time <- 2				
  n.rep.null.model <- 10				
  n.rep.scenario <- 2					
  
  # species pool
  n.species.pool <- 200					
  n.invader.pool <- 1						
  niche.breadth <- c(5, 20)
  evol.model <- 1					
  min.phyl.signal <- NA				
  
  # communities
  n.communities <- 5					
  K <- 30					
  env <- 50				
  beta.env <- c(0,1)					
  beta.comp <- c(0,1)					
  beta.abun <- 1					
  species.pool.abundance <- NA
  
  # create table with parameters
  param <- expand.grid(rep.scenario=1:n.rep.scenario, K=K, n.species.pool=n.species.pool, n.invader.pool=n.invader.pool, beta.env=beta.env, 
      beta.comp=beta.comp, beta.abun=beta.abun, niche.breadth=niche.breadth, evol.model=evol.model, years=years, invasion.time=invasion.time, 
      n.communities=n.communities, n.rep.null.model=n.rep.null.model, env=env, min.phyl.signal=min.phyl.signal, species.pool.abundance=species.pool.abundance)
  param$simu.id <- 1:nrow(param)

  #################################################################################################
  # run the experiment
  #################################################################################################
  # relatively slow versions:
  # only one line
  single.run <- simulation.experiment(param[1,]) # run the first line
  str(single.run) 
  # more lines with lapply
  experiment.run <- lapply(1:2, function(x) simulation.experiment(param[x,])) 
  str(experiment.run)

}

\keyword{ assembly }
\keyword{ invasion }
\keyword{ competition }
\keyword{ habitat filtering }