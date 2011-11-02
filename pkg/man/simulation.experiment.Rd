\name{simulation.experiment}
\alias{simulation.experiment}

\title{ Runs simulation experiments for community assembly }
\description{
  This function creates first a species pool (with phylogenetic and functional attributes), second assembles communities (see \code{\link{tamaure}}) and third computes functional and phylogenetic diversity indices. After natural communities are assembled invasion can occur.
}
\usage{
simulation.experiment(parameters)
}

\arguments{
  \item{parameters}{ vector of input parameter values (see details)}
  }


\value{
  Several objects are created and returned as a list of elements: 
  
  \item{parameter}{ a vector of the parameters used in the simulation }
  \item{pool}{ a list of different objects that define the species pool. "Traits" is a dataframe of the trait values for each species and the identity of the invaders. "phylo" is a phylo object of the phylogenetic relationship between the species of the species pool. "K" is the value of the phylogenetic signal for the evolved traits. "invader" is the list of element that characterize the invaders: Its identity (ID), trait value (opt), and its distance to the environment of the community (distToOpt)  }
  \item{native}{a list of two elements "communities" and "indices". communities is a matrix of species (columns) abundance (i.e. number of individuals) par sites (rows). indices is a list of matrices, each matrix represent the value of different indices (columns) per communities (see details for more precision on the indices) }
  \item{invasion}{ a list of three elements "communities", "process", and "indices". communities and indices are the equivalent of "communities" and "indices" for the native communities but applied to invaded communities. Process is a dataframe showing the presence of the invaders in the communities (rows) at each year (columns) }
 
  
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
  
  
}

\author{ Laure Gallien, Tamara Münkemüller }
\seealso{ \code{\link{tamaure}}, \code{\link{trait.evolution}}, \code{\link{get.attractors}} }

\examples{
  #################################################################################################
  # parameters (needed for all examples!)
  #################################################################################################
  # general
  years <- 10					
  invasion.time <- 0				
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
  simulation.experiment(param[1,]) # run the first line 
  # more lines with lapply
  results <- lapply(1:2, function(x) simulation.experiment(param[x,])) 
  

}

\keyword{ assembly }
\keyword{ invasion }
\keyword{ competition }
\keyword{ habitat filtering }