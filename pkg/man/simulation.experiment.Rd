\name{simulation.experiment}
\alias{simulation.experiment}

\title{ Runs simulation experiments with eco-evolutionary model }
\description{
  This function creates (1) a species pool (with phylogenetic and functional attributes), (2) assembles communities (see \code{\link{tamaure}}) and (3) computes functional and phylogenetic diversity indices. After natural communities are assembled invasion can occur.
}
\usage{
  simulation.experiment(parameters)
}

\arguments{
  The only argument is a named vector of input parameter values with the following entries:
  
  \item{niche.breadth}{ value of standard deviation of the Gaussian distribution that describes the niches (for now the same for all species)  } 
  \item{env}{ value of the environmental optimum of the habitat the community is occupying  }
  \item{beta.env}{ value of the strength of the environmental filter (see details) }
  \item{beta.comp}{ value of the strength of the competition filter (see details) }
  \item{beta.abun}{ value of the strength of the local abundance filter, i.e. the advantage of already being present in the community (see also \code{\link{tamaure}}) }   
  \item{years}{ number of simulated years; If plots == TRUE equilibrium conditions can be checked by eye  }  
  \item{K}{ value of carrying capacity, i.e. number of individuals in the community  }  
  \item{species.pool.abundance}{ vector of abundances in species pool; if NA than all are set to 1 }  
  \item{invasion.time}{ number of years simulated for invasion (after native community has established); if 0 than there is no invasion   }  
  \item{n.rep.null.model}{ umber of repetitios for null model     }
  \item{n.species.pool}{ size of species pool  }
  \item{n.invader.pool}{ number of invaders in species pool  }                                                                                                                                 
  \item{evol.model}{ choice of evolutionary model which determines phylogenetic signal; 1: BM, 2: OU, 3: delta        }
  \item{min.phyl.signal}{ minimum value of phylogenetic signal in species pool (NA means that the first evolved value is taken)}   
  \item{n.communities}{ number of simulated communities}  
  
}


\value{
  Information on input parameters, species pool, native assembly structure and invader assembly structure are returned as a list : 
  
  \item{parameter}{ a vector of the input parameters }
  \item{pool}{ a list of different objects that define the species pool. "func" is a dataframe of trait values and invader identity. "phylo" is a phylo object. "invader" is a list of the invaders ID, its niche optimum and its performance given the environmental optimum. }
  \item{natives}{a list of the native commuity structure and diversity measures. "all.abundances" is a site X species matrix of abundances. "indices" is a list of matrices containing functional and phylohgenetic diversity information (see details) }
  \item{invaders}{ a list of commuity structure and diversity measures after the invasion processes, see "natives" for structure }
}

\details{
  The output lists for natives and invaders consist of 5 matrices describing the diversity structur. 
  They contain information on (1) observed diversity, (2) z-scores, (3) ranks, (4) mean of the null distribution and (5) standard deviation of the null distribution. Each matrix gives these information for all sites and for a list of diversity measures (mpd, mntd, ...).  

}

\author{ Laure Gallien, Tamara Munkemuller }
\seealso{ \code{\link{tamaure}} for the community assembly function}

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
  # only one line
  single.run <- simulation.experiment(param[1,]) # run the first line
  str(single.run) 
  
  # more lines 
  # preparation 
  wrapper <- function(a, path){
  	library(VirtualCom)
    load(paste(path,"param", sep=""))
	  return(simulation.experiment(param[a,]))	
  }
  path <- getwd()  
  save(param, file=paste(path,"param", sep="")) # save parameter table in file  
  
  # more lines with lapply
  experiment.run <- lapply(1:2, wrapper, path=path) 
  str(experiment.run)
  
  # more lines with snowfall for parallel computing (here for 2 CPUs)
  library(snowfall)
  sfInit(parallel=TRUE, cpus=2)
  stopifnot(sfCpus() == 2)
  stopifnot(sfParallel() == TRUE)
  sfStop()
  experiment.snow <- sfLapply(1:2, wrapper, path=path)
  str(experiment.snow)
  
}

\keyword{ assembly }
\keyword{ invasion }
\keyword{ competition }
\keyword{ habitat filtering }