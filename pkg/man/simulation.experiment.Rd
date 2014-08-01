\name{simulation.experiment}
\alias{simulation.experiment}

\title{ Runs simulation experiments with VirtualCom }
\description{
  This function creates (1) a species pool (see \code{\link{create.pool}}), (2) assembles communities (see \code{\link{tamaure}}) and (3) computes functional and phylogenetic diversity indices. After natural communities are assembled invasion can occur.
}
\usage{
  simulation.experiment(parameters)
}

\arguments{
  The only argument is a named vector of input parameter values with the following entries:
 
  \item{n.species.pool}{ number of species in the species pool (including native and invasive species) }
  \item{evol.model}{ choice of the trait evolution model (see also \code{\link{create.pool}}) }
  \item{evol.model.param}{ parameterization of the trait evolution model (see also \code{\link{create.pool}}) }
  \item{species.pool.abundance}{ vector of abundances in the species pool; if NA then all species are equally abundant }  
  \item{min.phyl.signal}{ minimum value of phylogenetic signal in species pool (see also \code{\link{create.pool}}) }   
  \item{years}{ number of simulated timesteps; If plots=TRUE equilibrium conditions can be checked by eye  }  
  \item{n.communities}{ number of simulated communities}
  \item{env}{ value of the environment in the simulated community (e.g. temperature) }
  \item{K}{ value of carrying capacity, i.e. number of individuals in the community  }  
  \item{niche.breadth}{ value of standard deviation of the Gaussian distribution that describes the niche traits (identical for all species)  } 
  \item{beta.env}{ value of the strength of the environmental filter (see details) }
  \item{beta.comp}{ value of the strength of the competition filter (see details) }
  \item{beta.abun}{ value of the strength of the recruitment filter, i.e. the advantage of already being present in the community (see also \code{\link{tamaure}}) }   
  \item{invasion.time}{ number of time-steps to simulate invasion (after native community has established); if 0 than there is no invasion }   
  \item{n.invader.pool}{ number of invaders in species pool  }                                                                                                                                 
  \item{InvDistrib}{ The invader's distribution in the phylogney can be either clusterd ("1"), random ("2") or  overdispersed ("3)  }                                                                                                                                 
  \item{n.rep.null.model}{ number of repetitions for null models to test the diversity indices describing community composition  }
  \item{null.model}{ different null models can be choosen, see details  }   
}


\value{
  Information on input parameters, species pool, native assembly structure and invader assembly structure are returned as a list:   
  \item{parameter}{ a vector of the input parameters }
  \item{pool}{ a list of different objects that define the species pool. "func" is a dataframe of trait values and invader identity. "phylo" is a phylo object. "invader" is a list of the invader's ID, their niche optima and their performances given the environment in the community }
  \item{natives}{a list of the native commuity structure and diversity measures. "all.abundances" is a site-by-species matrix of abundances. "indices" is a list of matrices containing functional and phylohgenetic diversity information (see details) }
  \item{invaders}{ a list of community structure and diversity measures after the invasion processes, see "natives" for the structure }
}

\details{
  The following null models can be choosen to test for the structure of the native community: "1" = taxa.labels,"2" = richness, "3" = frequency, "4" = independentswap, "5" = trialswap. These models are explained in more detail in the function ses.mpd in the package picante. 
  The output lists for natives and invaders consist of 5 matrices describing the diversity structure. 
  They contain information on (1) observed diversity, and based on the chosen null models: (2) z-scores, (3) ranks, (4) mean of the null distribution and (5) standard deviation of the null distribution. Each matrix gives this information for all sites and for a set of diversity indices (mpd, mntd, ...).  

}

\author{ Laure Gallien, Tamara Münkemüller }
\seealso{ \code{\link{tamaure}} for the community assembly function}

\examples{
	# load pre-prepared parameter table
	data(simple_param)	
	# run the first line of the parameter table
	single.run <- simulation.experiment(simple_param[1,]) # run the first line
	str(single.run)
	
	# setting up a full experiment (this may take a few minutes) 
	wrapper <- function(a){
	  library(VirtualCom)
	  data(simple_param)
	  print(a)
	  return(try(simulation.experiment(simple_param[a,])))	
	}  
	
	# running the experiment either with lapply
	output <- lapply(1:nrow(simple_param), wrapper) 
	
	# # or running the experiment with sfLapply using snowfall to allow for parallel computing
	# require(snowfall)
	# sfInit(parallel=TRUE, cpus=11)
	# output <- sfLapply(1:nrow(simple_param), wrapper)
	# sfStop()
	
	# extract results from list
    result.table <- get.results(output=output, myvar="obs", invader="FALSE")
	result.table$process <- ifelse(result.table$beta.env==0 & result.table$beta.comp==0, "Random", ifelse(result.table$beta.env!=0 & result.table$beta.comp==0, "Env", ifelse(result.table$beta.env==0 & result.table$beta.comp!=0, "Comp", "Both") ))

	# plot the funtional diversity (mpd) of communities in dependence on assembly processes
	require(ggplot2)
	ggplot(data=result.table, aes(x=process, y=FD_ab_mpd)) + geom_boxplot(width=0.8) +  xlab("Assembly rules") +  ylab("Functional diversity")
	
}

\keyword{ assembly }
\keyword{ invasion }
\keyword{ competition }
\keyword{ habitat filtering }