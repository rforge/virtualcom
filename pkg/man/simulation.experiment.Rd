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
 
  \item{n.species.pool}{ size of species pool  }
  \item{evol.model}{ choice of evolutionary model (see also \code{\link{create.pool}}) }
  \item{evol.model.param}{ parameterization of evolutionary model (see also \code{\link{create.pool}}) }
  \item{species.pool.abundance}{ vector of abundances in species pool; if NA than all are equally abundant }  
  \item{min.phyl.signal}{ minimum value of phylogenetic signal in species pool (see also \code{\link{create.pool}}) }   
  \item{years}{ number of simulated timesteps; If plots == TRUE equilibrium conditions can be checked by eye  }  
  \item{n.communities}{ number of simulated communities}
  \item{env}{ value of the environmental optimum of the plot the community is occupying  }
  \item{K}{ value of carrying capacity, i.e. number of individuals in the community  }  
  \item{niche.breadth}{ value of standard deviation of the Gaussian distribution that describes the niche traits (for now the same for all species)  } 
  \item{beta.env}{ value of the strength of the environmental filter (see details) }
  \item{beta.comp}{ value of the strength of the competition filter (see details) }
  \item{beta.abun}{ value of the strength of the local abundance filter, i.e. the advantage of already being present in the community (see also \code{\link{tamaure}}) }   
  \item{invasion.time}{ number of time-steps simulated for invasion (after native community has established); if 0 than there is no invasion   }   
  \item{n.invader.pool}{ number of invaders in species pool  }                                                                                                                                 
  \item{InvDistrib}{ The invader's distribution in the phylogney can be either clusterd ("1"), random ("2") or  overdispersed ("3)  }                                                                                                                                 
  \item{n.rep.null.model}{ number of repetitions for null models for the diversity indices calculated to describe community composition  }
  \item{null.model}{ different null models can be choosen, see details  }   
}


\value{
  Information on input parameters, species pool, native assembly structure and invader assembly structure are returned as a list:   
  \item{parameter}{ a vector of the input parameters }
  \item{pool}{ a list of different objects that define the species pool. "func" is a dataframe of trait values and invader identity. "phylo" is a phylo object. "invader" is a list of the invaders ID, its niche optimum and its performance given the environmental optimum. }
  \item{natives}{a list of the native commuity structure and diversity measures. "all.abundances" is a site X species matrix of abundances. "indices" is a list of matrices containing functional and phylohgenetic diversity information (see details) }
  \item{invaders}{ a list of commuity structure and diversity measures after the invasion processes, see "natives" for structure }
}

\details{
  The following null models can be choosen to test for the structure of the native community: "1" = taxa.labels,"2" = richness, "3" = frequency, "4" = independentswap, "5" = trialswap. These models are explained in more detail in the function ses.mpd in the package picante. 
  The output lists for natives and invaders consist of 5 matrices describing the diversity structur. 
  They contain information on (1) observed diversity, and based on the chosen null models: (2) z-scores, (3) ranks, (4) mean of the null distribution and (5) standard deviation of the null distribution. Each matrix gives these information for all sites and for a list of diversity measures (mpd, mntd, ...).  

}

\author{ Laure Gallien, Tamara Münkemüller }
\seealso{ \code{\link{tamaure}} for the community assembly function}

\examples{
	# load pre-prepared parameter table
	data(param1)	
	# run the first line of the parameter table
	single.run <- simulation.experiment(param1[1,]) # run the first line
	str(single.run)
	
	# setting up a full experiment (this may take a few minutes) 
	# In this example traits evolve either on an untransformed or on a delta transformed tree. The latter should result in higher phylogenetic signal. Communities are assembled randomly, under a pure environmental filter or under a competition filter.
	wrapper <- function(a){
	  library(VirtualCom)
	  data(param1)
	  print(a)
	  return(try(simulation.experiment(param1[a,])))	
	}  
	# running the experiment either with lapply ...
	output1 <- lapply(1:nrow(param1), wrapper) 	
	# ... or running the experiment with sfLapply using snowfall to allow for parallel computing
	require(snowfall)
	sfInit(parallel=TRUE, cpus=11)
	output1 <- sfLapply(1:nrow(param1), wrapper)
	sfStop()
	
	# results put together in one table
	result.table1 <- getResults(output=output1, myvar="obs", invader="FALSE")
	result.table1$evol.model.param <- ifelse(result.table1$evol.model.param==1, "no", "yes")
	result.table1$process <- ifelse(result.table1$beta.env==0 & result.table1$beta.comp==0, "Random", ifelse(result.table1$beta.env!=0 & result.table1$beta.comp==0, "Env", ifelse(result.table1$beta.env==0 & result.table1$beta.comp!=0, "Comp", "Both") ))
	result.table1 <- result.table1[result.table1$process != "Both",]
	head(result.table1)
	Ks <- aggregate(result.table1$Blom_K, list(result.table1$evol.model.param), mean)
	
	# plot results for functional and phylogenetic mpd
	require(ggplot2)
	plot1 <- ggplot(data=result.table1, aes(x=process, y=FD_ab_mpd))+ geom_boxplot(width=0.8) +  xlab("Assembly rules") +  ylab("Functional diversity")   	
	plot2 <- ggplot(data=result.table1, aes(x=process, y=PD_ab_mpd)) + geom_boxplot(aes(fill= evol.model.param), width=0.8) + xlab("Assembly rules") +  ylab("Phylogenetic diversity") + theme(legend.justification=c(1,0), legend.position=c(0.42,0)) + scale_fill_discrete(name = "Tree transformed?", labels=c(paste("no; K~", round(Ks[1,2],2),sep=""), paste("yes; K~", round(Ks[2,2],2),sep=""))) 
	multiplot(plot1, plot2, cols=2)

  
}

\keyword{ assembly }
\keyword{ invasion }
\keyword{ competition }
\keyword{ habitat filtering }