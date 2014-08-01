\name{ get.results }
\alias{ get.results }

\title{ result collection  }

\description{
  extracts results from the list that is given by simulation.experiment
}

\usage{
get.results(output, myvar, invader="FALSE")
}

\arguments{
  \item{output}{ list of results as it is given by simulation experiment }
  \item{myvar}{ variables that should be extracted. can be either "obs" for observed values, "zNULL" for z-values or "rankNULL" for ranks}
  \item{invader}{ should variables be extracted for the native or the invasive community? Default is FALSE and thus extraction for native community  }
  }


\value{
  Dataframe
  }

\author{ Tamara Münkemüller }
\seealso{ \code{\link{simulation.experiment}} for running simulation experiments with an integrated simulation of species pools and diversity analyses of final community structures}

\examples{
	# setting up a full experiment (this may take a few minutes) 
	wrapper <- function(a){
	  library(VirtualCom)
	  data(simple_param)
	  print(a)
	  return(try(simulation.experiment(simple_param[a,])))	
	}  
	
	# running the experiment either with lapply
	output <- lapply(1:nrow(simple_param), wrapper) 

	## running the experiment with sfLapply using snowfall to allow for parallel computing
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


