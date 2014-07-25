\name{ create.pool }
\alias{ create.pool }

\title{ simulates a species pools }

\description{
  simulates a species pools containing a phylogeny, species trait niches and information on invasion status
}

\usage{
create.pool(n.species.pool, n.invader.pool, evol.model, evol.model.param, min.phyl.signal, nrep)
}

\arguments{
  \item{n.species.pool}{ number of species (natives and invaders) }
  \item{n.invader.pool}{ number of invaders }
  \item{evol.model}{ choice of evolutionary model which determines phylogenetic signal, see details  }
  \item{evol.model.param}{ when evol.model is different from Brownian motion, see details }
  \item{min.phyl.signal}{ minimal phylogenetic signal }
  \item{n.rep}{ number of repetition for the null models }
  }


\value{
  List of objects: 
  \item{func}{ Functional traits for each species of the pool }
  \item{phylo}{ Phylogenetic tree of the species pool  }
  \item{indices}{ Indices of the phylogenetic signal } 
  \item{invader}{ List of the invader's identities and trait values } 
  \item{parameters}{ The parameters used for creating the species pool } 
  }

\details{
   The different evolutionary models (argument evol.mod) are based on tree transformations see (\code{\link{rescale} for more details): 1= 	
	"1" = BM (Brownian motion), "2" = deltaTree, "3" = kappaTree, "4" = ouTree, 5 = "lambdaTree". The tree transformation can be parameterized using the argument evol.model.param. 
}

\author{ Laure Gallien, Tamara Münkemüller }
\seealso{ \code{\link{simulation.experiment}} for running simulation experiments with an integrated simulation of species pools and diversity analyses of final community structures}

\examples{


}

\keyword{ species pool }
\keyword{ invasion }
\keyword{ trait evolution }
