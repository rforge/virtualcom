\name{create.pool.invasion}
\alias{create.pool.invasion}

\title{ Create several species pools containing more than 1 invader }

\description{
  Create several species pools containing several invaders.
}

\usage{
create.pool.invasion(n.species.pool, n.invader.pool, env, evol.model, evol.model.param, min.phyl.signal, nreps, ...)
}

\arguments{
  \item{n.species.pool}{ number of species in the species pool }
  \item{n.invader.pool}{ number of invaders in the species pool }
  \item{env}{ trait value of the phylogenetic root, used when evolving the trait on the phylogeny }
  \item{evol.model}{ choice of evolutionary model which determines phylogenetic signal; 1: BM, 2: delta, 3:Kappa, 4: OU }
  \item{evol.model.param}{ only when evol.model different from brownian motion. It is the alpha parameter }
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


\author{ Laure Gallien }
\seealso{ \code{\link{simulation.experiment}} for complete simulation experiments with simulation of species pool and diversity analyses of final community structures}

\examples{


}

\keyword{ species pool }
\keyword{ invasion }
\keyword{ trait evolution }
