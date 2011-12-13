\name{create.pool.invasions}
\alias{create.pool.invasions}

\title{ Create several species pools containing more than 1 invader }
\description{
  Creates several species pools containing several invaders
}
\usage{
create.pool.invasions(paramPool, ...)

}

\arguments{
  \item{n.species.pool}{ number of species in the species pool }
  \item{n.invader.pool}{ number of invaders in the species pool }
  \item{evol.model}{ hoice of evolutionary model which determines phylogenetic signal; 1: BM, 2: OU, 3: delta   }
  \item{min.phyl.signal}{ minimal phylogenetic signal }
  \item{evol.model.param}{ For evol.model different from brownian motion }
  \item{n.rep.null.model}{ number of repetition for the null models }
  }


\value{
  List of objects: 
  \item{func}{ Functional traits for each species of the pool }
  \item{phylo}{ Phylogenetic tree of the species pool  }
  \item{indices}{ indices of the phylogenetic signal } 
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
