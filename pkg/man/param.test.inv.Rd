\name{param.test.inv}
\alias{param.test.inv}

\title{ Wrapping function for creating the species pool, assembling the community and invading it }
\description{
  This function creates a species pool (with phylogenetic and functional attributes) and assembles the communities according to a given set of ecological rules
}
\usage{
param.test.inv(timesteps, invasionTime, NulRep, N, Ninv, niche.breath, evolModel, NicheVal, physig, nbCom, capacity, envOptimum, cstE, cstB, cstSR, comp, NN)
}

\arguments{
  \item{timesteps}{ timesteps in each simulation to get equilibrium conditions }
  \item{invasionTime}{ how many years can invaders try to invade }
  \item{NulRep}{ how many repetitions in null models }
  \item{N}{  size of species pool}
  \item{Ninv}{ number of invaders in species pool}
  \item{niche.breath}{ widths of species niches }
  \item{evolModel}{ choice of evolutionary model which determines phylogenetic signal. Use 1 for Brownian motion, 2 for OU, and 3 for delta }
  \item{NicheVal}{ take the evolved traits or equidistant traits. Use 1 for evolved traits, and 2 for equidistant ones }
  \item{physig}{ min value of phylogenetic signal in species pool }
  \item{nbCom}{ number of sinulated communities }
  \item{capacity}{  carrying capacity of communities (number of individuals) }
  \item{envOptimum}{ optimal value of environmental conditions for communities }
  \item{cstE}{ strength of environmental filter }
  \item{cstB}{ strength of biotic filter  }
  \item{cstSR}{ strength of local seed reproduction }
  \item{comp}{different descriptions of biotic filter (competition). Use 1 for niche overlap, 2 for competition with the nearest neighbour only, and 3 for niche overlap among a given number of nearest neighbours (use NN to set the nuber of neighbours for the competition) }
  \item{NN}{ if comp=3, how many neighbours included for biotic filter (comeptition)}
  }


\value{
  Several objects are created and returned as a list of elements: 
  
  \item{parameter}{ a vector of the parameters used in the simulation }
  \item{pool}{ a list of different objects that define the species pool. "Traits" is a dataframe of the trait values for each species and the identity of the invaders. "phylo" is a phylo object of the phylogenetic relationship between the species of the species pool. "K" is the value of the phylogenetic signal for the evolved traits. "invader" is the list of element that characterize the invaders: Its identity (ID), trait value (opt), and its distance to the environment of the community (distToOpt)  }
  \item{native}{a list of two elements "communities" and "indices". communities is a matrix of species (columns) abundance (i.e. number of individuals) par sites (rows). indices is a list of matrices, each matrix represent the value of different indices (columns) per communities (see details for more precision on the indices) }
  \item{invasion}{ a list of three elements "communities", "process", and "indices". communities and indices are the equivalent of "communities" and "indices" for the native communities but applied to invaded communities. Process is a dataframe showing the presence of the invaders in the communities (rows) at each year (columns) }
 
  
}

\details{
  Several indices are computed for describing the native species communities:
  \item{mpd}{ mean pairwise distances}
  \item{mntd}{ mean of the nearest taxon distances}
  \item{mntd_med}{ median of the nearest taxon distances}
  
  \item{fun or phy}{ when using functional or phylogenetic distances}
  \item{ab}{ when the distances are weighted by the abundances}
 
  
}

\author{ Laure Gallien, Tamara Münkemüller }
\seealso{ \code{\link{tamaure}}, \code{\link{trait.evolution}}, \code{\link{get.attractors}} }

\examples{
}

\keyword{ assembly }
\keyword{ invasion }
\keyword{ competition }
\keyword{ habitat filtering }