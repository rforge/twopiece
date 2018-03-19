\name{clusterprobs}
\alias{clusterprobs}
\alias{clusterprobs-methods}
\alias{clusterprobs,normFit-method}
\alias{clusterprobs,skewtFit-method}

\title{
Posterior cluster probabilities based on mixture of skew-t distributions
}
\description{
Posterior cluster probabilities from a mixture of multivariate
skew-t distributions, obtained by averaging over Gibbs draws.
}
\usage{
  clusterprobs(fit, x, iter)
}
\arguments{
  \item{fit}{Object of class normFit or skewtFit containing Gibbs draws,
  as returned by mixnormGibbs or mixskewtGibbs}
  \item{x}{Data points at which to evaluate the probabilities.
  \code{ncol(x)} must match \code{ncol(x)} when calling \code{mixskewtGibbs}}
  \item{iter}{To avoid excessive time only a subset of the Gibbs iterations are
  used. \code{iter} can be used to indicate the iteration indices, by default
  1000 equally space iterations are used}
}
\value{
  Matrix with cluster probabilities (individuals in rows, clusters in columns)
}
\author{
David Rossell
}
\seealso{
\code{\link{mixnormGibbs}}, \code{\link{mixskewtGibbs}}
}
\examples{
#See help(mixnormGibbs)
#See help(mixskewtGibbs)
}
\keyword{ cluster }
\keyword{ models }
\keyword{ distribution }
