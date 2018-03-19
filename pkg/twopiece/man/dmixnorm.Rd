\name{dnorm}
\alias{dmixnorm}
\alias{rmixnorm}
\title{ Density and random number generation for
Normal mixtures}
\description{
  \code{dmixnorm} evaluate the pdf for a Normal mixture,
  \code{rmixnorm} draws random numbers
}
\usage{
dmixnorm(x,mu,Sigma,probs,logscale=FALSE)

rmixnorm(n,mu,Sigma,probs)
}
\arguments{
  \item{x}{Vector, matrix or data.frame with values at which to evaluate the density}

  \item{n}{Number of random draws to be obtained}

  \item{mu}{A list of location parameters for each component}

  \item{Sigma}{A list of scale matrices for each component}

  \item{probs}{Vector with mixture component weights}

  \item{logscale}{If set to \code{TRUE} log-density is returned}
}
\value{
  \code{dmixnorm} returns a vector with the probability density
  function evaluated at the given \code{x}.

  \code{rmixnorm} return a matrix with the generated random draws.
}
\examples{
library(twopiece)

n <- 100; probs <- c(2/3,1/3); param <- 'eps'
mu1 <- c(0,0); S1 <- matrix(c(1,0,0,1),nrow=2)
mu2 <- c(3,3); S2 <- matrix(c(1,.5,.5,1),nrow=2)
mu <- list(mu1,mu2); Sigma <- list(S1,S2)
xsim <- rmixnorm(n,mu=mu,Sigma=Sigma,probs=probs)

head(xsim$x) #simulated values

head(xsim$cluster) #true cluster
}
\keyword{ distribution }
\keyword{ cluster }
