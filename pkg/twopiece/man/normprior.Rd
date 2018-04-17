\name{normprior}
\alias{normprior}
\title{Prior parameters for a multivariate Normal mixture
}
\description{
normprior creates a list with prior parameters required by mixnormGibbs.
The prior is given by

 (mu,Sigma) ~ N(mu; m,g*Sigma) * IWishart(Sigma; Q, q)

 probs ~ Symmetric Dirichlet(r), where r is a scalar
}
\usage{
normprior(p,G,m=rep(0,p),g=1,Q=diag(p),q=p+1,r=1/G)
}
\arguments{
  \item{p}{Number of variables in the observed data}
  \item{G}{Number of mixture components}
  \item{m}{Location parameter for prior on mu given Sigma}
  \item{g}{Scale parameter for prior on mu given Sigma}
  \item{Q}{Scale matrix for Inverse Wishart prior on Sigma}
  \item{q}{Degrees of freedom of Inverse Wishart prior on Sigma}
  \item{r}{Prior on mixture component weights is Dirichlet(r)}
}
\value{
List with elements mu, Sigma and probs containing the prior parameters.
}
\examples{
##See help(mixnormGibbs)
}
\keyword{ distribution }
