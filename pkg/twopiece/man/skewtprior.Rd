\name{skewtprior}
\alias{skewtprior}
\title{Prior parameters for mixture of multivariate skew-t's.
}
\description{
skewtprior creates a list with prior parameters required by mixskewtGibbs.
The prior is given by

 (mu,Sigma) ~ N(mu; m,g*Sigma) * IWishart(Sigma; Q, q)

 0.5*(alpha+1) ~ Beta(a,b)

 nu is discrete with P(nu=j)=nuprobs[j]

 probs ~ Symmetric Dirichlet(r), where r is a scalar
}
\usage{
skewtprior(p, G, m=rep(0,p), g=1, Q=diag(p), q=p+1, a=2, b=2, r=1/G,
nuprobs=priornuJS(1:30,k=2.78,numax=30))
}
\arguments{
  \item{p}{Number of variables in the observed data}
  \item{G}{Number of mixture components}
  \item{m}{Location parameter for prior on mu given Sigma}
  \item{g}{Scale parameter for prior on mu given Sigma}
  \item{Q}{Scale matrix for Inverse Wishart prior on Sigma}
  \item{q}{Degrees of freedom of Inverse Wishart prior on Sigma}
  \item{a}{Prior on (alpha+1)/2 is Beta(a,b)}
  \item{b}{Prior on (alpha+1)/2 is Beta(a,b)}
  \item{r}{Prior on mixture component weights is Dirichlet(r)}
  \item{nuprobs}{vector where P(nu=j)=nuprobs[j]. Names can be provided,
  else it is assumed that support of nu ranges from 1 to length(nuprobs)}
}
\details{
The Juarez-Steel prior is proportional to k*nu/(nu+k)^3.

The Villa-Walker prior is proportional to

(log(nu+1)-log(nu))/2 + lbeta(0.5,(nu+1)/2) - lbeta(0.5,nu/2)
- 0.5*(nu+1)*int1 + 0.5*(nu+2)*int2

where f1= E(log(1+x^2/nu)), f2= E(log(1+x^2/(nu+1))) and
x follows a standard t with nu degrees of freedom.
}
\value{
List with elements mu, Sigma, alpha, nu and probs containing the prior parameters.
}
\references{
Rossell D., Steel M.F.J. Continuous mixtures with skewness and heavy
tails. Handbook of Mixture Analysis, Chapter 10, CRC Press 2018.
}
\examples{
##See help(mixskewtGibbs)
}
\keyword{ distribution }
