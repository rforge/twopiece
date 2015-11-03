\name{dskewt}
\alias{dskewt}
\alias{dmixskewt}
\alias{rskewt}
\alias{rmixskewt}
\title{ Density and random number generation for
multivariate skew-t and mixture of multivariate skew-t's}
\description{
  \code{dskewt} gives the multivariate skew-t density,
  \code{dmixskewt} for a mixture of multivariate skew-t densities.
  \code{rskewt} draws random values from the multivariate skew-t density,
  \code{rmixskewt} from a mixture of multiviariate skew-t's.
}
\usage{
dskewt(x,mu,Sigma,A,D,alpha,nu,param='eps',logscale=FALSE,ttype='independent')

dmixskewt(x,mu,Sigma,alpha,nu,probs,param='eps',logscale=FALSE,ttype='independent')

rskewt(n,mu,Sigma,alpha,nu,param='eps',ttype='independent')

rmixskewt(n,mu,Sigma,alpha,nu,probs,param='eps',ttype='independent')
}
\arguments{
  \item{x}{Vector, matrix or data.frame with values at which to evaluate the density}

  \item{n}{Number of random draws to be obtained}

  \item{mu}{For dskewt and rskewt mu is a vector indicating the location parameter.
  For dmixskewt and rmixskewt mu is a list of location parameters for each component}

  \item{Sigma}{For dskewt and rskewt Sigma is the scale matrix.
  For dmixskewt and rmixskewt it is a list of scale matrices for each component}

  \item{A}{Optional matrix containing eigenvectores of Sigma, i.e. Sigma= t(A) D A.
  For identifiability all elements in A[,1] must be positive, if any entry
  is negative then the corresponding row is changed sign}

  \item{D}{Optional diagonal matrix with eigenvalues of Sigma, i.e. Sigma= t(A) D A.
  For identifiability these are assumed to be in decreasing order.}

  \item{alpha}{For dskewt and rskewt alpha is the vector or asymmetry parameters.
  For dmixskewt and rmixskewt it is a list of asymmetry parameters for each component}

  \item{nu}{For \code{dskewt} and \code{rskewt} nu is a scalar indicating the degrees
  of freedom. For \code{dmixskewt} and \code{rmixskewt} nu is a list with the degrees
  of freedom for each component}

  \item{probs}{Vector with mixture component weights}

  \item{param}{Set to 'eps' for the epsilon-skew parametization, 'isf' for the inverse
  scale factor parameterization. See help(dtp3) for details}

  \item{ttype}{Set to 'independent' for iskew-t, to 'dependent' for dskew-t}

  \item{logscale}{If set to \code{TRUE} log-density is returned}
}
\details{
A draw from the iskew-t distribution is obtained as
y= sqrt(D) A x, where x is a vector with independent draws
from univariate t distributions.

A draw from the dskew-t distribution is obtained as
y= sqrt(D) A x, where x is a draw from a multivariate t distribution
with identity scale matrix. That is, the elements in x are dependent
but uncorrelated, in constrast with the iskew-t where they are independent.
}
\value{
  \code{dskewt} and \code{dmixskewt} return a vector with the probability density
  function evaluated at the given \code{x}.

  \code{rskewt} and \code{rmixskewt} return a matrix with the generated random draws.
}
\references{
Rossell D., Steel M.F.J. Continuous non-Gassian mixtures. In Handbook of Cluster Analysis, CRC Press.
}

\examples{
#Single skew-t
library(twopiece)
n <- 100; mu1 <- c(0,0); S1 <- matrix(c(1,0,0,1),nrow=2)
alpha1 <- c(0,0); nu1 <- 100
y1 <- rskewt(n,mu=mu1,Sigma=S1,alpha=alpha1,nu=nu1,param='eps')
xseq <- seq(-3,3,length=10)
xgrid <- expand.grid(xseq,xseq)
dy1 <- dskewt(xgrid,mu=mu1,Sigma=S1,alpha=alpha1,nu=nu1)

#Mixture of skew-t's
n <- 100; probs <- c(2/3,1/3); param <- 'eps'
mu1 <- c(0,0); S1 <- matrix(c(1,0,0,1),nrow=2)
alpha1 <- c(0,0); nu1 <- 100
mu2 <- c(3,3); S2 <- matrix(c(1,.5,.5,1),nrow=2)
alpha2 <- c(-.5,0); nu2 <- 100
mu <- list(mu1,mu2); Sigma <- list(S1,S2); alpha <- list(alpha1,alpha2); nu <- list(nu1,nu2)
xsim <- rmixskewt(n,mu=mu,Sigma=Sigma,alpha=alpha,nu=nu,probs=probs,param=param)

head(xsim$x) #simulated values

head(xsim$cluster) #true cluster
}
\keyword{ distribution }
\keyword{ cluster }
