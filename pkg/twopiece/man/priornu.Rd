\name{priornuJS}
\alias{priornuJS}
\alias{priornuVW}
\title{ Discretized Juarez-Steel and Villa-Walker priors on degrees of freedom}
\description{
  \code{priornuJS} evaluates the Juarez-Steel prior,
  \code{priornuVW} the Villa-Walker prior (in the particular case of location
  parameter equal to 0).
  In both cases the prior is discretized,
  i.e. only consider integer values for the degrees of freedom.
}
\usage{
priornuJS(nu,k=2.78,numax=30)

priornuVW(nu,numax=30)
}
\arguments{
  \item{nu}{Values of the degrees of freedom}

  \item{k}{Parameter of the Juazer-Steel prior, see details}

  \item{numax}{Maximum value of nu, i.e. nu is restricted to take
  values between 1 and numax}
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
Vector with prior probabilities evaluated at nu
}
\references{
M.A. Juarez and M.F.J. Steel. Model-based clustering of non-Gaussian
panel data based on skew-t distributions.
Journal of Business and Economic Statistics, 28(1):52-66, 2010

C. Villa and S. Walker. Objective prior for the number of degrees
of freedom of a t distribution. Bayesian analysis, 9(1):197-200, 2014.
}

\examples{
library(twopiece)

nuseq <- 1:30
plot(nuseq, priornuJS(nuseq,k=1.2,numax=30), ylab='Prior density', xlab=expression(nu))
points(nuseq, priornuVW(nuseq,numax=30), col='grey', pch=16)
legend('topright',c('Juarez-Steel (k=1.2)','Villa-Walker'),pch=c(1,16),col=c('black','grey'))
}
\keyword{ distribution }
