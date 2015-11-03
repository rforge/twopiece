\encoding{latin1}
\name{dtp4}
\alias{dtp4}
\alias{ptp4}
\alias{qtp4}
\alias{rtp4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The 4-Parameter Two Piece Distribution
}
\description{
Density, distribution function, quantile function and random generation for the 4-parameter two piece distribution with 3 parameterizations: two-piece (tp), epsilon-skew (eps), and inverse scale factors (isf).
}
\usage{
dtp4(x, mu, par1, par2, delta, FUN, param = "tp", log = FALSE)
ptp4(x, mu, par1, par2, delta, FUN, param = "tp", log.p = FALSE)
qtp4(p, mu, par1, par2, delta, FUN, param = "tp")
rtp4(n, mu, par1, par2, delta, FUN, param = "tp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of quantiles.
  }
  \item{p}{ vector of probabilities.
  }
  \item{n}{ number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
  }
  \item{mu}{ location parameter, \eqn{\mu}{%
  mu}.}
  \item{par1}{ scale parameter 1, \eqn{\sigma_1}{%
  par1}.}
  \item{par2}{ scale parameter 2, \eqn{\sigma_2}{%
  par2}.}
  \item{delta}{ shape parameter, \eqn{\delta}{%
  delta}.}
  \item{FUN}{ a symmetric density \code{f}.}
  \item{param}{ parameterizations used.}
  \item{log, log.p}{ logical; if TRUE, probabilities p are given as log(p).}
}
\details{
The 4-parameter two piece distribution with parameters \eqn{\mu}{%
mu}, \eqn{\sigma_1}{%
par1}, \eqn{\sigma_2}{%
par2} and \eqn{\delta}{%
delta} has the following density:
\deqn{s(x) = \frac{2}{\sigma_1 + \sigma_2}f( ( x -\mu )/\sigma_1, \delta ) \quad for \quad x<\mu }{%
s(x) = ( 2/(par1 + par2) ) f( (x - mu)/par1, delta ) for x < mu}
and
\deqn{s(x) = \frac{2}{\sigma_1 + \sigma_2}f( ( x -\mu )/\sigma_2, \delta ) \quad for \quad x\geq \mu }{%
s(x) = ( 2/(par1 + par2) ) f( (x - mu)/par2, delta ) for x >= mu}


where \code{f(x,delta)} is a symmetric density about zero.

More details

If param is not specified, it assumes the default value of "tp". Information about the "eps" and "isf" parameterizations can be found in the References.
}
\value{
   \code{dtp4} gives the density, \code{ptp4} gives the distribution function, \code{qtp4} gives the quantile function and \code{rtp4} generates random deviates.

}

\author{
F. J. Rubio \email{fxrubio@gmail.com} and A. M. López \email{amonteslop@gmail.com}.
}
\references{
Arellano-Valle, R. B Gómez, H. W. and Quintana, F. A. (2005). Statistical inference for general class of
asymmetric distributions. \emph{Journal of Statistical Planning and Inference}, \bold{128}: 427-443.

Fernández, C. and Steel, M. F. J. (1998). On Bayesian modeling of fat tails and skewness. 
\emph{Journal of the American Statistical Asociation}, \bold{93}, 359-371.

Mudholkar, G. S. and Hutson, A. D. (2000). The epsilon-skew-normal distribution for analyzing
 near-normal data. \emph{Journal of Statistical Planning and Inference}, \bold{83}: 291-309.

Rubio, F. J. and Steel, M. F. J. (2014). Inference in Two-Piece Location-Scale models with 
Jeffreys Priors, with discussion. \emph{Bayesian Analysis}, \bold{9}: 1-22.
}  
\seealso{
\link{dnorm} for the normal distribution and \link{dt} for the Student t distribution. %%\code{\link{dnorm}}

\code{\link{dtp3}} for the 3-parameter two piece distribution.
}
\examples{

## 4-parameter two piece Student-t density with parameterization 'tp'
tempf = function(x) dtp4(x,0,3,1,4,dt,param="tp")
curve(tempf,-10,5)
         
## 4-parameter two piece Student-t distribution with parameterization 'tp'
tempf = function(x) ptp4(x,0,3,1,4,pt,param="tp")
curve(tempf,-10,5)

## random number generation for 4-parameter two piece Student-t distribution
## with parameterization 'tp'
sim <- rtp4(1000,0,1,1,10,rt)
hist(sim, probability=TRUE, xlim=c(-10,10),ylim=c(0,dt(0,4)))

## quantile function for the 4-parameter two piece Student-t distribution
## with parameterization 'tp'
qtp4(0.5, 0, 1, 1, 4, qt ,param = "tp")

}