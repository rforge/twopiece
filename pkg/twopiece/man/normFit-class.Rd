\name{normFit-class}
\Rdversion{1.1}
\docType{class}
\alias{normFit}
\alias{normFit-class}
\alias{coef,normFit-method}
\alias{show,normFit-method}
\alias{[,normFit,ANY,ANY,ANY-method}
\alias{[[,normFit-method}

\title{Class "normFit"}
\description{
Stores the Markov Chain Monte Carlo output from fitting a mixture of
multivariate skew-t distributions.
This is essentially a list with a few methods that were added
for convenience (see below).
}
\section{Objects from the Class}{
Objects from this class are automatically created
as the output of mixnormGibbs.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Single slot containing all the information. Since
    the class is built on a list one can access as usual its elements
    (mu, Sigma, alpha, nu, probs, probcluster, cluster, G, ttype).
    See details.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "normFit")}: Displays
      general information about the object. }

    \item{"["}{Selects an element of the list}

    \item{"[["}{Selects an element of the list}

    \item{coef}{\code{signature(object = "normFit")}:
    Computes the posterior mean of all model parameters}

    \item{clusterprobs}{\code{signature(fit = "normFit")}:
    Computes posterior cluster probabilities on a specified grid of value.}
  }
}
\details{

Below is a description of what each element in the list contains:

    - mu: List with G elements. Each contains Gibbs draws for the location
    parameter of each mixture component.

    - Sigma: List with G elements. Each contains Gibbs draws for the scale
    matrix of each mixture component. Only diagonal and upper diagonal elements
    are stored, since the matrix is symmetric.

    - alpha: List with G elements. Each contains Gibbs draws for the asymmetry
    parameters of each mixture component. Asymmetry parameters are in the epsilon-skew
    parameterization, see help(dtp3) for details.

    - nu: Matrix with G columns, each row contains Gibbs draws for the degrees
    of freedom of each mixture component.

    - probs: Matrix with G columns, each row contains Gibbs draws for the
    mixture weights.

    - probcluster: Matrix with the estimated posterior probabilities that
    each observation belongs to each cluster (Rao-Blackwellized)

    - cluster: if mixnormGibbs was called with \code{returnCluster=TRUE}
    cluster contains the latent cluster allocations at each Gibbs iteration
    (iterations in rows, observations in columns). Else it contains NA.

    - G: Assumed number of mixture components

    - ttype: Value of argument \code{ttype} passed to \code{mixnormGibbs},
    'independent' indicates that the iskew-t was used, 'dependent' the dskew-t.
}
\references{
Rossell D., Steel M.F.J. Continuous mixtures with skewness and heavy
tails. Handbook of Mixture Analysis, Chapter 10, CRC Press 2018.
}
\author{
David Rossell
}
\seealso{
See also \code{\link{mixnormGibbs}}.
}
\examples{
showClass("normFit")
}
\keyword{classes}
