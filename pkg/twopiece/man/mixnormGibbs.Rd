\name{mixnormGibbs}
\alias{mixnormGibbs}
\title{Gibbs sampling for multivariate Normal mixture}
\description{
Gibbs sampling to obtain posterior samples for
a multivariate Normal mixture
}
\usage{
mixnormGibbs(x, G, clusini='kmedians', priorParam=normprior(ncol(x),G),
niter, burnin=round(niter/10), returnCluster=FALSE, relabel='ECR',
verbose=TRUE)
}
\arguments{
\item{x}{Observed data matrix (observations in rows, variables in columns)}

\item{G}{Assumed number of components}

\item{clusini}{It can be a vector with initial cluster allocations, alternatively
set to 'kmedians' to use K-medians (as in cclust from package flexclust),
'kmeans' for K-means, or 'em' for EM algorithm-based on mixture of normals.
The latter is based on Mclust, which is initialized by hierarchical clustering
and can hence be slow}

\item{priorParam}{List with named elements 'mu', 'Sigma', 'alpha', 'nu', 'probs'
containing prior parameters. See help(normprior)}

\item{niter}{Number of Gibbs iterations}

\item{burnin}{Number of burn-in iterations (i.e. initial iterations to be discarded)}

\item{returnCluster}{If set to TRUE the allocated cluster at each MCMC
iteration is returned. Note that this can be memory-consuming if nrow(x) is large}

\item{relabel}{Cluster relabelling algorithm. Set to 'none' to do no relabelling.
'ECR' corresponds to the Papastimouis-Iliopoulos 2010 algorithm to make simulated
cluster allocations similar to a pivot (taken to be the last MCMC iteration).
Set to 'RW' for Rodriguez-Walker's (2014) proposal based on a loss function
to preserve cluster means across iterations.
Set to 'PC' for an artificial identifiability constraint based on projection
of location parameters on first principal component.
See help(label.switching) for details.}

\item{verbose}{Set to TRUE to output iteration progress}
}
\value{
Object of class normFit, i.e. essentially a list with elements
mu, Sigma and probs storing the posterior draws.
Posterior means can be extracted with coef() and cluster probabilities
with clusterprobs().
}
\seealso{help("normFit-class"), help("clusterprobs")}
\references{
Fuquene J., Steel M.F.J., Rossell D. On choosing mixture components via
non-local priors. arxiv:1604.00314, 2018

Rossell D., Steel M.F.J. Continuous mixtures with skewness and heavy
tails. Handbook of Mixture Analysis, Chapter 10, CRC Press 2018.
}
\examples{

library(twopiece)
set.seed(1)
n <- 1000; probs <- c(2/3,1/3)
mu1 <- c(0,0); S1 <- matrix(c(1,0,0,1),nrow=2)
mu2 <- c(3,3); S2 <- matrix(c(1,.5,.5,1),nrow=2)
mu <- list(mu1,mu2)
Sigma <- list(S1,S2)
xsim <- rmixnorm(n,mu=mu,Sigma=Sigma,probs=probs)
fit <- mixnormGibbs(xsim$x,G=2,clusini='kmedians',niter=100)

#Posterior means
coef(fit)

#Trace plots
plot(fit$mu[[1]][,1],type='l',xlab='Iteration',ylab='',main='Cluster 1')
lines(fit$mu[[1]][,2], col=2)
abline(h=c(mu1[1],mu1[2]),col=1:2,lty=2)
legend('topright',c(expression(mu[11]),expression(mu[12])),col=1:2,lty=1)

#Cluster probabilities
x1seq <- seq(-3,8,length=20); x2seq <- seq(-4,7,length=20)
xgrid <- expand.grid(x1seq,x2seq)
probgrid <- clusterprobs(fit, x=xgrid)
filled.contour(x=x1seq,y=x2seq,z=matrix(probgrid[,1],nrow=length(x1seq)),
plot.axes={ axis(1); axis(2); points(xsim$x,pch=xsim$cluster)})
}
\keyword{ cluster }
\keyword{ distribution }
\keyword{ model }