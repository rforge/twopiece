\name{mixskewtGibbs}
\alias{mixskewtGibbs}
\title{MH-within-Gibbs sampling for mixture of multivariate skew t's}
\description{
Metropolis-Hastings within Gibbs algorithm to obtain posterior samples for
a mixture of multivariate skew-t's.
}
\usage{
mixskewtGibbs(x, G, clusini='kmedians', priorParam=skewtprior(ncol(x),G), niter, burnin=round(niter/10), ttype='independent', returnCluster=FALSE, verbose=TRUE)
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
containing prior parameters. See help(skewtprior)}

\item{niter}{Number of Gibbs iterations}

\item{burnin}{Number of burn-in iterations (i.e. initial iterations to be discarded)}

\item{ttype}{Set to 'independent' for iskew-t, to 'dependent' for dskew-t.
See help(dskewt) for details}

\item{returnCluster}{If set to TRUE the allocated cluster at each MCMC
iteration is returned. Note that this can be memory-consuming if nrow(x) is large}

\item{verbose}{Set to TRUE to output iteration progress}
}
\value{
Object of class skewtFit, i.e. essentially a list with elements
mu, Sigma, alpha, nu and probs storing the posterior draws.
Posterior means can be extracted with coef() and cluster probabilities
with clusterprobs().
}
\seealso{help("skewtFit-class"), help("clusterprobs")}
\references{
Rossell D., Steel M.F.J. Continuous non-Gassian mixtures. In Handbook of Cluster Analysis, CRC Press.
}
\examples{

library(twopiece)
set.seed(1)
n <- 1000; probs <- c(2/3,1/3)
mu1 <- c(0,0); S1 <- matrix(c(1,0,0,1),nrow=2)
alpha1 <- c(0,0); nu1 <- 100
mu2 <- c(3,3); S2 <- matrix(c(1,.5,.5,1),nrow=2)
alpha2 <- c(-.5,0); nu2 <- 100
mu <- list(mu1,mu2)
Sigma <- list(S1,S2)
alpha <- list(alpha1,alpha2)
nu <- list(nu1,nu2)
xsim <- rmixskewt(n,mu=mu,Sigma=Sigma,alpha=alpha,nu=nu,probs=probs,param='eps')
fit <- mixskewtGibbs(xsim$x,G=2,clusini='kmedians',niter=100)

#Posterior means
coef(fit)

#Trace plots
plot(fit$mu[[1]][,1],type='l',xlab='Iteration',ylab='',main='Cluster 1',ylim=c(-1,1),cex.lab=1.25,cex.axis=1.2)
lines(fit$mu[[1]][,2], col=2)
abline(h=c(mu1[1],mu1[2]),col=1:2,lty=2)
legend('topright',c(expression(mu[11]),expression(mu[12])),col=1:2,lty=1,cex=1.25)

#Cluster probabilities
x1seq <- seq(-3,8,length=20); x2seq <- seq(-4,7,length=20)
xgrid <- expand.grid(x1seq,x2seq)
probgrid <- clusterprobs(fit, x=xgrid)
filled.contour(x=x1seq,y=x2seq,z=matrix(probgrid[,1],nrow=length(x1seq)),plot.axes={ axis(1); axis(2); points(xsim$x,pch=xsim$cluster)})
}
\keyword{ cluster }
\keyword{ distribution }
\keyword{ model }