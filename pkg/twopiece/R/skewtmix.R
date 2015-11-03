################################################################################################################
## METHODS FOR OBJECTS OF CLASS skewtFit
################################################################################################################

setMethod("show", signature(object="skewtFit"), function(object) {
  cat("skew-t mixture with ",object$G," components\n")
  cat("Use coef() to obtain posterior means. clusterprobs() for cluster probabilities\n")
}
)

setMethod("coef", signature(object="skewtFit"), function(object, ...) {
  mu <- lapply(object$mu,colMeans)
  nondiag <- (object$G+1):ncol(object$Sigma[[1]])
  Sigma <- lapply(object$Sigma, function(z) { colm <- colMeans(z); S <- diag(colm[1:object$G]); S[upper.tri(S)] <- S[lower.tri(S)] <- colm[nondiag]; return(S) })
  alpha <- lapply(object$alpha,colMeans)
  nu <- colMeans(object$nu)
  probs <- colMeans(object$probs)
  return(list(mu=mu,Sigma=Sigma,alpha=alpha,nu=nu,probs=probs))
}
)


#Compute cluster probabilities under skew-t mixture, averaging over posterior draws contained in x
# - fit: object of type clusterprobs
# - x: data points at which to evaluate the probabilities
# - iter: iterations of posterior draws over which to average. Defaults to using 1,000 equally spaced iterations.
setMethod("clusterprobs", signature(fit='skewtFit'), function(fit, x, iter) {
    if (missing(iter)) iter <- as.integer(seq(1,nrow(fit$probs),length=1000))
    if (!is.integer(iter)) stop("iter must contain integer iteration numbers")
    if ((min(iter)<1) | max(iter)>nrow(fit$mu[[1]])) stop("Specified iter values are out of valid range")
    if (missing(x)) stop("x values must be provided")
    G <- fit$G; p <- ncol(fit$mu[[1]])
    diag <- 1:G; nondiag <- (G+1):ncol(fit$Sigma[[1]])
    S <- lapply(1:G, function(i) { ans <- apply(fit$Sigma[[i]][iter,], 1, vec2matrix, diag=diag, nondiag=nondiag); array(as.vector(ans),dim=c(p,p,ncol(ans))) } )
    proboneiter <- function(i) {
        dx <- sapply(1:G,function(g) dskewt(x,mu=fit$mu[[g]][i,],Sigma=S[[g]][,,i],alpha=fit$alpha[[g]][i,],nu=fit$nu[i,g],param='eps',logscale=TRUE,ttype=fit$ttype))
        dx <- t(t(dx)+log(fit$probs[i,]))
        dx <- exp(dx - apply(dx,1,max))
        dx/rowSums(dx)
    }
    ans <- lapply(iter,function(i) proboneiter(i))
    ans <- Reduce('+',ans)/length(ans)
    return(ans)
}
)

vec2matrix <- function(vec,diag,nondiag) {
    #Format input vector as symmetric matrix. Indices of diagonal elements are in diag, those on the upper triangle in nondiag
    S <- diag(vec[diag])
    S[upper.tri(S)] <- S[lower.tri(S)] <- vec[nondiag]
    return(S)
}


################################################################################################################
## ROUTINES TO EVALUATE DENSITY, RANDOM NUMBER GENERATION
################################################################################################################

#Univariate skew T density
# When param=='eps': dt(x;mu,s/(1-alpha),nu) I(x>=mu) + dt(x;mu,s/(1+alpha),nu) I(x<mu)
# When param=='isf': dt(x;mu,s*alpha,nu) I(x>=mu) + dt(x;mu,s/alpha,nu) I(x<mu)
dskewtuniv <- function(x,mu,s,alpha,nu,param='eps',logscale=FALSE) {
    if (param=='eps') {
        ans <- dtp3(x,mu=mu,par1=s,par2=alpha,FUN=function(z,log) dt(z,df=nu,log=log),param="eps",log=logscale)
    } else if (param=='isf') {
        ans <- dtp3(x,mu=mu,par1=s,par2=1/alpha,FUN=function(z,log) dt(z,df=nu,log=log),param="isf",log=logscale)
    } else stop("Invalid value was specified for 'param'")
    return(ans)
}



check.skewtpars <- function(mu,Sigma,A,D,alpha,nu,param) {
    #Check that input parameters for multivariate t are correct
    p <- length(mu)
    if (missing(A) | missing(D)) {
        if (ncol(Sigma) != p) stop("Dimension of mu and Sigma doesn't match")
    } else {
        if ((ncol(A)!=p) | (ncol(D)!=p)) stop("Dimensions of mu and A,D don't match")
    }
    if (param=='eps') {
        if (any(abs(alpha)>1)) stop("When param=='eps' alpha must have entries in [-1,1]")
    } else if (param=='isf') {
        if (any(alpha<=0)) stop("When param=='isf' alpha must be a vector with entries >0")
    }
    if (length(nu)>1) stop("nu must have length 1")
    if (nu<=0) stop("nu must be >=0")
}

#Multivariate skew T density (p dimensions)
# Input
# - x: vector of length(p), alternatively matrix or data.frame with arbitrary number of rows and ncol(x)==p
# - mu: mean vector of length(p)
# - Sigma: p*p positive-definite matrix
# - A: optional parameter containing pre-computed A=t(eigen(Sigma)). For identifiability all elements in A[,1] must be positive, if any entry is negative then the corresponding row is changed sign
# - D: optional parameter containing pre-computed A=diag(eigen(Sigma)$values)
# - alpha: asymmetry parameters. Either a vector of length p or length 1, in the latter case all asymmetry parameters are assumed equal
# - nu: scalar indicating degrees of freedom
# - param: set to 'eps' for the epsilon-skew parametization, 'isf' for the inverse scale factor parameterization
# - logscale: set to TRUE to obtain log-pdf
# - ttype: set to 'independent' for independent skew-t, to 'dependent' for dependent skew-t
# Output: multivariate skew T pdf evaluated at x. It has length 1 if x was a vector, length equal to nrow(x) if x was a matrix or data.frame
dskewt <- function(x,mu,Sigma,A,D,alpha,nu,param='eps',logscale=FALSE,ttype='independent') {
    p <- length(mu)
    if (length(alpha) != p) { if (length(alpha)==1) alpha <- rep(alpha,p) else stop("length(alpha) must be either 1 or length(mu)") }
    check.skewtpars(mu=mu,Sigma=Sigma,A=A,D=D,alpha=alpha,nu=nu,param=param)
    #Eigendecomposition
    if (missing(A) | missing(D)) {
        e <- eigen(Sigma,symmetric=TRUE); A <- t(e$vectors); D <- diag(e$values); A <- A * ifelse(sign(A[,1])== -1, -1, 1)
    }
    if (any(diag(D)<=0)) stop("Sigma is not positive definite")
    if (is.vector(x)) {
        x <- as.numeric(x)
        if (length(x) != p) stop("Dimensions of x and mu don't match")
        z <- diag(1/sqrt(diag(D))) %*% A %*% matrix(x - as.vector(as.numeric(mu)),ncol=1)
    } else if (is.matrix(x) | is.data.frame(x)) {
        if (ncol(x) != p) stop("Dimensions of x and mu don't match")
        z <- diag(1/sqrt(diag(D))) %*% A %*% (t(as.matrix(x)) - as.vector(as.numeric(mu)))
    } else stop("x must be a vector, matrix or data.frame")
    if (ttype=='independent') {
        ans <-  -0.5*sum(log(diag(D)))
        for (i in 1:p) { ans <- ans + dskewtuniv(z[i,],mu=0,s=1,alpha=alpha[i],nu=nu,param=param,logscale=TRUE) }
    } else if (ttype=='dependent') {
        ans <-  -0.5*sum(log(diag(D))) + lgamma((p+nu)/2) - lgamma(nu/2) - 0.5*p*log(pi*nu)
        xi <- (z>=0) / (1-alpha)^2 + (z<0) / (1+alpha)^2
        ans <- ans - 0.5*(p+nu)*log(1 + rowSums(xi * z^2)/nu)
    } else { stop("Invalid ttype") }
    if (!logscale) ans <- exp(ans)
    return(ans)
}

#Density for mixture of multivariate skew t's (arguments as in dskewt)
dmixskewt <- function(x,mu,Sigma,alpha,nu,probs,param='eps',logscale=FALSE,ttype='independent') {
    K <- length(probs)
    if (length(mu) != K) stop("mu must be a list with length = length(probs)")
    if (length(Sigma) != K) stop("Sigma must be a list with length = length(probs)")
    if (length(alpha) != K) stop("alpha must be a list with length = length(probs)")
    if (is.vector(x)) x <- matrix(x,nrow=1)
    logd <- matrix(NA,nrow=nrow(x),ncol=K)
    for (i in 1:K) { logd[,i] <- dskewt(x,mu=mu[[i]],Sigma=Sigma[[i]],alpha=alpha[[i]],nu=nu[[i]],param=param,logscale=TRUE,ttype=ttype) + log(probs[i]) }
    ct <- apply(logd,1,max)
    ans <- exp(ct) * rowSums(exp(logd-ct))
    if (logscale) ans <- log(ans)
    return(ans)
}



#Draws from multivariate skew-t. Set ttype='independent' for independent skew-t, ttype='dependent for dependent skew-t
rskewt <- function(n,mu,Sigma,alpha,nu,param='eps',ttype='independent') {
    p <- length(mu)
    check.skewtpars(mu=mu,Sigma=Sigma,alpha=alpha,nu=nu,param=param)
    if (length(alpha) != p) { if (length(alpha)==1) alpha <- rep(alpha,p) else stop("length(alpha) must be either 1 or length(mu)") }
    #Eigendecomposition
    e <- eigen(Sigma,symmetric=TRUE); A <- t(e$vectors); D <- diag(e$values); A <- A * ifelse(sign(A[,1])== -1, -1, 1)
    if (any(e$values<=0)) stop("Sigma is not positive definite")
    #Draw from underlying uncorrelated variables
    ans <- matrix(NA,nrow=n,ncol=p)
    if (ttype=='independent') {
        if (param=='eps') {
            for (i in 1:p) ans[,i] <- rtp3(n,mu=0,par1=1,par2=alpha[i],FUN=function(n) rt(n=n,df=nu), param="eps")
        } else if (param=='isf') {
            for (i in 1:p) ans[,i] <- rtp3(n,mu=0,par1=1,par2=1/alpha[i],FUN=function(n) rt(n=n,df=nu), param="isf")
        } else stop("Invalid value was specified for 'param'")
    } else if (ttype=='dependent') {
        if (param=='eps') {
            omega <- 1/rgamma(n,nu/2,nu/2)
            for (i in 1:p) ans[,i] <- sqrt(omega) * rtp3(n,mu=0,par1=1,par2=alpha[i],FUN=function(n) rnorm(n=n), param="eps")
        } else if (param=='isf') {
            for (i in 1:p) ans[,i] <- sqrt(omega) * rtp3(n,mu=0,par1=1,par2=1/alpha[i],FUN=function(n) rnorm(n=n), param="isf")
        } else stop("Invalid value was specified for 'param'")
    } else {
        stop("Invalid ttype")
    }
    #Perform linear transform
    ans <- t(t(A) %*% (sqrt(e$values) * t(ans)) + mu)
    return(ans)
}



#Multivariate draws from mixture of skew-t's
rmixskewt <- function(n,mu,Sigma,alpha,nu,probs,param='eps',ttype='independent') {
    #Checks
    K <- length(probs)
    if (!is.list(mu)) { mu <- lapply(1:K,function(z) mu); warning("Formatted mu as list") } else { if (length(mu) != K) stop("mu has the wrong length") }
    if (!is.list(Sigma)) { Sigma <- lapply(1:K,function(z) Sigma); warning("Formatted Sigma as list") } else { if (length(Sigma) != K) stop("Sigma has the wrong length") }
    if (!is.list(alpha)) { alpha <- lapply(1:K,function(z) alpha); warning("Formatted alpha as list") } else { if (length(alpha) != K) stop("alpha has the wrong length") }
    p <- length(mu[[1]])
    for (i in 1:K) {
        check.skewtpars(mu=mu[[i]],Sigma=Sigma[[i]],alpha=alpha[[i]],nu=nu[[i]],param=param)
        if (length(alpha[[i]]) != p) { if (length(alpha[[i]])==1) alpha[[i]] <- rep(alpha[[i]],p) else stop("length(alpha[[i]]) must be either 1 or length(mu[[i]])") }
    }
    #Obtain draws
    ngroup <- as.vector(rmultinom(1, size=n, prob=probs))
    groupidx <- c(0,cumsum(ngroup))
    x <- matrix(NA,nrow=n,ncol=p)
    for (i in 1:K) { x[(groupidx[i]+1):groupidx[i+1],] <- rskewt(ngroup[i],mu=mu[[i]],Sigma=Sigma[[i]],alpha=alpha[[i]],nu=nu[[i]],param=param,ttype=ttype) }
    idx <- sample(1:nrow(x), size=nrow(x), replace=FALSE)
    x <- x[idx,]
    cluster <- rep(1:K,ngroup)[idx]
    return(list(x=x,cluster=cluster))
}


################################################################################################################
## ROUTINES TO EVALUATE PRIOR DISTRIBUTIONS
################################################################################################################

#Discretized Juarez-Steel gamma-gamma prior on degrees of freedom truncated to [1,numax]
priornuJS <- function(nu,k=1.2,numax=30) {
    if (any(nu>numax)) stop("All elements in nu must be <=numax")
    if (any((nu %% 1) != 0)) stop("nu cannot have non-integer values")
    if (any(nu<=0)) stop("All elements in nu must be strictly positive")
    ans <- k*(1:numax)/((1:numax)+k)^3
    ans <- ans/sum(ans)
    return(ans[as.integer(nu)])
}

#Villa-Walker objective prior on degrees of freedom truncated to [1,numax]
priornuVW <- function(nu,numax=30) {
    if (any(nu>numax)) stop("All elements in nu must be <=numax")
    if (any((nu %% 1) != 0)) stop("nu cannot have non-integer values")
    if (any(nu<=0)) stop("All elements in nu must be strictly positive")
    f1 <- function(x,nufix) { log(1+x^2/nufix) * dt(x,df=nufix) }
    f2 <- function(x,nufix) { log(1+x^2/(nufix+1)) * dt(x,df=nufix) }
    #f3 <- function(x,nufix) { log(1+x^2/(nufix-1)) * dt(x,df=nufix) }
    #f4 <- function(x,nufix) { log(1+x^2/nufix) * dnorm(x) }
    nuseq <- 1:numax
    ans <- 0.5*(log(nuseq+1)-log(nuseq)) + lbeta(0.5,0.5*(nuseq+1)) - lbeta(0.5,0.5*nuseq)
    for (i in 1:numax) {
        int1 <- integrate(f1,-Inf,Inf,nufix=i)$value
        int2 <- integrate(f2,-Inf,Inf,nufix=i)$value
        ans[i] <- ans[i] - 0.5*(i+1)*int1 + 0.5*(i+2)*int2
    }
    #int1 <- integrate(f1,-Inf,Inf,nufix=numax-1)$value
    #int2 <- integrate(f3,-Inf,Inf,nufix=numax-1)$value
    #ans[numax-1] <- ans[numax-1] - 0.5*numax*int1 + 0.5*(numax-1)*int2
    #int2 <- integrate(f4,-Inf,Inf,nufix=numax-1)$value
    #ans[numax] <- ans[numax] - 0.5 + 0.5*(numax-1)*int2
    ans <- exp(ans-max(ans))
    ans <- ans/sum(ans)
    return(ans[as.integer(nu)])
}




################################################################################################################
## AUXILIARY FUNCTIONS
################################################################################################################

rdirichlet <- function (n, alpha) {
    #Obtain n draws from Dirichlet(alpha)
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}

diwishfast <- function (eigenvecW, eigenvalW, v, S, logscale=TRUE) {
    #Density of inverse Wishart ignoring constant terms depending on degrees of freedom (v), i.e. gives same output proportional to dinvwishart (below)
    k <- nrow(S)
    logdetS <- as.numeric(determinant(S,logarithm=TRUE)$modulus)
    logdetW <- sum(log(eigenvalW))
    hold <- S %*% t(eigenvecW) %*% diag(1/eigenvalW) %*% eigenvecW    #S %*% solve(W)
    ans <- 0.5*v*logdetS - 0.5*(v+k+1)*logdetW  - 0.5*sum(diag(hold))
    if (!logscale) ans <- exp(ans)
    return(ans)
}


dinvwishart <- function (Sigma, nu, S, log = FALSE) {
    k <- nrow(Sigma)
    gamsum <- 0
    for (i in 1:k) gamsum <- gamsum + lgamma((nu + 1 - i)/2)
    dens <- -(nu * k/2) * log(2) - ((k * (k - 1))/4) * log(pi) - gamsum + (nu/2) * log(det(S)) - ((nu + k + 1)/2) * log(det(Sigma)) - 0.5 * sum(diag(S %*% solve(Sigma)))
    if (log == FALSE) dens <- exp(dens)
    return(dens)
}

rinvwishart <- function(nu, S, Sinv) {
    if (missing(Sinv)) Sinv <- solve(S)
    solve(rWishart(1, df=nu, Sigma=Sinv)[,,1])
}


################################################################################################################
## GIBBS SAMPLING
################################################################################################################

#Create list containing prior parameters for mixture of multivariate skew-t's. The following prior functional form is used:
#  (mu,Sigma) ~ N(mu; m,g*Sigma) * IWishart(Sigma; Q, q)
#  0.5*(alpha+1) ~ Beta(a,b)
#  nu is discrete with P(nu=j)=nuprobs[j] (or proportional to nuprobs[j])
#  probs ~ Symmetric Dirichlet(r), where r is a scalar
# Input
# - p: dimension of the observed data
# - m, g: parameters for prior on mu given Sigma
# - Q, q: parameters for prior on Sigma
# - a, b: parameter for prior on asymetry parameter
# - r: parameter for symmetric Dirichlet prior on mixing probabilities
# - nuprobs: vector where P(nu=j)=nuprobs[j]. Names can be provided, else it is assumed that support of nu ranges from 1 to length(nuprobs)
skewtprior <- function(p,m=rep(0,p),g=1,Q=diag(p),q=p+1,a=2,b=2,r=1/p,nuprobs=priornuJS(1:30,k=1.2,numax=30)) {
    if (length(m) != p) stop("m has the wrong length")
    if (length(g)>1) stop("g must have length 1")
    if (!is.matrix(Q)) stop("Q must be a matrix")
    if ((nrow(Q) != ncol(Q)) | (nrow(Q) != p)) stop("Q must be a square matrix with p rows")
    if (any(Q[upper.tri(Q)] != Q[lower.tri(Q)])) stop("Q must be symmetric")
    if (any(eigen(Q,symmetric=TRUE)$values <= 0)) stop("Q is not positive-definite!")
    if (q < p) stop("q must be >= p")
    if (length(r)>1) stop("r must have length 1")
    ans <- vector("list",5)
    names(ans) <- c('mu','Sigma','alpha','nu','probs')
    ans[['mu']] <- list(m=m,g=g)
    ans[['Sigma']] <- list(Q=Q,q=q)
    ans[['alpha']] <- list(a=a,b=b)
    ans[['nu']] <- nuprobs; if (is.null(names(nuprobs))) names(ans[['nu']]) <- 1:length(nuprobs)
    ans[['probs']] <- c(r=r)
    return(ans)
}



#MH-within-Gibbs sampling for mixture of multivariate skew t's
# Input
# - x: observed data (observations in rows, variables in columns)
# - G: number of components
# - clusini: either vector with initial cluster allocation, 'kmedians' for K-medians (as in cclust from package flexclust) 'kmeans' for K-means, 'em' for EM algorithm based on mixture of normals (as in Mclust package, initialized by hierarchical clustering so it can be slow)
# - priorParam: list with named elements 'mu', 'Sigma', 'alpha', 'nu', 'probs' containing prior parameters. See help(skewtprior).
# - niter: number of Gibbs iterations
# - ttype: set to 'independent' for independent skew-t, to 'dependent' for dependent skew-t
# - returnCluster: if set to TRUE the allocated cluster at each MCMC iteration is returned. This can be memory-consuming if nrow(x) is large.
# - verbose: set to TRUE to output iteration progress
# Output:
# - mu: list of length G, where mu[[i]] is a matrix with posterior draws (niter-burnin rows)
# - Sigma: list of length G, where Sigma[[i]] is a matrix with posterior draws (niter-burnin rows). Each row contains diag(S),S[upper.tri(S)]
# - alpha: list of length G, where alpha[[i]] is a matrix with posterior draws (niter-burnin rows)
# - nu: matrix with niter-burnin rows and G columns with posterior draws for the degrees of freedom
# - probs: matrix with niter-burnin rows and G columns with posterior draws for the mixing probabilities
# - probcluster: matrix with nrow(x) rows and G columns with posterior probabilities that each observation belongs to each cluster (Rao-Blackwellized)
# - cluster: if returnCluster==TRUE, a matrix with niter-burnin rows and nrow(x) columns with latent cluster allocations at each MCMC iteration. If returnCluster==FALSE, NA is returned.
mixskewtGibbs <- function(x, G, clusini='kmedians', priorParam=skewtprior(ncol(x)), niter, burnin=round(niter/10), ttype='independent', returnCluster=FALSE, verbose=TRUE) {
    #require(mvtnorm)
    if (!(ttype %in% c('independent','dependent'))) stop("ttype must be equal to 'independent' or 'dependent'")
    p <- ncol(x); n <- nrow(x)
    m <- priorParam[['mu']]$m; gprior <- priorParam[['mu']]$g
    Q <- priorParam[['Sigma']]$Q; q <- priorParam[['Sigma']]$q
    a <- priorParam[['alpha']]$a; b <- priorParam[['alpha']]$b
    nuprobs <- priorParam[['nu']]
    if (is.null(names(nuprobs))) { nuvals <- 1:length(nuprobs) } else { nuvals <- as.integer(names(nuprobs)) }
    r <- priorParam[['probs']]['r']
    ## Initial cluster allocations ##
    if (is.character(clusini)) {
        if (clusini=='kmedians') {
            #require(flexclust)
            z <- predict(cclust(x, k=G, dist='manhattan', method='kmeans'))
        } else if (clusini=='kmeans') {
            z <- kmeans(x, centers=G)$cluster
        } else if (clusini=='em') {
            #require(mclust)
            z <- Mclust(x, G=G, modelNames='VVV')$classification
        } else stop("Invalid value for 'clusini'")
    } else {
        if (length(clusini) != nrow(x)) stop("length(clusini) must be equal to nrow(x)")
        z <- as.integer(factor(clusini))
    }
    zcount <- rep(0,G); names(zcount) <- (1:G)
    tab <- table(z); zcount[names(tab)] <- tab
    ## Initialize parameters ##
    probscur <- as.vector(rdirichlet(1,r+zcount))
    txstd <- matrix(NA,nrow=ncol(x),ncol=nrow(x))
    tx <- t(x); xstd <- t(txstd)
    xi <- matrix(1,nrow=nrow(x),ncol=ncol(x))
    mucur <- alphacur <- lapply(1:G, function(i) double(p))
    Scur <- Acur <- Dcur <- lapply(1:G, function(i) matrix(NA,nrow=p,ncol=p)); nucur <- rep(30,G)
    if (ttype=='independent') wcur <- matrix(1,nrow=n,ncol=p) else wcur <- matrix(1,nrow=n,ncol=1)
    for (i in 1:G) {
        sel <- which(z==i)
        #Initialize alpha
        alphacur[[i]] <- rep(0,p)
        #Initialize mu, Sigma
        if (length(sel)>0) {
            mucur[[i]] <- apply(x[sel,],2,median); Scur[[i]] <- ((tx[,sel] - mucur[[i]]) %*% t(tx[,sel] - mucur[[i]]))/length(sel)
            #mucur[[i]] <- colMeans(x[sel,]); Scur[[i]] <- var(x[sel,]) * (nucur[[i]]-2) / nucur[[i]]
        } else {
            mucur[[i]] <- m; Scur[[i]] <- Q / (q+p+1)
        }
        #Initialize w
        e <- eigen(Scur[[i]],symmetric=TRUE); Acur[[i]] <- t(e$vectors); Dcur[[i]] <- diag(e$values); Acur[[i]] <- Acur[[i]] * ifelse(sign(Acur[[i]][,1])== -1, -1, 1)
        if (length(sel)>0) {
            txstd[,sel] <- sqrt(Dcur[[i]]) %*% Acur[[i]] %*% (tx[,sel] - mucur[[i]])
            xstd[sel,] <- t(txstd[,sel])
            xi[sel,] <- t((txstd[,sel]>=0) / (1-alphacur[[i]])^2 + (txstd[,sel]<0) / (1+alphacur[[i]])^2)
            wcur[sel,] <- rwSkewtGibbs(xstd=xstd[sel,],xi=xi[sel,],nu=nucur[[i]],ttype=ttype)
        }
    }
    ##Gibbs sampling
    mu <- lapply(1:G,function(i) matrix(NA,nrow=niter-burnin,ncol=p))
    Sigma <- lapply(1:G,function(i) matrix(NA,nrow=niter-burnin,ncol=p*(p+1)/2))
    alpha <- lapply(1:G,function(i) matrix(NA,nrow=niter-burnin,ncol=p))
    nu <- probs <- matrix(NA,nrow=niter-burnin,ncol=G)
    probclus <- matrix(0,nrow=nrow(x),ncol=G); rownames(probclus) <- paste('indiv',1:n,sep='')
    if (returnCluster) { cluster <- matrix(NA,nrow=niter-burnin,ncol=nrow(x)); colnames(cluster) <- rownames(probclus) }
    if (verbose) { niter10 <- round(niter/10); cat("Running MCMC") }
    for (l in 1:niter) {
        #Sample latent cluster indicators z
        dx <- sapply(1:G,function(g) dskewt(x,mu=mucur[[g]],A=Acur[[g]],D=Dcur[[g]],alpha=alphacur[[g]],nu=nucur[[g]],param='eps',logscale=TRUE,ttype=ttype))
        dx <- t(t(dx)+log(probscur))
        dx <- exp(dx - apply(dx,1,max))
        dx <- dx/rowSums(dx)
        z <- apply(dx,1,function(pp) match(TRUE,runif(1)<cumsum(pp)))
        tab <- table(z); zcount[1:G] <- 0; zcount[names(tab)] <- tab
        #Sample component weights
        probscur <- as.vector(rdirichlet(1,alpha=r+zcount))
        for (g in 1:G) {
            sel <- which(z==g)
            if (length(sel)>0) {
                #Update precomputed stuff
                txstd[,sel] <- sqrt(Dcur[[g]]) %*% Acur[[g]] %*% (tx[,sel] - mucur[[g]])
                xstd[sel,] <- t(txstd[,sel])
                xi[sel,] <- t((txstd[,sel]>=0) / (1-alphacur[[g]])^2 + (txstd[,sel]<0) / (1+alphacur[[g]])^2)
                #Sample w
                wcur[sel,] <- rwSkewtGibbs(xstd=xstd[sel,],xi=xi[sel,],nu=nucur[[g]],ttype=ttype)
                #Sample mu, Sigma
                munew <- rmuSkewtGibbs(tx[,sel],mucur=mucur[[g]],Sigma=Scur[[g]],A=Acur[[g]],D=Dcur[[g]],w=wcur[sel,],alpha=alphacur[[g]],xi=xi[sel,],nu=nucur[[g]],ttype=ttype,m=m,g=gprior)
                if (munew$accept) {
                    mucur[[g]] <- munew$mu
                    txstd[,sel] <- sqrt(Dcur[[g]]) %*% Acur[[g]] %*% (tx[,sel] - mucur[[g]])
                    xi[sel,] <- t((txstd[,sel]>=0) / (1-alphacur[[g]])^2 + (txstd[,sel]<0) / (1+alphacur[[g]])^2)
                }
                newS <- rSigmaSkewtGibbs(tx[,sel],mu=mucur[[g]],Sigmacur=Scur[[g]],Acur=Acur[[g]],Dcur=Dcur[[g]],w=wcur[sel,],alpha=alphacur[[g]],xi=xi[sel,],ttype=ttype,Q=Q,q=q)
                if (newS$accept) {
                    Scur[[g]] <- newS$Sigmanew; Acur[[g]] <- newS$Anew; Dcur[[g]] <- newS$Dnew
                    txstd[,sel] <- sqrt(Dcur[[g]]) %*% Acur[[g]] %*% (tx[,sel] - mucur[[g]])
                    #xi[sel,] <- t((txstd[,sel]>=0) / (1-alphacur[[g]])^2 + (txstd[,sel]<0) / (1+alphacur[[g]])^2)  #xi not needed for alpha or nu
                }
                if (munew$accept | newS$accept) xstd[sel,] <- t(txstd[,sel])
                #Sample alpha
                alphanew <- ralphaSkewtGibbs(alphacur=alphacur[[g]],xstd=xstd[sel,],w=wcur[sel,],a=a,b=b,ttype=ttype)
                alphacur[[g]] <- alphanew$alphanew
                #xi[sel,] <- t((txstd[,sel]>=0) / (1-alphacur[[g]])^2 + (txstd[,sel]<0) / (1+alphacur[[g]])^2) #xi not needed for nu
                #Sample nu
                nucur[[g]] <- rnuSkewtGibbs(x=x[sel,],mu=mucur[[g]],A=Acur[[g]],D=Dcur[[g]],alpha=alphacur[[g]],nuprobs=nuprobs,ttype=ttype)
            } else {
                #If no observations in cluster, sample from the prior
                Scur[[g]] <- rinvwishart(q,Q)
                e <- eigen(Scur[[g]],symmetric=TRUE); Acur[[g]] <- t(e$vectors); Dcur[[g]] <- diag(e$values); Acur[[g]] <- Acur[[g]] * ifelse(sign(Acur[[g]][,1])== -1, -1, 1)
                mucur[[g]] <- dmvnorm(1,m,gprior * Scur[[g]])
                alphacur[[g]] <- 1 - 2*rbeta(p,a,b)
                nucur[[g]] <- as.numeric(names(nuprobs)[rmultinom(n=1,size=1,prob=nuprobs)])
            }
        }
        #Save output
        if (l>burnin) {
            idx <- l-burnin
            for (g in 1:G) {
                mu[[g]][idx,] <- mucur[[g]]
                Sigma[[g]][idx,] <- c(diag(Scur[[g]]),Scur[[g]][upper.tri(Scur[[g]])])
                alpha[[g]][idx,] <- alphacur[[g]]
            }
            nu[idx,] <- nucur
            probs[idx,] <- probscur
            probclus <- probclus + dx
            if (returnCluster) cluster[idx,] <- z
        }
        if (verbose & ((l %% niter10)==0)) cat('.')
    }
    if (verbose) cat('\n')
    names(mu) <- names(Sigma) <- names(alpha) <- colnames(nu) <- colnames(probs) <- colnames(probclus) <- paste('cluster',1:ncol(nu),sep='')
    probclus <- probclus/(niter-burnin)
    ans <- list(mu=mu,Sigma=Sigma,alpha=alpha,nu=nu,probs=probs,probcluster=probclus,cluster=NA,G=G,ttype=ttype)
    if (returnCluster) ans$cluster <- cluster
    return(new("skewtFit",ans))
}



rmuSkewtGibbs <- function(tx,mucur,Sigma,A,D,w,alpha,xi,nu,ttype,m,g) {
    #Metropolis-Hastings sampler from the conditional posterior of the skew-t location parameter mu given all other parameters
    #Input
    # - tx: t(x), where x are original observations
    # - mucur: current value of mu
    # - Sigma: scale matrix
    # - A: pre-computed t(eigen(Sigma)$vectors). For identifiability all elements in A[,1] must be positive, if any entry is negative then the corresponding row is changed sign
    # - D: pre-computed diag(eigen(Sigma)$values)
    # - w: skew-t latent scale parameter (it can either be a vector for dependent t or a matrix for dependent t)
    # - alpha: vector of length ncol(x) with asymmetry parameters
    # - xi: (x>=0)/(1-alpha)^2 + (x<0)/(1+alpha)^2
    # - nu: degrees of freedom
    # - ttype: independent or dependent skew-t
    # - m,g: prior on mu is N(m,g*Sigma)
    #Output: list with 2 components
    # - munew: sampled value of mu (it may be equal to mucur)
    # - accept: indicates if munew==mucur
    #
    #Pre-compute useful quantities
    sqDA <- diag(1/sqrt(diag(D))) %*% A
    Sigmainv <- solve(Sigma)
    tphi <- t(xi/w)
    xiwinv <- diag(rowSums(tphi))
    #Parameters of proposal
    Vinv <- Sigmainv/g + t(sqDA) %*% xiwinv %*% sqDA
    term1 <- t(sqDA) %*% (tphi * (sqDA %*% tx)) %*% matrix(1,nrow=ncol(tx))  #possible problem here
    term2 <- g * Sigma %*% m
    V <- solve(Vinv)
    muprop <- V %*% (term1 + term2)
    munew <- as.vector(rmvnorm(1,muprop,sigma=V))
    #Parameters of reverse proposal
    txstdnew <- sqDA %*% (tx - munew)
    xinew <- t((txstdnew>=0) / (1-alpha)^2 + (txstdnew<0) / (1+alpha)^2)
    tphinew <- t(xinew/w)
    xiwinvnew <- diag(rowSums(tphinew))
    Vinvnew <- Sigmainv/g + t(sqDA) %*% xiwinvnew %*% sqDA
    term1 <- t(sqDA) %*% (tphinew * (sqDA %*% tx)) %*% matrix(1,nrow=ncol(tx))
    Vnew <- solve(Vinvnew)
    mupropnew <- Vnew %*% (term1 + term2)
    #Acceptance probability (likelihood marginalized wrt w)
    u <- dmvnorm(mucur,mupropnew,solve(Vinvnew),log=TRUE) - dmvnorm(munew,muprop,V,log=TRUE)
    #num <- sum(dskewt(t(tx),mu=munew,Sigma=Sigma,A=A,D=D,alpha=alpha,nu=nu,param='eps',logscale=TRUE,ttype=ttype)) + dmvnorm(munew,m,g*Sigma,log=TRUE)
    #den <- sum(dskewt(t(tx),mu=mucur,Sigma=Sigma,A=A,D=D,alpha=alpha,nu=nu,param='eps',logscale=TRUE,ttype=ttype)) + dmvnorm(mucur,m,g*Sigma,log=TRUE)
    #Acceptance probability (conditional on w)
    num <- sqrt(tphinew) * (sqDA %*% (tx-munew))
    num <- -0.5 * sum(num^2) - 0.5 * (t(munew-m) %*% Sigmainv %*% (munew-m)) / g
    den <- sqrt(tphi) * (sqDA %*% (tx-mucur))
    den <- -0.5 * sum(den^2) - 0.5 * (t(mucur-m) %*% Sigmainv %*% (mucur-m)) / g
    u <- exp(u + num - den)
    if (runif(1)<u) { accept <- TRUE } else { accept <- FALSE; munew <- mucur }
    return(list(munew=munew,accept=accept))
}


rSigmaSkewtGibbs <- function(tx,mu,Sigmacur,Acur,Dcur,w,alpha,xi,ttype,Q=Q,q=q) {
    #MH draws from conditional posterior of Sigma= t(A) %*% D %*% A given all other parameters
    # Input
    # - mu: mean vector
    # - Sigmacur: current value of Sigma
    # - Acur, Dcur: current value of A,D, i.e. Sigmacur= t(Acur) %*% Dcur %*% Acur
    # - w: current value of latent scale parameters
    # - alpha: current value of asymmetry parameters
    # - xi: (x>=0)/(1-alpha)^2 + (x<0)/(1+alpha)^2
    # - ttype: independent or dependent skew-t
    # - Q, q: prior on Sigma is InvWishart(Q,q), where q are the degrees of freedom
    # Output: list with 2 components
    # - Sigmanew: sampled value of Sigma (it may be equal to Sigmacur)
    # - Anew: sampled value of A (it may be equal to Acur)
    # - Dnew: sampled value of D (it may be equal to Dcur)
    # - accept: indicates if Sigmanew==Sigmacur
    #
    n <- nrow(xi); p <- ncol(xi)
    sqphi <- sqrt(xi/w)
    txmu <- tx-mu
    Dcurinv <- diag(1/diag(Dcur))
    #Proposal parameters
    xortho <- t(Acur %*% txmu) * sqphi
    Sprop <- Q + matrix(mu,ncol=1) %*% matrix(mu,nrow=1) + t(Acur) %*% t(xortho) %*% xortho %*% Acur
    Sigmanew <- rinvwishart(n+q+p, Sprop)  #used to be riwish
    #Parameters for reverse proposal
    e <- eigen(Sigmanew,symmetric=TRUE); Anew <- t(e$vectors); Dnew <- diag(e$values); Dnewinv <- diag(1/e$values); Anew <- Anew * ifelse(sign(Anew[,1])== -1, -1, 1)
    xorthonew <- t(Anew %*% txmu) * sqphi
    Spropnew <- Q + matrix(mu,ncol=1) %*% matrix(mu,nrow=1) + t(Anew) %*% t(xorthonew) %*% xorthonew %*% Anew
    #Acceptance probability
    u <- diwishfast(Acur,diag(Dcur),n+q+p,Spropnew,logscale=TRUE) - diwishfast(Anew,diag(Dnew),n+q+p,Sprop,logscale=TRUE)
    #u <- dinvwishart(Sigmacur,n+q+p,Spropnew,log=TRUE) - dinvwishart(Sigmanew,n+q+p,Sprop,log=TRUE) #same but slower
    txorthonew <- t(sqphi) * (diag(sqrt(diag(Dnewinv))) %*% Anew %*% txmu)
    txorthocur <- t(sqphi) * (diag(sqrt(diag(Dcurinv))) %*% Acur %*% txmu)
    smunew <- matrix(mu,nrow=1) %*% t(Anew) %*% Dnewinv %*% Anew %*% matrix(mu,ncol=1)
    smucur <- matrix(mu,nrow=1) %*% t(Acur) %*% Dcurinv %*% Acur %*% matrix(mu,ncol=1)
    logposnew <- -0.5 * (sum(txorthonew^2) + smunew) - 0.5*n*sum(log(diag(Dnew))) + diwishfast(Anew,diag(Dnew),v=q+p,S=Q,logscale=TRUE)
    logposcur <- -0.5 * (sum(txorthocur^2) + smucur) - 0.5*n*sum(log(diag(Dcur))) + diwishfast(Acur,diag(Dcur),v=q+p,S=Q,logscale=TRUE)
    u <- u + logposnew - logposcur
    #Check: using normal pdf with sqphi=1 we get acceptance prob=1
    #logposnorm <- sum(dmvnorm(t(tx),mean=mu,sigma=Sigmanew,log=TRUE)) + diwishfast(t(eigen(Sigmanew)$vectors),eigen(Sigmanew)$values,v=q,S=Q,logscale=TRUE)
    #logposnorm2 <- sum(dmvnorm(t(tx),mean=mu,sigma=Sigmacur,log=TRUE)) + diwishfast(t(eigen(Sigmacur)$vectors),eigen(Sigmacur)$values,v=q,S=Q,logscale=TRUE)
    #u + logposnorm-logposnorm2
    if (runif(1)<exp(u)) {
        ans <- list(Sigmanew=Sigmanew,Anew=Anew,Dnew=Dnew,accept=TRUE)
    } else {
        ans <- list(Sigmanew=Sigmacur,Anew=Acur,Dnew=Dcur,accept=FALSE)
    }
    return(ans)
}


rwSkewtGibbs <- function(xstd,xi,nu,ttype) {
    #Sample from the conditional posterior of the skew-t latent weights w given all other parameters
    #Input
    # - xstd: matrix with standardized & uncorrelated observations, i.e. xstd= Sigma^(-1/2) (y-mu)
    # - xi: matrix with (xstd>=0)/(1-alpha)^2 + (xstd<0)/(1+alpha)^2, where alpha is the vector with asymmetry parameters
    # - nu: degrees of freedom
    # - ttype: independent skew-t or dependent skew-t
    #Output: if ttype=='independent' a matrix with dim(xstd), else a matrix with nrow(xstd) rows and 1 column
    if (ttype=='independent') {
        par1 <- (nu+1) / 2; par2 <- (nu + xstd^2) / 2
        ans <- matrix(1/rgamma(nrow(xstd)*ncol(xstd),shape=par1,rate=par2), nrow=nrow(xstd), ncol=ncol(xstd))
    } else {
        par1 <- (nu+ncol(xstd)) / 2; par2 <- (nu + rowSums(xstd^2)) / 2
        ans <- matrix(1/rgamma(nrow(xstd),shape=par1,rate=par2),ncol=1)
    }
    return(ans)
}


ralphaSkewtGibbs <- function(alphacur,xstd,w,a,b,ttype) {
    #Metropolis-Hastings sampler for vector of asymmetry parameters
    # Input
    # - xstd: standardized & orthogonalized x values, i.e. sqrt(D) A (x-mu)
    # - w: latent scale parameters
    # - a, b: prior on .5 * (alpha+1) ~ Beta(a,b)
    # - ttype: independent or dependent skew-t
    if (ttype=='independent') {
        s1 <- colSums((xstd>=0) * xstd^2 / w)
        s2 <- colSums((xstd< 0) * xstd^2 / w)
    } else {
        s1 <- colSums((xstd>=0) * xstd^2 / w)
        s2 <- colSums((xstd< 0) * xstd^2 / w)
    }
    alphaprop <- ralphaPropT(s1=s1,s2=s2,n=max(3,nrow(xstd)),a=a,b=b)
    u <- dt((alphacur-alphaprop$alpha.mode)/alphaprop$sd.mode,df=max(3,nrow(xstd)),log=TRUE) - alphaprop$logpdfProp
    #u <- dnorm(alphacur,alphaprop$alpha.mode,sd=alphaprop$sd.mode,log=TRUE) - alphaprop$logpdfProp
    u <- u + falpha(alphaprop$alpha,s1=s1,s2=s2,a=a,b=b,logscale=TRUE) - falpha(alphacur,s1=s1,s2=s2,a=a,b=b,logscale=TRUE)
    accept <- (runif(length(u)) < exp(u))
    alphanew <- alphaprop$alpha
    alphanew[!accept] <- alphacur[!accept]
    return(list(alphanew=alphanew,accept=accept))
}

#Full conditional posterior density of alpha
falpha <- function(alpha,s1,s2,a=1,b=1,logscale=FALSE) {
    ans <- -0.5*(s1/(1-alpha)^2 + s2/(1+alpha)^2) + (a-1)*log(1+alpha) + (b-1)*log(1-alpha)
    if (!logscale) ans <- exp(ans)
    return(ans)
}

#Density of Normal approximation to full posterior
falphaProp <- function(alpha,s1,s2,a=1,b=1,alpha.mode,alpha.sd,logscale=FALSE) {
    m <- polyroot(c(s2-s1,-3*(s1+s2),3*(s2-s1),-(s1+s2)))
    alpha.mode <- Re(m[which.min(abs(m))])
    vinv <- 3*s2/(1+alpha.mode)^4 + 3*s1/(1-alpha.mode)^4 + (a-1)/(1+alpha.mode)^2 + (b-1)/(1-alpha.mode)^2
    sd.mode <- sqrt(1/vinv)
    dnorm(alpha,mean=alpha.mode,sd=sd.mode,log=logscale)
}

ralphaPropT <- function(s1,s2,n,a,b) {
    #Draw from truncated T approximation to full posterior of alpha (s1,s2 can be vectors)
    # - s1,s2: parameters from the likelihood term
    # - n: sample size
    # - a,b: prior parameters 0.5(alpha+1) ~ Beta(a,b)
    # Output: list with following elements
    # - alpha: proposed values for alpha ~ N(alpha.mode,sd.mode)
    # - alpha.mode, sd.mode: vector mean and standard deviation of approximating Normal
    # - logpdfProp: log proposal density evaluated at the returned alpha
    alpha.mode <- double(length(s1))
    for (i in 1:length(s1)) {
        m <- polyroot(c(s2[i]-s1[i],-3*(s1[i]+s2[i]),3*(s2[i]-s1[i]),-(s1[i]+s2[i])))
        alpha.mode[i] <- Re(m[which.min(abs(m))])
    }
    vinv <- 3*s2/(1+alpha.mode)^4 + 3*s1/(1-alpha.mode)^4 + (a-1)/(1+alpha.mode)^2 + (b-1)/(1-alpha.mode)^2
    sd.mode <- sqrt(1/vinv)
    u <- runif(length(s1),pt((-1-alpha.mode)/sd.mode,df=n),pt((1-alpha.mode)/sd.mode,df=n))
    alpha <- qt(u,df=n)*sd.mode + alpha.mode
    logpdfProp <- dt((alpha-alpha.mode)/sd.mode,df=n,log=TRUE)
    #plot(dt((aseq-alpha.mode[1])/sd.mode[1],df=n,log=TRUE), (-(n+1)/2) * log(1+(1/n)*((aseq-alpha.mode[1])/sd.mode[1])^2))
    return(list(alpha=alpha,alpha.mode=alpha.mode,sd.mode=sd.mode,logpdfProp=logpdfProp))
}

ralphaPropNorm <- function(s1,s2,a,b) {
    #Draw from truncated Normal approximation to full posterior of alpha (s1,s2 can be vectors)
    # - s1,s2: parameters from the likelihood term
    # - a,b: prior parameters 0.5(alpha+1) ~ Beta(a,b)
    # Output: list with following elements
    # - alpha: proposed values for alpha ~ N(alpha.mode,sd.mode)
    # - alpha.mode, sd.mode: vector mean and standard deviation of approximating Normal
    # - logpdfProp: log proposal density evaluated at the returned alpha
    alpha.mode <- double(length(s1))
    for (i in 1:length(s1)) {
        m <- polyroot(c(s2[i]-s1[i],-3*(s1[i]+s2[i]),3*(s2[i]-s1[i]),-(s1[i]+s2[i])))
        alpha.mode[i] <- Re(m[which.min(abs(m))])
    }
    vinv <- 3*s2/(1+alpha.mode)^4 + 3*s1/(1-alpha.mode)^4 + (a-1)/(1+alpha.mode)^2 + (b-1)/(1-alpha.mode)^2
    sd.mode <- sqrt(1/vinv)
    u <- runif(length(s1),pnorm(-1,alpha.mode,sd=sd.mode),pnorm(1,alpha.mode,sd.mode))
    alpha <- qnorm(u,alpha.mode,sd=sd.mode)
    logpdfProp <- dnorm(alpha,mean=alpha.mode,sd=sd.mode,log=TRUE)
    return(list(alpha=alpha,alpha.mode=alpha.mode,sd.mode=sd.mode,logpdfProp=logpdfProp))
}



rnuSkewtGibbs <- function(x, mu, Sigma, A, D, alpha, nuprobs, ttype) {
#Sample from the conditional posterior of the skew-t degrees of freedom parameter nu (assumed to take values in names(nuprobs)) given mu, Sigma, alpha
# Input
# - x: matrix with observed data. Observations in rows, variables in columns
# - mu: location parameter
# - Sigma: scale matrix parameter
# - A: optional parameter containing pre-computed A=t(eigen(Sigma)). For identifiability all elements in A[,1] must be positive, if any entry is negative then the corresponding row is changed sign
# - D: optional parameter containing pre-computed A=diag(eigen(Sigma)$values)
# - alpha: vector with asymmetry coefficients
# Return: value of nu
    nuvals <- as.numeric(names(nuprobs))
    pp <- sapply(nuvals, function(v) sum(dskewt(x=x,mu=mu,A=A,D=D,alpha=alpha,nu=v,param='eps',logscale=TRUE,ttype=ttype))) + nuprobs
    pp <- exp(pp - max(pp))
    pp <- pp/sum(pp)
    nuvals[match(TRUE,runif(1)<cumsum(pp))]
}
