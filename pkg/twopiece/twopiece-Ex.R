pkgname <- "twopiece"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('twopiece')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dtp3")
### * dtp3

flush(stderr()); flush(stdout())

### Name: dtp3
### Title: The 3-Parameter Two Piece Distribution
### Aliases: dtp3 ptp3 qtp3 rtp3

### ** Examples


## 3-parameter two piece normal density with parameterization "tp"
tempf = function(x) dtp3(x,0,3,1,dnorm,param="tp")
curve(tempf,-10,5)

## 3-parameter two piece normal distribution with parameterization "tp"
tempf = function(x) ptp3(x,0,1,3,pnorm,param="tp")
curve(tempf,-10,10)

## random number generation for 3-parameter two piece normal distribution
## with parameterization "tp"
sim <- rtp3(1000,0,1,1,rnorm)
hist(sim,probability=TRUE)

## quantile function for the 3-parameter two piece normal distribution
## with parameterization "tp"
qtp3(0.5, 0, 1, 1, qnorm ,param = "tp")



cleanEx()
nameEx("dtp4")
### * dtp4

flush(stderr()); flush(stdout())

### Name: dtp4
### Title: The 4-Parameter Two Piece Distribution
### Aliases: dtp4 ptp4 qtp4 rtp4

### ** Examples


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




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
