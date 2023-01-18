setwd()
library(coda)
library(graphicalVAR)
library(nimble)
library(tidyverse)
source('fit_mlMSVAR_ranRes_PARALLEL.R')
source('simulate_mlMSVAR.R')
source('makeMat.R')
source('makeSingleMat.R')
source('lse.R')
source('init_mlMSVAR.R')
source('lse.R')
source('dForProbsMSVARnew.R')


##### simulate data ######
#model with 2 latent states and 3 variables

M <- 2 # number of latent state
probs1 <- c(.7,.3) #initial state probailities at time 1
n <- 10 # number of individuals
nTime <- 50 # number of time points
nVars <- 3 #number of variables

#generate average intercepts
b1.means <- array(c(rnorm(nVars,-3,1),rnorm(nVars, 3,1)), c(nVars,M))
#variances of average intercepts
b1.sd1 <-diag(3)
b1.sd2 <-diag(3)
#store in an array
b1.sd <- array(c(as.vector(b1.sd1), as.vector(b1.sd2)), c(nVars, nVars, M))

#generate initial average VAR matrices
var1 <- matrix(runif(nVars*nVars, -1,1),3,3)
var2 <- matrix(runif(nVars*nVars, -1,1),3,3)
b.means <- array(c(as.vector(var1), as.vector(var2)), c(nVars*nVars,M))

#generate initial covariance matrices of VAR effects
b.sd1 <- LaplacesDemon::rwishart(9,diag(9)/80)
b.sd2 <-  LaplacesDemon::rwishart(9,diag(9)/80)
b.sd <- array(c(as.vector(b.sd1),as.vector(b.sd2)), c(nVars*nVars, nVars*nVars, M))

#generate average transition logits and their covariance matrix
v.means <- c(-1,-.5)
v.cov <- diag(2)/3

#generate scale matrix to simulate random residual matrices
scaleResCov <- array(c(rwish_chol(1, diag(nVars)*5, nVars),
                      rwish_chol(1, diag(nVars)*5, nVars)),
                c(nVars,nVars,M))

#simulate data with above parameters
data <- simulate_mlMSVAR(probs1 = probs1,
                         n = n,
                         nTime = nTime,
                         nItems = nVars,
                         b1.means = b1.means,
                         b.means = b.means,
                         b1.sd = b1.sd,
                         b.sd = b.sd,
                         v.means = v.means,
                         v.cov = v.cov,
                         dfResCov = c(n/2,n/2),
                         scaleResCov=scaleResCov)
# note that, most likely, parameters needed to be shrinked to ensure stationarity
# in this case, the model is not generated with parameters specified above
# but with the shrinked parameters stored in the data object (e.g. data$b.means)

# with this small data set a relatively loose convergence criteria (neff=10 and rhat=1.1), 
#estimation will take maybe 5-30 minutes 
#(I guess it mainly depends on how good the initialization was)

# with larger data sets and stricter convergence criteria, 
#it can take up to a few days

res <- fit_mlMSVAR(data$y, M = 2, niter = 10000, nVars=nVars, 
                   ncores = 2, neff=10, rhat=1.1)

#view summary of results
res$summary

