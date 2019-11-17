rm(list=ls())

library(rstudioapi) # to automatically set the working directory to this file's path.

#set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in data
x <- c(-1.16225338, -0.01236762, 0.23144468, -2.08263805, 2.64870304, -0.52868938,
       -1.52280636, 0.03085357, -1.54249244, -0.23164183, 2.23750583, -0.64678326,
       -0.82817693, -1.50508448, 0.87125402, -1.67430203, 1.63702338, -0.85729792,
       -0.67855079, -1.36315013, -1.70412209, -0.05967576, -0.88579241, -0.33469221,
       -0.74940615)

n <- length(x)
epsil = 0.000001 # relative convergence criterion as a stop rule 

# Likelihood function
lik <- function(lambda){
  lambda^n * e^(-lambda*sum(x))
}

# loglikelihood function
loglik <- function(lambda){
  n*log(lambda)-lambda*sum(x)
}

# score function
s <- function(lambda){
  n/lambda + sum(x)
}

# score derivative
s.prime <- function(lambda){
  -n/lambda^2
}

# (a) Use the observed data and the method of Newton-Raphson to approximate the ML 
# estimate of ??. Justify your choice of the initial estimate and all the intermediate steps.

NR      <- function(lambda0,eps){
  # x      : observed sample
  # [lambda0: initial estimate of lambda
  # eps    : is the stopping rule
  
  # the NR method iterates until the stopping criterion is validated
  # we will use as stopping rule diff = |lambda_(t+1)-lambda_(t)| <= eps
  lambda.it    = vector()
  lambda.it[1] = alpha0
  k           = 1
  diff        = 1
  while(diff>eps){
    lambda.it[k+1] = lambda.it[k]-s(lambda.it[k])/s.prime(lambda.it[k])
    diff          = abs(alpha.it[k+1]-alpha.it[k])
    k             = k+1
  }
  result = as.matrix(lambda.it)
  colnames(result)<-"iterates"
  rownames(result)<-1:length(lambda.it)
  result
}

# (b) Consider the following reparametrization b = log(??) in the pdf above. Implement the 
# Newton-Rapshon algorithm to approximate the ML estimate of ?? using this reparametrization

# (c) Which approach, (a) or (b) is more sensitive to the initial values?