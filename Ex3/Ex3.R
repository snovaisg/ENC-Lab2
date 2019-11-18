rm(list=ls())

library(rstudioapi) # to automatically set the working directory to this file's path.

#set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


x <- c(0.5409477, 0.8184872, 0.7848854, 0.9850439, 0.8963032, 0.6089008, 0.9549606,
       0.6795304, 0.8451902, 0.5613979, 0.4029634, 0.2741569, 0.3996693, 0.6371445, 0.7521881)
n <- length(x)


# likelihood and log-likelihood
lik    <- function(alpha){
  alpha^n*prod((x^(alpha-1)))
  #((gamma(alpha+1)^n)/(gamma(alpha)^n))*prod(x^(alpha-1))
}

loglik <- function(alpha){
  n*log(alpha) + (alpha-1)*sum(log(x))
} 

s <- function(alpha){
  n/alpha + sum(log(x))
}


#---- plot of Beta function
alpha <- seq(.1,3,0.05)
plot(alpha, beta(alpha, 1))
plot(alpha, lik(alpha))

#---------------


# (b) Display the likelihood, log-likelihood and score functions graphically
# in order to locate the ML estimate of ??.

library(Cairo)
CairoPDF("lik_loglik_score.pdf",width=9,height=3)
par(mfrow=c(1,3))
alpha <- seq(1.5,2.9,0.1)
plot(alpha,lik(alpha),ylab="likelihood",xlab=expression(alpha),lwd=2,type="l")
box(lwd=2)
plot(alpha,loglik(alpha),ylab="log-likelihood",xlab=expression(alpha),lwd=2,type="l")
box(lwd=2)
plot(alpha,s(alpha),ylab="score",xlab=expression(alpha),lwd=2,type="l")
abline(h=0,v=2.235083,lty=3)
box(lwd=2)
dev.off()

#(c) Use the R function maxLik() from library maxLik to approximate the ML estimate of ??.
library(maxLik)
maxLik(logLik=loglik, start=mme)
#Maximum Likelihood estimation
#Newton-Raphson maximisation, 3 iterations
#Return code 1: gradient close to zero
#Log-Likelihood: 3.775335 (1 free parameter(s))
#Estimate(s): 2.235083 


# (d) Derive the algorithms of bisection, Newton-Raphson and Fisher scoring that 
# enable the approximation of the ML estimate of ??. Implement those in R and use the
# sample above to estimate ??. Justify your choice of the initial estimates. 

# BISECTION
# programming the bisection method, which needs only the score function
bisection <- function(a,b,eps){
  # x    : observed sample
  # [a,b]: interval where s verifies Bolzano's theorem
  # eps  : is the stopping rule
  
  # the bisection method iterates until the stopping criterion is validated
  # we will use as stopping rule diff = |alpha_(t+1)-alpha_(t)| < eps
  # alpha is initialized as the mid point of [a,b]
  alpha.it    = vector()
  alpha.it[1] = (a+b)/2
  k           = 1
  diff        = 1
  while(diff>eps){
    if(s(alpha.it[k])*s(a)<0){
      b             = alpha.it[k]
      alpha.it[k+1] = (a+b)/2
    }
    else{if(s(alpha.it[k])*s(a)>0){
      a             = alpha.it[k]
      alpha.it[k+1] = (a+b)/2
    }else{alpha.it[k+1]=alpha.it[k]}
    }
    diff = abs(alpha.it[k+1]-alpha.it[k])
    k = k+1
  }
  result = as.matrix(alpha.it)
  colnames(result)<-"iterates"
  rownames(result)<-1:length(alpha.it)
  result
}

epsil = 0.000001
t(bisection(2,2.5,epsil))
#             1     2      3       4        5        6        7        8        9       10
# iterates 2.25 2.125 2.1875 2.21875 2.234375 2.242188 2.238281 2.236328 2.235352 2.234863
#               11       12       13       14       15       16       17       18       19
# iterates 2.235107 2.234985 2.235046 2.235077 2.235092 2.235085 2.235081 2.235083 2.235084


##
# NEWTON-RAPHSON
# programming the NR method, which needs both the score and the score derivative functions
s.prime <- function(alpha){
  -n/alpha^2
}


NR      <- function(alpha0,eps){
  # x      : observed sample
  # [alpha0: initial estimate of alpha
  # eps    : is the stopping rule
  
  # the NR method iterates until the stopping criterion is validated
  # we will use as stopping rule diff = |alpha_(t+1)-alpha_(t)| <= eps
  alpha.it    = vector()
  alpha.it[1] = alpha0
  k           = 1
  diff        = 1
  while(diff>eps){
    alpha.it[k+1] = alpha.it[k]-s(alpha.it[k])/s.prime(alpha.it[k])
    diff          = abs(alpha.it[k+1]-alpha.it[k])
    k             = k+1
  }
  result = as.matrix(alpha.it)
  colnames(result)<-"iterates"
  rownames(result)<-1:length(alpha.it)
  result
}


# pesquisa de grelha pela estimativa dos momentos
mean(x)
# 0.676118

alpha <- seq(1.5,2.9,0.1)
k = 1
for(i in alpha){
  print(c(k,alpha[k], gamma(alpha[k]+1)/(gamma(alpha[k])*(alpha[k]+1))))
  k= k+1
}
# observamos que o valor mais proximo de 0.676118 e dado por alpha=2.1
mme.mean=2.1

mme.exp = mean(x)/(1-mean(x))

t(NR(0.6,epsil))

t(NR(mme.mean,epsil))
#            1        2        3        4        5
# iterates 2.1 2.226919 2.235053 2.235083 2.235083

t(NR(mme.exp,epsil))
#                 1        2       3        4        5
# iterates 2.087544 2.225344 2.23504 2.235083 2.235083



# FISHER SCORING
# programming the NR method, which needs both the score and the fisher information functions    
I <- function(alpha){
  -s.prime(alpha)}

SF      <- function(alpha0,eps){
  # x      : observed sample
  # alpha0: initial estimate of alpha
  # eps    : is the stopping rule
  alpha.it    = vector()
  alpha.it[1] = alpha0
  k           = 1
  diff        = 1
  while(diff>eps){
    alpha.it[k+1] = alpha.it[k]+s(alpha.it[k])/I(alpha.it[k])
    diff          = abs(alpha.it[k+1]-alpha.it[k])
    k             = k+1
  }
  result = as.matrix(alpha.it)
  colnames(result)<-"iterates"
  rownames(result)<-1:length(alpha.it)
  result
}

t(SF(0.6,epsil))
t(SF(mme.exp,0.000001))
#            1        2        3        4        5
# iterates 2.087544 2.225344 2.23504 2.235083 2.235083
