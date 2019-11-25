rm(list=ls())

library(rstudioapi) # to automatically set the working directory to this file's path.

#set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in data
x <- c(0.58169466, 0.39766226, 0.91658956, 0.02177749, 0.03876619, 0.21827670, 0.21670848,
       0.10018400, 1.88016308, 0.02010583, 0.69298625, 0.33118257, 0.19380250, 0.26008126,
       0.12985106, 0.58605940, 1.07807141, 0.33017960, 0.40754127, 2.78690463, 0.58148257,
       0.66611808, 1.02432813, 0.92968585, 0.80588206)

n <- length(x)

# compute alpha analytically for max lik
analytically <- n/sum(x); analytically # 1.645161

epsil = 0.000001 # relative convergence criterion as a stop rule 


# Likelihood function
lik <- function(lambda){
  lambda^n * exp(-lambda * sum(x))
}

# loglikelihood function
loglik <- function(lambda){
  n*log(lambda)-lambda*sum(x)
}

# score function
s <- function(lambda){
  n/lambda - sum(x)
}

# score derivative
s.prime <- function(lambda){
  -n/lambda^2
}

# max value at:
maxval <- lik(analytically); maxval #3.530575e-06


# Plot of the Likelihood, loglikelihood and score function to look for an estimator
library(Cairo)
CairoPDF("lik_loglik_score_Exp.pdf",width=9,height=5)
par(mfrow=c(1,3), oma = c(2, 5, 2, 0), mar = c(3.1, 1.5, 2.1, 2.1))
lambda <- seq(0,3.5,0.05)
plot(lambda,lik(lambda),ylab="likelihood",xlab=expression(lambda),lwd=2,type="l", 
     main = "Likelihood")
box(lwd=2)
plot(lambda,loglik(lambda),ylab="log-likelihood",xlab=expression(lambda),lwd=2,type="l", 
     main = "Loglikelihood")
box(lwd=2)
plot(lambda,s(lambda),ylab="score",xlab=expression(lambda),lwd=2,type="l", 
     main = "Score")
abline(h=0,lty=3)
box(lwd=2)
# mtext(expression(paste("f(x;", lambda , ") =", lambda, "e^(",lambda,"x)")), outer = TRUE, cex = 1.5)
dev.off()


# iteratively searching for an approximation using the estimation of moments
mean(x)
# 0.6078434

lambda <- seq(1,2.4,0.1)
k = 1
for(i in lambda){
  print(c(k,lambda[k],
          (1/(lambda[k]))))
  k= k+1
}

# 1           1          1
# 2.0000000   1.1000000  0.9090909
# 3.0000000   1.2000000  0.8333333
# 4.0000000   1.3000000  0.7692308
# 5.0000000   1.4000000  0.7142857
# 6.0000000   1.5000000  0.6666667
# 7.000       1.600      0.625
# 8.0000000   1.7000000  0.5882353  
# 9.0000000   1.8000000  0.5555556
# 10.0000000  1.9000000  0.5263158
# 11.0        2.0        0.5
# 12.0000000  2.1000000  0.4761905
# 13.0000000  2.2000000  0.4545455
# 14.0000000  2.3000000  0.4347826
# 15.0000000  2.4000000  0.4166667

# we see that the value closest to the mean 0.676118 is achieved by lambda = 1.7
mme.mean=1.7

# using the maxlik() function to determine a ML estimator
library(maxLik)
maxLik(logLik=loglik,start=mme.mean)
# Maximum Likelihood estimation
# Newton-Raphson maximisation, 3 iterations
# Return code 1: gradient close to zero
# Log-Likelihood: -12.55405 (1 free parameter(s))
# Estimate(s): 1.645161



# (a) Use the observed data and the method of Newton-Raphson to approximate the ML 
# estimate of ??. Justify your choice of the initial estimate and all the intermediate steps.

NR      <- function(lambda0,eps){
  # x      : observed sample
  # lambda0: initial estimate of lambda
  # eps    : is the stopping rule
  
  # the NR method iterates until the stopping criterion is validated
  # we will use as relative convergence rule 
  # diff = |lambda_(t+1)-lambda_(t)|/|lambda_(t)| <= eps
  lambda.it    = vector()
  lambda.it[1] = lambda0
  k           = 1
  diff        = 1
  broke = FALSE
  while(!broke && diff>eps){ # to see whether method deverges
    lambda.it[k+1] = lambda.it[k]-s(lambda.it[k])/s.prime(lambda.it[k])
    if (lambda.it[k+1] > 0){
      diff          = abs(lambda.it[k+1]-lambda.it[k])/abs(lambda.it[k])
    }else{
      broke = TRUE
    }
    k             = k+1
  }
  result = as.matrix(lambda.it)
  colnames(result)<-"iterates"
  rownames(result)<-1:length(lambda.it)
  result
}

# Newton Raphson with graphical estimation
mme.graph = 2
t(NR(mme.graph,epsil))
#          1        2      3        4        5        6
# iterates 2 1.568626 1.6416 1.645153 1.645161 1.645161

# # Newton Raphson with approximation of the mean
t(NR(mme.mean,epsil))
#            1        2        3        4        5
# iterates 1.7 1.643333 1.645159 1.645161 1.645161


mme.E = 1/mean(x)
t(NR(mme.E, epsil))
#                 1        2
# iterates 1.645161 1.645161

# # Newton Raphson starting at almost 0 (epsilon = 0.000001)
t(NR(epsil,epsil)) 
#              1            2            3            4            5           6    
# iterates 1e-06 1.999999e-06 3.999996e-06 7.999983e-06 1.599993e-05 3.19997e-05 
#                    7        8            9           10          11          12       
# iterates 6.399877e-05 0.0001279951 0.0002559802 0.0005119205 0.001023682 0.002046726 
#                13        14         15         16         17        18        19       
# iterates 0.004090907  0.00817164 0.01630269 0.03244383 0.06424785 0.1259866 0.2423252 
#               20        21       22       23      24       25       26       27
# iterates 0.448957 0.7753956  1.185332 1.516637 1.63512 1.645099 1.645161 1.645161


# Newton Raphson starting at almost 0 (epsilon = 0.000001)
t(NR(5,epsil))
#          1         2
# iterates 5 -5.196085


CairoPDF("score_exp.pdf",width=9,height=5)
lambda <- seq(1.5,10,0.1)
plot(lambda,s(lambda),ylab="score",xlab=expression(lambda),lwd=2,type="l", 
     main = "Score"); abline(h=0,lty=3); box(lwd=2)
dev.off()


# (b) Consider the following reparametrization b = log(??) in the pdf above. Implement the 
# Newton-Rapshon algorithm to approximate the ML estimate of ?? using this reparametrization

# loglikelihood function
loglik.b <- function(b){
  n*b-exp(b)*sum(x)
}

# score function
s.b <- function(b){
  n-exp(b)*sum(x)
  }

# derivative of the (pro???le) score function 
s.b.prime <- function(b){-exp(b)*sum(x)}


# Compute back
lambda.b.hat <- exp(0.497838); lambda.b.hat # 1.645161

library(Cairo)
CairoPDF("lik_loglik_score_Exp_reparam.pdf",width=9,height=5)
# par(mfrow = c(1,2))
par(mfrow=c(1,2), oma = c(2, 5, 2, 0), mar = c(3.1, 1.5, 2.1, 2.1))
b <- seq(-2,4,0.05)
plot(b,loglik.b(b),ylab="loglikelihood",xlab="b",lwd=2,type="l", 
     main = "Loglikelihood");box(lwd=2)

plot(b,s.b(b),ylab="Score",xlab="b",lwd=2,type="l", 
     main = "Score"); abline(h=0,lty=3); box(lwd=2)
dev.off()


mme.graph.b = 1
maxLik(logLik=loglik.b,start=mme.graph.b)
# Maximum Likelihood estimation
# Newton-Raphson maximisation, 4 iterations
# Return code 1: gradient close to zero
# Log-Likelihood: -12.55405 (1 free parameter(s))
# Estimate(s): 0.497838 

library(Cairo)
install.packages("autoimage")
library(autoimage)
reset.par() # reset parameters to default

# example using alpha = 5
slope <- s.b.prime(log(epsil)) # slope of tangent in x0 = 5 of score function
p.epsil <- s.b(log(epsil)) # point of score function
inter <- p.epsil - slope*log(epsil) # intercept of tangent
y0 <- (-inter/slope) # cut of y-axis

CairoPDF("score_exp_reparam_tangent.pdf",width=9,height=5)
b <- seq(-15,2,0.1)
plot(b,s.b(b),ylab="score",xlab=expression(b),lwd=2,type="l", 
     main = "Score tangent at lambda = epsilon"); abline(h=0,lty=3); box(lwd=2)
abline(h=0, v = 0, a = inter, b = slope,lty=c(1,3, 3),
       col=c("red","black", "black"))
box(lwd=2)
points(c(log(epsil), y0), c(p.epsil, 0), pch = 19, col = "red")
dev.off()

NR <- function(lambda0,eps){ 
  # x : observed sample 
  # lambda0 : initial estimate of lambda
  # eps : is the stopping rule (using relative convergence rule)
  b0 = log(lambda0)
  b.it = vector(); b.it[1] = b0
  lambda.it = vector(); lambda.it[1] = lambda0 
  k = 1; diff = 1 
  broke = FALSE
  while(!broke && diff>eps){ # to see whether method deverges
    b.it[k+1] = b.it[k]-s(b.it[k])/s.prime(b.it[k]) 
    lambda.it[k+1] = exp(b.it[k+1]) 
    if (lambda.it[k+1] > 0){
      diff = abs(b.it[k+1]-b.it[k])/abs(b.it[k])
    }else{
      broke = TRUE
    }
    k = k+1 
    } 
  result = as.matrix(cbind(b.it,lambda.it)) 
  rownames(result)<-1:length(lambda.it);
  result
}

# using graphical estimation
t(NR(mme.graph,epsil))
#                   1        2        3        4        5        6        7
# b.it      0.6931472 1.094254 1.460681 1.624474 1.644900 1.645161 1.645161
# lambda.it 2.0000000 2.986954 4.308895 5.075749 5.180494 5.181842 5.181842

# using approximation of the mean
t(NR(mme.mean,epsil))
#                   1         2        3        4        5        6        7        8
# b.it      0.5306283 0.8901083 1.298627 1.572167 1.641922 1.645154 1.645161 1.645161
# lambda.it 1.7000000 2.4353933 3.664261 4.817077 5.165087 5.181809 5.181842 5.181842

# using method of the moments
t(NR(mme.E,epsil))
#                  1         2        3        4        5        6        7        8
# b.it      0.497838 0.8450265 1.256010 1.553110 1.640010 1.645144 1.645161 1.645161
# lambda.it 1.645161 2.3280394 3.511384 4.726146 5.155222 5.181758 5.181842 5.181842


# using epsil
t(NR(epsil,epsil))
#                   1             2         3
# b.it      -13.815511 -1.436491e+02 -12830.18
# lambda.it   0.000001  4.111487e-63      0.00

# using 20
t(NR(20,epsil))
#                   1         2         3        4        5        6        7        8        9
# b.it       2.995732 0.5364276 0.8979454 1.305783 1.575151 1.642181 1.645155 1.645161 1.645161
# lambda.it 20.000000 1.7098875 2.4545549 3.690578 4.831471 5.166427 5.181814 5.181842 5.181842

# using 20
t(NR(5,epsil))
#                  1        2       3        4
# b.it      1.609438 1.644385 1.64516 1.645161
# lambda.it 5.000000 5.177824 5.18184 5.181842


## CrasherPlots
library(Cairo)
install.packages("autoimage")
library(autoimage)
reset.par() # reset parameters to default



# example using lambda = 5
slope <- s.prime(5) # slope of tangent in x0 = 5 of score function
p5 <- s(5) # point of score function
inter <- p5 - slope*5 # intercept of tangent
y0 <- (-inter/slope) # cut of y-axis

CairoPDF("score_exp_tangent.pdf",width=9,height=5)
lambda <- seq(-8,10,0.1)
plot(lambda,s(lambda),ylab="score",xlab=expression(lambda),lwd=2,type="l", 
     main = "Score tangent at lambda = 5", ylim = c(-40,40)); abline(h=0,lty=3); box(lwd=2)
abline(h=0, v = 0, a = inter, b = slope,lty=c(1,3, 3),
       col=c("red","black", "black"))
box(lwd=2)
points(c(5, y0), c(p5, 0), pch = 19, col = "red")
dev.off()


# Reparamtrization: example using epsilon
slope <- s.b.prime(log(epsil)) # slope of tangent in x0 = 5 of score function
p.epsil <- s.b(log(epsil)) # point of score function
inter <- p.epsil - slope*log(epsil) # intercept of tangent
y0 <- (-inter/slope) # cut of y-axis

CairoPDF("score_exp_reparam_tangent.pdf",width=9,height=5)
b <- seq(-15,2,0.1)
plot(b,s.b(b),ylab="score",xlab=expression(b),lwd=2,type="l", 
     main = "Score tangent at lambda = epsilon"); abline(h=0,lty=3); box(lwd=2)
abline(h=0, v = 0, a = inter, b = slope,lty=c(1,3, 3),
       col=c("red","black", "black"))
box(lwd=2)
points(c(log(epsil), y0), c(p.epsil, 0), pch = 19, col = "red")
dev.off()
