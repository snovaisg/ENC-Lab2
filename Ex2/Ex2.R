## Exercise 2
rm(list=ls())

library(rstudioapi) # to automatically set the working directory to this file's path.

#set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

x = c(0, 2, 3, 0, 0, 0, 0, 2, 3, 1, 0, 1, 0, 1, 0, 2, 5, 2, 5, 1, 2, 2, 1, 0, 0)
n = length(x)

# Maximum Likelihood Estimator of p
p.mle = 1/(1+mean(x))
p.mle
# 0.4310345

# likelihood function
lik <- function(p){
  (p^n)*((1-p)^sum(x))
}

#log-likelihood function
loglik <- function(p){
  n*log(p) + log(1-p)*sum(x)
} 

# score function
s <- function(p){
  (n/p)-(1/(1-p))*sum(x)
}


# Display the likelihood, log-likelihood and score functions graphically
# in order to locate the ML estimate of p
library(Cairo)
CairoPDF("lik_loglik_score_Geom.pdf",width=9,height=5)
par(mfrow=c(1,3), oma = c(2, 5, 2, 0), mar = c(3.1, 1.5, 2.1, 2.1))
p <- seq(0,1,0.01)
plot(p,lik(p),ylab="likelihood",xlab=expression(p),lwd=2,type="l", 
     main = "Likelihood")
box(lwd=2)
plot(p,loglik(p),ylab="log-likelihood",xlab=expression(p),lwd=2,type="l", 
     main = "Loglikelihood", ylim = c(-130,-40))
box(lwd=2)
plot(p,s(p),ylab="score",xlab=expression(p),lwd=2,type="l", 
     main = "Score", ylim = c(-2000, 2000))
abline(h=0,lty=3)
box(lwd=2)
dev.off()


# Jackknife
p.jack=numeric(n)
for(i in 1:n){
  p.jack[i] = 1/(1+mean(x[-i]))
}

# Estimative p 
Jack.p = mean(p.jack); Jack.p
#  0.4313412

# Standard error of the MLE
se.jack = sqrt((n-1)*mean((p.jack-Jack.p)^2)); se.jack
# 0.05726061

# Bias of the MLE
bias.jack = (n-1)*(Jack.p-p.mle); bias.jack
# 0.007362259



# Bootstrap
B=5000
p.boot=numeric(B)
set.seed(123)

for(i in 1:B){
  x.boot=sample(x,n,replace=T)
  p.boot[i]=(1/(1+mean(x.boot)))
}

# Estimative p 
boot.p = mean(p.boot); boot.p
# 0.4384257

# Standard error of the MLE
se.boot=sd(p.boot); se.boot
# 0.05590814

# Bias of the MLE
bias.boot=mean(p.boot)-p.mle; bias.boot
# 0.007391179