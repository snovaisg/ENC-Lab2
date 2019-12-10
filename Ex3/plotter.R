# Ex 3 generator
rm(list=ls())

library(rstudioapi) # to automatically set the working directory to this file's path.

#set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


x <- c(0.5409477, 0.8184872, 0.7848854, 0.9850439, 0.8963032, 0.6089008, 0.9549606,
       0.6795304, 0.8451902, 0.5613979, 0.4029634, 0.2741569, 0.3996693, 0.6371445, 
       0.7521881)
n <- length(x)

# (a) Derive the likelihood, log-likelihood and score functions.

# likelihood 
lik    <- function(alpha){
  alpha^n*prod(x)^(alpha-1)
}

# log-likelihood
loglik <- function(alpha){
  n*log(alpha) + (alpha-1)*sum(log(x))
} 

# score function
s <- function(alpha){
  n/alpha + sum(log(x))
}

# score derivative function
s.prime <- function(alpha){
  -n/alpha^2
}



# stopping criteria
epsil = 0.000001

# graphical choice
a.init = 2
b.init = 2.5




library(Cairo)
install.packages("autoimage")
library(autoimage)
reset.par() # reset parameters to default


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
  broke = FALSE
  while(!broke && diff>eps){
    alpha.it[k+1] = alpha.it[k]+s(alpha.it[k])/I(alpha.it[k])
    if (alpha.it[k+1] > 0){
      diff          = abs(alpha.it[k+1]-alpha.it[k])
    }else{
      broke = TRUE
    }
    k             = k+1
  }
  result = as.matrix(alpha.it)
  colnames(result)<-"iterates"
  rownames(result)<-1:length(alpha.it)
  result
}





### newton plotter
NR      <- function(alpha0,eps){
  # alpha0: initial estimate of alpha
  # eps    : is the stopping rule
  
  # the NR method iterates until the stopping criterion is validated
  # we will use as the absolute convergence rule diff = |alpha_(t+1)-alpha_(t)| <= eps
  alpha.it    = vector()
  alpha.it[1] = alpha0
  k           = 1
  diff        = 1
  broke = FALSE
  while(!broke && diff>eps){ # to see whether method diverges
    alpha.it[k+1] = alpha.it[k]-s(alpha.it[k])/s.prime(alpha.it[k])
    if (alpha.it[k+1] > 0){
      diff          = abs(alpha.it[k+1]-alpha.it[k])
    }else{
      broke = TRUE
    }
    k             = k+1
  }
  result = as.matrix(alpha.it)
  colnames(result)<-"iterates"
  rownames(result)<-1:length(alpha.it)
  result
}

length(NR(2.5,epsil))
#            1        2        3        4
# iterates 2.2 2.234532 2.235083 2.235083
library(plotrix)
epsil <- 0.000001
init <- 5 
cols = c("#FF0000", "#FF8000", "#00CC00", "#00FFFF", "#0080FF", "#0000FF", "#7F00FF",
         "#FF33FF", "#CC0066", "#660033", "#330019", "#003366", "#006666", "#000000")
newton.plotter <- function(init){
  ralph <- NR(init, epsil)
  alpha <- seq(-4,7,0.1)
  plot(alpha,s(alpha),ylab=expression("score"),
       xlab=expression(alpha),lwd=2,type="l", ylim = c(-40,40), 
       main = expression("Tangents depending on initial value"))
  abline(h=0, v = 0, lty=c(3, 3), col=c("black", "black"))
  box(lwd=2)
  
  # adding tangentes in for loop
  for (i in (1:length(ralph))){
    slope <- s.prime(ralph[i]) # slope of tangent in x0 = 5 of score function
    p.i <- s(ralph[i]) # point of score function
    inter <- p.i - slope*ralph[i] # intercept of tangent
    y0 <- (-inter/slope) # cut of y-axis
    abline(a = inter, b = slope,lty=1, col=cols[i])
    points(c(ralph[i], y0), c(p.i, 0), pch = 19, col = cols[i])
    print(y0)
  }
  
}

library(Cairo)
CairoPDF("NRepsil.pdf",width=9,height=5)
newton.plotter(epsil+0.01)
dev.off()

CairoPDF("NR2.pdf",width=9,height=5)
newton.plotter(2)
dev.off()

CairoPDF("NR05.pdf",width=9,height=5)
newton.plotter(0.5)
dev.off()

