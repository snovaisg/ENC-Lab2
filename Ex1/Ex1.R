library(rstudioapi) # to automatically set the working directory to this file's path.
library(e1071)  # to compute skewness of a distribution

#set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Ok let's begin
A = as.vector(read.csv("data/data.csv",header = F)); A
n = length(A); n

# Ex1 a)
X = length(A[A > 1200]); X
P = X/n;P

# Ex1 b)

B = 10000;
boot_samples = matrix(nrow=B,ncol=20)
bootPvector = c()
sd_boot_sample = c()
set.seed(123)
# gen 10000 bootstrap samples
for (i in 1:B){
  # generated 1 bootstrap sample
  bi_sample = sample(as.numeric(A),n,replace = TRUE)
  # save the sample
  boot_samples[i,] = bi_sample
  # estimator P applied to each bootstrap sample
  pi = length(bi_sample[bi_sample > 1200]) / n
  # save the estimator of this sample
  bootPvector = c(bootPvector, pi)
}
# bootstrap estimate of P
bootP = mean(bootPvector); bootP
# bootstrap estimate of the bias of P
boot_bias_P = bootP - P; boot_bias_P
# bootstrap estimate of the variance of P - theta
boot_var_P = 1/(B-1)*sum((bootPvector - bootP)^2); boot_var_P

# Now let's correct the estimates with the bias
bias_corrected_P = P - boot_bias_P; bias_corrected_P

# optional - visualizing the distribution of bootP - P
hist(bootPvector - P,xlim=c(-0.5,0.5))
# Conclusions
# it looks normal, centered at zero. As seen above, 
# the expected value of this bias is boot_bias_P = 0.00023.
# also, the skew looks to be about zero with a slight inclination to the left (-> negative)
skew_boot_bias_P = skewness(bootPvector - P); skew_boot_bias_P
# -0.1725537

###################################################
############ Confidence intervals #################
###################################################

# A) Percentile confidence interval
ci_percentile = quantile(bootPvector, c(.025,.975)); ci_percentile

# B) Basic(Pivot) confidence interval
deltastar = bootPvector - P
d = quantile(deltastar, c(.025,.975));
ci_basic = P - c(d[2], d[1]); ci_basic
# Percentile was exactly the same as the basic!

# C) Studentized confidence interval
# Steps
# 1º - for each boot_sample already generated, generate M bootstrap samples, 
#      for each of them (the 1,...,M samples), compute the statistic P and finally get its standard deviation
#      which represents the sd of the bootstrap estimator of the initial sample. Repeat for the B samples
#      and save the results in a vector called sd_boot_P.
# 2º - compute the delta_star's 
# 3º - get the quantiles of the delta_star's.
# 4º- compute the confidence interval through the formula

# 1º
sd_boot_P = c()
set.seed(123)
M = 100
for (boot_sample_idx in 1:nrow(boot_samples)){
  boot_sample = boot_samples[boot_sample_idx,]
  
  # vector will contain the estimated P's
  boot_P = c()
  for (m in 1:M){
    boot_sample_ception = sample(as.numeric(boot_sample), length(boot_sample), replace = TRUE)
    pi = length(boot_sample_ception[boot_sample_ception > 1200]) / length(boot_sample)
    boot_P = c(boot_P, pi)
  }
  # compute sd
  sd_P = sd(boot_P)
  sd_boot_P = c(sd_boot_P, sd_P)
}
# 2º
deltastar = (bootPvector - P) / sd_boot_P; 
# 3º
d = quantile(deltastar, c(.025,.975));
# 4ª
ci = P - c(d[2],d[1])*sd(bootPvector); ci


####################################
############ Plots #################
####################################
png(filename="media/histogram_boot_ps.png")
hist(bootPvector,main='Histogram of boostrap estimated p\'s')
dev.off()

