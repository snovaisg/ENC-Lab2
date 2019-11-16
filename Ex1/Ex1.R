library(rstudioapi) # to automatically set the working directory to this file's path.

#set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

A = as.vector(read.csv("data/data.csv",header = F)); A
n = length(A); n

# Ex1 a)
X = length(A[A > 1200]); X
P = X/n;P

# Ex1 b)

B = 10000;
boot_samples = matrix(nrow=B,ncol=20)
bootPvector = c()
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
# bootstrap estimate of the variance of P to the actual parameter p
boot_var_P = 1/(B-1)*sum((bootPvector - bootP)^2); boot_var_P

# Now let's correct the estimates with the bias
bias_corrected_P = P - boot_bias_P; bias_corrected_P

####################################
############ Plots #################
####################################
png(filename="media/histogram_boot_ps.png")
hist(bootPvector,main='Histogram of boostrap estimated p\'s')
dev.off()

