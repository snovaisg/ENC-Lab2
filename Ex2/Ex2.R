## Exercise 2
rm(list=ls())

x = c(0, 2, 3, 0, 0, 0, 0, 2, 3, 1, 0, 1, 0, 1, 0, 2, 5, 2, 5, 1, 2, 2, 1, 0, 0)
n = length(x)

p.mle = 1/(1+mean(x))

# Jackknife
p.jack=numeric(n)
for(i in 1:n){
  p.jack[i] = 1/(1+mean(x[-i]))
}

Jack.p = mean(p.jack); Jack.p 
se.jack = sqrt((n-1)*mean((p.jack-Jack.p)^2)); se.jack
bias.jack = (n-1)*(Jack.p-p.mle); bias.jack



# Bootstrap
B=5000
p.boot=numeric(B)
set.seed(123)

for(i in 1:B){
  x.boot=sample(x,n,replace=T)
  p.boot[i]=(1/(1+mean(x.boot)))
}

boot.p = mean(p.boot); boot.p
se.boot=sd(p.boot); se.boot
bias.boot=mean(p.boot)-p.mle; bias.boot




