library(rstudioapi) # to automatically set the working directory to this file's path.

#set the working directory to this file's path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


A = as.vector(read.csv("data/data.csv",header = F)); A

