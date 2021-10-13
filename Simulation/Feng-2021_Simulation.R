###################################################
#### This file replicates the simulation results ##
############## in Feng (2021)  ####################
########### Author: Yingjie Feng ##################
########## Last updated: 10/13/2021 ################

rm(list = ls())
library(Rfast)
library(foreach)
library(doParallel)
library(doRNG)
library(irlba)
library(RcppNumerical)
source("Feng-2021_SimFuns.R")

# param
rep <- 5000
par <- read.csv("Feng-2021_SimModels.csv", header = T, colClasses=c(rep("numeric", 6), "logical"))

j   <- 1
n   <- par$n[j]
p   <- par$p[j]
model <- funlist[[par$hdmodel[j]]]
nlam   <- par$nlam[j]
K      <- n^(2*nlam/(2*nlam+1))
Kseq   <- floor(K * c(0.5, 1, 1.5))
const  <- par$const[j]
if (nlam%%2 == 0) {
  odd.delta <- 1
} else {
  odd.delta <- 2
}

###############################################
# calculate true value
p1 <- integrate(px, 0, 1)$value
integrand <- function(alpha) {
  return(mu0(alpha)*px(alpha))
}
theta0 <- (integrate(integrand, 0, 1)$value)/p1
#####################################################


## simulation
cl <- makeCluster(32)
registerDoParallel(cl)

output <- foreach (i = 1:rep, .options.RNG=1234, .packages=c('Rfast','irlba','RcppNumerical'),
                   .combine=rbind) %dorng% {
                     output <- sim(i, n, p, model, Kseq, nlam, theta0, const, odd.delta)
                     output   # (5*rep) by (length(Kseq)+1) matrix
                   }

stopCluster(cl)
###################
write.table(output, paste("rawoutput_par", j, "txt", sep = "."), sep = ",", 
            row.names = F, col.names = F)


