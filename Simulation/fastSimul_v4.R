###################################################
#### This file replicates the simulation results ##
############## in Feng (2020)  ####################
########### Author: Yingjie Feng ##################
########## Last updated: 5/16/2020 ################

rm(list = ls())
setwd("E:/Dropbox/000000_JMP/Replication/Simulation/")
library(Rfast)
library(foreach)
library(doParallel)
library(doRNG)
#library(gpuR)
library(irlba)
library(RcppNumerical)
source("SuppFuns.R")

# param
rep <- 5000
par <- read.csv("SimulModel.csv", header = T, colClasses=c(rep("numeric", 6), "logical"))
j   <- 2
n   <- par$n[j]
p   <- par$p[j]
model <- funlist[[par$hdmodel[j]]]
#family.y <- gaussian(link="identity")
#family.d <- binomial(link="logit")    
nlam <- par$nlam[j]
K    <- n^(2*nlam/(2*nlam+1))
Kseq <- floor(K * c(0.5, 1, 1.5))
const <- par$const[j]

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
                     output <- sim(i, n, p, model, Kseq, nlam, theta0, const)
                     output   # (4*rep) by length(Kseq) matrix
                   }

stopCluster(cl)
###################
write.table(output, paste("rawoutput_par", j, "txt", sep = "."), sep = ",", 
            row.names = F, col.names = F)

