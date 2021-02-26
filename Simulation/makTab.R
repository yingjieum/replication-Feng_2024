################################################
##### This make Table 3 in the main paper#######
################################################
library(Hmisc)
library(Rfast)
setwd("E:/Dropbox/000000_JMP/Replication/Simulation/")
source("SuppFuns.R")

############################################
# Table used in the main paper
par <- read.csv("SimulModel.csv", header = T, colClasses=c(rep("numeric", 6), "logical"))

mat <- matrix(NA, 6, 11)
for (j in 1:4) {
n   <- par$n[j]
p   <- par$p[j]
model <- funlist[[par$hdmodel[j]]]
nlam <- par$nlam[j]
K    <- n^(2*nlam/(2*nlam+1))
Kseq <- floor(K * c(0.5, 1, 1.5))
###############################################
# calculate true value
p1 <- integrate(px, 0, 1)$value
integrand <- function(alpha) {
  return(mu0(alpha)*px(alpha))
}
theta0 <- (integrate(integrand, 0, 1)$value)/p1
#####################################################

output <- as.matrix(read.table(paste("rawoutput_par", j ,"txt", sep="."), sep = ","))

rep <- 5000
rowname <- rep(1:4, rep)
est   <- output[rowname==1,] 
se    <- output[rowname==2,]
rej   <- output[rowname==3,]
ci    <- output[rowname==4,]

bias  <- colMeans(est) - theta0
sd    <- sqrt(colVars(est))
rmse  <- sqrt(colMeans((est - theta0)^2))
CR    <- 1-colMeans(rej)
AL    <- colMeans(ci)
table <- round(cbind(bias, sd, rmse, CR, AL),3)

if (j==1) {
   table <- cbind(Kseq, table)
   mat[1:3, 1:6] <- table  
} else if (j==2) {
   table <- cbind(Kseq, table)
   mat[4:6, 1:6] <- table
} else if (j==3) {
   mat[1:3, 7:11] <- table
} else {
   mat[4:6, 7:11] <- table
}

}


n.cgroup <- c(1, 5, 5)
colheads <- c("$K$", rep(c("BIAS", "SD", "RMSE", "CR", "AL"),2))
cgroup   <- c("", "Model 1", "Model 2")

n.rgroup <- rep(3, 2)
rgroup   <- c("local linear", "local constant")
rowname  <- rep("", 6)
latex(mat, file=paste("Table_Pointwise_MainPaper", ".txt", sep = ""), 
      append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
      n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
      )