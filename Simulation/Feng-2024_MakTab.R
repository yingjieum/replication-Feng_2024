#########################################################
##### This file produces Table 3 in the main paper#######
########################################################
library(Hmisc)
library(Rfast)
source("Feng-2024_SimFuns.R")

############################################
# Table used in the main paper
par <- read.csv("Feng-2024_SimModels.csv", header = T, colClasses=c(rep("numeric", 6), "logical"))

mat <- matrix(NA, 15, 8)
sum.K <- NULL
for (j in 1:6) {
   
n   <- par$n[j]
p   <- par$p[j]
model <- funlist[[par$hdmodel[j]]]
nlam <- par$nlam[j]
K    <- n^(4/5)
Kseq <- floor(K * c(0.5, 1, 1.5))
###############################################
# calculate true value
p1 <- integrate(px, 0, 1)$value
integrand <- function(alpha) {
  return(mu0(alpha)*px(alpha))
}
theta0 <- (integrate(integrand, 0, 1)$value)/p1
#####################################################

output <- as.matrix(read.table(paste("output/rawoutput_par", j ,"txt", sep="."), sep = ","))

rep <- 2000
rownum <- rep(1:5, rep)
est   <- output[rownum==2,] 
se    <- output[rownum==3,]
rej   <- output[rownum==4,]
ci    <- output[rownum==5,]
k.cv  <- output[rownum==1,4]


bias  <- colMeans(est) - theta0
sd    <- sqrt(colVars(est))
rmse  <- sqrt(colMeans((est - theta0)^2))
CR    <- 1-colMeans(rej)
AL    <- colMeans(ci)
table <- round(rbind(bias, sd, rmse, CR, AL), 3)

sum.K <- rbind(sum.K, c(summary(k.cv), sd(k.cv)))

if (j==1) {
   mat[6:10, 1:4] <- table  
} else if (j==2) {
   mat[6:10, 5:8] <- table
} else if (j==3) {
   mat[11:15, 1:4] <- table
} else if (j==4) {
   mat[11:15, 5:8] <- table
} else if (j==5) {
   mat[1:5, 1:4] <- table
} else if (j==6) {
   mat[1:5, 5:8] <- table
}

}

# Report main simulation results
n.rgroup <- c(5, 5, 5)
rowname <- rep(c("BIAS", "SD", "RMSE", "CR", "AL"),3)
rgroup   <- c("Model 1", "Model 2", "Model 3")

n.cgroup <- rep(4, 2)
cgroup   <- c("Local PCA, $K=$", "Local average, $K=$")
colheads <- rep(c(paste(Kseq), "$\\widehat{K}_{\\mathtt{CV}}$"), 2)
latex(mat, file=paste("Table_Pointwise_MainPaper", ".txt", sep = ""), 
      append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
      n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
      )


