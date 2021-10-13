#########################################################
##### This file produces Table 3 in the main paper#######
########################################################
library(Hmisc)
library(Rfast)
source("Feng-2021_SimFuns.R")

############################################
# Table used in the main paper
par <- read.csv("Feng-2021_SimModels.csv", header = T, colClasses=c(rep("numeric", 6), "logical"))

mat <- matrix(NA, 8, 10)
sum.K <- NULL
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
rownum <- rep(1:5, rep)
est   <- output[rownum==2,] 
se    <- output[rownum==3,]
rej   <- output[rownum==4,]
ci    <- output[rownum==5,]
k.dpi <- output[rownum==1,4]


bias  <- colMeans(est) - theta0
sd    <- sqrt(colVars(est))
rmse  <- sqrt(colMeans((est - theta0)^2))
CR    <- 1-colMeans(rej)
AL    <- colMeans(ci)
table <- round(cbind(bias, sd, rmse, CR, AL),3)

sum.K <- rbind(sum.K, c(summary(k.dpi), sd(k.dpi)))

if (j==1) {
   rowname <- c(paste(Kseq), "$\\widehat{K}_{\\mathtt{DPI}}$")
   mat[1:4, 1:5] <- table  
} else if (j==2) {
   rowname <- c(rowname, paste(Kseq), "$\\widehat{K}_{\\mathtt{DPI}}$")
   mat[5:8, 1:5] <- table
} else if (j==3) {
   mat[1:4, 6:10] <- table
} else {
   mat[5:8, 6:10] <- table
}

}

# Report main simulation results
n.cgroup <- c(5, 5)
colheads <- rep(c("BIAS", "SD", "RMSE", "CR", "AL"),2)
cgroup   <- c("Model 1", "Model 2")

n.rgroup <- rep(4, 2)
rgroup   <- c("Local linear, $K=$", "Local constant, $K=$")
#rowname  <- rep("", 6)
latex(mat, file=paste("Table_Pointwise_MainPaper", ".txt", sep = ""), 
      append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
      n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
      )


# Report the summary of DPI selection (NOT reported in the paper)
sum.K[,c(4,7)] <- round(sum.K[,c(4,7)],2)
n.cgroup <- 7
colheads <- rep(c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Sd."), 1)
cgroup   <- NULL

n.rgroup <- c(2,2)
rgroup   <- c("Model 1", "Model 2")
rowname  <- rep(c("Local linear", "Local constant"), 2)
latex(sum.K, file=paste("Table_Pointwise_MainPaper", ".txt", sep = ""), 
      append=FALSE, table.env=FALSE, center="none", title="",
      n.cgroup=n.cgroup, cgroup=cgroup, colheads=colheads,
      n.rgroup=n.rgroup, rgroup=rgroup, rowname=rowname
)


