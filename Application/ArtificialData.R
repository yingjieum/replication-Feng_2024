#############################################################
##### This file replicates Figure 1 and 2 in Feng (2020) ####
############## Last updated: 08/30/2020 #####################
#############################################################

rm(list=ls())
set.seed(1234)
library(ggplot2)
library(plotly)
library(processx)
library(R.matlab)
library(Rfast)
library(irlba)
library(RcppNumerical)
library(reshape2)
library(Hmisc)
pdf.width <- 6
pdf.height <- 4.5
setwd("E:/Dropbox/000000_JMP/replication/Application")
source("E:/Dropbox/000000_JMP/replication/Application/suppfn.R")

# Generate illustrative pictures
#################################################
n=200
Lambda1 <- runif(n); Lambda2 <- runif(n)
L1 <- Lambda1; L2 <- Lambda2; L3 <- 5*(Lambda1-.5)^2+5*Lambda2^3
L  <- cbind(L1, L2, L3)
Y  <- L + matrix(rnorm(n*3,0,0.1), n, 3)
data <- data.frame(x1=L[,1], x2=L[,2], x3=L[,3], y1=Y[,1], y2=Y[,2], y3=Y[,3])
grid <- seq(-0.1,1.1,0.1); surf <- expand.grid(grid, grid)
z <- matrix(5*(surf[,2]-.5)^2+5*surf[,1]^3, 13, 13)

## tangent space
F1 <- c(1,0,-2); F2 <- c(0,1,7.35)
grid.tan <- seq(-.2,.2, 0.1); surf.tan <- expand.grid(grid.tan, grid.tan)
z.tan <- matrix(-2*surf.tan[,2]+7.35*surf.tan[,1]+1.915, 
                5, 5)

### Plotting
m <- list(
  l = 0,
  r = 0,
  b = 0,
  t = 0,
  pad = 0
)

fig <- plot_ly(width = 800, height = 597) %>%
       layout(title="", margin=list(top=0),
         scene = list(xaxis = list(title = 't=1', showticklabels=F, showgrid=F, zeroline=F, showline=T),
                      yaxis = list(title = 't=2', showticklabels=F, showgrid=F, zeroline=F, showline=T),
                      zaxis = list(title = 't=3', showticklabels=F, showgrid=F, zeroline=F, showline=T),
                      camera = list(eye = list(x = 2.1, y = -.9, z = 0.2))))
fig$sizingPolicy$padding <- "0"


fig3 <- fig %>%
  add_surface(x=grid, y=grid, z=z, type="surface", opacity = .5, 
              colorscale = list(c(0,1),c("grey","grey")),
              #colors=c("grey", "blue"),
              showscale=FALSE) %>%
  add_trace(x=data$y1, y=data$y2, z=data$y3, type="scatter3d", mode="markers", 
            marker=list(color="black", size=2), opacity=0.5, showlegend=FALSE)

# Figure 1 in the main paper
a0 <- c(0.3, 0.7, 1.915)
fig4 <- fig3 %>%
  add_trace(x=a0[1], y=a0[2], z=a0[3], type="scatter3d", mode="markers", 
            marker=list(color="red", size=4), opacity=1, showlegend=FALSE)


# Figure 2 in the main paper
fig5 <- fig4 %>%
  add_surface(x=grid.tan+0.3, y=grid.tan+0.7, z=z.tan, type="surface", opacity = .5,
              colorscale = list(c(0,1),c("navy","navy")),
              #colors=c("blue"),
              showscale=FALSE)

