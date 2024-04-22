###########################################
### Define supporting functions ###########
####### for simulation ####################
####### Last updated: Apr 22, 2024 ########
###########################################

# treatment model
mu0 <- function(alpha) {
    return(alpha + alpha^2)
}
mu1 <- function(alpha) {
    return(2*alpha + alpha^2 + 1)
}
px <- function(alpha) {
  alpha <- alpha - 0.5
  return(exp(alpha+alpha^2)/(1 + exp(alpha+alpha^2)))
}

# large-dim x
# 1 and 2 are old and not used anymore.
hdfun1 <- function(u,v) {
    out <- sin(pi*(u+v))
    return(out)
}

hdfun2 <- function(u,v) {
    out <- (u-v)^2
    return(out)
}

hdfun3 <- function(u, v, sig=.1) {
  out <- (1/sqrt(2*pi)/sig) * exp(-(u-v)^2/sig^2)
  return(out)
}

hdfun4 <- function(u, v, sig=.1) {
  out <- exp(-abs(u-v)/sig)
  return(out)
}

funlist <- list(hdfun1=hdfun1, hdfun2=hdfun2, hdfun3=hdfun3, hdfun4=hdfun4)


# DGP generator
# sigma: s.d. for error in x
dgp <- function(n, p, hdmodel, addz=FALSE, sigma=1) {
    # latent variable
    alpha <- runif(n, 0, 1)
    
    z <- NULL
    if (addz) {
        z <- runif(n, 0, 1)
    }
    
    # outcome
    Y0 <- mu0(alpha) + rnorm(n)
    Y1 <- mu1(alpha) + rnorm(n)
    
    # logistic assignment
    prob  <- px(alpha)
    runis <- runif(n,0,1)
    D     <- ifelse(runis < prob,1,0)
    
    y <- D*Y1+(1-D)*Y0
    
    # HD covariates
    eta <- runif(p, 0, 1)
    x <- outer(alpha, eta, FUN = hdmodel) + matrix(rnorm(n*p, 0, sigma), n, p)   # n by p matrix
    
    return(list(y=y, y0=Y0, y1=Y1, x=x, z=z, alpha=alpha, eta=eta, d=D))
}

################################################
# Note the following KNN can accept a vector of Ks and report a list of matrices
# Searching KNN
findknn <- function(v, A2, Kseq=NULL) {
     tmp   <- abs(A2 - v[-1])
     diag(tmp)  <- -1        # implicitly delete diagonal elements
     tmp[v[1],] <- -1        # implicitly delete the v[1]th element
     tmp.d <- colMaxs(tmp, value = T)
     out   <- sapply(Kseq, function(K) ifelse(tmp.d <= nth(tmp.d, K), TRUE, FALSE)) # n by length(Kseq) matrix
     return(out)
}

knn.index <- function(A, Kseq=NULL) {    # A: n by p; # Kseq: vector of tuning param.
    p    <- ncol(A)
    A2   <- tcrossprod(A)/p
    Kmat <- apply(cbind(1:nrow(A), A2), 1, function(v) findknn(v=v, A2=A2, Kseq=Kseq))
    return(Kmat)           # a matrix: neighborhood for each unit saved in each column, each n rows -> one choice of K in Kseq 
}

# select # of local factors, n should >=2
sel.r <- function(sval, n) {
  d <- length(sval)
  ratio <- (sval[-d] / sval[-1]) < log(log(n))
  if (any(ratio)) {
    r <- max(which.max(ratio)-1, 1)   # at least 1 factor
  } else {
    r <- d - 1
  }
  return(r)
}

# Local PCA (for each point)
# A: n by p, input matrix; index: n by 1, KNN index; nlam: Maximum number of vectors to be extracted
lpca <- function(A, index, nlam, n, K) {      
    A <- cbind(1:n, A)   # add a row number
    A <- A[index,]
    svd <- irlba(A[,-1], nu=nlam, nv=nlam)
    r <- sel.r(svd$d, K)
    return(list(Lam=svd$u[,1:r,drop=F], no.neigh=A[,1], sv=svd$d, d.i=r))    #Lam: K by r
}


###### Estimation ####################
# i: evaluate at ith obs. (order as in the whole sample); y: outcome of interest; 
# LAM: lambda matrix around i (NOT depend on d)!!; z: additional controls
# fullkmat: index matrix, relative to the whole sample, contains results for all K's
# subset: additional subsetting, relative to the whole sample
# i.Lam: PCs around i; i.no.neigh: no. of units around i; BOTH from lpca, of length K
linkinv.d <- binomial(link="logit")$linkinv

# local regression
locreg <- function(i, y, d, x, subset=NULL, kmat, nlam, n, p, shutdown.d=F, K) {
  y.fitval <- d.fitval <- NA
  index <- kmat[,i]   # length n
  PC    <- lpca(A=x[,(p/2+1):p], index=index, nlam=nlam, n=n, K=K)
  
  # prepare design
  if (is.null(subset)) sub <- index
  else                 sub <- index & subset   # of length n
  y.sub  <- y[sub]
  design <- Lam.sub <- PC$Lam[sub[index],,drop=F]
  
  # eval
  design.d <- PC$Lam
  eval <- design.d[which(PC$no.neigh==i),]
  
  # run a local GLM of y at the i-th obs.
  y.fitval <- sum(.lm.fit(design, y.sub)$coeff * eval)
  
  if (!shutdown.d) {
    # run a local GLM of d at ith obs.
    d.fitval <- linkinv.d(sum(fastLR(design.d, d[index])$coeff * eval))
  }
  
  return(c(y.fitval, d.fitval))   # two scalars
}

# local regression, only for local constant approx.
locreg.cons <- function(i, y, d, subset=NULL, kmat, shutdown.d=F) {   
  y.fitval <- d.fitval <- NA  
  index  <- kmat[,i]   # length n
  if (is.null(subset)) sub <- index
  else                 sub <- index & subset   # of length n
  y.fitval <- mean(y[sub])
  if (!shutdown.d) {
    d.fitval <- mean(d[index])
  }
  
  return(c(y.fitval, d.fitval))   # two scalars
}

# Prediction, for a sequence of K's
pred <- function(i, y, d, x, subset=NULL, fullkmat, nlam, n, p, shutdown.d=F, const=F, Kseq) {   
    L <- nrow(fullkmat)/n
    
    if (const) {
       out <- sapply(c(1:L), function(l) locreg.cons(i, y, d, subset, kmat=fullkmat[((l-1)*n+1):(l*n),], shutdown.d))
    } else {
       out <- sapply(c(1:L), function(l) locreg(i, y, d, x, subset, kmat=fullkmat[((l-1)*n+1):(l*n),], nlam, n, p, shutdown.d, Kseq[l]))
    }
    
    return(out)   # 2 by Kseq matrix or a vector
}


# Computation, for all units in range
compute <- function(range, y, d, x, subset=NULL, Kseq, nlam, n, p, const=F, shutdown.d=F) {
    fullkmat <- knn.index(A=x[,1:(p/2)], Kseq=Kseq)
    fit  <- sapply(range, function(i) pred(i=i, y=y, d=d, x=x, subset=subset,
                                           fullkmat=fullkmat, nlam=nlam, n=n, p=p, 
                                           shutdown.d=shutdown.d, const=const, Kseq=Kseq))
    no <- rep(1:2, length(Kseq))
    
    return(list(yfit=fit[no==1,,drop=F], ps=fit[no==2,,drop=F]))  # each is a length(Kseq) by length(range) matrix
}

#################################################

# calculate stats related to treatment effect
te.stats <- function(y, d, pr, theta0, result, l) {
    yfit <- result$yfit[l,]; ps <- result$ps[l,]
    
    psi <- (d * yfit + (1-d)*(y-yfit)*ps/(1-ps))/pr # NOT influence fun
    theta <- mean(psi)   # point estimate
    
    varphi <- (d * (yfit-theta) + (1-d)*(y-yfit)*ps/(1-ps))/pr
    se    <- sqrt(mean(varphi^2))/sqrt(n)
    
    rej <- (((theta-theta0)/se > qnorm(0.975)) | ((theta-theta0)/se < qnorm(0.025))) * 1
    ci  <- se*qnorm(0.975)*2
    
    return(c(theta, se, rej, ci))
}

# apply across K's, calculate necessary quantities
fitting <- function(range, y, d, x, subset=NULL, Kseq, nlam, n, p, pr, theta0, const=F) {
  result <- compute(range=range, y=y, d=d, x=x, subset=subset, 
                    Kseq=Kseq, nlam=nlam, n=n, p=p, const=const)
  
  # calculate counterfactual mean
  out <- sapply(1:length(Kseq), function(l) te.stats(y, d, pr, theta0, result, l=l))

  return(out)   # 4 by Kseq matrix
}

# sim function
sim <- function(i, n, p, model, Kseq, nlam, theta0, const) {
  data   <- dgp(n=n, p=p, hdmodel=model)
  subset <- (data$d==0)   # so compute mu.0 in the following
  pr     <- mean(data$d)
  range  <- 1:n
  
  # 20-fold CV choice
  Klist <- seq(Kseq[1], Kseq[length(Kseq)], 5)
  K0  <- Klist[which.min(knn.cv(nfolds=20, y=data$y, x=data$x[,p,drop=F], k=Klist, type="R")$crit)]
  
  # fitting
  Kseq <- c(Kseq, K0)   # add the DPI choice
  result <- fitting(range=range, y=data$y, d=data$d, x=data$x, subset=subset, 
                    Kseq=Kseq, nlam=nlam, n=n, p=p, pr=pr, theta0=theta0, const=const)
  
  result <- rbind(Kseq, result)
  
  return(result)  # (4+1) by length(Kseq)+1 matrix
}

