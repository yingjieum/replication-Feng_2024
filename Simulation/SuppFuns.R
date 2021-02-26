###########################################
### Define supporting functions ###########
######## for simulation ###################
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
hdfun1 <- function(u,v) {
    #out <- (u^2+v^2)/(3*cos(1/(u^2+v^2)))+0.15
    #out <- sin(5*pi*(u+v-1)+1)/2+0.5
    out <- sin(pi*(u+v))
  return(out)
}

hdfun2 <- function(u,v) {
    out <- (u-v)^2
  return(out)
}
funlist <- list(hdfun1=hdfun1, hdfun2=hdfun2)


# DGP generator
dgp <- function(n, p, hdmodel, addz=FALSE) {
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
    x <- outer(alpha, eta, FUN = hdmodel) + matrix(rnorm(n*p),n,p)   # n by p matrix
    
    return(list(y=y, y0=Y0, y1=Y1, x=x, z=z, alpha=alpha, eta=eta, d=D))
}

################################################
# Searching KNN
findknn <- function(v, A2, K) {
     tmp   <- abs(A2 - v)
     tmp.d <- colMaxs(tmp, value = T)
     out   <- ifelse(tmp.d <= nth(tmp.d, K), TRUE, FALSE) 
     return(out)
     #Kmat[i, which(tmp.d <= nth(tmp.d, K))] <- TRUE
}

knn.index <- function(A, K=NULL) {    # A: n by p; # K: tuning param.
    p    <- ncol(A)
    A2   <- tcrossprod(A)/p
    Kmat <- apply(A2, 2, function(v) findknn(v=v, A2=A2, K=K))
    return(Kmat)           # neighborhood saved in each column
}

# Local PCA (for each point)
# A: n by p, input matrix; index: n by 1, KNN index; nlam: number of vectors to be extracted
lpca <- function(A, index, nlam, n) {      
    A <- cbind(1:n, A)   # add a row number
    A <- A[index,]
    eigvec <- irlba(A[,-1], nu=nlam)$u
    return(list(Lam=eigvec, no.neigh=A[,1]))    #Lam: K by nlam
}


# Estimation
# i: evaluate at ith obs. (order as in the whole sample); y: outcome of interest; 
# LAM: lambda matrix around i (NOT depend on d)!!; z: additional controls
# kmat: index matrix, relative to the whole sample
# subset: additional subsetting, relative to the whole sample
# i.Lam: PCs around i; i.no.neigh: no. of units around i; BOTH from lpca, of length K
linkinv.d <- binomial(link="logit")$linkinv

pred <- function(i, y, d, x, subset=NULL, kmat, nlam, n, p) {   
    index <- kmat[,i]   # length n
    PC    <- lpca(A=x[,(p/2+1):p], index=index, nlam=nlam, n=n)

    # prepare design
    sub    <- index & subset   # of length n
    y.sub  <- y[sub]
    design <- Lam.sub <- PC$Lam[sub[index],,drop=F]
        
    # eval
    design.d <- PC$Lam
    eval <- design.d[which(PC$no.neigh==i),]
    
    # # run a local GLM of y at ith obs.
    # model    <- glm.fit(x=design, y=y.sub, family=family.y)
    # y.fitval <- model$family$linkinv(sum(model$coeff * eval))
    y.fitval <- sum(.lm.fit(design, y.sub)$coeff * eval)
    
    # run a local GLM of d at ith obs.
    # model    <- glm.fit(x=design.d, y=d[index], family=family.d)
    # d.fitval <- model$family$linkinv(sum(model$coeff * eval))
    d.fitval <- linkinv.d(sum(fastLR(design.d, d[index])$coeff * eval))
    
    return(c(y.fitval, d.fitval))   # two scalars
}

# Pred, only for local constant approx.
pred.cons <- function(i, y, d, subset=NULL, kmat) {   
    index  <- kmat[,i]   # length n
    sub    <- index & subset   # of length n
    y.fitval <- mean(y[sub])
    d.fitval <- mean(d[index])
    
    return(c(y.fitval, d.fitval))   # two scalars
}


# Computation
compute <- function(range, y, d, x, subset=NULL, K, nlam, n, p, const=F) {
    if (const) {
      kmat <- knn.index(A=x, K=K)
      fit  <- sapply(range, function(i) pred.cons(i=i, y=y, d=d, subset=subset,
                                           kmat=kmat)) 
    } else {
      kmat <- knn.index(A=x[,1:(p/2)], K=K)
      fit  <- sapply(range, function(i) pred(i=i, y=y, d=d, x=x, subset=subset,
                                           kmat=kmat, nlam=nlam, n=n, p=p))
    }
    return(list(yfit=fit[1,], ps=fit[2,]))
}

# apply across K
fitting <- function(range, y, d, x, subset=NULL, K, nlam, n, p, pr, theta0, const=F) {
    result <- compute(range=range, y=y, d=d, x=x, subset=subset, 
                      K=K, nlam=nlam, n=n, p=p, const=const)
    # calculate counterfactual mean
    psi <- (d * result$yfit + (1-d)*(y-result$yfit)*result$ps/(1-result$ps))/pr #NOT influence fun
    theta <- mean(psi)   # point estimate
    
    varphi <- (d * (result$yfit-theta) + (1-d)*(y-result$yfit)*result$ps/(1-result$ps))/pr
    se    <- sqrt(mean(varphi^2))/sqrt(n)
    
    rej <- (((theta-theta0)/se > qnorm(0.975)) | ((theta-theta0)/se < qnorm(0.025))) * 1
    ci  <- se*qnorm(0.975)*2
    return(c(theta, se, rej, ci))
}

# sim function
sim <- function(i, n, p, model, Kseq, nlam, theta0, const=const) {
    data   <- dgp(n=n, p=p, hdmodel=model)
    subset <- (data$d==0)   # so compute mu.0 in the following
    pr     <- mean(data$d)
    range  <- 1:n
    
    result <- sapply(Kseq, function(K) fitting(range=range, y=data$y, d=data$d, x=data$x, subset=subset, 
                                               K=K, 
                                               nlam=nlam, n=n, p=p, pr=pr, theta0=theta0, const=const))

    return(result)  # 4 by length(Kseq)
}