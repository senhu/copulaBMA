utdc.empcop.1 <- function(X){
  n <- nrow(X)
  t <- sqrt(n)
  c.approx <- copula::C.n(matrix(c(1-t/n, 1-t/n),nrow=1), X)
  utdc <- (1 - 2*(1-t/n) + c.approx)/(t/n)
  return(utdc)
}

# utdc.empcop.1(pobs(X))

ltdc.empcop.1 <- function(X){
  n <- nrow(X)
  t <- sqrt(n)
  c.approx <- copula::C.n(matrix(c(t/n, t/n),nrow=1), X, smoothing = "none")
  ltdc <- c.approx/(t/n)
  return(ltdc)
}

# ltdc.empcop.1(pobs(X))

#--------------------------------------------

utdc.empcop.2 <- function(X){
  n <- nrow(X)
  t <- sqrt(n)
  c.approx <- copula::C.n(matrix(c(1-t/n, 1-t/n),nrow=1), X, smoothing = "none")
  utdc <- (2 - log(c.approx)/log(1-t/n))
  return(utdc)
}
# utdc.empcop.2(pobs(X))

ltdc.empcop.2 <- function(X){
  n <- nrow(X)
  t <- sqrt(n)
  c.approx <- copula::C.n(matrix(c(t/n, t/n),nrow=1), X, smoothing = "none")
  ltdc <- 2 - log(1 - 2*(t/n) + c.approx)/log(1-t/n)
  return(ltdc)
}

# ltdc.empcop.2(pobs(X))

#------------------------------------

utdc.empcop.3 <- function(X, plot=FALSE){
  n <- nrow(X)
  k <- round(sqrt(n))
  Uemp <- function(i) {1 - 2*(1-i/n) + copula::C.n(matrix(c(1-i/n, 1-i/n),nrow=1), X, smoothing = "none")}
  emp.lambda <- Vectorize(Uemp)(k:1)
  xx <- (k:1)/n
  lmod <- lm(emp.lambda~(xx-1))
  if (plot) {
    plot(xx, emp.lambda, type="l", ylab="upper tail dependence")
    lines(xx, lmod$fitted.values, col=2)
  }
  return(sum(coef(lmod)))
}

ltdc.empcop.3 <- function(X, plot=FALSE){
  n <- nrow(X)
  k <- round(sqrt(n))
  Lemp <- function(i) {copula::C.n(matrix(c(i/n, i/n),nrow=1), X, smoothing = "none")}
  emp.lambda <- Vectorize(Lemp)(1:k)
  xx <- (1:k)/n
  lmod <- lm(emp.lambda~(xx-1))
  if (plot) {
    plot(xx, emp.lambda, type="l", ylab="lower tail dependence")
    lines(xx, lmod$fitted.values, col=2)
  }
  return(sum(coef(lmod)))
}

#------------------------------------

utdc.empcop.4 <- function(X, 
                          bootstrap = FALSE, bootstrap.num = 100, 
                          smooth = TRUE, plot = TRUE,
                          cutoff = 0.9, smoothing.span = 0.1) {
  n <- nrow(X)
  
  if (!bootstrap){
    U=X[order(X[,1]),1]
    V=X[order(X[,2]),2]
    Emp <- function(i, j) {sum(X[,1]<=U[i] & X[,2]<=V[j])/n}
    Uemp <- function(i) (1-(2*i/n)+Emp(i, i))/(1-(i/n))
    emp.lamb <- Vectorize(Uemp)(1:(n-1))
    if (!smooth){
      wm <- which.min(emp.lamb[1:((n-1)*cutoff)])
      m  <- min(emp.lamb[1:((n-1)*cutoff)]) 
    }
    if (smooth) {
      con <- 1:(n-1)
      emp.lamb.ss <- predict(loess(emp.lamb ~ con, span = smoothing.span))
      wm <- which.min(emp.lamb.ss[1:((n-1)*cutoff)])
      m <- min(emp.lamb.ss[1:((n-1)*cutoff)])
    }
    if (plot){
      plot(emp.lamb, ylim=c(0,1), type="l", ylab="upper tail dependence")
      if (smooth) lines(emp.lamb.ss, col=2)
      abline(h=m, v=wm)
    }
  }
  
  if (bootstrap){
    boot.res <- NULL
    for (k in 1:bootstrap.num){
      Z <- X[sample(1:n, n, replace=TRUE),]
      U.z=Z[order(Z[,1]),1]
      V.z=Z[order(Z[,2]),2]
      Emp.boot <- function(i, j) {sum(Z[,1]<=U.z[i] & Z[,2]<=V.z[j])/n}
      Uemp.boot <- function(i) (1-(2*i/n)+Emp.boot(i, i))/(1-(i/n))
      boot.res <- rbind(boot.res, Vectorize(Uemp.boot)(1:(n-1)))
    }
    lm <- apply(boot.res, 2, mean)
    lq <- apply(boot.res[,1:(n-1)], 2, quantile, .025)
    uq <- apply(boot.res[,1:(n-1)], 2, quantile, .975)
    if (!smooth){
      wm <- which.min(lm[1:((n-1)*cutoff)])
      m  <- min(lm[1:((n-1)*cutoff)]) 
    }
    if (smooth) {
      con <- 1:(n-1)
      emp.lamb.ss <- predict(loess(lm ~ con, span = smoothing.span))
      wm <- which.min(emp.lamb.ss[1:((n-1)*cutoff)])
      m <- min(emp.lamb.ss[1:((n-1)*cutoff)])
    }
    if (plot){
      plot(boot.res[1,], ylim=c(0,1), type="l", ylab="upper tail dependence")
      for (k in 2:100) lines(boot.res[k,])
      lines(lm, col=2, lwd=2)
      lines(lq, col=2)
      lines(uq, col=2)
      if (smooth) lines(emp.lamb.ss, col=4, lwd=2)
      abline(h=m, v=wm)
    }
  }
  
  return(c(m, wm))
}

# utdc.empcop.4(X, cutoff = .99, bootstrap = T)

ltdc.empcop.4 <- function(X, 
                          bootstrap = FALSE, bootstrap.num = 100,
                          smooth = TRUE, plot = TRUE,
                          cutoff = 0.9, smoothing.span = 0.1) {
  
  n <- dim(X)[1]
  
  if (!bootstrap){
    U=X[order(X[,1]),1]
    V=X[order(X[,2]),2]
    Emp <- function(i, j) {sum(X[,1]<=U[i] & X[,2]<=V[j])/n}
    Lemp <- function(i) Emp(i, i)/(i/n)
    emp.lamb <- Vectorize(Lemp)(1:(n-1))
    if (!smooth){
      wm <- which.min(emp.lamb[((n-1)*(1-cutoff)):(n-1)]) + (n-1)*(1-cutoff)
      m  <- min(emp.lamb[((n-1)*(1-cutoff)):(n-1)])
    }
    if (smooth) {
      con <- 1:(n-1)
      emp.lamb.ss <- predict(loess(emp.lamb ~ con, span = smoothing.span))
      wm <- which.min(emp.lamb.ss[((n-1)*(1-cutoff)):(n-1)]) + (n-1)*(1-cutoff)
      m <- min(emp.lamb.ss[((n-1)*(1-cutoff)):(n-1)])
    }
    if (plot){
      plot(emp.lamb, ylim=c(0,1), type="l", ylab="upper tail dependence")
      if (smooth) lines(emp.lamb.ss, col=2)
      abline(h=m, v=wm)
    }
  }
  
  if (bootstrap){
    boot.res <- NULL
    for (k in 1:bootstrap.num){
      Z <- X[sample(1:n, n, replace=TRUE),]
      U.z=Z[order(Z[,1]),1]
      V.z=Z[order(Z[,2]),2]
      Emp.boot <- function(i, j) {sum(Z[,1]<=U.z[i] & Z[,2]<=V.z[j])/n}
      Lemp.boot <- function(i) Emp.boot(i, i)/(i/n)
      boot.res <- rbind(boot.res, Vectorize(Lemp.boot)(1:(n-1)))
    }
    lm <- apply(boot.res, 2, mean)
    lq <- apply(boot.res[,1:(n-1)], 2, quantile, .025)
    uq <- apply(boot.res[,1:(n-1)], 2, quantile, .975)
    if (!smooth){
      wm <- which.min(lm[((n-1)*(1-cutoff)):(n-1)]) + (n-1)*(1-cutoff)
      m  <- min(lm[((n-1)*(1-cutoff)):(n-1)])
    }
    if (smooth) {
      con <- 1:(n-1)
      emp.lamb.ss <- predict(loess(lm ~ con, span = smoothing.span))
      wm <- which.min(emp.lamb.ss[((n-1)*(1-cutoff)):(n-1)]) + (n-1)*(1-cutoff)
      m <- min(emp.lamb.ss[((n-1)*(1-cutoff)):(n-1)])
    }
    if (plot){
      plot(boot.res[1,], ylim=c(0,1), type="l", ylab="lower tail dependence")
      for (k in 2:bootstrap.num) lines(boot.res[k,])
      lines(lm, col=2, lwd=2)
      lines(lq, col=2)
      lines(uq, col=2)
      if (smooth) lines(emp.lamb.ss, col=4, lwd=2)
      abline(h=m, v=wm)
    }
  }
  
  return(c(m, wm))
}

# ltdc.empcop.4(X, cutoff = .99, bootstrap = T)
