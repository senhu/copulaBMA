L2dist <- function(X, Y){
  if (nrow(X)!=nrow(Y)) {stop("X and Y do not have the same dimensions")}
  n <- nrow(X)
  dis <- rep(0,n)
  for (t in (1:n)){
    c1 <- copula::C.n(matrix(c(t/n, t/n),nrow=1), copula::pobs(X), smoothing = "none")
    c2 <- copula::C.n(matrix(c(t/n, t/n),nrow=1), copula::pobs(Y), smoothing = "none")
    dis[t] <- (c1-c2)^2
  }
  return(sqrt(sum(dis)))
}
# L2dist(h1, h2)


# dis <- matrix(0,nrow=n, ncol=n)
# for (t in (1:n)){
#   for (k in (1:n)){
#     c1 <- copula::C.n(matrix(c(t/n, k/n),nrow=1), h1, smoothing = "none")
#     c2 <- copula::C.n(matrix(c(t/n, k/n),nrow=1), h2, smoothing = "none")
#     dis[t,k] <- (c1-c2)^2
#     print(c(t,k))
#   }
# }
# sum(dis)
