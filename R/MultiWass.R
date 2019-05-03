require(transport)
MultiWass <- function(act, pred, prob=FALSE){
  x     <- pp(as.matrix(act))
  y     <- pp(as.matrix(pred))
  # match <- transport(x,y,p=1)
  # wass  <- wasserstein(x,y,p=1,match)
  if (prob==FALSE) wass <- wasserstein(x,y,p=1, prob=FALSE)
  if (prob) wass <- wasserstein(x,y,p=1, prob=TRUE)
  return(wass)
}
