# This software is under license BSD 3.

measure_conv <- function(X,N){
  # This function needs the "nnls" package in R, you can use the command: 
  # install.packages("nnls")
  # in R interface to install that package if you don't have it.
  require("nnls")
  #X=t(cluster$centers[corner,]);N=K
  M <- nrow(X); L <- ncol(X)
  total_index <- t(combn(1:L,N))
  comb <- nrow(total_index)
  
  error <- matrix(0,comb,1)
  #browser()
  for (p in 1:comb) {
    A <- X[,total_index[p,]]
    Others <- X[,-total_index[p,]]
    Ae <- rbind(1e-5*A, matrix(1,1,N))
    Oe <- rbind(1e-5*Others, matrix(1,1,L-N))
    #browser()
    alpha <- matrix(0,N,L-N)
    for (i in 1:(L-N))
      alpha[,i] <- nnls(Ae, Oe[,i])$x
    error[p] <- norm(Others - A %*% alpha, 'f')^2    
  }
  #browser()
  val <- min(error)
  ind <- which.min(error)
  
  eA <- X[,total_index[ind,]]
  cornerind <- total_index[ind,]
  return(list(eA,cornerind))
}