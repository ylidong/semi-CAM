# This software is under license BSD 3.

measure_conv.miss <- function(X,N,mks.known){
  
  #X=t(cluster$centers[corner,]);N=K
  # This function needs the "nnls" package in R, you can use the command: 
  # install.packages("nnls")
  # in R interface to install that package if you don't have it.
  require("nnls")
  
  #X=t(cluster$centers[corner,]);N=K;mks.known=vertice.known.nm
  
  M <- nrow(X); L <- ncol(X)
  total_index <- t(combn(1:L,N))
  
  mks.known.ind=match(mks.known,colnames(X))
  comb.select=which(apply(total_index,1,function(x) sum(x%in%mks.known.ind))==length(mks.known))
  
  total_index_sub=total_index[comb.select,]
  error <- matrix(0,length(comb.select),1)
  #browser()
  for (p in 1:nrow(total_index_sub)) {
    #p=1
    A <- X[,total_index_sub[p,]]
    Others <- X[,-total_index_sub[p,]]
    Ae <- rbind(1e-5*A, matrix(1,1,N))
    Oe <- rbind(1e-5*Others, matrix(1,1,L-N))
    #Oe <- rbind(as.matrix(1e-5*Others), matrix(1,1,L-N)) ##add 1's row to maker sure sumof alpha=1##
    #browser()
    alpha <- matrix(0,N,L-N)
    for (i in 1:(L-N))
      alpha[,i] <- nnls(Ae, Oe[,i])$x
    error[p] <- norm(Others - A %*% alpha, 'f')^2    
  }
  #browser()
  val <- min(error)
  ind <- which.min(error)
  
  eA <- X[,total_index_sub[ind,]]
  #eA <- eA[,rank(colnames(eA))]
  cornerind <- total_index_sub[ind,]
  return(list(eA,cornerind))
}
