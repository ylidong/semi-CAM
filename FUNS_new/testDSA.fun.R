
CAM.nWCA_force.new<-function(data, K, cluster_num,cluster_arg,MKS.initial) {
  #data=Data.sus;K=k;cluster_num=50;cluster_arg="Hartigan-Wong";MKS.initial=MKS.list.given
  #data=data_cl;K=k;cluster_num=50;cluster_arg="Hartigan-Wong";MKS.initial=mk.list_in
  
  #data_ssCAM=Data.cls;data_est=Data;ncell;cluster_num=50;mks_in=MKS.initial
  #data=data_ssCAM
  
  #data=data_ssCAM;K=ncell;cluster_num=cluster_num;cluster_arg="Hartigan-Wong";mks_in=MKS.initial
  
  data=as.matrix(data)
  
  # data=as.matrix(data.mix.filter);K=ncell;cluster_num=50;cluster_arg="Hartigan-Wong";MKS.initial=MKS.initial
  
  #data=data_filter;K=ncell;cluster_num=20;cluster_arg="Hartigan-Wong";MKS.initial=MKS.initial
  #data=data_filter_all;K=cn_ni;cluster_num=20;cluster_arg="Hartigan-Wong";MKS.initial=mks$MKS.cam
  
  # load the library "MASS" for the function "ginv"
  packageExist <- require("MASS")
  if (!packageExist) {
    install.packages("MASS")
    library("MASS")
  }
  
  # load the library "geometry" for the function "convhulln"
  packageExist <- require("geometry")
  if (!packageExist) {
    install.packages("geometry")
    library("geometry")
  }
  
  # load the library "nnls" for the function "nnls"
  packageExist <- require("nnls")
  if (!packageExist) {
    install.packages("nnls")
    library("nnls")
  }
  
  data_size<- dim(data)[1]
  L <- dim(data)[2]
  
  ##### use kmeans to cluster the observation #####
  cat("\nPerforming kmeans to cluster data ... \n")
  
  data_norm=data/rowSums(data)
  
  ##calculate given initial centers given markers##
  
  if(length(MKS.initial)!=0){
    center.mks=do.call(rbind,lapply(MKS.initial,function(mks){ if(length(mks)>1) out=colMeans(data_norm[mks,]);if(length(mks)==1) out=data_norm[mks,];out} ) )
    rownames(center.mks)=paste0("MKS_kown",1:nrow(center.mks))
  }else{
    center.mks=NULL
  }
  
  tt = lapply(MKS.initial, length)
  mm = length(tt)
  
  data_new = data_norm
  
  if(mm>0){
    for(k in 1:mm){
      data_new[rownames(data_new)%in%unlist(MKS.initial[k]),] <- matrix(rep(center.mks[k,], times= length(unlist(MKS.initial[k]))), nrow=length(unlist(MKS.initial[k])), ncol=length(center.mks[k,]), byrow = T)
    }
  }
  ##call kmeans##
  cluster <- kmeans(data_new,cluster_num,iter.max=100,algorithm =cluster_arg)
  
  for (i in 1:50){
    
    tmp <- kmeans(data_new,cluster_num,iter.max=100,algorithm =cluster_arg)
    
    if (cluster$tot.withinss>tmp$tot.withinss){
      cluster <- tmp
    }
    #print(cluster$tot.withinss)
  }
  
  ##exclude small clusters##
  small_cluster <- matrix(numeric(0),0,0)
  for (k in 1:cluster_num){
    if (cluster$size[k]<0.1*dim(data)[1]/cluster_num)
      small_cluster <- c(small_cluster,k)
  }
  
  if (length(small_cluster)==0){
    cluster <- cluster
  } else {
    cluster$centers <- cluster$centers[-small_cluster,]
  }
  
  Vs=lapply(rownames(cluster$centers),function(x) names(cluster$cluster)[which(cluster$cluster==x)])
  names(Vs) <- rownames(cluster$centers)
  
  ##remember we use location here not the real names!##
  vertice.known.nm=unlist(lapply(MKS.initial,function(mks) names(Vs)[which(unlist(lapply(Vs,function(x) sum(mks%in%x) ))==length(mks))] ) )
  vertice.known=unlist(lapply(MKS.initial,function(mks) which(unlist(lapply(Vs,function(x) sum(mks%in%x) ))==length(mks))  ))
  
  if(length(vertice.known)!=length(MKS.initial)) warning("Cannot find correct known vertice!")
  
  if (K==length(MKS.initial)){##if we known K markers,the K cluster will be the final vertice##
    #A_est <- t(cluster$centers[as.character(vertice.known),])
    #V.pick=lapply(vertice.known,function(kk) Vs[[as.character(kk)]])
    A_est <- t(cluster$centers[vertice.known,])
    V.pick=lapply(vertice.known,function(kk) Vs[[kk]])
    
  }else{
    
    J <- dim(cluster$centers)[1]
    ##convex hull method to choose potential vertice##
    convex <- convhulln(rbind(cluster$centers,0))
    
    corner <- matrix(numeric(0), 0,0)
    for (i in 1:L){
      corner <- union(corner,convex[,i])
    }
    
    for (j in 1:length(corner)){
      if (corner[j]==(J+1)){
        break
      }
    }
    corner <- corner[-j]  # throw away the origin point
    
    ##force the known cell type to be a corner##
    corner=union(corner,vertice.known)
    
    J_out <- length(corner)
    ##### estimate A and S ###########
    cat("\nEstimating A and S ... \n")  
    
    ##get K final vertices##
    if (K==J_out){
      A_est <- t(cluster$centers[corner,])
      vertice=colnames(A_est)
    }else{
      
      if(length(MKS.initial)!=0){# this step force the given markers have to be part of final vertices#
        #mks.known=as.character(seq(1,length(MKS.initial)))
        cornerResult <- measure_conv.miss(t(cluster$centers[corner,]),K,vertice.known.nm)
      }else{
        cornerResult <- measure_conv(t(cluster$centers[corner,]),K)
      }
      
      A_est <- cornerResult[[1]]
      col.nm=colnames(A_est)
      ind <- cornerResult[[2]] 
      vertice=colnames(A_est)
    }
    
    V.pick=lapply(1:K,function(x) names(cluster$cluster[cluster$cluster== vertice[x] ]) )
  }
  
  
  ##calculate proportions##
  #scale <- ginv(A_est)%*%matrix(1,L,1)
  solve=nnls(A_est, matrix(1,nrow=dim(A_est)[1],ncol=1))
  scale <- solve$x
  err=solve$deviance
  
  scale <- as.vector(scale)
  A_est <- A_est%*%diag(scale)
  colnames(A_est)=col.nm
  
  S_est <- matrix(0,nrow=dim(A_est)[2],ncol=dim(data)[1])
  for (i in 1:nrow(data)){
    S_est[,i] <- coef(nnls(A_est,t(data)[,i]))
  }
  colnames(S_est)=rownames(data)
  #V.pick=lapply(1:K,function(x) names(cluster$cluster[cluster$cluster== vertice[x] ]) )
  
  return(list(V.pick=V.pick,Vs=Vs,A_est=A_est,S_est=t(S_est),err=err,vertice.known.nm=vertice.known.nm))
  #V.pick: final K vertices;Vs: cluster_num clusters after Kmeans;A_est: proportion estimates;S_est: cell type specific expression estimates
}


ssCAM.main.NMF.1run.new<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;
  
  ssCAM=CAM.nWCA_force.new(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  ct.nm=which(colnames(ssCAM$A_est)%in%ssCAM$vertice.known.nm)
  
  NMF.KL=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssKL" ,log = FALSE, nrun = 1)
  P.NMF.KL=NMF.KL@fit@H
  X.NMF.KL=NMF.KL@fit@W
  P.NMF.KL.know=P.NMF.KL[ct.nm,]
  
  NMF.Fro=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssFrobenius" ,log = FALSE, nrun = 1)
  P.NMF.Fro=NMF.Fro@fit@H
  X.NMF.Fro=NMF.Fro@fit@W
  P.NMF.Fro.know=P.NMF.Fro[ct.nm,]
  
  
  return(list("P.KL"=P.NMF.KL,"P.KL.know"=P.NMF.KL.know,"X.KL"=X.NMF.KL,"P.Fro"=P.NMF.Fro,"P.Fro.know"=P.NMF.Fro.know,"X.Fro"=X.NMF.Fro,"MKS"=ssCAM$V.pick))
  
}