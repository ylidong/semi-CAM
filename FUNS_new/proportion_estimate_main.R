CAM.main.est<-function(data_mix,ncell,cluster_num=20){

  CAM=CAM.nWCA(X=as.matrix(t(data_mix)), K=ncell,cluster_num=cluster_num)
  P.CAM=t(CAM$A_est)
  P.CAM=scale(P.CAM, center=F, scale= apply(P.CAM, 2, sum))
  
  
  return(list("P"=P.CAM,"MKS"=CAM$V.pick))
  
}


CAM.main.est.DSA<-function(data_mix,data_est,ncell,cluster_num=20){
  
  CAM=CAM.nWCA(X=as.matrix(t(data_mix)), K=ncell,cluster_num=cluster_num)
  
  DSA.res= ged(as.matrix(data_est), MarkerList(CAM$V.pick), "DSA")
  P.DSA=DSA.res@fit@H
  
  P.CAM=scale(P.DSA, center=F, scale= apply(P.DSA, 2, sum))
  
  
  return(list("P"=P.CAM,"MKS"=CAM$V.pick))
  
}


ssCAM.main.est<-function(data_mix,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10){

  ssCAM=CAM.nWCA_force(data_mix, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_mix), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H

  return(list("P"=P.NMF,"MKS"=ssCAM$V.pick))
  
}



##ssCAM: run two NMF method together##
ssCAM.main.est.comb<-function(data_mix,ncell,cluster_num=20,mks_in,NMF.run=10){
  
  ssCAM=CAM.nWCA_force(data_mix, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.KL=ged(as.matrix(data_mix), MarkerList(ssCAM$V.pick), method="ssKL" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.KL=NMF.KL@fit@H
  
  NMF.Fro=ged(as.matrix(data_mix), MarkerList(ssCAM$V.pick), method="ssFrobenius" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.Fro=NMF.Fro@fit@H
  
  return(list("P.KL"=P.NMF.KL,"P.Fro"=P.NMF.Fro,"MKS"=ssCAM$V.pick))
  
}

ssCAM.main.est.comb.2data<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.run=10){
  #data_ssCAM=Data.sus.hc;data_est=Data.sus;ncell=k;cluster_num=50;mks_in=MKS.list.given;NMF.run=5
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.KL=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssKL" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.KL=NMF.KL@fit@H
  
  NMF.Fro=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssFrobenius" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.Fro=NMF.Fro@fit@H
  
  return(list("P.KL"=P.NMF.KL,"P.Fro"=P.NMF.Fro,"MKS"=ssCAM$V.pick))
  
}

##NMF: run two NMF method together##
NMF.main.est.comb<-function(data_mix,mks_in,NMF.run=10){
  

  NMF.KL=ged(as.matrix(data_mix), MarkerList(mks_in), method="ssKL" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.KL=NMF.KL@fit@H
  
  NMF.Fro=ged(as.matrix(data_mix), MarkerList(mks_in), method="ssFrobenius" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.Fro=NMF.Fro@fit@H
  
  return(list("P.KL"=P.NMF.KL,"P.Fro"=P.NMF.Fro))
  
}

NMF.main.est<-function(data_mix,mks_in,NMF.method,NMF.run=10){
  
  
  NMF.res=ged(as.matrix(data_mix), MarkerList(mks_in), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  return(list("P"=P.NMF))
  
}


#####################################
####only run once :no seed needed####
#####################################

semiCAM.main<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in){
  
  #data_ssCAM=data.mix;data_est=data.mix;ncell=2;cluster_num=50;mks_in=Marker.list;
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)

  NMF.KL=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssKL" ,log = FALSE, nrun = 1)
  P.NMF.KL=NMF.KL@fit@H
  X.NMF.KL=NMF.KL@fit@W
  
  
  
  return(list("P"=P.NMF.KL,"X"=X.NMF.KL,"MKS"=ssCAM$V.pick))
  
}



