CAM.main.est<-function(data_mix,ncell,cluster_num=20){

  CAM=CAM.nWCA(X=as.matrix(t(data_mix)), K=ncell,cluster_num=cluster_num)
  P.CAM=t(CAM$A_est)
  P.CAM=scale(P.CAM, center=F, scale= apply(P.CAM, 2, sum))
  
  
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

##NMF: run two NMF method together##
NMF.main.est.comb<-function(data_mix,mks_in,NMF.run=10){
  

  NMF.KL=ged(as.matrix(data_mix), MarkerList(mks_in), method="ssKL" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.KL=NMF.KL@fit@H
  
  NMF.Fro=ged(as.matrix(data_mix), MarkerList(mks_in), method="ssFrobenius" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.Fro=NMF.Fro@fit@H
  
  return(list("P.KL"=P.NMF.KL,"P.Fro"=P.NMF.Fro))
  
}

ssCAM.main.update.est<-function(data_mix,ncell,Ptrue,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10){

  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
 #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  ssCAM=CAM.nWCA_force(data_mix, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_mix), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  P.map=match.cor(P.NMF,Ptrue)
  
  P.initial.order=P.map$P.order
  Markers.initial.order=ssCAM$V.pick[P.map$cell]
  
  cat("Start markers updating procedure \n")
  NMF.update.res=umm_mks_NMF(data=data_mix,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  cat("done! \n")
  COR.map.update=cor(as.vector(NMF.update.res$P),as.vector(Ptrue))
  
  res=list("P"=P.initial.order,"MKS.initial"=Markers.initial.order,"COR.initial"=P.map$COR,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS,"COR.update"=COR.map.update)
  
  return(res)
}

if(F){
ssCAM.main.update.est.2data<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10){
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  if(F){
  NMF.KL=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssKL" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.KL=NMF.KL@fit@H
  
  NMF.Fro=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssFrobenius" ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF.Fro=NMF.Fro@fit@H
  }
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$MKS), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.NMF,MKS.known=mks_in,Markers.ls_initial=ssCAM$MKS,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y")
  
  res=list("P.initial"=P.NMF,"MKS.initial"=ssCAM$V.pick,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}
}

ssCAM.main.update.est.2data<-function(data_ssCAM,data_est,Ptrue,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;Ptrue=P.true;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  P.map=match.cor(P.NMF,Ptrue)
  
  P.initial.order=P.map$P.order
  Markers.initial.order=ssCAM$V.pick[P.map$cell]
  
  cat("Start markers updating procedure \n")
  NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  cat("done! \n")
  COR.map.update=cor(as.vector(NMF.update.res$P),as.vector(Ptrue))
  
  res=list("P"=P.initial.order,"MKS.initial"=Markers.initial.order,"COR.initial"=P.map$COR,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS,"COR.update"=COR.map.update)
  
  return(res)
}