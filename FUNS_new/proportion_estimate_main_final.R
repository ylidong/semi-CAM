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

ssCAM.main.est.2data<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10){
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
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


##markers updating methods: only add markers##

ssCAM.main.update.final<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=1000){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=2^9;
  
  #data_ssCAM=Data.sus.hc;data_est=Data.sus;ncell=k;cluster_num=50;mks_in=MKS.list.given;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  #data_ssCAM=data.mix.filter;data_est=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)

  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  if(NMF.method=="ssKL"){X.NMF=NMF.res@fit@W}
  if(NMF.method=="ssFrobenius"){X.NMF=2^(NMF.res@fit@W)}
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.final(data_est,P.initial=P.NMF,X.initial=X.NMF,Markers.ls_initial=ssCAM$V.pick,maxInter,NMF.method,NMF.run,f.cutoff,sh.cutoff)
  
  cat("Done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  res=list("P.initial"=P.NMF,"MKS.initial"=ssCAM$V.pick,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}

ssCAM.main.update.final.last<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=1000){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=2^9;
  
  #data_ssCAM=Data.sus.hc;data_est=Data.sus;ncell=k;cluster_num=50;mks_in=MKS.list.given;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  #data_ssCAM=data.mix.filter;data_est=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  if(NMF.method=="ssKL"){X.NMF=NMF.res@fit@W}
  if(NMF.method=="ssFrobenius"){X.NMF=2^(NMF.res@fit@W)}
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.final.last(data_est,P.initial=P.NMF,X.initial=X.NMF,Markers.ls_initial=ssCAM$V.pick,maxInter,NMF.method,NMF.run,f.cutoff,sh.cutoff)
  
  cat("Done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  res=list("P.initial"=P.NMF,"MKS.initial"=ssCAM$V.pick,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}


ssCAM.main.update.final.last2<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=1000){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=2^9;
  
  #data_ssCAM=Data.sus.hc;data_est=Data.sus;ncell=k;cluster_num=50;mks_in=MKS.list.given;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  #data_ssCAM=data.mix.filter;data_est=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  if(NMF.method=="ssKL"){X.NMF=NMF.res@fit@W}
  if(NMF.method=="ssFrobenius"){X.NMF=2^(NMF.res@fit@W)}
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.final.last2(data_est,P.initial=P.NMF,X.initial=X.NMF,Markers.ls_initial=ssCAM$V.pick,maxInter,NMF.method,NMF.run,f.cutoff,sh.cutoff)
  
  cat("Done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  res=list("P.initial"=P.NMF,"MKS.initial"=ssCAM$V.pick,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}
##markers updating methods: only use marker genes for NMF##

ssCAM.main.update.final.MKS<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=1000){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=2^8;
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  MKS.order=vector("list", length = ncell);cell.ind.rcd=NULL
  for (ele in 1:length(mks_in)){
    num.mks_in=unlist(lapply(ssCAM$V.pick,function(x) sum(x%in%mks_in[[ele]])))
    cell.ind.temp=which(num.mks_in>0)
    cell.ind.rcd=c(cell.ind.rcd,cell.ind.temp)
    MKS.order[[ele]]=ssCAM$V.pick[[cell.ind.temp]]
    
  }
  MKS.order[-c(1:length(mks_in))]=ssCAM$V.pick[-cell.ind.rcd]
  
  NMF.res=ged(as.matrix(data_est[unlist(MKS.order),]), MarkerList(MKS.order), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  #if(NMF.method=="ssKL"){X.NMF=NMF.res@fit@W}
  #if(NMF.method=="ssFrobenius"){X.NMF=2^(NMF.res@fit@W)}
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.final.MKS(data_est,P.initial=P.NMF,Markers.known=mks_in,Markers.ls_initial=MKS.order,maxInter,NMF.method,NMF.run,f.cutoff,sh.cutoff)
  
  cat("Done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  res=list("P.initial"=P.NMF,"MKS.initial"=ssCAM$V.pick,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}


##markers updating methods: markers can be in and out CSE is re-estimated##

ssCAM.main.update.inout<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=2^8){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  MKS.order=vector("list", length = ncell);cell.ind.rcd=NULL
  for (ele in 1:length(mks_in)){
    num.mks_in=unlist(lapply(ssCAM$V.pick,function(x) sum(x%in%mks_in[[ele]])))
    cell.ind.temp=which(num.mks_in>0)
    cell.ind.rcd=c(cell.ind.rcd,cell.ind.temp)
    MKS.order[[ele]]=ssCAM$V.pick[[cell.ind.temp]]
    
  }
  MKS.order[-c(1:length(mks_in))]=ssCAM$V.pick[-cell.ind.rcd]
  
  NMF.res=ged(as.matrix(data_est), MarkerList(MKS.order), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.inout(data_est,P.initial=P.NMF,MKS.known=mks_in,Markers.ls_initial=MKS.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  
  
  cat("done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  
  res=list("P.initial"=P.NMF,"MKS.initial"=MKS.order,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}

##updating methods: run more NMF##

ssCAM.main.update.final.NMF<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=1000){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=2^9;
  
  #data_ssCAM=Data.sus.hc;data_est=Data.sus;ncell=k;cluster_num=50;mks_in=MKS.list.given;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  #data_ssCAM=data.mix.filter;data_est=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=2;sh.cutoff=2^9
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = 2)
  P.NMF=NMF.res@fit@H
  
  cat("Start markers updating procedure \n")
  
  NMF.res.update=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = 10)
  P.NMF.update=NMF.res.update@fit@H
  
  cat("Done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  res=list("P.initial"=P.NMF,"MKS.initial"=ssCAM$V.pick,"P.update"=P.NMF.update,"MKS.update"=ssCAM$V.pick)
  
  return(res)
}

##updating methods: run more NMF##

ssCAM.main.update.final.NMF.nrun<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2
  
 
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE, nrun = 1,seed="rprop")
  P.NMF.initial=NMF.res@fit@H
  print(summary(NMF.res)[6])
  print(NMF.res@fit@H)
  showRNG(NMF.res)
  
  for (r in 1:NMF.run){
    set.seed(r)
    tmp=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE, nrun = 1,seed="rprop")
    print(summary(tmp)[6])
    showRNG(tmp)
    print(tmp@fit@H)
    if (summary(NMF.res)[6]>summary(tmp)[6]){
    NMF.res <- tmp
    }
  }
  P.NMF.update=NMF.res@fit@H
  
  
  cat("Done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  res=list("P.initial"=P.NMF.initial,"MKS.initial"=ssCAM$V.pick,"P.update"=P.NMF.update,"MKS.update"=ssCAM$V.pick)
  
  return(res)
}


#####################################
####only run once :no seed needed####
#####################################

ssCAM.main.NMF.1run<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in){
  
  #data_ssCAM=data.mix.use;data_est=data.mix.use;ncell;cluster_num=50;mks_in=MKS.initial;
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)

  NMF.KL=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssKL" ,log = FALSE, nrun = 1)
  P.NMF.KL=NMF.KL@fit@H
  X.NMF.KL=NMF.KL@fit@W
  
  
  NMF.Fro=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method="ssFrobenius" ,log = FALSE, nrun = 1)
  P.NMF.Fro=NMF.Fro@fit@H
  X.NMF.Fro=NMF.Fro@fit@W
  
  
  return(list("P.KL"=P.NMF.KL,"X.KL"=X.NMF.KL,"P.Fro"=P.NMF.Fro,"X.Fro"=X.NMF.Fro,"MKS"=ssCAM$V.pick))
  
}



main.NMF.1run<-function(data_est,mks_in){
  
 
  NMF.KL=ged(as.matrix(data_est), MarkerList(mks_in), method="ssKL" ,log = FALSE, nrun = 1)
  P.NMF.KL=NMF.KL@fit@H
  X.NMF.KL=NMF.KL@fit@W
  
  NMF.Fro=ged(as.matrix(data_est), MarkerList(mks_in), method="ssFrobenius" ,log = FALSE, nrun = 1)
  P.NMF.Fro=NMF.Fro@fit@H
  X.NMF.Fro=NMF.Fro@fit@W
  
  
  return(list("P.KL"=P.NMF.KL,"X.KL"=X.NMF.KL,"P.Fro"=P.NMF.Fro,"X.Fro"=X.NMF.Fro))
  
}
