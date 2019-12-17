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


##markers updating methods##
if(F){
ssCAM.main.update.est<-function(data_mix,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10){

  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  
  ssCAM=CAM.nWCA_force(data_mix, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_mix), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  NMF.update.res=umm_mks_NMF(data=data_mix,P.initial=P.NMF,MKS.known=mks_in,Markers.ls_initial=ssCAM$V.pick,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y")
  
  res=list("P"=P.NMF,"MKS.initial"=ssCAM$V.pick,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}
}

ssCAM.main.update.est<-function(data_mix,ncell,Ptrue,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=10,sh.cutoff=0.1){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  ssCAM=CAM.nWCA_force(data_mix, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_mix), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  P.map=match.cor(P.NMF,Ptrue)
  
  P.initial.order=P.map$P.order
  Markers.initial.order=ssCAM$V.pick[P.map$cell]
  
  cat("Start markers updating procedure")
  NMF.update.res=umm_mks_NMF(data=data_mix,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  cat("done!")
  COR.map.update=cor(as.vector(NMF.update.res$P),as.vector(Ptrue))
  
  res=list("P"=P.initial.order,"MKS.initial"=Markers.initial.order,"COR.initial"=P.map$COR,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS,"COR.update"=COR.map.update)
  
  return(res)
}

##data 1: identify markers; data 2: estimate proportions##

ssCAM.main.update.est.2data<-function(data_ssCAM,data_est,ncell,Ptrue,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=10,sh.cutoff=0.1){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;Ptrue=P.true.all;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10
  #data_ssCAM=data.mix;data_est=data.mix;ncell;Ptrue=P.true;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10;f.cutoff=5;sh.cutoff=0.1
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  P.map=match.cor(P.NMF,Ptrue)
  P.initial.order=P.map$P.order
  Markers.initial.order=ssCAM$V.pick[P.map$cell]
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.add.nodup(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff)
  
  cat("done! \n")
  
  P.map.update=match.cor(NMF.update.res$P,Ptrue)
  #COR.map.update=cor(as.vector(NMF.update.res$P),as.vector(Ptrue))
  
  res=list("P.initial"=P.initial.order,"MKS.initial"=Markers.initial.order,"COR.initial"=P.map$COR,"P.update"=P.map.update$P.order,"MKS.update"=NMF.update.res$MKS[P.map.update$cell],"COR.update"=P.map.update$COR)
  
  return(res)
}


##only update, no map##
ssCAM.main.update.est.2data.nomap<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=10,sh.cutoff=0.1){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;Ptrue=P.true.all;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10
  #data_ssCAM=data.mix;data_est=data.mix;ncell;Ptrue=P.true;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10;f.cutoff=5;sh.cutoff=0.1
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  NMF.res=ged(as.matrix(data_est), MarkerList(ssCAM$V.pick), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.add.nodup.simple(data_est,P.initial=P.NMF,Markers.ls_initial=ssCAM$V.pick,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff)
  
  cat("done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  
  res=list("P.initial"=P.NMF,"MKS.initial"=ssCAM$V.pick,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res) 
}


##update markers in and out, no map##

ssCAM.main.update.est.2data.inout<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=2^8){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;Ptrue=P.true.all;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssFrobenius";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=0.1
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3
  
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

ssCAM.main.update.est.2data.inout.new<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=1000){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;Ptrue=P.true.all;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssFrobenius";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=0.1
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=10000
  
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
  if(NMF.method=="ssKL"){X.NMF=NMF.res@fit@W}
  if(NMF.method=="ssFrobenius"){X.NMF=2^(NMF.res@fit@W)}
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.inout.new(data_est,P.initial=P.NMF,X.initial=X.NMF,MKS.known=mks_in,Markers.ls_initial=MKS.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  
  
  cat("done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  
  res=list("P.initial"=P.NMF,"MKS.initial"=MKS.order,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}

if(F){
ssCAM.main.update.est.2data.inout.new.MKS<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=1000){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;Ptrue=P.true.all;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssFrobenius";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=0.1
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=10000
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=10000
  
  
  ssCAM=CAM.nWCA_force(data_ssCAM, K=ncell,cluster_num=cluster_num,cluster_arg="Hartigan-Wong",mks_in)
  
  MKS.order=vector("list", length = ncell);cell.ind.rcd=NULL
  for (ele in 1:length(mks_in)){
    num.mks_in=unlist(lapply(ssCAM$V.pick,function(x) sum(x%in%mks_in[[ele]])))
    cell.ind.temp=which(num.mks_in>0)
    cell.ind.rcd=c(cell.ind.rcd,cell.ind.temp)
    MKS.order[[ele]]=ssCAM$V.pick[[cell.ind.temp]]
    
  }
  MKS.order[-c(1:length(mks_in))]=ssCAM$V.pick[-cell.ind.rcd]
  
  NMF.res=ged(as.matrix(data_est)[unlist(MKS.order),], MarkerList(MKS.order), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
  P.NMF=NMF.res@fit@H
  if(NMF.method=="ssKL"){X.NMF=NMF.res@fit@W}
  if(NMF.method=="ssFrobenius"){X.NMF=2^(NMF.res@fit@W)}
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.inout.new(data_est,P.initial=P.NMF,X.initial=X.NMF,MKS.known=mks_in,Markers.ls_initial=MKS.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  
  
  cat("done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  
  res=list("P.initial"=P.NMF,"MKS.initial"=MKS.order,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}
}
#########################################
##update markers by only adding, no map##
#########################################

ssCAM.main.update.est.2data.onlyadd<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=2^8){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;Ptrue=P.true.all;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssFrobenius";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=0.1
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3
  
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
  NMF.update.res=umm_mks_NMF.onlyadd(data_est,P.initial=P.NMF,MKS.known=mks_in,Markers.ls_initial=MKS.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  
  
  cat("done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  
  res=list("P.initial"=P.NMF,"MKS.initial"=MKS.order,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}


ssCAM.main.update.est.2data.onlyadd.new<-function(data_ssCAM,data_est,ncell,cluster_num=20,mks_in,NMF.method,NMF.run=10,maxInter=10,f.cutoff=3,sh.cutoff=2^8){
  
  #data_mix=data.mix.filter;ncell;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL"
  #data_mix=data.mix;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=10;maxInter=10
  
  #data_ssCAM=data.mix.filter.select;data_est=data.mix.filter;Ptrue=P.true.all;cluster_num=50;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=10
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssFrobenius";NMF.run=2;maxInter=5;f.cutoff=3;sh.cutoff=0.1
  #data_ssCAM=data.mix;data_est=data.mix;ncell;cluster_num=20;mks_in=MKS.initial;NMF.method="ssKL";NMF.run=2;maxInter=5;f.cutoff=3
  
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
  if(NMF.method=="ssKL"){X.NMF=NMF.res@fit@W}
  if(NMF.method=="ssFrobenius"){X.NMF=2^(NMF.res@fit@W)}
  
  
  cat("Start markers updating procedure \n")
  #NMF.update.res=umm_mks_NMF(data=data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=0.7,sh.cutoff=0.1)
  #NMF.update.res=umm_mks_NMF.mks(data_est,P.initial=P.initial.order,MKS.known=mks_in,Markers.ls_initial=Markers.initial.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  NMF.update.res=umm_mks_NMF.onlyadd.new(data_est,P.initial=P.NMF,X.initial=X.NMF,MKS.known=mks_in,Markers.ls_initial=MKS.order,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff)
  
  
  cat("done! \n")
  
  #.map.update=match.cor(NMF.update.res$P,Ptrue)
  
  
  res=list("P.initial"=P.NMF,"MKS.initial"=MKS.order,"P.update"=NMF.update.res$P,"MKS.update"=NMF.update.res$MKS)
  
  return(res)
}