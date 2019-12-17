Est.cse<-function(prop,Amat){
  
  #prop=P.initial;Amat=data #col(prop)=mixture row(prop)=cell types
  
  S_est=apply(Amat,1,function(x) coef(nnls(t(prop),x))) ##after transformation, rows are samples
  
  S_est
}


GET.newmarkers<-function(x,f.cutoff,sh.cutoff){
  
  #x=X.initial;f.cutoff=2;sh.cutoff=2^9
  
  K=ncol(x)
  
  ##markers only express on a single cell##
  MKS_0=vector("list", length = K)
  for(k in 1:K){
    mks.ind=apply(x, 1, function(m) (m[k]!=0) & (sum(m[-k]==0)==(K-1)))
    MKS_0[[k]]=names(mks.ind)[mks.ind==T]
  }
  
  ##markers highly express on a single cell##
  x.which.max.ind=apply(x[!rownames(x)%in%unlist(MKS_0),],1,which.max)
  x.which.max.list=lapply(1:K,function(k) names(x.which.max.ind)[which((x.which.max.ind==k)==T)])
  
  MKS_1=vector("list", length = K)
  for(i in 1:K){
    mks1.ind=apply(x[x.which.max.list[[i]],], 1, function(m){ 
      h.exp=m[i]
      sh.exp=sort(m)[K-1]
      f=h.exp/sh.exp
      
      if((f>f.cutoff)& (sh.exp<sh.cutoff)){ind=T
      }else{ind=F}
      })
    MKS_1[[i]]=names(mks1.ind)[mks1.ind==T]
  } 
    
  MKS_new=mapply(c,MKS_0,MKS_1, SIMPLIFY=FALSE)
  
  return(MKS_new)
}

GET.newmarkers2<-function(x,f.cutoff,sh.cutoff){
  
  #x=X.initial;f.cutoff=2;sh.cutoff=2^9
  
  K=ncol(x)
  
  ##markers only express on a single cell##
  MKS_0=vector("list", length = K)
  for(k in 1:K){
    mks.ind=apply(x, 1, function(m) (m[k]!=0) & (sum(m[-k]==0)==(K-1)))
    MKS_0[[k]]=names(mks.ind)[mks.ind==T]
  }
  
  ##markers highly express on a single cell##
  x.which.max.ind=apply(x[!rownames(x)%in%unlist(MKS_0),],1,which.max)
  x.which.max.list=lapply(1:K,function(k) names(x.which.max.ind)[which((x.which.max.ind==k)==T)])
  
  MKS_1=vector("list", length = K)
  for(i in 1:K){
    mks1.ind=apply(x[x.which.max.list[[i]],], 1, function(m){ 
      h.exp=m[i]
      sh.exp=sort(m)[K-1]
      f=h.exp/sh.exp
      
      if((f>f.cutoff)& (sh.exp<sh.cutoff) & (h.exp>2^9)){ind=T
      }else{ind=F}
    })
    MKS_1[[i]]=names(mks1.ind)[mks1.ind==T]
  } 
  
  MKS_new=mapply(c,MKS_0,MKS_1, SIMPLIFY=FALSE)
  
  return(MKS_new)
} 
#################################################################################
##update markers based on the initial give markers: markers can in      ##
#################################################################################

umm_mks_NMF.final<-function(data,P.initial,X.initial,Markers.ls_initial,maxInter,NMF.method,NMF.run,f.cutoff,sh.cutoff){
  
  #data=data_est;P.initial=P.NMF;X.initial=X.NMF;Markers.ls_initial=ssCAM$V.pick;maxInter;NMF.method;NMF.run;criterion="MSE.Y";f.cutoff=2;sh.cutoff=2^9

  
  MSE.Y.initial=sum((X.initial%*%P.initial-data)^2)/(ncol(data)*nrow(data))
  
  
  inter=1;res.list=list()
  
  res.list[[inter]]=list("MKS"=Markers.ls_initial,"P"=P.initial,"MSE.Y"=MSE.Y.initial)
  
  
  while(inter<maxInter){
    
    
    inter=inter+1
    #print(inter)
    
    Markers.ls_update=GET.newmarkers(x=X.initial,f.cutoff,sh.cutoff)
    #print(Markers.ls_update[[1]])
    
    NMF.res<- ged(as.matrix(data), MarkerList(Markers.ls_update), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
    P.update=NMF.res@fit@H
    
    if(NMF.method=="ssKL"){X.update=NMF.res@fit@W}
    if(NMF.method=="ssFrobenius"){X.update=2^(NMF.res@fit@W)}
    
    MSE.Y.update=sum((X.update%*%P.update-data)^2)/(ncol(data)*nrow(data))
    X.initial=X.update
    
    res.list[[inter]]=list("MKS"=Markers.ls_update,"P"=P.update,"MSE.Y"=MSE.Y.update)
    
  }
  
  optim=which.min(unlist(lapply(res.list,function(x) x$MSE.Y)))
  
  res.out=res.list[[optim]]
  res.out[["P.initial"]]=P.initial
  
  return(res.out) 
  
}

umm_mks_NMF.final.last<-function(data,P.initial,X.initial,Markers.ls_initial,maxInter,NMF.method,NMF.run,f.cutoff,sh.cutoff){
  
  #data=data_est;P.initial=P.NMF;X.initial=X.NMF;Markers.ls_initial=ssCAM$V.pick;maxInter;NMF.method;NMF.run;criterion="MSE.Y";f.cutoff=2;sh.cutoff=2^9
  
  
  MSE.Y.initial=sum((X.initial%*%P.initial-data)^2)/(ncol(data)*nrow(data))
  
  
  inter=1;res.list=list()
  
  res.list[[inter]]=list("MKS"=Markers.ls_initial,"P"=P.initial,"MSE.Y"=MSE.Y.initial)
  
  
  while(inter<maxInter){
    
    
    inter=inter+1
    #print(inter)
    
    Markers.ls_update=GET.newmarkers(x=X.initial,f.cutoff,sh.cutoff)
    #print(Markers.ls_update[[1]])
    
    NMF.res<- ged(as.matrix(data), MarkerList(Markers.ls_update), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
    P.update=NMF.res@fit@H
    
    if(NMF.method=="ssKL"){X.update=NMF.res@fit@W}
    if(NMF.method=="ssFrobenius"){X.update=2^(NMF.res@fit@W)}
    
    MSE.Y.update=sum((X.update%*%P.update-data)^2)/(ncol(data)*nrow(data))
    X.initial=X.update
    
    res.list[[inter]]=list("MKS"=Markers.ls_update,"P"=P.update,"MSE.Y"=MSE.Y.update)
    
  }
  
  #optim=which.min(unlist(lapply(res.list,function(x) x$MSE.Y)))
  
  res.out=res.list[[maxInter]]
  res.out[["P.initial"]]=P.initial
  
  return(res.out) 
  
}

umm_mks_NMF.final.last2<-function(data,P.initial,X.initial,Markers.ls_initial,maxInter,NMF.method,NMF.run,f.cutoff,sh.cutoff){
  
  #data=data_est;P.initial=P.NMF;X.initial=X.NMF;Markers.ls_initial=ssCAM$V.pick;maxInter;NMF.method;NMF.run;criterion="MSE.Y";f.cutoff=2;sh.cutoff=2^9
  
  
  MSE.Y.initial=sum((X.initial%*%P.initial-data)^2)/(ncol(data)*nrow(data))
  
  
  inter=1;res.list=list()
  
  res.list[[inter]]=list("MKS"=Markers.ls_initial,"P"=P.initial,"MSE.Y"=MSE.Y.initial)
  
  
  while(inter<maxInter){
    
    
    inter=inter+1
    #print(inter)
    
    Markers.ls_update=GET.newmarkers2(x=X.initial,f.cutoff,sh.cutoff)
    #print(Markers.ls_update[[1]])
    
    NMF.res<- ged(as.matrix(data), MarkerList(Markers.ls_update), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
    P.update=NMF.res@fit@H
    
    if(NMF.method=="ssKL"){X.update=NMF.res@fit@W}
    if(NMF.method=="ssFrobenius"){X.update=2^(NMF.res@fit@W)}
    
    MSE.Y.update=sum((X.update%*%P.update-data)^2)/(ncol(data)*nrow(data))
    X.initial=X.update
    
    res.list[[inter]]=list("MKS"=Markers.ls_update,"P"=P.update,"MSE.Y"=MSE.Y.update)
    
  }
  
  #optim=which.min(unlist(lapply(res.list,function(x) x$MSE.Y)))
  
  res.out=res.list[[maxInter]]
  res.out[["P.initial"]]=P.initial
  
  return(res.out) 
  
}
#################################################################################
##update markers based on the initial give markers: only use markers to estimate proportion    ##
#################################################################################

umm_mks_NMF.final.MKS<-function(data,P.initial,Markers.known,Markers.ls_initial,maxInter,NMF.method,NMF.run,f.cutoff,sh.cutoff){
  
  #data=data_est;P.initial=P.NMF;Markers.known=mks_in;Markers.ls_initial=MKS.order;maxInter;NMF.method;NMF.run;criterion="MSE.Y";f.cutoff=2;sh.cutoff=2^8
  
  X.initial=t(Est.cse(prop=P.initial,Amat=data))
  MSE.Y.initial=sum((X.initial%*%P.initial-data)^2)/(ncol(data)*nrow(data))
  
  
  inter=1;res.list=list()
  
  res.list[[inter]]=list("MKS"=Markers.ls_initial,"P"=P.initial,"MSE.Y"=MSE.Y.initial)
  
  
  while(inter<maxInter){
    
    
    inter=inter+1
    #print(inter)
    
    #find new markers#
    Markers.ls_update=GET.newmarkers(x=X.initial[!rownames(X.initial)%in%unlist(Markers.known),],f.cutoff,sh.cutoff)
    #add Markers.known
    for(ii in 1:length(Markers.known)){ ##original known markers + new finding markers##
      Markers.ls_update[[ii]]=union(Markers.ls_update[[ii]],Markers.known[[ii]])
    }
    NMF.res.update<- ged(as.matrix(data[unlist(Markers.ls_update),]), MarkerList(Markers.ls_update), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
    P.update=NMF.res.update@fit@H
    
    #if(NMF.method=="ssKL"){X.update=NMF.res@fit@W}
    #if(NMF.method=="ssFrobenius"){X.update=2^(NMF.res@fit@W)}
    X.update=t(Est.cse(prop=P.update,Amat=data))
    
    MSE.Y.update=sum((X.update%*%P.update-data)^2)/(ncol(data)*nrow(data))
    X.initial=X.update
    
    res.list[[inter]]=list("MKS"=Markers.ls_update,"P"=P.update,"MSE.Y"=MSE.Y.update)
    
  }
  
  optim=which.min(unlist(lapply(res.list,function(x) x$MSE.Y)))
  
  res.out=res.list[[optim]]
  res.out[["P.initial"]]=P.initial
  
  return(res.out) 
  
}


#################################################################################
##update markers based on the initial give markers: markers can in and out?? Cell type specific expression is re-estimated      ##
#################################################################################

umm_mks_NMF.inout<-function(data,P.initial,MKS.known,Markers.ls_initial,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff,sh.cutoff){
  
  #data=data_est;P.initial=P.NMF;MKS.known=mks_in;Markers.ls_initial=MKS.order;maxInter;ncell;NMF.method;NMF.run;criterion="MSE.Y";f.cutoff;sh.cutoff
  
  #1. estimate cell type specific expression
  A_cse=Est.cse(prop=P.initial,Amat=data)
  
  MSE.Y.initial=sum((t(A_cse)%*%P.initial-data)^2)/(ncol(data)*nrow(data))
  
  
  inter=1;res.list=list()
  
  res.list[[inter]]=list("MKS"=Markers.ls_initial,"P"=P.initial,"MSE.Y"=MSE.Y.initial)
  
  
  while(inter<maxInter){
    
    
    inter=inter+1
    #print(inter)
    
    Markers.ls_update=GET.newmarkers(x=t(A_cse[,!colnames(A_cse)%in%unlist(MKS.known)]),f.cutoff,sh.cutoff)
    
    #add initial known markers#
    for(ele in 1:length(MKS.known)){
      Markers.ls_update[[ele]]=union(Markers.ls_update[[ele]],MKS.known[[ele]])
    }
    
    ##check if there is zero markers could be added to the unkown markers cell type
    cell.zero=which(sapply(Markers.ls_update,length)==0)
    if(length(cell.zero)!=0){
      for(i in cell.zero){
        Markers.ls_update[[i]]=Markers.ls_initial[[i]]}
    }
    #print(Markers.ls_update)
    NMF.res<- ged(as.matrix(data), MarkerList(Markers.ls_update), method=NMF.method ,log = FALSE,  rng = 1234, nrun = NMF.run)
    P.update=NMF.res@fit@H
    
    A_cse=Est.cse(prop=P.update,Amat=data)
    
    MSE.Y.update=sum((t(A_cse)%*%P.update-data)^2)/(ncol(data)*nrow(data))
    
    res.list[[inter]]=list("MKS"=Markers.ls_update,"P"=P.update,"MSE.Y"=MSE.Y.update)
    
  }
  
  if(criterion=="MSE.Y") optim=which.min(unlist(lapply(res.list,function(x) x$MSE.Y)))
  
  res.out=res.list[[optim]]
  res.out[["P.initial"]]=P.initial
  
  return(res.out) 
  
}
