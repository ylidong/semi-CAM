Get.12mks_new3<-function(mat_in,num.cell,f.cutoff,sh.cutoff){
  #mat_in=est.cse_mks2;fold_cutoff=10;num.cell=4
  mat.rank=apply(mat_in,2,function(x) rank(x)) ###get rank of each cell type
  
  ##get highest expression values and the high expresssion cell type indicator##
  high=apply(mat_in,2,max) ##highest expression
  high.ind=apply(mat.rank,2,function(xx) which(xx==num.cell)) ##which cell type has the highest expression 
  
  ##get the second highest expression values##
  #sec.high.ind=apply(mat.rank,2,function(xx) which(xx==(num.cell-1)))
  sec.high=apply(mat_in,2,function(x) x[rank(x)==(num.cell-1)])
  fold=high/sec.high ##fold high/sec.high
  
  #second highest expression is ignorable#
  
  
  ##markers for each cell type##
  mk.ind=high.ind[(fold> f.cutoff)&(sec.high<quantile(sec.high,sh.cutoff))]
  print(summary(fold))
  mk.add=lapply(1:dim(mat_in)[1],function(i) names(mk.ind)[which(mk.ind==i)]) 
  mk.add 
}



GET.newmarkers_new3<-function(est.cse_in,num.cell,f.cutoff,sh.cutoff){
  
  #est.cse_in=A_cse;num.cell=4
  
  ###markers only expressed in one cell type###
  
  mks=apply(est.cse_in, 2, function(c) sum(c!=0))
  
  est.cse_mks=est.cse_in[,which(mks==1),drop=F]
  #colnames(est.cse_mks)=colnames(est.cse_in)[which(mks==1)]
  
  est.cse_mks_ind=apply(est.cse_mks, 2, function(c) which(c!=0))
  
  est.cse_mks_add=lapply(1:dim(est.cse_in)[1],function(i) names(est.cse_mks_ind)[which(est.cse_mks_ind==i)]) 
  
  ###markers expressed in >=2 cell types###
  est.cse_mks2=est.cse_in[,which(mks>=2),drop=F]
  #colnames(est.cse_mks2)=colnames(est.cse_in)[which(mks>=2)]
  est.cse_mks2_add=Get.12mks_new3(mat_in=est.cse_mks2,num.cell,f.cutoff,sh.cutoff) ##add markers detected by fold change high/sec.high
  
  #lapply(1:4,function(i) length(est.cse_mks2_add[[i]]))
  
  markers_update=mapply(c,est.cse_mks_add,est.cse_mks2_add, SIMPLIFY=FALSE)
  
  markers_update
  
  #lapply(1:4,function(i) length(markers_update[[i]]))
}


Est.cse<-function(prop,Amat){
  
  #prop=P.initial;Amat=data #col(prop)=mixture row(prop)=cell types
  
  S_est=apply(Amat,1,function(x) coef(nnls(t(prop),x))) ##after transformation, rows are samples
  
  S_est
}


umm_mks_NMF<-function(data,P.initial,MKS.known,Markers.ls_initial,maxInter,ncell,NMF.method,NMF.run,criterion="MSE.Y",f.cutoff=10,sh.cutoff=0.1){
  
  #data=data.mix.filter;P.initial=prop.est;MKS.known=NULL;Markers.ls_initial= CAM.org.out$V.pick;maxInter=10;ncell=ncell;NMF.method="ssKL";criterion="MSE.Y"
  #data=data.mix.filter;P.initial=ssCAM$P.KL;MKS.known=MKS.initial;Markers.ls_initial=ssCAM$MKS;maxInter=10;NMF.method="ssKL" 
  #data=data_est;P.initial=P.NMF;MKS.known=mks_in;Markers.ls_initial=ssCAM$MKS;criterion="MSE.Y"
  
  #data=data_mix;P.initial=P.map$P.order;MKS.known=mks_in;Markers.ls_initial=ssCAM$V.pick[P.map$ncell];maxInter=10;NMF.run=10;criterion="MSE.Y"
  
  if(length(MKS.known)==0) Markers.ls_knwon=lapply(1:ncell,function(x) NULL)
  
  if( (length(MKS.known)>0) & (length(MKS.known)<ncell)){ 
    Markers.ls_knwon=lapply(1:ncell,function(x) NULL)
    for(x in 1:length(Markers.ls_knwon)){
      if(x<=length(MKS.known)){Markers.ls_knwon[[x]]=MKS.known[[x]]
      }
      
    }
  }
  
  if(length(MKS.known)==ncell) Markers.ls_knwon=MKS.known
  
  
  #P.initial=t(res.out$A_est)
  A_cse=Est.cse(prop=P.initial,Amat=data)
  
  MSE.Y.initial=sum((t(A_cse)%*%P.initial-data)^2)/(ncol(data)*nrow(data))
  
  
  inter=1;res.list=list()
  
  res.list[[inter]]=list("MKS"=Markers.ls_initial,"P"=P.initial,"MSE.Y"=MSE.Y.initial)
  
  
  while(inter<maxInter){
    
    
    inter=inter+1
    
    
    Markers.ls_update=GET.newmarkers_new3(est.cse_in=A_cse,num.cell=ncell,f.cutoff,sh.cutoff)
    
    Markers.ls_update=mapply(union,Markers.ls_update,Markers.ls_knwon,SIMPLIFY=FALSE) ##original known markers + new finding markers##
    
    ##check if there is zero markers could be added to the unkown markers cell type
    cell.zero=which(sapply(Markers.ls_update,length)==0)
    if(length(cell.zero)!=0){
      for(i in cell.zero){
        Markers.ls_update[[i]]=Markers.ls_initial[[i]]}
    }
    
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
