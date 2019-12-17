match.sse<-function(P.est,P.true){
  #P.est=prop.est;P.true=P.true
  #P.est=P
  packageExist <- require("combinat")
  if (!packageExist) {
    install.packages("combinat")
    library("combinat")
  }
  
  pos.ord=permn(1:nrow(P.true))
  sse=c(0,length(pos.ord))
  for(i in 1:length(pos.ord)){
    P.try=P.est
    sse[i]=sum((P.try[pos.ord[[i]],]-P.true)^2)
    
  }
  
  tt=pos.ord[[which.min(sse)]]
  COR=cor(as.vector(P.est[tt,]),as.vector(P.true))
  
  P.order=P.est[tt,]
  list(err=min(sse), cell=tt,P.order=P.order,COR=COR)
}


match.mse<-function(P.est,P.true){
  #P.est=prop.est;P.true=P.true
  #P.est=P
  packageExist <- require("combinat")
  if (!packageExist) {
    install.packages("combinat")
    library("combinat")
  }
  
  pos.ord=permn(1:nrow(P.true))
  sse=c(0,length(pos.ord))
  for(i in 1:length(pos.ord)){
    P.try=P.est
    sse[i]=sum((P.try[pos.ord[[i]],]-P.true)^2)/length(P.true)
    
  }
  
  tt=pos.ord[[which.min(sse)]]
  COR=cor(as.vector(P.est[tt,]),as.vector(P.true))
  
  P.order=P.est[tt,]
  list(mse=min(sse), cell=tt,P.order=P.order,COR=COR)
}

match.cor<-function(P.est,P.true){
  #P.est=prop.est1;P.true=P.true
  #P.est=P.rela.est;P.true=P.rela_in
  packageExist <- require("combinat")
  if (!packageExist) {
    install.packages("combinat")
    library("combinat")
  }
  
  pos.ord=permn(1:nrow(P.true))
  cor=c(0,length(pos.ord))
  for(i in 1:length(pos.ord)){
    P.try=P.est
    cor[i]=cor(as.vector(P.try[pos.ord[[i]],]),as.vector(P.true))
  }
  
  tt=pos.ord[[which.max(cor)]]
  COR=cor(as.vector(P.est[tt,]),as.vector(P.true))
  
  P.order=P.est[tt,]
  rownames(P.order)=rownames(P.true)
  sse=sum((P.order-P.true)^2)
  
  list(err=sse, cell=tt,P.order=P.order,COR=COR)
}

match.cor.mse<-function(P.est,P.true){
  #P.est=prop.est1;P.true=P.true
  #P.est=P.rela.est;P.true=P.rela_in
  packageExist <- require("combinat")
  if (!packageExist) {
    install.packages("combinat")
    library("combinat")
  }
  
  pos.ord=permn(1:nrow(P.true))
  cor=c(0,length(pos.ord))
  for(i in 1:length(pos.ord)){
    P.try=P.est
    cor[i]=cor(as.vector(P.try[pos.ord[[i]],]),as.vector(P.true))
  }
  
  tt=pos.ord[[which.max(cor)]]
  COR=cor(as.vector(P.est[tt,]),as.vector(P.true))
  
  P.order=P.est[tt,]
  rownames(P.order)=rownames(P.true)
  mse=sum((P.order-P.true)^2)/length(P.true)
  
  list(mse=mse, cell=tt,P.order=P.order,COR=COR)
}


match.cor.MKS<-function(P.est,P.true,mks_known,mks_find){
  #P.est=ssCAM.res$P.update;P.true=P.true;mks_known=MKS.list.order;mks_find=ssCAM.res$MKS.update
  #P.est=ssCAM.res$P;P.true=P.true;mks_known=MKS.list.order;mks_find=ssCAM.res$MKS.initial
  
  #P.est=CAM.res$P;P.true=P.true;mks_known=MKS.list.order;mks_find=CAM.res$MKS
  #P.est=ssCAM.res$P;P.true;mks_known=MKS.list.order;mks_find=ssCAM.res$MKS.initial
  
  rownames(P.est)=paste0(1:nrow(P.est))
  if(length(mks_known)==0){stop('no known markers, cannot map to cell types!')}
  
  if(length(mks_known)<nrow(P.est)){
  ##map the prop to cell type having known markers (using markers)##
  map.flag=0
  for(m in 1:length(mks_known) ){
    
    celltype.ind=which(unlist(lapply(mks_find, function(x) sum(mks_known[[m]]%in%x)==length(mks_known[[m]])) ))
    #print(celltype.ind)
    ##one to one##
    if(length(celltype.ind)==1){rownames(P.est)[celltype.ind]=names(mks_known)[m];map.flag=1}
    
    ##one to more## totally wrong
    if(length(celltype.ind)>1) {map.flag=2}
  }
  
  #using markers map some cell types
  if(map.flag==1){
    P.est.map=P.est[rownames(P.est)%in%names(mks_known),,drop=F]
    
    ##map the prop to cell type not known having markers (using P.true)##
    
    if(length(P.true)!=0){
      if(nrow(P.est.map)==(nrow(P.true)-1)){
        map.cell=rownames(P.est)[rownames(P.est)%in%names(mks_known)]
        rownames(P.est)[rownames(P.est)==setdiff(rownames(P.est),map.cell)]=setdiff(rownames(P.true),map.cell)
        
        P.est.back=P.est[rownames(P.true),]
        cor.aftermap=cor(as.vector(P.est.back),as.vector(P.true))
      }else{
        map.cell=rownames(P.est)[rownames(P.est)%in%names(mks_known)]
        
        P.est.nomap=P.est[!rownames(P.est)%in%map.cell,]
        P.true.nomap=P.true[!rownames(P.true)%in%map.cell,]
        
        temp.map.cor=match.cor(P.est=P.est.nomap,P.true=P.true.nomap)
        rownames(P.est.nomap)[temp.map.cor$cell]=rownames(P.true.nomap)
        
        P.est.back=rbind(P.est.map,P.est.nomap)
        P.est.back=P.est.back[rownames(P.true),]
        
        cor.aftermap=cor(as.vector(P.est.back),as.vector(P.true))
      }
    }else{
      P.est.back=P.est
      cor.aftermap=NULL
    }
  }
  
  #using markers map none cell types
  if(map.flag==0){
    P.est.map=NULL
    if(length(P.true)!=0){
      P.est.nomap=P.est
      P.true.nomap=P.true
      
      temp.map.cor=match.cor(P.est=P.est.nomap,P.true=P.true.nomap)
      rownames(P.est.nomap)[temp.map.cor$cell]=rownames(P.true.nomap)
      
      P.est.back=rbind(P.est.map,P.est.nomap)
      P.est.back=P.est.back[rownames(P.true),]
      
      cor.aftermap=cor(as.vector(P.est.back),as.vector(P.true))
    }else{
      P.est.back=P.est
      cor.aftermap=NULL
    }
  }
  
  #using markers totally wrong
  if(map.flag==2){
    P.est.back=P.est
    cor.aftermap=0
  }
  }
  if(length(mks_known)==nrow(P.est)){
  if(length(P.true)!=0){
    P.est.back=P.est
    cor.aftermap=cor(as.vector(P.est.back),as.vector(P.true))
  }else{
    P.est.back=P.est
    cor.aftermap=NULL
  }
  }
    
  out=list("P.map"=P.est.back,"cor.map"=cor.aftermap)
  
  return(out)
}

match.comb<-function(P.est,P.true,mks_known,mks_find){
  
  cor.m1=match.cor(P.est, P.true)
  cor.m2=match.cor.MKS(P.est,P.true,mks_known,mks_find)
  
  return(list("cor.m1"=cor.m1,"cor.m2"=cor.m2))
}



##map back cell types using markers for each cell type##
map.MKS<-function(mks_in,mks_true){

}


map.realdata.fun<-function(mks_known,mks_find,P_in){
  ##try to map back##
  
  #mks_known=MKS.list.all;mks_find=ssCAM.pro$MKS;P_in=ssCAM.pro$P.KL
  
  map.cell=NULL
  for(mks in mks_known){
    #mks=mks_known[[1]]
    a=which(unlist(lapply(mks_find,function(x) length(intersect(x,mks))==length(mks)) )) 
    print(a)
    #map.cell=c(map.cell,which(unlist(lapply(mks_find,function(x) length(intersect(x,mks))==length(mks)) )) )
    
  }
  
  
  for(i in 1:length(map.cell)){
    #i=1
    rownames(P_in)[map.cell[i]]=names(mks_known)[i]
  }
  
  return(P_in)
}


MAP.realdata.fun<-function(P.rela_in,P_in){
  
  #P.rela_in=P.rela.use; P_in=CAM.res$P;
  
  packageExist <- require("combinat")
  if (!packageExist) {
    install.packages("combinat")
    library("combinat")
  }
  
  j.temp=1
  P.rela.est=scale(P_in[-j.temp,],center=F,scale=colSums(P_in[-j.temp,]))
  P.rela.est.match=match.cor.mse(P.est=P.rela.est,P.rela_in)
  
  for(j in 2:ncol(P_in) ){
    #print(P.rela.est.match$COR)
    P.rela.est.temp=scale(P_in[-j,],center=F,scale=colSums(P_in[-j,]))
    P.rela.est.match.temp=match.cor.mse(P.est=P.rela.est.temp,P.rela_in)
    
    if(P.rela.est.match.temp$COR>P.rela.est.match$COR){
      P.rela.est.match=P.rela.est.match.temp
      j.temp=j
    }
    
    #print(P.rela.est.match$COR)
  }
  
  P.all.est.order=rbind(P_in[-j.temp,][P.rela.est.match$cell,],P_in[j.temp,])
  rownames(P.all.est.order)=c(rownames(P.rela_in),"Netro")
  P.rela.est.match[["P.all.est.order"]]=P.all.est.order
  
  return(P.rela.est.match)
}

MAP.realdata.mse.fun<-function(P.rela_in,P_in){
  
  #P.rela_in=P.rela.use; P_in=CAM.res$P;
  
  packageExist <- require("combinat")
  if (!packageExist) {
    install.packages("combinat")
    library("combinat")
  }
  
  P.rela.est=scale(P_in[-1,],center=F,scale=colSums(P_in[-1,]))
  P.rela.est.match=match.mse(P.est=P.rela.est,P.rela_in)
  #P.rela.est.order=cbind(P_in[-1,][P.rela.est.match$cell],P_in[-1,])
  #rownames(P.rela.est.order)
  for(j in 2:ncol(P_in) ){
    P.rela.est.temp=scale(P_in[-j,],center=F,scale=colSums(P_in[-j,]))
    P.rela.est.match.temp=match.mse(P.est=P.rela.est.temp,P.rela_in)
    
    
    if(P.rela.est.match.temp$COR>P.rela.est.match$COR){
      P.rela.est.match=P.rela.est.match.temp
    }
    
    #print(P.rela.est.match$COR)
  }
  
  P.rela.est.match
}
