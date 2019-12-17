sim.data.fun<- function(Ngene=500,M.gene=50,H2.gene=25,mk_type,P.true,ncell,alpha){
  ###simulate cell type specific expression###
  #Ngene=500   #number of total genes
  #M.gene=50;H2.gene=25;ncell=4;alpha=0.8;H=1;mk_type=1
  
  ##simulate back groud expressions of all the genes##
  mu=runif(Ngene,2^5,2^11)
  #CS.noise=t(sapply(mu,function(x) rnorm(ncell,x,sqrt(x))))
  
  ##Fei changed########################################
  CS1.noise=t(sapply(mu,function(x) rnorm(dim(P.true)[2],x,2^(alpha*log2(x)) )))
  rownames(CS1.noise)=paste0("N",1:nrow(CS1.noise))
  ##########################################################
  
  
  ##simulate marker genes##
  
  ##marker genes 25 for each cell type##
  MK.mat=diag(ncell)%x%rep(1,M.gene)
  CS.mk=runif(nrow(MK.mat),2^8,2^12)*MK.mat
  M.names=lapply(1:ncell,function(x) paste0(rep(paste0("M_",x),each=M.gene),"_",1:M.gene))
  rownames(CS.mk)=unlist(M.names)
  if(mk_type==1) CS.mk[CS.mk==0]=runif(sum(CS.mk==0),2^5,2^6) ##low expression not equal to 0 for markers
  Markers=M.names
  
  ##simulate the two H genes##
  
  if(H2.gene!=0){
    
    for(k in 1:ncell){
      
      for(ww in 1:H2.gene){
        CS.mk[(k-1)*M.gene+ww,c(k,sample(setdiff(c(1:ncell),k),1))]=runif(2,2^8,2^12)
        rownames(CS.mk)[(k-1)*M.gene+ww]=paste0("FM_",k,"_",ww)
      }
      Markers[[k]]=Markers[[k]][-c(1:H2.gene)]
    }
    
  }
  
  ##combine all genes##
  CS.all=CS.mk
  colnames(CS.all)=paste0("Cell",1:dim(CS.all)[2])
  
  ###mixed with proportions###
  M=CS.all%*%P.true
  M.mean=apply(M,1,mean)
  E=sapply(M.mean,function(x) rnorm(dim(M)[2],0,2^(alpha*log2(x)) ))
  Mix=M+t(E)
  
  ##Fei changed##############################
  ##last set of genes should all be noise
  
  Mix=rbind(Mix, CS1.noise)
  #####################################################
  #print(sum(Mix<0))
  
  Mix[Mix<0]=2
  colnames(Mix)=paste0("Mix",1:dim(Mix)[2])
  
  mymix=list(dat_pure=CS.all,dat_mix=Mix,P.true=P.true,Markers=Markers)
  mymix
}


sim.data.fun2<- function(Ngene=500,M.gene=50,H2.gene=25,mk_type,P.true,ncell,alpha){
  ###simulate cell type specific expression###
  #Ngene=500   #number of total genes
  #M.gene=50;H2.gene=25;ncell=4;alpha=0.8;H=1;mk_type=1
  
  ##simulate back groud expressions of all the genes##
  mu=runif(Ngene,2^7,2^12)
  #CS.noise=t(sapply(mu,function(x) rnorm(ncell,x,sqrt(x))))
  
  ##Fei changed########################################
  CS1.noise=t(sapply(mu,function(x) rnorm(dim(P.true)[2],x,2^(alpha*log2(x)) )))
  rownames(CS1.noise)=paste0("N",1:nrow(CS1.noise))
  ##########################################################
  
  
  ##simulate marker genes##
  
  ##marker genes 25 for each cell type##
  MK.mat=diag(ncell)%x%rep(1,M.gene)
  CS.mk=runif(nrow(MK.mat),2^8,2^13)*MK.mat
  M.names=lapply(1:ncell,function(x) paste0(rep(paste0("M_",x),each=M.gene),"_",1:M.gene))
  rownames(CS.mk)=unlist(M.names)
  if(mk_type==1) CS.mk[CS.mk==0]=runif(sum(CS.mk==0),2^5,2^7) ##low expression not equal to 0 for markers
  Markers=M.names
  
  ##simulate the two H genes##
  
  if(H2.gene!=0){
    
    for(k in 1:ncell){
      
      for(ww in 1:H2.gene){
        CS.mk[(k-1)*M.gene+ww,c(k,sample(setdiff(c(1:ncell),k),1))]=runif(2,2^8,2^13)
        rownames(CS.mk)[(k-1)*M.gene+ww]=paste0("FM_",k,"_",ww)
      }
      Markers[[k]]=Markers[[k]][-c(1:H2.gene)]
    }
    
  }
  
  ##combine all genes##
  CS.all=CS.mk
  colnames(CS.all)=paste0("Cell",1:dim(CS.all)[2])
  
  ###mixed with proportions###
  M=CS.all%*%P.true
  M.mean=apply(M,1,mean)
  E=sapply(M.mean,function(x) rnorm(dim(M)[2],0,2^(alpha*log2(x)) ))
  Mix=M+t(E)
  
  ##Fei changed##############################
  ##last set of genes should all be noise
  
  Mix=rbind(Mix, CS1.noise)
  #####################################################
  #print(sum(Mix<0))
  
  Mix[Mix<0]=2
  colnames(Mix)=paste0("Mix",1:dim(Mix)[2])
  
  mymix=list(dat_pure=CS.all,dat_mix=Mix,P.true=P.true,Markers=Markers)
  mymix
}


sim.data.fun3<- function(Ngene=500,M.gene=50,H2.gene=25,mk_type,P.true,ncell,alpha){
  ###simulate cell type specific expression###
  #Ngene=500   #number of total genes
  #M.gene=50;H2.gene=25;ncell=4;alpha=0.8;H=1;mk_type=1
  
  ##simulate back groud expressions of all the genes##
  mu=runif(Ngene,2^7,2^10)
  #CS.noise=t(sapply(mu,function(x) rnorm(ncell,x,sqrt(x))))
  
  ##Fei changed########################################
  #CS1.noise=t(sapply(mu,function(x) rnorm(dim(P.true)[2],x,2^(alpha*log2(x)) )))
  CS1.noise=t(sapply(mu,function(x) rnorm(ncell,x,15 )))
  rownames(CS1.noise)=paste0("N",1:nrow(CS1.noise))
  ##########################################################
  
  
  ##simulate marker genes##
  
  ##marker genes 25 for each cell type##
  MK.mat=diag(ncell)%x%rep(1,M.gene)
  CS.mk=runif(nrow(MK.mat),2^8,2^13)*MK.mat
  M.names=lapply(1:ncell,function(x) paste0(rep(paste0("M_",x),each=M.gene),"_",1:M.gene))
  rownames(CS.mk)=unlist(M.names)
  if(mk_type==1) CS.mk[CS.mk==0]=runif(sum(CS.mk==0),2^5,2^7) ##low expression not equal to 0 for markers
  Markers=M.names
  
  ##simulate the two H genes##
  
  if(H2.gene!=0){
    
    for(k in 1:ncell){
      
      for(ww in 1:H2.gene){
        CS.mk[(k-1)*M.gene+ww,c(k,sample(setdiff(c(1:ncell),k),1))]=runif(2,2^8,2^13)
        rownames(CS.mk)[(k-1)*M.gene+ww]=paste0("FM_",k,"_",ww)
      }
      Markers[[k]]=Markers[[k]][-c(1:H2.gene)]
    }
    
  }
  
  ##combine all genes##
  CS.all=CS.mk
  colnames(CS.all)=paste0("Cell",1:dim(CS.all)[2])
  
  ###mixed with proportions###
  CS.all.all=rbind(CS.all,CS1.noise)
  M=CS.all.all%*%P.true
  M.mean=apply(M,1,mean)
  E=sapply(M.mean,function(x) rnorm(dim(M)[2],0,2^(alpha*log2(x)) ))
  Mix=M+t(E)
  
  ##Fei changed##############################
  ##last set of genes should all be noise
  
  #Mix=rbind(Mix, CS1.noise)
  #####################################################
  #print(sum(Mix<0))
  
  Mix[Mix<0]=2
  colnames(Mix)=paste0("Mix",1:dim(Mix)[2])
  
  mymix=list(dat_pure=CS.all.all,dat_mix=Mix,P.true=P.true,Markers=Markers)
  mymix
}

if(F){
sim.data.real<-function(dat_pure,prop_in,alpha){
  #dat_pure=data.pure.avg;prop_in=P.true;alpha=myalpha
  
  M=dat_pure%*%prop_in
  M.mean=apply(M,1,mean)
  E=sapply(M.mean,function(x) rnorm(dim(M)[2],0,2^(alpha*log2(x)) ))
  Mix=M+t(E)
  
  Mix.new=Mix-min(Mix)+1
  
  mymix=list(dat_pure=dat_pure,dat_mix=Mix.new,P.true=prop_in)
  
  return(mymix)
}
}

if(F){
sim.data.real<-function(dat_pure,prop_in,alpha){
  #dat_pure=data.pure.avg;prop_in=P.true;alpha=myalpha
  
  M=dat_pure%*%prop_in
  #M.mean=apply(M,1,mean)
  E=t(apply(M,1,function(x) rnorm(dim(M)[2],0,2^(alpha*log2(x)) )) )
  #E=t(apply(M,1,function(x) rnorm(dim(M)[2],0,alpha*x )) ) #0.3
  
  Mix=M+E
  
  neg.id=apply(Mix,1,function(x) sum(x<0)>0)
  neg.nm=names(neg.id)[neg.id==T]
  neg.mean=rowMeans(M[neg.nm,])
  
  M[neg.nm[neg.mean>2^7],]#1370496_at
  
  Mix.new=Mix-min(Mix)+1
  
  mymix=list(dat_pure=dat_pure,dat_mix=Mix.new,P.true=prop_in)
  
  return(mymix)
}
}
sim.data.real<-function(dat_pure,prop_in,alpha){
  #dat_pure=data.pure.avg;prop_in=P.true;alpha=myalpha
  
  M=dat_pure%*%prop_in

  E=t(apply(M,1,function(x) rnorm(dim(M)[2],0,2^(alpha*log2(x)) )) )

  Mix=M+E
  
  Mix[Mix<0]=1

  
  mymix=list(dat_pure=dat_pure,dat_mix=Mix,P.true=prop_in)
  
  return(mymix)
}