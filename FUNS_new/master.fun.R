master.fun<-function(data.mix.use_in,MKS.list.use_in,P.use_in){
  ##################
  ###CAM original###
  ##################
  ##avg estimate##
  ncell=length(MKS.list.use)
  
  CAM.res=CAM.main.est(data_mix=data.mix.use,ncell,cluster_num=50)
  CAM.cor=match.cor(P.est=CAM.res$P,P.use)
  
  res.CAM.list=list("CAM.cv"=CAM.res)
  #RES.CAM=data.frame(known.vertice=0,known.markers.pct=0,vertice.select=0,res.KL.m1=CAM.cor.cv$cor.m1$COR,res.KL.m2=CAM.cor.cv$cor.m2$cor.map,res.Fro.m1=CAM.cor.nocv$cor.m1$COR,res.Fro.m2=CAM.cor.nocv$cor.m2$cor.map)
  RES.CAM=data.frame(method="CAM",known.vertice=0,known.markers.pct=0,vertice.select=0,COR.KL=CAM.cor$COR,COR.Fro=CAM.cor$COR)
  
  CAM.NMF.res=main.NMF.1run(data.mix.use,CAM.res$MKS)
  CAM.KL=match.cor(CAM.NMF.res$P.KL,P.use)
  CAM.Fro=match.cor(CAM.NMF.res$P.Fro,P.use)
  
  RES.CAM.NMF=data.frame(method="CAM.NMF",known.vertice=0,known.markers.pct=0,vertice.select=0,COR.KL=CAM.KL$COR,COR.Fro=CAM.Fro$COR)
  
  ######################
  ####MY method####
  ######################
  
  RES.MY=NULL
  res.ssCAM.list=list();num.run=0
  ##force estimate##
  for(known.vertice in 1:ncell){
    
    num.run=num.run+1
    vertice.select=matrix(combn(ncell,known.vertice),nrow=known.vertice)
    
    
    for(j in 1:ncol(vertice.select)){
      
      MKS.initial=MKS.list.use[vertice.select[,j]]
      print(vertice.select[,j])
      
      ssCAM.res=ssCAM.main.NMF.1run(data_ssCAM=data.mix.use,data_est=data.mix.use,ncell,cluster_num=50,mks_in=MKS.initial)
      #ssCAM.res=ssCAM.main.NMF.1run(data_ssCAM=data.mix.use,data_est=data.mix.use,ncell,cluster_num=50,mks_in=MKS.initial)
      
      ssCAM.KL=match.cor(ssCAM.res$P.KL,P.use)
      ssCAM.Fro=match.cor(ssCAM.res$P.Fro,P.use)
      
      res.ssCAM.list[[num.run]]=list("vertice.select"=paste0(vertice.select),"ssCAM"=ssCAM.res)
      out=data.frame(method="ssCAM",known.vertice=known.vertice,known.markers.pct=known.markers.pct,vertice.select=paste0(as.vector(vertice.select[,j]),collapse = ""),COR.KL=ssCAM.KL$COR,COR.Fro=ssCAM.Fro$COR)
      
      
      RES.MY=rbind(RES.MY,out)
    }
  }
  #}
  ##NMF##
  
  NMF.res=main.NMF.1run(data_est=data.mix.use,mks_in=MKS.list.use)
  NMF.COR.KL=match.cor(P.est=NMF.res$P.KL,P.use)
  NMF.COR.Fro=match.cor(P.est=NMF.res$P.Fro,P.use)
  
  res.NMF.list=list("NMF"=NMF.res)
  RES.NMF=data.frame(method="NMF",known.vertice=ncell,known.markers.pct=known.markers.pct,vertice.select="123",COR.KL=NMF.COR.KL$COR,COR.Fro=NMF.COR.Fro$COR)
  
  ##DSA##
  DSA.res= ged(as.matrix(data.mix.use), MarkerList(MKS.list.use), "DSA")
  P.DSA=DSA.res@fit@H
  DSA.COR=match.cor(P.est=P.DSA,P.use)
  RES.DSA=data.frame(method="DSA",known.vertice=ncell,known.markers.pct=known.markers.pct,vertice.select="123",COR.KL=DSA.COR$COR,COR.Fro=DSA.COR$COR)
  
  RES.all=rbind(RES.CAM,RES.CAM.NMF,RES.MY,RES.NMF,RES.DSA)
  return(RES.all)
}