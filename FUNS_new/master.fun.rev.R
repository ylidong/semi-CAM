master.fun.rev<-function(data.mix.use_in,MKS.list.use_in,P.use_in){
  ##################
  ###CAM original###
  ##################
  ##avg estimate##
  ncell=length(MKS.list.use)
  
  ##CAM original##
  CAM.res=CAM.main.est(data_mix=data.mix.use,ncell,cluster_num=50)
  CAM.cor=match.cor(P.est=CAM.res$P,P.use)
  RES.CAM.org=data.frame(method="CAM",known.vertice=0,known.markers.pct=0,vertice.select=0,COR.KL=CAM.cor$COR,COR.Fro=CAM.cor$COR)
  
  ##CAM+NMF##
  CAM.NMF.res=main.NMF.1run(data.mix.use,CAM.res$MKS)
  CAM.KL=match.cor(CAM.NMF.res$P.KL,P.use)
  CAM.Fro=match.cor(CAM.NMF.res$P.Fro,P.use)
  
  RES.CAM.NMF=data.frame(method="CAM.NMF",known.vertice=0,known.markers.pct=0,vertice.select=0,COR.KL=CAM.KL$COR,COR.Fro=CAM.Fro$COR)
  
  ##CAM+NMFMKS##
  CAM.NMFMKS.res=main.NMF.1run(data.mix.use[unlist(CAM.res$MKS),],CAM.res$MKS)
  CAM.NMFMKS.KL=match.cor(CAM.NMFMKS.res$P.KL,P.use)
  CAM.NMFMKS.Fro=match.cor(CAM.NMFMKS.res$P.Fro,P.use)
  
  RES.CAM.NMFMKS=data.frame(method="CAM.NMFMKS",known.vertice=0,known.markers.pct=0,vertice.select=0,COR.KL=CAM.NMFMKS.KL$COR,COR.Fro=CAM.NMFMKS.Fro$COR)
  
  ##CAM+DSA##
  CAM.DSA.res= ged(as.matrix(data.mix.use), MarkerList(CAM.res$MKS), "DSA")
  CAM.DSA=match.cor(CAM.DSA.res@fit@H,P.use)
  RES.CAM.DSA=data.frame(method="CAM.DSA",known.vertice=0,known.markers.pct=0,vertice.select=0,COR.KL=CAM.DSA$COR,COR.Fro=CAM.DSA$COR)
  
  RES.CAM=rbind(RES.CAM.org,RES.CAM.NMF,RES.CAM.NMFMKS,RES.CAM.DSA)  
  
  ######################
  ####MY method####
  ######################
  
  RES.MY=NULL
  ##force estimate##
  for(known.vertice in 1:ncell){
  #for(known.vertice in 3){  
  
    vertice.select=matrix(combn(ncell,known.vertice),nrow=known.vertice)
    
    
    for(j in 1:ncol(vertice.select)){
      
      MKS.initial=MKS.list.use[vertice.select[,j]]
      print(vertice.select[,j])
      
      ssCAM.res=ssCAM.main.NMF.1run(data_ssCAM=data.mix.use,data_est=data.mix.use,ncell,cluster_num=50,mks_in=MKS.initial)
      ssCAM.KL=match.cor(ssCAM.res$P.KL,P.use)
      ssCAM.Fro=match.cor(ssCAM.res$P.Fro,P.use)
      RES.ssCAM.NMF=data.frame(method="ssCAM.NMF",known.vertice=known.vertice,known.markers.pct=known.markers.pct,vertice.select=paste0(as.vector(vertice.select[,j]),collapse = ""),COR.KL=ssCAM.KL$COR,COR.Fro=ssCAM.Fro$COR)
      
      ##ssCAM+NMFMKS
      ssCAM.NMFMKS.res=main.NMF.1run(data.mix.use[unlist(ssCAM.res$MKS),],ssCAM.res$MKS)
      ssCAM.NMFMKS.KL=match.cor(ssCAM.NMFMKS.res$P.KL,P.use)
      ssCAM.NMFMKS.Fro=match.cor(ssCAM.NMFMKS.res$P.Fro,P.use)
      RES.ssCAM.NMFMKS=data.frame(method="ssCAM.NMFMKS",known.vertice=known.vertice,known.markers.pct=known.markers.pct,vertice.select=paste0(as.vector(vertice.select[,j]),collapse = ""),COR.KL=ssCAM.NMFMKS.KL$COR,COR.Fro=ssCAM.NMFMKS.Fro$COR)
      
      ##ssCAM+DSA##
      ssCAM.DSA.res= ged(as.matrix(data.mix.use), MarkerList(ssCAM.res$MKS), "DSA")
      ssCAM.DSA=match.cor(ssCAM.DSA.res@fit@H,P.use)
      RES.ssCAM.DSA=data.frame(method="ssCAM.DSA",known.vertice=known.vertice,known.markers.pct=known.markers.pct,vertice.select=paste0(as.vector(vertice.select[,j]),collapse = ""),COR.KL=ssCAM.DSA$COR,COR.Fro=ssCAM.DSA$COR)
      
      out=rbind(RES.ssCAM.NMF,RES.ssCAM.NMFMKS,RES.ssCAM.DSA)
      
      RES.MY=rbind(RES.MY,out)
    }
  }
  #}
  
  ##NMF##
  NMF.res=main.NMF.1run(data_est=data.mix.use,mks_in=MKS.list.use)
  NMF.COR.KL=match.cor(P.est=NMF.res$P.KL,P.use)
  NMF.COR.Fro=match.cor(P.est=NMF.res$P.Fro,P.use)
  RES.NMF=data.frame(method="NMF",known.vertice=ncell,known.markers.pct=known.markers.pct,vertice.select="123",COR.KL=NMF.COR.KL$COR,COR.Fro=NMF.COR.Fro$COR)
  
  ##NMFMKS##
  NMFMKS.res=main.NMF.1run(data_est=data.mix.use[unlist(MKS.list.use),],mks_in=MKS.list.use)
  NMFMKS.COR.KL=match.cor(P.est=NMFMKS.res$P.KL,P.use)
  NMFMKS.COR.Fro=match.cor(P.est=NMFMKS.res$P.Fro,P.use)
  RES.NMFMKS=data.frame(method="NMFMKS",known.vertice=ncell,known.markers.pct=known.markers.pct,vertice.select="123",COR.KL=NMF.COR.KL$COR,COR.Fro=NMF.COR.Fro$COR)
  
  
  ##DSA##
  DSA.res= ged(as.matrix(data.mix.use), MarkerList(MKS.list.use), "DSA")
  P.DSA=DSA.res@fit@H
  DSA.COR=match.cor(P.est=P.DSA,P.use)
  RES.DSA=data.frame(method="DSA",known.vertice=ncell,known.markers.pct=known.markers.pct,vertice.select="123",COR.KL=DSA.COR$COR,COR.Fro=DSA.COR$COR)
  
  RES.4MKS=rbind(RES.NMF,RES.NMFMKS,RES.DSA)

  ##summary of all results##
  RES.all=rbind(RES.CAM,RES.MY,RES.4MKS)  
  return(RES.all)
}