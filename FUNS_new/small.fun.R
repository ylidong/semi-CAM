samp.avg<-function(datdat,num.rep){
  #datdat=dat_mix;num.rep=3
  #datdat=puredata_in[mks.cam[[1]],];num.rep=num.purerep
  
  id_start=seq(1,dim(datdat)[2],by=num.rep)
  id_end=id_start+num.rep-1
  out=NULL
  
  if (num.rep==1){
    out=datdat
  }else{
    for (i in 1:(dim(datdat)[2]/num.rep)){
      out=cbind(out,rowMeans(datdat[,id_start[i]:id_end[i]]))
    }
  }
  out
}


check.markers<-function(find.mks,true.mks){
  #find.mks=CAM.res$MKS[CAM.cor$cell];true.mks=data_0$Markers
  true.mks.pct=NULL
  for(i in 1:length(find.mks)){
    true.mks.pct.temp=round(mean(find.mks[[i]]%in%true.mks[[i]]),3)
    true.mks.pct=c(true.mks.pct,true.mks.pct.temp)
  }
  
  return(paste(true.mks.pct,collapse=","))
}