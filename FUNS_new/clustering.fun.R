##cluster samples for real data##
cluster.samp<-function(data_in,num.k){
  #data_in=Data.sus;num.k=10;
  hc.samp=hclust(as.dist(1-cor(data_in, method="pearson")), method="average")
  hc.cut=cutree(tree=hc.samp,k=num.k)
  
  data.out=NULL
  for(i in unique(hc.cut)){
    data.temp=data_in[,names(hc.cut)[hc.cut==i],drop=F]
    data.temp.mean=apply(data.temp,1,mean)
    data.out=cbind(data.out,data.temp.mean)
  }
  colnames(data.out)=unique(hc.cut)
  return(data.out)
}

##cluster relative proportions for real data##
cluster.relaprop.samp<-function(prop.in,data_in,num.k){
  #prop.in=P.sus.relative.true;data_in=Data.sus;num.k=10;
  hc.samp=hclust(as.dist(1-cor(prop.in, method="pearson")), method="average")
  hc.cut=cutree(tree=hc.samp,k=num.k)
  
  data.out=NULL
  for(i in unique(hc.cut)){
    data.temp=data_in[,names(hc.cut)[hc.cut==i],drop=F]
    data.temp.mean=apply(data.temp,1,mean)
    data.out=cbind(data.out,data.temp.mean)
  }
  colnames(data.out)=unique(hc.cut)
  return(data.out)
}
