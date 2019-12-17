CV.calc<-function(x){
  
  cv.out=(sd(x)/mean(x))*100
  
  return(cv.out)
}

filter.fun<-function(data,low_pct,high_pct,cv,exclude){
  
  norm_data=apply(data,1,function(x) sqrt(sum(x^2)) )
  
  data_sub=data[(norm_data < quantile(norm_data , high_pct) & norm_data > quantile(norm_data,low_pct) ),]
  
  if(cv!=0) {
    dd_cv=apply(data_sub,1,CV.calc)
    data_out_sub=data_sub[(dd_cv>quantile(dd_cv,cv)),]
  }else{
    data_out_sub=data_sub
  }
  
  data_out=data[union(rownames(data_out_sub),exclude),]
   
  return(data_out)
}

filter.ENSG.fun<-function(data_in,mean.low,cv,exclude){
  #data_in=data.mix;cv=0.5;exclude=unlist(MKS.list.order);mean.low=2^7
  
  mean_data=apply(data_in,1,mean )
  
  data_sub=data_in[(mean_data > mean.low ),]
  
  if(cv!=0) {
    dd_cv=apply(data_sub,1,CV.calc)
    data_out_sub=data_sub[(dd_cv>quantile(dd_cv,cv)),]
  }else{
    data_out_sub=data_sub
  }
  
  data_out=data_in[union(rownames(data_out_sub),exclude),]
  
  return(data_out)
}