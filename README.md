# semi-CAM
A semi-supervised deconvolution method for bulk transcriptomic data with partial marker gene information


##Input data: \

1. Gene expression data of mixture samples. semi-CAM can take around 10 mixtures, if have more than 15 mixtures, dimension reduction methods such as hierarchical clustering methods need to be applied on the data first. 

2. Markers for partial/all cell types. The algorithm takes markers for at least one cell type, up to markers for all number of cell types in the mixtures.

3. Number of cell types in the mixtures. 


# Example code

Deconvolute two tissues from the mixture tissue data GSE29832

###1. Read in the gene expression data of the mixtures

```
dir="http://raw.github.com/ylidong/semi-CAM/master/Test data/"
data=read.csv(paste0(dir, "GSE29832_mixture.csv"),row.names=1)
```

###2. Read in the marker list

Markers list is a R list with each element contains the marker genes for one cell/tissue

```
RData_githubURL <- "http://raw.github.com/ylidong/semi-CAM/master/Test data/MKS_GSE29832_100.RData"
load(url(RData_githubURL))

print(MKS.list.order) ##print the markers for each blood and breast
```

###3. Source the R code
```
RSource_githubURL <- "http://raw.github.com/ylidong/semi-CAM/master/FUNS_new/"

file.sources = list.files(path=url(RSource_githubURL))
sapply(file.sources,function(x) source )

```

###4. True proportions##
```
P.true=cbind(c(0.67,0.33),c(0.33,0.67),c(0.33,0.67))
colnames(P.true)=paste0("P",1:ncol(P.true))
rownames(P.true)=paste0("C",1:nrow(P.true))
P.true=P.true[,rep(1:3,each=3)]
```
###5. Use semi-CAM to estimate cell proportions

```
ssCAM.res=ssCAM.main.NMF.1run(data_ssCAM=data.mix.use,data_est=data.mix.use,ncell,cluster_num=50,mks_in=MKS.initial)
  
  
```