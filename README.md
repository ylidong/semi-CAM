# semi-CAM
A semi-supervised deconvolution method for bulk transcriptomic data with partial marker gene information


#Input data 

1. Gene expression data of mixture samples. semi-CAM can take around 10 mixtures, if have more than 15 mixtures, dimension reduction methods such as hierarchical clustering methods need to be applied on the data first. 

2. Markers for partial/all cell types. The algorithm takes markers for at least one cell type, up to markers for all number of cell types in the mixtures.

3. Number of cell types in the mixtures. 


# Example code

Deconvolute two tissues from the mixture tissue data mixed from blood and breast tissue

###1. Read in the gene expression data of the mixtures

```
dir="http://raw.github.com/ylidong/semi-CAM/master/Test data/"
data=read.csv(paste0(dir, "Example_mix.csv"),row.names=1)
```

###2. Read in the marker list

Markers list is a R list with each element contains the marker genes for one cell/tissue

```
RData_githubURL <- "http://raw.github.com/ylidong/semi-CAM/master/Test data/Marker.list.RData"
load(url(RData_githubURL))

print(Markers.list) ##print the markers for each blood and breast
```

###3. Source the R code
```
R_githubURL <- "https://raw.github.com/ylidong/semi-CAM/master/FUNS_new/"


source(paste0(source.dir,"proportion_estimate_main.R"))
source(paste0(source.dir,"measure_conv_miss.R"))
source(paste0(source.dir,"identify_markers.R"))
source(paste0(source.dir,"small.fun.R"))
source(paste0(source.dir,"match_celltype.R"))


script <- getURL(paste0(R_githubURL,"match_celltype.R"), ssl.verifypeer = FALSE)
eval(parse(text = script))

```

###4. True proportions##
```
P.true=cbind(c(0.67,0.33),c(0.33,0.67),c(0.33,0.67))
colnames(P.true)=paste0("MIX",1:ncol(P.true))
rownames(P.true)=c("Blood","Breast")
P.true=P.true[,rep(1:3,each=3)]
```
###5. Use semi-CAM to estimate cell proportions

```
set.seed(1234)

###############################
##give markers for two tissue##
###############################
semiCAM.res=semiCAM.main(data_ssCAM=data.mix,data_est=data.mix,ncell=2,cluster_num=50,mks_in=Marker.list)
#output estimated proportions#
semiCAM.res$P

###############################
##only give markers for blood##
###############################

semiCAM.res2=semiCAM.main(data_ssCAM=data.mix,data_est=data.mix,ncell=2,cluster_num=50,mks_in=Marker.list[1])

#match the tissue type using true proportions
semiCAM.match=match.cor(semiCAM.res2$P,P.true)
#output the proportions after matching
semiCAM.match$P.order

```