get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

n_cores <- 1
os = get_os()
if(os!="windows"){
  # try to determine the number of processors and get 1 -less
  if(!require(doParallel)){
    install.packages("doParallel", repos = "http://cran.us.r-project.org", dependencies = TRUE)
    tryCatch({
      install.packages("doParallel", repos = "http://cran.us.r-project.org", dependencies = TRUE);
      library(doParallel)
      n_cores <- detectCores() # Determine the number of cores
    }, error=function(e){cat("doParallel package cannot be installed. The number of cores cannot be determined, hence it is assumed to be 1, n_cores = 1", "\n")})
  } 
  # To not exahust the number of cores.
  if (n_cores > 1)
      n_cores <- n_cores - 1

}

##### BEGIN DSFMIX
time.stamp.dsfmixrun<-system.time( {
source("feature dsfmix prerequisite.R")
###spermdata
priorg=NULL
spadefilename="dsfmixhormone.fcs"

# merge the two matrices
load("dat1_hormonedata.rdata")
load("dat2_hormonedata.rdata")

dat1 <- rbind(dat1,dat2)
load("timeclust_hormonedata.rdata")
load("seuratclust_hormonedata.rdata")
timeclust=as.numeric(as.factor(timeclust))
seuratclust=as.numeric(as.factor(seuratclust))

poptime=paste(timeclust,seuratclust,sep=" ")
timeseuratclust=as.numeric(as.factor(poptime))
save(timeseuratclust,file="timeseuratclust.rdata")

if (dim(dat1)[1]>15000){
    
    x=1:dim(dat1)[1]
    n = 2
    chunk=split(x, sort(x%%n))
D1=list()
D1[[1]]=dat1[chunk[[1]],]
D1[[2]]=dat1[chunk[[2]],]

sq<-1:2
bootfxn2<-function(x) {
    sqdata=D1[[x]]
    bootfxn1<-function(y) mad5(sqdata[y,])
    trials <- seq(1, dim(sqdata)[1])
    dd3 <- mclapply(trials, bootfxn1, mc.cores = n_cores)
    }
    dd1 <- mclapply(sq, bootfxn2, mc.cores = n_cores)
    dd2=unlist(dd1)
    names(dd2)=rownames(dat1)

} else {
    trials <- seq(1, dim(dat1)[1])
    bootfxn<-function(x) mad5(dat1[x,])
    dd1 <- mclapply(trials, bootfxn, mc.cores = n_cores)
    dd2=unlist(dd1)
    names(dd2)=rownames(dat1)
}


g=names(sort(dd2,decreasing = TRUE)[1:20])
save(dd2,file="dd2.rdata")
### plot top mad genes
tiff(paste("top variable genes by mad5",".tiff",sep=""), height = 15, width = 20,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
par(mfrow=c(4,5))
for ( i in 1:20) {
  x=dat1[grep(g[i],rownames(dat1))[1],]
  r1=which(x==0)
  if(length(r1)>0) {
  xx=x[-r1]
  xxx=xx
  hist(xxx,breaks=50,xlab=g[i], main=g[i])
  }else{
      xx=x
      xxx=xx
      hist(xxx,breaks=50,xlab=g[i], main=g[i])
  }
}
dev.off();

dd3=dd2[dd2>0]
sg=names(sort(dd3,decreasing = TRUE))
XX=dat1[sg,]
X=as.matrix(rownames(XX),1)
##### Complete analysis with all data #####
if (dim(XX)[1]>15000){
    xx=1:dim(XX)[1]
    n = 2
    chunk1=split(xx, sort(xx%%n))
X=list()
X[[1]]=as.matrix(rownames(XX[chunk1[[1]],]),1)
X[[2]]=as.matrix(rownames(XX[chunk1[[2]],]),1)
alphac=0.0001
res=testsymmodspread2(X,dat=XX,alpha=alphac,multimodset=NULL)
save(res,file="res.rdata")
}else{
    alphac=0.0001
    res=testsymmodspread2(X,dat=XX,alpha=alphac,multimodset=NULL)
    save(res,file="res.rdata")
}

### top genes
cutg=1000
symres2=unlist(res)
if(length(symres2)>cutg) {
    ll=unlist(lapply(res,function(x) length(x)))
    sortll=sort(ll,decreasing=T)
    cut= sortll[2]

  topg=NULL
  for ( k in 1:length(res)) {
      if(length(res[[k]])>0) {
      topmad=sort(dd2[names(res[[k]])],decreasing=TRUE)[1:min(length(res[[k]]),cut)]
  topg=c(topg,topmad)
      }else{
       topg=c(topg,NULL)
      }
  }
  save(topg,file="topg.rdata")
} else {
  cut=length(symres2)
  topg=NULL
  for ( k in 1:length(res)) {
      if(length(res[[k]])>0) {
      topmad=sort(dd2[names(res[[k]])],decreasing=TRUE)[1:min(length(res[[k]]),cut)]
  topg=c(topg,topmad)
      }else{
       topg=c(topg,NULL)
      }
  }
  save(topg,file="topg.rdata")
}
for ( k in 1:length(res)) {
gg=names(res[[k]])
if(length(gg)>0){
gg=names(res[[k]][1:min(length(res[[k]]),12)])
tiff(paste("top asymmultimodspread genes",gg[k],".tiff",sep=""), height = 15, width = 20,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
par(mfrow=c(3,4))
for ( i in 1:length(gg)) {
  x=XX[gg[i],]
  r1=which(x==0)
  if(length(r1)>0){
  xx=x[-r1]
  }else{
   xx=x
  }
  xxx=xx
  hist(xxx,breaks=50,xlab=gg[i], main=gg[i])
}
} else {
    print("No genes to plot")
}
dev.off();
}


##### Complete SPADE Forest algorithm ########

nrcolors = 100
half = 1 + nrcolors/2
colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
colorpalette = colorRampPalette(colpal)(nrcolors)

D=as.matrix(dat1[names(topg),])
D2=rbind(D,groups=timeclust)
outdir <- paste(getwd(),"/enrichmentplotstime/",sep="")
dir.create(outdir)
setwd(outdir)
save(D2,file="D2.rdata")
cc1=names(table(timeclust))
toptimeenrich=list()

for ( k in 1:max(timeclust)) {
   ##toptimeenrich[[k]]<- enrichment(D=D2,xx=k,removerow=NULL,cc1)
    tryCatch({
   toptimeenrich[[k]]<- enrichmentctree(D=D2,xx=k,removerow=NULL,cc1)
    }, error=function(e){cat("ERROR :","Model is just constant. No features used", "\n")})
}
save(toptimeenrich,file="toptimeenrich.rdata")
setwd("..")
outdir <- paste(getwd(),"/enrichmentplotscluster/",sep="")
dir.create(outdir)
setwd(outdir)
D3=rbind(D,groups=seuratclust)
save(D3,file="D3.rdata")
cc1=names(table(seuratclust))
topclustenrich=list()
for ( k in 1:max(seuratclust)) {
    tryCatch({
   topclustenrich[[k]]<- enrichmentctree(D=D3,xx=k,removerow=NULL,cc1)
   }, error=function(e){cat("ERROR :","Model is just constant. No features used", "\n")})
}
save(topclustenrich,file="topclustenrich.rdata")
setwd("..")
outdir <- paste(getwd(),"/enrichmentplotstimecluster/",sep="")
dir.create(outdir)
setwd(outdir)
D4=rbind(D,groups=timeseuratclust)
save(D4,file="D4.rdata")
cc1=names(table(timeseuratclust))
toptimeseuratclust=list()
for ( k in 1:max(timeseuratclust)) {
    tryCatch({
   toptimeseuratclust[[k]]<- enrichmentctree(D=D4,xx=k,removerow=NULL,cc1)
    }, error=function(e){cat("ERROR :","Model is just constant. No features used", "\n")})
}
save(toptimeseuratclust,file="toptimeseuratclust.rdata")
setwd("..")
###### Variability across classes ####
topvargenelist=list()
D22=as.matrix(t(D2))
D33=as.matrix(t(D3))
D44=as.matrix(t(D4))
set.seed(123)
topvargenelist[[1]] <-xgboostmulticlass(D22)
topvargenelist[[2]]<- xgboostmulticlass(D33)
topvargenelist[[3]]<- xgboostmulticlass(D44)

save(D22,file="D22_time.rdata")
save(D33,file="D33_groups.rdata")
save(D44,file="D44_group_time.rdata")

toggenelist=list()
toggenelist[[1]]=toptimeenrich
toggenelist[[2]]=topclustenrich
toggenelist[[3]]=toptimeseuratclust
save(toggenelist,file="toggenelist.rdata")
save(topvargenelist,file="topvargenelist.rdata")

spadegenes1=NULL
spadegenes2=NULL
for ( i in 1:3) {
    for ( j in 1:length(toggenelist[[i]])) {
       if(length(toggenelist[[i]][[j]])>0) {
           selectg1=unlist(toggenelist[[i]][[j]][1])
            selectg2=unlist(toggenelist[[i]][[j]][1:min(length(toggenelist[[i]][[j]]),20)])
        spadegenes1=c(spadegenes1,selectg1)
        spadegenes2=c(spadegenes2,selectg2)
            }else{
             spadegenes1=c(spadegenes1,NULL)
            spadegenes2=c(spadegenes2,NULL)
            }
    }
}
spadegenes11=NULL
spadegenes22=NULL
for ( i in 1:length(topvargenelist)) {
       if(length(topvargenelist[[i]])>0) {
           selectg1=unlist(topvargenelist[[i]][1])
            selectg2=unlist(topvargenelist[[i]][1:min(length(topvargenelist[[i]]),20)])
        spadegenes11=c(spadegenes11,selectg1)
        spadegenes22=c(spadegenes22,selectg2)
            }else{
             spadegenes11=c(spadegenes11,NULL)
            spadegenes22=c(spadegenes22,NULL)
            }
    }

spadegenes=unique(c(spadegenes1,spadegenes11))
spadedatagenes=unique(c(spadegenes2,spadegenes22))
naid=which(is.na(spadedatagenes))
if(length(naid)>0) {
spadedatagenes=spadedatagenes[-naid]
} else{
  spadedatagenes=spadedatagenes
}
if(is.null(priorg)) {
priorgenes=rownames(dat1)[!(rownames(dat1) %in% spadedatagenes)]
if(length(priorgenes)>500) {
    spadedatagenesf=spadedatagenes
    }else{
spadedatagenesf=unique(c(spadedatagenes,priorgenes))
}
}else{
    priorgenes=priorg
    spadedatagenesf=unique(c(spadedatagenes,priorgenes))
}
save(spadegenes,file="spadegenes.rdata")
save(spadedatagenesf,file="spadedatagenesf.rdata")
##### Generating spade forest FCS file ###
source("feature dsfmix spadefunctions.R")
spadedat=t(as.matrix(dat1[spadedatagenesf,]))
spadedata1=apply(spadedat,2,normx)

spadedata2=cbind(spadedata1,timeclust)
spadedata3=cbind(spadedata2,seuratclust)
spadedata=cbind(spadedata3,timeseuratclust)
spadedata0=spadedata
bintime=list()
bintimeid=paste("timepoint",1:max(timeclust),sep="")
for ( l in 1:max(timeclust)) {
 bintime[[l]]=ifelse(timeclust==l,1,0)
 spadedata=cbind(spadedata,bintime[[l]])
}
colnames(spadedata)=c(colnames(spadedata0),bintimeid)
spadedata00=spadedata
binclust=list()
binclustid=paste("cluster",1:max(seuratclust),sep="")
for ( l in 1:max(seuratclust)) {
 binclust[[l]]=ifelse(seuratclust==l,1,0)
 spadedata=cbind(spadedata,binclust[[l]])
}
colnames(spadedata)=c(colnames(spadedata00),binclustid)

###source("rna.seq.spade.R")
file1="fcstospadedatamatrix.fcs"
file2=spadedata
if(dim(file2)[1]>150000){
    
    ###### subsample for spade forest ######
    days=as.numeric(as.factor(timeclust))
    states=as.numeric(as.factor(seuratclust))
    
statematrix=matrix(0,nrow=max(days),ncol=max(states))
colnames(statematrix)=paste("State",1:dim(statematrix)[2])
rownames(statematrix)=names(table(timeclust))
for (j in 1:max(days) ) {
     for (i in 1:max(states)) {
         ##xid=which(emtdata[,41]==j);
         xid=which(days==j);
         ##yid=which(clustersforcompletedata==i);
         yid=which(states==i);
         cid=intersect(xid,yid);
        n1=length(cid);
        statematrix[j,i]=n1
        }
     }

data1=spadedata
subdat=bivariatesubsample(n=100000, pmatrix=statematrix,data1=data1,days,states)
save(subdat,file="subdat_sperm.rdata")
file3=subdat[[1]]
save(file3,file="file3.rdata")
ddd=matrixtofcs (file1,file3,filename=spadefilename)
}else{
    file3=file2
    save(file3,file="file3.rdata")
    ddd=matrixtofcs (file1,file3,filename=spadefilename)
}

### RUNNING SPADE FOREST #####
# load required R libraries for spade forest

###loading highly enriched and highly variable genes for spade forest

##Selecting the unique genes
sperm_ids=spadegenes
#### Identity transformation function. Could be changed by user
transforms1= function(x) x

### input file path
input <- paste(getwd(),spadefilename,sep="/")
###output file path
output_dir=(paste(getwd(),"output",sep="/"))

##### Run SPADE FOREST on Lung mouse data example
### See the description of variables from the integrated.tool.functionstree.R file
#### All the functions needed to run spade forest are in integrated.tool.functionstree.R file.
## set seed for reproducibility of figures
### time.stamp.spadeforestrun is to record running time
#k=200
k=100
set.seed(123)
dsf<-dsfmixdriverc(files = input, file_pattern = "*.fcs", out_dir = output_dir, cluster_cols = sperm_ids, comp = FALSE,transforms = transforms1, downsampling_target_number = NULL, downsampling_target_pctile = 1, downsampling_target_percent = 1, downsampling_exclude_pctile = 0, k = k, graph_cols = NULL,pctile_color = c(0, 1), run_spade=TRUE, do_real_filtering = FALSE, run_spade_only=FALSE, gated=FALSE, forestcluster=NULL, dcut = NULL,leafsort=TRUE, minsize=20, remove_single_cell_clusters=FALSE, perplexity = NULL, costtol=10,reg=NULL,forestlayout=igraph:::layout_in_circle)
}, gc=TRUE)
save(time.stamp.dsfmixrun, file="./time.stamp.dsfmixrun.rda")
##### THE END OF DSFMIX
