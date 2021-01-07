#setwd("/Users/anchang/Desktop/lung 5 patient data sets/spadeforestmanuscript/dsfpublication")
###setwd("/Users/anchang/Desktop/lung 5 patient data sets/spadeforestmanuscript/dsfpublication/ipsccells/")
###setwd("/Users/anchang/Desktop/lung 5 patient data sets/spadeforestmanuscript/dsfpublication/spermcells/")
###setwd("/Users/anchang/Desktop/lung 5 patient data sets/spadeforestmanuscript/dsfpublication/breastcells/")
####setwd("/Users/anchang/Desktop/lung 5 patient data sets/spadeforestmanuscript/dsfpublication/emtcells/")
###source("/Users/anchang/Desktop/lung 5 patient data sets/feature forest running code.R")
##### BEGIN DSFMIX
time.stamp.dsfmixrun<-system.time( {
source("feature dsfmix prerequisite.R")
###spermdata
#priorg=c("Glis1","Glis2","Glis3","Rorc","Ror1","Rora","Ror2","Rorb","Jazf1","Shbg" ,"Scfd2", "Scfd1","Gata1", "Gata4","Gdnf","Shbg","Scfd2","Scfd1","Gata1", "Gata4")
#priorg=c("Hsf2","Dppa3","Pou5f1","Dnmt3a","Dnmt3b","Sohlh1","Neurog3","Kit","Stra8","Ddx4","Hoxb2","Hmgb2","Lhx1","Id4","Ret","Etv5","Gfra1","Rbmxl2","Asrgl1","Dmrtb1","Catsper3","Prm3","Prm1","Acrv1","Prdm9","Nanog","Sox2","Magea4","Sycp3","Dmc1","Piwil1","Pck2","Acr","Gapdhs","Meioc","Smc3","Top2a","Spata25","Izumo1","Tssk6","Dnajb3","Ovol1","Wt1","Rhox5","Rhox13","Sox8","Foxi1","Egr1","Cyp21a1","Serpina6","Zbtb16","Kdm3a","Pygo2","Brdt","Taf4b","Ar","Vim","Cdh1","Muc1","Twist1","Cd24a","Cd44","Cd14","Igfbp7" , "Apoe" ,"Aard" ,"Fabp9", "Tmsb4x" , "Calm2", "Malat1" , "mt-Atp6","Ctsl" ,"Tmsb4x" ,"Ftl1","Rpl13", "Itm2b" , "Dcn" , "Tmsb4x", "Cst3", "Rps4x" , "Rps4x" ,  "Apoe","Meg3" ,"Tmsb4x" ,"Dcn" ,"Apoe" ,   "Meg3","Mlf1","Sparc" )
#ipscdata
#priorg=c("Sox2","Prrx1","Fbn1","Gata4","Pecam1","Zfp42","Sall4","Pdgfra","Gata6","Fgf4","Tmem92","Zscan4b","Zscan4c","Zscan4f","Foxa2","Zeb2","Gm4340","Gm21761")
priorg=NULL
spadefilename="dsfmixhormone.fcs"
#spadefilename="spadeforestsperm.fcs"
#spadefilename="spadeforestemt.fcs"
#spadefilename="spadeforestbreast.fcs"
#priorg=c("CD324","Vimentin","CD44","MUC1","Twist","CD24")
##load("/Users/anchang/Desktop/lung 5 patient data sets/spermdata/SpermGenesis_24samples_Seurat3_2SD_MT10_Harmony_tsne_UMAP2_rsv.3_SampleDate_20200204.Rdata")
##load("/Users/anchang/Desktop/lung 5 patient data sets/spermdata/Embroy_12samples_Seurat3_MT15_Harmony_tsne_UMAP2_rsv.5_20200324.Rdata")
#setwd("/Users/anchang/Desktop/lung 5 patient data sets/spermdata/external memory/enrichmenttimeplots2/enrichmenttimeplots1/testrun/")
#setwd("/Users/anchang/Desktop/lung 5 patient data sets/spadeforestmanuscript/R files/testemtfeature/testsperm/")
load("dat1.rdata")
load("timeclust.rdata")
load("seuratclust.rdata")
#load("timeseuratclust.rdata")
timeclust=as.numeric(as.factor(timeclust))
seuratclust=as.numeric(as.factor(seuratclust))

poptime=paste(timeclust,seuratclust,sep=" ")
timeseuratclust=as.numeric(as.factor(poptime))
save(timeseuratclust,file="timeseuratclust.rdata")

#system.time({
if (dim(dat1)[1]>15000){
    
    x=1:dim(dat1)[1]
    n = 2
    chunk=split(x, sort(x%%n))
D1=list()
D1[[1]]=dat1[chunk[[1]],]
D1[[2]]=dat1[chunk[[2]],]
#D1[[3]]=dat1[chunk[[3]],]
#D1[[4]]=dat1[chunk[[4]],]
#dd2=NULL
#for ( k in 1:2) {
#dd1=apply(D1[[k]],1,mad5)
#dd2=c(dd2,dd1)
#}
#names(dd2)=rownames(dat1)

##system.time({
sq<-1:2
bootfxn2<-function(x) {
    sqdata=D1[[x]]
    bootfxn1<-function(y) mad5(sqdata[y,])
    trials <- seq(1, dim(sqdata)[1])
    dd3 <- mclapply(trials, bootfxn1, mc.cores = 7)
    }
    dd1 <- mclapply(sq, bootfxn2, mc.cores = 2)
    dd2=unlist(dd1)
   ##}
    names(dd2)=rownames(dat1)

} else {
    #dd2=apply(dat1,1,mad5)
    trials <- seq(1, dim(dat1)[1])
    bootfxn<-function(x) mad5(dat1[x,])
    dd1 <- mclapply(trials, bootfxn, mc.cores = 7)
    dd2=unlist(dd1)
    names(dd2)=rownames(dat1)
}
#}

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
##res=testsymmodspread(X,dat=XX,alpha=alphac,multimodset=NULL)
res=testsymmodspread2(X,dat=XX,alpha=alphac,multimodset=NULL)
save(res,file="res.rdata")
}else{
    alphac=0.0001
    ###res=testsymmodspread(X=XX,dat=XX,alpha=alphac,multimodset=NULL)
    res=testsymmodspread2(X,dat=XX,alpha=alphac,multimodset=NULL)
    save(res,file="res.rdata")
}

### top genes
cutg=1000
symres2=unlist(res)
if(length(symres2)>cutg) {
    ll=unlist(lapply(res,function(x) length(x)))
    ##maxll=which.max(ll)
    sortll=sort(ll,decreasing=T)
    ##ll[maxll]=sortll[2]
        cut= sortll[2]
  #  load("unimod3.rdata")
  #res=testsymmodspread(XX,dat=XX,alpha=alphac,multimodset=unimod3)
  #symres2=unlist(res)
  
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
###res2=testsymmodspread(X=XX,alpha=0.05,symmetryset=symres2)
##gn=lapply(res,function(x) while length(x)>0 ; names(sort(x,decreasing=FALSE)[1:min(length(x),cut)]))
for ( k in 1:length(res)) {
gg=names(res[[k]])
#g=colnames(X[,1:20])
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
#trees=list()
#topvargenelist=list()
#D22=as.data.frame(t(D2))
#D33=as.data.frame(t(D3))
#D44=as.data.frame(t(D4))
# trees[[1]] <- ctree(formula=as.factor(groups) ~ ., data = D22,controls = ctree_control(maxdepth = 2,minbucket=500))
 #trees[[2]]<- ctree(formula=as.factor(groups) ~ ., data = D33,controls = ctree_control(maxdepth = 2,minbucket = 500))
#trees[[3]]<- ctree(formula=as.factor(groups) ~ ., data = D44,controls = ctree_control(maxdepth = 2,minbucket = 500))
##tree2 <- ctree(formula=as.factor(groups) ~ ., data = D2,controls = ctree_control(maxdepth = 6,minbucket = 100)) wrong tree

#for ( k in 1:3) {
#tree_nodes = nodes(trees[[k]], 1:max(unique(where(trees[[k]]))))
#out1 <- sapply(tree_nodes, function(i) {
#    splitstat <- i$psplit$splitstatistic
#    x <- i$psplit$variableName
 #   return(x)
#}
#)
#print(head(unique(unlist(out1))))

#   topvargenelist[[k]]<- unlist(out1)
#}

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
save(D44,file="D44_group*time.rdata")
##save(trees,file="trees.rdata")

#for (j in 1:3){
#tiff(paste("Tree",j,"with top var genes across all clusters.tiff",sep=""), height = 15, width = 20,units = 'cm', pointsize = 6,compression = "lzw",type="cairo", res=300)
#plot(trees[[j]],xlab="groups",main=topvargenelist[[j]][1])
#dev.off();
#}
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
#sperm_ids=unique(spadegenes,priorg)
sperm_ids=spadegenes
###sperm_ids=c("CD324","Vimentin","CD44","CD24","Twist","MUC1" )
#### Identity transformation function. Could be changed by user
transforms1= function(x) x

### input file path
input <- paste(getwd(),spadefilename,sep="/")
#input <- paste(getwd(),"spadeforestemt.fcs",sep="/")
###output file path
output_dir=(paste(getwd(),"output",sep="/"))
#forestclust=list(c(1,2,3,4),c(5,6,7))
##forestclust=list(c(1:12),c(13:17),c(18:25))
##forestclust=list(c(1:12),c(13:16),c(17:25))
##forestclust=list(c(1),c(2), c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11),c(12))
##forestclust=list(c(1,2),c(3,4), c(5,6),c(7,8),c(9),c(10,11),c(12),c(13),c(14),c(15,16),c(17),c(18,19),c(20,21),c(22),c(23),c(24,25))
####spermdays=c(5,5.05,25,15,35,30,20,10,67,80.05,80,14,18,18.05,25.05,30.05,6,8,0,0.5,3,3.05,6.05,7,8.05)
##### Run SPADE FOREST on Lung mouse data example
### See the description of variables from the integrated.tool.functionstree.R file
#### All the functions needed to run spade forest are in integrated.tool.functionstree.R file.
## set seed for reproducibility of figures
### time.stamp.spadeforestrun is to record running time
#source("/Users/anchang/Desktop/lung 5 patient data sets/spermdata/integrated.tool.functionstree4.R")
##igraph:::layout_in_circle
#igraph:::layout.kamada.kawai
#k=200
k=100
set.seed(123)
dsf<-dsfmixmain(files = input, file_pattern = "*.fcs", out_dir = output_dir, cluster_cols = sperm_ids, comp = FALSE,transforms = transforms1, downsampling_target_number = NULL, downsampling_target_pctile = 1, downsampling_target_percent = 1, downsampling_exclude_pctile = 0, k = k, graph_cols = NULL,pctile_color = c(0, 1), run_spade=TRUE, do_real_filtering = FALSE, run_spade_only=FALSE, gated=FALSE, forestcluster=NULL,dcut = NULL,leafsort=TRUE, minsize=20, remove_single_cell_clusters=FALSE, perplexity = NULL, costtol=10,reg=NULL,forestlayout=igraph:::layout_in_circle)
}, gc=TRUE)
save(time.stamp.dsfmixrun, file="./time.stamp.dsfmixrun.rda")
#)
##### THE END OF DSFMIX


