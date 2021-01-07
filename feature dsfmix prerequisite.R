library(Seurat)
library(xgboost)
library(spade) ##install from the source spade_1.0.0.tar R package
library(car)
library(scales)
library(flowCore)
library(RColorBrewer)
library(cluster)
library(fastcluster)
library(sp)
library(igraph)
library(dendextend)
library(ggplot2)
library(plot3D)
library(lawstat)
library(diptest)
library(mutoss)
library(robustbase)
library(uwot)
library(Rtsne)
library(phateR)
library(party)
library(genie)
library(dendsort)
library(networkD3)
library(magrittr)
library(htmlwidgets)
library(plyr)
library(tidyr)
library(dynamicTreeCut)
library(gTests)
library(parallel)
###COMMENT on installing spade (spade_1.0.0.tar)##
#installing spade algorithm (spade_1.0.0.tar) might produce C++ clang related errors. I recommend open terminal and run the following:
#cat <<- EOF > ~/.R/Makevars
#C=/usr/bin/clang
#CXX=/usr/bin/clang++
#EOF
#R CMD INSTALL spade_1.0.0.tar
#### Load all functions
DoubleMAD <- function(x, zero.mad.action="warn"){
   # The zero.mad.action determines the action in the event of an MAD of zero.
   # Possible values: "stop", "warn", "na" and "warn and na".
   x         <- x[!is.na(x)]
   m         <- median(x)
   abs.dev   <- abs(x - m)
   left.mad  <- median(abs.dev[x<=m])
   right.mad <- median(abs.dev[x>=m])
   if (left.mad == 0 || right.mad == 0){
      if (zero.mad.action == "stop") stop("MAD is 0")
      if (zero.mad.action %in% c("warn", "warn and na")) warning("MAD is 0")
      if (zero.mad.action %in% c(  "na", "warn and na")){
         if (left.mad  == 0) left.mad  <- NA
         if (right.mad == 0) right.mad <- NA
      }
   }
   return(c(left.mad, right.mad))
}

DoubleMADsFromMedian <- function(x, zero.mad.action="warn"){
   # The zero.mad.action determines the action in the event of an MAD of zero.
   # Possible values: "stop", "warn", "na" and "warn and na".
   two.sided.mad <- DoubleMAD(x, zero.mad.action)
   print(two.sided.mad)
   m <- median(x, na.rm=TRUE)
   x.mad <- rep(two.sided.mad[1], length(x))
   x.mad[x > m] <- two.sided.mad[2]
   mad.distance <- abs(x - m) / x.mad
   mad.distance[x==m] <- 0
   return(mad.distance)
}

mad5<-function(x) {
    if(all(x==0)) {
    c=0
    } else {
        r1=which(x==0)
        if(length(r1)>0) {
        lr1=length(x[-r1])
         x1=x[-r1]
         if(lr1>=30) {
            a=DoubleMAD(x1, zero.mad.action="warn")
            if(all(a==0)) {
                c=0
            }else{
                c=Qn(x1, constant = 1,finite.corr = FALSE)
            }
         }else{
          c=0
         }
        } else {
            x1=x
               a=DoubleMAD(x1, zero.mad.action="warn")
               if(all(a==0)) {
                   c=0
               }else{
                   c=Qn(x1,constant = 1,finite.corr = FALSE)
               }
        }
    }
    y=c
    print(paste("dispersion=",c))
    return(y)

}

symtestp<-function(x,side="both"){
r1=which(x==0)
if(length(r1)>0){
r=symmetry.test(x[-r1], side=side, boot = FALSE)[[4]]
} else{
  r=symmetry.test(x, side=side, boot = FALSE)[[4]]
}
return(r)
}
symteststat<-function(x,side="both"){
r1=which(x==0)
if(length(r1)>0){
r=symmetry.test(x[-r1], side=side, boot = FALSE)[[3]]
} else{
  r=symmetry.test(x, side=side, boot = FALSE)[[3]]
}
return(r)
}
multimodtest<-function(x, stat="pval"){
r1=which(x==0)
if(length(r1)>0){
    if(stat=="pval"){
r=dip.test(x[-r1])$p.value
    }else{
  r=dip.test(x[-r1])$statistic
    }
} else{
    if(stat=="pval"){
       r=dip.test(x)$p.value
    }else{
     r=dip.test(x)$statistic
    }
}
return(r)
}
#### Simes test function
simes<-function (x, returnstat = FALSE,alpha=0.05)
{
    r = rank(x)
    T = min(length(x) * x/r)
    ###cutoff=which(c(length(x) * x/r)==T)
    if (returnstat)
    c(T, alpha/length(x))
    else T
}

###Bonferonni test function
bonfer<-function(x,returnstat = FALSE,alpha=0.05){
    require(graphics)
    p1=p.adjust(x, method = "bonferroni")
    T = length(x)* min(p1)
    if (returnstat)
    c(T, alpha/length(x))
    else T
    
}

make.R.matrixtop<-function (dat, wt, pi1 = 0.01,cutoff)
{
    ###require(limma)
    cn <- colnames(dat)
    cn <- setdiff(cn, wt)
    colnames(dat) <- make.names(colnames(dat))
    exps <- factor(colnames(dat))
    design <- model.matrix(~0 + factor(colnames(dat)))
    fit1 <- lmFit(dat, design)
    contrast.matrix <- makeContrasts(contrasts = paste(setdiff(exps,
    wt), "-", wt, sep = ""), levels = exps)
    fit2 <- contrasts.fit(fit1, contrast.matrix)
    fit3 <- eBayes(fit2,proportion=pi1)
    lods <- fit3$lods
    tscore<-fit3$t
    pscore<-fit3$p.value
    pdiff <- exp(lods)/(1+exp(lods))
    out=topTable(fit3,p.value=cutoff,number=Inf)
    ##cut=grep("TRUE",(out$adj.P.Value>cutoff))
    #out=topTable(fit3,number=15)
    print(dim(out))
    colnames(lods) <- cn
    res<-list(logodds=lods,tstat=tscore,pvalue=pscore,top=out,prob=pdiff)
    return(res)
}

enrichment<-function(D,xx,removerow=NULL,cc1,cutoff=0.05) {
##Enrichment analysis###
cid=which(rownames(D)=="groups")
s1=which(D[cid,] %in% cc1[xx])
if(!is.null(removerow)) {
    dat3=D[-removerow,s1]
    dat3b=D[-removerow,]
}else {
    dat3=D[,s1]
    dat3b=D
}

require(limma)
dat00=dat3b[-cid,]
colnames(dat00)=paste(D[cid,])
cc=ifelse(colnames(dat00)!=cc1[xx],paste("Non",cc1[xx],sep=""),cc1[xx])
colnames(dat00)=cc
control=paste("Non",cc1[xx],sep="")
toptable=make.R.matrixtop(dat00, wt= control, pi1 = 0.01,cutoff=cutoff)
print(dim(toptable[[4]]))
toptable1=toptable[[4]]
toptable2=toptable[[2]]
posttid=which(toptable2>0)
toptable3=intersect(rownames(toptable1),names(toptable2[,1])[posttid])
print(length(toptable3))
save(toptable,file=paste("toptable_",cc1[xx],".rdata",sep=""))

pdf(file=paste("Top enriched genes  for celltype",xx,"in sample ",cc1[xx],".pdf",sep=""),width=20, height=15)
if(length(toptable3)<10) {
    posttid2=which(toptable2<0)
    toptable4=intersect(rownames(toptable1),names(toptable2[,1])[posttid2])
    topgenes=toptable4
}else{
topgenes=toptable3
}
if(length(topgenes)>2) {
if(dim(dat00)[2]<1000) {
    heatmap(dat00[topgenes,,drop=FALSE],Rowv = NA,Colv =NA, scale="none",col= colorpalette,xlab="CELLS",ylab="GENES",cex.lab=1.8, main=paste("Top enriched genes for cluster", cc1[xx]),margins = c(8, 8))
    dev.off();
} else {
    rr1=sample(s1, min(length(s1),1000),replace=FALSE)
    rr2=sample((1:dim(dat00)[2])[-s1], 1000,replace=FALSE)
    rr=c(rr1,rr2)
    dat50=dat00[topgenes,rr,drop=FALSE]
    heatmap(dat50,Rowv = NA,Colv =NA,scale="none",col= colorpalette,xlab="CELLS",ylab="GENES",cex.lab=1.8, main=paste("Top enriched genes for cluster", cc1[xx]),margins = c(8, 8))
    dev.off();
   }
}else{
    print(paste("No genes to plot"))
}
return(topgenes)
    
}

##Decision tree Enrichment analysis###
enrichmentctree<-function(D,xx,removerow=NULL,cc1) {

if(!is.null(removerow)) {
    dat3b=D[-removerow,]
}else {
    dat3b=D
}

cid=which(rownames(dat3b)=="groups")
s1=which(dat3b[cid,] %in% cc1[xx])
#dat3b[cid,-s1]=rep(paste("Non",cc1[xx],sep=""),(dim(dat3b)[2]-length(s1)))
dat3b[cid,-s1]=rep(0,(dim(dat3b)[2]-length(s1)))
dat3b[cid,s1]=rep(1,length(s1))
#require(party)
##dat3b[cid,s1]=rep(paste(cc1[xx]),
dat00=as.matrix(t(dat3b))

##cc=ifelse(colnames(dat00)!=cc1[xx],paste("Non",cc1[xx],sep=""),cc1[xx])
##colnames(dat00)=cc
##control=paste("Non",cc1[xx],sep="")
Y=dat00[,cid]
set.seed(123)
if(dim(dat00[,-cid])[1]>30) {
##treeres<- ctree(formula=as.factor(groups) ~ ., data = dat00,controls = ctree_control(maxdepth = 4,minbucket=100))
test <- xgboost(dat00[,-cid], label = Y, objective = "binary:logistic", nrounds = 15, max.depth=3, subsample=0.9)
tree_nodes=xgb.importance(colnames(dat00[,-cid]), model = test )
#tree_nodes = nodes(treeres, 1:max(unique(where(treeres))))
#out1 <- sapply(tree_nodes, function(i) {
 #   splitstat <- i$psplit$splitstatistic
#    x <- i$psplit$variableName
#    return(x)
#}
#)
#print(head(unique(unlist(out1))))
print(head(tree_nodes))
toptable<- tree_nodes
save(toptable,file=paste("toptable_",cc1[xx],".rdata",sep=""))
topgenes=toptable$Feature
pdf(file=paste("Top enriched genes  for celltype",xx,"in sample ",cc1[xx],".pdf",sep=""),width=20, height=15)
if(length(topgenes)>=2) {
if(dim(dat00)[1]<1000) {
    heatmap(dat00[,topgenes,drop=FALSE],Rowv = NA,Colv =NA, scale="none",col= colorpalette,ylab="CELLS",xlab="GENES",cex.lab=1.8, main=paste("Top enriched genes for cluster", cc1[xx]),margins = c(8, 8))
    dev.off();
} else {
    rr1=sample(s1, min(length(s1),1000),replace=FALSE)
    rr2=sample((1:dim(dat00)[1])[-s1], 1000,replace=FALSE)
    rr=c(rr1,rr2)
    dat50=dat00[rr,topgenes,drop=FALSE]
    heatmap(dat50,Rowv = NA,Colv =NA,scale="none",col= colorpalette,ylab="CELLS",xlab="GENES",cex.lab=1.8, main=paste("Top enriched genes for cluster", cc1[xx]),margins = c(8, 8))
    dev.off();
   }
}else{
    print(paste("No genes to plot"))
}
}else{
    topgenes=NULL
  print(paste("Not enough data to fit model and plot genes"))
}
return(topgenes)
    
}

xgboostmulticlass<-function(DD) {
    cid=which(colnames(DD)=="groups")
    groups=DD[,cid]
    dat00=DD[,-cid]
    if(dim(dat00)[1]>30) {
    test <- xgboost(dat00, as.numeric(groups) - 1, objective = "multi:softmax", num_class=length(table(groups)), nrounds = 30, max.depth=3, subsample=0.9)
    ##test <- xgboost(dat00[,-cid], label = Y, objective = "binary:logistic", nrounds = 15, max.depth=3, subsample=0.9)
    tree_nodes=xgb.importance(colnames(dat00), model = test )
    print(head(tree_nodes))
    toptable<- tree_nodes
    ###save(toptable,file=paste("toptable_",cc1[xx],".rdata",sep=""))
    topgenes=toptable$Feature
    }else{
        topgenes=NULL
        print(paste("Not enough data to fit model and plot genes"))
    }
    return(topgenes)
}
    #}
count5test2<-function(testx,medpopx=NULL, medpopy=NULL,alpha=0.05,side="both") {
    sortx=sort(testx,decreasing=FALSE)
    xx=DoubleMADsFromMedian(sortx)
    ###print(table(xx))
     print(range(xx))
    ###zerox=which(xx==0)
    zerox= which.min(xx)
    print(length(zerox))
    x=xx[1:(min(zerox)-1)]
    y=xx[(max(zerox)+1):length(xx)]
    lx=length(x)
    ly=length(y)
    x2=sort(x,decreasing=TRUE)
    ##y1=mad(y,constant = 1)
    print(x2[1])
    y2=sort(y,decreasing=TRUE)
    print(y2[1])
    ##x1=mad(x,constant = 1)
    c1=length(which(x2>y2[1]))
    c2=length(which(y2>x2[1]))
    if(lx==ly) {
            if (side == "both") {
                ALTERNATIVE <- "dispersion left is different from dispersion right"
                c=max(c1,c2)
                p=(1-(phyper(c-1,lx,ly,k=c)))
                    P.VALUE <- p * 2
            }
            if (side == "left") {
                c=c1
                p=(1-(phyper(c-1,lx,ly,k=c)))
                ALTERNATIVE <- "dispersion left is bigger than dispersion right."
                P.VALUE <- p
            }
            if (side == "right") {
                 c=c2
                 p=(1-(phyper(c-1,lx,ly,k=c)))
                ALTERNATIVE <- "dispersion right is bigger than dispersion left."
                P.VALUE <- p
            }
            
        } else {
            
            if (side == "both") {
                ALTERNATIVE <- "dispersion left is different from dispersion right"
                mx=round((log(alpha/2)/(log(lx/(lx+ly)))),0)
                px=(lx/(lx+ly))^mx
                my=round((log(alpha/2)/(log(ly/(lx+ly)))),0)
                py=(ly/(lx+ly))^my
                    c=max(mx,my)
                    if(mx==c) p=(lx/(lx+ly))^c else p=(ly/(lx+ly))^c
                    P.VALUE <- 2*p
            }
            if (side == "left") {
                ALTERNATIVE <- "dispersion left is bigger than dispersion right."
                mx=round((log(alpha)/(log(lx/(lx+ly)))),0)
                px=(lx/(lx+ly))^mx
                c=mx
                P.VALUE <- px
            }
            if (side == "right") {
                ALTERNATIVE <- "dispersion right is bigger than dispersion left."
                my=round((log(alpha)/(log(ly/(lx+ly)))),0)
                py=(ly/(lx+ly))^my
                c=my
                P.VALUE <- py
            }

        }
        return(list(statistic = c, p.value = P.VALUE, alternative = ALTERNATIVE))
}
#######
spreadtest<-function(x,side="both",alpha=0.05) {
###xtest=round(x,0)
xtest=x
#print(table(xtest))
print(range(xtest))
r1=which(xtest==0)
if(length(r1)>0) {
xtest2=xtest[-r1]
x=xtest2
p=count5test2(x,alpha=alpha,side=side)
}else{
    xtest2=xtest
    x=xtest2
    p=count5test2(x,alpha=alpha,side=side)
}
print(paste(side, "stats =",p))
return(list(statistic = p[[1]], p.value = p[[2]], alternative = p[[3]]))
}
############### GENE SELECTION by multimodality,symetry and skewness and CODE #######
testsymmodspread<-function(X,dat,alpha=0.05,multimodset=NULL) {
    
finalgeneset=list()

#### Testing multimodality using list of large data sets ###
if (is.null(multimodset)) {
if (is.list(X)) {
        lsymres1=lapply(X,function(x) apply(x,1,function(x) multimodtest(x,"pval")))
        lsymres2=lapply(X,function(x) apply(x,1,function(x) multimodtest(x,"stats")))
        unimod1=unlist(lsymres1)
        unimod2=unlist(lsymres2)
        unimod3=list(pvalue=unimod1,statistics=unimod2)
        save(unimod3,file="unimod3.rdata")
    }else{
       unimod1=apply(X,1,function(x) multimodtest(x,"pval"))
       unimod2=apply(X,1,function(x) multimodtest(x,"stats"))
       unimod3=list(pvalue=unimod1,statistics=unimod2)
       save(unimod3,file="unimod3.rdata")
    }
} else {
 unimod1=multimodset[[1]]
 unimod2=multimodset[[2]]
}
tiff("Histograms of top genes %01d.tiff", height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
#names(symres)=gnames
#symgsimes=simes(symres,alpha=alpha)
ssa=names(sort(abs(unimod2),decreasing=TRUE))
ssb=names(sort(abs(unimod2),decreasing=FALSE))
unimodgsimes=hochberg(unimod1,alpha=alpha,silent=TRUE)
print(head(unimodgsimes[[1]]))
nounimodgene=intersect(ssa,names(unimod1)[(unimodgsimes$rejected)])
yesunimodgene= intersect(ssb,names(unimod1)[-match(nounimodgene,names(unimod1))])
if (length(nounimodgene)>0) hist(dat[nounimodgene[1],],breaks=50,xlab=nounimodgene[1], main=nounimodgene[1])
if (length(yesunimodgene)>0) hist(dat[yesunimodgene[1],],breaks=50,xlab=yesunimodgene[1], main=yesunimodgene[1])


####Testing multimodality spread ####
#if(length(nounimodgene)>0) {
#Y1=X[nounimodgene,,drop=FALSE]
#nonunimodspreadres=apply(Y1,1,function(x) spreadtest(x,alpha=alpha,side="both")[[2]])
#nonunimodspreadres2=apply(Y1,1,function(x) spreadtest(x,alpha=alpha,side="both")[[1]])

#nonunimodspreadsimes=hochberg(nonunimodspreadres,alpha=alpha,silent=TRUE)
#ss=names(sort(abs(nonunimodspreadres2),decreasing=TRUE))
#print(head(nonunimodspreadsimes[[1]][ss]))
#finalnonunimodspreadgene=intersect(ss,names(nonunimodspreadres)[nonunimodspreadsimes$rejected])
#finalunimodspreadgene2=unimodspreadsimes[finalunimodspreadgene]
#if (length(finalnonunimodspreadgene)>0) hist(X[finalnonunimodspreadgene[1],],breaks=50,xlab=finalnonunimodspreadgene[1], main=finalnonunimodspreadgene[1])
if(length(nounimodgene)>0) {
finalnounimodspreadgene=nounimodgene
}else{
    finalnounimodspreadgene=NULL
}

####Testing assymmetrical unimodal distributions ####
if(length(yesunimodgene)>0) {
Y2=dat[yesunimodgene,,drop=FALSE]
if (dim(Y2)[1]>15000){
    xx=1:dim(Y2)[1]
    n = 2
    chunk1=split(xx, sort(xx%%n))
D=list()
D[[1]]=Y2[chunk1[[1]],]
D[[2]]=Y2[chunk1[[2]],]
lmultimod1=lapply(D,function(x) apply(x,1,function(x) symtestp(x,side="both")))
lmultimod2=lapply(D,function(x) apply(x,1,function(x) symteststat(x,side="both")))
symres=unlist(lmultimod1)
symres2=unlist(lmultimod2)
}else{
symres=apply(Y2,1,function(x) symtestp(x,side="both"))
symres2=apply(Y2,1,function(x) symteststat(x,side="both"))
}
symgsimes=hochberg(symres,alpha=alpha,silent=TRUE)
ss1=names(sort(abs(symres2),decreasing=TRUE))

print(head(symgsimes[[1]][ss1]))
nosymgene=names(symres)[(symgsimes$rejected)]
yesymgene= names(symres)[-match(nosymgene,names(symres))]
if (length(nosymgene)>0) hist(dat[ss1[1],],breaks=50,xlab=ss1[1], main=ss1[1])
if (length(yesymgene)>0) hist(dat[ss1[length(ss1)],],breaks=50,xlab=ss1[length(ss1)], main=ss1[length(ss1)])

}else{
    nosymgene=NULL
    yesymgene=NULL
}

#### testing left skewness unimodality
if(length(nosymgene)>0) {
X1=dat[nosymgene,,drop=FALSE]
leftsymres=apply(X1,1,function(x) symtestp(x,side="left"))
leftsymres2=apply(X1,1,function(x) symteststat(x,side="left"))

leftsymsimes=hochberg(leftsymres,alpha=alpha,silent=TRUE)
ss2=names(sort(leftsymres2,decreasing=FALSE))
yesleftsymgene=intersect(ss2,names(leftsymres)[leftsymsimes$rejected])
rightskewgene=names(leftsymres)[-match(yesleftsymgene,names(leftsymres))]
if (length(yesleftsymgene)>0) hist(dat[yesleftsymgene[1],],breaks=50,xlab=yesleftsymgene[1], main=yesleftsymgene[1])
if (length(rightskewgene)>0) hist(dat[rightskewgene[1],],breaks=50,xlab=rightskewgene[1], main=rightskewgene[1])
}else{
  yesleftsymgene=NULL
  rightskewgene=NULL
}
#### testing left spread unimodality
if(length(yesleftsymgene)>0) {
X11=dat[yesleftsymgene,,drop=FALSE]
leftspreadres=apply(X11,1,function(x) spreadtest(x,alpha=alpha,side="left")[[2]])
leftspreadres2=apply(X11,1,function(x) spreadtest(x,alpha=alpha,side="left")[[1]])

leftspreadsimes=hochberg(leftspreadres,alpha=alpha,silent=TRUE)
ss3=names(sort(leftspreadres2,decreasing=TRUE))
finalleftspreadgene=intersect(ss3,names(leftspreadres)[leftspreadsimes$rejected])
##finalleftspreadgene2=leftspreadsimes$adjPValues[finalleftspreadgene]
if (length(finalleftspreadgene)>0) hist(dat[finalleftspreadgene[1],],breaks=50,xlab=finalleftspreadgene[1], main=finalleftspreadgene[1])
}else{
  finalleftspreadgene=NULL

}
#### testing right skewness unimodality
if(length(rightskewgene)>0) {
Z1=dat[rightskewgene,,drop=FALSE]
rightsymres=apply(Z1,1,function(x) symtestp(x,side="right"))
rightsymres2=apply(Z1,1,function(x) symteststat(x,side="right"))

rightsymsimes=hochberg(rightsymres,alpha=alpha,silent=TRUE)
ss4=names(sort(rightsymres2,decreasing=TRUE))
yesrightsymgene=intersect(ss4,names(rightsymres)[rightsymsimes$rejected])
if (length(yesrightsymgene)>0) hist(dat[yesrightsymgene[1],],breaks=50,xlab=yesrightsymgene[1], main=yesrightsymgene[1])

}else{
  yesrightsymgene=NULL
}
#### testing right skew spread unimodality
if(length(yesrightsymgene)>0) {
X22=dat[yesrightsymgene,,drop=FALSE]
rightspreadres=apply(X22,1,function(x) spreadtest(x,alpha=alpha,side="right")[[2]])
rightspreadres2=apply(X22,1,function(x) spreadtest(x,alpha=alpha,side="right")[[1]])

rightspreadsimes=hochberg(rightspreadres,alpha=alpha,silent=TRUE)
ss5=names(sort(rightspreadres2,decreasing=TRUE))
finalrightspreadgene=intersect(ss5,names(rightspreadres)[rightspreadsimes$rejected])
##finalleftspreadgene2=leftspreadsimes$adjPValues[finalleftspreadgene]
if (length(finalrightspreadgene)>0) hist(dat[finalrightspreadgene[1],],breaks=50,xlab=finalrightspreadgene[1], main=finalrightspreadgene[1])
}else{
  finalrightspreadgene=NULL

}

#### Testing unimodal symetry spread
if(length(yesymgene)>0) {
Y22=dat[yesymgene,,drop=FALSE]
unimodspreadres=apply(Y22,1,function(x) spreadtest(x,alpha=alpha,side="both")[[2]])
unimodspreadres2=apply(Y22,1,function(x) spreadtest(x,alpha=alpha,side="both")[[1]])

unimodspreadsimes=hochberg(unimodspreadres,alpha=alpha,silent=TRUE)
ss6=names(sort(abs(unimodspreadres2),decreasing=TRUE))

print(head(unimodspreadsimes[[1]][ss6]))

finalunimodspreadgene=intersect(ss6,names(unimodspreadres)[unimodspreadsimes$rejected])
#finalunimodspreadgene2=unimodspreadsimes[finalunimodspreadgene]
if (length(finalunimodspreadgene)>0) hist(dat[finalunimodspreadgene[1],],breaks=50,xlab=finalunimodspreadgene[1], main=finalunimodspreadgene[1])

}else{
    finalunimodspreadgene=NULL
}

##### final output
dev.off();
finalgeneset$multimodalspread=finalnounimodspreadgene
finalgeneset$leftskewspread=finalleftspreadgene
finalgeneset$rightskewspread=finalrightspreadgene
finalgeneset$symunimodspread=finalunimodspreadgene
#print(ssa)
#print(ssb)
#return(list(finalgeneset=finalgeneset,finalgenesetpvalues=finalgenesetp))
return(finalgeneset)
}

############### TOP GENE SELECTION by multimodality,symetry and skewness CODE #######
genetestsymmodspread<-function(nn,dat,alpha=0.05) {
#### Testing multimodality ###
finalgeneset=list()
       x=dat[nn,]
       multimod1=multimodtest(x,"pval")
       multimod2=multimodtest(x,"stats")
       if (multimod1<alpha) {
        finalgeneset$shape="multimodal"
        finalgeneset$name=nn
        finalgeneset$pvalue=multimod1
        finalgeneset$stat=multimod2
       } else{
           #### Testing unimodal symmetry ###
           symunimod1=symtestp(x,side="both")
           symunimod2=symteststat(x,side="both")
           if (symunimod1>alpha) {
               symunimodspreadres1=spreadtest(x,alpha=alpha,side="both")[[2]]
               symunimodspreadres2=spreadtest(x,alpha=alpha,side="both")[[1]]
            finalgeneset$shape="symunimodal"
            finalgeneset$name=nn
            finalgeneset$pvalue=symunimodspreadres1
            finalgeneset$stat=symunimodspreadres1
           } else {
                #### Testing unimodal skewness ###
               leftsymres1= symtestp(x,side="left")
               leftsymres2= symteststat(x,side="left")
               if (leftsymres1<alpha) {
                   #### Testing unimodal left skewness ###
                  leftspreadres1=spreadtest(x,alpha=alpha,side="left")[[2]]
                  leftspreadres2=spreadtest(x,alpha=alpha,side="left")[[1]]
               finalgeneset$shape="leftskewed"
               finalgeneset$name=nn
               finalgeneset$pvalue=leftspreadres1
               finalgeneset$stat=leftspreadres2
               } else {
                   #### Testing unimodal right skewness ###
                   rightsymres1= symtestp(x,side="right")
                   rightsymres2= symteststat(x,side="right")
                   if (rightsymres1<alpha) {
                   rightspreadres1=spreadtest(x,alpha=alpha,side="right")[[2]]
                   rightspreadres2=spreadtest(x,alpha=alpha,side="right")[[1]]
                   finalgeneset$shape="rightskewed"
                   finalgeneset$name=nn
                   finalgeneset$pvalue=rightspreadres1
                   finalgeneset$stat=rightspreadres2
               } else {
                       finalgeneset$shape="NULL"
                       finalgeneset$name=nn
                       finalgeneset$pvalue="NULL"
                       finalgeneset$stat="NULL"
                     }
               }
           }
       }
       return(finalgeneset)
}
       
############### GENE SELECTION by multimodality,symetry and skewness and CODE #######
testsymmodspread2<-function(X,dat,alpha=0.05,multimodset=NULL) {
finalgeneset=list()
#### Testing multimodality using list of large data sets ###
if (is.null(multimodset)) {
if (is.list(X)) {
    finalgeneset1=list()
        lsymres=lapply(X,function(x) apply(x,1,function(x) genetestsymmodspread(x,dat,alpha=alpha)))
        finalgeneset1$pvalue= unlist(lapply(lsymres,function(x) unlist(lapply(x,function(x) x$pvalue))))
        finalgeneset1$stat= unlist(lapply(lsymres,function(x) unlist(lapply(x,function(x) x$stat))))
        finalgeneset1$shape= unlist(lapply(lsymres,function(x) unlist(lapply(x,function(x) x$shape))))
        finalgeneset1$name= unlist(lapply(lsymres,function(x) unlist(lapply(x,function(x) x$name))))
        save(finalgeneset1,file="finalgeneset1.rdata")
    }else{
    finalgeneset1=list()
       lsymres=apply(X,1,function(x) genetestsymmodspread(x,dat,alpha=alpha))
       finalgeneset1$pvalue=unlist(lapply(lsymres,function(x) x$pvalue))
       finalgeneset1$stat=unlist(lapply(lsymres,function(x) x$stat))
       finalgeneset1$shape=unlist(lapply(lsymres,function(x) x$shape))
       finalgeneset1$name=unlist(lapply(lsymres,function(x) x$name))
       save(finalgeneset1,file="finalgeneset1.rdata")
    }
} else {
 finalgeneset1=multimodset
}
shapetype=c("multimodal","symunimodal","leftskewed","rightskewed")
tiff("Histograms of top genes %01d.tiff", height = 10, width = 10,units = 'cm',pointsize = 6,compression = "lzw", type="cairo", res=300)
for ( k in 1:4) {
id1=which(finalgeneset1$shape==shapetype[k])
if(length(id1)>0){
pv1=finalgeneset1$pvalue[id1]
names(pv1)=finalgeneset1$name[id1]
pv2=finalgeneset1$stat[id1]
names(pv2)=names(pv1)
ssa=names(sort(abs(pv2),decreasing=TRUE))
ssb=names(sort(abs(pv2),decreasing=FALSE))
adjpvalues=hochberg(pv1,alpha=alpha,silent=TRUE)
print(head(adjpvalues[[1]]))
sig=intersect(ssa,names(pv1[(adjpvalues$rejected)]))
nonsig= intersect(ssb,names(pv1[-match(sig,names(pv1))]))
if (length(nonsig)>0) hist(dat[nonsig[1],],breaks=50,xlab=nonsig[1], main=paste(shapetype[k],"nonsignicant ",nonsig[1]))
if (length(sig)>0) hist(dat[sig[1],],breaks=50,xlab=sig[1], main=paste(shapetype[k], "significant ",sig[1]))
##### final outut
finalgeneset[[shapetype[k]]]=pv1[sig]
}else{
    print(paste("No significant",shapetype[k],"genes"))
}
}
dev.off();
return(finalgeneset)
}

### bivariate non-parametric sampling
bivariatesubsample <- function(n, pmatrix,data1,days,states) {
dt <- expand.grid(X=1:dim(pmatrix)[1], Y=1:dim(pmatrix)[2])
dt$p <- as.numeric(pmatrix/ sum(pmatrix))  # get probabilities
idx <- sample(1:nrow(dt), size=n, replace=TRUE, prob=dt$p)
sampled.x <- dt$X[idx]
sampled.y <- dt$Y[idx]

dat1=NULL
dat2=NULL
newcluster=NULL

statematrix2=matrix(0,nrow=max(days),ncol=max(states))
colnames(statematrix2)=paste("State",1:dim(statematrix2)[2])
rownames(statematrix2)=names(table(days))

if(dim(data1)[1]==length(states)) {
data1=cbind(data1,days)
data1=cbind(data1,states)
rr=rownames(data1)
rownames(data1)=paste(1:dim(data1)[1])
#data2=cbind(data2,days)
#data2=cbind(data2,states)
} else {
    stop("Length of clusters must equal length of time points")
}
rr2=NULL
for ( i in 1:max(days)) {
          for ( j in 1:max(states)) {
         xid=which(sampled.x==i)
         yid=which(sampled.y==j)
         cid=intersect(xid,yid)
         n1=length(cid)
         print(paste("n1=",n1))
         xid2=which(days==i)
         yid2=which(states==j)
         cid2=intersect(xid2,yid2)
         subpop=data1[cid2,,drop=FALSE]
         print(paste("n2=",dim(subpop)[1]))
         if (dim(subpop)[1]>0){
         n2=ifelse(dim(subpop)[1]<n1,dim(subpop)[1],n1)
         subpopid=as.numeric(rownames(subpop)[sample(1:dim(subpop)[1],n2,replace=FALSE)])
         statematrix2[i,j]= n2
         dat=data1[subpopid,]
         rr2=c(rr2,rownames(dat))
         }else{
         statematrix2[i,j]=0
         dat=NULL
         rr2=c(rr2,NULL)
         }
         dat1=rbind(dat1,dat)
         #dat0=data2[subpopid,]
         #dat2=rbind(dat2,dat0)
         print(paste("day-cluster",i,j,dim(dat)))
         
}
}
print(paste("dimension of subsample data",dim(dat1)))

return(list(subdata=dat1,statemat1=statematrix,statemat2=statematrix2,rowids=unique(rr2)))
}

normx<-function(x) (x-min(x))/(max(x) - min(x))

testdispersion<-function(alpha,reg,file,geneids) {
    dat1=read.FCS(file)
    dat2=exprs(dat1)
    cid=which(colnames(dat2)=="cluster")
    cid2=dat2[,cid]
    spreadpval=list()
    spreadstat=list()
    adjspreadpval=list()
    adjspreadstat=list()
    for (i in 1:length(reg)) {
      cid3=cid2%in%reg[[i]]
      X=dat2[cid3,geneids]
    zeros= which(apply(X,2,function(x) all(x==0)))
    if(length(zeros)>0)  {
       X2=X[,-zeros,drop=FALSE]
    }else{
        X2=X
    }
    spreadpval[[i]]=apply(X2,2,function(x) spreadtest(x,alpha=alpha,side="both")[[2]])
    spreadstat[[i]]=apply(X2,2,function(x) spreadtest(x,alpha=alpha,side="both")[[1]])

    pv1=spreadpval[[i]]
    pv2=spreadstat[[i]]
    ssa=names(sort(abs(pv2),decreasing=TRUE))
    adjpvalues=hochberg(pv1,alpha=alpha,silent=TRUE)
    print(head(adjpvalues[[1]]))
    sig=intersect(ssa,names(pv1[(adjpvalues$rejected)]))
    
    adjspreadpval[[i]]=pv1[sig]
    adjspreadstat[[i]]=pv2[sig]
    mat=matrix(c(adjspreadpval[[i]],adjspreadstat[[i]]),ncol=2)
    colnames(mat)=c("AdjPvalues","Stats")
    rownames(mat)=names(adjspreadpval[[i]])
    write.csv(mat,file=paste("pavluestree",i,".csv",sep=""))
    }
   
    return(list(pvalue=adjspreadpval,stat=adjspreadstat))
}

splitegde2<-function(x) {
    ss=unlist(strsplit(x, "|"))
    idx=which("|" == ss)
    before=as.numeric(paste(ss[1:(idx-1)],collapse=""))
    after=as.numeric(paste(ss[(idx+1):length(ss)],collapse=""))
      res1=c(before,after)
    print(paste("edge",res1))
    return(res1)
}

g.tests2<-function (E, sample1ID, sample2ID, test.type = "all", maxtype.kappa = 1.14,
    perm = 0,N=200)
{
    temp = gTests:::getR1R2(E, sample1ID)
    R1 = temp$R1
    R2 = temp$R2
    n = length(sample1ID)
    m = length(sample2ID)
    ##N = n + m
    #N=200
    Ebynode = vector("list", N)
    for (i in 1:N) Ebynode[[i]] = rep(0, 0)
    for (i in 1:nrow(E)) {
        Ebynode[[E[i, 1]]] = c(Ebynode[[E[i, 1]]], E[i, 2])
        Ebynode[[E[i, 2]]] = c(Ebynode[[E[i, 2]]], E[i, 1])
    }
    nE = nrow(E)
    nodedeg = rep(0, N)
    for (i in 1:N) nodedeg[i] = length(Ebynode[[i]])
    nEi = sum(nodedeg * (nodedeg - 1))
    mu0 = nE * 2 * n * m/N/(N - 1)
    mu1 = nE * n * (n - 1)/N/(N - 1)
    mu2 = nE * m * (m - 1)/N/(N - 1)
    V0 = nEi * n * m/N/(N - 1) + (nE * (nE - 1) - nEi) * 4 *
        n * m * (n - 1) * (m - 1)/N/(N - 1)/(N - 2)/(N - 3) +
        mu0 - mu0^2
    V1 = nEi * n * (n - 1) * (n - 2)/N/(N - 1)/(N - 2) + (nE *
        (nE - 1) - nEi) * n * (n - 1) * (n - 2) * (n - 3)/N/(N -
        1)/(N - 2)/(N - 3) + mu1 - mu1^2
    V2 = nEi * m * (m - 1) * (m - 2)/N/(N - 1)/(N - 2) + (nE *
        (nE - 1) - nEi) * m * (m - 1) * (m - 2) * (m - 3)/N/(N -
        1)/(N - 2)/(N - 3) + mu2 - mu2^2
    V12 = (nE * (nE - 1) - nEi) * m * n * (m - 1) * (n - 1)/N/(N -
        1)/(N - 2)/(N - 3) - mu1 * mu2
    S = matrix(c(V1, V12, V12, V2), nrow = 2)
    Zw = ((m) * (R1 - mu1) + (n) * (R2 - mu2))/sqrt((m)^2 * V1 +
        (n)^2 * V2 + 2 * (m) * (n) * V12)
    Zd = (R1 - R2 - (mu1 - mu2))/sqrt(V1 + V2 - 2 * V12)
    if (is.na(match(test.type, c("all", "original", "o", "generalized",
        "g", "weighted", "w", "maxtype", "m")))) {
        cat("Wrong test.type input! All tests are performed!\n")
        test.type = "all"
    }
    if (test.type == "all" || test.type == "original" || test.type ==
        "o") {
        Zo = (nE - R1 - R2 - mu0)/sqrt(V0)
        po.approx = pnorm(Zo)
        ro = list(test.statistic = Zo, pval.approx = po.approx)
    }
    if (test.type == "all" || test.type == "generalized" || test.type ==
        "g") {
        Sinv = solve(S)
        Rmv = c(R1 - mu1, R2 - mu2)
        Zg = Rmv %*% Sinv %*% Rmv
        Zg = Zg[1]
        pg.approx = pchisq(Zg, df = 2, lower.tail = F)
        rg = list(test.statistic = Zg, pval.approx = pg.approx)
    }
    if (test.type == "all" || test.type == "weighted" || test.type ==
        "w") {
        pw.approx = pnorm(-Zw)
        rw = list(test.statistic = Zw, pval.approx = pw.approx)
    }
    if (test.type == "all" || test.type == "maxtype" || test.type ==
        "m") {
        M = max(maxtype.kappa * Zw, abs(Zd))
        pm.approx = 1 - pnorm(M/maxtype.kappa) * (2 * pnorm(M) -
            1)
        rmax = list(test.statistic = M, pval.approx = pm.approx)
    }
    if (perm > 0) {
        Zov = Zgv = Zwv = Mv = rep(0, perm)
        for (k in 1:perm) {
            g = sample(c(sample1ID, sample2ID), n)
            temp.p = gTests:::getR1R2(E, g)
            R1.p = temp.p$R1
            R2.p = temp.p$R2
            Zwv[k] = ((m) * (R1.p - mu1) + (n) * (R2.p - mu2))/sqrt((m)^2 *
                V1 + (n)^2 * V2 + 2 * (m) * (n) * V12)
            Zd.p = (R1.p - R2.p - (mu1 - mu2))/sqrt(V1 + V2 -
                2 * V12)
            if (test.type == "all" || test.type == "original" ||
                test.type == "o") {
                Zov[k] = (nE - R1.p - R2.p - mu0)/sqrt(V0)
            }
            if (test.type == "all" || test.type == "generalized" ||
                test.type == "g") {
                Rmv.p = c(R1.p - mu1, R2.p - mu2)
                Zgv[k] = Rmv.p %*% Sinv %*% Rmv.p
            }
            if (test.type == "all" || test.type == "maxtype" ||
                test.type == "m") {
                Mv[k] = max(maxtype.kappa * Zwv[k], abs(Zd.p))
            }
        }
        if (test.type == "all" || test.type == "original" ||
            test.type == "o") {
            po.perm = length(which(Zov <= Zo))/perm
            ro = c(ro, list(pval.perm = po.perm))
        }
        if (test.type == "all" || test.type == "generalized" ||
            test.type == "g") {
            pg.perm = length(which(Zgv >= Zg[1]))/perm
            rg = c(rg, list(pval.perm = pg.perm))
        }
        if (test.type == "all" || test.type == "weighted" ||
            test.type == "w") {
            pw.perm = length(which(Zwv >= Zw))/perm
            rw = c(rw, list(pval.perm = pw.perm))
        }
        if (test.type == "all" || test.type == "maxtype" || test.type ==
            "m") {
            pm.perm = length(which(Mv >= M))/perm
            rmax = c(rmax, list(pval.perm = pm.perm))
        }
    }
    r = list()
    if (test.type == "all" || test.type == "original" || test.type ==
        "o") {
        r = c(r, list(original = ro))
    }
    if (test.type == "all" || test.type == "generalized" || test.type ==
        "g") {
        r = c(r, list(generalized = rg))
    }
    if (test.type == "all" || test.type == "weighted" || test.type ==
        "w") {
        r = c(r, list(weighted = rw))
    }
    if (test.type == "all" || test.type == "maxtype" || test.type ==
        "m") {
        r = c(r, list(maxtype = rmax))
    }
    return(r)
}


#### END of dependencies ####
