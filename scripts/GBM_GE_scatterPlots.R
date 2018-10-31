# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
#setwd('~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/clean_july12/scripts/')
# loading results of DESeq models (made by scripts deseq_GE.R and deseq_GBM.R) :

ll=load("../datasets/GBM_splitModels.RData")
ll=load("../datasets/GE_splitModels.RData")
ll=load('../datasets/GBM_simpleModels.RData')
ll=load('../datasets/GE_simpleModels.RData')
# ll=load('../datasets/deseq_exons.RData')
# ll=load('../datasets/deseq_introns.RData')

# MBD-scores (made by MBDscore_calc.R)
ll=load('../datasets/capturedVflowthroughResults.Rdata')
classes=data.frame(cbind("gene"=row.names(mbd)))
classes$mbd.score= mbd$log2FoldChange
plot(density(na.omit(classes$mbd.score)))
abline(v=-1.05)
class.split=-1.05
head(classes)

str(vsd.ge)
goods=intersect(row.names(vsd.ge),row.names(mbd))
length(goods)
ge.means=apply(vsd.ge[goods,],1,mean)
plot(ge.means~mbd[goods,"log2FoldChange"],pch=16, cex=0.3,col=rgb(0,0,0,alpha=0.1))
abline(v=class.split)

#---------- log fold change against MBD score

# un-remark the dataset you want to plot:
 dset=o2k.gbm;YL=c(-1,1.3);YLAB="GBM change";DSCALE=0.6 # Fig 2A, Fig. S8A. 853 DMGs, 65 up, 788 down
# dset=k2o.gbm;YL=c(-1,1);YLAB="GBM change";DSCALE=0.6 # Fig 2B, , Fig. S8B. no DMGs
# dset=natives.gbm;YL=c(-1,1);YLAB="GBM difference";DSCALE=0.6 # Fig S8 C
# dset=ori.s;YL=c(-1,1);YLAB="GBM difference";DSCALE=1;dset$log2FoldChange=(-1)*dset$log2FoldChange # note: inverting the original comparison sign to make it K vs O ; Fig 3a. 379 DMGs, 196 up, 183 down
# dset=o2k.ge;YL=c(-1.5,1.5);YLAB="GE change";DSCALE=1.2 # Fig 2A, Fig. 2F. 239 DEGs, 87 up, 152 down.
# dset=k2o.ge;YL=c(-1.5,1.5);YLAB="GE change";DSCALE=1.8 # Fig 2A, Fig. 2G. 578 DEGs, 344 up, 234 down.
# dset=natives.ge;YL=c(-1.5,1.5);YLAB="GBM difference";DSCALE=1 
# dset=ori.ge;YL=c(-1.5,1.5);YLAB="GE difference";DSCALE=2 # 662 DEGs, O vs K: 321 up, 341 down
summary(dset)
head(dset)
gg=intersect(classes$gene,row.names(dset))

# plotting scatters
change=data.frame(cbind("gene"=row.names(dset[gg,]),"lfc"=dset[gg,]$log2FoldChange,"padj"=dset[gg,]$padj))
change$lfc=as.numeric(as.character(change$lfc))
change$padj=as.numeric(as.character(change$padj))
x=merge(classes,change,by="gene")
plot(lfc~mbd.score,x,col=rgb(0,0,0,alpha=0.01),pch=16,cex=0.5,ylim=YL,xlim=c(-7,5),mgp=c(2.3,1,0),xlab="MBD score",ylab=YLAB)
ll=lm(lfc~mbd.score,x)
summary(ll) #k2o: 0.13; o2k: 0.54 ; ori: 0.023 ; natives: 0.07 ; k2o.ge: 0.013 ; o2k.ge: 0.017
abline(ll,col="red")
abline(v=0,lty=3,col="grey80",lwd=1.5)
abline(h=0,lty=3,col="grey80",lwd=1.5)

# adding density curves for significant genes
siglimit=0.1
xn=na.omit(x)
sigs=(xn$padj<siglimit)
sum(sigs) 
sigsup=(xn$padj<siglimit & xn$lfc>0); sum(sigsup) 
sigsdn=(xn$padj<siglimit & xn$lfc<0); sum(sigsdn) 
#points(lfc~mbd.score,x[sigs,],col=rgb(0.8,0.2,0.2,alpha=0.25),pch=16,cex=0.3)
lines(density(xn[sigsup,"mbd.score"])$x,DSCALE*density(xn[sigsup,"mbd.score"])$y+YL[1],col="red")
lines(density(xn[sigsdn,"mbd.score"])$x, DSCALE*density(xn[sigsdn,"mbd.score"])$y+YL[1],col="blue")

#------------------- scatter o2k ~ k2o

# GBM (Fig. 2C)

plot(k2o.gbm$log2FoldChange~o2k.gbm$log2FoldChange,col=rgb(0,0,0,alpha=0.03),pch=16,cex=0.5,xlim=c(-1,1),ylim=c(-1,1),mgp=c(2.3,1,0),xlab="O > K",ylab="K > O")
abline(v=0,lty=3,col="grey80")
abline(h=0,lty=3,col="grey80")
ll=lm(k2o.gbm$log2FoldChange~o2k.gbm$log2FoldChange)
summary(ll) # 0.09 
abline(ll,col="red")
siglimit=0.1
table(k2o.gbm$padj<=siglimit | o2k.gbm$padj<=siglimit)
sigs=(k2o.gbm$padj<=siglimit | o2k.gbm$padj<=siglimit)

# GE (Fig. 2D)

plot(k2o.ge$log2FoldChange~o2k.ge$log2FoldChange,col=rgb(0,0,0,alpha=0.03),pch=16,cex=0.5,xlim=c(-1,1),ylim=c(-1,1),mgp=c(2.3,1,0),xlab="O > K",ylab="K > O")
abline(v=0,lty=3,col="grey80")
abline(h=0,lty=3,col="grey80")
ll=lm(k2o.ge$log2FoldChange~o2k.ge$log2FoldChange)
summary(ll) # 0.1 
abline(ll,col="red")
siglimit=0.1
table(k2o.ge$padj<=siglimit | o2k.ge$padj<=siglimit)
sigs=(k2o.ge$padj<=siglimit | o2k.ge$padj<=siglimit)

#------------------- GBM ~ GE

# o2k (Fig 2F)

goods=intersect(row.names(o2k.ge),row.names(o2k.gbm))
length(goods)
o2k.ge=o2k.ge[goods,]
o2k.gbm=o2k.gbm[goods,]
summary(o2k.gbm)
siglimit=1
sigs=(o2k.gbm$padj<=siglimit)
plot(o2k.ge$log2FoldChange[sigs]~o2k.gbm$log2FoldChange[sigs],col=rgb(0,0,0,alpha=0.05),pch=16,cex=0.5,xlab="GBM",ylab="GE",mgp=c(2.3,1,0),xlim=c(-0.35,0.5),ylim=c(-1.5,1.5))
abline(v=0,col="grey70",lty=3)
abline(h=0,col="grey70",lty=3)
lmm=lm(o2k.ge$log2FoldChange~o2k.gbm$log2FoldChange)
abline(lmm,col="red")
summary(lmm) # R2=0.009
siglimit=0.1
summary(o2k.gbm)
table(o2k.gbm$padj<=siglimit)
sigs=(o2k.gbm$padj<=siglimit)
sum(sigs)
points(o2k.ge$log2FoldChange[sigs]~o2k.gbm$log2FoldChange[sigs],col=rgb(0.8,0.2,0.2,alpha=0.25),pch=16,cex=0.5)

# k2o (Fig 2G)

goods=intersect(row.names(k2o.ge),row.names(k2o.gbm))
k2o.ge=k2o.ge[goods,]
k2o.gbm=k2o.gbm[goods,]
summary(k2o.gbm)
siglimit=1
sigs=(k2o.gbm$padj<=siglimit)
plot(k2o.ge$log2FoldChange[sigs]~k2o.gbm$log2FoldChange[sigs],col=rgb(0,0,0,alpha=0.01),pch=16,cex=0.5,xlab="GBM",ylab="GE",mgp=c(2.3,1,0),xlim=c(-0.35,0.5),ylim=c(-1.5,1.5))
abline(v=0,col="grey70",lty=3)
abline(h=0,col="grey70",lty=3)
lmm=lm(k2o.ge$log2FoldChange~k2o.gbm$log2FoldChange)
abline(lmm,col="red")
summary(lmm) # R2=0.00155
siglimit=0.1
sigs=(k2o.gbm$padj<=siglimit)
points(k2o.ge$log2FoldChange[sigs]~k2o.gbm$log2FoldChange[sigs],col=rgb(0.8,0.2,0.2,alpha=0.5),pch=16,cex=0.5)

# O vs K by origin (Fig 3B)

goods=intersect(row.names(ori.ge),row.names(ori.s))
length(goods)
ori.ge=ori.ge[goods,]
ori.gbm=ori.s[goods,]
summary(ori.gbm)
siglimit=1
sigs=(ori.gbm$padj<=siglimit)
plot(ori.ge$log2FoldChange[sigs]~ori.gbm$log2FoldChange[sigs],col=rgb(0,0,0,alpha=0.01),pch=16,cex=0.5,xlab="GBM",ylab="GE",mgp=c(2.3,1,0),xlim=c(-2,2),ylim=c(-3,3))
abline(v=0,col="grey70",lty=3)
abline(h=0,col="grey70",lty=3)
lmm=lm(ori.ge$log2FoldChange~ori.gbm$log2FoldChange)
abline(lmm,col="red")
summary(lmm) # R2=0.025
siglimit=0.1
sigs=(ori.gbm$padj<=siglimit)
points(ori.ge$log2FoldChange[sigs]~ori.gbm$log2FoldChange[sigs],col=rgb(0.8,0.2,0.2,alpha=0.5),pch=16,cex=0.5)
