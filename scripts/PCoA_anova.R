
#========= PCoA and multivariate anova
setwd('~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/clean_july12/scripts/')
#---------  GBM

ll=load("../datasets/GBM_splitmodels.RData")
colnames(vsd)=sub("_2m","",colnames(vsd))
library(vegan)

# multivariate anova:

ad=adonis(t(vsd)~ori+ind+tra,conds,method="manhattan")
ad 
          # Df  SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# ori        1   44768030 44768030  3.9300 0.02957  0.003 ** 
# ind       20 1203049617 60152481  5.2805 0.79463  0.001 ***
# tra        1   26934455 26934455  2.3644 0.01779  0.029 *  
# Residuals 21  239220955 11391474         0.15801           

# plotting pie-chart of variance explained
labs=c("origin","colony.id","transplant","residuals")
cols=c("skyblue","green2","coral","grey80")
pie(ad$aov.tab$R2[1:4],labels=labs,col=cols,main="GBM")

# PCoA

fitpc=capscale(dist(t(vsd),method="manhattan")~1)
pc1expl=100*round(summary(eigenvals(fitpc))[2,1],2)
pc2expl=100*round(summary(eigenvals(fitpc))[2,2],2)
grp=rep("red",ncol(vsd))
grp[grep("OK",colnames(vsd))]="gold"
grp[grep("KK",colnames(vsd))]="blue"
grp[grep("KO",colnames(vsd))]="skyblue"
conds$oritra=paste(conds$ori,conds$tra,sep="")
plot(fitpc$CA$u,pch=16, col="white",xlab=paste("PC1 ( ",pc1expl,"% )",sep=""),ylab=paste("PC2 ( ",pc2expl,"% )",sep=""),mgp=c(2.3,1,0),main="GBM")
ordiellipse(fitpc$CA$u,conds$oritra,draw="polygon",col=c("blue","skyblue","gold","red"))
ordispider(fitpc$CA$u,conds$ind,col="grey50")
points(fitpc$CA$u,pch=16, col=grp,xlab=paste("PC1 ( ",pc1expl,"% )",sep=""))

# % variation explained
plot(summary(eigenvals(fitpc))[2,]*100,cex=0.8, ylab="% var explained",xlab="PC",mgp=c(2.1,1,0))

# PCoA colored by GBM class delta

ll=load("../datasets/MBDclasses_q30.RData")
vsd1=data.frame(vsd)
vsd1$gene=row.names(vsd)
vsd1=merge(vsd1,classes,by="gene")
head(vsd1)

library(Rmisc)
class.delta=c()
for (s in colnames(vsd)){
	see=summarySE(vsd1,measurevar=s,groupvars="class")
	class.delta=append(class.delta,see[2,3]-see[1,3])
}

rbPal <- colorRampPalette(c('white','black'))
Cols <- rbPal(10)[as.numeric(cut(class.delta[order(class.delta)],breaks = 10))]
reo=fitpc$CA$u[order(class.delta),]
plot(reo,col="grey50",xlab=paste("PC1 ( ",pc1expl,"% )",sep=""),ylab=paste("PC2 ( ",pc2expl,"% )",sep=""),mgp=c(2.3,1,0),main="GBM")
points(reo,pch=16, col=Cols)
ge.cd=data.frame("GEcd"=class.delta)
ge.cd$sample=row.names(fitpc$CA$u)
ef=envfit(fitpc,ge.cd,permu=9999)
ef # p<1e-5, R2=0.36
plot(ef,cex=0.001,col="red")

summary(lm(class.delta~conds$ind))

ges=data.frame("ge.pc1"=fitpc$CA$u[,1],"GEcd"=class.delta)
ges$Colony.ID=row.names(fitpc$CA$u)
str(ges)

x=merge(traits,ges,by="Colony.ID",all.x=T)
nrow(x)
nrow(traits)

#---------  GE

ll=load("../datasets/GE_splitModels.RData")
vsd=vsd.ge
conds=conds.ge
colnames(vsd)=sub("_2m","",colnames(vsd))

# gene expression vs MBD score:
ll=load('../datasets/mbd_classes_FULL_GENES.Rdata')
head(classes)
goods=which(row.names(vsd) %in% classes$gene)
vsd1=vsd[goods,]
gene=row.names(vsd1)
ge=apply(vsd1,1,mean)
ge=data.frame(cbind(gene,ge))
ge$ge=as.numeric(as.character(ge$ge))
ge2mbd=merge(classes,ge,by="gene")

plot(ge~mbd.score,ge2mbd,pch=16, col=rgb(0,0,0,alpha=0.1),cex=0.2,mgp=c(2.1,1,0),xlim=c(-5.5,5))
abline(lm(ge~mbd.score,ge2mbd),col="red")
summary(lm(ge~mbd.score,ge2mbd)) # R2 = 0.0044,  p < 2e-16

# multivariate anova:

ad=adonis(t(vsd)~ori+ind+tra,conds,method="manhattan")
ad 
          # Df  SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# ori        1  111027912 111027912  2.8536 0.03640  0.001 ***
# ind       28 1943236567  69401306  1.7837 0.63705  0.001 ***
# tra        1  101183324 101183324  2.6005 0.03317  0.001 ***
# Residuals 23  894898303  38908622         0.29338           

# plotting pie-chart of variance explained
labs=c("origin","colony.id","transplant","residuals")
cols=c("skyblue","green2","coral","grey80")
pie(ad$aov.tab$R2[1:4],labels=labs,col=cols,main="GBM")

# PCoA

fitpc=capscale(dist(t(vsd),method="manhattan")~1)
pc1expl=100*round(summary(eigenvals(fitpc))[2,1],2)
pc2expl=100*round(summary(eigenvals(fitpc))[2,2],2)
grp=rep("red",ncol(vsd))
grp[grep("OK",colnames(vsd))]="gold"
grp[grep("KK",colnames(vsd))]="blue"
grp[grep("KO",colnames(vsd))]="skyblue"
conds$oritra=paste(conds$ori,conds$tra,sep="")
plot(fitpc$CA$u,pch=16, col="white",xlab=paste("PC1 ( ",pc1expl,"% )",sep=""),ylab=paste("PC2 ( ",pc2expl,"% )",sep=""),mgp=c(2.3,1,0),main="GE")
ordiellipse(fitpc$CA$u,conds$oritra,draw="polygon",col=c("blue","skyblue","gold","red"))
ordispider(fitpc$CA$u,conds$ind,col="grey50")
points(fitpc$CA$u,pch=16, col=grp,xlab=paste("PC1 ( ",pc1expl,"% )",sep=""))

# % variation explained
plot(summary(eigenvals(fitpc))[2,]*100,cex=0.8, ylab="% var explained",xlab="PC",mgp=c(2.1,1,0))

ll=load('~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/datasets/GE_3mo.RData')
allcs=capscale(dist(t(vsd),method="manhattan")~1)
ad=adonis(t(vsd)~origin+colony.id+transplant,conditions,method="manhattan")
           # Df  SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# origin      1  111031433 111031433  2.8535 0.03640  0.001 ***
# colony.id  28 1943255943  69401998  1.7837 0.63705  0.001 ***
# transplant  1  101186089 101186089  2.6005 0.03317  0.001 ***
# Residuals  23  894927964  38909911         0.29338           
# Total      53 3050401430                   1.00000           

labs=c("origin","colony.id","transplant","residuals")
cols=c("skyblue","green2","coral","grey80")
pie(ad$aov.tab$R2[1:4],labels=labs,col=cols,main="GE")

#---------- CpG-OE vs MBD score

ll=load('../datasets/mbd_classes_FULL_GENES.Rdata')
head(classes)

hist(classes$mbd.score,breaks=40,yaxt="n",ylab="",xlim=c(-5.5,5.5))
plot(density(classes$mbd.score),yaxt="n",ylab="",xlim=c(-6,6),bty="n")
abline(v=max(subset(classes,mbd.class==1)$mbd.score),col="red")
load("../datasets/cpgoe.Rdata")

plot(res2[,'cpgOE'] ~ res2[,'log2FoldChange'], xlab ="MBD-score", ylab =expression("CpG"["o/e"]), col=rgb(0,0,0,alpha=0.01),pch=16,cex=0.5,xlim=c(-7,7), ylim = c(0, 1.5),mgp=c(2.3,1,0),bty="n")
abline(v=max(subset(classes,mbd.class==1)$mbd.score)+0.5,col="red")
mbdcut=max(subset(classes,mbd.class==1)$mbd.score)+0.5
classes$mbd.class2=1
classes$mbd.class2[classes$mbd.score>mbdcut]=2
save(classes,file="datasets/mbd_classes_FULL_GENES.Rdata")


###--------- G B M PCA

setwd("~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2")
# ll=load("datasets/GBM_htseq_may26_deseq.RData")
# htseq.genes=row.names(o2k.gbm)

ll=load("datasets/GBM_q30_deseq_min20.RData")
# load("~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/datasets/GBM_q30_deseq_exons.RData")
# load("~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/datasets/GBM_q30_deseq_introns.RData")

#ll=load("datasets/mbd_3mo_FULLGENES.RData")
# mns=apply(vsd,1,mean)
# lim=quantile(mns,0)
# vsd=vsd[mns>lim,]
# sds=apply(vsd,1,sd)
# limsd=quantile(sds,0.75)
# vsd=vsd[sds>limsd,]

vsd=vsd[grep("repeat|intergen|bs",row.names(vsd),invert=T),]
dim(vsd)

library(vegan)
adonis(dist(t(vsd),method="manhattan")~ori+ind+tra,conditions)
          # Df  SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# ori        1   41747850 41747850  3.9249 0.02937  0.003 ** 
# ind       20 1130982324 56549116  5.3165 0.79552  0.001 ***
# tra        1   25586727 25586727  2.4055 0.01800  0.021 *  
# Residuals 21  223367455 10636545         0.15711           
# Total     43 1421684357                  1.00000           



#fitpc=capscale(dist(t(vsd),method="manhattan")~colony.id,conditions)
fitpc=capscale(dist(t(vsd),method="manhattan")~1)
pc1expl=100*round(summary(eigenvals(fitpc))$importance[2,1],2)
pc2expl=100*round(summary(eigenvals(fitpc))$importance[2,2],2)
grp=rep("red",ncol(vsd))
grp[grep("OK",colnames(vsd))]="gold"
grp[grep("KK",colnames(vsd))]="blue"
grp[grep("KO",colnames(vsd))]="skyblue"
#fitpc$CA$u=fitpc$CA$u*(-1)
#fitpc$CA$v=fitpc$CA$v*(-1)
conditions=conds
conditions$oritra=paste(conditions$ori,conditions$tra,sep="")
plot(fitpc$CA$u,pch=16, col="white",xlab=paste("PC1 ( ",pc1expl,"% )",sep=""),ylab=paste("PC2 ( ",pc2expl,"% )",sep=""),mgp=c(2.3,1,0),main="GBM")
ordiellipse(fitpc$CA$u,conditions$oritra,draw="polygon",col=c("blue","skyblue","gold","red"))
ordispider(fitpc$CA$u,conditions$ind,col="grey50")
points(fitpc$CA$u,pch=16, col=grp,xlab=paste("PC1 ( ",pc1expl,"% )",sep=""))
varexpl=100*round(summary(eigenvals(fitpc))$importance[2,],3)
#plot(varexpl[1:10],cex=0.8,mgp=c(2.3,1,0),main="GBM",xlab="PC",ylab="% variance")
lines(varexpl)
#str(fitpc)

varexpl=100*round(summary(eigenvals(fitpc))$importance[2,],3)
plot(varexpl[1:10],cex=0.8,mgp=c(2.3,1,0),main="GBM",xlab="PC",ylab="% variance")
lines(varexpl)

#------ GBMcd
ll=load('~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/datasets/capturedVflowthroughResults_may24_q30.Rdata')
classes=data.frame(cbind("gene"=row.names(res)))
classes$mbd.score=res$log2FoldChange
head(classes)
hist(classes$mbd.score,breaks=50)

# printing mean GBM vs mbdscore
ll=load("~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/datasets/GBM_q30_deseq_min20.RData")
dset=o2k.gbm;YL=c(-0.6,0.6)
gg=grep("repeat\\.|interg",row.names(dset),invert=T)
bmean=data.frame(cbind("gene"=row.names(dset[gg,]),"baseMean"=log(dset[gg,]$baseMean,10),"padj"=dset[gg,]$padj))
bmean$baseMean =as.numeric(as.character(bmean$baseMean))
bmean$padj=as.numeric(as.character(bmean$padj))
x=merge(classes,bmean,by="gene")
plot(baseMean~mbd.score,x,col=rgb(0,0,0,alpha=0.01),pch=16,cex=0.5,mgp=c(2.3,1,0),xlab="MBD score",ylab="GBM baseMean")
abline(v=-1,col="red")
classes$class=1
classes$class[classes$mbd.score>-1]=2
head(classes)
save(classes,file="datasets/MBDclasses_q30.RData")

vsd1=data.frame(vsd)
vsd1$gene=row.names(vsd)
vsd1=merge(vsd1,classes,by="gene")
head(vsd1)
s=colnames(vsd)[2]
head(vsd1)
plot(OO8_2m~mbd.score,vsd1,pch=16,cex=0.3,col=rgb(0,0,0,alpha=0.01))

install.packages("Rmisc")
library(Rmisc)
class.delta=c()
for (s in colnames(vsd)){
	see=summarySE(vsd1,measurevar=s,groupvars="class")
	class.delta=append(class.delta,see[2,3]-see[1,3])
}
hi=colnames(vsd)[class.delta==max(class.delta)]
lo=colnames(vsd)[class.delta==min(class.delta)]

plot(density(na.omit(vsd1[vsd1$class==1,hi])),col="blue",xlim=c(2,11),bty="n",yaxt="n",ylab="",xlab="",mgp=c(2.3,1,0),main="")
lines(density(na.omit(vsd1[vsd1$class==2,hi])),col="red")
abline(v=mean(na.omit(vsd1[vsd1$class==1,hi])),lty=3, col="blue")
abline(v=mean(na.omit(vsd1[vsd1$class==2,hi])),lty=3, col="red")

plot(density(na.omit(vsd1[vsd1$class==1,lo])),col="blue",xlim=c(2,11),bty="n",yaxt="n",ylab="",xlab="",mgp=c(2.3,1,0),main="")
lines(density(na.omit(vsd1[vsd1$class==2,lo])),col="red")
abline(v=mean(na.omit(vsd1[vsd1$class==1,lo])),lty=3, col="blue")
abline(v=mean(na.omit(vsd1[vsd1$class==2,lo])),lty=3, col="red")

plot(class.delta~fitpc$CA$u[,1],cex=0.8,mgp=c(2.3,1,0),main="",xlab="GBM PC1",ylab="GBM class delta")
ll=lm(class.delta~fitpc$CA$u[,1])
summary(ll) # GBM: R2 0.97 ; GE: R2 0.23 ; p = 1.5e-4
abline(ll,col="red")
adonis(dist(t(vsd),method="manhattan")~class.delta)
            # Df  SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# class.delta  1  498911227 498911227  20.643 0.32954  0.001 ***
# Residuals   42 1015062002  24168143         0.67046           
# Total       43 1513973229                   1.00000           

gbms=data.frame("gbm.pc1"=fitpc$CA$u[,1],"GBMcd"=class.delta)
gbms$Colony.ID=sub("_2m","",row.names(fitpc$CA$u))
str(gbms)

# x[,"mbd.pco1"]=NULL
# x[,"ge.pco1"]=NULL
# xx=merge(x,gbms,by="Colony.ID",all.x=T)
# names(xx)
# traits=xxx
# save(traits,file="datasets/traits_april9.RData")

#------ PCA colored by MBD class delta
rbPal <- colorRampPalette(c('white','black'))
Cols <- rbPal(10)[as.numeric(cut(class.delta[order(class.delta)],breaks = 10))]
reo=fitpc$CA$u[order(class.delta),]
plot(reo,col="grey25",xlab=paste("PC1 ( ",pc1expl,"% )",sep=""),ylab=paste("PC2 ( ",pc2expl,"% )",sep=""),mgp=c(2.3,1,0),main="GBM")
points(reo,pch=16, col=Cols,)
gbm.cd=data.frame("GBMcd"=class.delta)
gbm.cd$sample=sub("_2m","",row.names(fitpc$CA$u))
ef=envfit(fitpc,gbm.cd,permu=9999)
ef
plot(ef,cex=0.001,col="red")

class.delta=merge(ge.cd,gbm.cd,by="sample")

plot(GEcd~GBMcd,class.delta,mgp=c(2.3,1,0))

save(class.delta,file="datasets/GBM_class_delta.RData")

dim(o2k.gbm)[1]-1013
#----------