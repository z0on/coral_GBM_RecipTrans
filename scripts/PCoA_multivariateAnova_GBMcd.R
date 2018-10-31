
#========= PCoA and multivariate anova
#setwd('~/Dropbox/Documents/dixon_2017_RT-GBM/reciprocal_transplant_methylationV2/clean_july12/scripts/')
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

# ------- PCoA colored by GBM class difference

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

# correlation of GBMcd with coral's genotype
summary(lm(class.delta~conds$ind))


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
#save(classes,file="datasets/mbd_classes_FULL_GENES.Rdata")

