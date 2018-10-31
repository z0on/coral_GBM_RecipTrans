#--------- plotting PCA of fitness traits
library(vegan)

load('../datasets/traits.Rdata')

names(traits)
fits=traits[,c(5:7,10)]
row.names(fits)=traits$sample
goods=which(!is.na(apply(fits,1,mean)))
fits=fits[goods,]

fitpc=rda(fits,scale=T)
axes2plot=c(1,2)

pc1expl=100*round(summary(eigenvals(fitpc))[2, axes2plot[1]],2)
pc2expl=100*round(summary(eigenvals(fitpc))[2, axes2plot[2]],2)
grp=rep("red",nrow(fits))
grp[grep("OK",row.names(fits))]="gold"
grp[grep("KK",row.names(fits))]="blue"
grp[grep("KO",row.names(fits))]="skyblue"
fitpc$CA$u=fitpc$CA$u*(-1)
fitpc$CA$v=fitpc$CA$v*(-1)

plot(fitpc$CA$u,pch=16, cex=1,col=grp,xlim=c(-0.4,0.6),ylim=c(-0.7,0.7),xlab=paste("PC1 ( ",pc1expl,"% )",sep=""),ylab=paste("PC2 ( ",pc2expl,"% )",sep=""),mgp=c(2.3,1,0),main="Fitness proxies")
arrows(rep(0,4),rep(0,4), fitpc$CA$v[,1]/1.5, fitpc $CA$v[,2]/1.5,length=0.05)
text(fitpc$CA$v[,1]/1.15, fitpc $CA$v[,2]/1.15+c(-0.035,0,0,0),labels=row.names(fitpc$CA$v)
,cex=0.8,col=c(rep(1,4),rep(4,3)))
legend("topleft",pch=16,col= c("blue","skyblue","gold","red"),legend=c("KK","KO","OK","OO"),cex=0.8)

ordiellipse(fitpc$CA$u,grp,draw="polygon",col=c("blue","gold","red","skyblue"))

#--------- plotting PCA of fitness traits in transplanted frags, correlation with "similaroty to locals"

load('../datasets/traits.Rdata')

fits=traits[,c(5:7,10)]
row.names(fits)=traits$Colony.ID
envs=traits[,c(17,45:46)] # GE
#envs=traits[,c(16,47:48)] # GBM



# ori=traits$ori[goods]
# fitpc=rda(fits~ori,scale=T)
#plot(fitpc,choices=c(2,3))
ef=envfit(fitpc,choices=c(2:5),envs)
ef
# GE:
            # PC1      PC2     r2 Pr(>r)    
# ge.ld1 -0.80305 -0.59591 0.4955  0.001 ***
# ge.pc1  0.98802 -0.15430 0.0719  0.238    
# GEcd    0.98365  0.18007 0.2087  0.010 ** 

# GE after ori fitting:
# ge.ld1 -0.84499  0.18832  0.46895 -0.17498 0.3002  0.021 *
# ge.pc1  0.67096  0.49736  0.26176 -0.48365 0.1189  0.297  
# GEcd    0.91095 -0.15867 -0.08417 -0.37137 0.2332  0.053 .


# GBM:
# gbm.ld1 -0.68328 -0.73015 0.2575  0.013 *
# gbm.pc1 -0.58161  0.81346 0.0629  0.378  
# GBMcd   -0.46164  0.88707 0.0594  0.406  

# GBM, no origin effect
# gbm.ld1 -0.887790  0.138430  0.145919  0.413970 0.0573  0.779
# gbm.pc1 -0.581310 -0.737300  0.050789  0.340420 0.0769  0.677
# GBMcd   -0.492570 -0.818510  0.122036  0.269310 0.0679  0.724


plot(ef,col="purple",cex=0.7)




plot(gbm.ld1~GAIN,traits)
plot(ge.ld1~GAIN,traits)
plot(gbm.ld1~fitpc1,traits)
plot(ge.ld1~fitpc1,traits)

