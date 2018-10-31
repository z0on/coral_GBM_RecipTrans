# ----- GBM:
ll=load("../datasets/GBM_simpleModels.RData")
traits.gbm=conds
names(traits.gbm)[1]="sample"
vsd.gbm=vsd

# ----- GE:
ll=load("../datasets/GE_simpleModels.RData")
traits.ge=conds.ge
names(traits.ge)[1]="sample"

# ----- GT (identity-by-state from ANGSD):
bams=read.table("../datasets/bams.ibs")[,1] # list of bam files
bams=sub(".bam","",bams)
pop=rep("O",length(bams))
pop[K]="K"
popt=rep("K",length(bams))
popt[K]="O"
ind=sub("[OK]","",bams)
sample=c()
for(i in 1:length(bams)) {
	sample=append(sample,paste(pop[i],pop[i],ind[i],sep=""))
}
for(i in 1:length(bams)) {
	sample=append(sample,paste(pop[i],popt[i],ind[i],sep=""))
}
ma = as.matrix(read.table("../datasets/gbmDAll.ibsMat"))
dimnames(ma)=list(bams,bams)
traits.gt=data.frame("sample"=sample,"ori"=c(pop,pop),"tra"=c(pop,popt))
vsd.gt=cbind(ma,ma)
colnames(vsd.gt)=paste(traits$ori,traits$tra,sep="")


# ========== DAPC, similarity to locals ========== 

dapc.coral=function(vsd,traits){

	require(adegenet)
	
	# making masks for extracting natives
	traits$oritra=as.factor(paste(traits$ori,traits$tra,sep=""))
	vsd.natives=vsd[,traits$ori==traits$tra]
	tr.natives=traits[traits$ori==traits$tra,]
	
	# DAPC analysis: creating discriminant function to tell natives apart
	
	# establishing discriminant function
	dp=dapc(t(vsd.natives),tr.natives$ori,perc.pca=80, n.da=1)
	
	# scoring every sample according to discriminant function
	pred=predict.dapc(dp,newdata=(t(vsd))) 
	
	# plotting "ghost plots"
	
	plot(density(pred$ind.scores),col="white",bty="n",yaxt="n",ylab="",xlab="discriminant function",ylim=c(0,0.8),main="",mgp=c(2.3,1,0))
	polygon(density(pred$ind.scores[traits$oritra=="KK"]),col=rgb(0,0,1,alpha=0.5),border="blue")
	polygon(density(pred$ind.scores[traits$oritra=="OO"]),col=rgb(1,0,0,alpha=0.5),border="red")
	# don't print "transplant ghosts" for GT data :
	if (nrow(vsd)>100){ 
		polygon(density(pred$ind.scores[traits$oritra=="OK"]),col=rgb(252/255,228/255,10/255,alpha=0.5),border="gold")
		polygon(density(pred$ind.scores[traits$oritra=="KO"]),col=rgb(10/255,195/255,252/255,alpha=0.5),border="skyblue")	
	}

	# similarity to locals	

	Kloc=mean(pred$ind.scores[which(traits$ori=="K" & traits$tra=="K")])
	Oloc=mean(pred$ind.scores[which(traits$ori=="O" & traits$tra=="O")])
	
	sim2loc=c()
	for (i in 1:nrow(traits)) {
		if (traits$tra[i]=="K") {
			sim2loc=append(sim2loc,(-1)*abs(Kloc-pred$ind.scores[i]))
		} else {
			sim2loc=append(sim2loc,(-1)*abs(Oloc-pred$ind.scores[i]))
		}
	}
	traits$sim2loc=(sim2loc-mean(sim2loc))/sd(sim2loc)
	return(traits)
}

gbm=dapc.coral(vsd.gbm,traits.gbm)[,c("sample","sim2loc")]
ge=dapc.coral(vsd.ge,traits.ge)[,c("sample","sim2loc")]
gt=dapc.coral(vsd.gt,traits.gt)[,c("sample","sim2loc")]

# ---adding similarity to locals to the main traits dataset

load('../datasets/traits.Rdata')
traits=merge(traits,gbm,by="sample",all.x=TRUE)
names(traits)[ncol(traits)]="gbm.sim2loc"
traits=merge(traits,ge,by="sample",all.x=TRUE)
names(traits)[ncol(traits)]="ge.sim2loc"
traits=merge(traits,gt,by="sample",all.x=TRUE)
names(traits)[ncol(traits)]="gt.sim2loc"
head(traits)

save(traits,file='../datasets/traits_sim2loc.Rdata')

#--------- plotting PCA of fitness traits in transplanted frags, correlation with "similarity to locals"

load('../datasets/traits_sim2loc.Rdata')
head(traits)
# remove zoox and area
trans=traits[traits$ori!=traits$tra,-c(8,9)]
trans=na.omit(trans)
row.names(trans)=trans$sample
fits= trans[,c(5:8)]
ori=trans$ori
envs= trans[,c(9:11)]
names(envs)=c("GBM","GE","GT")

library(vegan)

# computing PCA of transplanted frags, removing effect of origin:
fitpc=rda(fits~ori,scale=T)
axes2plot=c(2,3)
plot(fitpc,choices=axes2plot)
str(fitpc)
pc1expl=100*round(summary(eigenvals(fitpc))[2, axes2plot[1]],2)
pc2expl=100*round(summary(eigenvals(fitpc))[2, axes2plot[2]],2)
grp=rep("red",nrow(fits))
grp[grep("KO",row.names(fits))]="blue"

# flipping axes to match figure in the paper
fitpc$CA$u=fitpc$CA$u*(-1)
fitpc$CA$u[,2]=fitpc$CA$u[,2]*(-1)
fitpc$CA$v=fitpc$CA$v*(-1)
fitpc$CA$v[,2]=fitpc$CA$v[,2]*(-1)

# plotting PCA of transplanted frags, after effect of origin is removed:
plot(fitpc$CA$u,pch=16, cex=1,col=grp,xlim=c(-0.5,0.7),ylim=c(-0.7,0.7),xlab=paste("PC1 ( ",pc1expl,"% )",sep=""),ylab=paste("PC2 ( ",pc2expl,"% )",sep=""),mgp=c(2.3,1,0),main="Fitness proxies")
arrows(rep(0,4),rep(0,4), fitpc$CA$v[,1]/1.5, fitpc $CA$v[,2]/1.5,length=0.05)
text(fitpc$CA$v[,1]/1.15, fitpc$CA$v[,2]/1.15+c(-0.035,0,0,0),labels=row.names(fitpc$CA$v)
,cex=0.8,col=c(rep(1,4),rep(4,3)))
legend("topleft",pch=16,col= c("blue","red"),legend=c("K","O"),title="origin",cex=0.8)

# fittning similarities-to-locals to the fitness ordination:
ef=envfit(fitpc,choices=c(2:3),envs,permutations = 9999)
ef
plot(ef,col="purple")

# correlation of GBM-similarity-to-locals with fitness PC1:
summary(lm(envs$GBM~fitpc$CA$u[,1]))

# correlation of GE-similarity-to-locals with fitness PC1:
summary(lm(envs$GE~fitpc$CA$u[,1]))

# correlation of GT-similarity-to-locals with fitness PC1:
summary(lm(envs$GT~fitpc$CA$u[,1]))
