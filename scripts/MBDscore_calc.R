#------- assembling genic counts table 
library(DESeq2)

counts=read.table('../datasets/genics.counts',sep="\t")
ll=load("../datasets/coord_tables.RData")
bams=read.table("../datasets/bams")[,1]
bams=gsub("allbams/|\\.trim.+","",bams)
head(counts)
table(genes$start == counts[,2])
counts=counts[,-c(1:3)]
row.names(counts)=genes[,1]
colnames(counts)=bams


#-------- selecting samples sequenced both ways (captured and flow-through fractions)

s12=c("KK1_2m","KK1_2ub","KK2_2m","KK2_2ub","KK3_2m","KK3_2ub","KO1_2m","KO1_2ub","KO2_2m","KO2_2ub","KO3_2m","KO3_2ub","OK1_2m", "OK1_2ub","OK2_2m","OK2_2ub","OK3_2m","OK3_2ub","OO1_2m","OO1_2ub","OO2_2m","OO2_2ub","OO3_2m","OO3_2ub")
counts12=counts[,s12]
dim(counts12)

#set up boolean variable for catpured vs flowthrough
sample=colnames(counts12)
captured = sample
captured[grep('m', captured)] <- TRUE
captured[grep('ub', captured)] <- FALSE
sample
captured

colony.id = sub("_.+","",sample)

#build the dataframe with all the varialbes
coldata <- data.frame(sample, captured, colony.id)
coldata

#build DESeq input table
ddsHTSeq<-DESeqDataSetFromMatrix(counts12,
	colData = coldata,
	design = formula(~colony.id + captured))

# run DESeq
dds=DESeq(ddsHTSeq,fitType="local")
mbd = results(dds, contrast=c('captured', 'TRUE', 'FALSE'))
head(mbd) # MBD score is res$log2FoldChange


# plotting Fig. S1 A
hist(mbd$log2FoldChange, breaks = 100)  

# save the result
save(mbd,file='../datasets/capturedVflowthroughResults.Rdata')

#---------------- comparing power of two-way and captured-only GBM change detection

coldata$ori=sub("..$","",coldata$colony.id)
coldata$tra=sub(".$","",coldata$colony.id)
coldata$tra=sub("^.","",coldata$tra)
coldata$ind=paste(coldata$ori,gsub("[KO]","",coldata$colony.id),sep="")

o2k.c=counts12[,coldata$ori=="O"]
o2k.m=coldata[coldata$ori=="O",]

k2o.c=counts12[,coldata$ori=="K", ]
k2o.m=coldata[coldata$ori=="K" ,]

# full model, o2k
o2kf=DESeqDataSetFromMatrix(o2k.c,
	colData = o2k.m,
	design = formula(~(ind+tra)*captured))
o2kf=DESeq(o2kf)
resultsNames(o2kf)
o2kf.res=results(o2kf,name="traO.capturedTRUE")

# captured-only, o2k
o2k.c0=o2k.c[,o2k.m$captured==T]
o2k.m0=o2k.m[o2k.m$captured==T,]
o2k0=DESeqDataSetFromMatrix(o2k.c0,
	colData = o2k.m0,
	design = formula(~ind+tra))
o2k0=DESeq(o2k0)
resultsNames(o2k0)
o2k0.res=results(o2k0,name="tra_O_vs_K")

plot(o2k0.res$log2FoldChange~o2kf.res$log2FoldChange,pch=16,cex=0.5,col=rgb(0,0,0,alpha=0.05),mgp=c(2.3,1,0),xlab="two-way",ylab="captured-only",main="O to K")
abline(0,1,col="red")

# full model, k2o
k2of=DESeqDataSetFromMatrix(k2o.c,
	colData = k2o.m,
	design = formula(~(ind+tra)*captured))
k2of=DESeq(k2of)
resultsNames(k2of)
k2of.res=results(k2of,name="traO.capturedTRUE")

# captured-only, k2o
k2o.c0=k2o.c[,k2o.m$captured==T]
k2o.m0=k2o.m[k2o.m$captured==T,]
k2o0=DESeqDataSetFromMatrix(k2o.c0,
	colData = k2o.m0,
	design = formula(~ind+tra))
k2o0=DESeq(k2o0)
resultsNames(k2o0)
k2o0.res=results(k2o0,name="tra_O_vs_K")

plot(k2o0.res$log2FoldChange~k2of.res$log2FoldChange,pch=16,cex=0.5,col=rgb(0,0,0,alpha=0.05),mgp=c(2.3,1,0),xlab="two-way",ylab="captured-only",main="K to O")
abline(0,1,col="red")

# ----- for shits and giggles: do we have reciprocal change? 
# (note: here we use only N=3 per group instead of N=11 as in main paper)

# full model:
plot(o2kf.res$log2FoldChange~k2of.res$log2FoldChange,pch=16,cex=0.5,col=rgb(0,0,0,alpha=0.05),mgp=c(2.3,1,0),xlab="k2o",ylab="o2k",main="full model",xlim=c(-1,1),ylim=c(-1,1))
summary(lm(o2kf.res$log2FoldChange~k2of.res$log2FoldChange))

# captured-only
plot(o2k0.res$log2FoldChange~k2o0.res$log2FoldChange,pch=16,cex=0.5,col=rgb(0,0,0,alpha=0.05),mgp=c(2.3,1,0),xlab="k2o",ylab="o2k",main="captured-only",xlim=c(-1,1),ylim=c(-1,1))
summary(lm(o2k0.res$log2FoldChange~k2o0.res$log2FoldChange))

	
	

