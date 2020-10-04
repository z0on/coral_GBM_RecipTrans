library(DESeq2)

#------- assembling genic counts table 

counts=read.table('../datasets/genics.counts',sep="\t")
ll=load("../datasets/coord_tables.RData")
bams=read.table("../datasets/bams")[,1]
bams=gsub("allbams/|\\.trim.+","",bams)
head(counts)
table(genes$start == counts[,2])
counts=counts[,-c(1:3)]
row.names(counts)=genes[,1]
colnames(counts)=bams

#-------------- retaining genes with mean count>=20

mns=rowMeans(counts)
table(mns>=20)
counts=counts[mns>=20,]

#-------------- adding bisulfite amplicons, repeats and intergenics

bs=read.table('../datasets/bisulfite.counts',sep="\t")[,-c(1:3)]
row.names(bs)=paste("bs",c(13:25),sep="")
colnames(bs)=bams
ll=load("../datasets/MBD_intergenics.RData")
ll=load("../datasets/MBD_repeats.RData")
row.names(inter1)=paste("intergen",row.names(inter1),sep=".")
row.names(tes1)=paste("repeat",row.names(tes1),sep=".")
counts=rbind(bs,inter1,tes1,counts)
dim(counts)

#-------------- subsetting for only captured, 2 month time point

sample = colnames(counts)
keep = sample[grep('ub|3m', sample,invert=T)]
counts=counts[,keep]
head(counts)
dim(counts)

sam=sub("_2m","",colnames(counts))
num=sub("[KO]+","",sam)
ori=sub("[KO][0-9]+","",sam)
tra=gsub("^[KO]|[0-9]+$","",sam)
conds=data.frame(cbind(sam,ori,tra))
conds$ind=paste(ori,num,sep="")

dim(counts)
dim(conds)

#---------- Generating variance stabilized data, fitting simple models

dds=DESeqDataSetFromMatrix(counts,
	colData = conds,
	design = ~ori+tra)

vsd=assay(vst(dds))


dds=DESeq(dds)

ori.s=results(dds,contrast=c("ori","O","K"))
tra.s=results(dds,contrast=c("tra","O","K"))

summary(ori.s)
summary(tra.s)
dim(vsd)
save(vsd,conds,dds,ori.s,tra.s,file="../datasets/GBM_simpleModels.RData")

#------------- split models (paired by genet,except for natives)

o2k.c=counts[,conds$ori=="O"]
o2k.meta=conds[conds$ori=="O",]
k2o.c=counts[,conds$ori=="K"]
k2o.meta=conds[conds$ori=="K",]
natives.c=counts[,conds$ori==conds$tra]
natives.meta=conds[conds$ori==conds$tra,]

o2k.dd=DESeqDataSetFromMatrix(o2k.c,
	colData = o2k.meta,
	design = ~ind+tra)
o2k.dd=DESeq(o2k.dd)
o2k.gbm=results(o2k.dd,contrast=c("tra","K","O"))
summary(o2k.gbm)

k2o.dd=DESeqDataSetFromMatrix(k2o.c,
	colData = k2o.meta,
	design = ~ind+tra)
k2o.dd=DESeq(k2o.dd)
k2o.gbm=results(k2o.dd,contrast=c("tra","O","K"))
summary(k2o.gbm)

natives.dd=DESeqDataSetFromMatrix(natives.c,
	colData = natives.meta,
	design = ~tra)
natives.dd=DESeq(natives.dd)
natives.gbm=results(natives.dd,contrast=c("tra","O","K"))
summary(natives.gbm)

save(vsd,conds,counts,o2k.gbm,k2o.gbm,natives.gbm,file="../datasets/GBM_splitModels.RData")

