#-------- aggregating counts per gene

counts=read.table("../datasets/introns.counts")
# un-remark next line to analyze exons
# counts=read.table("../datasets/exons.counts")
head(counts)
bams=read.table('../datasets/bams')[,1]
bams=gsub("allbams/|\\.trim.+$","",bams)
ll=load("../datasets/coord_tables.RData")

coords=introns
# un-remark next line to analyze exons
# coords=exons
table(counts[,1]==coords[,2])
counts=counts[,-c(1:3)]
dim(counts)
length(bams)
names(counts)=bams

counts.aggr=c()
pb=txtProgressBar(0,length(unique(coords$gene)))
j=0
for (i in unique(coords$gene)) {
	j=j+1
	setTxtProgressBar(pb,j)
	goods=which(coords $gene==i)
	gc=apply(counts[goods,],2,sum)
	counts.aggr=data.frame(rbind(counts.aggr,gc))
}
row.names(counts.aggr)=unique(coords$gene)
nrow(counts.aggr)
head(counts.aggr)
names(counts.aggr)=bams

#-------------- DESeq split models

# assembling captured-only dataset
library(DESeq2)

counts=counts.aggr
mns=rowMeans(counts)
table(mns>=5)
counts=counts[mns>=5,]

dim(counts)
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
#dim(vsd)

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

# natives.dd=DESeqDataSetFromMatrix(natives.c,
	# colData = natives.meta,
	# design = ~tra)
# natives.dd=DESeq(natives.dd)
# natives.gbm=results(natives.dd,contrast=c("tra","O","K"))
# summary(natives.gbm)

save(conds,counts,o2k.gbm,k2o.gbm,file="../datasets/deseq_exons.RData")
