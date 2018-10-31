library(DESeq2)

ll=load("../datasets/GE_data.RData")
names(conds.ge)=c("sample","ind","ori","tra")

mns=apply(counts.ge,1,mean)
min(mns)

dds=DESeqDataSetFromMatrix(counts.ge,
	colData = conds.ge,
	design = ~ori+tra)

vsd.ge=assay(vst(dds))

# outlier detection
# rl=vst(dds)
# library(Biobase)
# e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
# library(arrayQualityMetrics)
# arrayQualityMetrics(e,intgroup=c("ori","tra"),force=T)
# dev.off()
# double-click index.html

dds=DESeq(dds)

ori.ge=results(dds,contrast=c("ori","O","K"))
tra.ge=results(dds,contrast=c("tra","O","K"))

summary(ori.ge)
summary(tra.ge)
save(vsd.ge,conds.ge,ori.ge,tra.ge,file="../datasets/GE_simpleModels.RData")

#------------- split models (paired by genet,except for natives)

counts=counts.ge
conds=conds.ge
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
o2k.ge=results(o2k.dd,contrast=c("tra","K","O"))
summary(o2k.ge)

k2o.dd=DESeqDataSetFromMatrix(k2o.c,
	colData = k2o.meta,
	design = ~ind+tra)
k2o.dd=DESeq(k2o.dd)
k2o.ge=results(k2o.dd,contrast=c("tra","O","K"))
summary(k2o.ge)

natives.dd=DESeqDataSetFromMatrix(natives.c,
	colData = natives.meta,
	design = ~tra)
natives.dd=DESeq(natives.dd)
natives.ge=results(natives.dd,contrast=c("tra","O","K"))
summary(natives.ge)

save(vsd.ge,conds.ge,counts.ge,o2k.ge,k2o.ge,natives.ge,file="../datasets/GE_splitModels.RData")
