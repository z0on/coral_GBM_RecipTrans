# sampling 500 intergenic and 500 repeated regions

tes=read.table("../datasets/intergenics.counts")
bams=read.table("../datasets/bams")[,1]
bams=gsub("allbams/|\\..+","",bams)
row.names(tes)=paste(tes[,1],tes[,2],tes[,3],sep=".")
tes.lng=tes[,3]-tes[,2]
tes=tes[,-c(1:3)]
names(tes)=bams
sds=apply(log(tes+0.1),1,sd)
hist(sds)
tes=tes[sds>quantile(sds,0.05) & sds<quantile(sds,0.95),]
sds=apply(log(tes+0.1),1,sd)
hist(sds)
head(tes)
tes1=tes[sample(1:nrow(tes),500),]
nrow(tes1)
inter=tes
inter1=tes1
save(inter,inter1,file="datasets/MBD_intergenics.RData")

tes=read.table("../datasets/repeats.counts")
bams=read.table("../datasets/bams")[,1]
bams=gsub("allbams/|\\..+","",bams)
row.names(tes)=paste(tes[,1],tes[,2],tes[,3],sep=".")
tes.lng=tes[,3]-tes[,2]
tes=tes[,-c(1:3)]
names(tes)=bams
sds=apply(log(tes+0.1),1,sd)
hist(sds)
tes=tes[sds>quantile(sds,0.05) & sds<quantile(sds,0.95),]
head(tes)
tes1=tes[sample(1:nrow(tes),500),]
nrow(tes1)
save(tes1,tes,file="../datasets/MBD_repeats.RData")
head(tes1)

