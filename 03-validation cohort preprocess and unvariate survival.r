pheno<-read.table("E-MTAB-10767.sdrf.txt",h=T,sep="\t",na.strings="NA",row.names=1)

data<-read.table("tpm.genes.reorder.txt",h=T,sep="\t")

bm<-read.table("biomart112HG38clean.tsv",h=T,sep="\t")


library(dplyr)

bm%>%inner_join(data,by=c("ensembl"=="X"))->data

data$ensembl<-NULL

library(transpipe15)

ok<-filtermatrix(data)

data<-log2(ok+1)

library(preprocessCore)

edata = data[rowMeans(data) > 3, ]

colramp = colorRampPalette(c(3,"white",2))(50)
edata = data[rowMeans(data) > 2, ]
plot(density(data[,1]),col=colramp[1],lwd=3,ylim=c(0,.25))
	for(i in 1:331){lines(density(edata[,i]),lwd=1,col=colramp[i])}






plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 1:331){lines(density(edata[,i]),lwd=3,col=colramp[i])}

norm_edata = normalize.quantiles(as.matrix(edata))
rownames(norm_edata)<-rownames(edata)
colnames(norm_edata)<-colnames(edata)

plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 1:331){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}

data<-norm_edata



all(row.names(pheno)==colnames(data))

data<-data[,row.names(pheno)]

save(data,file="quantile.rda")

load("ferrdbhs.rda")

ferrdbhs%>%select(hs.gene)%>%distinct()%>%pull(hs.gene)->vector

sel<-data[row.names(data)%in%vector,]

pcatrans(data,pheno, group="Characteristics.subgroup.")


## kmeans
library(stats)
wss <- function(k) {
  kmeans(t(sel), centers = k, nstart = 25)$tot.withinss
}

# Compute WSS for k = 1 to 10
k.values <- 1:10
wss.values <- sapply(k.values, wss)

# Plot to find the "elbow"
plot(k.values, wss.values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster Sum of Squares (WSS)")

kmeans_result <- kmeans(t(sel), centers = 2, nstart = 25)

pheno$kmeans<-kmeans_result$cluster

pheno$kmeans<-as.factor(pheno$kmeans)




pcatrans(sel,pheno,group="Characteristics.subgroup.",pal="Set1",alpha=1,names=F)


library(sva)
batch = pheno$kmeans
mod = model.matrix(~1, data=pheno)

combat_edata = ComBat(dat=sel,  batch=batch, mod=mod,par.prior=TRUE, prior.plots=TRUE)
pcatrans(combat_edata,pheno,group="Characteristics.subgroup.",pal="Paired",alpha=1,names=F)
pcatrans(combat_edata,pheno,group="kmeans",pal="Set1",alpha=1,names=F)


data<-combat_edata
save(data,file="ferroptosisSuperdataCombat.rda")

pheno%>%select(Factor.Value.os.status.,Factor.Value.pfs.status.,Characteristics.subgroup.,Characteristics.sex.)->annot
annot$Factor.Value.os.status. <- as.factor (annot$Factor.Value.os.status.)
annot$Factor.Value.pfs.status. <- as.factor (annot$Factor.Value.pfs.status.)
colnames(annot)<-c("os_status","pfs_status","subgroups","gender")

bestheat(data,annot,font=10,rownames=F)

matrix<-t(data)

library(loopcolcox)
annot$os_time<-pheno$Factor.Value.os.time.
annot$pfs_time<-pheno$Factor.Value.pfs.time.

annot$os_status<-as.numeric(annot$os_status)
df<-coxbycol(annot$os_time ,annot$os_status ,matrix)
head(df)



df%>%filter(coef.beta>-1)->df2
coxvolcano(df,font.size=18)


library(patchwork)
## nb : number of covariates to put on the graph
p1<-plotbeta(df,nb=30,title="",size=16)
p2<-plotnlphr(df,nb=30,title="",size=16)
p1+p2

write.table(df,file="OS_loopcolcox.tsv",row.names=F,sep="\t")



