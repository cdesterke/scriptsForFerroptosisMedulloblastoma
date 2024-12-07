load("datamb.rda")
load("phenomb.rda")

load("ferrdbhs.rda")

library(preprocessCore)


data<-sel


colramp = colorRampPalette(c(3,"white",2))(50)
edata = data[rowMeans(data) > 0.5, ]
edata<-log(edata+1,10)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,2))
	for(i in 1:257){lines(density(edata[,i]),lwd=1,col=colramp[i])}


norm_edata = normalize.quantiles(as.matrix(edata))
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,1))
for(i in 1:257){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}
row.names(norm_edata)<-row.names(edata)
colnames(norm_edata)<-colnames(edata)





ferrdbhs%>%select(hs.gene)%>%distinct()%>%pull(hs.gene)->vector
sel<-norm_edata[row.names(edata)%in%vector,]


df<-read.table("OS_loopcolcox.tsv",h=T,sep="\t")
 
library(dplyr)

ferrdbhs%>%right_join(df,by=c("hs.gene"="identifiers"))%>%arrange(desc(NLP))%>%filter(pvalues<0.05)->results

table(results$class,results$prognosis)

write.table(results,file="oswithferroptosis.tsv",row.names=F)

library(transpipe14)

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

annot$kmeans<-kmeans_result$cluster

annot$kmeans<-as.factor(annot$kmeans)




pcatrans(sel,annot,group="kmeans",pal="Set1",alpha=1,names=F)


library(sva)
batch = annot$kmeans
mod = model.matrix(~1, data=pheno)

combat_edata = ComBat(dat=sel,  batch=batch, mod=mod,par.prior=TRUE, prior.plots=TRUE)

pcatrans(combat_edata,annot,group="TUMOR_TYPE",pal="Set1",alpha=1,names=F)


data<-combat_edata
save(data,file="ferroptosisSuperdataCombat.rda")






annot%>%select(CANCER_TYPE_DETAILED,TUMOR_TYPE,CNS_REGION)->pheno
row.names(pheno)<-colnames(data)

bestheat(data,pheno,font=10,rownames=F,scale="row")


matrix<-t(data)

library(loopcolcox)
load("phenomb.rda")
annot$OS_STATUS<-as.numeric(annot$OS_STATUS)
df<-coxbycol(annot$OS_MONTHS ,annot$OS_STATUS ,matrix)
head(df)



df%>%filter(coef.beta>-1)->df2
coxvolcano(df,font.size=18)


library(patchwork)
## nb : number of covariates to put on the graph
p1<-plotbeta(df,nb=30,title="",size=16)
p2<-plotnlphr(df,nb=30,title="",size=16)
p1+p2

write.table(df,file="OS_loopcolcox.tsv",row.names=F,sep="\t")

