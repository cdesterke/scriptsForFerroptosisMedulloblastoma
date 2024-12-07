pheno<-read.table("E-MTAB-10767.sdrf.txt",h=T,sep="\t",na.strings="NA",row.names=1)
library(dplyr)
pheno%>%select(Factor.Value.os.status.,Factor.Value.pfs.status.,Characteristics.subgroup.,Characteristics.sex.)->annot
annot<-pheno
annot$Factor.Value.os.status. <- as.factor (annot$Factor.Value.os.status.)
annot$Factor.Value.pfs.status. <- as.factor (annot$Factor.Value.pfs.status.)
colnames(annot)<-c("os_status","pfs_status","subgroups","gender")
load("ferroptosisSuperdataCombat.rda")
cross<-read.csv("cross_sig.csv",h=T)
library(dplyr)

trans<-as.data.frame(t(data))
all(row.names(trans)==row.names(annot))


cross%>%select(identifiers)%>%distinct()%>%pull(identifiers)->vector

sel<-trans[,colnames(trans)%in%vector]
dim(sel)

sel%>%mutate(score=(USP11*-2.18002474518665)+(MCF2L*-1.04530832774745)+(SOX15*-0.992224884476285)+(BEX1*-0.907587628768039)+
(CDO1*-0.871822526242265)+(TRPV1*-1.46645982768831)+(COQ10A*-1.68798125350889)+(TFAP2C*0.967832606083715)+
(CAV1*0.722841430780254)+(SNX5*1.35635732231724)+(PARP6*-1.36262350508992)+(SQOR*1.33174127872047)+(SLC39A14*1.16349465313119)+
(MAP1LC3A*-0.860689113835478)+(CEMIP*-0.59005755740291)+(TXN*0.94085270257141)+(PCDHB14*-0.811693695123737)+(TFR2*-0.798241687760241)+
(USF2*-1.47385325409416)+(ATF4*0.956404022419562)+(PHGDH*0.483156894698922)+(PPARA*-0.923685758091893)+(G3BP1*1.17800387032181)+
(PRR5*-0.539547152725295)+(CARS1*1.17046928959313)+(PRDX4*0.703561532331052)+(CCT3*1.36836908075737)+(TSC1*-0.710935716498839)+
(CBS*0.553452384967907)+(CLOCK*-1.20946536046247)+(VAMP2*-0.757846020766605)+(IL6*0.577234292149136)+(TIMM9*0.736046933599977)+
(TCF4*-0.754520322329335)+(STC1*0.589174240198105)+(ADAM23*-0.615428764368122)+(HPX*-0.716302803991269)+(KIF20A*0.717871874916388)+
(FAM98A*0.999125573497597)+(SMG9*-1.15882491323182)+(CAPRIN2*-0.852103896784302)+(MTDH*0.90965683321584)+(FXR1*0.978687677829785)+
(AHCY*0.787735924986664)+(CUL9*-0.710689105355085))->sel


annot$score<-sel$score
sel$score<-NULL
library(transpipe15)

library(survminer)
annot%>%select(Characteristics.os.time.,Characteristics.os.status.,Characteristics.subgroup.,Characteristics.sex.,score.cat,score)->pheno
colnames(pheno)<-c("os.time","os.status","groups","gender","score.cat","score")
pheno$os.status<-as.factor(pheno$os.status)
trans<-as.data.frame(t(sel))
bestheat(trans,pheno,rownames=T,font=8)
pcatrans(trans,pheno,"score.cat",alpha=1)

write.table(pheno,file="annotation_score.tsv",row.names=T,sep="\t")

res.cut <- surv_cutpoint(annot, time = "Characteristics.os.time.", event = "Characteristics.os.status.",
                         variables = c("score"))

summary(res.cut)

plot(res.cut, "score", palette = "Set1")

res.cat <- surv_categorize(res.cut)
head(res.cat)
library(survival)
fit <- survfit(Surv(Characteristics.os.time., Characteristics.os.status.) ~score, data = res.cat)

ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE,palette="Set1",pval = TRUE,
risk.table.y.text=F,ggtheme=theme_classic2(base_size=16),conf.int.style="step")

annot$score.cat<-res.cat$score
annot$score.cat<-as.factor(annot$score.cat)
annot$score.cat<-relevel(annot$score.cat,ref="low")
 

library(Publish)

u<-univariateTable(score.cat~Characteristics.developmental.stage.+ Characteristics.sex.+
Characteristics.subgroup.+Characteristics.os.status.+Characteristics.os.time.,data=annot)
res<-summary(u)
res
write.table(res,file="validation_publish_results.tsv",sep="\t",row.names=F)

