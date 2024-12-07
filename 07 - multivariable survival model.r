

library(Publish)

u<-univariateTable(score.cat~SEX+OS_STATUS+OS_MONTHS+CNS_REGION+CANCER_TYPE_DETAILED+EFS_STATUS,data=annot)
res<-summary(u)

write.table(res,file="pbta_publish_results.tsv",sep="\t",row.names=F)


res.cut <- surv_cutpoint(annot, time = "OS_MONTHS", event = "OS_STATUS",
                         variables = c("AGE"))
## age cut
summary(res.cut)

plot(res.cut, "AGE", palette = "Set1")

res.cat <- surv_categorize(res.cut)
head(res.cat)
library(survival)
fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~AGE, data = res.cat)

ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE,palette="Set1",pval = TRUE,
risk.table.y.text=F,ggtheme=theme_classic2(base_size=16),conf.int.style="step")

annot$age.cat<-res.cat$AGE
annot$age.cat<-as.factor(annot$age.cat)
annot$age.cat<-relevel(annot$age.cat,ref="low")

annot$CANCER_TYPE_DETAILED<-gsub("Medulloblastoma, ","",annot$CANCER_TYPE_DETAILED)
annot$CANCER_TYPE_DETAILED<-gsub(" ",".",annot$CANCER_TYPE_DETAILED)
annot%>%dplyr::rename(groups="CANCER_TYPE_DETAILED")->annot

fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~groups, data = annot)

ggsurvplot(fit, data = annot, risk.table = TRUE, conf.int = TRUE,palette="Set1",pval = TRUE,legend = "right",
risk.table.y.text=F,ggtheme=theme_classic2(base_size=12),conf.int.style="step")


##remove wnt subtypes for multivariable model
annot%>%filter(CANCER_TYPE_DETAILED!="Medulloblastoma, WNT-activated")->pheno
pheno$CANCER_TYPE_DETAILED<-gsub("Medulloblastoma, ","",pheno$CANCER_TYPE_DETAILED)
pheno$CANCER_TYPE_DETAILED<-gsub(" ",".",pheno$CANCER_TYPE_DETAILED)
pheno%>%dplyr::rename(groups="CANCER_TYPE_DETAILED")->pheno
pheno%>%dplyr::rename(gender="SEX")->pheno
pheno%>%dplyr::rename(older.15yo="age.cat")->pheno

library(survival)
m<-coxph(formula=Surv(OS_MONTHS, OS_STATUS)~score.cat+older.15yo+gender+groups,data=pheno)
m
summary(m)
+antpbh.cat
library(broom)
library(broom.helpers)
library(GGally)
###mvo$antpbh.cat

test<-cox.zph(m)
test
ggcoxzph(test,font.main = 8,ggtheme = theme_classic2(base_size=14))

ggcoef_model(m,exponentiate=T)+scale_color_brewer(palette="Dark2")+
	theme(text=element_text(size=18),legend.position="bottom")

out<-tidy(m,exponentiate=T, conf.int=T)
write.table(out,file="multivariateOK.tsv",row.names=F,sep="\t")


pheno%>%select(OS_STATUS,OS_MONTHS,score.cat,score,older.15yo,gender,groups)->df


library(rms)
ddist <- datadist(df)
oldoption <- options(datadist='ddist')


f<-cph(formula=Surv(OS_MONTHS, OS_STATUS)~score.cat+older.15yo+gender+groups,x=TRUE,y=TRUE,surv=TRUE,data=df)
surv<-Survival(f)

nomo <- nomogram(f,
lp=TRUE,
fun=list(function(x) surv(12,x),
function(x) surv(24,x)),
funlabel=c("1-Year Survival Prob",
"2-Year Survival Prob"))



plot(nomo)

set.seed(136879)
##m number of patient by groups, u=1 time 1 years
cal<-calibrate(f,B=500,method="boot",cmethod="KM",m=100,u=24)
plot(cal,xlab="Predicted probability at 24 months",ylab="Actual OS proportion" )


nomo <- nomogram(f,lp=TRUE,fun=list(function(x) surv(24,x)), funlabel=c("24-Months OS Prob"))

plot(nomo)
ggplot(Predict(f,risk.score))+theme_minimal(base_size=16)+ggtitle("GSE42568 risk.score prediction")
library(regplot)
regplot(m,points=TRUE, rank="sd",failtime = 24) 

save(m,file="os_model.rda")






