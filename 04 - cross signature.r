list.files()

train<-read.table("OS_loopcolcox_training.tsv",h=T,sep="\t")
train%>%filter(significance=="YES")->train

val<-read.table("OS_loopcolcox_validation.tsv",h=T,sep="\t")
val%>%filter(significance=="YES")->val

library(dplyr)
train%>%inner_join(val,by="identifiers")->cross

cross$conca<-paste(cross$prognosis.x,cross$prognosis.y,sep="_")

cross%>%mutate(nlp=(NLP.x+NLP.y)/2)->cross

library(ggplot2)
library(ggrepel)
cross%>%filter(conca!="favorable_unfavorable")->cross
   
ggplot(cross,aes(coef.beta.x,coef.beta.y))+
	geom_smooth(method = "lm", color = "green", se = F)+
	geom_point(aes(size=nlp,color=conca))+
	geom_text_repel(aes(label = identifiers),size=4)+theme_classic(base_size=18)+theme(legend.position = "none")+
xlab("beta coefficients training cohort")+ylab("beta coefficients validation cohort")+
	scale_color_manual(values = c("steelblue", "tomato")) 

cross%>%mutate(beta.coef=(coef.beta.x+coef.beta.y)/2)->cross


id<-cross$identifiers
beta<-cross$beta.coef

equation<-paste(id,beta,sep="*",collapse=")+(")
equation<-paste("(",equation,")",sep="")
equation


(USP11*-2.18002474518665)+(MCF2L*-1.04530832774745)+(SOX15*-0.992224884476285)+(BEX1*-0.907587628768039)+
(CDO1*-0.871822526242265)+(TRPV1*-1.46645982768831)+(COQ10A*-1.68798125350889)+(TFAP2C*0.967832606083715)+
(CAV1*0.722841430780254)+(SNX5*1.35635732231724)+(PARP6*-1.36262350508992)+(SQOR*1.33174127872047)+(SLC39A14*1.16349465313119)+
(MAP1LC3A*-0.860689113835478)+(CEMIP*-0.59005755740291)+(TXN*0.94085270257141)+(PCDHB14*-0.811693695123737)+(TFR2*-0.798241687760241)+
(USF2*-1.47385325409416)+(ATF4*0.956404022419562)+(PHGDH*0.483156894698922)+(PPARA*-0.923685758091893)+(G3BP1*1.17800387032181)+
(PRR5*-0.539547152725295)+(CARS1*1.17046928959313)+(PRDX4*0.703561532331052)+(CCT3*1.36836908075737)+(TSC1*-0.710935716498839)+
(CBS*0.553452384967907)+(CLOCK*-1.20946536046247)+(VAMP2*-0.757846020766605)+(IL6*0.577234292149136)+(TIMM9*0.736046933599977)+
(TCF4*-0.754520322329335)+(STC1*0.589174240198105)+(ADAM23*-0.615428764368122)+(HPX*-0.716302803991269)+(KIF20A*0.717871874916388)+
(FAM98A*0.999125573497597)+(SMG9*-1.15882491323182)+(CAPRIN2*-0.852103896784302)+(MTDH*0.90965683321584)+(FXR1*0.978687677829785)+
(AHCY*0.787735924986664)+(CUL9*-0.710689105355085)



write.csv(cross,file="cross_sig.csv",row.names=F) 