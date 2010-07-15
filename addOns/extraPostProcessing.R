# Load cMonkey data from a run
load('cmonkey-run-hsa.RData')
library(cMonkey)
cm.attach()

# Calculate the association scores
ks = cmonkey.env$cmonkey.params$k.clust
k.cols = sort(cMonkey:::get.cols(1))
k.cols = k.cols[-which(k.cols=="HF1708.CEL")]
m1 = matrix(ncol=length(k.cols),nrow=ks,dimnames=list(1:ks,k.cols))
resids = c(1:ks)
for(k in 1:ks) {
    k.rows = cMonkey:::get.rows(k)
    mu.1 = apply(ratios[k.rows,k.cols],2,median)
    m1[k,] = mu.1
    resids[k] = cMonkey:::cluster.resid(k,varNorm=T)
}

scores = cMonkey:::cluster.summary(e.cutoff=NA,nrow.cutoff=NA)
scores = scores[order(as.numeric(rownames(scores))),2]

p1 = read.csv('phenotypes.csv',header=T,row.names=1)
useEm = c('SEX.bi','SEX.INF.bi','AGE','GRADE.NUM','MODEL.1','MODEL.2','MODEL.3','SURVIVAL','OLIGO','ASTRO','GBM','CANCER')

plot(hclust(dist(t(m1))))

library(gplots)
#md1 = ifelse(is.na(p1[k.cols,'MODEL.1']),'white',ifelse(p1[k.cols,'MODEL.1']==0,'green',ifelse(p1[k.cols,'MODEL.1']==1,'red',ifelse(p1[k.cols,'MODEL.1']==2,'black','white'))))
md1 = ifelse(is.na(p1[k.cols,'GRADE.NUM']),'white',ifelse(p1[k.cols,'GRADE.NUM']==0,'green',ifelse(p1[k.cols,'GRADE.NUM']==2,'yellow',ifelse(p1[k.cols,'GRADE.NUM']==3,'red',ifelse(p1[k.cols,'GRADE.NUM']==4,'black',ifelse(p1[k.cols,'MODEL.1']==1,'orange',ifelse(p1[k.cols,'MODEL.1']==2,'darkgrey','white')))))))
na1 = which(is.na(p1[k.cols,'MODEL.1']))
pdf(height=8.5,width=11,file='/local/www/html/motDifLeng/SeS/medExpHeatMap_1Kbp_SES_nobi_norc.pdf')
heatmap.2(m1,dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1)
na2 = which(resids<0.4)
heatmap.2(-m1[-na2,],dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1)
dev.off()

par(mfrow=c(1,2))
cors = c(1:ks)
pVals = c(1:ks)
for(i in 1:ks) {
    c1 = cor.test(m1[i,],p1[k.cols,'MODEL.1'])
    cors[i] = c1$estimate
    pVals[i] = c1$p.value
}
which(pVals<=(0.05/150))
plot(-log10(pVals))
abline(h=-log10(0.05/150))

pValsSLM1 = c(1:ks)
for(i in 1:ks) {
    slm1 = summary.lm(aov(p1[k.cols,'MODEL.1'] ~ m1[i,]))
    pValsSLM1[i] = slm1$coef[2,4]
}
plot(-log10(pValsSLM1))
abline(h=-log10(0.05/150))

plot(-log10(pVals),-log10(pValsSLM1))

library(survival)
m2 = matrix(ncol=6,nrow=ks,dimnames=list(colnames(data.frame(t(m1))),c("Surv.e","Surv.p","Surv.age.e","Surv.age.p","Surv.age.sex.e","Surv.age.sex.p")))
for( i in colnames(data.frame(t(m1)))) {
    d2 = data.frame(t(m1),p1[k.cols,])
    cph1 = summary(coxph(formula(paste("Surv(SURVIVAL,DEAD) ~ ",i,sep="")), d2))
    cph2 = summary(coxph(formula(paste("Surv(SURVIVAL,DEAD) ~ ",i," + AGE",sep="")), d2))
    cph3 = summary(coxph(formula(paste("Surv(SURVIVAL,DEAD) ~ ",i," + AGE + SEX.INF.bi",sep="")), d2))
    m2[i,] = c(cph1$coef[1,4],cph1$coef[1,5],cph2$coef[1,4],cph2$coef[1,5],cph3$coef[1,4],cph3$coef[1,5])
}

plot(-log10(m2[,'Surv.p']),-log10(m2[,'Surv.age.p']))
plot(-log10(pVals),-log10(m2[,'Surv.p']))
points(-log10(pVals),-log10(m2[,'Surv.age.p']),col='blue')
points(-log10(pVals),-log10(m2[,'Surv.age.sex.p']),col='red')
abline(h=-log10(0.05/150),col='red')
abline(v=-log10(0.05/150),col='red')

pdf(height=8.5,width=11,file='/local/www/html/motDifLeng/SeS/medExpHeatMap_500bp_SES_nobi_norc_c_all.pdf')
heatmap.2(m1,dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1,main='All')
na2 = which(resids>0.4)
heatmap.2(-m1[-na2,],dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1,main='Residual <= 0.4')
na3 = which(m2[,'Surv.age.p']>=(0.05/150))
heatmap.2(-m1[-na3,],dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1,main='Surv.age.p < (0.05/150)')
na4 = union(na2,na3)
heatmap.2(-m1[-na4,],dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1,main='Resid. <=0.4 and Surv.age.p < (0.05/150)')
na5 = which(pVals>=(0.05/150))
heatmap.2(-m1[-na5,],dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1,main='Model.1 < (0.05/150)')
na6 = union(na2,na5)
heatmap.2(-m1[-na6,],dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1,main='Resid. <=0.4 and Model.1 < (0.05/150)')
na7 = which(scores>0)
heatmap.2(-m1[-na7,],dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1,main='Score <= 0')
na7 = union(na2,na7)
heatmap.2(-m1[-na7,],dendrogram='both',col=redgreen,trace='none',scale='row',ColSideColors=md1,main='Resid. <=0.4 and Score <= 0')
dev.off()


