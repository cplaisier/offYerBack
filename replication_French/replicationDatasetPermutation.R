# Read in the cMonkey run
#load('../e.resplit_vAll_0.1.RData')
#library(cMonkey)
#update.cmonkey.env(e)
library(survival)
library(hgu133plus2.db)
convertMe = as.list(hgu133plus2ALIAS2PROBE)

residual <- function( rats ) {
  d.rows <- rowMeans( rats, na.rm=T )
  d.cols <- colMeans( rats, na.rm=T )
  d.all <- mean( d.rows, na.rm=T )
  rij <- rats + d.all
  rij <- rij - matrix( d.cols, nrow=nrow( rij ), ncol=ncol( rij ), byrow=T )
  rij <- rij - matrix( d.rows, nrow=nrow( rij ), ncol=ncol( rij ), byrow=F )
  average.r <- mean( abs( rij ), na.rm = TRUE )
  average.r
}

residual.norm <- function( rats, maxRowVar ) {
  d.rows <- rowMeans( rats, na.rm=T )
  d.cols <- colMeans( rats, na.rm=T )
  d.all <- mean( d.rows, na.rm=T )
  rij <- rats + d.all
  rij <- rij - matrix( d.cols, nrow=nrow( rij ), ncol=ncol( rij ), byrow=T )
  rij <- rij - matrix( d.rows, nrow=nrow( rij ), ncol=ncol( rij ), byrow=F )
  average.r <- mean( abs( rij ), na.rm = TRUE )
  row.var <- mean( apply( rats, 1, var, use = "pairwise.complete.obs" ), na.rm=T )
  if ( is.na( row.var ) || row.var > maxRowVar ) row.var <- maxRowVar
  average.r <- average.r / row.var
  average.r
}

# Read in genes for each cluster
d1 = read.csv('../output/cluster.members.genes.txt',header=F)
biclustMembership = list()
allGenes = c()
for(j in 1:length(d1[,1])) {
    biclustMembership[[j]] = strsplit(as.character(d1[j,]),split=' ')[[1]][-1]
}

# Read in a second dataset
ratSec <- read.delim( file=gzfile('mergedFrench_groups.csv.gz'), sep=",", as.is=T, header=T,row.names=1 )
#ratSec <- read.delim( file='normalized_AltCDFmsk.csv', sep=",", as.is=T, header=T,row.names=1 )
#ratSec <- read.delim( file='normalized_AltCDFmsk_calls_present_groups.csv', sep=",", as.is=T, header=T,row.names=1 )
rownames( ratSec ) <- toupper( rownames(ratSec) )
biclustMembership.sec = list()
for(j in 1:length(biclustMembership)) {
    biclustMembership.sec[[j]] = intersect(toupper(unlist(convertMe[biclustMembership[[j]]])), rownames(ratSec))
}

ratSec <- as.matrix(ratSec)
controls2 <- read.csv('controls.txt',header=F,as.is=T)[[1]]
rug2 <- read.csv('rug.csv',header=F,as.is=T,row.names=1)
rug2 <- rug2[intersect(rownames(rug2),colnames(ratSec)),]
ratSec.gbm <- ratSec[,c(rownames(rug2)[c(which(rug2[,1]=='GBM (grade IV)'))],controls2)]
ratSec.all <- ratSec[,c(rownames(rug2)[c(which(rug2[,1]=='GBM (grade IV)'),which(rug2[,1]=='OD (grade III)'),which(rug2[,1]=='OA (grade III)'),which(rug2[,1]=='A (grade III)'))],controls2)]
maxRowVar.sec.gbm = mean( apply( ratSec.gbm, 1, var, use="pair" ), na.rm=T )
maxRowVar.sec.all = mean( apply( ratSec.all, 1, var, use="pair" ), na.rm=T )

# Calculate the residuals for all clusters in the second dataset
ks = length(biclustMembership)
outNames = c('n.rows','overlap.rows','new.resid.norm.gbm','avg.norm.perm.resid.gbm','norm.perm.p.gbm','new.resid.norm.all','avg.norm.perm.resid.all','norm.perm.p.all','pc1.var.exp.gbm','avg.pc1.var.exp.gbm','pc1.perm.p.gbm','pc1.var.exp.all','avg.pc1.var.exp.all','pc1.perm.p.all','survival.gbm','survival.p.gbm','survival.age.gbm','survival.age.p.gbm','survival.all','survival.p.all','survival.age.all','survival.age.p.all')
m1 = matrix(ncol=length(outNames),nrow=ks,dimnames=list(1:ks,outNames))
permutations = 1000
p1 = read.csv('phenotypes_French.csv',header=T,row.names=1)
for(k in 1:ks) {
    # Get and add number of rows and columns
    k.rows.sec = biclustMembership.sec[[k]]
    print(k.rows.sec)
    if(length(k.rows.sec)>1) {
        m1[k,1] = length(biclustMembership[[k]])
        m1[k,2] = length(k.rows.sec)
        m1[k,3] = residual.norm(ratSec.gbm[k.rows.sec,],maxRowVar.sec.gbm)
        sub = sapply(1:permutations, function(i) { residual.norm(ratSec.gbm[sample(rownames(ratSec.gbm),m1[k,2]),],maxRowVar.sec.gbm) })
        m1[k,4] = mean(sub)
        m1[k,5] = length(which(sub <= m1[k,3]))/permutations
        m1[k,6] = residual.norm(ratSec.all[k.rows.sec,],maxRowVar.sec.all)
        sub = sapply(1:permutations, function(i) { residual.norm(ratSec.all[sample(rownames(ratSec.all),m1[k,2]),],maxRowVar.sec.all) })
        m1[k,7] = mean(sub)
        m1[k,8] = length(which(sub <= m1[k,6]))/permutations
        tmp.pc.1.gbm = princomp(as.matrix(t(ratSec.gbm[k.rows.sec,])))
        m1[k,9] = ((tmp.pc.1.gbm$sdev^2)/sum(tmp.pc.1.gbm$sdev^2))[1]
        sub = sapply(1:permutations, function(i) { pc1 = princomp(as.matrix(t(ratSec.gbm[sample(rownames(ratSec.gbm),m1[k,2]),]))); return(((pc1$sdev^2)/sum(pc1$sdev^2))[1]) })
        m1[k,10] = mean(sub)
        m1[k,11] = length(which(sub >= m1[k,9]))/permutations
        tmp.pc.1.all = princomp(as.matrix(t(ratSec.all[k.rows.sec,])))
        m1[k,12] = ((tmp.pc.1.all$sdev^2)/sum(tmp.pc.1.all$sdev^2))[1]
        sub = sapply(1:permutations, function(i) { pc1 = princomp(as.matrix(t(ratSec.all[sample(rownames(ratSec.all),m1[k,2]),]))); return(((pc1$sdev^2)/sum(pc1$sdev^2))[1]) })
        m1[k,13] = mean(sub)
        m1[k,14] = length(which(sub >= m1[k,12]))/permutations
        # Survival analysis
        pc.1.gbm = tmp.pc.1.gbm$scores[,1]
        pc.1.all = tmp.pc.1.all$scores[,1]
        d2 = data.frame(p1[colnames(ratSec.gbm),],pc.1.gbm)
        scph1 = summary(coxph(Surv(SURVIVAL.YEARS,DEAD=='Dead') ~ pc.1.gbm, data=d2))
        m1[k,15] = scph1$coef[1,4]
        m1[k,16] = scph1$coef[1,5]
        scph2 = summary(coxph(Surv(SURVIVAL.YEARS,DEAD=='Dead') ~ pc.1.gbm + AGE, data=d2))
        m1[k,17] = scph2$coef[1,4]
        m1[k,18] = scph2$coef[1,5]
        d2 = data.frame(p1[colnames(ratSec.all),],pc.1.all)
        scph1 = summary(coxph(Surv(SURVIVAL.YEARS,DEAD=='Dead') ~ pc.1.all, data=d2))
        m1[k,19] = scph1$coef[1,4]
        m1[k,20] = scph1$coef[1,5]
        scph2 = summary(coxph(Surv(SURVIVAL.YEARS,DEAD=='Dead') ~ pc.1.all + AGE, data=d2))
        m1[k,21] = scph2$coef[1,4]
        m1[k,22] = scph2$coef[1,5]
    }
}
write.csv(m1,file='replicationPvalues.csv')

