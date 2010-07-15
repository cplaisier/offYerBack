# Read in the cMonkey run
load('../iter3000.RData')
library(cMonkey)
update.cmonkey.env(e)
library(survival)

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

# Read in a second dataset
ratSec <- read.delim( file=gzfile( 'merged_French.csv.gz' ), sep=",", as.is=T, header=T,row.names=1 )
rownames( ratSec ) <- toupper( rownames(ratSec) )
ratSec <- as.matrix(ratSec)
controls2 <- read.csv('controls.txt',header=F,as.is=T)[[1]]
rug2 <- read.csv('rug.csv',header=F,as.is=T,row.names=1)
rug2 <- rug2[intersect(rownames(rug2),colnames(ratSec)),]
ratSec <- ratSec[,c(rownames(rug2)[c(which(rug2[,1]=='GBM (grade IV)'),which(rug2[,1]=='OD (grade III)'),which(rug2[,1]=='OA (grade III)'),which(rug2[,1]=='A (grade III)'))],controls2)]
groups2 <- rep(1, length(colnames(ratSec)))
maxRowVar.sec = mean( apply( ratSec, 1, var, use="pair" ), na.rm=T )

# Calculate the residuals for all clusters in the second dataset
ks = e$cmonkey.params$k.clust
outNames = c('n.rows','orig.resid','orig.resid.norm','overlap.rows','new.resid','avg.perm.resid','perm.p','new.resid.norm','avg.norm.perm.resid','norm.perm.p','survival','survival.p','survival.age','survival.age.p')
m1 = matrix(ncol=length(outNames),nrow=ks,dimnames=list(1:ks,outNames))
permutations = 1000
p1 = read.csv('phenotypes_French.csv',header=T,row.names=1)
for(k in 1:ks) {
    # Get and add number of rows and columns
    k.rows = e$get.rows(k)
    k.cols = e$get.cols(k)
    if(length(k.rows)>1) {
        m1[k,1] = length(k.rows)
        m1[k,2] = e$cluster.resid(k,varNorm=F)
        m1[k,3] = e$cluster.resid(k,varNorm=T)
        k.rows.sec = k.rows[c(k.rows %in% rownames(ratSec))]
        m1[k,4] = length(k.rows.sec)
        m1[k,5] = residual(ratSec[k.rows.sec,])
        sub = sapply(1:permutations, function(i) { residual(ratSec[sample(rownames(ratSec),m1[k,4]),]) })
        m1[k,6] = mean(sub)
        m1[k,7] = length(which(sub <= m1[k,5]))/permutations
        m1[k,8] = residual.norm(ratSec[k.rows.sec,],maxRowVar.sec)
        sub = sapply(1:permutations, function(i) { residual.norm(ratSec[sample(rownames(ratSec),m1[k,4]),],maxRowVar.sec) })
        m1[k,9] = mean(sub)
        m1[k,10] = length(which(sub <= m1[k,8]))/permutations
        mu.1 = apply(ratSec[k.rows.sec,],2,median)
        # Survival analysis
        d2 = data.frame(p1[names(mu.1),],mu.1)
        scph1 = summary(coxph(Surv(SURVIVAL.YEARS,DEAD=='Dead') ~ mu.1, data=d2))
        m1[k,11] = scph1$coef[1,4]
        m1[k,12] = scph1$coef[1,5]
        scph2 = summary(coxph(Surv(SURVIVAL.YEARS,DEAD=='Dead') ~ mu.1 + AGE, data=d2))
        m1[k,13] = scph2$coef[1,4]
        m1[k,14] = scph2$coef[1,5]
    }
}
write.csv(m1,file='replicationPvalues.csv')

