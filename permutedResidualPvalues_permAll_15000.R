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
d1 = read.csv('cluster.members.genes.txt',header=F)
biclustMembership.gene = list()
allGenes = c()
for(j in 1:length(d1[,1])) {
    biclustMembership.gene[[j]] = strsplit(as.character(d1[j,]),split=' ')[[1]][-1]
}

# Read in genes for each cluster
d1 = read.csv('cluster.members.conditions.txt',header=F)
biclustMembership.cond = list()
allGenes = c()
for(j in 1:length(d1[,1])) {
    biclustMembership.cond[[j]] = strsplit(as.character(d1[j,]),split=' ')[[1]][-1]
}

# Read in expression ratios file
ratios <- read.delim( file=gzfile('../ratios.tsv.gz'), sep="\t", as.is=T, header=T,row.names=1 )
maxRowVar = mean( apply( ratios, 1, var, use="pair" ), na.rm=T )

# Calculate the residuals for all clusters in the second dataset
ks = length(biclustMembership.gene)
outNames = c('n.rows','n.cols','orig.resid','avg.perm.resid','perm.p','orig.resid.norm','avg.norm.perm.resid','norm.perm.p','pc1.var.exp','avg.pc1.var.exp','pc1.perm.p')
m1 = matrix(ncol=length(outNames),nrow=ks,dimnames=list(1:ks,outNames))
permutations = 1/0.05*ks
print(paste('Running ',permutations,' permutations...',sep=''))
p1 = read.csv('phenotypes.csv',header=T,row.names=1)
for(k in 1:ks) {
    # Get and add number of rows and columns
    k.rows = biclustMembership.gene[[k]]
    k.cols = biclustMembership.cond[[k]]
    if(length(k.rows)>1 && length(k.cols)>1) {
        m1[k,1] = length(k.rows)
        m1[k,2] = length(k.cols)
        m1[k,3] = residual(as.matrix(ratios[k.rows,k.cols]))
        sub = sapply(1:permutations, function(i) { residual(as.matrix(ratios[sample(rownames(ratios),m1[k,1]), sample(colnames(ratios),m1[k,2])])) })
        m1[k,4] = mean(sub)
        m1[k,5] = length(which(sub <= m1[k,3]))/permutations
        m1[k,6] = residual.norm(as.matrix(ratios[k.rows,k.cols]),maxRowVar)
        sub = sapply(1:permutations, function(i) { residual.norm(as.matrix(ratios[sample(rownames(ratios),m1[k,1]), sample(colnames(ratios),m1[k,2])]), maxRowVar) })
        m1[k,7] = mean(sub)
        m1[k,8] = length(which(sub <= m1[k,4]))/permutations
        pc1 = princomp(as.matrix(t(ratios[k.rows,k.cols])))
        m1[k,9] = ((pc1$sdev^2)/sum(pc1$sdev^2))[1]
        sub = sapply(1:permutations, function(i) { pc1 = princomp(as.matrix(t(ratios[sample(rownames(ratios),m1[k,1]), sample(colnames(ratios),m1[k,2])]))); return(((pc1$sdev^2)/sum(pc1$sdev^2))[1]) })
        m1[k,10] = mean(sub)
        m1[k,11] = length(which(sub >= m1[k,3]))/permutations
    }
    print(c(k,as.numeric(m1[k,])))
}
write.csv(m1,file='residualPermutedPvalues_permAll.csv')

