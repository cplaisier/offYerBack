#################################################################
# @Program: cMonkeyWrapper.py                                   #
# @Version: 2 (python-cMonkey)                                  #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  2/18/2013                      #
#################################################################

import os
import sqlite3 as lite
import gzip
from bicluster import bicluster
from subprocess import *
from sys import stdout

# A class designed to hold the information from a cMonkey RData object
# to facilitate downstream analyses.
#
# Variables:
# biclusters = dictionary
# seqsUpstream = dictionary
# seqs3pUTR = dictionary
# nucFreqsUpstream = dictionary
# nucFreqs3pUTR = dictionary
#
# Functions:
# getBicluster(k)
# getBiclusters()
# getBiclusterNames()
# getPssmsUpstream(maxScore='NA',maxNormResid='NA',maxEValue='NA',maxSurv='NA')
# getPssms3pUTR(maxScore='NA',maxNormResid='NA',maxEValue='NA',maxSurv='NA')
# getSeqsUpstream() - returns the upstream sequences
# getSeqs3pUTR()
# getNucFreqsUpstream()
# getNucFreqs3pUTR()
# getBiclusterRCode()
# getPssmRCode(maxEValue)
# getBiclusterSeqsUpstream(k)
# getBiclusterSeqs3pUTR(k)
class cMonkeyWrapper:
    # Initialize the pssm
    def __init__(self, sqliteDb, maxEValue='NA', meme_upstream=0, weeder_upstream=0, weeder_3pUTR=0, pita_3pUTR=0, targetscan_3pUTR=0):
        # What has been run on this cMonkey run, legend = [0: not run, 1: run]
        de_novo_method_upstream = None
        de_novo_method_3pUTR = None
        if meme_upstream==1 and weeder_upstream==1:
            raise RuntimeError('You trained the same run on both MEME and Weeder! Are you stupid or something?')
        elif meme_upstream==1:
            de_novo_method_upstream = 'meme'
        elif weeder_upstream==1:
            de_novo_method_upstream = 'weeder'
        self.meme_upstream = meme_upstream
        self.weeder_upstream = weeder_upstream
        if weeder_3pUTR==1:
            de_novo_method_3pUTR = 'weeder'
        self.weeder_3pUTR = weeder_3pUTR
        self.pita_3pUTR = pita_3pUTR
        self.targetscan_3pUTR = targetscan_3pUTR
        # Attach to the database
        con = lite.connect(sqliteDb)
        con.row_factory = lite.Row
        cur = con.cursor()
        # Get the number of biclusters in run
        cur.execute('SELECT * FROM run_infos')
        data = cur.fetchall()
        con.close()
        ks = data[0]['num_clusters']
        print 'Found '+str(ks)+' clusters.'
        self.biclusters = {}
        for k in range(1,ks+1):
            self.biclusters[k] = bicluster(k, de_novo_method_upstream=de_novo_method_upstream, de_novo_method_3pUTR=de_novo_method_3pUTR, sqliteDb=sqliteDb)
            if k%10==0:
                stdout.write(str(k))
            else:
                stdout.write('.')
            stdout.flush()
        # Now read in the upstream sequences
        upstreamSeqsFile = gzip.open('seqs/promoterSeqs_Homo_sapiens_gs_fp.csv.gz','rb')
        upstreamSeqsFile.readline() # Skip header
        self.seqsUpstream = {}
        for line in upstreamSeqsFile.readlines():
            splitUp = line.strip().split(',')
            self.seqsUpstream[splitUp[0].strip('"')] = splitUp[1].strip('"')
        upstreamSeqsFile.close()
        # Now read in the 3' UTR sequences
        p3utrSeqsFile = gzip.open('seqs/p3utrSeqs_Homo_sapiens_gs.csv.gz','rb')
        p3utrSeqsFile.readline() # Skip header
        self.seqs3pUTR = {}
        for line in p3utrSeqsFile.readlines():
            splitUp = line.strip().split(',')
            self.seqs3pUTR[splitUp[0].strip('"')] = splitUp[1].strip('"')
        p3utrSeqsFile.close()
        # Now read in nucleotide frequencies
        nucFreqsFile = open('seqs/nucFreqs.csv','r')
        nucFreqsFile.readline() # Skip header
        upFreq = nucFreqsFile.readline().strip().split(',')
        self.nucFreqsUpstream = {'A':upFreq[1],'C':upFreq[2],'G':upFreq[2],'T':upFreq[1]}
        p3utrFreq = nucFreqsFile.readline().strip().split(',')
        self.nucFreqs3pUTR = {'A':p3utrFreq[1],'C':p3utrFreq[2],'G':p3utrFreq[2],'T':p3utrFreq[1]}
        nucFreqsFile.close()
        # Close database connection
        con.close()
        print '\nDone loading.\n'

    # Return a particular bicluster
    def getBicluster(self,k):
        return self.biclusters[k]

    # Return a dictionary of all biclusters
    def getBiclusters(self):
        return self.biclusters

    # Return a list of bicluster names
    def getBiclusterNames(self):
        return self.biclusters.keys()

    # Get all Upstream pssms
    def getPssmsUpstream(self,maxNormResid='NA',maxEValue='NA',maxSurv='NA',de_novo_method='NA'):
        pssmsNames = []
        pssms = []
        for bi in self.biclusters.keys():
            # Temporarily store the PSSMs
            biOk = 0
            if maxNormResid=='NA' or float(self.biclusters[bi].getNormResidual())<=float(maxNormResid):
                if maxSurv=='NA' or float(self.biclusters[bi].getSurvival()['"Survival"']['pValue'])<=float(maxSurv):
                    biOk = 1
            if biOk==1:
                tmpPssms = self.biclusters[bi].getPssmsUpstream()
                for pssm in tmpPssms:
                    if de_novo_method=='NA' or de_novo_method==pssm.getMethod():
                        # Only add it if it is less than an E-Value threshold
                        if maxEValue=='NA' or float(pssm.getEValue())<=float(maxEValue):
                            pssms.append(pssm)
                            pssmsNames.append(pssm.getName())
        return dict(zip(pssmsNames,pssms))

    # Get all 3' UTR pssms
    def getPssms3pUTR(self,maxNormResid='NA',maxEValue='NA',maxSurv='NA',de_novo_method='NA'):
        pssmsNames = []
        pssms = []
        for bi in self.biclusters.keys():
            # Temporarily store the PSSMs
            biOk = 0
            if maxNormResid=='NA' or float(self.biclusters[bi].getNormResidual())<=float(maxNormResid):
                if maxSurv=='NA' or float(self.biclusters[bi].getSurvival()['"Survival"']['pValue'])<=float(maxSurv):
                    biOk = 1
            if biOk==1:
                tmpPssms = self.biclusters[bi].getPssms3pUTR()
                for pssm in tmpPssms:
                    # Only add it if it is less than an E-Value threshold
                    if de_novo_method=='NA' or de_novo_method==pssm.getMethod():
                        if maxEValue=='NA' or float(pssm.getEValue())<=float(maxEValue):
                            pssms.append(pssm)
                            pssmsNames.append(pssm.getName())
        return dict(zip(pssmsNames,pssms))

    # getSeqsUpstream() - returns the upstream sequences
    def getSeqsUpstream(self):
        return self.seqsUpstream

    # getSeqs3pUTR() - returns the 3' UTR sequences
    def getSeqs3pUTR(self):
        return self.seqs3pUTR

    # getBiclusterSeqsUpstream() - returns the upstream sequences for a bicluster as a dictionary of {<gene_name>: <seqeunce>, ...}
    def getBiclusterSeqsUpstream(self, k):
        genes = self.biclusters[k].getGenes()
        seqs = dict(zip([gene for gene in genes if gene in self.seqsUpstream], [self.seqsUpstream[gene] for gene in genes if gene in self.seqsUpstream]))
        return seqs

    # getBiclusterSeqs3pUTR() - returns the 3' UTR sequences for a bicluster as a dictionary of {<gene_name>: <seqeunce>, ...}
    def getBiclusterSeqs3pUTR(self, k):
        genes = self.biclusters[k].getGenes()
        seqs = dict(zip([gene for gene in genes if gene in self.seqs3pUTR], [self.seqs3pUTR[gene] for gene in genes if gene in self.seqs3pUTR]))
        return seqs

    # getNucFreqsUpstream() - retunres the 
    def getNucFreqsUpstream(self):
        return self.nucFreqsUpstream

    # getNucFreqs3pUTR() - returns the nucleotide frequencies for the 3pUTR
    def getNucFreqs3pUTR(self):
        return self.nucFreqs3pUTR

    # The R code to get what is needed for a pssm needs to be run in a loop for k, where k = a cluster number
    # useEm = c('SEX.INF.bi','AGE','MODEL.1','MODEL.2','MODEL.3','MODEL.4','MODEL.5','SURVIVAL')
    def getBiclusterRCode(self):
        # First upstream motif
        sendToR = []
        sendToR.append('k.rows = e$get.rows(k)')
        sendToR.append('write.csv(k.rows,paste(\'biclust/\',k,\'/genes.csv\',sep=\'\'))')
        sendToR.append('k.cols = e$get.cols(k)')
        sendToR.append('write.csv(k.cols,paste(\'biclust/\',k,\'/conditions.csv\',sep=\'\'))')
        sendToR.append('resid = e$cluster.resid(k,varNorm=F)')
        sendToR.append('residNorm = e$cluster.resid(k,varNorm=T)')
        sendToR.append('write.csv(c(resid,residNorm),paste(\'biclust/\',k,\'/resid.csv\',sep=\'\'))')
        # Correlate with traits
        sendToR.append('if(length(k.rows)>1) {')
        sendToR.append('    mu.1 = apply(e$ratios$ratios[k.rows,k.cols],2,median)')
        sendToR.append('} else {')
        sendToR.append('    mu.1 = e$ratios$ratios[k.rows,k.cols]')
        sendToR.append('}')
        sendToR.append('d1 = p1[k.cols,useEm]')
        sendToR.append('sex1 = try(cor.test(mu.1,d1[,\'SEX.bi\']), TRUE)')
        sendToR.append('if(class(sex1)==\'try-error\') { sex1 = c(); sex1$estimate = NA; sex1$p.value = NA }')
        sendToR.append('age1 = try(cor.test(mu.1,d1[,\'AGE\']), TRUE)')
        sendToR.append('if(class(age1)==\'try-error\') { age1 = c(); age1$estimate = NA; age1$p.value = NA }')
        sendToR.append('c1 = rbind(c(sex1$estimate,sex1$p.value),c(age1$estimate,age1$p.value))')
        sendToR.append('rownames(c1) = c(\'Sex\',\'Age\')')
        sendToR.append('library(survival)')
        sendToR.append('d2 = data.frame(mu.1,d1)')
        sendToR.append('cph1 = try(summary(coxph(Surv(SURVIVAL,DEAD==\'DEAD\') ~ mu.1, d2)), TRUE)')
        sendToR.append('if(class(cph1)==\'try-error\') { cph1 = c(); cph1$coef = matrix(nrow=2,ncol=5); }')
        sendToR.append('cph2 = try(summary(coxph(Surv(SURVIVAL,DEAD==\'DEAD\') ~ mu.1 + AGE, d2)),TRUE)')
        sendToR.append('if(class(cph2)==\'try-error\') { cph2 = c(); cph2$coef = matrix(nrow=2,ncol=5); }')
        sendToR.append('cph3 = try(summary(coxph(Surv(SURVIVAL,DEAD==\'DEAD\') ~ mu.1 + AGE + SEX.bi, d2)),TRUE)')
        sendToR.append('if(class(cph3)==\'try-error\') { cph3 = c(); cph3$coef = matrix(nrow=2,ncol=5); }')
        sendToR.append('s1 = rbind(c(cph1$coef[1,4],cph1$coef[1,5]),c(cph2$coef[1,4],cph2$coef[1,5]),c(cph3$coef[1,4],cph3$coef[1,5]))')
        sendToR.append('rownames(s1) = c(\'Survival\',\'Survival.Age\',\'Survival.Age.Sex\')')
        sendToR.append('write.csv(c1,paste(\'biclust/\',k,\'/correlation.csv\',sep=\'\'))')
        sendToR.append('write.csv(s1,paste(\'biclust/\',k,\'/survival.csv\',sep=\'\'))')
        return sendToR

    # The R code to get what is needed for a pssm needs to be run in a loop for k, where k = a cluster number
    def getPssmRCode(self,maxEValue):
        # First upstream motif
        sendToR = []
        sendToR.append('e1 = try(e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$e.value,TRUE)')
        if maxEValue=='NA':
            sendToR.append('if(!class(e1)==\'try-error\' && !is.null(e1)) {')
        else:
            sendToR.append('if(!class(e1)==\'try-error\' && !is.null(e1) && e1<='+str(maxEValue)+') {')
        sendToR.append('tmp1 = unique(e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$posns[,1])')
        sendToR.append('write.table(rbind(c(e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$e,e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$sites,NA,NA),cbind(c(tmp1),rep(NA,length(tmp1)),rep(NA,length(tmp1)),rep(NA,length(tmp1))),as.matrix(e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$pssm)),file=paste(\'biclust/\',k,\'/upstream/motif1.csv\',sep=\'\'),col.names=F,row.names=F,sep=\',\')')
        sendToR.append('}')
        # Second upstream motif
        sendToR.append('e2 = try(e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$e.value,TRUE)')
        if maxEValue=='NA':
            sendToR.append('if(!class(e2)==\'try-error\' && !is.null(e2)) {')
        else:
            sendToR.append('if(!class(e2)==\'try-error\' && !is.null(e2) && e2<='+str(maxEValue)+') {')
        sendToR.append('tmp1 = unique(e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$posns[,1])')
        sendToR.append('write.table(rbind(c(e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$e,e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$sites,NA,NA),cbind(c(tmp1),rep(NA,length(tmp1)),rep(NA,length(tmp1)),rep(NA,length(tmp1))),as.matrix(e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$pssm)),file=paste(\'biclust/\',k,\'/upstream/motif2.csv\',sep=\'\'),col.names=F,row.names=F,sep=\',\')')
        sendToR.append('}')
        # First 3' UTR motif
        sendToR.append('e3 = try(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$e.value,TRUE)')
        if maxEValue=='NA':
            sendToR.append('if(!class(e3)==\'try-error\' && !is.null(e3)) {')
        else:
            sendToR.append('if(!class(e3)==\'try-error\' && !is.null(e3) && e3<='+str(maxEValue)+') {')
        sendToR.append('write.table(rbind(c(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$e,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$sites,NA,NA),cbind(c(unique(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$posns[,1])),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$sites),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$sites),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$sites)),as.matrix(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$pssm)),file=paste(\'biclust/\',k,\'/3pUTR/motif1.csv\',sep=\'\'),col.names=F,row.names=F,sep=\',\')')
        sendToR.append('}')
        # Second 3' UTR motif
        sendToR.append('e4 = try(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[2]]$e.value,TRUE)')
        if maxEValue=='NA':
            sendToR.append('if(!class(e4)==\'try-error\' && !is.null(e4)) {')
        else:
            sendToR.append('if(!class(e4)==\'try-error\' && !is.null(e4) && e4<='+str(maxEValue)+') {')
        sendToR.append('write.table(rbind(c(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[2]]$e,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[2]]$sites,NA,NA),cbind(c(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[2]]$posns[,1]),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[2]]$sites),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[2]]$sites),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[2]]$sites)),as.matrix(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[2]]$pssm)),file=paste(\'biclust/\',k,\'/3pUTR/motif2.csv\',sep=\'\'),col.names=F,row.names=F,sep=\',\')')
        sendToR.append('}')
        return sendToR

