#################################################################
# @Program: cMonkeyWrapper.py                                   #
# @Version: 1                                                   #
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
# Copyrighted by Chris Plaisier  12/4/2009                      #
#################################################################

from bicluster import bicluster
import os
from subprocess import *
from copy import deepcopy
from sys import stdout

# A class designed to hold the information from a cMonkey RData object
# to facilitate downstream analyses.
#
# Variables:
# RDataFile - where the cmonkey-run-hsa.RData object is stored
# biclusts - the data for all the biclusters (a list of biclust objects)
#
# Functions:
#
class cMonkeyWrapper:
    # Initialize the pssm
    def __init__(self, RDataFile,maxEValue='NA'):
        #
        print 'Beginning to load '+RDataFile+'...'
        # Make the directory to store the data
        if not os.path.exists('biclust'):
            os.mkdir('biclust')
        if not os.path.exists('tmp/meme'):
            os.makedirs('tmp/meme')
        # Build an R script to gather the data from the cMonkey run
        sendToR = []
        # First load up the cmonkey run
        sendToR.append('load(\''+RDataFile+'\')')
        sendToR.append('library(cMonkey)')
        #sendToR.append('cm.attach()')
        sendToR.append('write.csv(e$ratios$raitos, file=\'ratios.csv\')')
        # Load phenotypes for correlations
        sendToR.append('p1 = read.csv(\'phenotypes.csv\',header=T,row.names=1)')
        sendToR.append('useEm = c(\'SEX.bi\',\'AGE\',\'MODEL.1\',\'SURVIVAL\',\'DEAD\')')
        # Then dump everything needed on a bicluster by bicluster basis
        sendToR.append('ks = e$cmonkey.params$k.clust')
        sendToR.append('write.csv(e$cluster.summary(e.cutoff=NA,nrow.cutoff=NA),file=\'cluster.summary.csv\')')
        sendToR.append('write.csv(ks,\'biclust/k.csv\')')
        # Code to get each bicluster information
        sendToR.append('for(k in 1:ks) {')
        sendToR.append('dir.create(paste(\'biclust/\',k,sep=\'\'))') #,showWarnings=FALSE)')
        sendToR += self.getBiclusterRCode()
        # Code to get each PSSM
        sendToR.append('dir.create(paste(\'biclust/\',k,\'/upstream\',sep=\'\'))') #,showWarnings=FALSE)')
        sendToR.append('dir.create(paste(\'biclust/\',k,\'/3pUTR\',sep=\'\'))') #,showWarnings=FALSE)')
        sendToR += self.getPssmRCode(maxEValue)
        sendToR.append('}')
        # Then dump the sequence data for the upstream regions
        sendToR.append('write.csv(e$genome.info$all.upstream.seqs,\'all_upstream_seqs.csv\')')
        # Then dump the sequence data for the 3' UTRs
        sendToR.append('write.csv(e$genome.info$genome.seqs.p3utr,\'all_p3utr_seqs.csv\')')
        # Then dump the nucleotide frequencies for meme 3 file formats
        sendToR.append('write.csv(rbind(c(e$genome.info$bg.list[[\'upstream\']]$A,e$genome.info$bg.list[[\'upstream\']]$G),c(e$genome.info$bg.list[[\'p3utr\']]$A,e$genome.info$bg.list[[\'p3utr\']]$G)),\'nucFreqs.csv\')')
        # Get bg.lists
        sendToR.append('write.table(as.numeric(e$genome.info$bg.list[[\'upstream\']])[2:length(e$genome.info$bg.list[[\'upstream\']])], row.names=names(e$genome.info$bg.list[[\'upstream\']])[2:length(e$genome.info$bg.list[[\'upstream\']])], col.names=\'# 4th order Markov background model\', quote=F, file=\'tmp/meme/bgFile.meme\')')
        sendToR.append('write.table(as.numeric(e$genome.info$bg.list[[\'p3utr\']])[2:length(e$genome.info$bg.list[[\'p3utr\']])], row.names=names(e$genome.info$bg.list[[\'p3utr\']])[2:length(e$genome.info$bg.list[[\'p3utr\']])], col.names=\'# 4th order Markov background model\', quote=F, file=\'tmp/meme/bgFile3pUTR.meme\')')
        ######################
        #!  Run the R code  !#
        print 'Accessing RData object.'
        stderrFile = open('stderr.out','w')
        rProcess = Popen('R --no-save',shell=True,stdin=PIPE,stdout=PIPE,stderr=stderrFile)
        stdout_val = rProcess.communicate('\n'.join(sendToR))[0]
        stderrFile.close()
        #!  Done!           !#
        ######################
        # Build the list of biclusters
        kFile = open('biclust/k.csv','r')
        kFile.readline() # Skip header
        ks = int(kFile.readline().strip().split(',')[1])
        print 'Found '+str(ks)+' clusters.'
        self.biclusters = {}
        for k in range(1,ks+1):
            self.biclusters[k] = bicluster(k)
            if k%10==0:
                stdout.write(str(k))
            else:
                stdout.write('.')
            stdout.flush()
        # Now read in the upstream sequences
        upstreamSeqsFile = open('all_upstream_seqs.csv','r')
        upstreamSeqsFile.readline() # Skip header
        self.seqsUpstream = {}
        for line in upstreamSeqsFile.readlines():
            splitUp = line.strip().split(',')
            self.seqsUpstream[splitUp[0].strip('"')] = splitUp[1].strip('"')
        upstreamSeqsFile.close()
        # Now read in the 3' UTR sequences
        p3utrSeqsFile = open('all_p3utr_seqs.csv','r')
        p3utrSeqsFile.readline() # Skip header
        self.seqs3pUTR = {}
        for line in p3utrSeqsFile.readlines():
            splitUp = line.strip().split(',')
            self.seqs3pUTR[splitUp[0].strip('"')] = splitUp[1].strip('"')
        p3utrSeqsFile.close()
        # Now read in nucleotide frequencies
        nucFreqsFile = open('nucFreqs.csv','r')
        nucFreqsFile.readline() # Skip header
        upFreq = nucFreqsFile.readline().strip().split(',')
        self.nucFreqsUpstream = {'A':upFreq[1],'C':upFreq[2],'G':upFreq[2],'T':upFreq[1]}
        p3utrFreq = nucFreqsFile.readline().strip().split(',')
        self.nucFreqs3pUTR = {'A':p3utrFreq[1],'C':p3utrFreq[2],'G':p3utrFreq[2],'T':p3utrFreq[1]}
        nucFreqsFile.close()
        print '\nDone loading.\n'
        
    # Return a particular bicluster
    def getBicluster(self,k):
        return deepcopy(self.biclusters[k])
    
    # Return a dictionary of all biclusters
    def getBiclusters(self):
        return deepcopy(self.biclusters)
    
    # Return a list of bicluster names
    def getBiclusterNames(self):
        return deepcopy(self.biclusters.keys())
    
    # Get all Upstream pssms
    def getPssmsUpstream(self,maxScore='NA',maxNormResid='NA',maxEValue='NA',maxSurv='NA'):
        pssmsNames = []
        pssms = []
        for bi in self.biclusters.keys():
            # Temporarily store the PSSMs
            biOk = 0
            if maxScore=='NA' or float(self.biclusters[bi].getScore())<=float(maxScore):
                if maxNormResid=='NA' or float(self.biclusters[bi].getNormResidual())<=float(maxNormResid):
                    if maxSurv=='NA' or float(self.biclusters[bi].getSurvival()['"Survival"']['pValue'])<=float(maxSurv):
                        biOk = 1
            if biOk==1:
                tmpPssms = self.biclusters[bi].getPssmsUpstream()
                for pssm in tmpPssms:
                    # Only add it if it is less than an E-Value threshold
                    if maxEValue=='NA' or float(pssm.getEValue())<=float(maxEValue):
                        pssms.append(deepcopy(pssm))
                        pssmsNames.append(pssm.getName())
        return dict(zip(pssmsNames,pssms))
    
    # Get all 3' UTR pssms
    def getPssms3pUTR(self,maxScore='NA',maxNormResid='NA',maxEValue='NA',maxSurv='NA'):
        pssmsNames = []
        pssms = []
        for bi in self.biclusters.keys():
            # Temporarily store the PSSMs
            biOk = 0
            if maxScore=='NA' or float(self.biclusters[bi].getScore())<=float(maxScore):
                if maxNormResid=='NA' or float(self.biclusters[bi].getNormResidual())<=float(maxNormResid):
                    if maxSurv=='NA' or float(self.biclusters[bi].getSurvival()['"Survival"']['pValue'])<=float(maxSurv):
                        biOk = 1
            if biOk==1:
                tmpPssms = self.biclusters[bi].getPssms3pUTR()
                for pssm in tmpPssms:
                    # Only add it if it is less than an E-Value threshold
                    if maxEValue=='NA' or float(pssm.getEValue())<=float(maxEValue):
                        pssms.append(deepcopy(pssm))
                        pssmsNames.append(pssm.getName())
        return dict(zip(pssmsNames,pssms))
    
    # getSeqsUpstream() - returns the upstream sequences
    def getSeqsUpstream(self):
        return deepcopy(self.seqsUpstream)
    
    # getSeqs3pUTR() - returns the 3' UTR sequences
    def getSeqs3pUTR(self):
        return deepcopy(self.seqs3pUTR)

    # getNucFreqsUpstream() - retunres the 
    def getNucFreqsUpstream(self):
        return deepcopy(self.nucFreqsUpstream)
    
    # getNucFreqs3pUTR() - returns the nucleotide frequencies for the 3pUTR
    def getNucFreqs3pUTR(self):
        return deepcopy(self.nucFreqs3pUTR)
    
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
        sendToR.append('sex1 = cor.test(mu.1,d1[,\'SEX.bi\'])')
        sendToR.append('age1 = cor.test(mu.1,d1[,\'AGE\'])')
        sendToR.append('model1 = cor.test(mu.1,d1[,\'MODEL.1\'])')
        sendToR.append('c1 = rbind(c(sex1$estimate,sex1$p.value),c(age1$estimate,age1$p.value),c(model1$estimate,model1$p.value))')
        sendToR.append('rownames(c1) = c(\'Sex\',\'Age\',\'Model.1\')')
        sendToR.append('library(survival)')
        sendToR.append('d2 = data.frame(mu.1,d1)')
        sendToR.append('cph1 = summary(coxph(Surv(SURVIVAL,DEAD) ~ mu.1, d2))')
        sendToR.append('cph2 = summary(coxph(Surv(SURVIVAL,DEAD) ~ mu.1 + AGE, d2))')
        sendToR.append('cph3 = summary(coxph(Surv(SURVIVAL,DEAD) ~ mu.1 + AGE + SEX.bi, d2))')
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
        sendToR.append('write.table(rbind(c(e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$e,e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$sites,NA,NA),cbind(c(e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$posns[,1]),rep(NA,e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$sites),rep(NA,e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$sites),rep(NA,e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$sites)),as.matrix(e$meme.scores[[\'upstream\']][[k]]$meme.out[[1]]$pssm)),file=paste(\'biclust/\',k,\'/upstream/motif1.csv\',sep=\'\'),col.names=F,row.names=F,sep=\',\')')
        sendToR.append('}')
        # Second upstream motif
        sendToR.append('e2 = try(e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$e.value,TRUE)')
        if maxEValue=='NA':
            sendToR.append('if(!class(e2)==\'try-error\' && !is.null(e2)) {')
        else:
            sendToR.append('if(!class(e2)==\'try-error\' && !is.null(e2) && e2<='+str(maxEValue)+') {')
        sendToR.append('write.table(rbind(c(e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$e,e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$sites,NA,NA),cbind(c(e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$posns[,1]),rep(NA,e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$sites),rep(NA,e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$sites),rep(NA,e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$sites)),as.matrix(e$meme.scores[[\'upstream\']][[k]]$meme.out[[2]]$pssm)),file=paste(\'biclust/\',k,\'/upstream/motif2.csv\',sep=\'\'),col.names=F,row.names=F,sep=\',\')')
        sendToR.append('}')
        # First 3' UTR motif
        sendToR.append('e3 = try(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$e.value,TRUE)')
        if maxEValue=='NA':
            sendToR.append('if(!class(e3)==\'try-error\' && !is.null(e3)) {')
        else:
            sendToR.append('if(!class(e3)==\'try-error\' && !is.null(e3) && e3<='+str(maxEValue)+') {')
        sendToR.append('write.table(rbind(c(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$e,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$sites,NA,NA),cbind(c(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$posns[,1]),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$sites),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$sites),rep(NA,e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$sites)),as.matrix(e$meme.scores[[\'p3utr\']][[k]]$meme.out[[1]]$pssm)),file=paste(\'biclust/\',k,\'/3pUTR/motif1.csv\',sep=\'\'),col.names=F,row.names=F,sep=\',\')')
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

